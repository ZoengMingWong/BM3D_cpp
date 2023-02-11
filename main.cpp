#include <iostream>
#include <string>
#include <math.h>
#include "cbm3d.h"
#include "cbm3d_wiener.h"
using namespace std;

double get_psnr(ImageType *img1, ImageType *img2, int pixels, ImageType vmax)
{
	double mse = 0;
	double diff;
	for (int i = 0; i < pixels; i++)
	{
		diff = (double)img1[i] - (double)img2[i];
		mse += diff * diff;
	}
	mse /= pixels;
	return (10 * log10((double)vmax * vmax / mse));
}

FILE *openfile(const char *fname, const char *mode)
{
	FILE *f = fopen(fname, mode);
	if (NULL == f)
	{
		cout << "Failed to open: " << fname << endl;
		exit(1);
	}
	return f;
}

int main()
{
	int w = 512, h = 512;
	int chnl = 3;			// YUV 4:0:0 or 4:4:4

	int en_bm3d_step2 = 1;	// enable step2 of bm3d
	int sigma_step1 = 36;	// here same for Y/U/V, can be different
	int sigma_step2 = 25;	// bigger for smoother, usually a little smaller than step1

	int frames = 1;		// frames to process

	// ground truth
	FILE *gtf = openfile("test/yuv444_512x512_lena_gt.yuv", "rb");
	ImageType *gt = new ImageType[w * h * chnl];
	fread(gt, sizeof(ImageType), w * h * chnl, gtf);

	// noisy input and denoised output
	FILE *inf = openfile("test/yuv444_512x512_lena.yuv", "rb");
	FILE *ouf = openfile("test/yuv444_512x512_lena_deno.yuv", "wb");
	
	ImageType *noisy = new ImageType[w * h * chnl];
	ImageType *clean = new ImageType[w * h * chnl];

	/* hard-thresholding denoiser
	 */ 
	BM3D *denoiser = NULL;
	if (chnl == 1) {
		// used for YUV 4:0:0
		denoiser = new BM3D(w, h, 16, 8, 3, 16, 1, 16, 1); // at present the psize must be 8
	} else {
		// used for YUV 4;4:4
		denoiser = new CBM3D(w, h, 16, 8, 3, 16, 1, 16, 1);
	}

	/* wiener-filtering denoiser
	 * Note that Step2 is independent of Step1, except the basic denoised image.
	 * So you can even just process the Y component in Step2, and reuse the U/V result from Step1.
	 */
	BM3D_WIE* denoiser_wie = NULL;
	if (chnl == 1) {
		// used for YUV 4:0:0
		denoiser_wie = new BM3D_WIE(w, h, 32, 8, 3, 16, 1, 16, 1); // at present the psize must be 8
	}
	else {
		// used for YUV 4;4:4
		denoiser_wie = new CBM3D_WIE(w, h, 32, 8, 3, 16, 1, 16, 1);
	}

	int frame = 0;
	while (frames < 0 ||frame < frames)
	{
		if (fread(noisy, sizeof(ImageType), w * h * chnl, inf) != w * h * chnl) break;
		cout << "Processing frame " << frame << "..." << endl;

		// hard thresholding
		denoiser->load(noisy, sigma_step1);
		denoiser->run(clean);

		cout << "noisy PSNR: "    << get_psnr(noisy, gt, w*h*chnl, 255) << "    "
			 << "denoised PSNR: " << get_psnr(clean, gt, w*h*chnl, 255) << endl;

		fwrite(clean, sizeof(ImageType), w * h * chnl, ouf);

		if (en_bm3d_step2)
		{
			// wiener filtering
			denoiser_wie->load(noisy, clean, sigma_step2);
			denoiser_wie->run(clean);

			cout << "noisy PSNR: " << get_psnr(noisy, gt, w * h * chnl, 255) << "    "
				<< "wiener denoised PSNR: " << get_psnr(clean, gt, w * h * chnl, 255) << endl;

			fwrite(clean, sizeof(ImageType), w * h * chnl, ouf);
		}

		cout << "Frame " << frame << " done!" << endl << endl;
		frame++;
	}

	delete denoiser;
	delete denoiser_wie;

	fclose(gtf);
	delete[] gt;

	fclose(inf);
	fclose(ouf);
	delete[] noisy;
	delete[] clean;

	return 0;
}
