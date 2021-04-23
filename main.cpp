#include <iostream>
#include <string>
#include <math.h>
#include "cbm3d.h"
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
	int chnl = 3;		// YUV 4:0:0 or 4:4:4

	int sigma = 36;		// same for Y/U/V, can be different
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

	BM3D *denoiser = NULL;
	if (chnl == 1) {
		// used for YUV 4:0:0
		denoiser = new BM3D(w, h, 16, 8, 3, 16, 1, 16, 1); // at present the psize must be 8
	} else {
		// used for YUV 4;4:4
		denoiser = new CBM3D(w, h, 16, 8, 3, 16, 1, 16, 1);
	}

	int frame = 0;
	while (frames < 0 ||frame < frames)
	{
		if (fread(noisy, sizeof(ImageType), w * h * chnl, inf) != w * h * chnl) break;
		cout << "Processing frame " << frame << "..." << endl;

		denoiser->load(noisy, sigma);
		denoiser->run(clean);

		cout << "noisy PSNR: "    << get_psnr(noisy, gt, w*h*chnl, 255) << "    "
			 << "denoised PSNR: " << get_psnr(clean, gt, w*h*chnl, 255) << endl;

		fwrite(clean, sizeof(ImageType), w * h * chnl, ouf);
		cout << "Frame " << frame << " done!" << endl << endl;
		frame++;
	}

	delete denoiser;

	fclose(gtf);
	delete[] gt;

	fclose(inf);
	fclose(ouf);
	delete[] noisy;
	delete[] clean;

	return 0;
}
