#ifndef __CBM3D_WIENER_H__
#define __CBM3D_WIENER_H__

#include "bm3d_wiener.h"

/* Implementation of the second step, i.e. Wiener-filtering, of the BM3D denoising method for YUV 4:4:4 images.
 * Note that the implementation only supports YUV 4:4:4 format, that the U/V component has the same size as Y.
 * The memory of the YUV frame is in a planar format, i.e. [(w*h) Y | (w*h) U | (w*h) V].
 * The sigma of the Y/U/V component is usually different from each other,
 * which is relative to the transform matrix from RGB sapce to the YUV space,
 * and the sigma of the R/G/B is usually thought as the same.
 * For example, if Y = a*R + b*G + c*B, we can compute that sigmaY = sqrt(a*a + b*b + c*c) * sigmaR/G/B.
 * However, it's not a big difference between sigmaY and sigmaU or sigmaV for the BT.601 or BT.709,
 * and we can approximately compute sigmaY/U/V = 0.6 * sigmaR/G/B.
 */
class CBM3D_WIE : public BM3D_WIE
{
public:
	CBM3D_WIE(
		int w_,						// width
		int h_,						// height
		int max_sim = 16,			// maximum similar patches
		int psize_ = 8,				// reference patch size
		int pstep_ = 3,				// reference patch step
		int swinrh_ = 16,			// horizontal search window radius
		int ssteph_ = 1,			// horizontal search step
		int swinrv_ = 16,			// vertical search window radius
		int sstepv_ = 1				// vertical search step
	);
	~CBM3D_WIE();

	/* Load a new noisy yuv444 frame and reset the buffers. */
	void load(
		ImageType* org_noisy_yuv,	// pointer of the input noisy yuv444 (planar) frame
		ImageType* org_basic_yuv,	// pointer of the input noisy yuv444 (planar) frame
		int sigmay,					// sigma of the Y component
		DistType max_mdist = 2500,	// maximum mean distance (L2/L1) between the reference and the candidate
		int sigmau = -1,			// sigma of the U component, same as Y if <0
		int sigmav = -1				// sigma of the V component, same as Y if <0
	);

	/* reset the buffers and redirect the processing reference patch to the first one of the noisy image */
	void reset();

	/* Denoise just a line of reference patches and write out the completed rows. */
	int next_line(
		ImageType* clean_yuv		// pointer of output denoised yuv444 (planar) frame
	);

protected:

	ImageType* noisy_yuv[3];
	ImageType* basic_yuv[3];
	PatchType* numerator_yuv[3];
	PatchType* denominator_yuv[3];

	ImageType* refer_noisy_yuv[3];
	ImageType* refer_basic_yuv[3];
	PatchType* numer_yuv[3];
	PatchType* denom_yuv[3];

	PatchType wie_thres[3];
};

#endif

