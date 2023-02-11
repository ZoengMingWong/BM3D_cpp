#ifndef __BM3D_WIENER_H__
#define __BM3D_WIENER_H__

#include <iostream>
#include <time.h>
#include <omp.h>

#include "global_define.h"
#include "patch_2d.h"
#include "group_3d.h"

// 8x8 Kaiser window
extern const PatchType Kaiser[64];

/* Implementation of the second step, i.e. Wiener-filtering, of the BM3D denoising method for graysacle images.
 * Note that the implementation only supports graysacle images, or YUV 4:0:0 format.
 * If given the sigma of RGB colr space, which is usually thought as the same for R/G/B component,
 * the sigma of Y is relative to the transform matrix from RGB sapce to the YUV space,
 * For example, if Y = a*R + b*G + c*B, we can compute that sigmaY = sqrt(a*a + b*b + c*c) * sigmaR/G/B.
 * We can approximately compute sigmaY = 0.6 * sigmaR/G/B for the BT.601 or BT.709.
*/
class BM3D_WIE
{
public:
	BM3D_WIE(
		int w_,						// width
		int h_,						// height
		int max_sim = 16,			// maximum similar patches
		int psize_  = 8,			// reference patch size
		int pstep_  = 3,			// reference patch step
		int swinrh_ = 16,			// horizontal search window radius
		int ssteph_ = 1,			// horizontal search step
		int swinrv_ = 16,			// vertical search window radius
		int sstepv_ = 1				// vertical search step
		);
	virtual ~BM3D_WIE();

	/* Load a new grayscale image and reset the buffers. */
	virtual void load(
		ImageType *org_noisy,		// pointer of the input noisy grayscale image
		ImageType *org_basic,		// pointer of the input noisy grayscale image
		int sigma,					// sigma of the Y component
		DistType max_mdist = 2500,	// maximum mean distance (L2/L1) between the reference and the candidate
		int sigmau = -1,			// sigma of the U component, has no use for YUV 4:0:0
		int sigmav = -1				// sigma of the V component, has no use for YUV 4:0:0
		);

	/* reset the buffers and redirect the processing reference patch to the first one of the noisy image */
	virtual void reset();

	/* Denoise just a line of reference patches and write out the completed rows. */
	virtual int next_line(
		ImageType *clean			// pointer of the output denoised grayscale image
		);

	/* Denoise a whole grayscale image and write out the result. */
	void run(
		ImageType *clean			// pointer of the output denoised grayscale image
		);

	/* grouping step of a single patch */
	void grouping();

	/* filtering step of a single patch */
	void filtering();

	/* aggregation step of a single patch */
	void aggregation();

	/* discard the first (pstep) rows and shift the remains of the numerator/denominator buffer */
	void shift_numer_denom();

protected:
	int orig_w;	// original image width
	int orig_h;			// original image height
	int w;				// padded image width
	int h;				// padded image height
	ImageType *noisy;	// padded noisy image
	ImageType *basic;	// padded basic image

	int psize;			// patch size
	int pstep;			// reference patch step

	int swinrh;			// horizontal search window radius
	int ssteph;			// horizontal search step

	int swinrv;			// vertical search window radius
	int sstepv;			// vertical search step

	Group3D *g3d_basic;		// 3d group containg the reference patch and all its similar ones
	Group3D *g3d_noisy;		// 3d group containg the reference patch and all its similar ones
	ImageType *refer_noisy;	// reference patch pointer (top-left) of noisy image
	ImageType *refer_basic;	// reference patch pointer (top-left) of basic image

	int row_cnt;		// counter of the processed rows of the original image

	PatchType *numerator;		// size: w * (2 * swinrv + psize)
	PatchType *denominator;		// size: w * (2 * swinrv + psize)
	PatchType *numer;			// template pointer
	PatchType *denom;			// template pointer
	double wie_wgt_sum;

	DistType *dist_buf;		// sliding buffer to record the distances step by step
	DistType *dist_sum;		// distances buffer of each candidate patch

	int nbuf;				// number of steps in a single patch, ceil(psize / pstep)
	int ncnt;				// counter of the steps

	int nsh;				// number of horizontal candidate patches in a searching window
	int nsv;				// number of vertical candidate patches in a searching window

	clock_t gtime;			// timer of the grouping step
	clock_t ftime;			// timer of the filtering step
	clock_t atime;			// timer of the aggregation step
};

#endif
