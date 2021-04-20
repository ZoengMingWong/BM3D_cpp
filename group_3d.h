#ifndef __PATCH_3D_H__
#define __PATCH_3D_H__

#include <iostream>
#include "patch_2d.h"

struct Group3D
{
	int w;				// patch width
	int h;				// patch height
	int max_patches;	// maximum similar patches

	int num;			// number of patches (will be truncated to power of 2)
	int log_num;		// log2(num)

	DistType max_dist;	// maximum sum of distances (L2/L1) between two patches

	PatchType thres;	// hard threshold of the filtering
	int nonzeros;		// number of nonzero coefficients

	Patch2D **patch;	// array of pointers of 2D patches
	Patch2D **buf;		// array of pointers used as buffer (in Hadamard transform)

	static const PatchType sqrt_powN_x32[8];	// integer of (sqrt(1<<n) * 32)

	Group3D(int w_, int h_, int maxp);
	~Group3D();

	void set_thresholds(int sigma, DistType maxd);

	void set_reference(ImageType *refer, int stride);
	void set_reference();

	// find the index to insert with given distance
	int find_idx(DistType d);

	void insert_patch(int x, int y, DistType d);
	void fill_patches_values(ImageType *refer, int stride);

	// forward and backward are the same except the scaling
	void hadamard_1d();

	void transform_3d();
	void inv_transform_3d();

	void hard_thresholding();

	PatchType get_weight();
};

#endif