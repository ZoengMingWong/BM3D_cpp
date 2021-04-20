#ifndef __PATCH_2D_H__
#define __PATCH_2D_H__

#include <iostream>
#include "global_define.h"
#include "transform.h"

struct Patch2D
{
	int w;				// patch width
	int h;				// patch height
	int x;				// patch horizontal offset (top-left corner) to the reference patch
	int y;				// patch vertical offset
	PatchType *values;	// patch pixels' values
	DistType dist;		// L2/L1 distance between the patch and its reference one

	Patch2D(int w_, int h_);
	Patch2D(ImageType *image, int x_, int y_, DistType d, int w_, int h_, int stride);
	~Patch2D();

	void update(ImageType *image, int x_, int y_, DistType d, int stride);
	void update(ImageType *image, int stride);
	void update(int x_, int y_, DistType d);

	void transform_2d();
	void inv_transform_2d();
};










#endif

