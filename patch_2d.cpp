#include <iostream>
#include "patch_2d.h"

Patch2D::Patch2D(int w_, int h_)
	:w(w_), h(h_)
{
	values = new PatchType[w * h];
}

Patch2D::Patch2D(ImageType *image, int x_, int y_, DistType d, int w_, int h_, int stride)
	: w(w_), h(h_), x(x_), y(y_), dist(d)
{
	values = new PatchType[w * h];
	for (int i = 0, r = 0; r < h; r++)
	{
		for (int c = 0; c < w; c++, i++)
		{
			values[i] = (PatchType)image[(r + y) * stride + c + x];
		}
	}
}

Patch2D::~Patch2D()
{
	delete[] values;
}

void Patch2D::update(ImageType *image, int x_, int y_, DistType d, int stride)
{
	x = x_, y = y_, dist = d;
	for (int i = 0, r = 0; r < h; r++) 
	{
		for (int c = 0; c < w; c++, i++) 
		{
			values[i] = (PatchType)image[(r + y) * stride + c + x];
		}
	}
}

void Patch2D::update(ImageType *image, int stride)
{
	for (int i = 0, r = 0; r < h; r++)
	{
		for (int c = 0; c < w; c++, i++)
		{
			values[i] = (PatchType)image[(r + y) * stride + c + x];
		}
	}
}

void Patch2D::update(int x_, int y_, DistType d)
{
	x = x_, y = y_;
	dist = d;
}

void Patch2D::transform_2d()
{
	inplace_forward_bior15_2d_8x8(values);
}

void Patch2D::inv_transform_2d()
{
	inplace_backward_bior15_2d_8x8(values);
}


