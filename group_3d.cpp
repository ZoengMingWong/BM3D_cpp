#include <iostream>
#include "group_3d.h"

const PatchType Group3D::sqrt_powN_x32[8] = {32, 45, 64, 90, 128, 180, 256, 360};

Group3D::Group3D(int w_, int h_, int maxp)
	: w(w_), h(h_), max_patches(maxp)
{
	patch = new Patch2D *[max_patches];
	buf   = new Patch2D *[max_patches];
	for (int i = 0; i < max_patches; i++) 
	{
		patch[i] = new Patch2D(w, h);
	}
}

Group3D::~Group3D()
{
	for (int i = 0; i < max_patches; i++)
	{
		delete patch[i];
	}
	delete[] patch;
	delete[] buf;
}

void Group3D::set_thresholds(int sigma, DistType maxd)
{
	max_dist = maxd;
	thres = (PatchType)(HARD_THRES_MULTIPLIER * sigma) * (1 << COEFF_DICI_BITS);
}

void Group3D::set_reference(ImageType *refer, int stride)
{
	patch[0]->update(refer, 0, 0, 0, stride);
	num = 1;
}

void Group3D::set_reference()
{
	patch[0]->update(0, 0, 0);
	num = 1;
}

int Group3D::find_idx(DistType d)
{
	int i = 0;
	int n = num;
	while (n > 0)
	{
		if (d >= patch[i + (n - 1) / 2]->dist)
		{
			i += (n + 1) / 2;
		}
		n /= 2;
	}
	return i;
}

void Group3D::insert_patch(int x, int y, DistType d)
{
	if (x == 0 && y == 0) return;
	if (d > max_dist) return;

	int idx = find_idx(d);
	if (idx >= max_patches) return;
	if (num >= max_patches) num--;

	Patch2D *tmp = patch[num];
	for (int i = num; i > idx; i--)
	{
		patch[i] = patch[i - 1];
	}
	patch[idx] = tmp;
	patch[idx]->update(x, y, d);
	num++;
}

void Group3D::fill_patches_values(ImageType *refer, int stride)
{
	// truncate the number of patches to power of 2
	log_num = 0;
	while (num > 1) {
		num >>= 1;
		log_num++;
	}
	num = 1 << log_num;

	for (int p = 0; p < num; p++)
	{
		patch[p]->update(refer, stride);
	}
}

void Group3D::transform_3d()
{
	for (int p = 0; p < num; p++) 
	{
		patch[p]->transform_2d();
	}
	hadamard_1d();
}

void Group3D::inv_transform_3d()
{
	if (num == 1) 
	{
		patch[0]->inv_transform_2d();
		return;
	}
	
	hadamard_1d();
	for (int p = 0; p < num; p++) 
	{
		// normalization of the Hadamard transform
		for (int i = 0; i < w*h; i++) 
		{
#if USE_INTEGER
			patch[p]->values[i] = (patch[p]->values[i] + (1 << (log_num - 1))) >> log_num;
#else
			patch[p]->values[i] /= num;
#endif
		}
		patch[p]->inv_transform_2d();
	}
}

void Group3D::hard_thresholding()
{
	nonzeros = 0;
	PatchType tmp_thres = thres * sqrt_powN_x32[log_num] / 32;
	for (int p = 0; p < num; p++)
	{
		for (int i = 0; i < w*h; i++)
		{
			if (patch[p]->values[i] >= tmp_thres || patch[p]->values[i] <= -tmp_thres)
				nonzeros++;
			else
				patch[p]->values[i] = 0;
		}
	}
}

/* Get the weight of the 3D group.
 * All the pixels (regardless of the pixel location) in the 3D group share the same weight.
 * The weight is usually inversely proportional to the number of nonzero coefficients after hard-threshold filtering.
 * You can modify the function between the weight and the nonzero coefficients to obtain a better result.
 */
PatchType Group3D::get_weight()
{
#if USE_INTEGER
	return 1;
#else
	return 1;
	//return (PatchType)(nonzeros ? (1.f / nonzeros) : 1.f);
#endif
}

/* Inplace implementation of 1D Hadamard transform for the 3D group.
 * The length of the 1D transform is (this->num), 
 * and there are totally (w * h) transforms that one for each pixel location independently.
 * For a 1D array with length 4, for example, [a1, a2, a3, a4],
 * we first construct a new array [b1, b2, b3, b4] = [a1+a2, a3+a4, a1-a2, a3-a4],
 * and then we construct a new array [c1, c2, c3, c4] = [b1+b2, b3+b4, b1-b2, b3-b4].
 * Now, we have (c1=a1+a2+a3+a4), and the array [c] is the 1D Hadamard transform of array [a].
 * For other lengths, we just need to repeatedly construct the new array as above, 
 * until the fisrt element of the new array is the sum of all elements of the original array.
 * For a length (n), we need to repeat (log2 n) times to construct the new array.
 */
void Group3D::hadamard_1d()
{
	PatchType tmp;
	for (int n = 0; n < log_num; n++)
	{
		for (int p = 0; p < num / 2; p++)
		{
			buf[p] = patch[2 * p];
			buf[p + num / 2] = patch[2 * p + 1];
		}
		for (int p = 0; p < num; p += 2)
		{
			for (int i = 0; i < w*h; i++)
			{
				tmp = patch[p]->values[i] - patch[p + 1]->values[i];
				patch[p]->values[i] = patch[p]->values[i] + patch[p + 1]->values[i];
				patch[p + 1]->values[i] = tmp;
			}
		}
		for (int p = 0; p < num; p++)
		{
			patch[p] = buf[p];
		}
	}
}
