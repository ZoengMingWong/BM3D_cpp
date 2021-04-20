#include "transform.h"

/* Inplace implementation of the forward 2D 8x8 Bior-1.5 wavelet transform.
 * Firstly, the 8x8 matrix below is applied to each row and then each column (vice versa) of the input 8x8 patch. 
 *					[ 64,  64,  11, -11,   0,   0, -11,  11] [x0]
 *					[-11,  11,  64,  64,  11, -11,   0,   0] [x1]
 *		   1		[  0,   0, -11,  11,  64,  64,  11, -11] [x2]
 *	  ----------- *	[ 11, -11,   0,   0, -11,  11,  64,  64] [x3]
 *	   64*sqrt(2)	[ 64, -64,   0,   0,   0,   0,   0,   0] [x4]
 *					[  0,   0,  64, -64,   0,   0,   0,   0] [x5]
 *					[  0,   0,   0,   0,  64, -64,   0,   0] [x6]
 *					[  0,   0,   0,   0,   0,   0,  64, -64] [x7]
 *
 * Secondly, the 4x4 matrix below is applied to each row and then each column of the top-left 4x4 patch 
 * of the output of the fisrt step, and the remains are under *unchanged* (except scaling if necessary).
 *					[1,  1,  0,  0] [x0]
 *			 		[0,  0,  1,  1] [x1]
 *		1/sqrt(2) * [1, -1,  0,  0] [x2]
 *					[0,  0,  1, -1] [x3]
 *
 * Finally, the 2x2 matrix below is applied to each row and then each column of the top-left 2x2 patch 
 * of the output of the second step, and the remains are under *unchanged* (except scaling if necessary).
 *					[1,  1] [x0]
 *		1/sqrt(2) * [1, -1] [x1]
 */
void inplace_forward_bior15_2d_8x8(float *src)
{
	static float buf[4];
	float *org_src = src;

	// horizontal transform of the 1st step (8x8)
	for (int i = 0; i < 8; i++)
	{
		// row (i)
		buf[0] = src[0] - src[1];
		buf[1] = src[2] - src[3];
		buf[2] = src[4] - src[5];
		buf[3] = src[6] - src[7];

		src[0] = 64 * (src[0] + src[1]) + 11 * (buf[1] - buf[3]);
		src[1] = 64 * (src[2] + src[3]) + 11 * (buf[2] - buf[0]);
		src[2] = 64 * (src[4] + src[5]) + 11 * (buf[3] - buf[1]);
		src[3] = 64 * (src[6] + src[7]) + 11 * (buf[0] - buf[2]);

		src[4] = buf[0] * 64;
		src[5] = buf[1] * 64;
		src[6] = buf[2] * 64;
		src[7] = buf[3] * 64;
		src += 8;
	}

	// vertical transform of the 1st step (8x8)
	src = org_src;
	for (int j = 0; j < 8; j++)
	{
		// row (i)
		buf[0] = src[0 * 8 + j] - src[1 * 8 + j];
		buf[1] = src[2 * 8 + j] - src[3 * 8 + j];
		buf[2] = src[4 * 8 + j] - src[5 * 8 + j];
		buf[3] = src[6 * 8 + j] - src[7 * 8 + j];

		src[0 * 8 + j] = 64 * (src[0 * 8 + j] + src[1 * 8 + j]) + 11 * (buf[1] - buf[3]);
		src[1 * 8 + j] = 64 * (src[2 * 8 + j] + src[3 * 8 + j]) + 11 * (buf[2] - buf[0]);
		src[2 * 8 + j] = 64 * (src[4 * 8 + j] + src[5 * 8 + j]) + 11 * (buf[3] - buf[1]);
		src[3 * 8 + j] = 64 * (src[6 * 8 + j] + src[7 * 8 + j]) + 11 * (buf[0] - buf[2]);

		src[4 * 8 + j] = buf[0] * 64;
		src[5 * 8 + j] = buf[1] * 64;
		src[6 * 8 + j] = buf[2] * 64;
		src[7 * 8 + j] = buf[3] * 64;
	}

	for (int i = 0; i < 8; i++) 
	{
		for (int j = 0; j < 8; j++) {
			src[8 * i + j] /= 8192.f;	// (64 * 2^0.5)^2, normalization
		}
	}

	// horizontal transform of the 2nd step (4x4)
	for (int i = 0; i < 4; i++)
	{
		buf[0] = src[0] - src[1];
		buf[1] = src[2] - src[3];
		src[0] = src[0] + src[1];
		src[1] = src[2] + src[3];
		src[2] = buf[0];
		src[3] = buf[1];
		src += 8;
	}

	// vertical transform of the 2nd step (4x4)
	src = org_src;
	for (int j = 0; j < 4; j++)
	{
		buf[0] = src[0 * 8 + j] - src[1 * 8 + j];
		buf[1] = src[2 * 8 + j] - src[3 * 8 + j];
		src[0 * 8 + j] = src[0 * 8 + j] + src[1 * 8 + j];
		src[1 * 8 + j] = src[2 * 8 + j] + src[3 * 8 + j];
		src[2 * 8 + j] = buf[0];
		src[3 * 8 + j] = buf[1];
	}

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++) {
			src[8 * i + j] /= 2.f;	// (2^0.5)^2, normalization
		}
	}

	// horizontal transform of the 3rd step (2x2)
	for (int i = 0; i < 2; i++)
	{
		buf[0] = src[0] - src[1];
		src[0] = src[0] + src[1];
		src[1] = buf[0];
		src += 8;
	}

	// vertical transform of the 3rd step (2x2)
	src = org_src;
	for (int j = 0; j < 2; j++)
	{
		buf[0] = src[0 * 8 + j] - src[1 * 8 + j];
		src[0 * 8 + j] = src[0 * 8 + j] + src[1 * 8 + j];
		src[1 * 8 + j] = buf[0];
	}

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++) {
			src[8 * i + j] /= 2.f;	// (2^0.5)^2, normalization
		}
	}
}

/* Inplace implementation of the backward 2D 8x8 Bior-1.5 wavelet transform.
 * Firstly, the 2x2 matrix below is applied to each column and then each row (vice versa)
 * of the top-left 2x2 patch of the input, and the remains are under *unchanged* (except scaling if necessary).
 * Note that the output 2x2 patch will replace the original top-left 2x2 patch of the input.
 *					[1,  1] [x0]
 *		1/sqrt(2) * [1, -1] [x1]
 *
 * Secondly, the 4x4 matrix below is applied to each column and then each row of the top-left 4x4 patch 
 * of the input after the fisrt step, and the remains are under *unchanged* (except scaling if necessary).
 *					[1,  0,  1,  0] [x0]
 *			 		[1,  0, -1,  0] [x1]
 *		1/sqrt(2) * [0,  1,  0,  1] [x2]
 *					[0,  1,  0, -1] [x3]
 *
 * Finally, the 8x8 matrix below is applied to each column and then each row 
 * of the whole input patch after the second step.
 *					[ 64,   0,   0,   0,  64, -11,   0,  11] [x0]
 *					[ 64,   0,   0,   0, -64, -11,   0,  11] [x1]
 *         1		[  0,  64,   0,   0,  11,  64, -11,   0] [x2]
 *    ----------- *	[  0,  64,   0,   0,  11, -64, -11,   0] [x3]
 *	   64*sqrt(2)	[  0.,  0,  64,   0,   0,  11,  64, -11] [x4]
 *					[  0.,  0,  64,   0,   0,  11, -64, -11] [x5]
 *					[  0,   0,   0,  64, -11,   0,  11,  64] [x6]
 *					[  0,   0,   0,  64, -11,   0,  11, -64] [x7]
 */
void inplace_backward_bior15_2d_8x8(float *src)
{
	static float buf[4];
	float *org_src = src;

	// vertical transform of the 1st step (2x2)
	for (int j = 0; j < 2; j++)
	{
		buf[0] = src[0 * 8 + j] - src[1 * 8 + j];
		src[0 * 8 + j] = src[0 * 8 + j] + src[1 * 8 + j];
		src[1 * 8 + j] = buf[0];
	}

	// horizontal transform of the 1st step (2x2)
	for (int i = 0; i < 2; i++)
	{
		buf[0] = src[0] - src[1];
		src[0] = src[0] + src[1];
		src[1] = buf[0];
		src += 8;
	}
	src = org_src;

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++) {
			src[8 * i + j] /= 2.f;	// (2^0.5)^2, normalization
		}
	}

	// vertical transform of the 2nd step (4x4)
	for (int j = 0; j < 4; j++)
	{
		buf[0] = src[0 * 8 + j] - src[2 * 8 + j];
		buf[1] = src[1 * 8 + j] - src[3 * 8 + j];
		src[0 * 8 + j] = src[0 * 8 + j] + src[2 * 8 + j];
		src[2 * 8 + j] = src[1 * 8 + j] + src[3 * 8 + j];
		src[1 * 8 + j] = buf[0];
		src[3 * 8 + j] = buf[1];
	}

	// horizontal transform of the 2nd step (4x4)
	for (int i = 0; i < 4; i++)
	{
		buf[0] = src[0] - src[2];
		buf[1] = src[1] - src[3];
		src[0] = src[0] + src[2];
		src[2] = src[1] + src[3];
		src[1] = buf[0];
		src[3] = buf[1];
		src += 8;
	}
	src = org_src;

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			src[8 * i + j] /= 2.f;	// (2^0.5)^2, normalization

	// vertical transform of the 3rd step (8x8)
	for (int j = 0; j < 8; j++)
	{
		buf[0] = (src[0 * 8 + j] - src[4 * 8 + j]) * 64;
		buf[1] = (src[1 * 8 + j] - src[5 * 8 + j]) * 64;
		buf[2] = (src[2 * 8 + j] - src[6 * 8 + j]) * 64;
		buf[3] = (src[3 * 8 + j] - src[7 * 8 + j]) * 64;

		src[0 * 8 + j] = (src[0 * 8 + j] + src[4 * 8 + j]) * 64;
		src[1 * 8 + j] = (src[1 * 8 + j] + src[5 * 8 + j]) * 64;
		src[2 * 8 + j] = (src[2 * 8 + j] + src[6 * 8 + j]) * 64;
		src[3 * 8 + j] = (src[3 * 8 + j] + src[7 * 8 + j]) * 64;

		src[5 * 8 + j] = (src[5 * 8 + j] - src[7 * 8 + j]) * 11;
		src[6 * 8 + j] = (src[6 * 8 + j] - src[4 * 8 + j]) * 11;

		src[4 * 8 + j] = src[5 * 8 + j];
		src[5 * 8 + j] = buf[2] + src[5 * 8 + j];
		src[7 * 8 + j] = buf[3] + src[6 * 8 + j];

		src[0 * 8 + j] = src[0 * 8 + j] - src[4 * 8 + j];

		buf[2] = src[1 * 8 + j];
		buf[3] = src[2 * 8 + j];
		src[1 * 8 + j] = buf[0] - src[4 * 8 + j];
		src[2 * 8 + j] = buf[2] - src[6 * 8 + j];

		buf[0] = src[3 * 8 + j];
		src[3 * 8 + j] = buf[1] - src[6 * 8 + j];
		src[4 * 8 + j] = buf[3] + src[4 * 8 + j];
		src[6 * 8 + j] = buf[0] + src[6 * 8 + j];
	}

	// horizontal transform of the 3rd step (8x8)
	for (int i = 0; i < 8; i++)
	{
		buf[0] = (src[0] - src[4]) * 64;
		src[0] = (src[0] + src[4]) * 64;
		buf[1] = (src[1] - src[5]) * 64;
		src[1] = (src[1] + src[5]) * 64;
		buf[2] = (src[2] - src[6]) * 64;
		src[2] = (src[2] + src[6]) * 64;
		buf[3] = (src[3] - src[7]) * 64;
		src[3] = (src[3] + src[7]) * 64;

		src[5] = (src[5] - src[7]) * 11;
		src[6] = (src[6] - src[4]) * 11;
		src[4] = src[5];
		src[5] = buf[2] + src[5];
		src[7] = buf[3] + src[6];
		src[0] = src[0] - src[4];

		buf[2] = src[1];
		buf[3] = src[2];
		src[1] = buf[0] - src[4];
		src[2] = buf[2] - src[6];
		buf[0] = src[3];
		src[3] = buf[1] - src[6];
		src[4] = buf[3] + src[4];
		src[6] = buf[0] + src[6];
		src += 8;
	}
	src = org_src;

	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++) {
			src[8 * i + j] /= 8192.f;	// (64 * 2^0.5)^2, normalization
		}
	}
}


/* Integer version of the 8x8 Bior-1.5 forward transform.
 * It's almost the same as the floating-point version but only the scaling.
 * The output coefficients are multiplied by 2 in comparision with the floating-point ones.
 * There are 5 more bits of the ouput coefficients than the input values, and totally at most 
 * 18 bits are needed for the intermediate results compared to the input values, e.g. 26 bits for uint8.
 */
void inplace_forward_bior15_2d_8x8(int *src)
{
	static int buf[4];
	int *org_src = src;

	// horizontal transform of the 1st step (8x8)
	for (int i = 0; i < 8; i++)
	{
		// row (i)
		buf[0] = src[0] - src[1];
		buf[1] = src[2] - src[3];
		buf[2] = src[4] - src[5];
		buf[3] = src[6] - src[7];

		src[0] = 64 * (src[0] + src[1]) + 11 * (buf[1] - buf[3]);
		src[1] = 64 * (src[2] + src[3]) + 11 * (buf[2] - buf[0]);
		src[2] = 64 * (src[4] + src[5]) + 11 * (buf[3] - buf[1]);
		src[3] = 64 * (src[6] + src[7]) + 11 * (buf[0] - buf[2]);

		src[4] = buf[0] * 64;
		src[5] = buf[1] * 64;
		src[6] = buf[2] * 64;
		src[7] = buf[3] * 64;
		src += 8;
	}

	// vertical transform of the 1st step (8x8)
	src = org_src;
	for (int j = 0; j < 8; j++)
	{
		// row (i)
		buf[0] = src[0 * 8 + j] - src[1 * 8 + j];
		buf[1] = src[2 * 8 + j] - src[3 * 8 + j];
		buf[2] = src[4 * 8 + j] - src[5 * 8 + j];
		buf[3] = src[6 * 8 + j] - src[7 * 8 + j];

		src[0 * 8 + j] = 64 * (src[0 * 8 + j] + src[1 * 8 + j]) + 11 * (buf[1] - buf[3]);
		src[1 * 8 + j] = 64 * (src[2 * 8 + j] + src[3 * 8 + j]) + 11 * (buf[2] - buf[0]);
		src[2 * 8 + j] = 64 * (src[4 * 8 + j] + src[5 * 8 + j]) + 11 * (buf[3] - buf[1]);
		src[3 * 8 + j] = 64 * (src[6 * 8 + j] + src[7 * 8 + j]) + 11 * (buf[0] - buf[2]);

		src[4 * 8 + j] = buf[0] * 64;
		src[5 * 8 + j] = buf[1] * 64;
		src[6 * 8 + j] = buf[2] * 64;
		src[7 * 8 + j] = buf[3] * 64;
	}

	// horizontal transform of the 2nd step (4x4)
	for (int i = 0; i < 4; i++)
	{
		buf[0] = src[0] - src[1];
		buf[1] = src[2] - src[3];
		src[0] = src[0] + src[1];
		src[1] = src[2] + src[3];
		src[2] = buf[0];
		src[3] = buf[1];
		src += 8;
	}

	// vertical transform of the 2nd step (4x4)
	src = org_src;
	for (int j = 0; j < 4; j++)
	{
		buf[0] = src[0 * 8 + j] - src[1 * 8 + j];
		buf[1] = src[2 * 8 + j] - src[3 * 8 + j];
		src[0 * 8 + j] = (src[0 * 8 + j] + src[1 * 8 + j] + 1) >> 1;
		src[1 * 8 + j] = (src[2 * 8 + j] + src[3 * 8 + j] + 1) >> 1;
		src[2 * 8 + j] = (buf[0] + 1) >> 1;
		src[3 * 8 + j] = (buf[1] + 1) >> 1;
	}

	// horizontal transform of the 3rd step (2x2)
	for (int i = 0; i < 2; i++)
	{
		buf[0] = src[0] - src[1];
		src[0] = src[0] + src[1];
		src[1] = buf[0];
		src += 8;
	}

	// vertical transform of the 3rd step (2x2)
	src = org_src;
	for (int j = 0; j < 2; j++)
	{
		buf[0] = src[0 * 8 + j] - src[1 * 8 + j];
		src[0 * 8 + j] = (src[0 * 8 + j] + src[1 * 8 + j] + 1) >> 1;
		src[1 * 8 + j] = (buf[0] + 1) >> 1;
	}

	// normalization
	for (int i = 0; i < 64; i++) 
	{
		src[i] = (src[i] + (1 << 11)) >> 12;
	}
}

/* Integer version of the 8x8 Bior-1.5 backward transform.
 * It's almost the same as the floating-point version but only the scaling.
 * The output has the same bits as the original image values.
 * The bits needed for the intermediate results remain the same as the forward transform.
 */
void inplace_backward_bior15_2d_8x8(int *src)
{
	static int buf[4];
	int *org_src = src;

	// vertical transform of the 1st step (2x2)
	for (int j = 0; j < 2; j++)
	{
		buf[0] = src[0 * 8 + j] - src[1 * 8 + j];
		src[0 * 8 + j] = src[0 * 8 + j] + src[1 * 8 + j];
		src[1 * 8 + j] = buf[0];
	}

	// horizontal transform of the 1st step (2x2)
	for (int i = 0; i < 2; i++)
	{
		buf[0] = src[0] - src[1];
		src[0] = (src[0] + src[1] + 1) >> 1;
		src[1] = (buf[0] + 1) >> 1;
		src += 8;
	}
	src = org_src;

	// vertical transform of the 2nd step (4x4)
	for (int j = 0; j < 4; j++)
	{
		buf[0] = src[0 * 8 + j] - src[2 * 8 + j];
		buf[1] = src[1 * 8 + j] - src[3 * 8 + j];
		src[0 * 8 + j] = src[0 * 8 + j] + src[2 * 8 + j];
		src[2 * 8 + j] = src[1 * 8 + j] + src[3 * 8 + j];
		src[1 * 8 + j] = buf[0];
		src[3 * 8 + j] = buf[1];
	}

	// horizontal transform of the 2nd step (4x4)
	for (int i = 0; i < 4; i++)
	{
		buf[0] = src[0] - src[2];
		buf[1] = src[1] - src[3];
		src[0] = (src[0] + src[2] + 1) >> 1;
		src[2] = (src[1] + src[3] + 1) >> 1;
		src[1] = (buf[0] + 1) >> 1;
		src[3] = (buf[1] + 1) >> 1;
		src += 8;
	}
	src = org_src;

	// vertical transform of the 3rd step (8x8)
	for (int j = 0; j < 8; j++)
	{
		buf[0] = (src[0 * 8 + j] - src[4 * 8 + j]) * 64;
		buf[1] = (src[1 * 8 + j] - src[5 * 8 + j]) * 64;
		buf[2] = (src[2 * 8 + j] - src[6 * 8 + j]) * 64;
		buf[3] = (src[3 * 8 + j] - src[7 * 8 + j]) * 64;

		src[0 * 8 + j] = (src[0 * 8 + j] + src[4 * 8 + j]) * 64;
		src[1 * 8 + j] = (src[1 * 8 + j] + src[5 * 8 + j]) * 64;
		src[2 * 8 + j] = (src[2 * 8 + j] + src[6 * 8 + j]) * 64;
		src[3 * 8 + j] = (src[3 * 8 + j] + src[7 * 8 + j]) * 64;

		src[5 * 8 + j] = (src[5 * 8 + j] - src[7 * 8 + j]) * 11;
		src[6 * 8 + j] = (src[6 * 8 + j] - src[4 * 8 + j]) * 11;

		src[4 * 8 + j] = src[5 * 8 + j];
		src[5 * 8 + j] = buf[2] + src[5 * 8 + j];
		src[7 * 8 + j] = buf[3] + src[6 * 8 + j];

		src[0 * 8 + j] = src[0 * 8 + j] - src[4 * 8 + j];

		buf[2] = src[1 * 8 + j];
		buf[3] = src[2 * 8 + j];
		src[1 * 8 + j] = buf[0] - src[4 * 8 + j];
		src[2 * 8 + j] = buf[2] - src[6 * 8 + j];

		buf[0] = src[3 * 8 + j];
		src[3 * 8 + j] = buf[1] - src[6 * 8 + j];
		src[4 * 8 + j] = buf[3] + src[4 * 8 + j];
		src[6 * 8 + j] = buf[0] + src[6 * 8 + j];
	}

	// horizontal transform of the 3rd step (8x8)
	for (int i = 0; i < 8; i++)
	{
		buf[0] = (src[0] - src[4]) * 64;
		src[0] = (src[0] + src[4]) * 64;
		buf[1] = (src[1] - src[5]) * 64;
		src[1] = (src[1] + src[5]) * 64;
		buf[2] = (src[2] - src[6]) * 64;
		src[2] = (src[2] + src[6]) * 64;
		buf[3] = (src[3] - src[7]) * 64;
		src[3] = (src[3] + src[7]) * 64;

		src[5] = (src[5] - src[7]) * 11;
		src[6] = (src[6] - src[4]) * 11;
		src[4] = src[5];
		src[5] = buf[2] + src[5];
		src[7] = buf[3] + src[6];
		src[0] = src[0] - src[4];

		buf[2] = src[1];
		buf[3] = src[2];
		src[1] = buf[0] - src[4];
		src[2] = buf[2] - src[6];
		buf[0] = src[3];
		src[3] = buf[1] - src[6];
		src[4] = buf[3] + src[4];
		src[6] = buf[0] + src[6];
		src += 8;
	}
	src = org_src;

	// normalization
	for (int i = 0; i < 64; i++)
	{
		src[i] = (src[i] + (1 << 13)) >> 14;
	}
}
