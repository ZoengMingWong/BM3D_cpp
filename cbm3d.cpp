#include <iostream>
#include "cbm3d.h"

CBM3D::CBM3D(
	int w_,					// width
	int h_,					// height
	int max_sim,			// maximum similar patches
	int psize_,				// reference patch size
	int pstep_,				// reference patch step
	int swinrh_,			// horizontal search window radius
	int ssteph_,			// horizontal search step
	int swinrv_,			// vertical search window radius
	int sstepv_				// vertical search step
) : BM3D(w_, h_, max_sim, psize_, pstep_, swinrh_, ssteph_, swinrv_, sstepv_)
{
	noisy_yuv[0]       = noisy;
	numerator_yuv[0]   = numerator;
	denominator_yuv[0] = denominator;

	for (int i = 1; i < 3; i++)
	{
		noisy_yuv[i]       = new ImageType[w * h]();
		numerator_yuv[i]   = new PatchType[w * (psize + swinrv * 2)];
		denominator_yuv[i] = new PatchType[w * (psize + swinrv * 2)];
	}
}

CBM3D::~CBM3D()
{
	for (int i = 1; i < 3; i++)
	{
		delete[] noisy_yuv[i];
		delete[] numerator_yuv[i];
		delete[] denominator_yuv[i];
	}

	// the ~BM3D is called after the ~CBM3D
	noisy       = noisy_yuv[0];
	numerator   = numerator_yuv[0];
	denominator = denominator_yuv[0];
}

void CBM3D::reset()
{
	row_cnt = 0;
	for (int i = 0; i < 3; i++)
	{
		memset(numerator_yuv[i],   0, (psize + swinrv * 2) * w * sizeof(PatchType));
		memset(denominator_yuv[i], 0, (psize + swinrv * 2) * w * sizeof(PatchType));
	}
}

void CBM3D::load(ImageType *org_noisy_yuv, int sigmay, DistType max_mdist, int sigmau, int sigmav)
{
	for (int i = 0; i < 3; i++)
	{
		noisy       = noisy_yuv[i];
		numerator   = numerator_yuv[i];
		denominator = denominator_yuv[i];

		BM3D::load(org_noisy_yuv, sigmay, max_mdist);
		org_noisy_yuv += (orig_w * orig_h);
	}

	hard_thres[0] = (PatchType)(HARD_THRES_MULTIPLIER * sigmay) * (1 << COEFF_DICI_BITS);
	hard_thres[1] = sigmau < 0 ? hard_thres[0] : (PatchType)(HARD_THRES_MULTIPLIER * sigmau) * (1 << COEFF_DICI_BITS);
	hard_thres[2] = sigmav < 0 ? hard_thres[0] : (PatchType)(HARD_THRES_MULTIPLIER * sigmav) * (1 << COEFF_DICI_BITS);
}

int CBM3D::next_line(ImageType *clean)
{
	if (row_cnt >= orig_h + pstep - psize) return -1;	// beyond the last line of reference patches

	for (int i = 0; i < 3; i++)
	{
		refer_yuv[i] = noisy_yuv[i]       + swinrh + (row_cnt + swinrv) * w;
		numer_yuv[i] = numerator_yuv[i]   + swinrh + swinrv * w;
		denom_yuv[i] = denominator_yuv[i] + swinrh + swinrv * w;
	}

	memset(dist_buf, 0, nsh * nsv * nbuf * sizeof(DistType));
	memset(dist_sum, 0, nsh * nsv * sizeof(DistType));

	// initialize the distance buffer
#pragma omp parallel for num_threads(USE_THREADS_NUM)
	for (int sy = -swinrv; sy <= swinrv; sy += sstepv)
	{
		for (int sx = -swinrh; sx <= swinrh; sx += ssteph)
		{
			int idx = (sy + swinrv) / sstepv * nsh + (sx + swinrh) / ssteph;
			for (int y = 0; y < psize; y++) {
				for (int x = 0; x < psize - pstep; x++) {
					dist_buf[nbuf * idx + x / pstep] += get_dist(refer_yuv[0][y * w + x], refer_yuv[0][(y + sy) * w + x + sx]);
				}
			}
			for (int i = 0; i < nbuf - 2; i++) {
				dist_sum[idx] += dist_buf[nbuf * idx + i];
			}
		}
	}
	ncnt = nbuf;

	clock_t t;
	// proceesing the line
	for (int x = 0; x < orig_w + pstep - psize; x += pstep, ncnt++)
	{
		refer = refer_yuv[0];
		numer = numer_yuv[0];
		denom = denom_yuv[0];

		g3d->thres = hard_thres[0];

		t = clock();
		grouping();
		gtime += clock() - t;

		t = clock();
		filtering();
		ftime += clock() - t;

		t = clock();
		aggregation();
		atime += clock() - t;

		for (int i = 1; i < 3; i++)
		{
			refer = refer_yuv[i];
			numer = numer_yuv[i];
			denom = denom_yuv[i];

			g3d->thres = hard_thres[i];

			t = clock();
			g3d->fill_patches_values(refer, w);
			gtime += clock() - t;

			t = clock();
			filtering();
			ftime += clock() - t;

			t = clock();
			aggregation();
			atime += clock() - t;
		}

		for (int i = 0; i < 3; i++)
		{
			refer_yuv[i] += pstep;
			numer_yuv[i] += pstep;
			denom_yuv[i] += pstep;
		}
	}

	// output the completed rows
	for (int i = 0; i < 3; i++)
	{
		numer_yuv[i] = numerator_yuv[i] + swinrh;
		denom_yuv[i] = denominator_yuv[i] + swinrh;
	}

	int output_rows;
	if (row_cnt < swinrv)
	{
		if (row_cnt + pstep <= swinrv) {
			// no row is completed in the begining
			output_rows = 0;
		}
		else {
			output_rows = row_cnt + pstep - swinrv;
			for (int i = 0; i < 3; i++) {
				numer_yuv[i] += (pstep - output_rows) * w;
				denom_yuv[i] += (pstep - output_rows) * w;
			}
		}
	}
	else
	{
		if (row_cnt >= orig_h - psize)
			output_rows = orig_h - row_cnt + swinrv;
		else
			output_rows = pstep;
		clean += orig_w * (row_cnt - swinrv);
	}

	for (int i = 0; i < 3; i++)
	{
		numer = numer_yuv[i];
		denom = denom_yuv[i];
		for (int idx = 0, r = 0; r < output_rows; r++)
		{
			for (int c = 0; c < orig_w; c++, idx++)
			{
				clean[idx] = (ImageType)(numer[c] / denom[c]);
			}
			numer += w;
			denom += w;
		}
		clean += (orig_w * orig_h);
	}

	// remove the fisrt pstep rows of the numerator and denominator buffers
	// and insert pstep new rows to the end of the buffers
	for (int i = 0; i < 3; i++)
	{
		numerator = numerator_yuv[i];
		denominator = denominator_yuv[i];
		shift_numer_denom();
	}

	row_cnt += pstep;
	return output_rows;
}