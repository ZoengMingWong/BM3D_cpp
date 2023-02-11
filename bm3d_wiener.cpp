#include <iostream>
#include "bm3d_wiener.h"

static inline DistType get_dist(ImageType a, ImageType b)
{
	int diff = (int)a - b;
#if USE_L2_DIST
	return diff * diff;
#else
	return diff >= 0 ? diff : -diff;
#endif
}

BM3D_WIE::BM3D_WIE(
	int w_,					// width
	int h_,					// height
	int max_sim,			// maximum similar patches
	int psize_,				// reference patch size
	int pstep_,				// reference patch step
	int swinrh_,			// horizontal search window radius
	int ssteph_,			// horizontal search step
	int swinrv_,			// vertical search window radius
	int sstepv_			// vertical search step
) : orig_w(w_), orig_h(h_), psize(psize_), pstep(pstep_),
	swinrh(swinrh_), ssteph(ssteph_), swinrv(swinrv_), sstepv(sstepv_)
{
	// pad the last patch to ensure a completed step (copy the last row/column)
	int w_pad = (orig_w - psize + pstep - 1) / pstep * pstep + psize - orig_w;
	int h_pad = (orig_h - psize + pstep - 1) / pstep * pstep + psize - orig_h;

	// pad the surrounding with half search window (zero padding)
	w = orig_w + w_pad + swinrh * 2;
	h = orig_h + h_pad + swinrv * 2;

	g3d_noisy = new Group3D(psize, psize, max_sim);
	g3d_basic = new Group3D(psize, psize, max_sim);
	noisy = new ImageType[w * h]();
	basic = new ImageType[w * h]();

	numerator   = new PatchType[w * (psize + swinrv * 2)];
	denominator = new PatchType[w * (psize + swinrv * 2)];

	// the distances computed by the last patch can be partially reused when stepping forward
	nbuf = (psize + pstep - 1) / pstep;
	nsh  = (2 * swinrh + ssteph) / ssteph;
	nsv  = (2 * swinrv + sstepv) / sstepv;

	dist_buf = new DistType[nsh * nsv * nbuf];
	dist_sum = new DistType[nsh * nsv];

	row_cnt = h;	// avoid processing without the noisy image initialization
}

BM3D_WIE::~BM3D_WIE()
{
	delete g3d_noisy;
	delete g3d_basic;
	delete[] noisy;
	delete[] basic;
	delete[] numerator;
	delete[] denominator;
	delete[] dist_buf;
	delete[] dist_sum;
}

void BM3D_WIE::run(ImageType *clean)
{
	gtime = 0;
	ftime = 0;
	atime = 0;

	if (row_cnt > 0) 
		reset();
	while (next_line(clean) >= 0);

	clock_t stime = gtime + ftime + atime;
	std::cout << "time(s): " 
			  << (double)gtime / CLOCKS_PER_SEC << ' ' 
			  << (double)ftime / CLOCKS_PER_SEC << ' ' 
			  << (double)atime / CLOCKS_PER_SEC << ' ' 
			  << (double)stime / CLOCKS_PER_SEC << std::endl;

	std::cout << "percentage: " 
			  << (double)gtime / stime * 100 << ' ' 
			  << (double)ftime / stime * 100 << ' ' 
			  << (double)atime / stime * 100 << std::endl;
}

void BM3D_WIE::reset()
{
	row_cnt = 0;
	memset(numerator,   0, (psize + swinrv * 2) * w * sizeof(PatchType));
	memset(denominator, 0, (psize + swinrv * 2) * w * sizeof(PatchType));
}

void BM3D_WIE::load(ImageType *org_noisy, ImageType *org_basic, int sigma, DistType max_mdist, int sigmau, int sigmav)
{
	row_cnt = 0;
	g3d_basic->max_dist = max_mdist * psize * psize;
	g3d_basic->thres = sigma * sigma * (1 << (COEFF_DICI_BITS * 2));

	int w_pad = w - 2 * swinrh - orig_w;
	int h_pad = h - 2 * swinrv - orig_h;

	ImageType *tmp_noisy = noisy + swinrv * w + swinrh;
	ImageType *tmp_basic = basic + swinrv * w + swinrh;
	for (int i = 0; i < orig_h; i++)
	{
		memcpy(tmp_noisy, org_noisy, orig_w * sizeof(ImageType));
		memcpy(tmp_basic, org_basic, orig_w * sizeof(ImageType));
		for (int j = 0; j < w_pad; j++)
		{
			// edge padding
			tmp_noisy[orig_w + j] = tmp_noisy[orig_w + j - 1];
			tmp_basic[orig_w + j] = tmp_basic[orig_w + j - 1];
		}
		tmp_noisy += w;
		tmp_basic += w;
		org_noisy += orig_w;
		org_basic += orig_w;
	}
	for (int i = 0; i < h_pad; i++)
	{
		memcpy(tmp_noisy, tmp_noisy - w, (orig_w + w_pad) * sizeof(ImageType));
		memcpy(tmp_basic, tmp_basic - w, (orig_w + w_pad) * sizeof(ImageType));
		tmp_noisy += w;
		tmp_basic += w;
	}
	memset(numerator,   0, (psize + swinrv * 2) * w * sizeof(PatchType));
	memset(denominator, 0, (psize + swinrv * 2) * w * sizeof(PatchType));
}

/* Porcess a line of reference patches, the location of the line is recorded by (this->row_cnt).
 * The distances buffers used in the last reference patches line will be reset in the beginning, 
 * meaning that the grouping process of each reference patches line is independent.
 * When all reference patches in the line are processed, as we will step downward to next line, 
 * the beginning (this->pstep) rows of the buffers (numerator and denominator) will no longer be modified,
 * and we can write out the result to the clean (denoised) image.
 * Note that the number of rows we can write out is different in the beginning or the end, 
 * the function will return how exactly many rows we can write out.
 * The function is unidirectional that, you must call the function one by one 
 * so that the (this->row_cnt) increases step by step from 0 to the end of the image, 
 * as the numerator/denominator buffer records only partial information and updates progressively.
 */
int BM3D_WIE::next_line(ImageType *clean)
{
	if (row_cnt >= orig_h + pstep - psize) return -1;	// beyond the last line of reference patches

	refer_noisy = noisy + (row_cnt + swinrv) * w + swinrh;	// the first reference patch of the line
	refer_basic = basic + (row_cnt + swinrv) * w + swinrh;	// the first reference patch of the line

	numer = numerator   + swinrv * w + swinrh;
	denom = denominator + swinrv * w + swinrh;

	memset(dist_buf, 0, nsh * nsv * nbuf * sizeof(DistType));
	memset(dist_sum, 0, nsh * nsv * sizeof(DistType));

	// initialize the distance buffer
#pragma omp parallel for num_threads(USE_THREADS_NUM)
	for (int sy = -swinrv; sy <= swinrv; sy += sstepv)
	{
		for (int sx = -swinrh; sx <= swinrh; sx += ssteph)
		{
			int idx = (sy + swinrv) / sstepv * nsh + (sx + swinrh) / ssteph;
			for (int y = 0; y < psize; y++)
			{
				for (int x = 0; x < psize - pstep; x++)
				{
					dist_buf[nbuf * idx + x / pstep] += get_dist(refer_basic[y * w + x], refer_basic[(y + sy) * w + x + sx]);
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
		t = clock();
		grouping();
		gtime += clock() - t;

		t = clock();
		filtering();
		ftime += clock() - t;

		t = clock();
		aggregation();
		atime += clock() - t;

		refer_noisy += pstep;
		refer_basic += pstep;
		numer += pstep;
		denom += pstep;
	}

	// output the completed rows
	numer = numerator   + swinrh;
	denom = denominator + swinrh;

	int output_rows;
	if (row_cnt < swinrv) 
	{
		if (row_cnt + pstep <= swinrv) {
			// no row is completed in the begining
			output_rows = 0;
		}
		else {
			output_rows = row_cnt + pstep - swinrv;
			numer += (pstep - output_rows) * w;
			denom += (pstep - output_rows) * w;
		}
	} 
	else 
	{
		if (row_cnt >= orig_h - psize)
			output_rows = orig_h - row_cnt + swinrv;	// the last line of reference patches
		else
			output_rows = pstep;
		clean += orig_w * (row_cnt - swinrv);
	}

	for (int i = 0, r = 0; r < output_rows; r++)
	{
		for (int c = 0; c < orig_w; c++, i++)
		{
			clean[i] = (ImageType)(numer[c] / denom[c]);
		}
		numer += w;
		denom += w;
	}

	// remove the fisrt (pstep) rows of the numerator and denominator buffers
	// and insert (pstep) new rows to the end of the buffers
	shift_numer_denom();

	row_cnt += pstep;
	return output_rows;
}

void BM3D_WIE::grouping()
{
#pragma omp parallel for num_threads(USE_THREADS_NUM)
	for (int sy = -swinrv; sy <= swinrv; sy += sstepv)
	{
		for (int sx = -swinrh; sx <= swinrh; sx += ssteph)
		{
			int idx = (sy + swinrv) / sstepv * nsh + (sx + swinrh) / ssteph;
			for (int y = 0; y < psize; y++)
			{
				for (int x = psize - pstep; x < psize; x++)
				{
					dist_buf[nbuf * idx + (ncnt + x / pstep) % nbuf] += get_dist(refer_basic[y * w + x], refer_basic[(y + sy) * w + x + sx]);
				}
			}
			dist_sum[idx] += dist_buf[nbuf * idx + (ncnt - 2) % nbuf];
			dist_sum[idx] += dist_buf[nbuf * idx + (ncnt - 1) % nbuf];
		}
	}

	g3d_noisy->set_reference();
	g3d_basic->set_reference();

	for (int idx = 0, sy = -swinrv; sy <= swinrv; sy += sstepv)
	{
		for (int sx = -swinrh; sx <= swinrh; sx += ssteph, idx++)
		{
			g3d_basic->insert_patch(sx, sy, dist_sum[idx]);
		}
	}
	g3d_noisy->num = g3d_basic->num;
	for (int p = 0; p < g3d_basic->num; p++)
	{
		g3d_noisy->patch[p]->update(g3d_basic->patch[p]->x, g3d_basic->patch[p]->y, 0);
	}

#pragma omp parallel for num_threads(USE_THREADS_NUM)
	for (int i = 0; i < nsv; i++)
	{
		for (int j = 0; j < nsh; j++)
		{
			int idx = i * nsh + j;
			dist_sum[idx] -= dist_buf[nbuf * idx + (ncnt - 1) % nbuf];
			dist_sum[idx] -= dist_buf[nbuf * idx + (ncnt - 0) % nbuf];
			dist_buf[nbuf * idx + ncnt % nbuf] = 0;
		}
	}

	g3d_noisy->fill_patches_values(refer_noisy, w);
	g3d_basic->fill_patches_values(refer_basic, w);
}

void BM3D_WIE::filtering()
{
	g3d_noisy->transform_3d();
	g3d_basic->transform_3d();

	// wiener filtering
	double wie_wgt = 0.;
	for (int p = 0; p < g3d_basic->num; p++)
	{
		for (int i = 0; i < psize * psize; i++)
		{
			wie_wgt = (double)g3d_basic->patch[p]->values[i] * g3d_basic->patch[p]->values[i];
			wie_wgt = wie_wgt / (wie_wgt + g3d_basic->thres);
			wie_wgt_sum += wie_wgt;
			g3d_noisy->patch[p]->values[i] = (PatchType)((double)g3d_noisy->patch[p]->values[i] * wie_wgt);
		}
	}

	g3d_noisy->inv_transform_3d();
}

void BM3D_WIE::aggregation()
{
#if USE_INTEGER
	PatchType weight = 1;
#else
	PatchType weight = 64. / wie_wgt_sum;
#endif

	for (int p = 0; p < g3d_noisy->num; p++)
	{
		int x = g3d_noisy->patch[p]->x;
		int y = g3d_noisy->patch[p]->y;
		for (int i = 0, r = 0; r < psize; r++)
		{
			for (int c = 0; c < psize; c++, i++)
			{
				numer[(r + y) * w + c + x] += Kaiser[i] * weight * g3d_noisy->patch[p]->values[i];
				denom[(r + y) * w + c + x] += Kaiser[i] * weight;
			}
		}
	}
}

void BM3D_WIE::shift_numer_denom()
{
	numer = numerator;
	denom = denominator;
	for (int i = 0; i < 2 * swinrv + psize - pstep; i++)
	{
		memcpy(numer, numer + w * pstep, w * sizeof(PatchType));
		memcpy(denom, denom + w * pstep, w * sizeof(PatchType));
		numer += w;
		denom += w;
	}
	memset(numer, 0, w * pstep * sizeof(PatchType));
	memset(denom, 0, w * pstep * sizeof(PatchType));
}