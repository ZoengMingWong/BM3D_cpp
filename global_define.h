#ifndef __GLOBAL_DEFINE_H__
#define __GLOBAL_DEFINE_H__

#include <iostream>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

typedef uint8_t ImageType;				// data-type of the input/ouput image (up to 12 bits for integer version)
typedef uint32_t DistType;				// data-type of the distance between two patches

#define USE_INTEGER				1		// use integer or floating-point version
#define USE_L2_DIST				1		// use L2 (square) or L1 (absolute) distance

#define HARD_THRES_MULTIPLIER	2.7f	// multiply by the sigma is the threshold of the hard filtering 

#define USE_THREADS_NUM			4		// number of CPU threads can be used in the grouping step

#if USE_INTEGER

#define COEFF_DICI_BITS			1		// decimal bits of the interger coefficients (according to the transform implementation)
typedef int32_t PatchType;				// int32_t is enough for the intermediate results

#else

#define COEFF_DICI_BITS			0
typedef float PatchType;

#endif


#endif

