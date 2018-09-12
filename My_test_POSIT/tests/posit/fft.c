/* 
 * Free FFT and convolution (C)
 * 
 * Copyright (c) 2017 Project Nayuki. (MIT License)
 * https://www.nayuki.io/page/free-small-fft-in-multiple-languages
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 * - The above copyright notice and this permission notice shall be included in
 *   all copies or substantial portions of the Software.
 * - The Software is provided "as is", without warranty of any kind, express or
 *   implied, including but not limited to the warranties of merchantability,
 *   fitness for a particular purpose and noninfringement. In no event shall the
 *   authors or copyright holders be liable for any claim, damages or other
 *   liability, whether in an action of contract, tort or otherwise, arising from,
 *   out of or in connection with the Software or the use or other dealings in the
 *   Software.
 */
#include "/home/gps/Downloads/My_test_POSIT/tests/posit/common.hpp"
#include "../../posit/posit.hpp"
#include "../../posit/posit_manipulators.hpp"
#include "/home/gps/Downloads/My_test_POSIT/tests/test_helpers.hpp"
#include "/home/gps/Downloads/My_test_POSIT/tests/posit_test_helpers.hpp"
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "/home/gps/Downloads/My_test_POSIT/tests/fft.h"

sw::unum::posit<16, 1> pa[], pb[],pc[],pd[], pref[], pdif[];
// Private function prototypes
static size_t reverse_bits(size_t x, int n);
static void *memdup(const void *src, size_t n);


bool Fft_transform(pa, pb, size_t n) {
	if (n == 0)
		return true;
	else if ((n & (n - 1)) == 0)  // Is power of 2
		return Fft_transformRadix2(pa, pb, n);
	else  // More complicated algorithm for arbitrary sizes
		return Fft_transformBluestein(pa, pb, n);
}


bool Fft_inverseTransform(pa, pb, size_t n) {
	return Fft_transform(pb, pa, n);
}


bool Fft_transformRadix2(pa, pb, size_t n) {
	// Length variables
	bool status = false;
	int levels = 0;  // Compute levels = floor(log2(n))
	for (size_t temp = n; temp > 1U; temp >>= 1)
		levels++;
	if ((size_t)1U << levels != n)
		return false;  // n is not a power of 2
	
	// Trignometric tables
	if (SIZE_MAX / sizeof(double) < n / 2)
		return false;
	size_t size = (n / 2) * sizeof(double);
	double *cos_table = malloc(size);
	double *sin_table = malloc(size);
	if (cos_table == NULL || sin_table == NULL)
		goto cleanup;
	for (size_t i = 0; i < n / 2; i++) {
		cos_table[i] = cos(2 * M_PI * i / n);
		sin_table[i] = sin(2 * M_PI * i / n);
	}
	
	// Bit-reversed addressing permutation
	for (size_t i = 0; i < n; i++) {
		size_t j = reverse_bits(i, levels);
		if (j > i) {
			double temp = pa[i];
			pa[i] = pa[j];
			pa[j] = temp;
			temp = pb[i];
			pb[i] = pb[j];
			pb[j] = temp;
		}
	}
	
	// Cooley-Tukey decimation-in-time radix-2 FFT
	for (size_t size = 2; size <= n; size *= 2) {
		size_t halfsize = size / 2;
		size_t tablestep = n / size;
		for (size_t i = 0; i < n; i += size) {
			for (size_t j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
				size_t l = j + halfsize;
				double tpre =  pa[l] * cos_table[k] + pb[l] * sin_table[k];
				double tpim = -pa[l] * sin_table[k] + pb[l] * cos_table[k];
				pa[l] = pa[j] - tpre;
				pb[l] = pb[j] - tpim;
				pa[j] += tpre;
				pb[j] += tpim;
			}
		}
		if (size == n)  // Prevent overflow in 'size *= 2'
			break;
	}
	status = true;
	
cleanup:
	free(cos_table);
	free(sin_table);
	return status;
}


bool Fft_transformBluestein(pa, pb, size_t n) {
	bool status = false;
	
	// Find a power-of-2 convolution length m such that m >= n * 2 + 1
// Why do we need this part??
	size_t m = 1;
	while (m / 2 <= n) {
		if (m > SIZE_MAX / 2)
			return false;
		m *= 2;
	}
	
	// Allocate memory
	if (SIZE_MAX / sizeof(double) < n || SIZE_MAX / sizeof(double) < m)
		return false;
	size_t size_n = n * sizeof(double);
	size_t size_m = m * sizeof(double);
	double *cos_table = malloc(size_n);
	double *sin_table = malloc(size_n);
sw::unum::posit<16, 1> areal[m],aimag[m],breal[m],bimag[m],creal[size_m],cimag[size_m];
	/*double *areal = calloc(m, sizeof(double));
	double *aimag = calloc(m, sizeof(double));
	double *breal = calloc(m, sizeof(double));
	double *bimag = calloc(m, sizeof(double));
	double *creal = malloc(size_m);
	double *cimag = malloc(size_m);*/
	if (cos_table == NULL || sin_table == NULL
			|| areal == NULL || aimag == NULL
			|| breal == NULL || bimag == NULL
			|| creal == NULL || cimag == NULL)
		goto cleanup;
	
	// Trignometric tables
	for (size_t i = 0; i < n; i++) {
		unsigned long long temp = (unsigned long long)i * i;
		temp %= (unsigned long long)n * 2;
		double angle = M_PI * temp / n;
		// Less accurate version if long long is unavailable: double angle = M_PI * i * i / n;
		cos_table[i] = cos(angle);
		sin_table[i] = sin(angle);
	}
	
	// Temporary vectors and preprocessing
	for (size_t i = 0; i < n; i++) {
		areal[i] =  pa[i] * cos_table[i] + pb[i] * sin_table[i];
		aimag[i] = -pa[i] * sin_table[i] + pb[i] * cos_table[i];
	}
	breal[0] = cos_table[0];
	bimag[0] = sin_table[0];
	for (size_t i = 1; i < n; i++) {
		breal[i] = breal[m - i] = cos_table[i];
		bimag[i] = bimag[m - i] = sin_table[i];
	}
	
	// Convolution
	if (!Fft_convolveComplex(areal, aimag, breal, bimag, creal, cimag, m))
		goto cleanup;
	
	// Postprocessing
	for (size_t i = 0; i < n; i++) {
		pa[i] =  creal[i] * cos_table[i] + cimag[i] * sin_table[i];
		pb[i] = -creal[i] * sin_table[i] + cimag[i] * cos_table[i];
	}
	status = true;
	
	// Deallocation
cleanup:
	free(cimag);
	free(creal);
	free(bimag);
	free(breal);
	free(aimag);
	free(areal);
	free(sin_table);
	free(cos_table);
	return status;
}


bool Fft_convolveReal(const double x[], const double y[], double out[], size_t n) {
	bool status = false;
	double *ximag = calloc(n, sizeof(double));
	double *yimag = calloc(n, sizeof(double));
	double *zimag = calloc(n, sizeof(double));
	if (ximag == NULL || yimag == NULL || zimag == NULL)
		goto cleanup;
	
	status = Fft_convolveComplex(x, ximag, y, yimag, out, zimag, n);
cleanup:
	free(zimag);
	free(yimag);
	free(ximag);
	return status;
}


bool Fft_convolveComplex(
		pa, pb,
		pc, pd,
		pref, pdif, size_t n) {
	
	bool status = false;
	if (SIZE_MAX / sizeof(double) < n)
		return false;
	size_t size = n * sizeof(double);
	
	sw::unum::posit<16, 1>xr[size] = pa;
	sw::unum::posit<16, 1>xi[size] = pb;
	sw::unum::posit<16, 1>yr[size] = pc;
	sw::unum::posit<16, 1>yi[size] = pd;
	if (pa == NULL || pb == NULL || pc == NULL || pd == NULL)
		goto cleanup;
	
	if (!Fft_transform(pa, pb, n))
		goto cleanup;
	if (!Fft_transform(pc, pd, n))
		goto cleanup;
	
	for (size_t i = 0; i < n; i++) {
		double temp = pa[i] * pc[i] - pb[i] * pd[i];
		xi[i] = pb[i] * pc[i] + pa[i] * pd[i];
		xr[i] = temp;
	}
	if (!Fft_inverseTransform(pa, pb, n))
		goto cleanup;
	
	for (size_t i = 0; i < n; i++) {  // Scaling (because this FFT implementation omits it)
		pref[i] = pa[i] / n;
		pdif[i] = pb[i] / n;
	}
	status = true;
	
cleanup:
	free(yi);
	free(yr);
	free(xi);
	free(xr);
	return status;
}


static size_t reverse_bits(size_t x, int n) {
	size_t result = 0;
	for (int i = 0; i < n; i++, x >>= 1)
		result = (result << 1) | (x & 1U);
	return result;
}


static void *memdup(const void *src, size_t n) {
	void *dest = malloc(n);
	if (n > 0 && dest != NULL)
		memcpy(dest, src, n);
	return dest;
}
