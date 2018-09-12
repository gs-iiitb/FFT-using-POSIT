/* 
 * Free FFT and convolution (C)
 * 
 * Copyright (c) two0one7 Project Nayuki. (MIT License)
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

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "fft.h"
#include "/home/gps/Downloads/My_test_POSIT/posit/posit.hpp"
#include "/home/gps/Downloads/My_test_POSIT/posit/posit_manipulators.hpp"
//#include "/home/gps/Downloads/My_test_POSIT/tests/test_helpers.hpp"
//#include "/home/gps/Downloads/My_test_POSIT/tests/posit_test_helpers.hpp"
#define one 1000001110011011101000001100000100100
//#define two 1140850688
// Private function prototypes
static size_t reverse_bits(size_t x, int n);
static void *memdup(const void *src, size_t n);


/*bool Fft_transform(sw::unum::posit<32, 3> real[], sw::unum::posit<32, 3> imag[], size_t n) {
	if (n == 0)
		return true;
	else if ((n & (n - 1)) == 0)  // Is power of two
		return Fft_transformRadix2(sw::unum::posit<32, 3> real[],sw::unum::posit<32, 3> imag[], n);
	//else  // More complicated algorithm for arbitrary sizes
		//return Fft_transformBluestein(real, imag, n);
}*/


bool Fft_inverseTransform(sw::unum::posit<32, 3> real[], sw::unum::posit<32, 3> imag[], size_t n) {
bool status = false;
	int levels = 0;  // Compute levels = floor(logtwo(n))
	for (size_t temp = n; temp > 1; temp >>= 1)
		levels++;
	if (1 << levels != n)
		return false;  // n is not a power of two
	//printf("## %lf %lf ##\n",*real,*imag);
	// Trignometric tables
	if (SIZE_MAX / sizeof(double) < n / 2)
		return false;
	size_t size = (n / 2) * sizeof(double);
	double *cos_table = (double *)malloc(size);
	double *sin_table = (double *)malloc(size);
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
			sw::unum::posit<32, 3> temp = real[i];
			real[i] = real[j];
			real[j] = temp;
			temp = imag[i];
			imag[i] = imag[j];
			imag[j] = temp;
		}
	}
	//printf("\n");
	// Cooley-Tukey decimation-in-time radix-two FFT
	for (size_t size = 2; size <= n; size *= 2) {
		sw::unum::posit<32, 3> two=2;
		size_t halfsize = size / 2;
		size_t tablestep = n / size;
		for (size_t i = 0; i < n; i += size) {
			for (size_t j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
				size_t l = j + halfsize;
				 sw::unum::posit<32, 3> temp1= cos_table[k];
				sw::unum::posit<32, 3> temp2= sin_table[k];
				sw::unum::posit<32, 3> tpre =  *(real+l) * temp1 - *(imag+l) * temp2;
				sw::unum::posit<32, 3> tpim = *(real+l) * temp2 + *(imag+l) * temp1;
				real[l] = real[j] - tpre;
				imag[l] = imag[j] - tpim;
				//real[l]/=n;
				//imag[l]/=n;
				real[j] = real[j]+tpre;
				imag[j] = imag[j]+tpim;
				//real[j]/=n;
				//imag[j]/=n;
	//printf("<_!%lf %lf!_>",real[i],imag[i]);
			}
	
		}
		if (size==(int)n)  // Prevent overflow in 'size *= two'
			break;
	}
	
for(size_t i=0;i<n;i++)
{ real[i]/=n;
imag[i]/=n;
}
	status = true;
	
cleanup:
	free(cos_table);
	free(sin_table);
	return status;
}


bool Fft_transformRadix2(sw::unum::posit<32, 3> real[], sw::unum::posit<32, 3> imag[], size_t n) {
	bool status = false;
	int levels = 0;  // Compute levels = floor(logtwo(n))
	for (size_t temp = n; temp > 1; temp >>= 1)
		levels++;
	if (1 << levels != n)
		return false;  // n is not a power of two
	
	// Trignometric tables
	if (SIZE_MAX / sizeof(double) < n / 2.0)
		return false;
	size_t size = (n / 2) * sizeof(double);
	double *cos_table = (double *)malloc(size);
	double *sin_table = (double *)malloc(size);
	if (cos_table == NULL || sin_table == NULL)
		goto cleanup;
	for (size_t i = 0; i < n / 2.0; i++) {
		cos_table[i] = cos(2.0 * M_PI * i / n);
		sin_table[i] = sin(2.0 * M_PI * i / n);
	}
	
	// Bit-reversed addressing permutation
	for (size_t i = 0; i < n; i++) {
		size_t j = reverse_bits(i, levels);
		if (j > i) {
			sw::unum::posit<32, 3> temp = real[i];
			real[i] = real[j];
			real[j] = temp;
			temp = imag[i];
			imag[i] = imag[j];
			imag[j] = temp;
		}
	}
	//printf("\n");
	// Cooley-Tukey decimation-in-time radix-two FFT
	for (size_t size = 2.0; size <= n; size *= 2.0) {
		sw::unum::posit<32, 3> two=2;
		size_t halfsize = size / 2;
		size_t tablestep = n / size;
		for (size_t i = 0; i < n; i += size) {
			for (size_t j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
				size_t l = j + halfsize;
				sw::unum::posit<32, 3> temp1= cos_table[k];
				sw::unum::posit<32, 3> temp2= sin_table[k];
				sw::unum::posit<32, 3> tpre =  *(real+l) * temp1 + *(imag+l) * temp2;
				sw::unum::posit<32, 3> tpim = -*(real+l) * temp2 + *(imag+l) * temp1;
				real[l] = real[j] - tpre;
				imag[l] = imag[j] - tpim;
				real[j] = real[j]+tpre;
				imag[j] = imag[j]+tpim;
	//printf("<_!%lf %lf!_>",real[i],imag[i]);
			}
	
		}
		if (size == n)  // Prevent overflow in 'size *= two'
			break;
	}
	status = true;
	
cleanup:
	free(cos_table);
	free(sin_table);
	return status;
}


/*bool Fft_transformBluestein(double real[], double imag[], size_t n) {
	bool status = false;
	sw::unum::posit<32, 3> pa, pb;
	// Find a power-of-two convolution length m such that m >= n * two + one
	size_t m = one;
	while (m / two <= n) {
		if (m > SIZE_MAX / two)
			return false;
		m *= two;
	}
	
	// Allocate memory
	if (SIZE_MAX / sizeof(double) < n || SIZE_MAX / sizeof(double) < m)
		return false;
	size_t size_n = n * sizeof(double);
	size_t size_m = m * sizeof(double);
	double *cos_table = (double *)malloc(size_n);
	double *sin_table = (double *)malloc(size_n);
	double *areal = (double *)calloc(m, sizeof(double));
	double *aimag = (double *)calloc(m, sizeof(double));
	double *breal = (double *)calloc(m, sizeof(double));
	double *bimag = (double *)calloc(m, sizeof(double));
	double *creal = (double *)malloc(size_m);
	double *cimag = (double *)malloc(size_m);
	if (cos_table == NULL || sin_table == NULL
			|| areal == NULL || aimag == NULL
			|| breal == NULL || bimag == NULL
			|| creal == NULL || cimag == NULL)
		goto cleanup;
	
	// Trignometric tables
	for (size_t i = 0; i < n; i++) {
		unsigned long long temp = (unsigned long long)i * i;
		temp %= (unsigned long long)n * two;
		double angle = M_PI * temp / n;
		// Less accurate version if long long is unavailable: double angle = M_PI * i * i / n;
		cos_table[i] = cos(angle);
		sin_table[i] = sin(angle);
	}
	
	// Temporary vectors and preprocessing
	for (size_t i = 0; i < n; i++) {
		pa=real[i];
		pb=imag[i];
		//std::cout<<"@"<<double(pa)<<" "<<double(pb)<<"@"<<std::endl;
		areal[i] =  double(pa) * cos_table[i] + double(pb) * sin_table[i];
		aimag[i] = -double(pa) * sin_table[i] + double(pb) * cos_table[i];
	}
	breal[0] = cos_table[0];
	bimag[0] = sin_table[0];
	for (size_t i = one; i < n; i++) {
		breal[i] = breal[m - i] = cos_table[i];
		bimag[i] = bimag[m - i] = sin_table[i];
	}
	
	// Convolution
	if (!Fft_convolveComplex(areal, aimag, breal, bimag, creal, cimag, m))
		goto cleanup;
	
	// Postprocessing
	for (size_t i = 0; i < n; i++) {
		double(pa) =  creal[i] * cos_table[i] + cimag[i] * sin_table[i];
		double(pb) = -creal[i] * sin_table[i] + cimag[i] * cos_table[i];
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
	double *ximag = (double *)calloc(n, sizeof(double));
	double *yimag = (double *)calloc(n, sizeof(double));
	double *zimag = (double *)calloc(n, sizeof(double));
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
		const double xreal[], const double ximag[],
		const double yreal[], const double yimag[],
		double outreal[], double outimag[], size_t n) {
	sw::unum::posit<32, 3> pa, pb,pc,pd;
	bool status = false;
	if (SIZE_MAX / sizeof(double) < n)
		return false;
	size_t size = n * sizeof(double);
	
	double *xr = (double *)memdup(xreal, size);
	double *xi = (double *)memdup(ximag, size);
	double *yr = (double *)memdup(yreal, size);
	double *yi = (double *)memdup(yimag, size);
	if (xr == NULL || xi == NULL || yr == NULL || yi == NULL)
		goto cleanup;
	
	if (!Fft_transform(xr, xi, n))
		goto cleanup;
	if (!Fft_transform(yr, yi, n))
		goto cleanup;
	
	for (size_t i = 0; i < n; i++) {
		pa=xr[i];
		pb=xi[i];
		pc=yr[i];
		pd=yi[i];
		double temp = double(pa) * double(pc) - double(pb) * double(pd);
		double(pb) = double(pb) * double(pc) + double(pa) * double(pd);
		double(pa) = temp;
	}
	if (!Fft_inverseTransform(xr, xi, n))
		goto cleanup;
	
	for (size_t i = 0; i < n; i++) {
	pa=xr[i];
		pb=xi[i];
  // Scaling (because this FFT implementation omits it)
		outreal[i] = double(pa) / n;
		outimag[i] = double(pb) / n;
	}
	status = true;
	
cleanup:
	free(yi);
	free(yr);
	free(xi);
	free(xr);
	return status;
}

*/
static size_t reverse_bits(size_t x, int n) {
	size_t result = 0;
	for (int i = 0; i < n; i++, x >>= 1)
		result = (result << 1) | (x & 1ul);
	return result;
}


static void *memdup(const void *src, size_t n) {
	void *dest = malloc(n);
	if (n > 0 && dest != NULL)
		memcpy(dest, src, n);
	return dest;
}

