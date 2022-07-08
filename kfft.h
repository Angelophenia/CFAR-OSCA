#pragma once
#ifndef __KFFT_H__
#define __KFFT_H__
#include "complex_op.h"


void kfft(double* pr, double* pi, int n, double* fr, double* fi);
void fftshift(Complex *x, int length);
void fftshift_matrix(Complex* x, unsigned size1, unsigned size2, unsigned dim = 1);
void fft_matrix(Complex* x, unsigned size1, unsigned size2, unsigned dim = 1);
float abs_complex(Complex x);
#endif
