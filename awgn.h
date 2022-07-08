#pragma once
#ifndef __AWGN_H__
#define __AWGN_H__
#include "complex_op.h"

double generateGaussianNoise(double mu, double sigma);
Complex* awgn(Complex* x, unsigned NrNd, double snr, double sigpower);

#endif