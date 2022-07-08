#include "awgn.h"
//引入时间库 和 随机函数库 
#include<time.h>
#include<stdlib.h>
#include <math.h>

#define PI 3.1415926

//step1 生成服从U(0, 1)分布的u1, u2;
//step2 令 y = [-2 * ln(u1)] ^ 0.5*sin(2 * pi*u2);
//step3 令 x = miu + y * sigma，其中miu为均值，sigma为标准差
double generateGaussianNoise(double mu=0, double sigma=1)
{
	double X;
	double U1, U2;
	//均匀分布[0,1]
	U1 = (double)rand() / RAND_MAX;
	U2 = (double)rand() / RAND_MAX;
	X = sqrt(-2 * log(U1))*sin(2 * PI * U2);
	return mu + sigma * X;
}

Complex* awgn(Complex* x, unsigned NrNd, double snr, double sigpower=0) {
	// x:	信号，Nr*Nd
	// snr:	每个采样点的信噪比10lg(S/N)，dB
	// sigpower:	x的能量dBW，
		// 设为-1时自动计算： sigPower = sum(abs(sig(:)).^2)/length(sig(:));
	if (-1 == sigpower) {
		sigpower = 0;
		for (unsigned i = 0; i < NrNd; i++) {
			sigpower += pow(x[i].real, 2) + pow(x[i].img, 2);
		}
		sigpower = sigpower / NrNd;
	}
	else {
		sigpower = pow(10, sigpower / 10);	//dB转linear
	}
	
	snr = pow(10, snr / 10);	//dB转linear
	double noisePower = sigpower / snr;

	/*Complex noise;
	noise.real = sqrt(noisePower / 2)* generateGaussianNoise();
	noise.img = sqrt(noisePower / 2)* generateGaussianNoise();*/

	srand((unsigned)(time(NULL)));//调用srand是每次产生的随机数不同 
	for (unsigned i = 0; i < NrNd; i++) {
		x[i].real += sqrt(noisePower / 2)* generateGaussianNoise();
		x[i].img += sqrt(noisePower / 2)* generateGaussianNoise();
	}
	return x;
	
}
