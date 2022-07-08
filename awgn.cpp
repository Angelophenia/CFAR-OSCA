#include "awgn.h"
//����ʱ��� �� ��������� 
#include<time.h>
#include<stdlib.h>
#include <math.h>

#define PI 3.1415926

//step1 ���ɷ���U(0, 1)�ֲ���u1, u2;
//step2 �� y = [-2 * ln(u1)] ^ 0.5*sin(2 * pi*u2);
//step3 �� x = miu + y * sigma������miuΪ��ֵ��sigmaΪ��׼��
double generateGaussianNoise(double mu=0, double sigma=1)
{
	double X;
	double U1, U2;
	//���ȷֲ�[0,1]
	U1 = (double)rand() / RAND_MAX;
	U2 = (double)rand() / RAND_MAX;
	X = sqrt(-2 * log(U1))*sin(2 * PI * U2);
	return mu + sigma * X;
}

Complex* awgn(Complex* x, unsigned NrNd, double snr, double sigpower=0) {
	// x:	�źţ�Nr*Nd
	// snr:	ÿ��������������10lg(S/N)��dB
	// sigpower:	x������dBW��
		// ��Ϊ-1ʱ�Զ����㣺 sigPower = sum(abs(sig(:)).^2)/length(sig(:));
	if (-1 == sigpower) {
		sigpower = 0;
		for (unsigned i = 0; i < NrNd; i++) {
			sigpower += pow(x[i].real, 2) + pow(x[i].img, 2);
		}
		sigpower = sigpower / NrNd;
	}
	else {
		sigpower = pow(10, sigpower / 10);	//dBתlinear
	}
	
	snr = pow(10, snr / 10);	//dBתlinear
	double noisePower = sigpower / snr;

	/*Complex noise;
	noise.real = sqrt(noisePower / 2)* generateGaussianNoise();
	noise.img = sqrt(noisePower / 2)* generateGaussianNoise();*/

	srand((unsigned)(time(NULL)));//����srand��ÿ�β������������ͬ 
	for (unsigned i = 0; i < NrNd; i++) {
		x[i].real += sqrt(noisePower / 2)* generateGaussianNoise();
		x[i].img += sqrt(noisePower / 2)* generateGaussianNoise();
	}
	return x;
	
}
