#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "complex_op.h"
#include "eig.h"
#include "doa.h"
#include "tool.h"
#include "kfft.h"
#include "CFAR_OSCA.h"
#include "main.h"
#include "awgn.h"


int main()
{
	/*——CFAR——*/
	const float R_res = 1; //Range_Resolution
	const float V_max = 100;  //Max_Velocity
	const float c = 3e8;
	const float RCS = 100;
	int AA = 1;
	const double fc = 77e9;
	const float lambda = c / fc;

	const float B_sweep = c / (2 * R_res);
	const float fs = 2 * B_sweep;
	const float Ts = 1 / fs;
	
	const short Nd = 64; // slow - time dimension, num of chirps
	const short Nr = 128; // fast - time dimension, num of samples
	//range
	float T_c = Nr * Ts; //T of chirp
	float k = B_sweep / T_c;
	float R_max = c * fs / (2 * k);
	
	//doppler
	float PRF = 2 * V_max / lambda;
	float T_p = 1 / PRF;
	float V_res = lambda * PRF / (2.0 * Nd);
	float* tr = (float*)malloc(sizeof(float) * Nr); //快时间
	unsigned s_num = 3;

	Complex* Tx, * Rx, * Sr, * Sr1, * Mix, * fft1D;
	float* RDM;
	Tx = (Complex*)malloc(sizeof(Complex) * Nr * Nd);
	Rx = (Complex*)malloc(sizeof(Complex) * Nr * Nd * s_num);
	Sr = (Complex*)malloc(sizeof(Complex) * Nr * Nd);
	Sr1 = (Complex*)malloc(sizeof(Complex) * Nr * Nd);
	Mix = (Complex*)malloc(sizeof(Complex) * Nr * Nd);
	RDM = (float*)malloc(sizeof(float) * Nr * Nd); 
	

	if (!Tx || !Rx || !Sr || !Sr1 || !Mix || !tr ||!RDM) {
		printf("fail malloc");
		return 0;
	}
	for (int i = 0; i < Nr; i++) {
		tr[i] = (-Nr / 2 + i) / double(Nr / 2) * T_c / 2;
	}
	//生成CFAR测试信号
	int tarR[] = { 20,30,10 };
	int tarV[] = { 60,30,10 };
	init_complex(Sr, Nr * Nd);
	double t_all, R, delay;
	const unsigned P_received = 1;
	//加噪声
	//int SNR = 5;
	for (int n_tar = 0; n_tar < s_num; n_tar++) {
		for (int nd = 0; nd < Nd; nd++) {
			t_all = nd * T_p;
			R = tarR[n_tar] + t_all * tarV[n_tar];
			delay = 2 * R / c;
			
			for (int nr = 0; nr < Nr; nr++) {
				Tx[nr * Nd + nd].real = AA * cos(2 * pi * (fc * tr[nr] + 0.5 * k * tr[nr] * tr[nr]));
				Tx[nr * Nd + nd].img = AA * sin(2 * pi * (fc * tr[nr] + 0.5 * k * tr[nr] * tr[nr]));
				Rx[n_tar * Nr * Nd + nr * Nd + nd].real = AA * P_received * cos(2.0 * pi * (fc * (tr[nr] - delay) + 0.5 * k * (tr[nr] - delay) * (tr[nr] - delay)));
				Rx[n_tar * Nr * Nd + nr * Nd + nd].img = AA * P_received * sin(2.0 * pi * (fc * (tr[nr] - delay) + 0.5 * k * (tr[nr] - delay) * (tr[nr] - delay)));
				Sr[nr * Nd + nd].real += Rx[n_tar * Nr * Nd + nr * Nd + nd].real;
				Sr[nr * Nd + nd].img += Rx[n_tar * Nr * Nd + nr * Nd + nd].img;
				//printf("Rx: %.4lf+%.4lf \n", Rx[n_tar * Nr * Nd + nr * Nd + nd].real, Rx[n_tar * Nr * Nd + nr * Nd + nd].img);
			}
			//scanf_s("end");
		}
	}
	Sr = awgn(Sr, Nr*Nd, 5, 20);

	conj(Sr, Sr1, Nr * Nd);
	Point_Mul(Tx, Sr1, Mix, Nr * Nd);
	/*Range FFT*/
	fftshift_matrix(Mix, Nr, Nd, 1);
	fft_matrix(Mix, Nr, Nd, 1);
	fftshift_matrix(Mix, Nr, Nd, 1); //fft1D -> Mix
	/*Doppler FFT*/
	fftshift_matrix(Mix, Nr, Nd, 2);
	fft_matrix(Mix, Nr, Nd, 2);
	fftshift_matrix(Mix, Nr, Nd, 2); //fft2D -> Mix
	/*RDM*/
	fftshift_matrix(Mix, Nr, Nd, 1);
	fftshift_matrix(Mix, Nr, Nd, 2);
	free_n(Tx);
	free_n(Rx);
	free_n(Sr);
	free_n(Sr1);
	float sqrt_out;
	for (unsigned i = 0; i < Nr;i++) {
		//print("%d: ", i);
		for (unsigned j = 0; j < Nd; j++) {
			sqrt_out = abs_C(Mix[i * Nd + j]);
			RDM[i * Nd + j] = sqrt_out * sqrt_out;
			//print("%e  ", RDM[i * Nd + j]);
		}
		//print("\n");
	}
	//print_matrix(Mix, Nr, Nd);
	OSCA_OBJ obj = form_OSCA_obj();
	unsigned* Target_ind = NULL, Target_num;

	Target_ind = CFAR_OSCA_Range_2D(obj, RDM, Nr, Nd, &Target_num);
	printf("Target_num: %d\n", Target_num);
	for (unsigned i = 0; i < Target_num; i++) {
		printf("%d: %.4f %.4f\n", i, Target_ind[i * 2] * R_res, Target_ind[i * 2 + 1] * V_res);
	}
	free_n(Mix);
	free_n(RDM);
	

	/*——DOA——*/
	const unsigned NUM = 100; /*最大迭代次数*/
	int N = 8; /*阵元数*/
	int M = 3; /*信源数*/
	double theta[] = { 24, 21, 14 };
	int K = 512; /*快排数*/
	int snr = 10; /*信噪比*/
	double dd = 0.5; /*阵元间隔*/
	double derad = pi / 180; /*弧度*/
	double temp;
	int flag = 1;
	double* out;
	double* theta_doa = (double*)malloc(sizeof(double) * 360);
	double* d = (double*)malloc(sizeof(double) * N);
	Complex* A = (Complex*)malloc(sizeof(Complex) * N * M); /*方向矢量*/
	double* S = (double*)malloc(sizeof(double) * K * M);
	Complex* X = (Complex*)malloc(sizeof(Complex) * N * K);
	if (!d || !A || !S || !X || !theta_doa) {
		printf("fail malloc");
		return 0;
	}

	for (int i = 0; i < N; i++) {
		d[i] = i * dd;
		for (int j = 0; j < M; j++) {
			A[i * M + j].real = cos(d[i] * 2 * pi * sin(theta[j] * derad));
			A[i * M + j].img = -sin(d[i] * 2 * pi * sin(theta[j] * derad));
			//printf("%.2lf + %.2lfi\n", A[i * M + j].real,A[i * M + j].img);
		}
	}
	/*高斯分布*/
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < K; j++) {
			S[i * K + j] = randn();
		}
	}
	/* 构造信号X=A*S */
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < K; j++) {
			X[i * K + j] = { 0, 0 };
			for (int m = 0; m < M; m++) {
				X[i * K + j].real += A[i * M + m].real * S[m * K + j];
				X[i * K + j].img += A[i * M + m].img * S[m * K + j];
			}
			//printf("%.2lf + %.2lfi\n", X[i * K + j].real, X[i * K + j].img);
		}
	}
	
	out = DOA(X, N, M, K, dd);
	unsigned* local_out, local_num;  //find the local max location
	float* DOA_out;
	local_out = local_max(out, 360, &local_num);
	DOA_out = (float*)malloc(sizeof(float)*local_num);
	if (!DOA_out) {
		printf("malloc failed!\n");
	}
	else {
		for (unsigned i = 0; i < local_num; i++) {
			DOA_out[i] = (local_out[i] - 180) / 2;
			printf("DOA(%d): %.2lf\n", i, DOA_out[i]);
		}
		free_n(DOA_out);
	}
	free_n(local_out);
	free_n(d);
	free_n(A); 
	free_n(S);
	free_n(X);
	free_n(Target_ind);
	return 0;
}