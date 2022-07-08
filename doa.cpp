#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "complex_op.h"
#include "eig.h"
#include "tool.h"
#define pi 3.1415926

/*
	Function: doa()
	Input:
		X1:接收信号
		N :阵元数
		M :信源数
		K :快拍数
		dd:阵元间隔
*/
double* DOA(Complex *X , int N, int M, int K, double dd)
{
	unsigned i, j, k;
	const unsigned NUM = 100; //最大迭代次数
	double theta[] = { -30, 0, 60 };
	double* angle = (double*)malloc(sizeof(double) * 361);
	Complex* Pmusic = (Complex*)malloc(sizeof(Complex) * 361);
	double* Pmusic_abs = (double*)malloc(sizeof(double) * 361);
	double* d = (double*)malloc(sizeof(double) * N);
	Complex* a = (Complex*)malloc(sizeof(Complex) * N);
	double* RA = (double*)malloc(sizeof(double) * N * N);
	double* RB = (double*)malloc(sizeof(double) * N * N);
	double* Rxx = (double*)malloc(sizeof(double) * 2 * N * 2 * N);
	//double* Rt = (double*)malloc(sizeof(double) * 2 * N * 2 * N);
	double rr = rand();
	if (!angle || !Pmusic || !Pmusic_abs || !a || !d || !RA || !RB || !Rxx) {
		printf("Fail malloc!");
		return 0;
	}
	for (int k = 0; k < N; k++) {
		d[k] = k * dd;
	}

	Matrix A, Q, R, temp, eValue, eValue1, Q1;
	Matrix_C Q_C, Q_C1;
	//分配内存
	SetMatrixSize(&A, 2 * N, 2 * N, 1);
	SetMatrixSize(&Q, A.row, A.column, A.height);
	SetMatrixSize(&R, A.row, A.column, A.height);
	SetMatrixSize(&temp, A.row, A.column, A.height);
	SetMatrixSize(&eValue, A.row, 1, 1);
	SetMatrixSize(&eValue1, A.row / 2, 1, 1);
	SetMatrixSize(&Q1, 2 * N, N, 1);
	SetMatrixSize_C(&Q_C, N, N, 1);
	SetMatrixSize_C(&Q_C1, N, N - M, 1);

	double derad = pi / 180;
	double phim;
	Complex temp_1[5], temp_2[8], temp_3 = { 0,0 }, tt, tt1;

	/*初始化Q_C*/
	for (int i = 0; i < N * N; i++) {
		Q_C.array[i].real = Q_C.array[i].img = 0;
	}

	/* 求信号X的协方差矩阵: Rxx=X*X'/K */
	Complex X_t, *X1;
	double RA_sum=0, RB_sum=0;
	//complex *X_m = (complex*)malloc(sizeof(complex) * N);
	X1 = (Complex*)malloc(sizeof(Complex) * N * N);
	for (int i = 0; i < N;i++) {
		for (int j = 0; j < N;j++) {
			RA[i * N + j] = RB[i * N + j] = X1[i * N + j].real = X1[i * N + j].img = 0;
			for (int k = 0; k < K; k++) {
				X_t.real = X[j * K + k].real;
				X_t.img = -X[j * K + k].img;
				tt = ComplexMul(X[i * K + k], X_t); //转置虚部取逆！
				RA[i * N + j] += tt.real;
				RB[i * N + j] += tt.img;
				//printf("%.3lf + %.3lfi\n", tt.real, tt.img);
			}
			RA[i * N + j] /= K;
			RB[i * N + j] /= K;
			//printf("%d %d : %.3lf",i,j, RA[i * N + j]);
		}
		//printf("\n");
	}
	/*Toeplitz*/
	//for (int i = 0; i < N; i++) {
	//	X_m[i].real = X_m[i].img = 0;
	//	for (int j = 0; j < N - i; j++) {
	//		X_m[i].real += RA[j + i + j * N];
	//		X_m[i].img += RB [j + i + j * N];
	//	}
	//	X_m[i].real /= (N-i);
	//	X_m[i].img /= (N-i);
	//	//printf("X_m: %.2lf, %.2lf\n", X_m[i].real, X_m[i].img);
	//}
	//X1 = toeplitz(X_m, N);
	//for (int i = 0; i < N; i++) {
	//	for (int j = 0; j < N; j++) {
	//		printf("X1: %.2lf, %.2lf", X1[i * N + j].real, X1[i * N + j].img);
	//	}
	//	printf("\n");
	//}
	//for (int i = 0; i < N; i++) {
	//	for (int j = 0; j < N; j++) {
	//		Rt[i * 2 * N + j] = X1[i * N + j].real;
	//		Rt[i * 2 * N + (j + N)] = -X1[i * N + j].img;
	//		Rt[(i + N) * 2 * N + j] = X1[i * N + j].img;
	//		Rt[(i + N) * 2 * N + (j + N)] = X1[i * N + j].real;
	//	}
	//}

	/*Rxx = [RA, -RB;RB, RA]*/
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			Rxx[i * 2 * N + j] = RA[i * N + j];
			Rxx[i * 2 * N + (j + N)] = -RB[i * N + j];
			Rxx[(i + N) * 2 * N + j] = RB[i * N + j];
			Rxx[(i + N) * 2 * N + (j + N)] = RA[i * N + j];
		}
	}

	A.array = Rxx; //数组在栈上，不能主动销毁

	/*拷贝A矩阵元素至temp*/
	CopyMatrix(&temp, &A, 0);

	/*初始化Q、R所有元素为0*/
	SetMatrixZero(&Q);
	SetMatrixZero(&R);

	/*使用QR分解求矩阵特征值*/
	for (k = 0;k < NUM; ++k)
	{
		QR(&temp, &Q, &R);
		MatrixMulMatrix(&temp, &R, &Q);
	}

	/*获取特征值，将之存储于eValue*/
	for (k = 0;k < temp.column;++k)
	{
		eValue.array[k] = temp.array[k * temp.column + k];
	}

	/*对特征值按照降序排序*/
	SortVector(&eValue, 1);
	for (k = 0; k < temp.column / 2; ++k)
	{
		eValue1.array[k] = eValue.array[k * 2];
	}

	/*根据特征值eValue，原始矩阵求解矩阵特征向量Q*/
	SetMatrixZero(&Q);
	Eigenvectors(&Q, &A, &eValue);
	for (int i = 0; i < 16; i++) {
		for (k = 0; k < 8; ++k)
		{
			Q1.array[i * 8 + k] = Q.array[i * 16 + k * 2]; //去除重复特征的特征向量
		}
	}
	/* [u,v]->v+ui*/
	for (int i = 0; i < 8;i++) {
		for (int k = 0; k < 8; k++) {
			Q_C.array[i * 8 + k].img = Q1.array[i * 8 + k];
			Q_C.array[i * 8 + k].real = -Q1.array[(i+8) * 8 + k];
		}
	}
	//printf("特征值");
	//PrintMatrix(&eValue1);
	//printf("特征向量Q1");
	//PrintMatrix(&Q1);
	//printf("特征向量Q_C");
	//PrintMatrix_C(&Q_C);
	
	//取后(N-M)列做噪声子空间
	for (int i = 0; i < N; i++) {
		for (int j = M, k=0; j < N; j++, k++) {
			Q_C1.array[i * Q_C1.column + k] = Q_C.array[i * Q_C.column + j];
		}
	}
	//printf("噪声子空间");
	//PrintMatrix_C(&Q_C1);

	/*谱峰搜索*/
	//printf("谱峰搜索");
	double abs_out, max_Pmusic=0;
	for (int i = 0; i < 361;i++) {
		angle[i] = double(i - 180) / 2.0;
		phim = derad * angle[i];
		for (int j = 0; j < N; j++) {
			/*欧拉公式*/
			a[j].real = cos(2 * pi * d[j] * sin(phim)); 
			a[j].img = -sin(2 * pi * d[j] * sin(phim));
		}
		/*temp_1=a'*Q_C1*/
		for (int j = 0; j < N-M; j++) { //列
			temp_1[j].real = temp_1[j].img = 0;
			for (k = 0; k < N; k++) {  //行
				tt = ComplexMul(a[k], Q_C1.array[k * Q_C1.column + j]);
				temp_1[j].real += tt.real;
				temp_1[j].img += tt.img;
				
			}
			//printf("%d temp_1: %.2lf %.2lf\n", j, temp_1[j].real, temp_1[j].img);
		}
		//scanf_s("end");
		/*temp2=a'*Q_C1*Q_C1'*/
		for (int k = 0; k < N; k++) {
			temp_2[k].real = temp_2[k].img = 0;
			for (int j = 0; j < N - M; j++) {
				tt1.real = Q_C1.array[j + Q_C1.column * k].real;
				tt1.img = -Q_C1.array[j + Q_C1.column * k].img; //复矩阵转置虚部取反
				tt = ComplexMul(temp_1[j], tt1);
				temp_2[k].real += tt.real;
				temp_2[k].img += tt.img;
			}
		}
		/*temp3=a'*Q_C1*Q_C1'*a*/
		temp_3 = { 0,0 };
		for (int j = 0; j < N; j++) {
			tt1.real = a[j].real;
			tt1.img = -a[j].img;
			tt = ComplexMul(temp_2[j], tt1);
			temp_3.real += tt.real;
			temp_3.img += tt.img;
		}
		abs_out = abs_C(temp_3);
		Pmusic[i].real = temp_3.real / (abs_out * abs_out);
		Pmusic[i].img = -temp_3.img / (abs_out * abs_out);
		Pmusic_abs[360-i] = abs_C(Pmusic[i]);  //掉头
		if (Pmusic_abs[i] > max_Pmusic)max_Pmusic = Pmusic_abs[i];
	}
	/*归一化*/
	for (int i = 0; i < 361; i++) {
		Pmusic_abs[i] = 10 * log10(Pmusic_abs[i] / max_Pmusic);
		//printf("%d Pmusic_abs: %f\n", i, Pmusic_abs[i]);
	}

	DestroyMatrix(&R);
	DestroyMatrix(&Q);
	DestroyMatrix(&Q1);
	DestroyMatrix(&eValue);
	DestroyMatrix(&eValue1);
	DestroyMatrix(&temp);
	DestroyMatrix_C(&Q_C);
	DestroyMatrix_C(&Q_C1);
	return Pmusic_abs;
}