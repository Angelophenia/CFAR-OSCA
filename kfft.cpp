#include<math.h>
#include<stdlib.h>
#include "kfft.h"
#include "tool.h"
#include "complex_op.h"


/*void kfft(double* pr, double* pi, int n, double* fr, double* fi){ //k=log2(n)
	int it, m, is, i, j, nv, l0;
	double p, q, s, vr, vi, poddr, poddi;
	int k = log2(n);
	//蝶形运算倒序
	for (it = 0; it <= n - 1; it++)
	{
		m = it;
		is = 0;
		for (i = 0; i <= k - 1; i++)
		{
			j = m / 2;
			is = 2 * is + (m - 2 * j);
			m = j;
		}//(m - 2 * j)依次得到低位->高位
		fr[it] = pr[is];
		fi[it] = pi[is];
	}
	pr[0] = 1.0;
	pi[0] = 0.0;//W0N
	p = 6.283185306 / (1.0 * n);
	pr[1] = cos(p);
	pi[1] = -sin(p);//W1N

	for (i = 2; i <= n - 1; i++)//WiN
	{
		p = pr[i - 1] * pr[1];
		q = pi[i - 1] * pi[1];
		s = (pr[i - 1] + pi[i - 1]) * (pr[1] + pi[1]);
		pr[i] = p - q; pi[i] = s - p - q;
	}
	for (it = 0; it <= n - 2; it = it + 2)
	{
		vr = fr[it];
		vi = fi[it];
		fr[it] = vr + fr[it + 1];
		fi[it] = vi + fi[it + 1];
		fr[it + 1] = vr - fr[it + 1];
		fi[it + 1] = vi - fi[it + 1];
	}
	m = n / 2;
	nv = 2;
	for (l0 = k - 2; l0 >= 0; l0--)
	{
		m = m / 2;
		nv = 2 * nv;
		for (it = 0; it <= (m - 1) * nv; it = it + nv)
			for (j = 0; j <= (nv / 2) - 1; j++)
			{
				p = pr[m * j] * fr[it + j + nv / 2];
				q = pi[m * j] * fi[it + j + nv / 2];
				s = pr[m * j] + pi[m * j];
				s = s * (fr[it + j + nv / 2] + fi[it + j + nv / 2]);
				poddr = p - q;
				poddi = s - p - q;
				fr[it + j + nv / 2] = fr[it + j] - poddr;
				fi[it + j + nv / 2] = fi[it + j] - poddi;
				fr[it + j] = fr[it + j] + poddr;
				fi[it + j] = fi[it + j] + poddi;
			}
	}
	for (i = 0; i <= n - 1; i++)
	{
		pr[i] = sqrt(double(fr[i] * fr[i]) + double(fi[i] * fi[i]));  //幅值计算
	}
}*/
void kfft(double* Wr, double* Wi, int N, double* outr, double* outi) { //k=log2(n)
//void fft()
		int i = 0, j = 0, k = 0, l = 0;
		//蝶形运算倒序
		for (int it = 0; it < N; it++)
		{
			int m = it;
			int reverse = 0;
			for (i = 0; i < log(N) / log(2); i++)
			{
				j = m / 2;
				reverse = 2 * reverse + (m - 2 * j);
				m = j;
			}//(m - 2 * j)依次得到低位->高位
			outr[it] = Wr[reverse];
			outi[it] = Wi[reverse];
		}

		Wr[0] = 1.0;
		Wi[0] = 0.0;//W0N
		double p = 6.283185306 / (1.0 * N);//e(j*2*pi*n/N),p=j*2*pi*1/N
		Wr[1] = cos(p);
		Wi[1] = -sin(p);//W1N=e(jp)=cos(p)-j*sin(p)

		for (int i = 2; i < N; i++)//WiN
		{
			Wr[i] = Wr[i - 1] * Wr[1] - Wi[i - 1] * Wi[1];
			Wi[i] = Wr[i - 1] * Wi[1] + Wi[i - 1] * Wr[1];
		}

	
	double upr, downr, productr, upi, downi, producti; 
	for (i = 0; i < log(N) / log(2); i++)        /*一级蝶形运算 stage */
	{
		l = 1 << i;
		for (j = 0; j < N; j += 2 * l)     /*一组蝶形运算 group,每个group的蝶形因子乘数不同, j is the start of a group*/
		{
			for (k = 0; k < l; k++)        /*一个蝶形运算 每个group内的蝶形运算 共l个蝶形*/
			{
				//mul(x[j + k + l], W[size_x*k / 2 / l], &product);
				productr = outr[j + k + l] * Wr[k* N / 2 / l] - outi[j + k + l] * Wi[k* N / 2 / l];
				producti = outr[j + k + l] * Wi[k* N / 2 / l] + outi[j + k + l] * Wr[k* N / 2 / l];
				//add(x[j + k], product, &up);
				upr = outr[j + k] + productr;
				upi = outi[j + k] + producti;
				//sub(x[j + k], product, &down);
				downr = outr[j + k] - productr;
				downi = outi[j + k] - producti;
				//x[j + k] = up;
				outr[j + k] = upr;
				outi[j + k] = upi;
				//
				outr[j + k + l] = downr;
				outi[j + k + l] = downi;
			}
		}
	}
	for (i = 0; i < N; i++)
	{
		Wr[i] = sqrt(double(outr[i] * outr[i]) + double(outi[i] * outi[i]));  //幅值计算
	}
}


void fftshift(Complex *x, int length){
	int k = length / 2;
	Complex temp;
	for (int i = 0; i < length / 2; i++) {
		temp = x[k+i];
		x[k+i] = x[i];
		x[i] = temp;
	}
}

void fftshift_matrix(Complex* x, unsigned size1, unsigned size2, unsigned dim) {
	Complex temp;
	if (dim == 1) { //横对调
		//printf("fftshift_matrix, dim == 1\n");
		unsigned k = size1 / 2;
		for (unsigned i = 0; i < size2; i++) { //列

			for (unsigned j = 0; j < k; j++) {
				temp = x[j * size2 + i];
				x[j * size2 + i] = x[(j + k) * size2 + i];
				x[(j + k) * size2 + i] = temp;
			}
		}
	}
	else { //竖对调
		unsigned k = size2 / 2;
		for (unsigned i = 0; i < size1; i++) { //行
			for (unsigned j = 0; j < k; j++) {
				temp = x[i * size2 + j];
				x[i * size2 + j] = x[i * size2 + j + k];
				x[i * size2 + j + k] = temp;
			}
		}
	}
}


void fft_matrix(Complex *x, unsigned size1, unsigned size2, unsigned dim) {
	int n, k;
	Complex* temp;
	double* temp_r, * temp_i, * out_r, * out_i;
	if (dim == 1) {//每列做fft
		temp = (Complex*)malloc(sizeof(Complex) * size1);
		temp_r = (double*)malloc(sizeof(double) * size1);
		temp_i = (double*)malloc(sizeof(double) * size1);
		out_r = (double*)malloc(sizeof(double) * size1);
		out_i = (double*)malloc(sizeof(double) * size1);
		init_array(temp_r, size1);
		init_array(temp_i, size1);
		init_array(out_r, size1);
		init_array(out_i, size1);
		init_complex(temp, size1);
		n = size1;
		k = log2(n);
		for (unsigned k = 0; k < size2; k++) {
			get_column(x, temp, size1, size2, k);  //取第k列temp
			real(temp, temp_r, size1);  //取实部temp_r
			img(temp, temp_i, size1);  //取虚部temp_i
			/*for (int i = 0; i < 5; i++) {
				printf("%e + %ei | ", temp_r[i], temp_i[i]);
			}*/
			//printf("\n");
			kfft(temp_r, temp_i, size1, out_r, out_i); /* fft */
			for (unsigned i = 0; i < size1; i++) { //原地赋值
				x[i * size2 + k].real = out_r[i];
				x[i * size2 + k].img = out_i[i];
				//printf("%d %e + %ei \n",i+1, out_r[i], out_i[i]);

			}
			//scanf_s("end");
		}
	}
	else { //每行做fft
		temp = (Complex*)malloc(sizeof(Complex) * size2);
		temp_r = (double*)malloc(sizeof(double) * size2);
		temp_i = (double*)malloc(sizeof(double) * size2);
		out_r = (double*)malloc(sizeof(double) * size2);
		out_i = (double*)malloc(sizeof(double) * size2);
		init_array(out_r, size2);
		init_array(out_i, size2);
		init_complex(temp, size2);
		n = size2;
		k = log2(n);
		for (unsigned k = 0; k < size1; k++) {
			get_row(x, temp, size1, size2, k);
			real(temp, temp_r, size2);
			img(temp, temp_i, size2);
			kfft(temp_r, temp_i, size2, out_r, out_i); /* fft */
			for (unsigned i = 0; i < size2; i++) { //原地赋值
				x[k * size2 + i].real = out_r[i];
				x[k * size2 + i].img = out_i[i];
			}
		}
		
	}
	free_n(temp);
	free_n(temp_r);
	free_n(temp_i);
	free_n(out_r);
	free_n(out_i);
}


float abs_complex(Complex x){
	float y = sqrt(x.real*x.real + x.img*x.img);
	return y;
}

