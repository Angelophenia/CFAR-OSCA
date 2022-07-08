#pragma once
#ifndef COMPLEX_OP_H
#define COMPLEX_OP_H
#include <math.h>
#include <stdio.h>

typedef struct{
    double real;
    double img;
}Complex;

//�����˷�
Complex ComplexMul(Complex a, Complex b);

//��������
Complex ComplexDiv(Complex a, Complex b);

//����ƽ����(����)
Complex root(Complex x);

//����ȡģ
double abs_C(Complex x);

//ȡʵ��
void real(Complex* x, double* out, unsigned num);
//ȡ�鲿
void img(Complex* x, double* out, unsigned num);
//��ʼ��complex����
void init_complex(Complex* x, unsigned size);
//complex����ȡ����
void conj(Complex* x, Complex* y, unsigned size);
//complex������
void Point_Mul(Complex* x, Complex* y, Complex* z, unsigned size);
//��ӡcomplex����
void print_matrix(Complex* x, unsigned size1, unsigned size2);
//ȡ��
void get_column(Complex* x, Complex* y, unsigned size1, unsigned size2, unsigned k);  //��������
void get_column_f(float* x, float* y, unsigned size2, unsigned k, unsigned start, unsigned stop);  //ʵ������, ��ȡ
//ȡ��
void get_row(Complex* x, Complex* y, unsigned size1, unsigned size2, unsigned k);
void get_row_f(float* x, float* y, unsigned size2, unsigned k, unsigned start, unsigned stop);
//��ֵ
float mean(float* x, unsigned size);
#endif

