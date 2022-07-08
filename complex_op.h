#pragma once
#ifndef COMPLEX_OP_H
#define COMPLEX_OP_H
#include <math.h>
#include <stdio.h>

typedef struct{
    double real;
    double img;
}Complex;

//复数乘法
Complex ComplexMul(Complex a, Complex b);

//复数除法
Complex ComplexDiv(Complex a, Complex b);

//复数平方根(算数)
Complex root(Complex x);

//复数取模
double abs_C(Complex x);

//取实部
void real(Complex* x, double* out, unsigned num);
//取虚部
void img(Complex* x, double* out, unsigned num);
//初始化complex数组
void init_complex(Complex* x, unsigned size);
//complex数组取共轭
void conj(Complex* x, Complex* y, unsigned size);
//complex数组点乘
void Point_Mul(Complex* x, Complex* y, Complex* z, unsigned size);
//打印complex数组
void print_matrix(Complex* x, unsigned size1, unsigned size2);
//取列
void get_column(Complex* x, Complex* y, unsigned size1, unsigned size2, unsigned k);  //复数矩阵
void get_column_f(float* x, float* y, unsigned size2, unsigned k, unsigned start, unsigned stop);  //实数矩阵, 截取
//取行
void get_row(Complex* x, Complex* y, unsigned size1, unsigned size2, unsigned k);
void get_row_f(float* x, float* y, unsigned size2, unsigned k, unsigned start, unsigned stop);
//均值
float mean(float* x, unsigned size);
#endif

