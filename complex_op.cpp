#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "complex_op.h"

//复数乘法
Complex ComplexMul(Complex a, Complex b) {
    Complex c;
    c.real = a.real * b.real - a.img * b.img;
    c.img = a.img * b.real + a.real * b.img;
    return c;
}

//复数除法
Complex ComplexDiv(Complex a, Complex b) {
    Complex c;
    c.real = double(a.real * b.real + a.img * b.img) / double(b.real * b.real + b.img * b.img);
    c.img = double(a.img * b.real - a.real * b.img) / double(b.real * b.real + b.img * b.img);
    return c;
}

//复数平方根(算数)
Complex root(Complex x) {
    Complex r;
	//(a+ib)^2 = x
	//a2-b2=real;a2+b2=|x|
    double t = sqrt(x.real * x.real + x.img * x.img);
    r.real = sqrt((t + x.real) / 2.0);
    r.img = sqrt((t - x.real) / 2.0);
    return r;
}

double abs_C(Complex x) {
    return sqrt(x.real * x.real + x.img * x.img);
}

void real(Complex *x, double* out, unsigned num) {
    for (int i = 0; i < num; i++) {
        out[i] = double(x[i].real);
    }

}

void img(Complex* x, double* out, unsigned num) {
    for (int i = 0; i < num; i++) {
        out[i] = double(x[i].img);
    }
}

void init_complex(Complex *x, unsigned size) {
    for (unsigned i = 0; i < size; i++) {
        x[i].real = x[i].img = 0.0;
    }
}

void conj(Complex* x, Complex*y, unsigned size) {
    
    for (unsigned i = 0; i < size; i++) {
        y[i].real = x[i].real;
        y[i].img = -x[i].img;
    }
}

void Point_Mul(Complex* x, Complex* y, Complex* z, unsigned size) {
    for (unsigned i = 0; i < size; i++) {
        z[i] = ComplexMul(x[i], y[i]);
    }
}

void print_matrix(Complex* x, unsigned size1, unsigned size2) {
    printf("print matrix:\n");
    for (unsigned i = 0; i < size1; i++) {
        printf("%d: ", i+1);
        for (unsigned j = 0; j < size2; j++) {
            printf("%e+%e | ", x[i * size2 + j].real, x[i * size2 + j].img);
        }
        printf("\n\n");
    }
}

//取矩阵第k列
void get_column(Complex* x, Complex* y, unsigned size1, unsigned size2, unsigned k) {
    if (k > size2) {
        printf("out of size2\n");
    }
    for (unsigned i = 0; i < size1; i++) {
        y[i] = x[i * size2 + k];
    }
}

//取矩阵第k行
void get_row(Complex* x, Complex* y, unsigned size1, unsigned size2, unsigned k) {
    if (k > size1) {
        printf("out of size2\n");
    }
    for (unsigned i = 0; i < size2; i++) {
        y[i] = x[i + k*size2];
    }
}

//取矩阵第k列
void get_column_f(float* x, float* y, unsigned size2, unsigned k, unsigned start, unsigned stop) {
    for (unsigned i = start, j=0; i < stop; i++, j++) {
        y[j] = x[i * size2 + k];
    }
}

//取矩阵第k行
void get_row_f(float* x, float* y, unsigned size2, unsigned k, unsigned start, unsigned stop) {
    for (unsigned i = start, j=0; i < stop; i++, j++) {
        y[j] = x[i + k * size2];
    }
}

//均值
float mean(float *x, unsigned size) {
    float sum = 0;
    for (unsigned i = 0; i < size; i++) {
        sum += x[i];
    }
    return sum / float(size);
}
