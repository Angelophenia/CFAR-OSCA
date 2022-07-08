#pragma once
#ifndef EIG_H
#define EIG_H
#include <stdio.h>
#include "complex_op.h"

//��C���涨��Ĳ�������
typedef enum { False = 0, True = 1 }Bool;
//�������Ԫ�ص�����ΪmatrixType
typedef double matrixType;

//�˽ṹ��������ʾ��������rowΪ�У�columnΪ�У�heightΪ�ߣ�array������ž���Ԫ��(��һά��ģ���ά/��ά)
typedef struct
{
	unsigned  row, column, height;
	matrixType* array; //ʹ��ʱ�������*array���г�ʼ��
}Matrix;


typedef struct
{
	unsigned  row, column, height;
	Complex* array;
}Matrix_C;

//��������ڴ�
Bool SetMatrixSize(Matrix* matrix, const unsigned row, const unsigned column, const unsigned height);
//����Matrix�����е�����Ԫ�ش�СΪele
void SetMatrixEle(Matrix* matrix, matrixType ele);
//����Matrix�����е�����Ԫ�ش�СΪ0
void SetMatrixZero(Matrix* matrix);
//�жϾ����Ƿ�Ϊ�գ���Ϊ���򷵻�1�����򷵻�0
Bool IsNullMatrix(const Matrix* matrix);
//���پ��󣬼��ͷ�Ϊ����̬������ڴ�,��������ĳ������0
void DestroyMatrix(Matrix* matrix);
//������������Ԫ�ظ�������return row*column*height
unsigned MatrixCapacity(const Matrix* matrix);
//||matrix||_2  ��A�����2����
matrixType MatrixNorm2(const Matrix* matrix);
//matrixB = matrix(:,:,height)��������ά�����ĳ�㣬��matrixΪ��ά�����轫height����Ϊ0
Bool CopyMatrix(Matrix* matrixB, Matrix* matrix, unsigned height);
//matrixC = matrixA * matrixB
Bool MatrixMulMatrix(Matrix* matrixC, const Matrix* matrixA, const Matrix* matrixB);
//��vector������Ԫ������sign= 0 ʱΪ��������Ϊ����
void SortVector(Matrix* vector, int sign);
//��ӡ����
void PrintMatrix(const Matrix* matrix);
//��A�ֽ�ΪQ��R
void QR(Matrix* A, Matrix* Q, Matrix* R);
//��������ֵ����������
void Eigenvectors(Matrix* eigenVector, Matrix* A, Matrix* eigenValue);

Bool SetMatrixSize_C(Matrix_C* matrix, const unsigned row, const unsigned column, const unsigned height);
Bool IsNullMatrix_C(const Matrix_C* matrix);
void DestroyMatrix_C(Matrix_C* matrix);
Bool CopyMatrix_C(Matrix_C* matrixB, Matrix_C* matrix, unsigned height);
void PrintMatrix_C(const Matrix_C* matrix);
#endif


