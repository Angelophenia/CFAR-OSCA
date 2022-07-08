#pragma once
#ifndef EIG_H
#define EIG_H
#include <stdio.h>
#include "complex_op.h"

//纯C里面定义的布尔类型
typedef enum { False = 0, True = 1 }Bool;
//定义矩阵元素的类型为matrixType
typedef double matrixType;

//此结构体用来表示矩阵，其中row为行，column为列，height为高，array用来存放矩阵元素(用一维来模拟二维/三维)
typedef struct
{
	unsigned  row, column, height;
	matrixType* array; //使用时，必须对*array进行初始化
}Matrix;


typedef struct
{
	unsigned  row, column, height;
	Complex* array;
}Matrix_C;

//矩阵分配内存
Bool SetMatrixSize(Matrix* matrix, const unsigned row, const unsigned column, const unsigned height);
//设置Matrix矩阵中的所有元素大小为ele
void SetMatrixEle(Matrix* matrix, matrixType ele);
//设置Matrix矩阵中的所有元素大小为0
void SetMatrixZero(Matrix* matrix);
//判断矩阵是否为空，若为空则返回1，否则返回0
Bool IsNullMatrix(const Matrix* matrix);
//销毁矩阵，即释放为矩阵动态分配的内存,并将矩阵的长宽高置0
void DestroyMatrix(Matrix* matrix);
//计算矩阵可容纳元素个数，即return row*column*height
unsigned MatrixCapacity(const Matrix* matrix);
//||matrix||_2  求A矩阵的2范数
matrixType MatrixNorm2(const Matrix* matrix);
//matrixB = matrix(:,:,height)即拷贝三维矩阵的某层，若matrix为二维矩阵，需将height设置为0
Bool CopyMatrix(Matrix* matrixB, Matrix* matrix, unsigned height);
//matrixC = matrixA * matrixB
Bool MatrixMulMatrix(Matrix* matrixC, const Matrix* matrixA, const Matrix* matrixB);
//对vector中所有元素排序，sign= 0 时为升序，其余为降序
void SortVector(Matrix* vector, int sign);
//打印矩阵
void PrintMatrix(const Matrix* matrix);
//将A分解为Q和R
void QR(Matrix* A, Matrix* Q, Matrix* R);
//计算特征值和特征向量
void Eigenvectors(Matrix* eigenVector, Matrix* A, Matrix* eigenValue);

Bool SetMatrixSize_C(Matrix_C* matrix, const unsigned row, const unsigned column, const unsigned height);
Bool IsNullMatrix_C(const Matrix_C* matrix);
void DestroyMatrix_C(Matrix_C* matrix);
Bool CopyMatrix_C(Matrix_C* matrixB, Matrix_C* matrix, unsigned height);
void PrintMatrix_C(const Matrix_C* matrix);
#endif


