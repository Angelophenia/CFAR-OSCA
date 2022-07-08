//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include <time.h>
//#include "complex_op.h" 
//
//
////--------------------------������һЩ����Ľṹ�����������---------
////��C���涨��Ĳ�������
//typedef enum { False = 0, True = 1 }Bool;
//
////�������Ԫ�ص�����ΪmatrixType
//typedef complex matrixType;
//
////�˽ṹ��������ʾ��������rowΪ�У�columnΪ�У�heightΪ�ߣ�array������ž���Ԫ��(��һά��ģ���ά/��ά)
//typedef struct
//{
//	unsigned  row, column, height;
//	matrixType* array; //ʹ��ʱ�������*array���г�ʼ��
//}Matrix;
//
////��������ڴ�
//Bool SetMatrixSize(Matrix* matrix, const unsigned row, const unsigned column, const unsigned height);
////����Matrix�����е�����Ԫ�ش�СΪele
//void SetMatrixEle(Matrix* matrix, matrixType ele);
////����Matrix�����е�����Ԫ�ش�СΪ0
//void SetMatrixZero(Matrix* matrix);
////�жϾ����Ƿ�Ϊ�գ���Ϊ���򷵻�1�����򷵻�0
//Bool IsNullMatrix(const Matrix* matrix);
////���پ��󣬼��ͷ�Ϊ����̬������ڴ�,��������ĳ������0
//void DestroyMatrix(Matrix* matrix);
////������������Ԫ�ظ�������return row*column*height
//unsigned MatrixCapacity(const Matrix* matrix);
////||matrix||_2  ��A�����2����
//matrixType MatrixNorm2(const Matrix* matrix);
////matrixB = matrix(:,:,height)��������ά�����ĳ�㣬��matrixΪ��ά�����轫height����Ϊ0
//Bool CopyMatrix(Matrix* matrixB, Matrix* matrix, unsigned height);
////matrixC = matrixA * matrixB
//Bool MatrixMulMatrix(Matrix* matrixC, const Matrix* matrixA, const Matrix* matrixB);
////��vector������Ԫ������sign= 0 ʱΪ��������Ϊ����
//void SortVector(Matrix* vector, int sign);
////��ӡ����
//void PrintMatrix(const Matrix* matrix);
////��A�ֽ�ΪQ��R
//void QR(Matrix* A, Matrix* Q, Matrix* R);
////��������ֵ����������
//void Eigenvectors(Matrix* eigenVector, Matrix* A, Matrix* eigenValue);
//
////---------������QR�ֽ⣬������Է������õ���һЩ����-----------
///*
//matrixΪҪ���ô�С�������ڴ�ľ���row��column��height�ֱ�Ϊ�У��У��ߡ�
//�������óɹ����򷵻�true,���򷵻�false
//*/
//Bool SetMatrixSize(Matrix* matrix, const unsigned row, const unsigned column, const unsigned height)
//{
//	unsigned size = row * column * height * sizeof(matrixType);
//
//	if (size <= 0)
//	{
//		return False;
//	}
//
//	matrix->array = (matrixType*)malloc(size);
//
//	//��������ڴ�ɹ�
//	if (matrix->array)
//	{
//		matrix->row = row;
//		matrix->column = column;
//		matrix->height = height;
//
//		return True;
//	}
//	else
//	{
//		matrix->row = matrix->column = matrix->height = 0;
//
//		return False;
//	}
//}
//
////����Matrix�����е�����Ԫ�ش�СΪele
//void SetMatrixEle(Matrix* matrix, double ele)
//{
//	unsigned size = matrix->row * matrix->column * matrix->height;
//	unsigned i;
//
//	for (i = 0;i < size;++i)
//	{
//		matrix->array[i].real = ele;
//		matrix->array[i].img = 0;
//	}
//}
//
////����Matrix�����е�����Ԫ�ش�СΪ0
//void SetMatrixZero(Matrix* matrix)
//{
//	SetMatrixEle(matrix, 0);
//}
//
////�жϾ����Ƿ�Ϊ�գ���Ϊ���򷵻�1�����򷵻�0
//Bool IsNullMatrix(const Matrix* matrix)
//{
//	unsigned size = matrix->row * matrix->column * matrix->column;
//
//	if (size <= 0 || matrix->array == NULL)
//	{
//		return True;
//	}
//	return False;
//}
//
////���پ��󣬼��ͷ�Ϊ����̬������ڴ�,��������ĳ������0
//void DestroyMatrix(Matrix* matrix)
//{
//	if (!IsNullMatrix(matrix))
//	{
//		free(matrix->array);
//		matrix->array = NULL;
//	}
//	matrix->row = matrix->column = matrix->height = 0;
//}
//
////������������Ԫ�ظ�������return row*column*height
//unsigned MatrixCapacity(const Matrix* matrix)
//{
//	return matrix->row * matrix->column * matrix->height;
//}
//
////||matrix||_2  ��A�����2����
//matrixType MatrixNorm2(const Matrix* matrix)
//{
//	unsigned size = matrix->row * matrix->column * matrix->height;
//	unsigned i;
//
//	matrixType norm = { 0, 0 }, temp;
//
//	for (i = 0;i < size;++i)
//	{
//		temp = ComplexMul(matrix->array[i], matrix->array[i]);
//		//printf("%d A: %.2lf, B: %.2lf, temp: %.2lf\n",i, matrix->array[0].real, matrix->array[0].real, temp.real);
//		norm.real += temp.real;
//		norm.img += temp.img;
//
//	}
//	//printf("norm:%f + %f\n", norm.real, norm.img);
//	return (matrixType)root(norm);
//}
//
////matrixB = matrix(:,:,height)��������ά�����ĳ�㣬��matrixΪ��ά�����轫height����Ϊ0
//Bool CopyMatrix(Matrix* matrixB, Matrix* matrix, unsigned height)
//{
//	unsigned size, i;
//	Matrix matrixA;
//
//	//�ж�heightֵ�Ƿ���ȷ
//	if (height < 0 || height >= matrix->height)
//	{
//		printf("ERROR: CopyMatrix() parameter error��\n");
//		return False;
//	}
//
//	//��matrix(:,:,height) ת��Ϊ��ά����matrixA
//	matrixA.row = matrix->row;
//	matrixA.column = matrix->column;
//	matrixA.height = 1;
//	matrixA.array = matrix->array;
//
//	//�ж�������ָ����ڴ��Ƿ����
//	if (matrixA.array == matrixB->array)
//	{
//		return True;
//	}
//
//	//����matrixA������
//	size = MatrixCapacity(&matrixA);
//
//	//�ж�matrixB��������matrixA�������Ƿ����
//	if (MatrixCapacity(matrixB) != size)
//	{
//		DestroyMatrix(matrixB);
//		SetMatrixSize(matrixB, matrixA.row, matrixA.column, matrixA.height);
//	}
//	else
//	{
//		matrixB->row = matrixA.row;
//		matrixB->column = matrixA.column;
//		matrixB->height = matrixA.height;
//	}
//
//	for (i = 0;i < size;++i)
//	{
//		matrixB->array[i] = matrixA.array[i];
//	}
//
//	return True;
//}
//
////matrixC = matrixA * matrixB
//Bool MatrixMulMatrix(Matrix* matrixC, const Matrix* matrixA, const Matrix* matrixB)
//{
//	size_t row_i, column_i, i;
//	size_t indexA, indexB, indexC;
//
//	matrixType temp, tt;
//	Matrix tempC;
//
//	if (IsNullMatrix(matrixA) || IsNullMatrix(matrixB))
//	{
//		return False;
//	}
//
//	if (matrixA->column != matrixB->row)
//	{
//		return False;
//	}
//
//	if (MatrixCapacity(matrixC) != matrixA->row * matrixB->column)
//	{
//		SetMatrixSize(&tempC, matrixA->row, matrixB->column, 1);
//	}
//	else
//	{
//		tempC.array = matrixC->array;
//		tempC.row = matrixA->row;
//		tempC.column = matrixB->column;
//		tempC.height = 1;
//	}
//
//	for (row_i = 0;row_i < tempC.row;++row_i)
//	{
//		for (column_i = 0;column_i < tempC.column;++column_i)
//		{
//			temp = { 0,0 };
//			for (i = 0;i < matrixA->column;++i)
//			{
//				indexA = row_i * matrixA->column + i;
//				indexB = i * matrixB->column + column_i;
//				tt = ComplexMul(matrixA->array[indexA], matrixB->array[indexB]);
//				temp.real += tt.real;
//				temp.img += tt.img;
//			}
//			indexC = row_i * tempC.column + column_i;
//			tempC.array[indexC] = temp;
//		}
//	}
//
//	if (tempC.array != matrixC->array)
//	{
//		DestroyMatrix(matrixC);
//		matrixC->array = tempC.array;
//	}
//	matrixC->row = tempC.row;
//	matrixC->column = tempC.column;
//	matrixC->height = tempC.height;
//
//	return True;
//}
//
////��vector������Ԫ������sign= 0 ʱΪ��������Ϊ����
//void SortVector(Matrix* vector, int sign)
//{
//	double mid;
//
//	int midIndex;
//	int size = MatrixCapacity(vector);
//	int i, j;
//
//	if (0 == sign)
//	{
//		for (i = 0;i < size;++i)
//		{
//			mid = vector->array[i].real;
//			midIndex = i;
//			for (j = i + 1; j < size; ++j)
//			{
//				if (mid > vector->array[j].real)
//				{
//					mid = vector->array[j].real;
//					midIndex = j;
//				}
//			}
//			vector->array[midIndex] = vector->array[i];
//			vector->array[i].real = mid;
//		}
//	}
//	else
//	{
//		for (i = 0;i < size;++i)
//		{
//			mid = vector->array[i].real;
//			midIndex = i;
//
//			for (j = i + 1; j < size; ++j)
//			{
//				if (mid < vector->array[j].real)
//				{
//					mid = vector->array[j].real;
//					midIndex = j;
//				}
//			}
//
//			vector->array[midIndex] = vector->array[i];
//			vector->array[i].real = mid;
//		}
//	}
//}
//
////��ӡ����
//void PrintMatrix(const Matrix* matrix)
//{
//	size_t row_i, column_i, height_i, index;
//
//	for (height_i = 0;height_i < matrix->height;++height_i)
//	{
//		(matrix->height == 1) ? printf("[:,:] = \n") : printf("[%d,:,:] = \n", height_i);
//
//		for (row_i = 0;row_i < matrix->row;++row_i)
//		{
//			for (column_i = 0;column_i < matrix->column;++column_i)
//			{
//				index = height_i * matrix->row * matrix->column + row_i * matrix->column + column_i;
//				printf("%.2lf + %.2lf, ", matrix->array[index].real, matrix->array[index].img);
//			}
//			printf("\n");
//		}
//	}
//}
//
////----------------------QR�ֽ�-------------------------------------------
////��A�ֽ�ΪQ��R
//void QR(Matrix* A, Matrix* Q, Matrix* R)
//{
//	unsigned  i, j, k, m;
//	unsigned size;
//	const unsigned N = A->row;
//
//	matrixType temp, tt;
//	Matrix a, b;
//
//	//���A����һ����ά��������ʾ���󣬺����������
//	if (A->row != A->column || 1 != A->height)
//	{
//		printf("ERROE: QR() parameter A is not a square matrix!\n");
//		return;
//	}
//
//	size = MatrixCapacity(A);
//	//printf(" A->height: %d !!\n", A->height);
//
//	if (MatrixCapacity(Q) != size)
//	{
//		DestroyMatrix(Q);
//		SetMatrixSize(Q, A->row, A->column, A->height);
//		SetMatrixZero(Q);
//	}
//	else
//	{
//		Q->row = A->row;
//		Q->column = A->column;
//		Q->height = A->height;
//	}
//
//	if (MatrixCapacity(R) != size)
//	{
//		DestroyMatrix(R);
//		SetMatrixSize(R, A->row, A->column, A->height);
//		SetMatrixZero(R);
//	}
//	else
//	{
//		R->row = A->row;
//		R->column = A->column;
//		R->height = A->height;
//	}
//
//	SetMatrixSize(&a, N, 1, 1);
//	SetMatrixSize(&b, N, 1, 1);
//	for (j = 0; j < N;++j)
//	{
//		for (i = 0;i < N; ++i)
//		{
//			a.array[i] = b.array[i] = A->array[i * A->column + j];
//		}
//		for (k = 0; k < j; ++k)
//		{
//			R->array[k * R->column + j] = { 0,0 };
//
//			for (m = 0;m < N; ++m)
//			{
//				tt = ComplexMul(a.array[m], Q->array[m * Q->column + k]);
//				printf("%d %d %d A: %.2lf+%.2lfi, B: %.2lf+%.2lfi, temp: %.2lf+%.2lfi\n", j, k, m, a.array[m].real, a.array[m].img, Q->array[m * Q->column + k].real, Q->array[m * Q->column + k].img, tt.real, tt.img);
//
//				R->array[k * R->column + j].real += tt.real;
//				R->array[k * R->column + j].img += tt.img;
//			}
//
//			for (m = 0;m < N; ++m)
//			{
//				tt = ComplexMul(R->array[k * R->column + j], Q->array[m * Q->column + k]);
//				//printf("����%d %d %d A: %.2lf, B: %.2lf, temp: %.2lf\n", j, k, m, a.array[m].real, Q->array[m].real, tt.real);
//
//				b.array[m].real -= tt.real;
//				b.array[m].img -= tt.img;
//
//			}
//		}
//
//		temp = MatrixNorm2(&b);
//
//		R->array[j * R->column + j] = temp;
//
//		for (i = 0; i < N; ++i)
//		{
//			tt = ComplexMul(b.array[i], temp);
//			Q->array[i * Q->column + j] = tt;
//		}
//	}
//
//	DestroyMatrix(&a);
//	DestroyMatrix(&b);
//}
//
////----------------------ʹ������ֵ���������������-----------------
////eigenVectorΪ������������A����������
////eigenValueΪ����A����������ֵ��
////AΪҪ�������������ľ���
//void Eigenvectors(Matrix* eigenVector, Matrix* A, Matrix* eigenValue)
//{
//	unsigned i, j, q;
//	unsigned count;
//	int m;
//	const unsigned NUM = A->column;
//
//	matrixType eValue;
//	matrixType sum, midSum, mid, tt;
//	Matrix temp;
//
//	SetMatrixSize(&temp, A->row, A->column, A->height);
//
//	for (count = 0; count < NUM;++count)
//	{
//		//��������ֵΪeValue�������������ʱ��ϵ������
//		eValue = eigenValue->array[count];
//		CopyMatrix(&temp, A, 0);
//		for (i = 0; i < temp.column; ++i)
//		{
//			temp.array[i * temp.column + i].real -= eValue.real;
//			temp.array[i * temp.column + i].img -= eValue.img;
//
//		}
//
//		//��temp��Ϊ�����;���
//		for (i = 0; i < temp.row - 1; ++i)
//		{
//			mid = temp.array[i * temp.column + i];
//			for (j = i; j < temp.column; ++j)
//			{
//				tt = ComplexMul(temp.array[i * temp.column + j], mid);
//				temp.array[i * temp.column + j] = tt;
//			}
//
//			for (j = i + 1;j < temp.row;++j)
//			{
//				mid = temp.array[j * temp.column + i];
//				for (q = i; q < temp.column; ++q)
//				{
//					tt = ComplexMul(mid, temp.array[i * temp.column + q]);
//					temp.array[j * temp.column + q].real -= tt.real;
//					temp.array[j * temp.column + q].img -= tt.img;
//				}
//			}
//		}
//		midSum = eigenVector->array[(eigenVector->row - 1) * eigenVector->column + count] = { 1, 0 };
//		for (m = temp.row - 2; m >= 0; --m)
//		{
//			sum = { 0,0 };
//			for (j = m + 1;j < temp.column; ++j)
//			{
//				tt = ComplexMul(temp.array[m * temp.column + j], eigenVector->array[j * eigenVector->column + count]);
//				sum.real += tt.real;
//				sum.img += tt.img;
//			}
//
//			if (temp.array[m * temp.column + m].real == 0 && temp.array[m * temp.column + m].img == 0) {
//				printf("000000~!!!!!\n");
//			}
//			tt = ComplexDiv(sum, temp.array[m * temp.column + m]);
//			sum.real = -tt.real;
//			sum.img = -tt.img;
//			tt = ComplexMul(sum, sum);
//			midSum.real += tt.real;
//			midSum.img += tt.img;
//
//			eigenVector->array[m * eigenVector->column + count] = sum;
//		}
//		printf("midSum:%f + %f\n", midSum.real, midSum.img);
//		midSum = root(midSum);
//
//		for (i = 0; i < eigenVector->row; ++i)
//		{
//			if (midSum.real == 0 && midSum.img == 0) {
//				printf("000000~!!!!!\n");
//			}
//			tt = ComplexDiv(eigenVector->array[i * eigenVector->column + count], midSum);
//			eigenVector->array[i * eigenVector->column + count] = tt;
//		}
//	}
//	DestroyMatrix(&temp);
//}
//
//
//int main()
//{
//	const unsigned NUM = 50; //����������
//
//	unsigned N = 3;
//	unsigned k;
//
//	Matrix A, Q, R, temp;
//	Matrix eValue;
//
//	//�����ڴ�
//	SetMatrixSize(&A, N, N, 1);
//	SetMatrixSize(&Q, A.row, A.column, A.height);
//	SetMatrixSize(&R, A.row, A.column, A.height);
//	SetMatrixSize(&temp, A.row, A.column, A.height);
//	SetMatrixSize(&eValue, A.row, 1, 1);
//
//	//A����Ϊһ���򵥾���
//	complex B[] = {
//		{3.258, 0},{0.1214,-0.597},{0.8564,-0.757},
//		{0.1214,0.597},{3.215,0},{0.3058,-0.5341},
//		{0.8564,0.757},{0.3058,0.5341},{3.3409,0}
//	};
//	A.array = B;
//	//A.array[0] = 5.8751;
//	//A.array[1] = 6.0764;
//	//A.array[2] = 3.5669;
//	//A.array[3] = 2.7109;
//	//A.array[4] = 6.0764;
//	//A.array[5] = 6.8258;
//	//A.array[6] = 3.8662;
//	//A.array[7] = 3.0022;
//	//A.array[8] = 3.5669;
//	//A.array[9] = 3.8662;
//	//A.array[10] = 3.9534;
//	//A.array[11] = 2.9391;
//	//A.array[12] = 2.7109;
//	//A.array[13] = 3.0022;
//	//A.array[14] = 2.9391;
//	//A.array[15] = 2.2816;
//
//
//	PrintMatrix(&A);
//
//	//����A����Ԫ����temp
//	CopyMatrix(&temp, &A, 0);
//	//��ʼ��Q��R����Ԫ��Ϊ0
//	SetMatrixZero(&Q);
//	SetMatrixZero(&R);
//	//ʹ��QR�ֽ����������ֵ
//	for (k = 0;k < NUM; ++k)
//	{
//		QR(&temp, &Q, &R);
//		MatrixMulMatrix(&temp, &R, &Q);
//	}
//
//	//��ȡ����ֵ����֮�洢��eValue
//	for (k = 0;k < temp.column;++k)
//	{
//		eValue.array[k] = temp.array[k * temp.column + k];
//	}
//
//	//������ֵ���ս�������
//	SortVector(&eValue, 1);
//
//	//��������ֵeValue��ԭʼ������������������Q
//	Eigenvectors(&Q, &A, &eValue);
//
//	//��ӡ����ֵ
//	printf("����ֵ��");
//	PrintMatrix(&eValue);
//
//	//��ӡ��������
//	printf("��������");
//	PrintMatrix(&Q);
//
//	//DestroyMatrix(&A);
//	DestroyMatrix(&R);
//	DestroyMatrix(&Q);
//	DestroyMatrix(&eValue);
//	DestroyMatrix(&temp);
//
//	return 0;
//}