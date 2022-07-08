#include <stdio.h>
#include <stdlib.h>
#include "eig.h"


//---------������QR�ֽ⣬������Է������õ���һЩ����-----------
/*
matrixΪҪ���ô�С�������ڴ�ľ���row��column��height�ֱ�Ϊ�У��У��ߡ�
�������óɹ����򷵻�true,���򷵻�false
*/
Bool SetMatrixSize(Matrix* matrix, const unsigned row, const unsigned column, const unsigned height)
{
	unsigned size = row * column * height * sizeof(matrixType);

	if (size <= 0)
	{
		return False;
	}

	matrix->array = (matrixType*)malloc(size);

	//��������ڴ�ɹ�
	if (matrix->array)
	{
		matrix->row = row;
		matrix->column = column;
		matrix->height = height;

		return True;
	}
	else
	{
		matrix->row = matrix->column = matrix->height = 0;

		return False;
	}
}
Bool SetMatrixSize_C(Matrix_C* matrix, const unsigned row, const unsigned column, const unsigned height)
{
	unsigned size = row * column * height * sizeof(Complex);

	if (size <= 0)
	{
		return False;
	}

	matrix->array = (Complex*)malloc(size);

	//��������ڴ�ɹ�
	if (matrix->array)
	{
		matrix->row = row;
		matrix->column = column;
		matrix->height = height;

		return True;
	}
	else
	{
		matrix->row = matrix->column = matrix->height = 0;

		return False;
	}
}

//����Matrix�����е�����Ԫ�ش�СΪele
void SetMatrixEle(Matrix* matrix, matrixType ele)
{
	unsigned size = matrix->row * matrix->column * matrix->height;
	unsigned i;

	for (i = 0; i < size; ++i)
	{
		matrix->array[i] = ele;
	}
}

//����Matrix�����е�����Ԫ�ش�СΪ0
void SetMatrixZero(Matrix* matrix)
{
	SetMatrixEle(matrix, 0);
}

//�жϾ����Ƿ�Ϊ�գ���Ϊ���򷵻�1�����򷵻�0
Bool IsNullMatrix(const Matrix* matrix)
{
	unsigned size = matrix->row * matrix->column * matrix->column;

	if (size <= 0 || matrix->array == NULL)
	{
		return True;
	}
	return False;
}
Bool IsNullMatrix_C(const Matrix_C* matrix)
{
	unsigned size = matrix->row * matrix->column * matrix->column;

	if (size <= 0 || matrix->array == NULL)
	{
		return True;
	}
	return False;
}

//���پ��󣬼��ͷ�Ϊ����̬������ڴ�,��������ĳ������0
void DestroyMatrix(Matrix* matrix)
{
	if (!IsNullMatrix(matrix))
	{
		free(matrix->array);
		matrix->array = NULL;
	}
	matrix->row = matrix->column = matrix->height = 0;
}
void DestroyMatrix_C(Matrix_C* matrix)
{
	if (!IsNullMatrix_C(matrix))
	{
		free(matrix->array);
		matrix->array = NULL;
	}
	matrix->row = matrix->column = matrix->height = 0;
}

//������������Ԫ�ظ�������return row*column*height
unsigned MatrixCapacity(const Matrix* matrix)
{
	return matrix->row * matrix->column * matrix->height;
}
unsigned MatrixCapacity_C(const Matrix_C* matrix)
{
	return matrix->row * matrix->column * matrix->height;
}

//||matrix||_2  ��A�����2����
matrixType MatrixNorm2(const Matrix* matrix)
{
	unsigned size = matrix->row * matrix->column * matrix->height;
	unsigned i;

	matrixType norm = 0;

	for (i = 0; i < size; ++i)
	{
		norm += (matrix->array[i]) * (matrix->array[i]);
	}

	return (matrixType)sqrt(norm);
}

//matrixB = matrix(:,:,height)��������ά�����ĳ�㣬��matrixΪ��ά�����轫height����Ϊ0
Bool CopyMatrix(Matrix* matrixB, Matrix* matrix, unsigned height)
{
	unsigned size, i;
	Matrix matrixA;

	//�ж�heightֵ�Ƿ���ȷ
	if (height < 0 || height >= matrix->height)
	{
		printf("ERROR: CopyMatrix() parameter error��\n");
		return False;
	}

	//��matrix(:,:,height) ת��Ϊ��ά����matrixA
	matrixA.row = matrix->row;
	matrixA.column = matrix->column;
	matrixA.height = 1;
	matrixA.array = matrix->array + height * matrix->row * matrix->column;

	//�ж�������ָ����ڴ��Ƿ����
	if (matrixA.array == matrixB->array)
	{
		return True;
	}

	//����matrixA������
	size = MatrixCapacity(&matrixA);

	//�ж�matrixB��������matrixA�������Ƿ����
	if (MatrixCapacity(matrixB) != size)
	{
		DestroyMatrix(matrixB);
		SetMatrixSize(matrixB, matrixA.row, matrixA.column, matrixA.height);
	}
	else
	{
		matrixB->row = matrixA.row;
		matrixB->column = matrixA.column;
		matrixB->height = matrixA.height;
	}

	for (i = 0; i < size; ++i)
	{
		matrixB->array[i] = matrixA.array[i];
	}

	return True;
}

Bool CopyMatrix_C(Matrix_C* matrixB, Matrix_C* matrix, unsigned height)
{
	unsigned size, i;
	Matrix_C matrixA;

	//�ж�heightֵ�Ƿ���ȷ
	if (height < 0 || height >= matrix->height)
	{
		printf("ERROR: CopyMatrix() parameter error��\n");
		return False;
	}

	//��matrix(:,:,height) ת��Ϊ��ά����matrixA
	matrixA.row = matrix->row;
	matrixA.column = matrix->column;
	matrixA.height = 1;
	matrixA.array = matrix->array + height * matrix->row * matrix->column;

	//�ж�������ָ����ڴ��Ƿ����
	if (matrixA.array == matrixB->array)
	{
		return True;
	}

	//����matrixA������
	size = MatrixCapacity_C(&matrixA);

	//�ж�matrixB��������matrixA�������Ƿ����
	if (MatrixCapacity_C(matrixB) != size)
	{
		DestroyMatrix_C(matrixB);
		SetMatrixSize_C(matrixB, matrixA.row, matrixA.column, matrixA.height);
	}
	else
	{
		matrixB->row = matrixA.row;
		matrixB->column = matrixA.column;
		matrixB->height = matrixA.height;
	}

	for (i = 0; i < size; ++i)
	{
		matrixB->array[i] = matrixA.array[i];
	}

	return True;
}

//matrixC = matrixA * matrixB
Bool MatrixMulMatrix(Matrix* matrixC, const Matrix* matrixA, const Matrix* matrixB)
{
	size_t row_i, column_i, i;
	size_t indexA, indexB, indexC;

	matrixType temp;
	Matrix tempC;

	if (IsNullMatrix(matrixA) || IsNullMatrix(matrixB))
	{
		return False;
	}

	if (matrixA->column != matrixB->row)
	{
		return False;
	}

	if (MatrixCapacity(matrixC) != matrixA->row * matrixB->column)
	{
		SetMatrixSize(&tempC, matrixA->row, matrixB->column, 1);
	}
	else
	{
		tempC.array = matrixC->array;
		tempC.row = matrixA->row;
		tempC.column = matrixB->column;
		tempC.height = 1;
	}

	for (row_i = 0; row_i < tempC.row; ++row_i)
	{
		for (column_i = 0; column_i < tempC.column; ++column_i)
		{
			temp = 0;
			for (i = 0; i < matrixA->column; ++i)
			{
				indexA = row_i * matrixA->column + i;
				indexB = i * matrixB->column + column_i;

				temp += matrixA->array[indexA] * matrixB->array[indexB];
			}
			indexC = row_i * tempC.column + column_i;
			tempC.array[indexC] = temp;
		}
	}

	if (tempC.array != matrixC->array)
	{
		DestroyMatrix(matrixC);
		matrixC->array = tempC.array;
	}
	matrixC->row = tempC.row;
	matrixC->column = tempC.column;
	matrixC->height = tempC.height;

	return True;
}

//��vector������Ԫ������sign= 0 ʱΪ��������Ϊ����
void SortVector(Matrix* vector, int sign)
{
	matrixType mid;

	int midIndex;
	int size = MatrixCapacity(vector);
	int i, j;

	if (0 == sign)
	{
		for (i = 0; i < size; ++i)
		{
			mid = vector->array[i];
			midIndex = i;
			for (j = i + 1; j < size; ++j)
			{
				if (mid > vector->array[j])
				{
					mid = vector->array[j];
					midIndex = j;
				}
			}
			vector->array[midIndex] = vector->array[i];
			vector->array[i] = mid;
		}
	}
	else
	{
		for (i = 0; i < size; ++i)
		{
			mid = vector->array[i];
			midIndex = i;

			for (j = i + 1; j < size; ++j)
			{
				if (mid < vector->array[j])
				{
					mid = vector->array[j];
					midIndex = j;
				}
			}

			vector->array[midIndex] = vector->array[i];
			vector->array[i] = mid;
		}
	}
}

//��ӡ����
void PrintMatrix(const Matrix* matrix)
{
	size_t row_i, column_i, height_i, index;

	for (height_i = 0; height_i < matrix->height; ++height_i)
	{
		(matrix->height == 1) ? printf("[:,:] = \n") : printf("[%d,:,:] = \n", height_i);
		//������
		for (row_i = 0; row_i < matrix->row; ++row_i)
		{
			for (column_i = 0; column_i < matrix->column; ++column_i)
			{
				index = height_i * matrix->row * matrix->column + row_i * matrix->column + column_i;
				printf("%.3lf  ", matrix->array[index]);
			}
			printf("\n");
		}
		//������
		//for (column_i = 0; column_i < matrix->column; ++column_i)
		//{
		//	for (row_i = 0; row_i < matrix->row; ++row_i)
		//	{
		//		index = height_i * matrix->row * matrix->column + row_i * matrix->column + column_i;
		//		printf("%.3lf  ", matrix->array[index]);
		//	}
		//	printf("\n");
		//}
	}
}

void PrintMatrix_C(const Matrix_C* matrix)
{
	size_t row_i, column_i, height_i, index;

	for (height_i = 0; height_i < matrix->height; ++height_i)
	{
		(matrix->height == 1) ? printf("[:,:] = \n") : printf("[%d,:,:] = \n", height_i);
		//������
		for (row_i = 0; row_i < matrix->row; ++row_i)
		{
			for (column_i = 0; column_i < matrix->column; ++column_i)
			{
				index = height_i * matrix->row * matrix->column + row_i * matrix->column + column_i;
				printf("%.4lf + %.4lfi ", matrix->array[index].real, matrix->array[index].img);
			}
			printf("\n");
		}
		//������
		//for (column_i = 0; column_i < matrix->column; ++column_i)
		//{
		//	for (row_i = 0; row_i < matrix->row; ++row_i)
		//	{
		//		index = height_i * matrix->row * matrix->column + row_i * matrix->column + column_i;
		//		printf("%.3lf + %.3lfi ", matrix->array[index].real, matrix->array[index].img);
		//	}
		//	printf("\n");
		//}
	}
}

//----------------------QR�ֽ�-------------------------------------------
//��A�ֽ�ΪQ��R
void QR(Matrix* A, Matrix* Q, Matrix* R)
{
	unsigned  i, j, k, m;
	unsigned size;
	const unsigned N = A->row;

	matrixType temp;
	Matrix a, b;

	//���A����һ����ά��������ʾ���󣬺����������
	if (A->row != A->column || 1 != A->height)
	{
		printf("ERROE: QR() parameter A is not a square matrix!\n");
		return;
	}

	size = MatrixCapacity(A);

	if (MatrixCapacity(Q) != size)
	{
		DestroyMatrix(Q);
		SetMatrixSize(Q, A->row, A->column, A->height);
		SetMatrixZero(Q);
	}
	else
	{
		Q->row = A->row;
		Q->column = A->column;
		Q->height = A->height;
	}

	if (MatrixCapacity(R) != size)
	{
		DestroyMatrix(R);
		SetMatrixSize(R, A->row, A->column, A->height);
		SetMatrixZero(R);
	}
	else
	{
		R->row = A->row;
		R->column = A->column;
		R->height = A->height;
	}

	SetMatrixSize(&a, N, 1, 1);
	SetMatrixSize(&b, N, 1, 1);

	for (j = 0; j < N; ++j)
	{
		for (i = 0; i < N; ++i)
		{
			a.array[i] = b.array[i] = A->array[i * A->column + j];
		}
		for (k = 0; k < j; ++k)
		{
			R->array[k * R->column + j] = 0;

			for (m = 0; m < N; ++m)
			{
				R->array[k * R->column + j] += a.array[m] * Q->array[m * Q->column + k];
			}

			for (m = 0; m < N; ++m)
			{
				b.array[m] -= R->array[k * R->column + j] * Q->array[m * Q->column + k];
			}
		}

		temp = MatrixNorm2(&b);

		R->array[j * R->column + j] = temp;

		for (i = 0; i < N; ++i)
		{
			Q->array[i * Q->column + j] = b.array[i] / temp;
		}
	}

	DestroyMatrix(&a);
	DestroyMatrix(&b);
}

//----------------------ʹ������ֵ���������������-----------------
//eigenVectorΪ������������A����������
//eigenValueΪ����A����������ֵ��
//AΪҪ�������������ľ���
void Eigenvectors(Matrix* eigenVector, Matrix* A, Matrix* eigenValue)
{
	unsigned i, j, q;
	unsigned count;
	int m;
	const unsigned NUM = A->column;

	matrixType eValue;
	matrixType sum, midSum, mid;
	Matrix temp;

	SetMatrixSize(&temp, A->row, A->column, A->height);

	for (count = 0; count < NUM; ++count)
	{
		//if (count%2!=0) {	
		//	printf("%d repeat!\n", count);
		//	continue;//�ظ�������ִ��
		//}
		
		//��������ֵΪeValue�������������ʱ��ϵ������
		eValue = eigenValue->array[count];
		CopyMatrix(&temp, A, 0);
		for (i = 0; i < temp.column; ++i)
		{
			temp.array[i * temp.column + i] -= eValue;
		}

		//��temp��Ϊ�����;���
		for (i = 0; i < temp.row - 1; ++i)
		{
			mid = temp.array[i * temp.column + i];
			for (j = i; j < temp.column; ++j)
			{
				temp.array[i * temp.column + j] /= mid;
			}

			for (j = i + 1; j < temp.row; ++j)
			{
				mid = temp.array[j * temp.column + i];
				for (q = i; q < temp.column; ++q)
				{
					temp.array[j * temp.column + q] -= mid * temp.array[i * temp.column + q];
				}
			}
		}
		midSum = eigenVector->array[(eigenVector->row - 1) * eigenVector->column + count] = 1;
		for (m = temp.row - 2; m >= 0; --m)
		{
			sum = 0;
			for (j = m + 1; j < temp.column; ++j)
			{
				sum += temp.array[m * temp.column + j] * eigenVector->array[j * eigenVector->column + count];
			}
			sum = -sum / temp.array[m * temp.column + m];
			midSum += sum * sum;

			eigenVector->array[m * eigenVector->column + count] = sum;
		}
		midSum = sqrt(midSum);

		for (i = 0; i < eigenVector->row; ++i)
		{
			eigenVector->array[i * eigenVector->column + count] /= midSum;
		}
	}
	DestroyMatrix(&temp);
}
