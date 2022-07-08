#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "complex_op.h"
#define PI 3.1415926
#define Stride 0.005


float randn() {
	float u = ((float)rand() / (RAND_MAX)) * 2 - 1;
	float v = ((float)rand() / (RAND_MAX)) * 2 - 1;
	float r = u * u + v * v;
	if (r == 0 || r > 1) return randn();
	float c = sqrt(-2 * log(r) / r);
	return u * c;
}

void free_n(void *x) {
	free(x);
	x = NULL;
}

//Complex* toeplitz(Complex* X, unsigned num) {
//	Complex* X1 = (Complex*)malloc(sizeof(Complex) * num * num);
//	if (!X1) {
//		printf("malloc error");
//		return 0;
//	}
//	for (int i = 0; i < num; i++) {
//		for (int j = 0; j < (num-i); j++) {
//			X1[j + i + j * num] = X[i];
//			X1[j + (j + i) * num].real = X[i].real; //上下三角共轭
//			X1[j + (j + i) * num].img = -X[i].img;
//		}
//	}
//	return X1;
//}

void init_array(double* x, unsigned num) {
	for (unsigned i = 0; i < num; i++) {
		x[i] = 0.0;
	}
}

void quickSort(float* array, int left, int right)
{
	if (NULL == array)
	{
		return;
	}

	if (left < right)
	{
		float pivot = array[left];
		int low = left, high = right;
		while (low < high)
		{
			while (array[high] >= pivot && low < high)
				high--;
			array[low] = array[high];

			while (array[low] <= pivot && low < high)
				low++;
			array[high] = array[low];
		}
		array[low] = pivot;

		quickSort(array, left, low - 1);
		quickSort(array, low + 1, right);
	}
}


unsigned find(unsigned *x, unsigned size1, unsigned size2, unsigned target, unsigned *y) {
	unsigned num = 0, capacity=3;
	for (unsigned j = 0; j < size2; j++) {
		for (unsigned i = 0; i < size1; i++) {
			if (x[i * size2 + j] == target) {
				//if (capacity < num) {  //扩容
				//	capacity *= 2;
				//	//y = (unsigned*)realloc(y, capacity * 2 * sizeof(unsigned));

				//	tmp = (unsigned*)realloc(y, capacity * 2 * sizeof(unsigned));
				//	if (tmp != NULL)
				//	{
				//		y = tmp;
				//		tmp = NULL;
				//	}
				//	else
				//	{
				//		printf("reallocmemory failed\n");
				//	}
				//	printf("realloc %d\n", capacity);
				//}
				y[num * 2] = i;
				y[num * 2 + 1] = j;
				num++;
				//printf("%d: %d %d\n", num, i, j);
			}
		}
	}
	return num;
}


unsigned* local_max(double *arr, unsigned num, unsigned* out_count) {
	if (num < 3) {
		printf("too small!\n");
		return 0;
	}
	unsigned* out_ind = (unsigned*)malloc(sizeof(unsigned)*num);
	if (!out_ind) {
		printf("malloc failed\n");
	}
	unsigned count = 0;
	for (unsigned i = 1; i < num-1;i++) {
		if ((arr[i] - arr[i - 1])>3 && (arr[i] - arr[i + 1])>3) {
			out_ind[count] = i;
			count++;
			i++;  //next one more
		}
	}
	*out_count = count;
	return out_ind;
}




