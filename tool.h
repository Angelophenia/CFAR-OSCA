#pragma once
#ifndef TOOL_H
#define TOOL_H

float randn();
void free_n(void* x);
//Complex* toeplitz(Complex* X, unsigned num);
void init_array(double* x, unsigned num);
void quickSort(float* array, int left, int right);
unsigned find(unsigned* x, unsigned size1, unsigned size2, unsigned target, unsigned* y);
unsigned* local_max(double* arr, unsigned num, unsigned* out_count);
#endif
