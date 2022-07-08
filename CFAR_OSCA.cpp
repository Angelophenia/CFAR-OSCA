#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "tool.h"
#include "complex_op.h"
#include  "CFAR_OSCA.h"
#include "main.h"


float form_alpha_OSCA(OSCA_OBJ obj) {
	unsigned N = obj.ReferenceNumber;
	float K = ceil(obj.KO * double(N));  //向上取整
	float sum = 0;
	for (int i = N-K+1; i < N+1; i++) {
		sum += 1.0/i;
	}
	float alpha = -log(obj.Pfa) / sum;
	return alpha;
}


OSCA_OBJ form_OSCA_obj() {
	OSCA_OBJ obj;
	obj.alpha = form_alpha_OSCA(obj);
	return obj;
}


unsigned* CFAR_OSCA_Range_2D(OSCA_OBJ obj, float *sig, unsigned size1, unsigned size2, unsigned * Target_num) {
	unsigned cellNum = obj.ReferenceNumber / 2, gapNum = obj.ProtectNumber / 2;
	unsigned gaptot = gapNum + cellNum;
	float alpha = obj.alpha, K0=0.75;
	unsigned N_obj=0;
	unsigned discardCellLeft, discardCellRight, discardCellTop, discardCellBottom;
	discardCellLeft = obj.discardCellLR[0];
	discardCellRight = obj.discardCellLR[1];
	discardCellTop = obj.discardCellUD[0];
	discardCellBottom = obj.discardCellUD[1];
	unsigned size1_1 = (size1 - discardCellBottom - discardCellTop);
	unsigned size2_1 = (size2 - discardCellLeft - discardCellRight);
	unsigned size3 = (size2_1 + 2 * gaptot);
	unsigned size4 = (size1_1 + 2 * gaptot);
	float* sigv = (float*)malloc(sizeof(float) * (size2+size1));
	float* v = (float*)malloc(sizeof(float) * (size2_1+size1_1));
	float* vec = (float*)malloc(sizeof(float) * size1 * size3);
	float* vector = (float*)malloc(sizeof(float) * size4 * size3);
	float* ordered_vector = (float*)malloc(sizeof(float) * (2 * gaptot + 1) * (2 * gaptot + 1));
	float* temp_v = (float*)malloc(sizeof(float) * (2 * gaptot+1));
	float* temp_v1 = (float*)malloc(sizeof(float) * (2 * gaptot + 1));
	float* Clutter = (float*)malloc(sizeof(float) * size1_1 * size2_1);
	unsigned* result = (unsigned*)malloc(sizeof(unsigned) * size1_1 * size2_1);
	unsigned* TargetInd = NULL, * Ind_obj = NULL, capacity = 310;
	TargetInd = (unsigned*)malloc(sizeof(unsigned) * 2 * capacity); //目标的RD图坐标
	if (!sigv || !v || !vec || !vector || !ordered_vector || !temp_v || !Clutter || !result || !TargetInd) {
		printf("fail malloc");

	}
	for (int k = 0; k < size1; k++) {
		get_row_f(sig, sigv, size2, k, 0, size2);
		for (int i = discardCellLeft, j=0; i < size2 - discardCellRight; i++, j++) { //去头掐尾
			v[j] = sigv[i];
			//print("%e ", sigv[i]);
		}
		//print("\n");
		/* vec(k,:) = [vecLeft v vecRight] */
		for (int q = 0; q <  gaptot; q++) {
			vec[k * size3 + q] = v[q];
			vec[k * size3 + gaptot + size2_1 + q] = v[size2_1-gaptot + q];
		}
		for (int p = gaptot; p < gaptot + size2_1; p++) {
			vec[k * size3 + p] = v[p - gaptot];
		}
		
	}

	for (int k = 0; k < size3; k++) {
		get_column_f(vec, sigv, size3, k, 0, size1);
		for (int i = discardCellTop, j = 0; i < size1 - discardCellBottom; i++, j++) { //去头掐尾
			v[j] = sigv[i];
		}
		/* vector(:,k) = [vecLeft v vecRight]' */
		for (int q = 0; q < gaptot; q++) {
			vector[q * size3 + k] = v[q];
			vector[(gaptot + size1_1 + q) * size3 + k] = v[size1_1-gaptot + q];
		}
		for (int p = gaptot; p < gaptot + size1_1; p++) {
			vector[p * size3 + k] = v[p - gaptot];
		}
		
	}
	for (int nd = 0; nd < size2_1;nd++) {
		
		for (int nr = 0; nr < size1_1;nr++) {
			
			for (int c = nd, num=0; c < nd + 2 * gaptot+1; c++, num++) {  //取nr 至 nr + 2 * gaptot列 （共2 * gaptot+1列）
				get_column_f(vector, temp_v, size3, c, nr, nr + 2 * gaptot+1);  //截取每列的 nr 至 nr + 2 * gaptot行 （共2 * gaptot+1行）
				quickSort(temp_v, 0, 2 * gaptot); // 列排序，小到大 //dui????
				for (int j=0; j < 2 * gaptot + 1; j++) {  //一列
					ordered_vector[j * (2 * gaptot + 1) + num] = temp_v[j];
				} 
			}
			/* K0_vector = ordered_vector(ceil(K0*(cellNum*2+1)),:) */
			get_row_f(ordered_vector,  temp_v, 2 * gaptot + 1, unsigned(ceil(K0 * (cellNum * 2.0))), 0, 21);
			Clutter[nr * size2_1 + nd] = mean(temp_v, 2 * gaptot + 1);
			//print("\n——%d %d:\n", nr + 1, nd + 1);
			//print("%f \n", Clutter[nr * size2_1 + nd]);	
		}
	}

	for (int i = gaptot, ii=0; i < gaptot + size1_1; i++, ii++) {
		for (int j = gaptot, jj=0; j < gaptot + size2_1; j++, jj++) {
			result[ii * size2_1 + jj] = (alpha * Clutter[ii * size2_1 + jj]) < vector[i * size3 + j] ? 1 : 0;
			/*if (result[ii * size2_1 + jj]==1)
				print("%d %d %d", result[ii * size2_1 + jj],ii,jj);*/
		}
		//print("\n");
	}
	
	N_obj = find(result, size1_1, size2_1, 1, TargetInd);
	*Target_num = N_obj;
	if (N_obj == 0) {
		printf("Nothing!!\n");
	}
	else {
		Ind_obj = (unsigned*)malloc(sizeof(unsigned) * N_obj * 2); //坐标
		float *noise_obj = (float*)malloc(sizeof(float) * N_obj); //噪声估计
		float* CFAR_SNR = (float*)malloc(sizeof(float) * N_obj); //信噪比
		if (!Ind_obj || !noise_obj || !CFAR_SNR) {
			printf("malloc failed!\n");
		}
		else {
			//printf("CFAR_SNR\n");
			for (int i = 0; i < N_obj; i++) {
				Ind_obj[i * 2] = TargetInd[i * 2] + discardCellTop;
				Ind_obj[i * 2 + 1] = TargetInd[i * 2 + 1] + discardCellLeft;
				//print("%d %d", TargetInd[i * 2], TargetInd[i * 2 + 1]);
				//print("%d %d", Ind_obj[i * 2], Ind_obj[i * 2+1]);
				noise_obj[i] = Clutter[TargetInd[i * 2+1] + TargetInd[i * 2] * size2_1];
				CFAR_SNR[i] = sig[Ind_obj[i * 2+1] + Ind_obj[i * 2] * size2]/ noise_obj[i];
				//printf("%.3lf\n", noise_obj[i]);
				//printf("%.3lf\n", CFAR_SNR[i]);
			}
			free_n(noise_obj);
			//free_n(Ind_obj);
			free_n(CFAR_SNR);
		}
	}
	free_n(TargetInd);
	free_n(sigv);
	free_n(vec);
	free_n(v);
	free_n(vector);
	free_n(ordered_vector);
	free_n(temp_v);
	free_n(Clutter);
	free_n(result);
	return Ind_obj;
}


 