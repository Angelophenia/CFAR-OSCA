#pragma once
#ifndef CFAR_OSCA_H
#define CFAR_OSCA_H


typedef struct OSCA_OBJ {
	unsigned ProtectNumber = 0;
	unsigned ReferenceNumber = 20;
	float alpha = 0;
	float Pfa = 1e-6;
	float KO = 0.75;
	unsigned discardCellLR[2] = { 0,0 };
	unsigned discardCellUD[2] = { 0,0 };
};

float form_alpha_OSCA(OSCA_OBJ obj);
OSCA_OBJ form_OSCA_obj();
unsigned* CFAR_OSCA_Range_2D(OSCA_OBJ obj, float* sig, unsigned size1, unsigned size2, unsigned* Target_num);
#endif