#pragma once

#include <math.h>

#include "gauss.h"
#include "el.h"

struct Element4_2D {
	double nMrx[4][4], KsiMrx[4][4], EtaMrx[4][4];
	double dNdX[4][4], dNdY[4][4];
	double dKsi[4], dEta[4];
	double *w;
	double ksi[4], eta[4];
	double detJ = 0;
	int p = 4;
	Element4_2D();
	Element4_2D(int integralpts = 2);
	
	
	
	double* nKsiEta(double ksi, double eta);
	double* dNdEta(double ksi);
	double* dNdKsi(double eta);
	double dXYdEta(int pc, double xy);
	double dXYdKsi(int i, double xy);
	
	
};