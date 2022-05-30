#pragma once

#include <math.h> 
#include "op_on_matr.h"

typedef double (*vFunctionCall)(double arg);
typedef double (*vFunctionCall2)(double arg1, double arg2);

struct Gauss {
	
	
	
	double w[5][5] = {
		
		{ 2 },
		{ 1, 1 },
		{ 0.555556, 0.888889, 0.555556 },
		{ 0.347855, 0.6521445, 0.6521445, 0.347855 },
		{ 0.236927, 0.478629, 0.568889, 0.478629, 0.236927 }
	
	};;

	double x[5][5] = {

		{ 0 },
		{ -0.577350, 0.577350 },
		{ -0.774597, 0, 0.774597 },
		{ -0.861136, -0.339981, 0.339981, 0.861136 },
		{ -0.906180, -0.538469, 0, 0.538469, 0.906180 }

	};;

	double quad1D(vFunctionCall f, int n = 2);
	double quad1D(double a, double b, vFunctionCall f, int n = 2);
	double quad2D(vFunctionCall2 f, int n = 2);
};