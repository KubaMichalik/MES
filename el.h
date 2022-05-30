#pragma once

struct Element {

	double H[4][4] = { 0 };
	double C[4][4] = { 0 };
	double Hbc[4][4] = { 0 };
	double detJ = 0;
	double J[2][2] = { 0 };
	double invJ[2][2] = { 0 };
	double P[4] = { 0 };
	int index = 0;
	int id[4] = { 0, 0, 0, 0 };
};