#include "el4_2d.h"

Element4_2D::Element4_2D(int integralpts) {
	this->p = integralpts * integralpts;
	

	ksi[0] = -0.5774;
	eta[0] = -0.5774;
	ksi[1] = 0.5774;
	eta[1] = -0.5774;
	ksi[2] = 0.5774;
	eta[2] = 0.5774;
	ksi[3] = -0.5774;
	eta[3] = 0.5774;

	
	

	for (int i = 0; i < p; i++)
		for (int j = 0; j < 4; j++) {
			
			nMrx[i][j] = nKsiEta(ksi[i], eta[i])[j];
			EtaMrx[i][j] = dNdEta(ksi[i])[j];
			KsiMrx[i][j] = dNdKsi(eta[i])[j];
		
		}
}

// Zwraca tablicê funkcji kszta³tu.
double* Element4_2D::nKsiEta(double ksi, double eta) {
	double* N = new double[4];
	N[0] = 0.25 * (1 - ksi) * (1 - eta);
	N[1] = 0.25 * (1 + ksi) * (1 - eta);
	N[2] = 0.25 * (1 + ksi) * (1 + eta);
	N[3] = 0.25 * (1 - ksi) * (1 + eta);
	return N;
}

// Zwraca tablicê pochodnych funkcji kszta³tu po eta.
double* Element4_2D::dNdEta(double ksi) {
	dEta[0] = -0.25 * (1 - ksi);
	dEta[1] = -0.25 * (1 + ksi);
	dEta[2] = 0.25 * (1 + ksi);
	dEta[3] = 0.25 * (1 - ksi);
	return dEta;
}

// Zwraca tablicê pochodnych funkcji kszta³tu po ksi.
double* Element4_2D::dNdKsi(double eta) {
	dKsi[0] = -0.25 * (1 - eta);
	dKsi[1] = 0.25 * (1 - eta);
	dKsi[2] = 0.25 * (1 + eta);
	dKsi[3] = -0.25 * (1 + eta);

	return dKsi;
}

// Zwraca Jakobian przekszta³cenia po eta.
double Element4_2D::dXYdEta(int pc, double xy) {
	return dEta[pc] * xy;
}


// Zwraca Jakobian przekszta³cenia po ksi.
double Element4_2D::dXYdKsi(int i, double xy) {
	return dKsi[i] * xy;
}



