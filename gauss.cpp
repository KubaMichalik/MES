#include "gauss.h"



// Zwraca przybliżony wynik całkowanej funkcji w przedziale <-1, 1>.
double Gauss::quad1D(vFunctionCall f, int n) {
	double result = 0;

	for (size_t i = 0; i < n + 1; i++)
		result += w[n - 1][i] * f(x[n - 1][i]);

	return result;
}

// Zwraca przybliżony wynik całkowanej funkcji w przedziale ze zmianą przedziału.
double Gauss::quad1D(double a, double b, vFunctionCall f, int n) {
	double dm = (b - a) / 2;
	double dp = (a + b) / 2;
	double result = 0;

	for (size_t i = 0; i < n + 1; i++)
		result += w[n - 1][i] * f(dm * x[n - 1][i] + dp);

	return result;
}

// Zwraca przybliżony wynik całkowanej funkcji w przestrzeni 2D.
double Gauss::quad2D(vFunctionCall2 f, int n) {
	double result = 0;

	for (size_t i = 0; i < n + 1; i++)
		for (size_t j = 0; j < n + 1; j++)
			result += w[n - 1][i] * w[n - 1][j] * f(x[n - 1][i], x[n - 1][j]);	

	return result;
}