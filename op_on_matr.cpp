#include "op_on_matr.h"

using namespace std;

//Mno¿y macie¿ przez wektor i zapisuje potem w vectorC
double* multiply(double** A, double* B, const int N) {
	double *temp = new double[N];

	for (int i = 0; i < N; i++) {
		temp[i] = 0;

		for (int j = 0; j < N; j++) {
			temp[i] += A[i][j] * B[j];
		}
	}

	return temp;
}

double* multiply(double* A, double *B, const int N) {
	double* temp = new double[N];
	for (size_t i = 0; i < N; i++)
		temp[i] = A[i] * B[i];

	return temp;
}

double* multiply(double* A, double b, const int N) {
	double* temp = new double[N];
	for (size_t i = 0; i < N; i++)
		temp[i] = A[i] * b;

	return temp;
}

double* add(double* A, double* B, const int M, const int N) {
	double* temp = new double[M];
	for (size_t i = 0; i < M; i++) {
		for (size_t j = 0; j < N; j++) {
			temp[i] = A[i] + B[i];
		}
	}

	return temp;
}

double** add(double** A, double** B, const int M, const int N) {
	double** temp = new double* [M];
	for (size_t i = 0; i < M; i++) {
		temp[i] = new double[N];
		for (size_t j = 0; j < N; j++) {
			temp[i][j] = A[i][j] + B[i][j];
		}
	}

	return temp;
}

double** divide(double** A, double b, const int M, const int N) {
	double** temp = new double* [M];
	for (size_t i = 0; i < M; i++) {
		temp[i] = new double[N];
		for (size_t j = 0; j < N; j++) {
			temp[i][j] = A[i][j] / b;
		}
	}

	return temp;
}

// Scala macierz równañ z macierz¹ rozwi¹zañ w jedn¹.
double** merge(double** A, double* B, int N) {
	double** AB = new double* [N];

	for (int i = 0; i < N; i++)
		AB[i] = new double[N + 1];
	for (int i = 0; i < N; i++)	
		for (int j = 0; j < N; j++)		
			AB[i][j] = A[i][j];
	
	for (int j = 0; j < N; j++)
		AB[j][N] = B[j];
	
	return AB;
}

// Zwraca wyznacznik macierzy 2x2.
double detMrx(double matrix[WIDTH][HEIGHT]) {
	return (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]);
}

// Oblicza odwrotnoœæ przekazanej macierzy 2x2
void matrixInversion(double invM[WIDTH][HEIGHT], double M[WIDTH][HEIGHT], double detM) {
	
}

// Wype³nia zadan¹ macierz jej przetransponowan¹ wersj¹.
void transpose(double* src, double* dst, const int M, const int N) {
#pragma omp parallel for
	for (int n = 0; n < N * M; n++) {
		int i = n / N;
		int j = n % N;
		dst[n] = src[M * j + i];
	}
}



