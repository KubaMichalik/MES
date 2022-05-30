#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>

#include "gauss.h"
#include "op_on_matr.h"
#include "el4_2d.h"
#include "el.h"
#include "node.h"

static double alf = 300;
static double Tsimulation = 100;
static double stepTsimulation = 50;
static double conductance = 25;
static double initialTemp = 100;
static double aT = 1200;
static double scpecH = 700;
static double density = 7800;

struct Grid {
	Element* el;
	Node* nodes;
	double height;
	double width;
	int  nH, nE, nN, nW, p;

	double** C_aggregation;
	double **H_aggregation;
	double* P_aggregation;

	Grid();
	Grid(double height, double width, int nH, int nW);
	~Grid();
		
	void turnOn(Element4_2D &e4, std::string path);
	void turnOn(Element4_2D &e4, double Tsimulation, double stepTsimulation, double conductance, double alf, double density, double scpecH, double initialTemp, double aT);
	double* GaussianElimination(double** A, double* B, int N);
	void read(std::string path, double& simulationTime, double& simulationStepTime, double& conductance, double& alf, double& density, double& scpecH, double& initialTemp, double& aT);
	void collectData(Element4_2D e4);
	void aggregation(Element &presentElement);
	void HBCcalculate(Element4_2D &e4, Element &presentElement, double alf = alf, double aT = aT);
	void Hcalculate(Element4_2D &e4, Element &presentElement, double k = conductance);
	void Ccalculate(Element4_2D &e4, Element &presentElement, double c = scpecH, double ro = density);
	void Jcalculate(Element4_2D &e4, Element &presentElement);
	void jacobian(Element4_2D &e4, Element &presentElement);
	double distance(double x1, double y1, double x2, double y2);
	void initMatrices();
	void nodesCoordsSetting();
	void printNodesCoords();
	void nodesElemSetting();
	void printElementsNodes();
	void nodesBCSetting();
	void printBCNodes();
	void prepare(Element4_2D& e4);
};