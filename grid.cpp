#include "grid.h"

using namespace std;


Grid::Grid() 
{
	
	width = 0;
	height = 0;
	nW = 0;
	nH = 0;
	nN = 0;
	nE = 0;
	p = 0;
	nodes = nullptr;
	el = nullptr;
	H_aggregation = C_aggregation = nullptr;
	P_aggregation = nullptr;

	P_aggregation = new double[nN];
	H_aggregation = new double* [nN];
	C_aggregation = new double* [nN];

	for (int i = 0; i < nN; i++) {
		P_aggregation[i] = 0;
		H_aggregation[i] = new double[nN];
		C_aggregation[i] = new double[nN];

		for (int j = 0; j < nN; j++) {
			H_aggregation[i][j] = 0;
			C_aggregation[i][j] = 0;
		}
	}

	
}

Grid::Grid(double height, double width, int nH, int nW) 
{
	
	this->nE = 0;
	this->height = height;
	this->width = width;
	this->nN = nH * nW;
	this->nW = nW;
	this->nE = (nW - 1) * (nH - 1);
	this->nH = nH;
	this->nodes = new Node[nN];
	this->el = new Element[nE];

	
	P_aggregation = new double[nN];
	H_aggregation = new double* [nN];
	C_aggregation = new double* [nN];

	for (int i = 0; i < nN; i++) {
		P_aggregation[i] = 0;
		H_aggregation[i] = new double[nN];
		C_aggregation[i] = new double[nN];

		for (int j = 0; j < nN; j++) {
			H_aggregation[i][j] = 0;
			C_aggregation[i][j] = 0;
		}
	}



	nodesCoordsSetting();
	nodesElemSetting();
	nodesBCSetting();
}

Grid::~Grid() 
{	
	
	for (int i = 0; i < nN; i++) {
		delete[] C_aggregation[i];
		delete[] H_aggregation[i];
	}	
	delete[] C_aggregation;
	delete[] H_aggregation;	

}

void Grid::prepare(Element4_2D& e4) 
{
	
	for (int i = 0; i < nE; i++) {
		
		Element presentElement = el[i];
		jacobian(e4, presentElement);
		
		Hcalculate(e4, presentElement, conductance);
		Ccalculate(e4, presentElement, scpecH, density);
		HBCcalculate(e4, presentElement, alf, aT);
		
		aggregation(presentElement);
	}

}

double* Grid::GaussianElimination(double** A, double* B, int N) 
{
	double** AB = merge(A, B, N);
	const double accuracy = 1e-15;
	double* result = new double[N];
	int* vector = new int[N + 1];

	for (int i = 0; i < N + 1; i++)
		vector[i] = i;
	
	for (int i = 0; i < N - 1; i++) {
		bool hasChanged = false;
		int largest = i;

		for (int j = i + 1; j < N; j++)
			if (fabs(AB[i][vector[largest]]) < fabs(AB[i][vector[j]])) {
				hasChanged = true;
				largest = j;
			}

		if (hasChanged) {
			int pom = vector[i];
			vector[i] = vector[largest];
			vector[largest] = pom;
		}

		for (int j = i + 1; j < N; j++) {
			if (fabs(AB[i][vector[i]]) < accuracy)
				return NULL;
			
			double divisor = AB[j][vector[i]] / AB[i][vector[i]];

			for (int k = i + 1; k < N + 1; k++)
				AB[j][vector[k]] -= (AB[i][vector[k]] * divisor);
		}
	}

	for (int i = N - 1; i >= 0; i--) {
		if (fabs(AB[i][vector[i]]) < accuracy)
			return NULL;

		for (int j = N - 1; j > i; j--)
			AB[i][N] -= AB[i][vector[j]] * result[vector[j]];

		result[vector[i]] = AB[i][N] / AB[i][vector[i]];
	}

	return result;
}




//pomiary dla wczytanych danych z pliku
 
void Grid::turnOn(Element4_2D &e4, std::string path) 
{	
	
	read(path, Tsimulation, stepTsimulation, conductance, alf, density, scpecH, initialTemp, aT);
	
	P_aggregation = new double[nN];
	H_aggregation = new double* [nN];
	C_aggregation = new double* [nN];

	for (int i = 0; i < nN; i++) {
		P_aggregation[i] = 0;
		H_aggregation[i] = new double[nN];
		C_aggregation[i] = new double[nN];

		for (int j = 0; j < nN; j++) {
			H_aggregation[i][j] = 0;
			C_aggregation[i][j] = 0;
		}
	}

	turnOn(e4, Tsimulation, stepTsimulation, conductance, alf, density, scpecH, initialTemp, aT);

}

// Dokonywanie pomiarów dla wprowadzonych danych.
void Grid::turnOn(Element4_2D &e4, double Tsimulation, double stepTsimulation, double conductance, double alf, double density, double scpecH, double initialTemp, double aT) {
	
	this->prepare(e4);
	

	const int N = nN;
	double** primC = divide(C_aggregation, stepTsimulation, N, N);
	double** primH = add(H_aggregation, primC, N, N);
	double* primP, * vectorC, * t1;
	double* t0 = new double[N];
	
	for (int i = 0; i < N; i++)
	{
		t0[i] = initialTemp;
	}

	cout << "-Time[s]------Minimal Temperature[C]-----------Maximal Temperature[C]\n";
	
	for (double t = stepTsimulation, i = 0; t <= Tsimulation; t += stepTsimulation, i++) {
		
		vectorC = multiply(primC, t0, N);
		primP = add(P_aggregation, vectorC, N, N);
	
		t1 = GaussianElimination(primH, primP, N);
		t0 = t1;


		double min = 9999, max = 0;
		
		for (int i = 0; i < N; i++) {
			if (t1[i] < min) min = t1[i];
			if (t1[i] > max) max = t1[i];
		}

		
		cout << setw(3) << t << "   " << setw(18) << min << "   " << setw(35) << max << endl;
	
	}
	delete[] t0;
}

//Odczytuje plik pod wskazan¹ œcie¿k¹ i zapisuje dane do przekazanych referencji.
void Grid::read(string path, double& simulationTime, double& simulationStepTime, double& conductance, double& alf, double& density, double& scpecH, double& initialTemp, double& aT) {	
	int nodesNo = 0, elementsNo = 0;
	Element* elements = NULL;
	Node* nodes = NULL;

	string line, temp;
	ifstream file;

	
	file.open(path);
	
	if (!file.is_open())
	{
		return;
	}
	while (getline(file, line)) 
	{
		if (regex_match(line, std::regex("SimulationTime \\d+"))) 
		{
			temp = regex_replace(line, std::regex("SimulationTime "), "$1");
			simulationTime = stod(temp);
		}
		else if (regex_match(line, std::regex("(SimulationStepTime|dt) \\d+"))) 
		{
			temp = regex_replace(line, std::regex("SimulationStepTime|dt "), "$1");
			simulationStepTime = stod(temp);
		}
		else if (regex_match(line, std::regex("(conductance|K) \\d+"))) 
		{
			temp = regex_replace(line, std::regex("conductance|K "), "$1");
			conductance = stod(temp);
		}
		else if (regex_match(line, std::regex("Alfa \\d+"))) 
		{
			temp = regex_replace(line, std::regex("Alfa "), "$1");
			alf = stod(temp);
		}
		else if (regex_match(line, std::regex("Tot \\d+"))) 
		{
			temp = regex_replace(line, std::regex("Tot "), "$1");
			aT = stod(temp);
		}
		else if (regex_match(line, std::regex("InitialTemp \\d+"))) 
		{
			temp = regex_replace(line, std::regex("InitialTemp "), "$1");
			initialTemp = stod(temp);
		}
		else if (regex_match(line, std::regex("Density \\d+"))) 
		{
			temp = regex_replace(line, std::regex("Density "), "$1");
			density = stod(temp);
		}
		else if (regex_match(line, std::regex("scpecH \\d+"))) 
		{
			temp = regex_replace(line, std::regex("scpecH "), "$1");
			scpecH = stod(temp);
		}
		else if (regex_match(line, std::regex("Nodes number \\d+"))) 
		{
			temp = regex_replace(line, std::regex("Nodes number "), "$1");
			nodesNo = stoi(temp);
		}
		else if (regex_match(line, std::regex("Elements number \\d+"))) 
		{
			temp = regex_replace(line, std::regex("Elements number "), "$1");
			elementsNo = stoi(temp);
		}

		
		if (regex_match(line, std::regex("\\*Node"))) 
		{
			nodes = new Node[nodesNo];

			while (getline(file, line)) {
				if (regex_match(line, std::regex("\\s+\\d+,\\s+[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+),\\s+[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)"))) {
					istringstream ss(line);
					string token;
					string results[3];
					int i = 0;

					while (getline(ss, token, ',')) {
						results[i] = token;
						i = i > 3 ? 0 : i + 1;
					}

					int id = stoi(results[0]) - 1;
					double x = stod(results[1]);
					double y = stod(results[2]);

					nodes[id].id = id;
					nodes[id].x = x;
					nodes[id].y = y;
				}
				else
					break;
			}
		}

		
		if (regex_match(line, regex("\\*Element, type=DC2D4"))) 
		{ 
			elements = new Element[elementsNo];

			while (getline(file, line)) {
				if (regex_match(line, regex("\\s*\\d+,\\s+\\d+,\\s+\\d+,\\s+\\d+,\\s+\\d+"))) {
					istringstream ss(line);
					string token;
					string results[5];
					int i = 0;

					while (getline(ss, token, ',')) {
						results[i] = token;
						i = i > 5 ? 0 : i + 1;
					}

					int index = stoi(results[0]) - 1;
					int id1 = stoi(results[1]);
					int id2 = stoi(results[2]);
					int id3 = stoi(results[3]);
					int id4 = stoi(results[4]);

					elements[index].id[0] = id1;
					elements[index].id[1] = id2;
					elements[index].id[2] = id3;
					elements[index].id[3] = id4;
				}
				else
					break;
			}
		}

		

		if (regex_match(line, std::regex("\\*BC"))) {
			while (getline(file, line)) {
				if (nodes == NULL)
					return;

				string token;
				istringstream ss(line);
				while (getline(ss, token, ',')) {
					int IDnode = stoi(token) - 1;
					nodes[IDnode].bc = 1;
				}
			}
		}
	}
	
	
	file.close();

	this->nN = nodesNo;
	this->nE = elementsNo;
	this->nodes = nodes;
	this->el = elements;

}



void Grid::aggregation(Element &presentElement) {
	
	int H[4] = { presentElement.id[0] - 1, presentElement.id[1] - 1, presentElement.id[2] - 1, presentElement.id[3] - 1 };
	
	
	int auxiliaryMatrix[4][4][2] = {									
		
		{ { H[0], H[0] }, { H[0], H[1] }, { H[0], H[2] }, { H[0], H[3] } },
		{ { H[1], H[0] }, { H[1], H[1] }, { H[1], H[2] }, { H[1], H[3] } },
		{ { H[2], H[0] }, { H[2], H[1] }, { H[2], H[2] }, { H[2], H[3] } },
		{ { H[3], H[0] }, { H[3], H[1] }, { H[3], H[2] }, { H[3], H[3] } },

	};

	for (int i = 0; i < p; i++) {
		
		for (int j = 0; j < p; j++) {
			
			int IDrows = auxiliaryMatrix[i][j][0];
			int IDcolumns = auxiliaryMatrix[i][j][1];

			H_aggregation[IDrows][IDcolumns] += presentElement.H[i][j] + presentElement.Hbc[i][j];
			C_aggregation[IDrows][IDcolumns] += presentElement.C[i][j];

		}

		int IDrows = H[i];
		P_aggregation[IDrows] += presentElement.P[i];
	}
}




void Grid::HBCcalculate(Element4_2D &e4, Element &presentElement, double alf, double aT) {
	
	struct Side { double* pc1, * pc2; };
	double Npc1[4][4] = { 0 };
	double Npc2[4][4] = { 0 };
	double w[2] = { 1, 1 };
	double ksi;
	double eta;
	ksi = 0.5774;
	eta = 0.5774;

	double pts[8][2] = {
		{ -ksi, -1 }, { ksi, -1 }, 
		{ 1, -eta }, { 1, eta },   
		{ ksi, 1 }, { -ksi, 1 },   
		{ -1, eta }, { -1, -eta }, 
	};
	Side walls[4] = {
		{ pts[0], pts[1] }, 
		{ pts[2], pts[3] }, 
		{ pts[4], pts[5] }, 
		{ pts[6], pts[7] }, 
	};

	
	for (int i = 0; i < p; i++) {
		
		double *Npc1 = e4.nKsiEta(walls[i].pc1[0], walls[i].pc1[1]);
		double *Npc2 = e4.nKsiEta(walls[i].pc2[0], walls[i].pc2[1]);
		
		int idN1 = presentElement.id[i] - 1;
		int idN2 = presentElement.id[(i + 1) % 4] - 1;

		
		if (nodes[idN1].bc > 0 && nodes[idN2].bc > 0) {
			double detJ = sqrt(pow(nodes[idN2].x - nodes[idN1].x, 2) + pow(nodes[idN2].y - nodes[idN1].y, 2))/2;
			
			for (int j = 0; j < p; j++) {
				for (int k = 0; k < p; k++) {
					presentElement.Hbc[j][k] += (1 * (Npc1[j] * Npc1[k]) + 1 * (Npc2[j] * Npc2[k])) * alf * detJ;
				}

				presentElement.P[j] += (1 * Npc1[j] + 1 * Npc2[j]) * alf * aT * detJ;
			}
		}
	}
}



void Grid::Hcalculate(Element4_2D &e4, Element &presentElement, double k) {
	for (int i = 0; i < p; i++) {
		double rowsX[4] = { e4.dNdX[i][0],	e4.dNdX[i][1], e4.dNdX[i][2],	e4.dNdX[i][3] };
		double rowsY[4] = { e4.dNdY[i][0], e4.dNdY[i][1], e4.dNdY[i][2], e4.dNdY[i][3] };

		for (int j = 0; j < p; j++)
			for (int l = 0; l < p; l++)
				presentElement.H[j][l] += k * (rowsX[j] * rowsX[l] + rowsY[j] * rowsY[l]) * presentElement.detJ;
	}
}



void Grid::Ccalculate(Element4_2D &e4, Element &presentElement, double c, double ro) {	
	for (int i = 0; i < p; i++) {
		double* N = e4.nMrx[i];
		for (int j = 0; j < p; j++)
			for (int k = 0; k < p; k++)
				presentElement.C[j][k] += c * ro * (N[j] * N[k]) * presentElement.detJ;
	}
}




void Grid::jacobian(Element4_2D &e4, Element &presentElement) {
	this->p = e4.p;
	

	for (int i = 0; i < p; i++) {
		int IDnode = presentElement.id[i] - 1;
		double x = nodes[IDnode].x;
		double y = nodes[IDnode].y;

		presentElement.J[0][0] += e4.dXYdKsi(i, x);
		presentElement.J[0][1] += e4.dXYdEta(i, x);
		presentElement.J[1][0] += e4.dXYdKsi(i, y);
		presentElement.J[1][1] += e4.dXYdEta(i, y);
	}



	presentElement.detJ = (presentElement.J[0][0] * presentElement.J[1][1]) - (presentElement.J[0][1] * presentElement.J[1][0]);
	
	double invDetJ = 1 / presentElement.detJ;
	presentElement.invJ[0][0] = presentElement.J[1][1] * invDetJ;
	presentElement.invJ[0][1] = -presentElement.J[0][1] * invDetJ;
	presentElement.invJ[1][0] = -presentElement.J[1][0] * invDetJ;
	presentElement.invJ[1][1] = presentElement.J[0][0] * invDetJ;


	
	for (int j = 0; j < 4; j++) {
		for (int k = 0; k < 4; k++) {
			e4.dNdX[j][k] = presentElement.invJ[0][0] * e4.KsiMrx[j][k] + presentElement.invJ[0][1] * e4.EtaMrx[j][k];
			e4.dNdY[j][k] = presentElement.invJ[1][0] * e4.KsiMrx[j][k] + presentElement.invJ[1][1] * e4.EtaMrx[j][k];
		}
	}

}

void Grid::nodesBCSetting() {
	for (int i = 0; i < nN; i++) {
		if (nodes[i].x == 0 || nodes[i].x == width)  nodes[i].bc = 1;
		if (nodes[i].y == 0 || nodes[i].y == height) nodes[i].bc = 1;
	}
}



void Grid::nodesCoordsSetting() {
	int i = 0;

	for (int x = 0; x < nW; x++) {
		for (int y = 0; y < nH; y++) {
		
			nodes[i].id = i;
			++i;
			nodes[i].x = nodes[i - 1].x;
			nodes[i].y = nodes[i - 1].y + height / (static_cast<double>(nH) - 1);
		
		}

		nodes[i].x = nodes[i - 1].x + width / (static_cast<double>(nW) - 1);
		nodes[i].y = 0;
	}
}



void Grid::nodesElemSetting() {
	int n = 1, nhSum = nH;

	for (int i = 0; i < nE; i++) {
		el[i].index = i;
		el[i].id[0] = n;
		el[i].id[1] = el[i].id[0] + nH;
		el[i].id[2] = el[i].id[1] + 1;
		el[i].id[3] = el[i].id[0] + 1;
		n++;

		if (n == nhSum) {
			nhSum += nH;
			n++;
		}
	}
}





