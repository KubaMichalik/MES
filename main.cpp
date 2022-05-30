#include <iostream>
#include "grid.h"
#include "func_to_gaus.h"
#include "el4_2d.h"


using namespace std;


// Wpisuj¹c 4 wybierasz 4x4, wpisuj¹c natomiast 31, wybierasz 31x31
#define TYPE_MES 4


int main() {	
	Element4_2D e4(2);

	#if TYPE_MES == 31
	Grid grid(0.1, 0.1, 4, 4);
		grid.turnOn(e4, "data/31_31.txt");

	#elif TYPE_MES == 4
		Grid grid(0.1, 0.1, 4, 4);
		grid.turnOn(e4, 500, 50, 25, 300, 7800, 700, 100, 1200);
		
	#endif
	
	
	system("pause");
	return 0;
}