#ifndef EXAMPLEFUNCTION_H
#define EXAMPLEFUNCTION_H
#include "Placement.h"
#include "Module.h"
#include "Net.h"
#include "Pin.h"
#include "NumericalOptimizerInterface.h"
#include <cmath>
#include <algorithm>
#include <vector>

class BIN {
public:
	double xi;
	double yi;
};

class ExampleFunction : public NumericalOptimizerInterface
{
public:
    ExampleFunction(Placement &placement);

    void evaluateFG(const vector<double> &x, double &f, vector<double> &g);
    void evaluateF(const vector<double> &x, double &f);
    unsigned dimension();

	int NumModules;
    double Ci;
    Placement &_placement;

    static const int NumBinPerRow = 16;

    double Tb; // target bin density
    double wb; // bin width
    double hb; // bin height
    double eta; //coeff for LSE
	BIN bin_matrix[NumBinPerRow][NumBinPerRow]; //bin's central coord
};
#endif // EXAMPLEFUNCTION_H
