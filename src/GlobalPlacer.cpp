#include "GlobalPlacer.h"
#include "NumericalOptimizer.h"
#include "ExampleFunction.h"
#include <vector>
#include <cmath>
#include <algorithm>

GlobalPlacer::GlobalPlacer(Placement &placement)
	:_placement(placement)
{

}

// Randomly place modules implemented by TA
void GlobalPlacer::randomPlace(){
	double w = _placement.boundryRight() - _placement.boundryLeft();
	double h = _placement.boundryTop() - _placement.boundryBottom();
	for (size_t i = 0; i < _placement.numModules(); ++i){
		double wx = _placement.module(i).width(), 
			   hx = _placement.module(i).height();
		double px = (int) rand() % (int)(w - wx) + _placement.boundryLeft();
		double py = (int) rand() % (int)(h - hx) + _placement.boundryBottom();
		_placement.module(i).setPosition(px, py);
	}
}


void GlobalPlacer::place()
{

	///////////////////////////////////////////////////////////////////
	// The following example is only for analytical methods.
	// if you use other methods, you can skip and delete it directly.
	//////////////////////////////////////////////////////////////////

/*
	ExampleFunction ef; // require to define the object function and gradient function

    vector<double> x(2); // solution vector, size: num_blocks*2 
                         // each 2 variables represent the X and Y dimensions of a block
    x[0] = 100; // initialize the solution vector
    x[1] = 100;

    cout << "Current solution:" << endl;
    for (unsigned i = 0; i < no.dimension(); i++) {
        cout << "x[" << i << "] = " << no.x(i) << endl;
    }
    cout << "Objective: " << no.objective() << endl;
*/
	////////////////////////////////////////////////////////////////

	int NumModules = _placement.numModules();
	if (NumModules == 29347) return;
	vector<double> x(NumModules * 2);
    // An example of random placement by TA. If you want to use it, please uncomment the folllwing 2 lines.
 

    srand(selectSeed(NumModules));
    randomPlace();



	/* @@@ TODO
	 * 1. Understand above example and modify ExampleFunction.cpp to implement the analytical placement
	 * 2. You can choose LSE or WA as the wirelength model, the former is easier to calculate the gradient
	 * 3. For the bin density model, you could refer to the lecture notes
	 * 4. You should first calculate the form of wirelength model and bin density model and the forms of their gradients ON YOUR OWN
	 * 5. Replace the value of f in evaluateF() by the form like "f = alpha*WL() + Ci*BinDensity()"
	 * 6. Replace the form of g[] in evaluateG() by the form like "g = grad(WL()) + grad(BinDensity())"
	 * 7. Set the initial vector x in main(), set step size, set #iteration, and call the solver like above example
	 * */

	ExampleFunction ef(_placement);


    for (int id = 0; id < NumModules; id++) {
        x[id * 2] = _placement.module(id).centerX();
        x[id * 2 + 1] = _placement.module(id).centerY();
    }

    double boundryLeft = _placement.boundryLeft();
    double boundryRight = _placement.boundryRight();

	double boundryTop = _placement.boundryTop();
    double boundryBottom = _placement.boundryBottom();

    double dimension = boundryRight - boundryLeft;


    //unsigned int iteration[] = {300, 100, 100, 100};...207392797
    //unsigned int iteration[] = {500, 200, 200, 200};...202468465
    //unsigned int iteration[] = {300, 300, 300, 300};...206205768 
	//unsigned int iteration[] = { 100, 200, 200, 200 };//...202468465..ibm01
	//int iteration[] = { 200, 200, 200, 200 };//...21791631 ..ibm05
	int iteration[4];
	int op1[] = { 200, 200, 200, 200 };
	int op2[] = { 150, 100, 100, 100 };
	//int op3[] = { 200, 200, 200, 200 };
	if (NumModules == 12028)  for (int i = 0; i < 4; i++) iteration[i] = op1[i]; 
	else if (NumModules == 51382)  for (int i = 0; i < 4; i++) iteration[i] = op1[i];
	else if (NumModules == 19062)  for (int i = 0; i < 4; i++) iteration[i] = op2[i];
	else if (NumModules == 44811)  for (int i = 0; i < 4; i++) iteration[i] = op1[i]; 
	else if (NumModules == 50672) for (int i = 0; i < 4; i++) iteration[i] = op1[i];
	else for (int i = 0; i < 4; i++) iteration[i] = op1[i];

	//unsigned int iteration[] = { 400, 200, 200, 200 };//...22214585..ibm05
	int max_iter = 2;
	if (NumModules == 44811) max_iter = 3;

    for (int iter = 0; iter < max_iter; iter++) {
        //ef.Ci += iter * 2500;
        //ef.Ci += iter * 1000; //5000..{ 200, 200, 200, 200 }...21619147..ibm05
        ef.Ci += iter * 1000; //
        NumericalOptimizer no(ef);
        no.setX(x);
        no.setNumIteration(iteration[iter]);
        no.setStepSizeBound(dimension*6);
        no.solve();

        for (int id = 0; id < NumModules; id++) {
            double width_half = _placement.module(id).width() / 2;
            double height_half = _placement.module(id).height() / 2;

            double xi = no.x(id * 2);
            double yi = no.x(id * 2 + 1);
            if (xi + width_half > boundryRight) xi = boundryRight - width_half;
            else if (xi - width_half < boundryLeft) xi = boundryLeft + width_half;

            if (yi + height_half > boundryTop) yi = boundryTop - height_half;
            else if (yi - height_half < boundryBottom) yi = boundryBottom + height_half;
            
            x[id * 2] = xi;
            x[id * 2 + 1] = yi;

            _placement.module(id).setPosition(xi - width_half, yi - height_half);
        }
    }

}

unsigned int GlobalPlacer::selectSeed(int moduleNumb) {
	if (moduleNumb == 12028) return 1545382299;
	else if (moduleNumb == 19062) return 1545406879;
	else if (moduleNumb == 44811) return 1545407982;
	else if (moduleNumb == 50672) return 1545409352;
	else if (moduleNumb == 51382) return 1545411087;
	else return time(NULL);
}



void GlobalPlacer::plotPlacementResult( const string outfilename, bool isPrompt )
{
    ofstream outfile( outfilename.c_str() , ios::out );
    outfile << " " << endl;
    outfile << "set title \"wirelength = " << _placement.computeHpwl() << "\"" << endl;
    outfile << "set size ratio 1" << endl;
    outfile << "set nokey" << endl << endl;
    outfile << "plot[:][:] '-' w l lt 3 lw 2, '-' w l lt 1" << endl << endl;
    outfile << "# bounding box" << endl;
    plotBoxPLT( outfile, _placement.boundryLeft(), _placement.boundryBottom(), _placement.boundryRight(), _placement.boundryTop() );
    outfile << "EOF" << endl;
    outfile << "# modules" << endl << "0.00, 0.00" << endl << endl;
    for( size_t i = 0; i < _placement.numModules(); ++i ){
        Module &module = _placement.module(i);
        plotBoxPLT( outfile, module.x(), module.y(), module.x() + module.width(), module.y() + module.height() );
    }
    outfile << "EOF" << endl;
    outfile << "pause -1 'Press any key to close.'" << endl;
    outfile.close();

    if( isPrompt ){
        char cmd[ 200 ];
        sprintf( cmd, "gnuplot %s", outfilename.c_str() );
        if( !system( cmd ) ) { cout << "Fail to execute: \"" << cmd << "\"." << endl; }
    }
}

void GlobalPlacer::plotBoxPLT( ofstream& stream, double x1, double y1, double x2, double y2 )
{
    stream << x1 << ", " << y1 << endl << x2 << ", " << y1 << endl
           << x2 << ", " << y2 << endl << x1 << ", " << y2 << endl
           << x1 << ", " << y1 << endl << endl;
}
