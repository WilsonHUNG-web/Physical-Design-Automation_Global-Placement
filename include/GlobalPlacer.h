#ifndef GLOBALPLACER_H
#define GLOBALPLACER_H

#include "Placement.h"
#include <cstdlib>
#include <cstdio>
#include <fstream>

class GlobalPlacer 
{
public:
    GlobalPlacer(Placement &placement);
	void place();
    void plotPlacementResult( const string outfilename, bool isPrompt = false );
	void randomPlace(); // An example of random placement implemented by TA
	unsigned int selectSeed(int moduleNumb); // An example of random placement implemented by TA

private:
    Placement& _placement;
    void plotBoxPLT( ofstream& stream, double x1, double y1, double x2, double y2 );
};

#endif // GLOBALPLACER_H
