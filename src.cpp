#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <math.h>
#include <iomanip>
#include "src/defaultValues.h"
#include "src/commandLineParser.h"
#include "src/MassReservoir.h"





int main(int argc, char** argv)
{
	
	std::cout << "S-ChEM Routine Initialised" << std::endl;
	
	//Parse command line
	bool parseSucceeded = parseCommandLine(argc, argv);
	if (!(parseSucceeded))
	{
		std::cout << "Parse failed. Quitting" << std::endl;
		return 1;
	}
	
	
	double radius = 7;
	double width = 0.2;
	densityToMassCorrection = 2 * 3.141592654*radius*width;
	massToDensityCorrection = 1.0/massToDensityCorrection;
	totalToRingMassCorrection = radius*width/(galaxyScaleLength*galaxyScaleLength)*exp(-radius/galaxyScaleLength);

	//Annulus A = Annulus(radius,width,tauInf);
	

	//~A.Evolve();
	
	double tMax = 14;
	double deltaT = 0.001;
	std::vector<double> time(0);
	double t = 0;
	while (t <= tMax)
	{
		time.push_back(t);
		t+=deltaT;
	}
	
	
	MassReservoir ism = MassReservoir(time,true);


}
