#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <math.h>
#include <iomanip>
#include "src/ParameterPack.h"
#include "src/commandLineParser.h"
#include "src/MassReservoir.h"





int main(int argc, char** argv)
{
	
	std::cout << "S-ChEM Routine Initialised" << std::endl;
	
	//Parse command line
	ParameterPack pp = parseCommandLine(argc, argv);
	if (!(pp.InitialisedCorrectly))
	{
		std::cout << "Parse failed. Quitting" << std::endl;
		return 1;
	}
	std::cout << "Parsed command line arguments" << std::endl;
	
	double radius = 7;
	double width = 0.2;
	pp.densityToMassCorrection = 2 * 3.141592654*radius*width;
	pp.massToDensityCorrection = 1.0/pp.massToDensityCorrection;
	pp.totalToRingMassCorrection = radius*width/(pp.galaxyScaleLength*pp.galaxyScaleLength)*exp(-radius/pp.galaxyScaleLength);

	//Annulus A = Annulus(radius,width,tauInf);
	

	//~A.Evolve();
	
	double tMax = pp.tMax;
	double deltaT = pp.timeStep;
	std::vector<double> time(0);
	double t = 0;
	while (t <= tMax)
	{
		time.push_back(t);
		t+=deltaT;
	}
	

	
	MassReservoir ism = MassReservoir(time,ppCopy,true);


}
