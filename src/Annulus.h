#pragma once
#include <vector>
#include "MassReservoir.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
class PathAnnulus
{
	public:
		PathAnnulus(ParameterPack pp);
	
		void Evolve();
		
	private:
		ParameterPack PP;
		
		MassReservoir ISM;
		
		std::vector<double> Europium;
		std::vector<double> Iron;
		std::vector<double> Magnesium;
		
		//normalised constants
		double alpha;
		double beta;
		double gamma;
		double delta;
		double epsilon;
		double eta;
		
		std::vector<double> TimeVector;
		
		double Quick(double t, bool FMode);
		double SlowIntegrand(double t, double tau, double nu);
		double Slow(double t, double tau, double nu);
		void Calibrate();
		
		
};



