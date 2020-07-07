#pragma once
#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include "ProcessClasses.h"
class Annulus
{
	public:
		Annulus(ParameterPack pp);

		
		void Evolve();
		
		bool QuickAnalysis();
		void SaveAnnulus(std::string fileName);
		bool ValueAnalysis(bool printMode);
	private:
		ParameterPack PP;
		int NSteps;
		
		std::vector<double> Europium;
		
		std::vector<double> Europium_NSM;
		std::vector<double> Europium_Coll;
		std::vector<double> Europium_S;
		
		std::vector<double> Iron;
		std::vector<double> Magnesium;
	
		Accretion MassTracker;
		StarFormation SFRTracker;
		
		CCSN CCSNTracker;
		Collapsar CollapsarTracker;
		Decayer NSMTracker;
		Decayer SNIaTracker;
		
		//normalised constants
		double alpha;
		double beta;
		double gamma;
		double delta;
		double epsilon;
		double eta;
		
	
		void Calibrate();
		void PrintCalibration();
};




