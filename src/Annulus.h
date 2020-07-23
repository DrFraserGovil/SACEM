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
		Annulus(ParameterPack * pp);
		ParameterPack * PP;

		bool QuickAnalysis();
		void SaveAnnulus(std::string fileName);
		bool ValueAnalysis(bool printMode);
		void SaveDerivedParams();
		void Calibrate();
	private:
		
		int NSteps;
		
		
	
		GalaxyMass MassTracker;
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
		
	
		
		void PrintCalibration();
};




