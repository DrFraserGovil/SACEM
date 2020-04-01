#pragma once
#include <string>
#include <vector>
#include <math.h>
#include <random>
#include <iostream>
#include <fstream>
#include "Parameters.h"
#include <iomanip>

extern std::random_device global_rd;
extern std::mt19937 global_mt;

class ParameterPack
{
		public:
		
			ParameterPack();
			bool InitialisedCorrectly;
			
			// global variable store for command-line modification
			std::string FILEROOT;
			int IntegrationSteps;
			int IterationSteps;
			int NThreads;
			int Mode;
			int NRandomGalaxies;
			int SaveValue;
			int NGrid;
			double tMax;
			double timeStep; 
			double tauInf;
			//calibration data
			
			RandomisableParameter<double> FeH_SN;
			RandomisableParameter<double> MgFe_SN;
			RandomisableParameter<double> MgFe_Sat;
			RandomisableParameter<double> EuMg_SN;
			RandomisableParameter<double> sProcFrac;
			IterableParameter<double> collFrac;
			
			//constraining values
			double finalEuFe_Min;
			double finalEuFe_Max;
			double finalFe_Min;
			double finalFe_Max;
				
			//accretion/infall parameter
			RandomisableParameter<double> galaxyM0;
			RandomisableParameter<double> galaxyM1;
			RandomisableParameter<double> galaxyM2;
			RandomisableParameter<double> galaxyB1;
			RandomisableParameter<double> galaxyB2;
			RandomisableParameter<double> galaxyScaleLength;
			RandomisableParameter<double> nuSFR;
			RandomisableParameter<double> nuCool;
			RandomisableParameter<double> alphaKS;
			RandomisableParameter<double> StellarLifeTimeSlope;
			RandomisableParameter<double> hotFrac;
			
			double massToDensityCorrection;
			double densityToMassCorrection;
			double totalToRingMassCorrection;
					
			//uncalibrated stuff		
			IterableParameter<double> tauColls;
			RandomisableParameter<double> collWidth;
			RandomisableParameter<double> tauSNIa;
			RandomisableParameter<double> nuSNIa;
			RandomisableParameter<double> tauNSM;
			RandomisableParameter<double> nuNSM;
	
			void UpdateRadius(double r, double deltaR);

			void ScrambleAll();
			void PrintState(std::string saveFile);
			bool WasSuccessful;
		private:
			double Radius;
			double Width;
};

