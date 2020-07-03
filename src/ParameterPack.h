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
		
			int NThreads;
			int Mode;
			int NRandomGalaxies;
			int NGrid;
			double tMax;
			double timeStep; 
			double tauInf;
			//calibration data
			
			RandomisableParameter<double> HFrac;
			RandomisableParameter<double> FeH_Sat;
			RandomisableParameter<double> MgFe_SN;
			RandomisableParameter<double> MgFe_Sat;
			RandomisableParameter<double> EuFe_Sat;
			RandomisableParameter<double> sProcFrac;
			IterableParameter<double> collFrac;
			
			//constraining values
			double finalEuFe_Min;
			double finalEuFe_Max;
			double finalFe_Min;
			double finalFe_Max;
			
			double EuFeCeiling;
			double EuFeFloor;	
			double maxLoopBack;

			
			//accretion/infall parameter
			RandomisableParameter<double> galaxyM0;
			RandomisableParameter<double> galaxyM1;
			RandomisableParameter<double> galaxyM2;
			RandomisableParameter<double> galaxyB1;
			RandomisableParameter<double> galaxyB2;
			
			std::vector<double> Betas;
			std::vector<double> galaxyMs;
			
			RandomisableParameter<double> galaxyScaleLength;
			RandomisableParameter<double> nuSFR;
			RandomisableParameter<double> content_modified_nuSFR;
			RandomisableParameter<double> stellarDeathParameter;
			
			
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
	
			RandomisableParameter<double> CollapsarHotFrac;
			RandomisableParameter<double> CCSNHotFrac;
			RandomisableParameter<double> SNIaHotFrac;
			RandomisableParameter<double> NSMHotFrac;
			
			RandomisableParameter<double> CollapsarCool;
			RandomisableParameter<double> CCSNCool;
			RandomisableParameter<double> SNIaCool;
			RandomisableParameter<double> NSMCool;
	
			void UpdateRadius(double r, double deltaR);
			void UpdateInfall();
			void ValueChecks();

			void ScrambleAll();
			void PrintState(std::string saveFile);
			bool WasSuccessful;
			bool MeetsValueLimits;
			
			double defaultCollFrac;
			double defaultTauColls;
		private:
			double Radius;
			double Width;
};

