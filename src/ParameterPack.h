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
			
			bool allTimeFractionToggle;
			//calibration data
			
			RandomisableParameter<double> HFrac;
			RandomisableParameter<double> FeH_Sat;
			RandomisableParameter<double> MgFe_SN;
			RandomisableParameter<double> MgFe_Sat;
			RandomisableParameter<double> EuFe_Sat;
			RandomisableParameter<double> sProcFrac;
			IterableParameter<double> collFrac;
			
			
			
			double EuFeCeiling;	
			double maxLoopBack;
			double minGasFrac;
			double maxGasFrac;

			//parameters are { original_yVal, original_xLimit, final_yVal, final_xLimit} 
			std::vector<double> EuFeMax;
			std::vector<double> EuFeMin;
			std::vector<double> MgFeMax;
			std::vector<double> MgFeMin;
			std::vector<double> EuMgMax;
			std::vector<double> EuMgMin;
			
			//accretion/infall parameter
			RandomisableParameter<double> galaxyM0;
			RandomisableParameter<double> galaxyM1;
			RandomisableParameter<double> galaxyM2;
			RandomisableParameter<double> galaxyB1;
			RandomisableParameter<double> galaxyB2;
			
			std::vector<double> Betas;
			std::vector<double> galaxyMs;
			
			RandomisableParameter<double> galaxyScaleLength;
			RandomisableParameter<double> OutFlowFraction;
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
			
			RandomisableParameter<double> CoolingFrequency;
			RandomisableParameter<double> CollapsarCoolMod;
			RandomisableParameter<double> NSMCoolMod;
			RandomisableParameter<double> SNIaCoolMod;
	
	
			void UseLaxConstraints();
			void UseTightConstraints();
			void UseMediumConstraints();
			
			void UpdateRadius(double r, double deltaR);
			void UpdateInfall();
			void ValueChecks();

			void ScrambleAll();
			std::string PrintState();
			
			std::vector<std::string> PrinterHeaders();
			std::vector<double> PrinterValues();
			
			void SaveState(std::string saveFile);
			bool WasSuccessful;
			bool MeetsValueLimits;
			
			double OriginalTau;
			
			double nSuccess;
			std::vector<int> SuccessfulTaus;
			std::vector<int> SuccessfulFracs;
			
			std::vector<std::vector<double>> derivedParams;
			double defaultCollFrac;
			double defaultTauColls;
			
			
		private:
			double Radius;
			double Width;
};

