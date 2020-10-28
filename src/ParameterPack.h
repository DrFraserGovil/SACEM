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
			int constraintMode;
			// global variable store for command-line modification
			std::string FILEROOT;
		
			int NThreads;
			int Mode;
			int NRandomGalaxies;
			int NGrid;
			int IterationSaveInterval;
			double tMax;
			double timeStep; 
			double tauInf;
			
			int gradientSeverity;
			bool allTimeFractionToggle;
			
			//observables
			double EuFeCeiling;	
			double maxLoopBack;
			double minGasFrac;
			double maxGasFrac;
			double maxGradient;
			double minGradient;
			
			//calibration data
			
			RandomisableParameter<double> HFrac;
			RandomisableParameter<double> FeH_Sat;
			RandomisableParameter<double> MgFe_SN;
			RandomisableParameter<double> MgFe_Sat;
			RandomisableParameter<double> EuFe_Sat;
			RandomisableParameter<double> sProcFrac;
			IterableParameter<double> collFrac;
			
			double nominal_HFrac;
			double nominal_FeH_Sat;
			double nominal_MgFe_SN;
			double nominal_MgFe_Sat;
			double nominal_EuFe_Sat;
			double nominal_sProcFrac;
			double nominal_collFrac;
			

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
			
			double nominal_galaxyM0;
			double nominal_galaxyM1;
			double nominal_galaxyM2;
			double nominal_galaxyB1;
			double nominal_galaxyB2;
			
			std::vector<double> Betas;
			std::vector<double> galaxyMs;
			
			RandomisableParameter<double> galaxyScaleLength;
			RandomisableParameter<double> OutFlowFraction;
			RandomisableParameter<double> nuSFR;
			RandomisableParameter<double> sfrModifier;
			RandomisableParameter<double> stellarDeathParameter;
			
			double nominal_galaxyScaleLength;
			double nominal_OutFlowFraction;
			double nominal_nuSFR;
			double nominal_sfrModifier;
			double nominal_stellarDeathParameter;
			
			
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
	
			double nominal_tauColls;
			double nominal_collWidth;
			double nominal_tauSNIa;
			double nominal_nuSNIa;
			double nominal_tauNSM;
			double nominal_nuNSM;
			
			
			RandomisableParameter<double> CollapsarHotFrac;
			RandomisableParameter<double> CCSNHotFrac;
			RandomisableParameter<double> SNIaHotFrac;
			RandomisableParameter<double> NSMHotFrac;
			
			double nominal_CollapsarHotFrac;
			double nominal_CCSNHotFrac;
			double nominal_SNIaHotFrac;
			double nominal_NSMHotFrac;
			
			RandomisableParameter<double> CoolingFrequency;
			RandomisableParameter<double> CollapsarCoolMod;
			RandomisableParameter<double> NSMCoolMod;
			RandomisableParameter<double> SNIaCoolMod;
	
			double nominal_CoolingFrequency;
			double nominal_CollapsarCoolMod;
			double nominal_NSMCoolMod;
			double nominal_SNIaCoolMod;
	
			void UseLaxConstraints();
			void UseMixedConstraints();
			void UseTightConstraints();
			void UseMediumConstraints();
			
			void SetGradientBounds();
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
			
			void LoadVariables();
		private:
			double Radius;
			double Width;
};

