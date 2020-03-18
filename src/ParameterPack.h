#pragma once
#include <string>
#include <vector>
#include <math.h>
#include <random>

class VariedParameter
{
	public:
		VariedParameter();
		VariedParameter(double value, double min, double max, int NSteps);
		
		
		int NSteps;
		double Value;
		
		//stepthrough functions
		void UpdateValue(int stepIndex);
		double IntermediateValue(int stepIndex);
		
		//randomiser
		void RandomiseValue()
		
		
	private:
	
		double bruch;
		double MinValue;
		double MaxValue;
		int CurrentIndex;
		 std::uniform_real_distribution<double> dist;
};


class ParameterPack
{
		public:
		
			ParameterPack();
			void LosePresets();
			bool InitialisedCorrectly;
			
			// global variable store for command-line modification
			std::string FILEROOT;
			int IntegrationSteps;
			int NThreads;
			int Mode;
			double tMax;
			double timeStep; 
			double tauInf;
			//calibration data
			
			VariedParameter FeH_SN;
			VariedParameter MgFe_SN;
			VariedParameter MgFe_Sat;
			VariedParameter EuMg_SN;
			VariedParameter sProcFrac;
			VariedParameter collFrac;
			
			//constraining values
			double finalEuFe_Min;
			double finalEuFe_Max;
						
			//accretion/infall parameter
			VariedParameter galaxyM0;
			VariedParameter galaxyM1;
			VariedParameter galaxyM2;
			VariedParameter galaxyB1;
			VariedParameter galaxyB2;
			VariedParameter galaxyScaleLength;
			VariedParameter nuSFR;
			VariedParameter nuCool;
			VariedParameter alphaKS;
			VariedParameter hotFrac;
			
			double massToDensityCorrection;
			double densityToMassCorrection;
			double totalToRingMassCorrection;
					
			//uncalibrated stuff		
			VariedParameter tauColls;
			VariedParameter collWidth;
			VariedParameter tauSNIa;
			VariedParameter nuSNIa;
			VariedParameter tauNSM;
			VariedParameter nuNSM;
	
			void UpdateRadius(double r, double deltaR);
			int CountThreadLoops();
	
};

