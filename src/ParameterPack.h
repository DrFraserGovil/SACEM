#pragma once
#include <string>
#include <vector>
#include <math.h>
class ParameterPack
{
		public:
		
			ParameterPack();
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
			
			//Iron value at SNIa turn on
			double FeH_SN;
			double FeH_SN_Min;
			double FeH_SN_Max;
			int FeH_SN_N;
			
			//[Mg/Fe] plateu value
			double MgFe_SN;
			double MgFe_SN_Min;
			double MgFe_SN_Max;
			int MgFe_SN_N;
						
			//Final [Mg/Fe] value
			double MgFe_Sat;
			double MgFe_Sat_Min;
			double MgFe_Sat_Max;
			int MgFe_Sat_N;
			
			//[Eu/Mg] plateu
			double EuMg_SN;
			double EuMg_SN_Min;
			double EuMg_SN_Max;
			int EuMg_SN_N;
			
			//fraction originating from s process
			double sProcFrac;
			double sProcFrac_Min;
			double sProcFrac_Max;
			int sProcFrac_N;
			
			//fraction originating from collapsars
			double collFrac;
			double collFrac_Min;
			double collFrac_Max;
			double collFrac_N;
			
			//constraining values
			double finalEuFe_Min;
			double finalEuFe_Max;
						
			//accretion/infall parameters
			
			// initial mass (10^10 solar mass)
			double galaxyM0;
			double galaxyM0_Min;
			double galaxyM0_Max;
			int galaxyM0_N;
			
			// thick disk infall mass (10^10 solar mass)
			double galaxyM1;
			double galaxyM1_Min;
			double galaxyM1_Max;
			int galaxyM1_N;
			
			//thin disk infall mass (10^10 solar mass)
			double galaxyM2;
			double galaxyM2_Min;
			double galaxyM2_Max;
			int galaxyM2_N;
			
			//thin disk infall timescale (Gyr)
			double galaxyB1;
			double galaxyB1_Min;
			double galaxyB1_Max;
			int galaxyB1_N;
			
			//thick disk infall timescale (Gyr)
			double galaxyB2;
			double galaxyB2_Min;
			double galaxyB2_Max;
			int galaxyB2_N;
			
			//thin disk radial scale length (kpc)
			double galaxyScaleLength;
			double galaxyScaleLength_Min;
			double galaxyScaleLength_Max;
			int galaxyScaleLength_N;
			
			//star formation efficiency parameter (per gyr)
			double nuSFR;
			double nuSFR_Min;
			double nuSFR_Max;
			int nuSFR_N;
			
			//gas cooling efficiency (per gyr)
			double nuCool;
			double nuCool_Min;
			double nuCool_Max;
			int nuCool_N;
			
			//kennicut-schmidt parameter
			double alphaKS;
			double alphaKS_Min;
			double alphaKS_Max;
			int alphaKS_N;
			
			//hotgas fraction
			double hotFrac;
			double hotFrac_Min;
			double hotFrac_Max;
			int hotFrac_N;
			
			double massToDensityCorrection;
			double densityToMassCorrection;
			double totalToRingMassCorrection;
					
					
					
			//uncalibrated stuff
			
			//collapsar cutoff time (Gyr after start)
			double tauColls;
			double tauColls_Min;
			double tauColls_Max;
			int tauColls_N;
			
			//collapsar cutoff width
			double collWidth;
			double collWidth_Min;
			double collWidth_Max;
			int collWidth_N;
			
			//snIa timescale
			double tauSNIa;
			double tauSNIa_Min;
			double tauSNIa_Max;
			int tauSNIa_N;
			
			//snIa frequency
			double nuSNIa;
			double nuSNIa_Min;
			double nuSNIa_Max;
			int nuSNIa_N;
			
			
			//NSM timescale
			double tauNSM;
			double tauNSM_Min;
			double tauNSM_Max;
			int tauNSM_N;
			
			//NSM frequency
			double nuNSM;
			double nuNSM_Min;
			double nuNSM_Max;
			int nuNSM_N;
					

			void UpdateRadius(double r, double deltaR);
		
	
};
