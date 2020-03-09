#include "ParameterPack.h"

ParameterPack::ParameterPack()
{
	//default variables
	InitialisedCorrectly = true;
	
	FILEROOT = "Output/";
	IntegrationSteps = 5000;
	NThreads = 10;
	Mode = 0;
	tMax = 14;
	timeStep = 0.03;
	
	tauInf = 14;
	
	
	// global variable store for command-line modification
	
	
	//Iron value at SNIa turn on
	FeH_SN = -1.2;
	FeH_SN_Min = -1.5;
	FeH_SN_Max = -0.8;
	FeH_SN_N = 10;
	
	//[Mg/Fe] plateu value
	MgFe_SN = 0.35;
	MgFe_SN_Min = 0.3;
	MgFe_SN_Max = 0.4;
	MgFe_SN_N = 10;
	
	//Final [Mg/Fe] value
	MgFe_Sat = -0.05;
	MgFe_Sat_Min = -0.1;
	MgFe_Sat_Max = 0.1;
	MgFe_Sat_N = 10;
	
	//[Eu/Mg] plateu
	EuMg_SN = 0.05;
	EuMg_SN_Min = -0.1;
	EuMg_SN_Max = 0.1;
	EuMg_SN_N = 20;
	
	//fraction originating from s process
	sProcFrac = 0.02;
	sProcFrac_Min = 0;
	sProcFrac_Max = 0.1;
	sProcFrac_N = 10;
	
	//fraction originating from collapsars
	collFrac = 0.98;
	collFrac_Min = 0;
	collFrac_Max = 1.0;
	collFrac_N = 100;
	
	
	//constraining values
	finalEuFe_Min = -0.2;
	finalEuFe_Max = 0.1;
	
	//accretion/infall parameters
	
	// initial mass (10^10 solar mass)
	galaxyM0 = 8.5;
	galaxyM0_Min = 4.0;
	galaxyM0_Max = 10;
	galaxyM0_N = 3;
	
	// thick disk infall mass (10^10 solar mass)
	galaxyM1 = 4.5;
	galaxyM1_Min = 0.0;
	galaxyM1_Max = 10.0;
	galaxyM1_N = 3;
	
	//thin disk infall mass (10^10 solar mass)
	galaxyM2 = 46;
	galaxyM2_Min = 10.0;
	galaxyM2_Max = 50;
	galaxyM2_N = 3;
	
	//thin disk infall timescale (Gyr)
	galaxyB1 = 0.3;
	galaxyB1_Min = 0.1;
	galaxyB1_Max = 1;
	galaxyB1_N = 3;
	
	//thick disk infall timescale (Gyr)
	galaxyB2 = 14;
	galaxyB2_Min = 5;
	galaxyB2_Max = 25;
	galaxyB2_N = 3;
	
	//thin disk radial scale length (kpc)
	galaxyScaleLength = 3.0;
	galaxyScaleLength_Min = 2.0;
	galaxyScaleLength_Max = 4.0;
	galaxyScaleLength_N = 3;
	
	//star formation efficiency parameter (per gyr)
	nuSFR = 0.5;
	nuSFR_Min = 0.1;
	nuSFR_Max = 1.0;
	nuSFR_N = 3;
	
	//gas cooling efficiency (per gyr)
	nuCool = 1.0;
	nuCool_Min = 0;
	nuCool_Max = 2.0;
	nuCool_N = 5;
	
	//kennicut-schmidt parameter
	alphaKS = 2.3;
	alphaKS_Min = 2;
	alphaKS_Max = 2.6;
	alphaKS_N = 1;
	
	//hotFrac
	hotFrac = 1.0;
	hotFrac_Min = 1.0;
	hotFrac_Max = 1.0;
	hotFrac_N = 3;
	
	massToDensityCorrection = 1;
	densityToMassCorrection = 1;
	totalToRingMassCorrection = 1;
	
	//uncalibrated stuff
	
	//collapsar cutoff time (Gyr after start)
	tauColls = 200;
	tauColls_Min = 0;
	tauColls_Max = 20;
	tauColls_N = 100;
	
	//collapsar cutoff width
	collWidth = 2;
	collWidth_Min = 0.1;
	collWidth_Max = 10;
	collWidth_N = 3;
	
	//snIa timescale
	tauSNIa = 0.15;
	tauSNIa_Min = 0;
	tauSNIa_Max = 1;
	tauSNIa_N = 10;
	
	//snIa frequency
	nuSNIa = 0.15;
	nuSNIa_Min = 0.4;
	nuSNIa_Max = 1.2;
	nuSNIa_N = 4;
	
	
	//NSM timescale
	tauNSM = 0.05;
	tauNSM_Min = 0;
	tauNSM_Max = 1;
	tauNSM_N = 10;
	
	//NSM frequency
	nuNSM = 0.4;
	nuNSM_Min = 0.2;
	nuNSM_Max = 1;
	nuNSM_N = 5;

}

void ParameterPack::UpdateRadius(double radius, double width)
{
	densityToMassCorrection = 2 * 3.141592654*radius*width;
	massToDensityCorrection = 1.0/massToDensityCorrection;
	totalToRingMassCorrection = radius*width/(galaxyScaleLength * galaxyScaleLength) * exp(-radius/galaxyScaleLength);
}
