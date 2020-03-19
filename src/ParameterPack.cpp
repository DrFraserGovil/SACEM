#include "ParameterPack.h"



std::random_device global_rd;
std::mt19937 global_mt(global_rd());
ParameterPack::ParameterPack()
{
	//default variables
	InitialisedCorrectly = true;

	FILEROOT = "Output/";
	IntegrationSteps = 5000;
	IterationSteps = 10;
	NThreads = 12;
	Mode = 0;
	tMax = 14;
	timeStep = 0.1;
	
	tauInf = 14;
	
	
	// global variable store for command-line modification
	
	
	//Iron value at SNIa turn on
	
	FeH_SN = RandomisableParameter<double>(-1.2,-1.5,-0.8,&global_mt);
	MgFe_SN = RandomisableParameter<double>(0.35,-1.5,-0.8,&global_mt);
	MgFe_Sat = RandomisableParameter<double>(-0.05,-0.1,0.1,&global_mt);
	EuMg_SN = RandomisableParameter<double>(0.05,-0.1,0.1,&global_mt);
	sProcFrac = RandomisableParameter<double>(0.02,0.00001,0.1,&global_mt);
	collFrac = IterableParameter<double>(0.98,0,1.0,21);
		
	//constraining values
	finalEuFe_Min = -0.2;
	finalEuFe_Max = 0.1;
	
	//accretion/infall parameters
	
	// initial mass (10^10 solar mass)
	galaxyM0 = RandomisableParameter<double>(8.5,4.0,10.0,&global_mt);
	galaxyM1 = RandomisableParameter<double>(4.5,0.0,10.0,&global_mt);
	galaxyM2 = RandomisableParameter<double>(46,10.0,50.0,&global_mt);
	galaxyB1 = RandomisableParameter<double>(0.3,0.1,1.0,&global_mt);
	galaxyB2 = RandomisableParameter<double>(14.0,5,25,&global_mt);
	galaxyScaleLength = RandomisableParameter<double>(3.0,2.0,4.0,&global_mt);
	nuSFR = RandomisableParameter<double>(0.5,0.1,1.0,&global_mt);
	nuCool = RandomisableParameter<double>(1.0,0,2.0,&global_mt);
	alphaKS = RandomisableParameter<double>(2.3,2,2.6,&global_mt);
	hotFrac = RandomisableParameter<double>(1.0,0,1.0,&global_mt);
			
	massToDensityCorrection = 1;
	densityToMassCorrection = 1;
	totalToRingMassCorrection = 1;
	
	//uncalibrated stuff
	
	tauColls = IterableParameter<double>(900,0,20,21);
	collWidth = RandomisableParameter<double>(2,0.1,10,&global_mt);
	tauSNIa = RandomisableParameter<double>(0.15,0.01,1,&global_mt);
	nuSNIa = RandomisableParameter<double>(0.15,0,2,&global_mt);
	tauNSM = RandomisableParameter<double>(0.05,0,1,&global_mt);
	nuNSM = RandomisableParameter<double>(0.4,0,1,&global_mt);
}



void ParameterPack::ScrambleAll()
{
	FeH_SN.Scramble();
	MgFe_SN.Scramble();
	MgFe_Sat.Scramble();
	EuMg_SN.Scramble();
	sProcFrac.Scramble();
	galaxyM0.Scramble();
	galaxyM1.Scramble();
	galaxyM2.Scramble();
	galaxyB1.Scramble();
	galaxyB2.Scramble();
	galaxyScaleLength.Scramble();
	nuSFR.Scramble();
	nuCool.Scramble();
	hotFrac.Scramble();
	collWidth.Scramble();
	tauSNIa.Scramble();
	nuSNIa.Scramble();
	tauNSM.Scramble();
	nuNSM.Scramble();
}

void ParameterPack::UpdateRadius(double radius, double width)
{
	densityToMassCorrection = 2 * 3.141592654*radius*width;
	massToDensityCorrection = 1.0/massToDensityCorrection;
	totalToRingMassCorrection = radius*width/(galaxyScaleLength.Value * galaxyScaleLength.Value) * exp(-radius/galaxyScaleLength.Value);
}

