#include "ParameterPack.h"


VariedParameter::VariedParameter()
{
	Value = 0.0;
	MinValue = 0.0;
	MaxValue = 0.0;
	NSteps = 0;
	bruch = 1;
	CurrentIndex = 0;
}

VariedParameter::VariedParameter(double value, double min, double max, int nSteps)
{
	Value = value;
	MinValue= min;
	MaxValue = max;
	NSteps = nSteps;
	CurrentIndex = 0;
	if (NSteps > 1)
	{
		bruch = (MaxValue - MinValue)/(NSteps - 1);
	}
	else
	{
		bruch = 0;
		MinValue = (MaxValue + MinValue)/2;
	}
	
}

int VariedParameter::Index()
{
	return CurrentIndex;
}


void VariedParameter::UpdateValue(int stepIndex)
{
	Value = MinValue + bruch*stepIndex;
	CurrentIndex = stepIndex;
}
double VariedParameter::IntermediateValue(int stepIndex)
{
	return MinValue + bruch*stepIndex;
}

ParameterPack::ParameterPack()
{
	//default variables
	InitialisedCorrectly = true;
	
	FILEROOT = "Output/";
	IntegrationSteps = 5000;
	NThreads = 10;
	Mode = 0;
	tMax = 14;
	timeStep = 0.1;
	
	tauInf = 14;
	
	
	// global variable store for command-line modification
	
	
	//Iron value at SNIa turn on
	
	FeH_SN = VariedParameter(-1.2,-1.5,-0.8,2);
	MgFe_SN = VariedParameter(0.35,0.3,0.4,2);
	MgFe_Sat = VariedParameter(-0.05,-0.1,0.1,2);
	EuMg_SN = VariedParameter(0.05,-0.1,0.1,2);
	sProcFrac = VariedParameter(0.02,0.00001,0.1,2);
	collFrac = VariedParameter(0.98,0,1.0,1);
		
	//constraining values
	finalEuFe_Min = -0.2;
	finalEuFe_Max = 0.1;
	
	//accretion/infall parameters
	
	// initial mass (10^10 solar mass)
	galaxyM0 = VariedParameter(8.5,4.0,10.0,1);
	galaxyM1 = VariedParameter(4.5,0.0,10.0,1);
	galaxyM2 = VariedParameter(46,10.0,50.0,1);
	galaxyB1 = VariedParameter(0.3,0.1,1.0,1);
	galaxyB2 = VariedParameter(14.0,5,25,3);
	galaxyScaleLength = VariedParameter(3.0,2.0,4.0,1);
	nuSFR = VariedParameter(0.5,0.1,1.0,1);
	nuCool = VariedParameter(1.0,0,2.0,1);
	alphaKS = VariedParameter(2.3,2,2.6,1);
	hotFrac = VariedParameter(1.0,0,1.0,2);
			
	massToDensityCorrection = 1;
	densityToMassCorrection = 1;
	totalToRingMassCorrection = 1;
	
	//uncalibrated stuff
	
	tauColls = VariedParameter(900,0,20,31);
	collWidth = VariedParameter(2,0.1,10,3);
	tauSNIa = VariedParameter(0.15,0.01,1,4);
	nuSNIa = VariedParameter(0.15,0,2,2);
	tauNSM = VariedParameter(0.05,0,1,2);
	nuNSM = VariedParameter(0.4,0,1,3);
}

void ParameterPack::UpdateRadius(double radius, double width)
{
	densityToMassCorrection = 2 * 3.141592654*radius*width;
	massToDensityCorrection = 1.0/massToDensityCorrection;
	totalToRingMassCorrection = radius*width/(galaxyScaleLength.Value * galaxyScaleLength.Value) * exp(-radius/galaxyScaleLength.Value);
}

int ParameterPack::CountThreadLoops()
{
	std::vector<VariedParameter> loopers = {galaxyM0, galaxyM1, galaxyM2, galaxyB1, galaxyB2, galaxyScaleLength, nuSFR, nuCool, alphaKS, hotFrac};
	int N = 1;
	for (int i = 0; i < loopers.size(); ++i)
	{
		N*=loopers[i].NSteps;
	}
	return N;
}

void ParameterPack::LosePresets()
{
	MgFe_SN.UpdateValue(0);
	MgFe_Sat.UpdateValue(0);
	EuMg_SN.UpdateValue(0);
	sProcFrac .UpdateValue(0);
	collFrac.UpdateValue(0);
		
	
	// initial mass (10^10 solar mass)
	galaxyM0.UpdateValue(0);
	galaxyM1.UpdateValue(0);
	galaxyM2.UpdateValue(0);
	galaxyB1.UpdateValue(0);
	galaxyB2.UpdateValue(0);
	galaxyScaleLength.UpdateValue(0);
	nuSFR.UpdateValue(0);
	nuCool.UpdateValue(0);
	alphaKS.UpdateValue(0);
	hotFrac.UpdateValue(0);

	
	//uncalibrated stuff
	
	tauColls.UpdateValue(0);;
	collWidth.UpdateValue(0);;
	tauSNIa.UpdateValue(0);;
	nuSNIa.UpdateValue(0);;
	tauNSM.UpdateValue(0);;
	nuNSM.UpdateValue(0);
}
