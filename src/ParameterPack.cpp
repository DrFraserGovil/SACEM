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
	NThreads = 4;
	Mode = 0;
	NRandomGalaxies = 200000;
	SaveValue = 1000;
	tMax = 14;
	timeStep = 0.05;
	
	tauInf = 14;
	
	
	// global variable store for command-line modification
	

	NGrid = 101;
	
	HFrac = RandomisableParameter<double>(0.714,0.68,0.8,&global_mt);
	FeH_Sat = RandomisableParameter<double>(0.12,0,0.6,&global_mt);
	MgFe_SN = RandomisableParameter<double>(0.38,0.25,0.4,&global_mt);
	MgFe_Sat = RandomisableParameter<double>(0.015,-0.1,0.1,&global_mt);
	EuFe_SN = RandomisableParameter<double>(0.31,0.3,0.5,&global_mt);
	sProcFrac = RandomisableParameter<double>(0.06,0.0000001,0.1,&global_mt);
	collFrac = IterableParameter<double>(0.92,0,1.0,NGrid);
		
	//constraining values
	finalEuFe_Min = -0.1;
	finalEuFe_Max = 0.1;
	finalFe_Min = -0.1;
	finalFe_Max = 0.4;
	
	//limiting values
	maxEuFe = 0.7;
	maxFeH = 0.6;
	
	
	//accretion/infall parameters
	
	// initial mass (10^10 solar mass)
	galaxyM0 = RandomisableParameter<double>(0.82,0.1,1.0,&global_mt);
	galaxyM1 = RandomisableParameter<double>(8.79,1,10.0,&global_mt);
	galaxyM2 = RandomisableParameter<double>(43.16,1,100.0,&global_mt);
	galaxyB1 = RandomisableParameter<double>(0.68,0.01,1.0,&global_mt);
	galaxyB2 = RandomisableParameter<double>(20.2,5,25,&global_mt);
	galaxyScaleLength = RandomisableParameter<double>(3.0,1.0,5.0,&global_mt);
	nuSFR = RandomisableParameter<double>(0.99,0.1,1.0,&global_mt);
	nuCool = RandomisableParameter<double>(0.04,0.001,2.0,&global_mt);
	alphaKS = RandomisableParameter<double>(2.3,2,2.6,&global_mt);
	StellarLifeTimeSlope =  RandomisableParameter<double>(2.3,1.5,3.5,&global_mt);
	hotFrac = RandomisableParameter<double>(0.36,0.3,1.0,&global_mt);
	stellarDeathParameter = RandomisableParameter<double>(0.19,0.000001,0.1,&global_mt);
	UpdateInfall();
	massToDensityCorrection = 1;
	densityToMassCorrection = 1;
	totalToRingMassCorrection = 1;
	
	//uncalibrated stuff
	
	tauColls = IterableParameter<double>(0.8,0,20,NGrid);
	collWidth = RandomisableParameter<double>(9.47,0.001,15,&global_mt);
	tauSNIa = RandomisableParameter<double>(0.499,0.001,1,&global_mt);
	nuSNIa = RandomisableParameter<double>(0.37,0.0001,2,&global_mt);
	tauNSM = RandomisableParameter<double>(0.7,0.0001,0.8,&global_mt);
	nuNSM = RandomisableParameter<double>(2.01,0.01,3,&global_mt);
	
	double hotMin = 0.3;
	double hotMax = 0.99;
	CollapsarHotFrac = RandomisableParameter<double>(0.505,hotMin,hotMax,&global_mt);
	CCSNHotFrac = RandomisableParameter<double>(0.82,hotMin,hotMax,&global_mt);
	SNIaHotFrac = RandomisableParameter<double>(0.73,hotMin,hotMax,&global_mt);
	NSMHotFrac = RandomisableParameter<double>(0.78,hotMin,hotMax,&global_mt);
	
	
	double coolMin = 0.1;
	double coolMax = 2;
	CollapsarCool = RandomisableParameter<double>(1.7,coolMin,coolMax,&global_mt);
	CCSNCool = RandomisableParameter<double>(0.13,coolMin,coolMax,&global_mt);
	SNIaCool = RandomisableParameter<double>(0.93,coolMin,coolMax,&global_mt);
	NSMCool = RandomisableParameter<double>(1.81,coolMin,coolMax,&global_mt);
	
	WasSuccessful = false;
}



void ParameterPack::ScrambleAll()
{
	HFrac.Scramble();
	FeH_Sat.Scramble();
	MgFe_SN.Scramble();
	MgFe_Sat.Scramble();
	EuFe_SN.Scramble();
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
	alphaKS.Scramble();
	StellarLifeTimeSlope.Scramble();
	stellarDeathParameter.Scramble();
	
	CollapsarHotFrac.Scramble();
	CCSNHotFrac.Scramble();
	SNIaHotFrac.Scramble();
	NSMHotFrac.Scramble();
	
	CollapsarCool.Scramble();
	CCSNCool.Scramble();
	SNIaCool.Scramble();
	NSMCool.Scramble();
	
	UpdateRadius(Radius, Width);
	UpdateInfall();
}

void ParameterPack::UpdateInfall()
{
	Betas = {1.0/galaxyB1.Value, 1.0/galaxyB2.Value};
	galaxyMs = {galaxyM1.Value,galaxyM2.Value};
}

void ParameterPack::UpdateRadius(double radius, double width)
{
	Radius = radius;
	Width = width;
	densityToMassCorrection = 2 * 3.141592654*radius*width;
	massToDensityCorrection = 1.0/massToDensityCorrection;
	totalToRingMassCorrection = radius*width/(galaxyScaleLength.Value * galaxyScaleLength.Value) * exp(-radius/galaxyScaleLength.Value);
}

void ParameterPack::PrintState(std::string saveFileName)
{
	std::vector<std::string> keyTitles = {"Collapsar Fraction", "Collapsar Turn-off"};
	std::vector<IterableParameter<double>> keyParams = {collFrac, tauColls};
	
	std::vector<std::string> calibrationTitles = {"Hydrogen Fraction", "[Fe/H]_Inf", "[Mg/Fe]_0", "[Mg/Fe]_Inf", "[Eu/Fe]_0", "s-process Fraction"};
	std::vector<RandomisableParameter<double>> calibrationParams = {HFrac, FeH_Sat, MgFe_SN, MgFe_Sat, EuFe_SN, sProcFrac};
		
	std::vector<std::string> galaxyTitles = {"Initial mass", "Thick disk mass", "Thin disk mass", "Thick disk infall time", "Thin disk infall time", "SFR frequency", "Cooling Frequency", "CCSN Hot Fraction","Stellar Death Mu"};
	std::vector<RandomisableParameter<double>> galaxyParams = {galaxyM0, galaxyM1, galaxyM2, galaxyB1, galaxyB2, nuSFR, nuCool, hotFrac,stellarDeathParameter};
	
	std::vector<std::string> processTitles = {"SNIa delay time", "SNIa frequency", "NSM delay time", "NSM frequency", "Collapsar turnoff width"};
	std::vector<RandomisableParameter<double>> processParams = {tauSNIa, nuSNIa, tauNSM, nuNSM, collWidth};

	std::vector<std::string> coolingTitles = {"CCSN Hot Frac", "CCSN Cooling Frequency", "SNIa hot frac", "SNIa cooling frequency", "Collapsar Hot Frac", "Collapsar Cooling Frequency", "NSM Hot Frac", "NSM Cooling Frequency"};
	std::vector<RandomisableParameter<double>> coolingParams = {CCSNHotFrac, CCSNCool, SNIaHotFrac, SNIaCool, CollapsarHotFrac, CollapsarCool, NSMHotFrac, NSMCool};
		
		
	std::ofstream saveFile;
	saveFile.open(FILEROOT + saveFileName);
	
	if (!saveFile.is_open() )
	{
		std::cout << "Could not open file " << saveFileName << " for saving." << std::endl;
	}
	else
	{
		saveFile << "SIMULATION PARAMETERS:\n\nGrid Parameters:\n";
		int width = 22;
		for (int i = 0; i < keyParams.size(); ++i)
		{
			saveFile << std::setw(width) << std::right << keyTitles[i] << ":";
			saveFile << std::setw(width) << std::right  << keyParams[i].Value << "\n";
		}
		saveFile << "\n\n";
		
		std::vector<std::string> sections = {"Data Calibration Values", "Macro-Galactic Properties", "Astrophysical Process Properties","Cooling Properties"};
		std::vector<std::vector<std::string>> titles = {calibrationTitles, galaxyTitles, processTitles,coolingTitles};
		std::vector<std::vector<RandomisableParameter<double>>> params = {calibrationParams, galaxyParams, processParams,coolingParams};
		
		for (int i = 0; i < sections.size(); ++i)
		{
			saveFile << sections[i] << "\n\n";
			for (int j = 0; j < titles[i].size(); ++j)
			{
				saveFile << std::setw(width) << std::right << titles[i][j] << ":";
				saveFile << std::setw(width) << std::right <<  params[i][j].Value << "\n";
			}
			saveFile << "\n\n";
		}
		
		saveFile << "Model Success: \t";
		if (WasSuccessful)
		{
			saveFile << "SUCCESS";
		}
		else
		{
			saveFile << "FAILURE";
		}
	}
	saveFile.close();
}
