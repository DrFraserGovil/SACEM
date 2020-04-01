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
	NThreads = 15;
	Mode = 0;
	NRandomGalaxies = 100;
	SaveValue = 1000;
	tMax = 14;
	timeStep = 0.05;
	
	tauInf = 14;
	
	
	// global variable store for command-line modification
	
	NGrid = 41;
	
	FeH_SN = RandomisableParameter<double>(-1.2,-1.5,-0.8,&global_mt);
	MgFe_SN = RandomisableParameter<double>(0.35,0.25,0.45,&global_mt);
	MgFe_Sat = RandomisableParameter<double>(-0.05,-0.1,0.1,&global_mt);
	EuMg_SN = RandomisableParameter<double>(0.05,0,0.15,&global_mt);
	sProcFrac = RandomisableParameter<double>(0.02,0.0000001,0.1,&global_mt);
	collFrac = IterableParameter<double>(0.98,0,1.0,NGrid);
		
	//constraining values
	finalEuFe_Min = -0.2;
	finalEuFe_Max = 0.1;
	finalFe_Min = -0.1;
	finalFe_Max = 0.4;
	
	//accretion/infall parameters
	
	// initial mass (10^10 solar mass)
	galaxyM0 = RandomisableParameter<double>(8.5,4.0,10.0,&global_mt);
	galaxyM1 = RandomisableParameter<double>(4.5,0.0,10.0,&global_mt);
	galaxyM2 = RandomisableParameter<double>(46,1,100.0,&global_mt);
	galaxyB1 = RandomisableParameter<double>(0.3,0.01,1.0,&global_mt);
	galaxyB2 = RandomisableParameter<double>(14.0,5,25,&global_mt);
	galaxyScaleLength = RandomisableParameter<double>(3.0,1.0,5.0,&global_mt);
	nuSFR = RandomisableParameter<double>(0.5,0.001,1.0,&global_mt);
	nuCool = RandomisableParameter<double>(1.0,0.001,2.0,&global_mt);
	alphaKS = RandomisableParameter<double>(2.3,2,2.6,&global_mt);
	hotFrac = RandomisableParameter<double>(1.0,0.3,1.0,&global_mt);
			
	massToDensityCorrection = 1;
	densityToMassCorrection = 1;
	totalToRingMassCorrection = 1;
	
	//uncalibrated stuff
	
	tauColls = IterableParameter<double>(900,0,20,NGrid);
	collWidth = RandomisableParameter<double>(2,0.01,15,&global_mt);
	tauSNIa = RandomisableParameter<double>(0.15,0.001,1,&global_mt);
	nuSNIa = RandomisableParameter<double>(0.15,0.0001,2,&global_mt);
	tauNSM = RandomisableParameter<double>(0.05,0.0001,0.8,&global_mt);
	nuNSM = RandomisableParameter<double>(0.4,0.1,3,&global_mt);
	
	WasSuccessful = false;
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
	
	UpdateRadius(Radius, Width);
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
	
	std::vector<std::string> calibrationTitles = {"[Fe/H]_0", "[Mg/Fe]_0", "[Mg/Fe]_Inf", "[Eu/Mg]_0", "s-process Fraction"};
	std::vector<RandomisableParameter<double>> calibrationParams = {FeH_SN, MgFe_SN, MgFe_Sat, EuMg_SN, sProcFrac};
		
	std::vector<std::string> galaxyTitles = {"Initial mass", "Thick disk mass", "Thin disk mass", "Thick disk infall time", "Thin disk infall time", "Galaxy Scale Length", "SFR frequency", "Cooling Frequency", "CCSN Hot Fraction"};
	std::vector<RandomisableParameter<double>> galaxyParams = {galaxyM0, galaxyM1, galaxyM2, galaxyB1, galaxyB2, galaxyScaleLength, nuSFR, nuCool, hotFrac};
	
	std::vector<std::string> processTitles = {"SNIa delay time", "SNIa frequency", "NSM delay time", "NSM frequency", "Collapsar turnoff width"};
	std::vector<RandomisableParameter<double>> processParams = {tauSNIa, nuSNIa, tauNSM, nuNSM, collWidth};

		
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
		
		std::vector<std::string> sections = {"Data Calibration Values", "Macro-Galactic Properties", "Astrophysical Process Properties"};
		std::vector<std::vector<std::string>> titles = {calibrationTitles, galaxyTitles, processTitles};
		std::vector<std::vector<RandomisableParameter<double>>> params = {calibrationParams, galaxyParams, processParams};
		
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
