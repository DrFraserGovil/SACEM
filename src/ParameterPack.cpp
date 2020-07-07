#include "ParameterPack.h"



std::random_device global_rd;
std::mt19937 global_mt(global_rd());
ParameterPack::ParameterPack()
{
	//default variables
	InitialisedCorrectly = true;

	FILEROOT = "Output/";
	
	NThreads = 4;
	Mode = 0;
	NRandomGalaxies = 2000;
	tMax = 14;
	timeStep = 0.05;
	
	tauInf = 14;
	
	
	// global variable store for command-line modification
	

	NGrid = 101;
	
	HFrac = RandomisableParameter<double>(0.7,0.68,0.8,&global_mt);
	FeH_Sat = RandomisableParameter<double>(0.3,0.1,0.4,&global_mt);
	MgFe_SN = RandomisableParameter<double>(0.35,0.3,0.4,&global_mt);
	MgFe_Sat = RandomisableParameter<double>(-0.05,-0.1,0.1,&global_mt);
	EuFe_Sat = RandomisableParameter<double>(0,-0.1,0.05,&global_mt);
	sProcFrac = RandomisableParameter<double>(0.01,0.0000001,0.1,&global_mt);
	collFrac = IterableParameter<double>(0.4,0,1.0,NGrid);
		
	//constraining values
	finalEuFe_Min = -0.1;
	finalEuFe_Max = 0.1;
	finalFe_Min = -0.1;
	finalFe_Max = 0.4;
	
	//limiting values
	EuFeCeiling = 0.62;
	EuFeFloor = 0.25;
	maxLoopBack = 0.05;

	
	//accretion/infall parameters
	
	// initial mass (10^10 solar mass)
	galaxyM0 = RandomisableParameter<double>(0.8,0.1,1.0,&global_mt);
	galaxyM1 = RandomisableParameter<double>(4,1,10.0,&global_mt);
	galaxyM2 = RandomisableParameter<double>(46,10,50.0,&global_mt);
	galaxyB1 = RandomisableParameter<double>(0.3,0.01,1.0,&global_mt);
	galaxyB2 = RandomisableParameter<double>(14,5,25,&global_mt);
	galaxyScaleLength = RandomisableParameter<double>(3.0,1.0,5.0,&global_mt);
	nuSFR = RandomisableParameter<double>(1.1,0.01,2.0,&global_mt);
	content_modified_nuSFR = RandomisableParameter<double>(0.01,0.01,2.0,&global_mt,true);
	stellarDeathParameter = RandomisableParameter<double>(0.02,0.0001,0.1,&global_mt,true);
	UpdateInfall();
	massToDensityCorrection = 1;
	densityToMassCorrection = 1;
	totalToRingMassCorrection = 1;
	
	//uncalibrated stuff
	
	tauColls = IterableParameter<double>(300.333,0,20,NGrid);
	collWidth = RandomisableParameter<double>(7.333,0.01,10,&global_mt);
	tauSNIa = RandomisableParameter<double>(0.15,0.05,0.5,&global_mt);
	nuSNIa = RandomisableParameter<double>(30.01,0.01,30,&global_mt);
	tauNSM = RandomisableParameter<double>(0.0001,0.00001,0.3,&global_mt,true);
	nuNSM = RandomisableParameter<double>(2.3,0.01,30,&global_mt);
	
	double hotMin = 0.4;
	double hotMax = 0.99;
	CollapsarHotFrac = RandomisableParameter<double>(0.75,hotMin,hotMax,&global_mt);
	CCSNHotFrac = RandomisableParameter<double>(0.75,hotMin,hotMax,&global_mt);
	SNIaHotFrac = RandomisableParameter<double>(0.99,hotMin,hotMax,&global_mt);
	NSMHotFrac = RandomisableParameter<double>(0.5,hotMin,hotMax,&global_mt);
	
	
	double coolMin = 0.1;
	double coolMax = 2;
	CollapsarCool = RandomisableParameter<double>(1,coolMin,coolMax,&global_mt);
	CCSNCool = RandomisableParameter<double>(1,coolMin,coolMax,&global_mt);
	SNIaCool = RandomisableParameter<double>(1.3,coolMin,coolMax,&global_mt);
	NSMCool = RandomisableParameter<double>(1,coolMin,coolMax,&global_mt);
	
	WasSuccessful = false;
}



void ParameterPack::ScrambleAll()
{
	HFrac.Scramble();
	FeH_Sat.Scramble();
	MgFe_SN.Scramble();
	MgFe_Sat.Scramble();
	EuFe_Sat.Scramble();
	sProcFrac.Scramble();
	galaxyM0.Scramble();
	galaxyM1.Scramble();
	galaxyM2.Scramble();
	galaxyB1.Scramble();
	galaxyB2.Scramble();
	galaxyScaleLength.Scramble();
	nuSFR.Scramble();
	content_modified_nuSFR.Scramble();
	collWidth.Scramble();
	tauSNIa.Scramble();
	nuSNIa.Scramble();
	tauNSM.Scramble();
	nuNSM.Scramble();
	
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
	
	ValueChecks();
}

void ParameterPack::ValueChecks()
{
	if (collWidth.Value > tauColls.Value)
	{
		collWidth.Value = tauColls.Value*0.99;
	}

	if (sProcFrac.Value + collFrac.Value > 1.0)
	{
		//std::cout << "Emergency! Sprocess Reduced" << std::endl;
		sProcFrac.Value = 1.0 - collFrac.Value;
	}
	
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



std::string ParameterPack::PrintState()
{
	std::ostringstream output;
	
	std::vector<std::string> keyTitles = {"Collapsar Fraction", "Collapsar Turn-off"};
	std::vector<IterableParameter<double>> keyParams = {collFrac, tauColls};
	
	std::vector<std::string> calibrationTitles = {"Hydrogen Fraction", "[Fe/H]_Inf", "[Mg/Fe]_0", "[Mg/Fe]_Inf", "[Eu/Fe]_Inf", "s-process Fraction"};
	std::vector<RandomisableParameter<double>> calibrationParams = {HFrac, FeH_Sat, MgFe_SN, MgFe_Sat, EuFe_Sat, sProcFrac};
		
	std::vector<std::string> galaxyTitles = {"Initial mass", "Thick disk mass", "Thin disk mass", "Thick disk infall time", "Thin disk infall time", "SFR frequency", "Cooling Frequency", "CCSN Hot Fraction","Stellar Death Mu"};
	std::vector<RandomisableParameter<double>> galaxyParams = {galaxyM0, galaxyM1, galaxyM2, galaxyB1, galaxyB2, nuSFR, CCSNCool, CCSNHotFrac,stellarDeathParameter};
	
	std::vector<std::string> processTitles = {"SNIa delay time", "SNIa frequency", "NSM delay time", "NSM frequency", "Collapsar turnoff width","Modified Removal Frequency"};
	std::vector<RandomisableParameter<double>> processParams = {tauSNIa, nuSNIa, tauNSM, nuNSM, collWidth,content_modified_nuSFR};

	std::vector<std::string> coolingTitles = {"CCSN Hot Frac", "CCSN Cooling Frequency", "SNIa hot frac", "SNIa cooling frequency", "Collapsar Hot Frac", "Collapsar Cooling Frequency", "NSM Hot Frac", "NSM Cooling Frequency"};
	std::vector<RandomisableParameter<double>> coolingParams = {CCSNHotFrac, CCSNCool, SNIaHotFrac, SNIaCool, CollapsarHotFrac, CollapsarCool, NSMHotFrac, NSMCool};
		
		
	
	output << "SIMULATION PARAMETERS:\n\nGrid Parameters:\n";
	int width = 22;
	for (int i = 0; i < keyParams.size(); ++i)
	{
		output << std::setw(width) << std::right << keyTitles[i] << ":";
		output << std::setw(width) << std::right  << keyParams[i].Value << "\n";
	}
	output << "\n\n";
	
	std::vector<std::string> sections = {"Data Calibration Values", "Macro-Galactic Properties", "Astrophysical Process Properties","Cooling Properties"};
	std::vector<std::vector<std::string>> titles = {calibrationTitles, galaxyTitles, processTitles,coolingTitles};
	std::vector<std::vector<RandomisableParameter<double>>> params = {calibrationParams, galaxyParams, processParams,coolingParams};
	
	for (int i = 0; i < sections.size(); ++i)
	{
		output << sections[i] << "\n\n";
		for (int j = 0; j < titles[i].size(); ++j)
		{
			output << std::setw(width) << std::right << titles[i][j] << ":";
			output << std::setw(width) << std::right <<  params[i][j].Value << "\n";
		}
		output << "\n\n";
	}
	
	output << "Model Success: \t";
	if (WasSuccessful)
	{
		output << "SUCCESS";
	}
	else
	{
		output << "FAILURE";
	}
	
	return output.str();

}

void ParameterPack::SaveState(std::string saveFileName)
{
	std::string output = PrintState();
	std::string name = FILEROOT + saveFileName + ".dat";
	std::ofstream saveFile;
	saveFile.open(name);
	
	if (!saveFile.is_open() )
	{
		std::cout << "Could not open file " << saveFileName << " for saving." << std::endl;
	}
	else
	{
		saveFile << output;
	}
	
	saveFile.close();
}
