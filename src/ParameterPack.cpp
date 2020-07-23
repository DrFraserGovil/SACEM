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
	allTimeFractionToggle = true;
	
	//Grid Parameters
	collFrac = IterableParameter<double>(0.2,0,1.0,NGrid);
	tauColls = IterableParameter<double>(3,0,16,NGrid);
	
	
	//limiting values
	EuFeCeiling = 1.6;
	maxLoopBack = 0.25;

	//parameters are { original_yVal, original_xLimit, final_yVal, final_xLimit} 
	EuFeMax = {EuFeCeiling,-0.7,0.15,0};
	EuFeMin = {0.2,-0.9,-0.12, -0.4};

	MgFeMax = {0.45,-0.7,0.15,-0.1};
	MgFeMin = {0.2,-0.9,-0.15,-0.1};
	
	EuMgMax = {0.3,-0.7,0.2,0};
	EuMgMin = {-0.2,-0.7,-0.15,-0.4};
	
	//accretion/infall parameters
	
	// initial mass (10^10 solar mass)
	
	
	UseTightConstraints();
	
	//uncalibrated stuff
	
	
	
	
	UpdateInfall();
	double nSuccess = 0;
	WasSuccessful = false;
}



void ParameterPack::UseTightConstraints()
{
	// calibration params
	
	HFrac = RandomisableParameter<double>(0.7,0.68,0.75,&global_mt);
	FeH_Sat = RandomisableParameter<double>(0.3,0.05,0.3,&global_mt);
	MgFe_SN = RandomisableParameter<double>(0.35,0.3,0.4,&global_mt);
	MgFe_Sat = RandomisableParameter<double>(-0.05,-0.1,0,&global_mt);
	EuFe_Sat = RandomisableParameter<double>(0,-0.1,0.05,&global_mt);
	sProcFrac = RandomisableParameter<double>(0.01,0.0000001,0.05,&global_mt);
	
	
	// infall params
	galaxyM0 = RandomisableParameter<double>(0.8,0.1,1.0,&global_mt,true);
	galaxyM1 = RandomisableParameter<double>(5,1,10.0,&global_mt);
	galaxyM2 = RandomisableParameter<double>(46,20,70.0,&global_mt);
	galaxyB1 = RandomisableParameter<double>(0.3,0.1,3,&global_mt);
	galaxyB2 = RandomisableParameter<double>(14,10,30,&global_mt);

	nuSFR = RandomisableParameter<double>(1,0.05,5.0,&global_mt);
	content_modified_nuSFR = RandomisableParameter<double>(0.01,0.001,8.0,&global_mt,true);
	stellarDeathParameter = RandomisableParameter<double>(0.1,0.001,0.5,&global_mt,true);
	OutFlowFraction = RandomisableParameter<double>(2.5,0.01,1.5,&global_mt,true);
	
	
	//Process parameters
	collWidth = RandomisableParameter<double>(1,0.1,15,&global_mt);
	tauSNIa = RandomisableParameter<double>(0.15,0.1,0.3,&global_mt);
	nuSNIa = RandomisableParameter<double>(30.01,0.05,25,&global_mt);
	tauNSM = RandomisableParameter<double>(0.0001,0.0001,0.2,&global_mt,true);
	nuNSM = RandomisableParameter<double>(2.3,0.05,25,&global_mt);
	
	double hotMin = 0.7;
	double hotMax = 0.9999;
	CollapsarHotFrac = RandomisableParameter<double>(0.75,hotMin,hotMax,&global_mt);
	CCSNHotFrac = RandomisableParameter<double>(0.75,hotMin,hotMax,&global_mt);
	SNIaHotFrac = RandomisableParameter<double>(0.99,hotMin,hotMax,&global_mt);
	NSMHotFrac = RandomisableParameter<double>(0.5,0.3,hotMax,&global_mt);
	
	double mod = 0.1;
	double modMin =  - mod;
	double modMax =  + mod;
	CoolingFrequency = RandomisableParameter<double>(1,0.4,2.5,&global_mt);
	CollapsarCoolMod = RandomisableParameter<double>(0,modMin,modMax,&global_mt);
	SNIaCoolMod = RandomisableParameter<double>(0,modMin,modMax,&global_mt);
	NSMCoolMod = RandomisableParameter<double>(0,modMin,modMax,&global_mt);
	
	minGasFrac = 0.05;
	maxGasFrac = 0.15;
}

void ParameterPack::UseMediumConstraints()
{
	HFrac = RandomisableParameter<double>(0.7,0.65,0.75,&global_mt);
	FeH_Sat = RandomisableParameter<double>(0.3,0,0.5,&global_mt);
	MgFe_SN = RandomisableParameter<double>(0.35,0.3,0.5,&global_mt);
	MgFe_Sat = RandomisableParameter<double>(-0.05,-0.1,0.1,&global_mt);
	EuFe_Sat = RandomisableParameter<double>(0,-0.1,0.1,&global_mt);
	sProcFrac = RandomisableParameter<double>(0.01,0.0000001,0.1,&global_mt);
	
	
	// infall params
	galaxyM0 = RandomisableParameter<double>(0.8,0.1,10.0,&global_mt,true);
	galaxyM1 = RandomisableParameter<double>(5,1,10.0,&global_mt);
	galaxyM2 = RandomisableParameter<double>(46,10,70.0,&global_mt);
	galaxyB1 = RandomisableParameter<double>(0.3,0.1,10,&global_mt);
	galaxyB2 = RandomisableParameter<double>(14,10,50,&global_mt);

	nuSFR = RandomisableParameter<double>(1,0.05,5.0,&global_mt);
	content_modified_nuSFR = RandomisableParameter<double>(0.01,0.001,8.0,&global_mt,true);
	stellarDeathParameter = RandomisableParameter<double>(0.1,0.001,0.5,&global_mt,true);
	OutFlowFraction = RandomisableParameter<double>(2.5,0.01,2,&global_mt);
	
	
	//Process parameters
	collWidth = RandomisableParameter<double>(1,0.1,15,&global_mt);
	tauSNIa = RandomisableParameter<double>(0.15,0.05,1,&global_mt);
	nuSNIa = RandomisableParameter<double>(30.01,0.05,25,&global_mt);
	tauNSM = RandomisableParameter<double>(0.0001,0.0001,0.6,&global_mt,true);
	nuNSM = RandomisableParameter<double>(2.3,0.05,25,&global_mt);
	
	double hotMin = 0.6;
	double hotMax = 0.9999;
	CollapsarHotFrac = RandomisableParameter<double>(0.75,hotMin,hotMax,&global_mt);
	CCSNHotFrac = RandomisableParameter<double>(0.75,hotMin,hotMax,&global_mt);
	SNIaHotFrac = RandomisableParameter<double>(0.99,hotMin,hotMax,&global_mt);
	NSMHotFrac = RandomisableParameter<double>(0.5,0.3,hotMax,&global_mt);
	
	double mod = 0.2;
	double modMin =  - mod;
	double modMax =  + mod;
	CoolingFrequency = RandomisableParameter<double>(1,0.4,2.5,&global_mt);
	CollapsarCoolMod = RandomisableParameter<double>(0,modMin,modMax,&global_mt);
	SNIaCoolMod = RandomisableParameter<double>(0,modMin,modMax,&global_mt);
	NSMCoolMod = RandomisableParameter<double>(0,modMin,modMax,&global_mt);
	
	minGasFrac = 0.05;
	maxGasFrac = 0.25;
}

void ParameterPack::UseLaxConstraints()
{
	HFrac = RandomisableParameter<double>(0.7,0.1,0.9,&global_mt);
	FeH_Sat = RandomisableParameter<double>(0.3,-0.5,0.7,&global_mt);
	MgFe_SN = RandomisableParameter<double>(0.35,0.2,0.6,&global_mt);
	MgFe_Sat = RandomisableParameter<double>(-0.05,-0.2,0.2,&global_mt);
	EuFe_Sat = RandomisableParameter<double>(0,-0.3,0.2,&global_mt);
	sProcFrac = RandomisableParameter<double>(0.01,0.0000001,0.4,&global_mt);
	
	
	// infall params
	galaxyM0 = RandomisableParameter<double>(0.8,0.0001,200.0,&global_mt);
	galaxyM1 = RandomisableParameter<double>(5,0.000001,20.0,&global_mt);
	galaxyM2 = RandomisableParameter<double>(46,0.000001,200.0,&global_mt);
	galaxyB1 = RandomisableParameter<double>(0.3,0.0001,10,&global_mt);
	galaxyB2 = RandomisableParameter<double>(14,5,50,&global_mt);

	nuSFR = RandomisableParameter<double>(1,0.01,8.0,&global_mt);
	content_modified_nuSFR = RandomisableParameter<double>(0.01,0.001,8.0,&global_mt,true);
	stellarDeathParameter = RandomisableParameter<double>(0.1,0.001,2,&global_mt,true);
	OutFlowFraction = RandomisableParameter<double>(2.5,0.0001,5,&global_mt,true);
	
	
	//Process parameters
	collWidth = RandomisableParameter<double>(1,0.01,15,&global_mt);
	tauSNIa = RandomisableParameter<double>(0.15,0.001,2,&global_mt);
	nuSNIa = RandomisableParameter<double>(30.01,0.05,50,&global_mt);
	tauNSM = RandomisableParameter<double>(0.0001,0.00001,2,&global_mt,true);
	nuNSM = RandomisableParameter<double>(2.3,0.05,50,&global_mt);
	
	double hotMin = 0.001;
	double hotMax = 0.9999;
	CollapsarHotFrac = RandomisableParameter<double>(0.75,hotMin,hotMax,&global_mt);
	CCSNHotFrac = RandomisableParameter<double>(0.75,hotMin,hotMax,&global_mt);
	SNIaHotFrac = RandomisableParameter<double>(0.99,hotMin,hotMax,&global_mt);
	NSMHotFrac = RandomisableParameter<double>(0.5,0.3,hotMin,&global_mt);
	
	double mod = 0.999;
	double modMin =  - mod;
	double modMax =  + 8.0;
	CoolingFrequency = RandomisableParameter<double>(1,0.4,2.5,&global_mt);
	CollapsarCoolMod = RandomisableParameter<double>(0,modMin,modMax,&global_mt);
	SNIaCoolMod = RandomisableParameter<double>(0,modMin,modMax,&global_mt);
	NSMCoolMod = RandomisableParameter<double>(0,modMin,modMax,&global_mt);
	
	minGasFrac = 0.0000001;
	maxGasFrac = 10;
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

	nuSFR.Scramble();
	content_modified_nuSFR.Scramble();
	collWidth.Scramble();
	tauSNIa.Scramble();
	nuSNIa.Scramble();
	tauNSM.Scramble();
	nuNSM.Scramble();
	OutFlowFraction.Scramble();
	stellarDeathParameter.Scramble();
	
	CollapsarHotFrac.Scramble();
	CCSNHotFrac.Scramble();
	SNIaHotFrac.Scramble();
	NSMHotFrac.Scramble();
	
	CoolingFrequency.Scramble();
	SNIaCoolMod.Scramble();
	NSMCoolMod.Scramble();
	CollapsarCoolMod.Scramble();
	
	//UpdateRadius(Radius, Width);
	UpdateInfall();
	
	ValueChecks();
}

void ParameterPack::ValueChecks()
{
	if (collWidth.Value > tauColls.Value)
	{
		//collWidth.Value = tauColls.Value/2;
	}

	if (sProcFrac.Value + collFrac.Value > 1.0)
	{
		std::cout << "Emergency! Sprocess Reduced" << std::endl;
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


std::vector<std::string> ParameterPack::PrinterHeaders()
{
	std::vector<std::string> keyTitles = {"Collapsar Fraction", "Collapsar Turn-off"};
	std::vector<std::string> calibrationTitles = {"X", "FeH_Inf", "MgFe_0", "MgFe_Inf", "EuFe_Inf", "s-Fraction"};
	std::vector<std::string> galaxyTitles = {"M0", "M1", "M2", "b1", "b2", "nu_SFR","Mu_Stellar","OutflowFrac"};
	std::vector<std::string> processTitles = {"tau_SNIa", "nu_SNIa", "tau_NSM", "nu_NSM", "Delta_Colls","nu_modified"};
	std::vector<std::string> coolingTitles = {"f_CCSN", "f_SNIa", "f_Coll", "f_NSM", "BaseCooling", "SNIa_CoolFrac", "Coll_CoolFrac", "NSM_CoolFrac"};
	
	std::vector<std::string> output;
	std::vector<std::vector<std::string>> vecs = {calibrationTitles, galaxyTitles, processTitles,coolingTitles};
	
	for (int i = 0; i < vecs.size(); ++i)
	{
		output.insert(output.end(), vecs[i].begin(), vecs[i].end() );
	}
	return output;
}

std::vector<double> ParameterPack::PrinterValues()
{
	std::vector<RandomisableParameter<double>> calibrationParams = {HFrac, FeH_Sat, MgFe_SN, MgFe_Sat, EuFe_Sat, sProcFrac};
	std::vector<RandomisableParameter<double>> galaxyParams = {galaxyM0, galaxyM1, galaxyM2, galaxyB1, galaxyB2, nuSFR,stellarDeathParameter,OutFlowFraction};
	std::vector<RandomisableParameter<double>> processParams = {tauSNIa, nuSNIa, tauNSM, nuNSM, collWidth,content_modified_nuSFR};
	std::vector<RandomisableParameter<double>> coolingParams = {CCSNHotFrac, SNIaHotFrac, CollapsarHotFrac, NSMHotFrac, CoolingFrequency, SNIaCoolMod, CollapsarCoolMod, NSMCoolMod};

	std::vector<std::vector<RandomisableParameter<double>>> vals = {calibrationParams, galaxyParams, processParams,coolingParams};
	
	std::vector<double> output;
	
	for (int i = 0; i < vals.size(); ++i)
	{
		for(int j = 0; j < vals[i].size(); ++j)
		{
			output.push_back(vals[i][j].Value);
		}
	}

	return output;
}

std::string ParameterPack::PrintState()
{
	//~ std::ostringstream output;
	
	
	
		
	
	//~ output << "SIMULATION PARAMETERS:\n\nGrid Parameters:\n";
	//~ int width = 22;
	//~ for (int i = 0; i < keyParams.size(); ++i)
	//~ {
		//~ output << std::setw(width) << std::right << keyTitles[i] << ":";
		//~ output << std::setw(width) << std::right  << keyParams[i].Value << "\n";
	//~ }
	//~ output << "\n\n";
	
	//~ std::vector<std::string> sections = {"Data Calibration Values", "Macro-Galactic Properties", "Astrophysical Process Properties","Cooling Properties"};
	//~ std::vector<std::vector<std::string>> titles = {calibrationTitles, galaxyTitles, processTitles,coolingTitles};
	//~ std::vector<std::vector<RandomisableParameter<double>>> params = {calibrationParams, galaxyParams, processParams,coolingParams};
	
	//~ for (int i = 0; i < sections.size(); ++i)
	//~ {
		//~ output << sections[i] << "\n\n";
		//~ for (int j = 0; j < titles[i].size(); ++j)
		//~ {
			//~ output << std::setw(width) << std::right << titles[i][j] << ":";
			//~ output << std::setw(width) << std::right <<  params[i][j].Value << "\n";
		//~ }
		//~ output << "\n\n";
	//~ }
	
	//~ output << "Model Success: \t";
	//~ if (WasSuccessful)
	//~ {
		//~ output << "SUCCESS";
	//~ }
	//~ else
	//~ {
		//~ output << "FAILURE";
	//~ }
	
	//~ return output.str();

}

void ParameterPack::SaveState(std::string saveFileName)
{
	std::string output = PrintState();
	std::string name = FILEROOT + saveFileName + ".dat";
	std::ofstream saveFile;
	saveFile.open(name);
	
	if (!saveFile.is_open() )
	{
		std::cout << "Could not open file " << name << " for saving." << std::endl;
	}
	else
	{
		saveFile << output;
	}
	
	saveFile.close();
}
