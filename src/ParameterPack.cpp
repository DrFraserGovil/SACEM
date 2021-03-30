#include "ParameterPack.h"



std::random_device global_rd;
std::mt19937 global_mt(global_rd());
ParameterPack::ParameterPack()
{
	//default variables
	InitialisedCorrectly = true;

	FILEROOT = "Output/";
	constraintMode = 0;
	NThreads = 4;
	Mode = 0;
	NRandomGalaxies = 2000;
	tMax = 14;
	timeStep = 0.05;
	IterationSaveInterval = 100;
	tauInf = 14;
	
	
	// global variable store for command-line modification
	

	NGrid = 101;
	
	allTimeFractionToggle = true;
	gradientSeverity = 0;
	
	
	//Grid Parameters
	collFrac = IterableParameter<double>(1,0,1.0,NGrid);
	tauColls = IterableParameter<double>(30,0,16,NGrid);
	
	
	//limiting values
	EuFeCeiling = 0.6;
	maxLoopBack = 0.05;

	//parameters are { original_yVal, original_xLimit, final_yVal, final_xLimit} 
	EuFeMax = {EuFeCeiling,-0.7,0.15,0};
	EuFeMin = {0.2,-0.9,-0.12, -0.4};

	MgFeMax = {0.45,-0.7,0.15,-0.1};
	MgFeMin = {0.2,-0.9,-0.15,-0.1};
	
	EuMgMax = {0.3,-0.7,0.2,0};
	EuMgMin = {-0.2,-0.7,-0.15,-0.4};
	
	//accretion/infall parameters
	
	// initial mass (10^10 solar mass)
	
	
	
	//uncalibrated stuff
	nominal_HFrac = 0.72;
	nominal_FeH_Sat = 0.15;
	nominal_MgFe_SN = 0.36;
	nominal_MgFe_Sat = -0.05;
	nominal_EuFe_Sat = 0;
	nominal_sProcFrac = 0.04;
	nominal_collFrac = 0.96;
	
	nominal_galaxyM0 = 3;
	nominal_galaxyM1 = 4.5;
	nominal_galaxyM2 = 29;
	nominal_galaxyB1 = 0.3;
	nominal_galaxyB2 = 14;
	
	nominal_galaxyScaleLength = 3;
	nominal_OutFlowFraction = 0.3;
	nominal_nuSFR = 0.3;
	nominal_sfrModifier = 0.3;
	nominal_stellarDeathParameter = 0.001;
	
	nominal_tauColls = 30;
	nominal_collWidth = 0.3;
	nominal_tauSNIa = 0.2;
	nominal_nuSNIa = 0.5;
	nominal_tauNSM = 0.01;
	nominal_nuNSM = 0.5;
	
	nominal_CollapsarHotFrac = 0.75;
	nominal_CCSNHotFrac = 0.75;
	nominal_SNIaHotFrac = 0.99;
	nominal_NSMHotFrac = 0.4;
	

	nominal_CoolingFrequency = 1.0;
	nominal_CollapsarCoolMod = 0.01;
	nominal_NSMCoolMod = 0.02;
	nominal_SNIaCoolMod = 0.03;

	
	double nSuccess = 0;
	WasSuccessful = false;
}

void ParameterPack::SetGradientBounds()
{
	std::vector<double> mins = {-999, -0.01};
	std::vector<double> maxes = {999,  0.01};
	
	minGradient = mins[gradientSeverity];
	maxGradient = maxes[gradientSeverity];
	
}


void ParameterPack::UseViableConstraints()
{
	
	// calibration params
	
	HFrac = RandomisableParameter<double>(nominal_HFrac,0.68,0.75,&global_mt);
	FeH_Sat = RandomisableParameter<double>(nominal_FeH_Sat,0.05,0.3,&global_mt);
	MgFe_SN = RandomisableParameter<double>(nominal_MgFe_SN,0.3,0.4,&global_mt);
	MgFe_Sat = RandomisableParameter<double>(nominal_MgFe_Sat,-0.1,0,&global_mt);
	EuFe_Sat = RandomisableParameter<double>(nominal_EuFe_Sat,-0.1,0.05,&global_mt);
	sProcFrac = RandomisableParameter<double>(nominal_sProcFrac,0.0000001,0.05,&global_mt);
	
	
	// infall params
	galaxyM0 = RandomisableParameter<double>(nominal_galaxyM0,0.1,1.0,&global_mt,true);
	galaxyM1 = RandomisableParameter<double>(nominal_galaxyM1,1,10.0,&global_mt);
	galaxyM2 = RandomisableParameter<double>(nominal_galaxyM2,20,70.0,&global_mt);
	galaxyB1 = RandomisableParameter<double>(nominal_galaxyB1,0.1,3,&global_mt);
	galaxyB2 = RandomisableParameter<double>(nominal_galaxyB2,10,30,&global_mt);

	nuSFR = RandomisableParameter<double>(nominal_nuSFR,0.05,5.0,&global_mt);
	sfrModifier = RandomisableParameter<double>(nominal_sfrModifier,0.3,1.0,&global_mt,true);
	stellarDeathParameter = RandomisableParameter<double>(nominal_stellarDeathParameter,0.0001,0.1,&global_mt,true);
	OutFlowFraction = RandomisableParameter<double>(nominal_OutFlowFraction,0.01,1.5,&global_mt,true);
	
	
	//Process parameters
	collWidth = RandomisableParameter<double>(nominal_collWidth,0.1,15,&global_mt);
	tauSNIa = RandomisableParameter<double>(nominal_tauSNIa,0.1,0.6,&global_mt);
	nuSNIa = RandomisableParameter<double>(nominal_nuSNIa,0.1,15,&global_mt);
	tauNSM = RandomisableParameter<double>(nominal_tauNSM,0.001,0.1,&global_mt,true);
	nuNSM = RandomisableParameter<double>(nominal_nuNSM,0.05,25,&global_mt);
	
	double hotMin = 0.7;
	double hotMax = 0.9999;
	CollapsarHotFrac = RandomisableParameter<double>(nominal_CollapsarHotFrac,hotMin,hotMax,&global_mt);
	CCSNHotFrac = RandomisableParameter<double>(nominal_CCSNHotFrac,hotMin,hotMax,&global_mt);
	SNIaHotFrac = RandomisableParameter<double>(nominal_SNIaHotFrac,hotMin,hotMax,&global_mt);
	NSMHotFrac = RandomisableParameter<double>(nominal_NSMHotFrac,0.3,hotMax,&global_mt);
	
	double mod = 0.1;
	double modMin =  - mod;
	double modMax =  + mod;
	CoolingFrequency = RandomisableParameter<double>(nominal_CoolingFrequency,0.5,1.5,&global_mt);
	CollapsarCoolMod = RandomisableParameter<double>(nominal_CollapsarCoolMod,modMin,modMax,&global_mt);
	SNIaCoolMod = RandomisableParameter<double>(nominal_SNIaCoolMod,modMin,modMax,&global_mt);
	NSMCoolMod = RandomisableParameter<double>(nominal_NSMCoolMod,modMin,modMax,&global_mt);
	
	minGasFrac = 0.05;
	maxGasFrac = 0.15;

}



void ParameterPack::UseLaxSFR()
{
	galaxyM0 = RandomisableParameter<double>(0.8,0.001,200.0,&global_mt);
	galaxyM1 = RandomisableParameter<double>(5,1e-5,20.0,&global_mt);
	galaxyM2 = RandomisableParameter<double>(46,1e-5,200.0,&global_mt);
	galaxyB1 = RandomisableParameter<double>(0.3,0.001,10,&global_mt);
	galaxyB2 = RandomisableParameter<double>(14,1,1000,&global_mt);

	nuSFR = RandomisableParameter<double>(0,0.01,8.0,&global_mt);
	sfrModifier = RandomisableParameter<double>(0.5,0.3,1.0,&global_mt,true);
	stellarDeathParameter = RandomisableParameter<double>(0.1,0.0001,1,&global_mt,true);
	OutFlowFraction = RandomisableParameter<double>(2.5,0.00001,5,&global_mt,true);


	CoolingFrequency = RandomisableParameter<double>(1,0.04,10,&global_mt);
	minGasFrac = 0.0000001;
	maxGasFrac = 1;

}
void ParameterPack::UseMixedConstraints()
{
	//uses tight constraints on everything but SFR stuff
	UseViableConstraints();
	
	UseLaxSFR();
	
	
}

void ParameterPack::UseWeakConstraints()
{
	HFrac = RandomisableParameter<double>(0.7,0.65,0.75,&global_mt);
	FeH_Sat = RandomisableParameter<double>(0.3,0,0.5,&global_mt);
	MgFe_SN = RandomisableParameter<double>(0.35,0.3,0.5,&global_mt);
	MgFe_Sat = RandomisableParameter<double>(-0.05,-0.1,0.1,&global_mt);
	EuFe_Sat = RandomisableParameter<double>(0,-0.1,0.1,&global_mt);
	sProcFrac = RandomisableParameter<double>(0.01,0.0000001,0.1,&global_mt);
	
	
	// infall params
	galaxyM0 = RandomisableParameter<double>(0.8,0.01,10.0,&global_mt,true);
	galaxyM1 = RandomisableParameter<double>(5,1,10.0,&global_mt);
	galaxyM2 = RandomisableParameter<double>(46,10,70.0,&global_mt);
	galaxyB1 = RandomisableParameter<double>(0.3,0.1,10,&global_mt);
	galaxyB2 = RandomisableParameter<double>(14,10,50,&global_mt);

	nuSFR = RandomisableParameter<double>(1,0.05,5.0,&global_mt);
	sfrModifier = RandomisableParameter<double>(0.5,0.3,1,&global_mt,true);
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

void ParameterPack::UseNoConstraints()
{
	UseLaxSFR();
	
	HFrac = RandomisableParameter<double>(0.7,0.5,0.9,&global_mt);
	FeH_Sat = RandomisableParameter<double>(0.3,-0.5,1,&global_mt);
	MgFe_SN = RandomisableParameter<double>(0.35,0.2,0.6,&global_mt);
	MgFe_Sat = RandomisableParameter<double>(-0.05,-0.3,0.2,&global_mt);
	EuFe_Sat = RandomisableParameter<double>(0,-0.3,0.2,&global_mt);
	sProcFrac = RandomisableParameter<double>(0.02,0.0000001,0.2,&global_mt);
	
	
	//Process parameters
	collWidth = RandomisableParameter<double>(1,0.1,15,&global_mt);
	tauSNIa = RandomisableParameter<double>(0.15,0.001,2,&global_mt);
	nuSNIa = RandomisableParameter<double>(30.01,0.05,50,&global_mt);
	tauNSM = RandomisableParameter<double>(0.0001,0.00001,2,&global_mt,true);
	nuNSM = RandomisableParameter<double>(2.3,0.05,50,&global_mt);
	
	double hotMin = 0.001;
	double hotMax = 0.9999;
	CollapsarHotFrac = RandomisableParameter<double>(0.75,hotMin,hotMax,&global_mt);
	CCSNHotFrac = RandomisableParameter<double>(0.75,hotMin,hotMax,&global_mt);
	SNIaHotFrac = RandomisableParameter<double>(0.99,hotMin,hotMax,&global_mt);
	NSMHotFrac = RandomisableParameter<double>(0.5,hotMin,hotMax,&global_mt);
	
	double mod = 0.999;
	double modMin =  - mod;
	double modMax =  + mod;
	
	CollapsarCoolMod = RandomisableParameter<double>(0,modMin,modMax,&global_mt);
	SNIaCoolMod = RandomisableParameter<double>(0,modMin,modMax,&global_mt);
	NSMCoolMod = RandomisableParameter<double>(0,modMin,modMax,&global_mt);

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
	sfrModifier.Scramble();
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


std::vector<std::string> ParameterPack::PrinterHeaders()
{
	std::vector<std::string> keyTitles = {"Collapsar Fraction", "Collapsar Turn-off"};
	std::vector<std::string> calibrationTitles = {"X", "FeH_Inf", "MgFe_0", "MgFe_Inf", "EuFe_Inf", "s-Fraction"};
	std::vector<std::string> galaxyTitles = {"M0", "M1", "M2", "b1", "b2", "nu_SFR","Mu_Stellar","OutflowFrac"};
	std::vector<std::string> processTitles = {"tau_SNIa", "nu_SNIa", "tau_NSM", "nu_NSM", "Delta_Colls","nu_modifier"};
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
	std::vector<RandomisableParameter<double>> processParams = {tauSNIa, nuSNIa, tauNSM, nuNSM, collWidth,sfrModifier};
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

void ParameterPack::LoadVariables()
{
	switch (constraintMode)
	{
		case 0:
			UseViableConstraints();
			break;
		case 1:
			UseWeakConstraints();
			break;
		case 2:
			UseNoConstraints();
			break;
		case 3:
			UseMixedConstraints();
			break;
		default:
			UseViableConstraints();
			break;
	
	}
	collFrac = IterableParameter<double>(nominal_collFrac,0,1.0,NGrid);
	tauColls = IterableParameter<double>(nominal_tauColls,0,16,NGrid);
	SetGradientBounds();
	UpdateInfall();
}
