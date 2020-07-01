#pragma once
#include "Annulus.h"


Annulus::Annulus(ParameterPack pp)
{
	PP = pp;
	
	double tMax = pp.tauInf;
	double deltaT = pp.timeStep;
	
	NSteps = ceil(tMax/(deltaT));
	
	PP.timeStep = tMax/(NSteps);
	
	NSteps +=1;
	
	Europium.resize(NSteps);
	Europium_S.resize(NSteps);
	Europium_NSM.resize(NSteps);
	Europium_Coll.resize(NSteps);
	Iron.resize(NSteps);
	Magnesium.resize(NSteps);
	
	
	MassTracker = Accretion(pp);
	CCSNTracker = CCSN(pp);
	CollapsarTracker = Collapsar(pp);
	
	
	NSMTracker = Decayer(pp, pp.NSMCool.Value, pp.nuNSM.Value, pp.tauNSM.Value, pp.NSMHotFrac.Value);
	
	SFRTracker = StarFormation(pp,0,0);

	SNIaTracker = Decayer(pp, pp.SNIaCool.Value, pp.nuSNIa.Value, pp.tauSNIa.Value, pp.SNIaHotFrac.Value);
	
	Calibrate();
	
}


void Annulus::Calibrate()
{
		double tI = 15;
		double t0Cutoff = 0.002;
		//double t0 = std::max(std::min(std::min(PP.tauSNIa.Value, PP.tauColls.Value),PP.tauNSM),t0Cutoff);   // chooses the smallest time before weird stuff happens. Iff this owuld result in an answer below the cutoff, use cutoff instead
	
		double t0 = 0.002;
	
		double Hmass = PP.HFrac.Value * MassTracker.Mass(tI);
		
		//CCSN counts
		double S0 = CCSNTracker.Count(t0);
		double SInf = CCSNTracker.Count(tI);
		double St = CCSNTracker.Count(PP.tauSNIa.Value);
		double STotal = CCSNTracker.Total(tI);
		
		//collapsar counts
		double C0 = CollapsarTracker.Count(t0);
		double CInf = CollapsarTracker.Count(tI);
		double Ct = CollapsarTracker.Count(PP.tauSNIa.Value);
		double CTotal = CollapsarTracker.Total(tI);
		
		//nsm counts
		double N0 = NSMTracker.Count(t0);
		double NInf = NSMTracker.Count(tI);
		double Nt = NSMTracker.Count(PP.tauSNIa.Value);
		double NTotal = NSMTracker.Total(tI);
		
		//SNIa counts
		double WInf = SNIaTracker.Count(tI);


		double FInf = PP.FeH_Sat.Value;
		double M0 = PP.MgFe_SN.Value;
		double MInf = PP.MgFe_Sat.Value;
		double Et = PP.EuFe_SN.Value;
		double xi = PP.sProcFrac.Value;
		double Omega = PP.collFrac.Value;
		double nsmFrac = 1.0 - xi - Omega;
		
		alpha = Hmass/SInf * pow(10.0,FInf + MInf - M0);
		
		beta = alpha * SInf/WInf * (pow(10.0,M0 - MInf) - 1.0);
		
		eta = alpha * pow(10.0,M0);
		
		

		
		//RETROFIT TO USE TOTAL GENERATION

		double epsDenom = xi * St/STotal + Omega * Ct/CTotal + nsmFrac * Nt/NTotal;
	
		double epsFrac = St/NTotal * alpha * pow(10.0,Et) /epsDenom;
		
		
		epsilon = nsmFrac * epsFrac;
		delta = Omega * NTotal/CTotal * epsFrac;
		gamma = xi * NTotal/STotal * epsFrac;
	
		//PrintCalibration();
}


void Annulus::PrintCalibration()
{
	std::vector<std::string> names = {"alpha","beta", "gamma", "delta", "epsilon", "eta"};
	std::vector<double> vals = {alpha, beta, gamma, delta, epsilon, eta};
	
	for (int i = 0; i < names.size(); ++i)
	{
		std::cout << names[i] << ":\t" << vals[i] << "\n";
	}
	
}

void Annulus::Evolve()
{
	//clean vectors
	double t = 0;
	for (int i = 0; i < NSteps; ++i)
	{

		double qT = CCSNTracker.Count(t);
		double H = PP.HFrac.Value * MassTracker.Mass(t);
		
		double eu = gamma*qT + delta*CollapsarTracker.Count(t) + epsilon*NSMTracker.Count(t);
		double fe = alpha * qT + beta * SNIaTracker.Count(t);
		double mg = eta*qT;
		
		Europium[i] = log10(eu/H);
		Iron[i] = log10(fe/H);
		Magnesium[i] = log10(mg/H);

		Europium_S[i] = log10(gamma*qT/H);
		Europium_NSM[i] = log10(epsilon*NSMTracker.Count(t)/H);
		Europium_Coll[i] = log10(delta*CollapsarTracker.Count(t)/H);

		t+=PP.timeStep;
	}
}

bool Annulus::FinalStateEvaluate()
{
	double t = PP.tMax;
	
	
	double qT = CCSNTracker.Count(t);
	double H = PP.HFrac.Value * MassTracker.Mass(t);
	
	double eu = gamma*qT + delta*CollapsarTracker.Count(t) + epsilon*NSMTracker.Count(t);
	double fe = alpha * qT + beta * SNIaTracker.Count(t);
	double mg = eta*qT;
	
	
	double EuFe = log10(eu/fe);
	bool satisfiesEuropiumCondition = (EuFe >= PP.finalEuFe_Min && EuFe <= PP.finalEuFe_Max);
	bool successCondition = satisfiesEuropiumCondition ;


	return successCondition;	
}

void Annulus::SaveAnnulus(std::string fileName)
{
	std::ofstream saveFile;
	int width = 15;
	
	std::string saveFileName =PP.FILEROOT + fileName + ".dat";
	saveFile.open(saveFileName);
	std::vector<std::string> titles = {"Time","Fe/H", "Mg/H", "Eu/H","Collapsar","NSM","s-Process","SFR","CCSN (cold)", "CCSN (all)" , "Collapsar (cold)", "Collapsar (All)"};
	
	for (int i = 0; i < titles.size(); ++i)
	{
		saveFile << std::setw(width) << std::left << titles[i];  
	}
	saveFile << "\n";
	
	for (int i = 0 ; i < Iron.size(); ++i)
	{
		
		double t = i*PP.timeStep;
		double feh = Iron[i];
		double mgh = Magnesium[i];
		double euh = Europium[i];
		
		double sfr = SFRTracker.Rho(t);
		double ccsnCold = CCSNTracker.Count(t);
		double ccsnAll = CCSNTracker.Total(t);
		double collapsarCold = CollapsarTracker.Count(t);
		double collapsarAll = CollapsarTracker.Total(t);
		
		std::vector<double> vals = {i*PP.timeStep, feh,mgh,euh,Europium_Coll[i],Europium_NSM[i], Europium_S[i],sfr,ccsnCold,ccsnAll,collapsarCold,collapsarAll};
		for(int i = 0; i < vals.size(); ++i)
		{
			saveFile << std::setw(width) << vals[i];
		}
		saveFile << "\n";
	}
	
	saveFile.close();
}


bool Annulus::ValueAnalysis()
{
	//clean vectors
	double t = 0;
	for (int i = 0; i < NSteps; ++i)
	{

		double qT = CCSNTracker.Count(t);
		double H = PP.HFrac.Value * MassTracker.Mass(t);
		
		double eu = gamma*qT + delta*CollapsarTracker.Count(t) + epsilon*NSMTracker.Count(t);
		double fe = alpha * qT + beta * SNIaTracker.Count(t);
		double mg = eta*qT;
		
		double euH = log10(eu/H);
		double feH = log10(fe/H);
		double mgH = log10(mg/H);
		
		
		bool meetsEuFeCriteria = ((euH - feH) < PP.maxEuFe);
		bool meetsFeHCriteria = (feH < PP.maxFeH);
		
		bool criteriaMet = meetsEuFeCriteria & meetsFeHCriteria;
		
		if (criteriaMet == false)
		{
			return false;
		}
		
		
		t+=PP.timeStep;
	}
	
	return true;
	
	
}
