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
		double STotal = SInf; //CCSNTracker.Total(tI);
		
		//collapsar counts
		double C0 = CollapsarTracker.Count(t0);
		double CInf = CollapsarTracker.Count(tI);
		double CTotal = CInf; //CollapsarTracker.Total(tI);
		
		//nsm counts
		double N0 = NSMTracker.Count(t0);
		double NInf = NSMTracker.Count(tI);
		double NTotal = NInf; //NSMTracker.Total(tI);
		
		//SNIa counts
		double WInf = SNIaTracker.Count(tI);


		double FInf = PP.FeH_Sat.Value;
		double M0 = PP.MgFe_SN.Value;
		double MInf = PP.MgFe_Sat.Value;
		double Et = PP.EuFe_Sat.Value;
		double xi = PP.sProcFrac.Value;
		double Omega = PP.collFrac.Value;
		double nsmFrac = 1.0 - xi - Omega;
		
		alpha = Hmass/SInf * pow(10.0,FInf + MInf - M0);
		
		beta = alpha * SInf/WInf * (pow(10.0,M0 - MInf) - 1.0);
		
		eta = alpha * pow(10.0,M0);
		
		

		
		//RETROFIT TO USE TOTAL GENERATION

		double denominator = NTotal * (xi * SInf/STotal + Omega * CInf/CTotal + nsmFrac*NInf/NTotal);
		
		double epsFrac = (alpha * SInf + beta * WInf) * pow(10.0,Et) /denominator;
		
		epsilon = nsmFrac * epsFrac;
		if (epsilon < 0)
		{
			epsilon = 0;
		}
		gamma = xi * epsFrac * NTotal/STotal;
		delta = Omega * epsFrac * NTotal / CTotal;
	
		double euInf = gamma * SInf + delta * CInf + epsilon * NInf;
		double feInf = alpha * SInf + beta * WInf;
		
		//~ PrintCalibration();
		//~ std::vector<std::string> vars = {"Final [Eu/H]", "Final [Fe/H]", "Final [Eu/Fe]", "Final Collapsar Frac", "Final s-Process Frac"};
		//~ std::vector<double> vals = {log10(euInf/Hmass), log10(feInf/Hmass), log10(euInf/feInf), delta*CTotal/(gamma*STotal + delta*CTotal + epsilon * NTotal), gamma*STotal/(gamma*STotal + delta*CTotal + epsilon * NTotal)};
		
		//~ for (int i = 0; i < vals.size(); ++i)
		//~ {
			//~ std::cout << vars[i] << ":\t" << vals[i] << std::endl;
		//~ }
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

bool Annulus::QuickAnalysis()
{
	//does a trivial check to see if model is successful or not
	
	
	double t = PP.tauSNIa.Value;
	
	
	double qT = CCSNTracker.Count(t);
	//double H = PP.HFrac.Value * MassTracker.Mass(t);
	
	double eu = gamma*qT + delta*CollapsarTracker.Count(t) + epsilon*NSMTracker.Count(t);
	double fe = alpha * qT + beta * SNIaTracker.Count(t);
	//double mg = eta*qT;
	
	
	if ( log10(eu/fe) > PP.EuFeCeiling)
	{
		return false;
	}


	return true;
	//~ return successCondition;	
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
	
	double t = 0;
	while (t < PP.tMax)
	{
		
		double qT = CCSNTracker.Count(t);
		double H = PP.HFrac.Value * MassTracker.Mass(t);
		
		double eu = gamma*qT + delta*CollapsarTracker.Count(t) + epsilon*NSMTracker.Count(t);
		double fe = alpha * qT + beta * SNIaTracker.Count(t);
		double mg = eta*qT;
		
		double euH = log10(eu/H);
		double feH = log10(fe/H);
		double mgH = log10(mg/H);

		double eu_S = log10(gamma*qT/H);
		double eu_NSM = log10(epsilon*NSMTracker.Count(t)/H);
		double eu_Coll = log10(delta*CollapsarTracker.Count(t)/H);
		
		double sfr = SFRTracker.Rho(t);
		double ccsnCold = CCSNTracker.Count(t);
		double ccsnAll = CCSNTracker.Total(t);
		double collapsarCold = CollapsarTracker.Count(t);
		double collapsarAll = CollapsarTracker.Total(t);
		
		std::vector<double> vals = {t, feH,mgH,euH,eu_Coll,eu_NSM, eu_S,sfr,ccsnCold,ccsnAll,collapsarCold,collapsarAll};
		for(int i = 0; i < vals.size(); ++i)
		{
			saveFile << std::setw(width) << vals[i];
		}
		saveFile << "\n";
		
		if (t < PP.tauSNIa.Value)
		{
			t+=PP.timeStep/10;
		}
		else
		{
			t+=PP.timeStep;
		}
	}
	
	saveFile.close();
}


bool Annulus::ValueAnalysis()
{
	//clean vectors
	double t = 0.02;
	
	double maxReachEu = -999999;
	double maxReachFe = -999999;
	
	bool exceededFloor = false;

	while (t <= PP.tMax)
	{

		double qT = CCSNTracker.Count(t);
		double H = PP.HFrac.Value * MassTracker.Mass(t);
		
		double eu = gamma*qT + delta*CollapsarTracker.Count(t) + epsilon*NSMTracker.Count(t);
		double fe = alpha * qT + beta * SNIaTracker.Count(t);
		double mg = eta*qT;
		
		double euH = log10(eu/H);
		double feH = log10(fe/H);
		double mgH = log10(mg/H);
		
		double minimumCheckValue = -2.5;
		if (feH > minimumCheckValue)
		{
			double eufe = euH - feH;
			if (eufe > maxReachEu)
			{
				maxReachEu = eufe;
			}
			if (feH > maxReachFe)
			{
				maxReachFe = feH;
			}
			
			
			if (eufe > PP.EuFeFloor)
			{
				exceededFloor = true;
			}
			
			
			
			bool exceededEuFeCeiling = (eufe > PP.EuFeCeiling);
			bool noDrop = (eufe > (-6.0/5.0*feH + 0.3));
			bool mgThickDiscMissing = (mgH-feH < PP.MgFe_SN.Value*0.9 & feH < -1.5);
			bool loopedBack = (feH < maxReachFe - 0.05);
			
			bool autoFail = exceededEuFeCeiling || loopedBack || noDrop ||mgThickDiscMissing;
			
			if (autoFail)
			{
				//~ std::cout << "Autofailed at t = " << t << " for :\n";
				//~ std::vector<bool> fails = {exceededEuFeCeiling,noDrop,mgThickDiscMissing, loopedBack};
				//~ std::string ceilString = "Going above [Eu/Fe] ceiling: [Eu/Fe] = " + std::to_string(eufe) + ">" + std::to_string(PP.EuFeCeiling);  
				//~ std::string noDropString = "No [Eu/Fe] drop present: [Eu/Fe] = " + std::to_string(eufe) + " at [Fe/H] = " + std::to_string(feH);
				//~ std::string mgMissString = "No [Mg/Fe] thick disc present: [Mg/Fe] = " + std::to_string(mgH - feH) + " at [FeH] = " + std::to_string(feH);
				//~ std::string loopString = "Looped back too far: [Fe/H] previously reached " + std::to_string(maxReachFe) + ", now at [Fe/H] = " + std::to_string(feH);
				//~ std::vector<std::string> reasons = {ceilString, noDropString, mgMissString, loopString};
				
				//~ for (int i = 0; i < fails.size(); ++i)
				//~ {
					//~ std::cout << i << std::endl;
					//~ if (fails[i] == true)
					//~ {
						//~ std::cout << "\t-" << reasons[i] <<std::endl;
					//~ }
				//~ }
				
				return false;
				
				
			}
		}	
		
		if (t < PP.tauSNIa.Value)
		{
			t+=PP.timeStep/10;
		}
		else
		{
			t+=PP.timeStep;
		}
	}
	
	if (exceededFloor)
	{
		return true;
	}
	else
	{
		//std::cout << "Failed at end of analysis, as did not rise above the floor" << std::endl;
		return false;
	}
	
}
