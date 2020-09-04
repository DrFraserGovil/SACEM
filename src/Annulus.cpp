#pragma once
#include "Annulus.h"


Annulus::Annulus(ParameterPack * pp)
{
	PP = pp;
	
	PP->OriginalTau = PP->tauColls.Value;
	if (PP->collFrac.Value < 10e-6)
	{
		PP->tauColls.Value = 5;
	}
	
	double tMax = pp->tauInf;
	double deltaT = pp->timeStep;
	
		
	
	MassTracker = GalaxyMass(PP);
	CCSNTracker = CCSN(PP);
	CollapsarTracker = Collapsar(PP);
	
	double nsmCool = PP->CoolingFrequency.Value * (1.0 + PP->NSMCoolMod.Value);
	double snIaCool = PP->CoolingFrequency.Value * (1.0 + PP->SNIaCoolMod.Value);
	
	NSMTracker = Decayer(PP, nsmCool, PP->nuNSM.Value, PP->tauNSM.Value, PP->NSMHotFrac.Value);
	
	SFRTracker = StarFormation(PP,-0.1,-0.2);

	SNIaTracker = Decayer(PP, snIaCool, PP->nuSNIa.Value, PP->tauSNIa.Value, PP->SNIaHotFrac.Value);
	
	Calibrate();
	
}




void Annulus::Calibrate()
{
		//Collapsar must update as tau coll changes 
		CollapsarTracker = Collapsar(PP);
	
		double tI = PP->tMax;
		double t0Cutoff = 0.1;
		//double t0 = std::max(std::min(std::min(PP->tauSNIa.Value, PP->tauColls.Value),PP->tauNSM),t0Cutoff);   // chooses the smallest time before weird stuff happens. Iff this owuld result in an answer below the cutoff, use cutoff instead
	
		double t0 = std::min(t0Cutoff, PP->tauSNIa.Value);
	
		double Hmass = PP->HFrac.Value * MassTracker.ColdGasMass(tI);
		
		//CCSN counts
		double S0 = CCSNTracker.Count(t0);
		double SInf = CCSNTracker.Count(tI);
		double STotal = CCSNTracker.Total(tI);
		
		//collapsar counts
		double C0 = CollapsarTracker.Count(t0);
		double CInf = CollapsarTracker.Count(tI);
		double CTotal = CollapsarTracker.Total(tI);
		
		//nsm counts
		double N0 = NSMTracker.Count(t0);
		double NInf = NSMTracker.Count(tI);
		double NTotal = NSMTracker.Total(tI);
		
		//SNIa counts
		double WInf = SNIaTracker.Count(tI);

		if (PP->allTimeFractionToggle == false)
		{
			STotal = SInf;
			CTotal = CInf;
			NTotal = NInf;
		}

		double FInf = PP->FeH_Sat.Value;
		double M0 = PP->MgFe_SN.Value;
		double MInf = PP->MgFe_Sat.Value;
		double Et = PP->EuFe_Sat.Value;
		double xi = PP->sProcFrac.Value;
		double Omega = PP->collFrac.Value;
		double nsmFrac = 1.0 - xi - Omega;
		
		
		if (nsmFrac < 0)
		{
			nsmFrac = 0;
			xi = 1.0 - Omega;
		}
		
		alpha = Hmass/SInf * pow(10.0,FInf + MInf - M0);
		
		beta = alpha * SInf/WInf * (pow(10.0,M0 - MInf) - 1.0);
		
		eta = alpha * pow(10.0,M0);
		
		
		
		if (Omega == 0)
		{
			CTotal = 1;
			CInf = 1;
		}
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


void Annulus::SaveDerivedParams()
{
	double tI = PP->tMax;
	
	double cg_Inf = MassTracker.ColdGasMass(tI);
	double sg_Inf = MassTracker.StellarMass(tI);
	double Mt_Inf = MassTracker.TotalMass(tI);
	double hg_Inf = Mt_Inf - sg_Inf -cg_Inf ;
	
	

	std::vector<double> params = {cg_Inf/sg_Inf,hg_Inf/sg_Inf,sg_Inf/Mt_Inf, PP->nuSFR.Value*cg_Inf};
	
	PP->derivedParams.push_back(params);
	
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


bool Annulus::QuickAnalysis()
{
	//does a trivial check to see if model is successful or not
	
	
	double t = PP->tauSNIa.Value;
	
	
	double qT = CCSNTracker.Count(t);
	//double H = PP->HFrac.Value * MassTracker.Mass(t);
	
	double eu = gamma*qT + delta*CollapsarTracker.Count(t) + epsilon*NSMTracker.Count(t);
	double fe = alpha * qT + beta * SNIaTracker.Count(t);
	//double mg = eta*qT;
	
	
	double tI = PP->tMax;
	
	double cg_Inf = MassTracker.ColdGasMass(tI);
	double sg_Inf = MassTracker.StellarMass(tI);
	
	
	
	if ( log10(eu/fe) > PP->EuFeCeiling)
	{
		return false;
	}

	double cSRatio = cg_Inf/sg_Inf;
	if (cSRatio < PP->minGasFrac || cSRatio > PP->maxGasFrac)
	{
		return false;
	}

	//change gradient

	if (PP->gradientSeverity > 0)
	{
		double tBack = 2;
		
		double t2 = PP->tMax - tBack;
		
		double qTOld = CCSNTracker.Count(t2);
		double euOld = gamma*qTOld + delta*CollapsarTracker.Count(t2) + epsilon*NSMTracker.Count(t2);
		double feOld = alpha * qTOld + beta * SNIaTracker.Count(t2);
		
		double eufeOld = log10(euOld/feOld);
		double eufeEnd = PP->EuFe_Sat.Value;
		
		double dX = (eufeOld-eufeEnd);
		
		double gradient = dX/tBack;
		
		if (gradient < PP->minGradient || gradient > PP->maxGradient)
		{
			return false;
		}
	}
	

	return true;
	//~ return successCondition;	
}



double cutter(double x, double cut1, double cut2, double y1, double y2)
{
	if (x < cut1)
	{
		return y1;
	}
	if (x > cut2)
	{
		return y2;
	}
	
	return (y2-y1)/(cut2-cut1) * (x - cut2) + y2;
	
}

double cutter(double x, std::vector<double> cutParams)
{
	return cutter(x, cutParams[1], cutParams[3], cutParams[0], cutParams[2]);
}


void Annulus::SaveAnnulus(std::string fileName)
{
	std::ofstream saveFile;
	int width = 15;
	
	std::string saveFileName =PP->FILEROOT + fileName + ".dat";
	saveFile.open(saveFileName);
	std::vector<std::string> titles = {"Time","Fe/H", "Mg/H", "Eu/H","Collapsar","NSM","s-Process","SFR","Mcg","Mhg", "Ms" };
	
	for (int i = 0; i < titles.size(); ++i)
	{
		saveFile << std::setw(width) << std::left << titles[i] << "\t";  
	}
	saveFile << "\n";
	
	double t = 0;
	while (t < PP->tMax)
	{
		
		double qT = CCSNTracker.Count(t);
		double H = PP->HFrac.Value * MassTracker.ColdGasMass(t);
		
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
		//~ double ccsnCold = CCSNTracker.Count(t);
		//~ double ccsnAll = CCSNTracker.Total(t);
		//~ double collapsarCold = CollapsarTracker.Count(t);
		//~ double collapsarAll = CollapsarTracker.Total(t);
		
		//~ double eufeMin = cutter(feH, PP->EuFeMin);
		//~ double eufeMax = cutter(feH, PP->EuFeMax);
		
		//~ double eumgMin = cutter(feH, PP->EuMgMin);
		//~ double eumgMax = cutter(feH, PP->EuMgMax);
		
		//~ double mgfeMin = cutter(feH, PP->MgFeMin);
		//~ double mgfeMax = cutter(feH, PP->MgFeMax);
		double cgm = MassTracker.ColdGasMass(t);
		double sm = MassTracker.StellarMass(t);
		double hgm = MassTracker.TotalMass(t) - cgm - sm;
		std::vector<double> vals = {t, feH,mgH,euH,eu_Coll,eu_NSM, eu_S,sfr,cgm,  hgm, sm};
		for(int i = 0; i < vals.size(); ++i)
		{
			saveFile << std::setw(width) << vals[i] << "\t";
		}
		saveFile << "\n";
		
		if (t < PP->tauSNIa.Value)
		{
			t+=PP->timeStep/100;
		}
		else
		{
			t+=PP->timeStep;
		}
	}
	
	saveFile.close();
}



bool Annulus::ValueAnalysis(bool printMode)
{
	//clean vectors
	double t = 0.02;
	
	double maxReachEu = -999999;
	double maxReachFe = -999999;
	
	bool beganEuFeDescent = false;
	double previousEuFeValue =-999999;
	bool beganMgFeDescent = false;
	double previousMgFeValue = -9999;
	


	while (t <= PP->tMax)
	{

		double qT = CCSNTracker.Count(t);
		double H = PP->HFrac.Value * MassTracker.ColdGasMass(t);
		double fe = alpha * qT + beta * SNIaTracker.Count(t);
		double feH = log10(fe/H);
		
		
		double minimumCheckValue = -1.5;
		if (feH > minimumCheckValue)
		{
			double eu = gamma*qT + delta*CollapsarTracker.Count(t) + epsilon*NSMTracker.Count(t);
			double mg = eta*qT;
			
			double euH = log10(eu/H);
			double mgH = log10(mg/H);
			
			
			double eufe = euH - feH;
			double mgfe = mgH - feH;
			double eumg = euH - mgH;
			

			if (feH > maxReachFe)
			{
				maxReachFe = feH;
			}
						
			//check for disallowed double loops
			if ((previousEuFeValue > eufe)  & (feH > -1))
			{
				beganEuFeDescent = true;
			}
			if ( (previousMgFeValue > mgH - feH) & (feH > -1))
			{
				beganMgFeDescent = true;
			}
			
			bool mgFeRisingAgain = false;
			if (beganMgFeDescent == true)
			{
				if (previousMgFeValue < mgfe - 0.02)
				{
					mgFeRisingAgain = true;
				}
			}
			
			bool euFeRisingAgain = false;
			if (beganEuFeDescent == true)
			{
				if (previousEuFeValue < eufe-0.02)
				{
					euFeRisingAgain = true;
				}
			}
			previousEuFeValue = eufe;
			previousMgFeValue=  mgH - feH;
			
			
			
			bool outOfEuFeBounds = (eufe < cutter(feH,PP->EuFeMin) ) || (eufe > cutter(feH,PP->EuFeMax) ) || isnan(eufe);
			bool outOfMgFeBounds = (mgfe < cutter(feH,PP->MgFeMin) ) || (mgfe > cutter(feH,PP->MgFeMax) ) || isnan(mgfe);
			bool outOfEuMgBounds = (eumg < cutter(feH, PP->EuMgMin) ) || (eumg > cutter(feH, PP->EuMgMax) ) || isnan(eumg);
			

			bool loopedBack = (feH < maxReachFe - PP->maxLoopBack);
			
			
			bool autoFail = loopedBack || euFeRisingAgain || mgFeRisingAgain || outOfEuFeBounds || outOfMgFeBounds || outOfEuMgBounds;
			
			if (autoFail)
			{
				if (printMode == true)
				{
					std::cout << "Autofailed at t = " << t << " for :\n";
					std::vector<bool> fails = {loopedBack,euFeRisingAgain, mgFeRisingAgain, outOfEuFeBounds, outOfMgFeBounds, outOfEuMgBounds};
					std::string loopString = "Looped back too far: [Fe/H] previously reached " + std::to_string(maxReachFe) + ", now at [Fe/H] = " + std::to_string(feH);
					std::string euLoop = "[Eu/Fe] looped back upwards";
					std::string mgLoop = "[Mg/Fe] looped back upwards";
					std::string eufeString = "[Eu/Fe] went out of bounds: [Eu/Fe] = " + std::to_string(eufe) + " at [Fe/H] = " + std::to_string(feH) + " ( bounds = " + std::to_string(cutter(feH,PP->EuFeMin)) + " < [Eu/Fe] < " + std::to_string(cutter(feH,PP->EuFeMax) ) + ")";
					std::string mgfeString = "[Mg/Fe] went out of bounds: [Mg/Fe] = " + std::to_string(mgfe) + " at [Fe/H] = " + std::to_string(feH)+ " ( bounds = " + std::to_string(cutter(feH,PP->MgFeMin)) + " < [Mg/Fe] < " + std::to_string(cutter(feH,PP->MgFeMax) ) + ")";;
					std::string eumgString = "[Eu/Mg] went out of bounds: [Eu/Mg] = " + std::to_string(eumg) + " at [Fe/H] = " + std::to_string(feH)+ " ( bounds = " + std::to_string(cutter(feH,PP->EuMgMin)) + " < [Eu/Mg] < " + std::to_string(cutter(feH,PP->EuMgMax) ) + ")";;
					
					
					std::vector<std::string> reasons = {loopString, euLoop, mgLoop, eufeString, mgfeString, eumgString};
					
					for (int i = 0; i < fails.size(); ++i)
					{
						if (fails[i] == true)
						{
							std::cout << "\t-" << reasons[i] <<std::endl;
						}
					}
				}
				return false;
				
				
			}
		}	
		
		if (t < PP->tauSNIa.Value)
		{
			t+=PP->timeStep/5;
		}
		else
		{
			t+=PP->timeStep;
		}
	}
	
	return true;
	
}
