#pragma once
#include "Annulus.h"


PathAnnulus::PathAnnulus(ParameterPack pp)
{
	PP = pp;
	
	double tMax = pp.tauInf*2.02;
	double deltaT = pp.timeStep;
	double t = 0;
	while (t <= tMax)
	{
		TimeVector.push_back(t);
		t+=deltaT;
	}

	int N = TimeVector.size();
	Europium.resize(N);
	Iron.resize(N);
	Magnesium.resize(N);
	ISM = MassReservoir(TimeVector, pp, false);
	//std::cout << "SFR Initialised" << std::endl;
	
	
	Calibrate();
		//	std::cout << "Calibration completed" << std::endl;
}

PathAnnulus::PathAnnulus(ParameterPack pp, MassReservoir ism)
{
	PP = pp;
	
	double tMax = pp.tMax*1.05;
	double deltaT = pp.timeStep;
	double t = 0;
	while (t <= tMax)
	{
		TimeVector.push_back(t);
		t+=deltaT;
	}

	int N = TimeVector.size();
	Europium.resize(N);
	Iron.resize(N);
	Magnesium.resize(N);
	ISM = ism;
	Calibrate();
}



double cutoff(double t, double cutT, double wT)
{
	double lowT = cutT - wT/2;
	double upT = cutT + wT/2;
	
	if (t < lowT)
	{
		return 1;
	}
	if (t > upT)
	{
		return 0;
	}
	
	return (lowT - t)/wT + 1;
}


double PathAnnulus::Quick(double t, bool decayActive)
{
	int N = 700;
	double dt = t/(N+1);
	
	double start = ISM.ColdGas(0);
	double end = ISM.ColdGas(t);
	if(decayActive)
	{
		start*=cutoff(0,PP.tauColls.Value, PP.collWidth.Value);
		end*=cutoff(t,PP.tauColls.Value,PP.collWidth.Value);
	}
	
	double sum = 0.5*(start + end);
	
	for (double x= dt; x < t; x+=dt)
	{
		double temp = ISM.ColdGas(x);
		if (decayActive)
		{
			temp*=cutoff(x,PP.tauColls.Value,PP.collWidth.Value);
		}
		sum+=temp;
	}
	
	return sum*PP.nuSFR.Value*dt;
	
}


double PathAnnulus::SlowIntegrand(double t, double tau, double nu)
{
	int N = 300;
	double dt = t/(N+1);
	double sum = 0;
	for (double x = tau +dt; x < t; x+=dt)
	{
		sum += ISM.ColdGas(x) * exp(-nu*(x - tau));
	} 
	if (t > tau)
	{
		sum += 0.5*(ISM.ColdGas(tau)*exp(nu*tau) + ISM.ColdGas(t)*exp(-nu*(t - tau)));
	}
	return sum*dt*PP.nuSFR.Value; 
}

double PathAnnulus::Slow(double t, double tau, double nu)
{
	if (t < tau)
	{
		return 0;
	}
	int N = 300;
	double dt = t/(N+1);
	double sum = 0;
	for (double x = tau +dt; x < t; x+=dt)
	{
		sum+=SlowIntegrand(x,tau,nu);
	} 
	sum += 0.5*(SlowIntegrand(t,tau,nu) );
	return sum*dt;
}



void PathAnnulus::Calibrate()
{
		double mass = ISM.ColdGas(PP.tauSNIa.Value) + ISM.HotGas(PP.tauSNIa.Value) + ISM.Stars(PP.tauSNIa.Value);
		
		double EInf = Quick(PP.tauInf,false);
		double E0 = Quick(PP.tauSNIa.Value,false);
		
		double FInf = Quick(PP.tauInf,true);
		double F0 = Quick(PP.tauSNIa.Value,true);
		
		double HSNIa_Inf = Slow(PP.tauInf,PP.tauSNIa.Value,PP.nuSNIa.Value);
		double HNSM_0 = Slow(PP.tauSNIa.Value, PP.tauNSM.Value, PP.nuNSM.Value);	
		double HNSM_Inf = Slow(PP.tauInf, PP.tauNSM.Value, PP.nuNSM.Value);
		
		double Fcal = PP.FeH_SN.Value;
		double M0 = PP.MgFe_SN.Value;
		double MInf = PP.MgFe_Sat.Value;
		double Eps = PP.EuMg_SN.Value;
		double zeta = PP.sProcFrac.Value;
		double omega = PP.collFrac.Value;
		double nsmFrac = 1.0 - zeta - omega;
		
		alpha = mass/E0 * pow(10.0,Fcal);
		beta = alpha * EInf/HSNIa_Inf*(pow(10.0,M0 - MInf) - 1.0);
		eta = alpha * pow(10.0,M0);
		
		
		double gammaFrac = zeta*E0/EInf + omega * F0/FInf + nsmFrac*HNSM_0/HNSM_Inf;
		gamma = eta * pow(10.0,Eps) * zeta * E0/EInf / gammaFrac;
		delta = omega/zeta * EInf/FInf * gamma;
		epsilon = nsmFrac/zeta * gamma * EInf;
			
}


void PathAnnulus::Evolve()
{
	//clean vectors
	int N = TimeVector.size();
	
	int i = 0;
	MaxIndex = 0;
	double t =0;
	while (i < N && t < PP.tMax)
	{
		t = TimeVector[i];
		double qT = Quick(t,false);
		double H = 0.7 * (ISM.ColdGas(t) + ISM.HotGas(t) + ISM.Stars(t));
		double eu = gamma*qT + delta*Quick(t,true) + epsilon*Slow(t,PP.tauNSM.Value, PP.nuNSM.Value);
		double fe = alpha * qT + beta * Slow(t,PP.tauSNIa.Value, PP.nuSNIa.Value);
		double mg = eta*qT;
		
		Europium[i] = eu;
		Iron[i] = fe;
		Magnesium[i] = mg;
		
		++i;
		++MaxIndex;
	}
}

bool PathAnnulus::FinalStateEvaluate()
{
	double t = PP.tMax;
	double qT = Quick(t,false);
	double eu = gamma*qT + delta*Quick(t,true) + epsilon*Slow(t,PP.tauNSM.Value, PP.nuNSM.Value);
	double fe = alpha * qT + beta * Slow(t,PP.tauSNIa.Value, PP.nuSNIa.Value);
	
	double EuFe_Sat = log10(eu / fe);
	
	if (EuFe_Sat >= PP.finalEuFe_Min && EuFe_Sat <= PP.finalEuFe_Max)
	{
		return true;
	} 
	return false;	
}

void PathAnnulus::SaveAnnulus(std::string fileName)
{
	std::ofstream saveFile;
	int width = 15;
	
	std::string saveFileName =PP.FILEROOT + fileName + ".dat";
	saveFile.open(saveFileName);
	std::vector<std::string> titles = {"Time", "H","Fe/H", "Mg/H", "Eu/H",};
	
	for (int i = 0; i < titles.size(); ++i)
	{
		saveFile << std::setw(width) << std::left << titles[i];  
	}
	saveFile << "\n";
	
	for (int i = 0 ; i < MaxIndex; ++i)
	{
		double t = TimeVector[i];
		double H = 0.7 * (ISM.ColdGas(t) + ISM.HotGas(t) + ISM.Stars(t));
		
		double feh = log10(Iron[i]/H);
		double mgh = log10(Magnesium[i]/H);
		double euh = log10(Europium[i]/H);
		
		std::vector<double> vals = {t, H, feh,mgh,euh};
		for(int i = 0; i < vals.size(); ++i)
		{
			saveFile << std::setw(width) << vals[i];
		}
		saveFile << "\n";
	}
	
	saveFile.close();
}
