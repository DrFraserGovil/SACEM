#include "MassReservoir.h"
#include "ParameterPack.h"
#include <fstream>
#include <math.h>
#include <iomanip>
double IMF(double M,double alpha)
{
	double beta = 1.0 - 1.0/alpha;
	if (M <= 1.0)
	{
		return beta;
	}
	else
	{
		if (isinf(M))
		{
			return 0;
		}
		else
		{
			return beta*pow(M,-alpha);
		}
	}
}

double MassReservoir::accretionRate(double t)
{
	std::vector<double> Ms = {PP.galaxyM1.Value,PP.galaxyM2.Value};
	std::vector<double> Bs = {PP.galaxyB1.Value,PP.galaxyB2.Value};
	
	double mass = 0;
	for (int i = 0; i < Ms.size(); ++i)
	{
		mass += Ms[i]/Bs[i]*exp(-t/Bs[i]);
	}
	return mass*PP.totalToRingMassCorrection;
}

				
MassReservoir::MassReservoir()
{

}
MassReservoir::MassReservoir(std::vector<double> time,ParameterPack pp, bool saveFileActive)
{
	PP = pp;
	int N = time.size();
	Size = N;
	StarMass.resize(N);
	ColdGasMass.resize(N);
	HotGasMass.resize(N);
	deathRate.resize(N);
	//~ for (int i = 0; i < N ; ++i)
	//~ {	
		//~ StarMass[i] = 0;
		//~ ColdGasMass[i] = 0;
		//~ HotGasMass[i] = 0;
		//~ deathRate[i] = 0;
	//~ }
	
	StarMass[0] = 0;
	ColdGasMass[0] = 0;
	HotGasMass[0] = 0;
	deathRate[0] = 0;
	
	ColdGasMass[0] = PP.galaxyM0.Value*PP.totalToRingMassCorrection;
	
	timeStorage = time;
	
	Evolve(saveFileActive);
}

double MassReservoir::Stars(double t)
{
	int lower = invertTime(t);
	if (lower >=Size - 1)
	{
		return StarMass[Size-1];
	}
	int upper = lower + 1;

	
	double bruch = (1.0*t - 1.0*timeStorage[lower])/(timeStorage[upper] - timeStorage[lower]);
	return (StarMass[upper] - StarMass[lower])*bruch + StarMass[lower];
}
double MassReservoir::ColdGas(double t)
{
	int lower = invertTime(t);
	//~ std::cout << Size << std::endl;
	//~ std::cout << t << " " << lower << std::endl;
	
	if (lower >=Size - 2)
	{
		return ColdGasMass[Size-1];
	}
	
	int upper = lower + 1;
	
	
	double bruch = (1.0*t - 1.0*timeStorage[lower])/(timeStorage[upper] - timeStorage[lower]);
	return (ColdGasMass[upper] - ColdGasMass[lower])*bruch + ColdGasMass[lower];

}
double MassReservoir::HotGas(double t)
{
	int lower = invertTime(t);
	
	if (lower >=Size - 1)
	{
		return HotGasMass[Size-1];
	}
	
	int upper = lower + 1;
	
	double bruch = (1.0*t - 1.0*timeStorage[lower])/(timeStorage[upper] - timeStorage[lower]);
	return (HotGasMass[upper] - HotGasMass[lower])*bruch + HotGasMass[lower];
}

int MassReservoir::invertTime(double t)
{
	double min = timeStorage[0];
	double max = timeStorage[timeStorage.size() - 1];
	

	if ((t < min) || (t > max))
	{
		std::cout << "ERROR: A time " << t << " was passed to the gas mas inverter. This is outside of valid timebase. QUITTING." << std::endl;
		exit(1);
	}
	
	
	
	double deltaT =timeStorage[1] - timeStorage[0];
	
	
	
	int I = (int)(t/deltaT);
	return I;
}


void MassReservoir::Evolve(bool saveFileActive)
{
	double dt = timeStorage[1] - timeStorage[0];
	for (int i = 1; i < Size; ++i)
	{
		double t = timeStorage[i];
		double D = DeathFunction(timeStorage[i-1]);
		double dmStar = PP.nuSFR.Value*ColdGasMass[i-1] - D;
		double dmCold = accretionRate(t) - PP.nuSFR.Value*ColdGasMass[i-1] + PP.nuCool.Value * HotGasMass[i-1] + (1.0 - PP.hotFrac.Value)*D;
		double dmHot = - PP.nuCool.Value * HotGasMass[i-1] + PP.hotFrac.Value*D;
		
		ColdGasMass[i] = ColdGasMass[i-1] + dt*dmCold;
		StarMass[i] = StarMass[i-1] + dt*dmStar;
		HotGasMass[i] = HotGasMass[i-1] + dt*dmHot;
		deathRate[i] = D;
	}	
	
	if (saveFileActive)
	{
		std::ofstream saveFile;
		std::string saveFileName =PP.FILEROOT + "SFRhistory.dat";
		saveFile.open(saveFileName);
		std::vector<std::string> titles = {"Time", "ColdGas", "HotGas", "StarMass", "Death","SFR"};
		int width = 15;
		for (int i = 0; i < titles.size(); ++i)
		{
			saveFile << std::setw(width) << std::left << titles[i];  
		}
		saveFile << "\n";
		for (int j = 0; j < Size; ++j)
		{
			std::vector<double> vals = { timeStorage[j], ColdGasMass[j], HotGasMass[j], StarMass[j], deathRate[j], PP.nuSFR.Value*ColdGasMass[j]};
			for(int i = 0; i < vals.size(); ++i)
			{
				saveFile << std::setw(width) << vals[i];
			}
			saveFile << "\n";
		}
	
		saveFile.close();
	}
}


double MassReservoir::Integrand(double Q,double t,double alpha, double gamma,double lifetime)
{
	double mass = pow((t-Q)/lifetime, -1.0/gamma);
	
	double imf = IMF(mass,alpha);
	double density = ColdGas(Q);
	double jacobian = pow(t - Q,-1.0*(1+gamma)/gamma);
	
	double I = imf*density*jacobian;
	
	if (isnan(I))
	{
		I =0;
	}
	
	return I;
}

double MassReservoir::DeathFunction(double t)
{
	double alpha = PP.AlphaKS;
	double gamma = 2.5;
	double tauTilde = 10;
	double prefactor = pow(tauTilde,1.0/gamma)*PP.nuSFR.Value/gamma;
	int N = 200;
	double h = 0.01;
	double pi = 3.141592654;
	double sum = 0;
	for (int k = -N; k <= N; ++k)
	{
		double xk = tanh(0.5*pi*sinh(k*h));
		
		double coshVal = cosh(0.5*pi*sinh(k*h));
		double wk = 0.5*h*pi*cosh(k*h)/(coshVal*coshVal);
		
		double Qk = t*(xk+1)/2;
		
		double fk = Integrand(Qk,t,alpha,gamma,tauTilde);
		
		if (!isinf(fk))
		{
			sum +=fk*wk;
		}
	}
	sum*=t/2;
	
	if (isnan(sum) || isinf(sum))
	{
		sum = 0;
	}
	
	return sum*prefactor;
}
