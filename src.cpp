#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <math.h>
#include <iomanip>
#include "src/defaultValues.h"
#include "src/commandLineParser.h"



double GalaxyMass(double t)
{
	double initialMass = galaxyM0;
	double thickDisk = galaxyM1 * ( 1.0- exp(-t/galaxyB1));
	double thinDisk = galaxyM2 * (1.0 - exp(-t/galaxyB2));
	
	return initialMass + thickDisk + thinDisk;
}


class Annulus
{
	double Radius;
	double Width;
	double Age;
	
	double MFe;
	double MEu;
	double MMg;
	
	
	double alpha;
	double beta;
	double gamma;
	double delta;
	double epsilon;
	double eta;
	
	std::vector<double> EVec;
	std::vector<double> FVec;
	std::vector<double> GVec;
	std::vector<double> HVecNSM;
	std::vector<double> HVecSNIa;
	
	private:
		inline const double radialMassCorrection()
		{
			return Radius*Width/(galaxyScaleLength*galaxyScaleLength)*exp(-Radius/galaxyScaleLength);
		}
	
		inline double RingMass(double t)
		{
			double thickDisk = galaxyM1 * (1 - exp(-t/galaxyB1));
			double thinDisk = galaxyM2 * (1 - exp(-t/galaxyB2));
			
			return (galaxyM0 + thickDisk + thinDisk)*radialMassCorrection();
		}
	
		double SFR(double t)
		{
				double basePart = galaxyM0*exp(-nuSFR*t);
				double thickDiskBoost = galaxyM1/(nuSFR*galaxyB1-1) * (exp(-t/galaxyB1) - exp(-nuSFR*t));
				double thinDiskBoost = galaxyM2/(nuSFR*galaxyB2-1) * (exp(-t/galaxyB2) - exp(-nuSFR*t));
				
				return nuSFR*(basePart+thickDiskBoost+thinDiskBoost);
		}
		
		
		double E(double t)
		{
			//simple first-order integrator
			int Nsteps = 10000;
			double deltaT = t/Nsteps;
			double sum = 0;
			for (double tPrime = 0; tPrime < t; tPrime +=deltaT)
			{
				sum += SFR(tPrime)*deltaT;
			}
			return sum;
		}
		
		inline double decay(double t, double tau, double width)
		{
			if (t < (tau - width/2))
			{
				return 1;
			}
			if (t > (tau + width/2))
			{
				return 0;
			}
			return (tau - t)/width + 0.5;
		}
		
		double F(double t, double tau,double width)
		{
			//simple first-order integrator
			int Nsteps = 10000;
			double deltaT = t/Nsteps;
			double sum = 0;
			for (double tPrime = 0; tPrime < t; tPrime +=deltaT)
			{
				sum += SFR(tPrime)*decay(tPrime,tau,width)*deltaT;
			}
			return sum;
		}
		
		double G(double t,double tau, double nu)
		{
			if (t <= tau)
			{
				return 0;
			}
			//simple first-order integrator
			int Nsteps = 1000;
			double deltaT = t/Nsteps;
			double sum = 0;
			for (double tPrime = tau; tPrime < t; tPrime +=deltaT)
			{
				sum += SFR(t-tPrime)*exp(-nu*(tPrime - tau))*deltaT;
			}
			return sum;
		}
		
		double H(double t, double tau, double nu)
		{
			//simple first-order integrator
			int Nsteps = 1000;
			double deltaT = t/Nsteps;
			double sum = 0;
			for (double tPrime = tau; tPrime < t; tPrime +=deltaT)
			{
				sum += G(tPrime,tau,nu)*deltaT;
			}
			return sum;
		}
		
		double CalibrateYields()
		{
			
			
			double E0 = E(tauSNIa);
			double EInf = E(tauInf);
			double F0 = F(tauSNIa,tauColls,collWidth);
			double FInf = F(tauInf,tauColls,collWidth);
			double H0_NSM = H(tauSNIa,tauNSM,nuNSM);
			double HInf_NSM = H(tauInf,tauNSM,nuNSM);
			double HInf_SNIa = H(tauInf,tauSNIa,nuSNIa);
			
			alpha = RingMass(tauSNIa)/E0 * pow(10,FeH_SN);
			std::cout << "Fe CCSN Yields Calibrated " << alpha << std::endl;
			
			beta = alpha*EInf/HInf_SNIa* (pow(10,MgFe_SN-MgFe_Sat)-1);
			std::cout << "Fe SNIa Yields Calibrated" << beta << std::endl;
			
			
			double gammaFactor = sProcFrac*E0/EInf + collFrac*F0/FInf + (1 - collFrac - sProcFrac)*H0_NSM/HInf_NSM;
			gamma = alpha * sProcFrac* E0/EInf/gammaFactor*pow(10,MgFe_SN + EuMg_SN);
			std::cout << "Eu CCSN Yields Calibrated" << std::endl;
			
			delta = collFrac/sProcFrac * EInf/FInf * gamma;
			std::cout << "Eu Collapsar Yields Calibrated" << std::endl;
			
			epsilon = (1.0 - collFrac - sProcFrac)/sProcFrac * EInf/HInf_NSM * gamma;
			std::cout << "Eu NSM Yields Calibrated" << std::endl;
			
			eta  = alpha * pow(10,MgFe_SN);
			std::cout << "Mg CCSN Yields Calibrated" << std::endl;
		}
	
	
		double PopulateVectors()
		{
			double deltaT = 1.0/IntegrationSteps*Age;
			for (int i = 1; i <= IntegrationSteps;++i)
			{
				double t = (float)i/IntegrationSteps*Age;
				double sfr = SFR(t);
				EVec[i] = EVec[i-1] + sfr*deltaT;
				FVec[i] = FVec[i-1] + sfr*decay(t,tauColls,collWidth)*deltaT;
				
				if (t > tauNSM)
				{
					double sum = 0;
					double dtP = (t-tauNSM)/(IntegrationSteps - 1);
					for (int j = 0; j < IntegrationSteps; ++j)
					{
						double tau = (t-tauNSM)/(IntegrationSteps - 1)*j + tauNSM;
						sum += SFR(t - tau)*exp(-nuNSM*(tau - tauNSM))*dtP;
					}
					HVecNSM[i] = HVecNSM[i-1] + sum*deltaT;
				}
				if (t > tauSNIa)
				{
					double sum = 0;
					double dtP = (t-tauSNIa)/(IntegrationSteps - 1);
					for (int j = 0; j < IntegrationSteps; ++j)
					{
						double tau = (t-tauSNIa)/(IntegrationSteps - 1)*j + tauSNIa;
						sum += SFR(t - tau)*exp(-nuSNIa*(tau - tauSNIa))*dtP;
					}
					HVecSNIa[i] = HVecSNIa[i-1] + sum*deltaT;
				}
			
			}
			
		}
		
		double TimeStep(double t)
		{
			int index = (int)(IntegrationSteps*t/Age);
			
			MFe = alpha*EVec[index] + beta*HVecSNIa[index];
			
			MEu = gamma*EVec[index] + delta*FVec[index];
			if (epsilon > 0)
			{
				MEu += epsilon * HVecNSM[index];
			}
			
			MMg = eta * EVec[index];
		}
	
	public: 
		Annulus(double r,double deltaR,double tMax)
		{
			std::cout << "An annulus has just been initialised..." << std::endl;
			Age = tMax;
			Radius = r;
			Width = deltaR;
			MFe = 0;
			MEu = 0;
			MMg = 0;			
			
			EVec=std::vector(IntegrationSteps+1,0.0);
			FVec=std::vector(IntegrationSteps+1,0.0);
			GVec=std::vector(IntegrationSteps+1,0.0);
			HVecNSM=std::vector(IntegrationSteps+1,0.0);
			HVecSNIa=std::vector(IntegrationSteps+1,0.0);
			
			PopulateVectors();
			
			CalibrateYields();
			
			
		}
	
		void Evolve()
		{
			int NSteps = 100;
			double sharpness = 2;
			
			std::ofstream saveFile;
			std::string saveFileName =FILEROOT + ".dat";
			saveFile.open(saveFileName);
			saveFile << std::setw(15);
	
	
			for (int i = 0; i <= NSteps; ++i)
			{
				double t = pow((float)i/NSteps,sharpness)*Age;
				
				TimeStep(t);
				double MH = RingMass(t)*0.7;
				std::vector<double> saveValues = {t,SFR(t), RingMass(t),log10(MFe/MH),log10(MEu/MH),log10(MMg/MH)};
				
				for (int j = 0; j < saveValues.size(); ++j)
				{
					saveFile << std::setw(10) << std::left << saveValues[j] << "\t";
				}
				
				saveFile << "\n";
			} 
			
			saveFile.close();
		}
	
};




int main(int argc, char** argv)
{
	
	std::cout << "S-ChEM Routine Initialised" << std::endl;
	
	//Parse command line
	bool parseSucceeded = parseCommandLine(argc, argv);
	if (!(parseSucceeded))
	{
		std::cout << "Parse failed. Quitting" << std::endl;
		return 1;
	}
	
	
	double radius = 7;
	double width = 0.2;
	
	Annulus A = Annulus(radius,width,14);
	
	A.Evolve();
}
