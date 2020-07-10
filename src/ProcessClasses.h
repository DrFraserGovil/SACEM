#pragma once
#include "Process.h"


class GalaxyMass : public Process
{
	public:
		GalaxyMass(){};
		ComplexVector StarVector;
		double StarConstant;
		GalaxyMass(ParameterPack * pp) : Process{pp,-0.001,-0.002}
		{
			ComplexVector newVec = ComplexVector(Powers.size());

			for (int i = 0; i < Powers.size(); ++i)
			{
	
				double ar = PP->stellarDeathParameter.Value;
				double ai = 0;
				double br = Powers.Real[i];
				double bi = Powers.Imaginary[i];
				
				double denom = (ar - br)*(ar - br) + (ai - bi)*(ai - bi);
				
				newVec.Real[i] = (ar - br)/denom;
				newVec.Imaginary[i] = -(ai - bi)/denom;
			}
			
			StarVector = Hadamard(newVec, SFRVector);
			
			StarConstant = - RealDot(StarVector,E(0));
		};
		
		double TotalMass(double t);
		double StellarMass(double t);
		double ColdGasMass(double t);
};
class StarFormation : public Process
{
	public:
		StarFormation(){};
		StarFormation(ParameterPack * pp, double nuCool, double nuDelay) : Process {pp,nuCool, nuDelay}
		{
			//all automated?
		};
	
	
		double Rho(double t);
};



class CCSN : public Process
{
	
	private:
		ComplexVector CCSNVector;
		ComplexVector TotalVector;
		void GeneratePrefactors();
	
	public:
		CCSN(){};
		CCSN(ParameterPack * pp) : Process {pp,pp->CoolingFrequency.Value, -0.1}
		{
			hotFrac = pp->CCSNHotFrac.Value;
			GeneratePrefactors();
		};
	
	
		double Count(double t);
		double Total(double t);
};


class Collapsar : public Process
{
	private:
		ComplexVector UntouchedVector;
		ComplexVector IntermediaryVector_Constant;
		ComplexVector IntermediaryVector_Linear;
		ComplexVector FinalVector;
		
		ComplexVector InitialTotal;
		ComplexVector MidConstantTotal;
		ComplexVector MidLinearTotal;
		double TotalConst; 
		
		double T;
		double Delta;
		
		void GeneratePrefactors();
	
	public:
		Collapsar(){};
		Collapsar(ParameterPack * pp) : Process {pp,pp->CoolingFrequency.Value, -0.1}
		{
			hotFrac = pp->CollapsarHotFrac.Value;
			
			
			T = pp->tauColls.Value;
			Delta = pp->collWidth.Value;
			
			GeneratePrefactors();
		};
	
	
		double Count(double t);
		double Total(double t);
		double Initial(double t);
		double Middle(double t);
		double End(double t);
		double TotalMid(double t);
};

class Decayer : public Process
{
	private:
		ComplexVector DecayVector;
		ComplexVector TotalVector;
		
		double Tau;
		double nu;
		
		void GeneratePrefactors();
		
	public:
		Decayer(){};
		Decayer(ParameterPack * pp, double nuCool, double nuDelay, double tau, double decayHotFrac) : Process {pp, nuCool, nuDelay}
		{
			hotFrac = decayHotFrac;
			
			Tau = tau;
			nu = nuDelay;
			
			GeneratePrefactors();
		}
		
		double Total(double t);
		double Count(double t);
		
		double Yield(double t);
};
