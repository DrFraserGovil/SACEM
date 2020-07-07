#pragma once
#include "Process.h"


class Accretion : public Process
{
	public:
		Accretion(){};
		Accretion(ParameterPack pp) : Process{pp,-0.001,-0.002}
		{
			
		};
		
		double Mass(double t);
};
class StarFormation : public Process
{
	public:
		StarFormation(){};
		StarFormation(ParameterPack pp, double nuCool, double nuDelay) : Process {pp,nuCool, nuDelay}
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
		CCSN(ParameterPack pp) : Process {pp,pp.CCSNCool.Value, -0.1}
		{
			hotFrac = pp.CCSNHotFrac.Value;
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
		Collapsar(ParameterPack pp) : Process {pp,pp.CollapsarCool.Value, -0.1}
		{
			hotFrac = pp.CollapsarHotFrac.Value;
			
			
			T = pp.tauColls.Value;
			Delta = pp.collWidth.Value;
			
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
		Decayer(ParameterPack pp, double nuCool, double nuDelay, double tau, double decayHotFrac) : Process {pp, nuCool, nuDelay}
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
