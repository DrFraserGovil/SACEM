#pragma once

#include <iostream>
#include <math.h>
#include <iomanip>
#include <vector>
#include "ParameterPack.h"
#include "Vectors.h"
class Process
{
	public:
		ComplexVector Powers;
		double hotFrac;
		ParameterPack * PP;
		double Omega;
		
		int ConstantID;
		int WPlusID;
		int WMinusID;
		int FormationID;
		int CoolingID;
		int DelayID;
		
		double MBar;
		
		ComplexVector SFRVector;
		bool ComplexDomain; 
	
	
		Process(){};
		Process(ParameterPack * PP, double nuCool, double nuDelay);
		
		ComplexVector E(double t);
	
		ComplexVector Integrator(int freqID);
		ComplexVector CCSNStyleVector();
	
		ComplexVector GenericIntegrator(ComplexVector input, int id);
		
		ComplexVector F(double t);
		
		ComplexVector TotalIntegrator();
};


