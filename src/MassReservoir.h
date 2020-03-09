#pragma once
#include <vector>
#include <iostream>
#include "ParameterPack.h"

class MassReservoir
{
	public:
				
		MassReservoir();
		
		MassReservoir(std::vector<double> time,ParameterPack pp, bool saveFileActive);
		
		double Stars(double t);
		double ColdGas(double t);
		double HotGas(double t);
		int Size;
	private:
		std::vector<double> StarMass;
		std::vector<double> ColdGasMass;
		std::vector<double> HotGasMass;
		std::vector<double> timeStorage;
		std::vector<double> deathRate;
		
		ParameterPack PP;
		int invertTime(double t);
		
		
		
		void Evolve(bool saveFileActive);
		
		double DeathFunction(double t);
		double Integrand(double Q,double t, double alpha, double gamma,double lifetime);
		double accretionRate(double t);
};
