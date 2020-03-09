#pragma once
#include <vector>
#include <iostream>
#include "ParameterPack.h"
#include "MassReservoir.h"

class ISMIterator
{
	public:
	
		ISMIterator(ParameterPack pp,std::vector<double> & timeVector);
		MassReservoir ISM;
	
		void Iterate(std::vector<double> & TimeVector);
	private:
		ParameterPack PP;
		
		
	
};

class EFHIterator
{
	public:
		EFHIterator(ParameterPack pp, ISMIterator * parent);
	
	private:
		
		void LoadVectors();
		
		void Iterate();
		
		
		double EInf;
		std::vector<double> E0;
		std::vector<std::vector<std::vector<double>>> F0_List;
		std::vector<std::vector<double>> FInf_List;
		std::vector<std::vector<double>> HInf_NSM_List;
		std::vector<std::vector<std::vector<double>>> H0_NSM_List;
		std::vector<std::vector<double>> HInf_SNIa_List;
		
		ParameterPack PP;
		ISMIterator * Parent;
		
		double Evalues(double t, double prevVal, double prevT);
		double Fvalues(double t, double tau, double w, double prevVal, double prevT);
		double Gvalues(double t, double tau, double nu);
		double Hvalues(double t, double tau, double nu, double prevVal, double prevT);
		double tauSN(int n);
};

class CalibrationIterator
{
	public:
		CalibrationIterator(ParameterPack pp, EFHIterator * parent);
		
	private:
		
		void Calibrate();
		bool Evaluate();
		
		
		//normalised constants
		double alpha;
		double beta;
		double gamma;
		double delta;
		double epsilon;
		double eta;
	
};





//~ class IterationAnnulus
//~ {
	//~ public:
		//~ IterationAnnulus(ParameterPack pp, std::vector<std::vector<double>> FGH, std::vector<std:vector<double>> ResultsVector);
	
		//~ void Evaluate();
		
	//~ private:
		//~ ParameterPack PP;
		
		//~ std::vector<std::vector<double>> * ResultsVector;
		
		//~ std::vector<double> Quick;
		//~ std::vector<double> Slow;
		
		//~ std::vector<double> Europium;
		//~ std::vector<double> Iron;
		//~ std::vector<double> Magnesium;
		
		
		
		
		//~ void Calibrate();
//~ };
