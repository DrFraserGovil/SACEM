#pragma once
#include <vector>
#include <iostream>
#include <iomanip>
#include "ParameterPack.h"
#include "MassReservoir.h"

class ISMIterator
{
	public:
	
		ISMIterator(ParameterPack pp,std::vector<double> & timeVector, std::vector<std::vector<int>> * grid);
		MassReservoir ISM;
	
		std::vector<std::vector<int>> * SuccessCount;
	
		void Iterate(std::vector<double> & TimeVector);
	private:
		ParameterPack PP;
		
		
	
};

class EFHIterator
{
	public:
		EFHIterator(ParameterPack pp, ISMIterator * parent);
		ISMIterator * Parent;
	private:
		
		void LoadVectors();
		
		void Iterate();
		
		
		double EInf;
		std::vector<double> E0_List;
		std::vector<std::vector<std::vector<double>>> F0_List;
		std::vector<std::vector<double>> FInf_List;
		std::vector<std::vector<double>> HInf_NSM_List;
		std::vector<std::vector<std::vector<double>>> H0_NSM_List;
		std::vector<std::vector<double>> HInf_SNIa_List;
		
		ParameterPack PP;
		
		
		double Evalues(double t, double prevVal, double prevT);
		double Fvalues(double t, double tau, double w, double prevVal, double prevT);
		double Gvalues(double t, double tau, double nu);
		double Hvalues(double t, double tau, double nu, double prevVal, double prevT);
		
};

class CalibrationIterator
{
	public:
		CalibrationIterator(ParameterPack pp, double e0, double eInf, double f0, double fInf, double h0_nsm, double hInf_nsm, double hInf_SN, EFHIterator* parent);
		void Iterate();
		
		
	private:
		
		ParameterPack PP;
		
		
		
		void Calibrate(double mass);
		void Evaluate();
		
		//given values
		double E0;
		double EInf;
		double F0;
		double FInf;
		double H0_NSM;
		double HInf_NSM;
		double HInf_SN;
		
		//normalised constants
		double alpha;
		double beta;
		double gamma;
		double delta;
		double epsilon;
		double eta;
	
		EFHIterator * Parent;
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
