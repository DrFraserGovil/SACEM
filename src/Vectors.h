#pragma once

#include <math.h>
#include <iomanip>
#include <vector>
#include <iostream>
class ComplexVector
{
	public:
	
		std::vector<double> Real;
		std::vector<double> Imaginary;
	
		ComplexVector();
		ComplexVector(std::vector<double> a, std::vector<double> b);
		ComplexVector(std::vector<double> a, bool realOnly);
		ComplexVector(int n);
		
		void Resize(int n);
		int size();
		void Print();
};



double Dot(std::vector<double> a, std::vector<double> b);
std::vector<double> Hadamard(std::vector<double> a, std::vector<double> b);
std::vector<double> Add(std::vector<double> a, std::vector<double> b);
std::vector<double> Scale(std::vector<double> a, double b);

std::vector<double> operator+(std::vector<double> a, std::vector<double> b);
std::vector<double> operator*(std::vector<double> a, double b);
std::vector<double> operator*(double a, std::vector<double> b);
std::vector<double> operator-(std::vector<double> a, std::vector<double> b);

ComplexVector operator+(ComplexVector a, ComplexVector b);
ComplexVector operator-(ComplexVector a, ComplexVector b);
ComplexVector operator*(double a, ComplexVector b);
ComplexVector operator*(ComplexVector a, double b);

std::vector<double> Dot(ComplexVector a, ComplexVector b);
double RealDot(ComplexVector a, ComplexVector b);
ComplexVector Hadamard(ComplexVector a, ComplexVector b);
//~ ComplexVector operator*(double a, double b, ComplexVector c);
//~ ComplexVector operator*(ComplexVector a, double b, double c);
