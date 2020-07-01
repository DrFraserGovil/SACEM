#include "Vectors.h"


ComplexVector::ComplexVector()
{
	Real.resize(1);
	Real[0] = 0;
	Imaginary.resize(1);
	Imaginary[0] = 0;
}

ComplexVector::ComplexVector(std::vector<double> a, std::vector<double> b)
{
	Real = a;
	Imaginary = b;
}

ComplexVector::ComplexVector(std::vector<double> a, bool realOnly)
{
	if (realOnly)
	{
		Real = a;
		Imaginary = std::vector(Real.size(), 0.0);
	}
	else
	{
		Imaginary = a;
		Real = std::vector(Imaginary.size(), 0.0);
	}
}
ComplexVector::ComplexVector(int n)
{
	Resize(n);
}


void ComplexVector::Resize(int n)
{
	Real.resize(n,0.0);
	Imaginary.resize(n,0.0);
}

int ComplexVector::size()
{
	int n1 = Real.size();
	int n2 = Imaginary.size();
	if (n1 != n2)
	{
		std::cout << "Error: Complex vector has differing internal dimensions" << std::endl;
		exit(99);
	}
	return n1;
}

void ComplexVector::Print()
{
	for (int i = 0; i < size(); ++i)
	{
		std::cout << Real[i] << " + i" << Imaginary[i] << "\n";
	}
	std::cout<<std::endl;
}
double Dot(std::vector<double> a, std::vector<double> b)
{
	int n1 = a.size();
	int n2 = b.size();
	if (n1 != n2)
	{
		std::cout << n1 << "  " << n2 <<std::endl;
		std::cout <<  "The Dot product only accepts vectors of a uniform size" << std::endl;
		exit(5);
	}
	
	double sum = 0;
	for (int i = 0; i < n1; ++i)
	{
		sum+=a[i] * b[i];
	}
	
	return sum;
}

std::vector<double> Hadamard(std::vector<double> a, std::vector<double> b)
{
	int n1 = a.size();
	int n2 = b.size();
	
	if (n1 != n2)
	{
		std::cout <<  "The Hadamard product only accepts vectors of a uniform size" << std::endl;
		exit(5);
	}
	
	std::vector<double> newVec = std::vector(n1,0.0);
	
	for (int i = 0; i < n1; ++i)
	{
		newVec[i] = a[i]*b[i];
	}
	return newVec;
}

std::vector<double> Add(std::vector<double> a, std::vector<double> b)
{
	int n1 = a.size();
	int n2 = b.size();
	
	if (n1 != n2)
	{
		std::cout << n1 << "  " << n2 <<std::endl;
		std::cout <<  "Elementwise Vector Addition only accepts vectors of a uniform size" << std::endl;
		exit(5);
	}
	
	std::vector<double> newVec = std::vector(n1,0.0);
	
	for (int i = 0; i < n1; ++i)
	{
		newVec[i] = a[i]+b[i];
	}
	return newVec;
}
std::vector<double> Scale(std::vector<double> a, double b)
{
	int n1 = a.size();

	
	std::vector<double> newVec = std::vector(n1,0.0);
	
	for (int i = 0; i < n1; ++i)
	{
		newVec[i] = a[i]*b;
	}
	return newVec;
}

std::vector<double> operator+(std::vector<double> a, std::vector<double> b)
{
	return Add(a,b);
}
std::vector<double> operator*(std::vector<double> a, double b)
{
	return Scale(a,b);
}
std::vector<double> operator*(double a, std::vector<double> b)
{
	return Scale(b,a);
}
std::vector<double> operator-(std::vector<double> a, std::vector<double> b)
{
	return Add(a,Scale(b,-1.0));
}


ComplexVector operator+(ComplexVector a, ComplexVector b)
{
	return ComplexVector(a.Real + b.Real, a.Imaginary + b.Imaginary);
}
ComplexVector operator-(ComplexVector a, ComplexVector b)
{
	return ComplexVector(a.Real - b.Real, a.Imaginary - b.Imaginary);
}
ComplexVector operator*(double a, ComplexVector b)
{
	return ComplexVector(a * b.Real, a *  b.Imaginary);
}
ComplexVector operator*(ComplexVector b, double a)
{
	return ComplexVector(a * b.Real, a *  b.Imaginary);
}

//~ ComplexVector operator*(ComplexVector a, double b, double c)
//~ {
	//~ return ComplexVector(b * a.Real - c*a.Imaginary, b *  a.Imaginary + c*a.Real);
//~ }

//~ ComplexVector operator*(double b, double c, ComplexVector a)
//~ {
	//~ return ComplexVector(b * a.Real - c*a.Imaginary, b *  a.Imaginary + c*a.Real);
//~ }

std::vector<double> Dot(ComplexVector a, ComplexVector b)
{
	std::vector<double> output = std::vector(2,0.0);
	
	output[0] = Dot(a.Real, b.Real) - Dot(a.Imaginary,b.Imaginary);
	output[1] = Dot(a.Real, b.Imaginary) + Dot(a.Imaginary, b.Real);
	return output;
}

ComplexVector Hadamard(ComplexVector a, ComplexVector b)
{
	std::vector<double> re = Hadamard(a.Real, b.Real) - Hadamard(a.Imaginary, b.Imaginary);
	std::vector<double> im = Hadamard(a.Real, b.Imaginary) + Hadamard(a.Imaginary, b.Real);
	
	return ComplexVector(re,im);
}
double RealDot(ComplexVector a, ComplexVector b)
{
	return Dot(a.Real, b.Real) - Dot(a.Imaginary,b.Imaginary);	
}
