#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <math.h>
#include <iomanip>

#include <thread>
#include <future>

#include "src/ParameterPack.h"
#include "src/commandLineParser.h"
#include "src/MassReservoir.h"
#include "src/Annulus.h"
#include "src/IteratorClasses.h"
#include "src/TimeSystems.h"


std::vector<bool> threadActive(1,false);
std::vector<double> TimeVector;

ParameterPack Dive(ParameterPack pp, std::vector<double> &minValues, std::vector<double> &maxValues, std::vector<int> &steps, std::vector<int> &pos)
{
	bool found = false;
	int N = pos.size();
	for (int j= 0; j < N; ++j)
	{
		pos[j] = pos[j] + 1;
		
		if (pos[j] >= steps[j] )
		{
			pos[j] = 0;
			found = false;
		}
		else
		{
			j = N;
			found = true;
		}
	}
	ParameterPack copy = pp;
	copy.InitialisedCorrectly = found;
	
	std::vector<double * > currentValues =  {&copy.galaxyM0, &copy.galaxyM1, &copy.galaxyM2, &copy.galaxyB1, &copy.galaxyB2, &copy.galaxyScaleLength, &copy.nuSFR, &copy.nuCool,&copy.alphaKS, &copy.hotFrac};

	if (found)
	{

		for (int j = 0; j < N; ++j)
		{
			double bruch = (maxValues[j] - minValues[j])/(N - 1);
		
			currentValues[j][0] = minValues[j] + pos[j] * bruch;
		}
	}

	
	return copy;
}


void LaunchProcess(ParameterPack copy, int id)
{
	ISMIterator iterator = ISMIterator(copy,TimeVector);
	threadActive[id] = false;
}

void generateTimeVector(ParameterPack PP)
{
	TimeVector = std::vector(0,0.0);
	double tMax = PP.tauInf*1.02;
	double deltaT = PP.timeStep;
	double t = 0;
	while (t <= tMax)
	{
		TimeVector.push_back(t);
		t+=deltaT;
	}
}

void IterationMode(ParameterPack pp)
{
	auto start = std::chrono::high_resolution_clock::now();
	
	
	generateTimeVector(pp);
	
	std::vector<double> minValues = {pp.galaxyM0_Min, pp.galaxyM1_Min, pp.galaxyM2_Min, pp.galaxyB1_Min, pp.galaxyB2_Min, pp.galaxyScaleLength_Min, pp.nuSFR_Min, pp.nuCool_Min,pp.alphaKS_Min, pp.hotFrac_Min};
	std::vector<double> maxValues = {pp.galaxyM0_Max, pp.galaxyM1_Max, pp.galaxyM2_Max, pp.galaxyB1_Max, pp.galaxyB2_Max, pp.galaxyScaleLength_Max, pp.nuSFR_Max, pp.nuCool_Max,pp.alphaKS_Max, pp.hotFrac_Max};
	
	
	std::vector<int>    valueSteps = {pp.galaxyM0_N, pp.galaxyM1_N, pp.galaxyM1_N, pp.galaxyB1_N, pp.galaxyB2_N, pp.galaxyScaleLength_N, pp.nuSFR_N, pp.nuCool_N,pp.alphaKS_N, pp.hotFrac_N};
	std::vector<int> pos = std::vector(minValues.size(),0);
	
	int NLoops = 1;
	for (int j = 0; j < minValues.size(); ++j)
	{
		NLoops *= valueSteps[j];
	}
	
	int i = 0;
	ParameterPack copy = pp;
	
	int currentThread = 0;
	std::vector<std::thread> persei(pp.NThreads);
	threadActive.resize(pp.NThreads);
	for (int i = 0; i < pp.NThreads; ++i)
	{
		threadActive[i] = false;
	}
	bool noThreadAssigned = false;
	
	while (copy.InitialisedCorrectly)
	{	
		copy  = Dive(pp,minValues, maxValues, valueSteps, pos);
		while (noThreadAssigned == true)
		{
			for (int j = 0; j < pp.NThreads; ++j)
			{

				if (threadActive[j]==false)
				{
					if (persei[j].joinable())
					{
						persei[j].join();
						
					}
					currentThread = j;
					j = pp.NThreads;
					noThreadAssigned = false;
				}
			}
			
		}
			
		threadActive[currentThread] = true;
		persei[currentThread] = std::thread(LaunchProcess,copy,currentThread);
		noThreadAssigned = true;	
		
		++i;
		if(i%100 == 0)
		{
			printTimeSince(start,i,NLoops);
		}
		
		
	}
	
	
	for (int j = 0; j < pp.NThreads; ++j)
	{
		if (persei[j].joinable())
		{
			persei[j].join();
		}
	}
	
}


int main(int argc, char** argv)
{
	
	std::cout << "S-ChEM Routine Initialised" << std::endl;
	
	//Parse command line
	ParameterPack pp = parseCommandLine(argc, argv);
	if (!(pp.InitialisedCorrectly))
	{
		std::cout << "Parse failed. Quitting" << std::endl;
		return 1;
	}
	std::cout << "Parsed command line arguments" << std::endl;
	
	double radius = 7;
	double width = 0.2;
	
	pp.UpdateRadius(radius,width);
	//Annulus A = Annulus(radius,width,tauInf);
	

	//~A.Evolve();
	
	if (pp.Mode == 0)
	{
		PathAnnulus A = PathAnnulus(pp);
		A.Evolve();
	}
	else
	{
		IterationMode(pp);
	}

}


