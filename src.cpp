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

#include "src/TimeSystems.h"


std::vector<bool> threadActive(1,false);
std::vector<double> TimeVector;

std::vector<std::vector<int>> BigGrid;


ParameterPack RandomiseGalaxy(ParameterPack pp)
{
	
	
}

void SaveGrid(ParameterPack copy)
{
	std::ofstream saveFile;
	std::string saveFileName =copy.FILEROOT + "SuccessGrid.dat";
	saveFile.open(saveFileName);
	int width = 15;
	for (int i = 0; i < copy.tauColls.NSteps; ++i)
	{
		for (int j = 0; j < copy.collFrac.NSteps; ++j)
		{
			saveFile << std::setw(width) << std::left << BigGrid[i][j];
		}
		saveFile << "\n";
	}
	
	saveFile.close();
}

void LaunchProcess(ParameterPack copy, int id)
{
	std::vector<std::vector<int>> grid = std::vector(copy.tauColls.NSteps,std::vector(copy.collFrac.NSteps,0));
//	ISMIterator iterator = ISMIterator(copy,TimeVector, &grid);
	
	for (int i = 0; i < copy.tauColls.NSteps; ++i)
	{
		for (int j = 0; j < copy.collFrac.NSteps; ++j)
		{
			BigGrid[i][j] += grid[i][j]; 
		}
	}

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
	int NLoops = pp.CountThreadLoops();
	
	//preparing iterator values
	ParameterPack copy = pp;
		
	//initialising threads
	int currentThread = 0;
	std::vector<std::thread> persei(pp.NThreads);
	threadActive.resize(pp.NThreads);
	for (int i = 0; i < pp.NThreads; ++i)
	{
		threadActive[i] = false;
	}
	bool noThreadAssigned = false;
	int i = 0;
	while (copy.InitialisedCorrectly)
	{	
		//perform the recursive search through the provided parameters, and save them to a ParameterPack object. 
		copy  = RandomiseGalaxy(pp);
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
		if(i%10 == 0)
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

	
	if (pp.Mode == 0)
	{
		PathAnnulus A = PathAnnulus(pp);
		A.Evolve();
	}
	else
	{
		pp.LosePresets();		
		BigGrid = std::vector(pp.tauColls.NSteps,std::vector(pp.collFrac.NSteps,0));
		IterationMode(pp);
		SaveGrid(pp);
	}

}


