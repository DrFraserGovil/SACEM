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
#include "src/Annulus.h"
#include "src/Parameters.h"
#include "src/TimeSystems.h"


std::vector<bool> threadActive(1,false);
std::vector<double> TimeVector;

std::vector<std::vector<int>> BigGrid;


ParameterPack RandomiseGalaxy(ParameterPack pp)
{
	ParameterPack copy = pp;
	
	copy.ScrambleAll();
	
	return copy;
}

void SaveGrid(ParameterPack copy)
{
	std::ofstream saveFile;
	std::string saveFileName =copy.FILEROOT + "/SuccessGrid.dat";
	saveFile.open(saveFileName);
	int width = 15;
	for (int i = 0; i < copy.collFrac.NSteps; ++i)
	{
		for (int j = 0; j < copy.tauColls.NSteps; ++j)
		{
			if (j > 0)
			{
				saveFile << ", ";
			}
			saveFile << std::setw(width) << std::left << BigGrid[i][j];
			
		}
		saveFile << "\n";
	}
	
	saveFile.close();
}



void LaunchProcess(ParameterPack state, std::vector<std::vector<int>> * miniGrid, int loopNumber, int id)
{	
	bool folderMade = false;
	for (int i = 0; i < state.collFrac.NSteps; ++i)
	{
		state.collFrac.IterateValue(i);
		for (int j = 0; j < state.tauColls.NSteps; ++j)
		{
			state.tauColls.IterateValue(j);
			
			ParameterPack safeState = state;
			safeState.ValueChecks();
			//create + evaluate an annulus using the given parameters
			Annulus A = Annulus(safeState);
			safeState.WasSuccessful = A.QuickAnalysis();
			
		
			if (state.WasSuccessful)
			{
				safeState.MeetsValueLimits = A.ValueAnalysis();
				
				if (safeState.MeetsValueLimits)
				{
					
					if (safeState.collFrac.Value > 0.7 & safeState.tauColls.Value < 4)
					{
						if (folderMade == false)
						{
							std::string commandStr= "mkdir FullPaths/Iteration" +std::to_string(loopNumber);
							const char * command = commandStr.c_str();
							
							system(command); 
							
							folderMade = true;
						}
						
						ostringstream simFileName;
						simFileName << "FullPaths/Iteration " << loopNumber << "/Grid_" << i << "_" << j; 
						
						
						A.SaveAnnulus(simFileName.str());
						
						ostringstream stateSave;
						stateSave << "ParamSaves/Iteration" << loopNumber << "/Param_"  << i << "_" << j;
						safeState.SaveState(stateSave.str());
					}
					
					++miniGrid[0][i][j];
				}
				
			}
			
			//~ bool isBeingSaved = evaluateSaveConditions(state,i,j);
			
			//~ if (isBeingSaved)
			//~ {
				//~ ostringstream simFileName;
				//~ simFileName << "FullPaths/Grid_" << loopNumber << "_" << i << "_" << j << "_" << state.WasSuccessful << ".dat"; 
				//~ A.Evolve();
				//~ A.SaveAnnulus(simFileName.str());
				
				//~ ostringstream stateSave;
				//~ stateSave << "ParamSaves/Param_" << loopNumber << "_"  << i << "_" << j << "_" << state.WasSuccessful << ".dat";
				//~ state.PrintState(stateSave.str());
			//~ }
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
	

	int NLoops = pp.NRandomGalaxies;
	
	//preparing iterator values
	ParameterPack copy = pp;
		
	//initialising threads
	int currentThread = 0;
	
	std::vector<std::vector<std::vector<int>>> miniGrids(pp.NThreads,std::vector(copy.collFrac.NSteps, std::vector<int>(copy.tauColls.NSteps,0)));
	std::vector<std::thread> persei(pp.NThreads);
	threadActive.resize(pp.NThreads);
	for (int i = 0; i < pp.NThreads; ++i)
	{
		threadActive[i] = false;
	}
	bool noThreadAssigned = false;


	

	for (int k = 0; k < NLoops; ++k)
	{	
		//std::cout << "Randomiser " << k << std::endl;
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
						for (int i = 0; i < copy.collFrac.NSteps; ++i)
						{
							for (int k = 0; k < copy.tauColls.NSteps; ++k)
							{
								BigGrid[i][k] +=miniGrids[j][i][k];
								miniGrids[j][i][k] = 0;
							}
						}
						
					}
					currentThread = j;
					j = pp.NThreads;
					noThreadAssigned = false;
				}
			}
			
		}

		threadActive[currentThread] = true;
		persei[currentThread] = std::thread(LaunchProcess,copy,&miniGrids[currentThread],k,currentThread);
		noThreadAssigned = true;	
		
		int nPrints = pp.NThreads + 1;
		if(k%nPrints == 0)
		{
			printTimeSince(start,k,NLoops);
		}
		
		
	}
	
	for (int j = 0; j < pp.NThreads; ++j)
	{
		if (persei[j].joinable())
		{
			persei[j].join();
			for (int i = 0; i < copy.collFrac.NSteps; ++i)
						{
							copy.collFrac.IterateValue(i);
							for (int k = 0; k< copy.tauColls.NSteps; ++k)
							{
								BigGrid[i][k] +=miniGrids[j][i][k];
							}
						}
		}
	}
	std::cout << "All threads closed" << std::endl;
	
}

void VerificationMode(ParameterPack pp)
{
	int N = 100;
	ParameterPack currentState;
	for (int i = 0; i < N; ++i)
	{
		currentState = pp;//RandomiseGalaxy(pp);
		std::string name = "verificationRun_" + std::to_string(i);
		Annulus A = Annulus(currentState);
		A.Evolve();
		A.SaveAnnulus(name);
		
		bool success = A.ValueAnalysis();
		currentState.WasSuccessful = success;
		std::cout << currentState.PrintState() <<std::endl;
		
		std::string command = "gnuplot -p -e \"filename = '" + pp.FILEROOT + name + ".dat'\" plotter.gp";
		const char *com = command.c_str();
		system(com);
		
		//~ command = "gedit " + pp.FILEROOT + name + "_state.dat";
		//~ const char *com2 = command.c_str();
		//~ system(com2);
		
		
		std::cin.get();
		std::cin.clear();
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
		
		Annulus A = Annulus(pp);
		A.Evolve();
		A.SaveAnnulus("SingleEvaluation");
		
	}
	
	if (pp.Mode == 1)
	{
		BigGrid = std::vector(pp.collFrac.NSteps,std::vector(pp.tauColls.NSteps,0));
		IterationMode(pp);
		SaveGrid(pp);
	}

	if (pp.Mode == 2)
	{
		VerificationMode(pp);
	}


}


