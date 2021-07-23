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
std::vector<ParameterPack> ActiveGalaxies;

std::vector<std::string> paramGridNames;
std::vector<std::vector<std::vector<double>>> ParamGrids;

std::vector<ParameterPack> SolvedGalaxies;

ParameterPack RandomiseGalaxy(ParameterPack pp)
{
	ParameterPack copy = pp;
	
	copy.ScrambleAll();
	
	return copy;
}

void InitialiseSaveGrids(ParameterPack copy)
{
	BigGrid = std::vector(copy.collFrac.NSteps,std::vector(copy.tauColls.NSteps,0));
	
	
	//save file headers
	std::vector<std::string> basicHeader = {"Iteration", "collFrac", "tauColl"};
	std::vector<std::string> ppHeader = copy.PrinterHeaders();
	std::vector<std::string> derivedHeader = {"ColdStarRatio","HotStarRatio","StarMassRatio","SFR_start","SFR_mid","SFR_end"};
	std::vector<std::vector<std::string>> allHeaders = {basicHeader, ppHeader,derivedHeader};
	
	int n = ppHeader.size() + derivedHeader.size();

	ParamGrids.resize(n);
	paramGridNames = ppHeader;
	paramGridNames.insert(paramGridNames.end(), derivedHeader.begin(), derivedHeader.end() );
	
	for (int i = 0; i < n; ++i)
	{
		double c = NAN;
		ParamGrids[i] = std::vector(copy.NGrid, std::vector(copy.NGrid, c));
	}
}

void StoreGalaxies(ParameterPack copy)
{

	std::vector<std::string> basicHeader = {"Iteration", "collFrac", "tauColl"};
	std::vector<std::string> ppHeader = copy.PrinterHeaders();
	std::vector<std::string> derivedHeader = {"ColdStarRatio","HotStarRatio","StarMassRatop","SFR_start","SFR_mid","SFR_end"};
	std::vector<std::vector<std::string>> allHeaders = {basicHeader, ppHeader,derivedHeader};
	
	
	
	

	//save file contents
	for (int i = 0; i < SolvedGalaxies.size(); ++i)
	{
		ParameterPack * gal = &SolvedGalaxies[i];
		
		
		std::vector<double> ppVals = gal->PrinterValues();
		
		//galaxyDerivedParams
		
		
		for (int j = 0; j < gal->nSuccess; ++j)
		{
			int fracIndex = gal->SuccessfulFracs[j];
			int tauIndex = gal->SuccessfulTaus[j];
			std::vector<double> basicInfo = {i, gal->collFrac.getStepValue(fracIndex), gal->tauColls.getStepValue(tauIndex)};
			
			
			
			std::vector<std::vector<double>> allVals = {basicInfo, ppVals, gal->derivedParams[j]};
			
			int totalIdx = 0;
			for (int k = 0; k < allHeaders.size(); ++k)
			{
				for (int l = 0;  l < allHeaders[k].size(); ++l)
				{
					//saveFile << std::setw(width) << std::left << allVals[k][l] << "\t";
					if (k > 0)
					{
						
						
						double v = ParamGrids[totalIdx][fracIndex][tauIndex];
						if (isnan(v))
						{
							ParamGrids[totalIdx][fracIndex][tauIndex] = 0;
						}
						ParamGrids[totalIdx][fracIndex][tauIndex] += allVals[k][l];
						++totalIdx;
					}
				}
			}
			//saveFile << "\n";
		}
	}
	
	

}

void SaveGrid(ParameterPack copy,bool trueFinish,int storedModels)
{
	if (trueFinish)
	{
		std::cout << "Simulation finished.\n\tSaving Success-Count Grid" << std::endl;
	}

	//Save number of models that have been searched
	std::ofstream progress;
	std::string progressName =copy.FILEROOT + "Progress.dat";
	progress.open(progressName);
	progress << storedModels << std::endl;
	progress.close();
	
	//load in the saved grids for analysis
	StoreGalaxies(copy);
	
	//save to file
	std::ofstream saveFile;
	std::string saveFileName =copy.FILEROOT + "SuccessCounts.dat";
	saveFile.open(saveFileName);
	
	std::vector<std::ofstream> paramFiles(paramGridNames.size() );
	for (int i = 0; i < paramGridNames.size(); ++i)
	{
		
		std::string name= copy.FILEROOT  + paramGridNames[i] + ".dat";

		paramFiles[i].open(name);
	}
	
	
	
	int width = 8;
	for (int i = 0; i < copy.collFrac.NSteps; ++i)
	{
		for (int j = 0; j < copy.tauColls.NSteps; ++j)
		{
			if (j > 0)
			{
				saveFile << ", ";
			}
			saveFile << std::setw(width) << std::left << BigGrid[i][j] << "\t";
			
			for (int k = 0; k < paramGridNames.size(); ++k)
			{
				if (j > 0)
				{
					paramFiles[k] << ", ";
				}
				paramFiles[k] << std::setw(width) << std::left << ParamGrids[k][i][j]/BigGrid[i][j] << "\t";
			}
		}
		
		//endlines
		saveFile << "\n";
		for (int k = 0; k < paramGridNames.size(); ++k)
		{
			
			paramFiles[k] << "\n";
		}
	}
	
	saveFile.close();
	for (int i=0; i < paramGridNames.size(); ++i)
	{
		paramFiles[i].close();
	}
	
	
}



void LaunchProcess(std::vector<std::vector<int>> * miniGrid, int loopNumber, int id)
{	
	ParameterPack * state = &ActiveGalaxies[id];
	Annulus A = Annulus(state);
	
	for (int i = 0; i < state->collFrac.NSteps; ++i)
	{
		state->collFrac.IterateValue(i);
		for (int j = 0; j < state->tauColls.NSteps; ++j)
		{
			state->tauColls.IterateValue(j);
			
			A.Calibrate();
			
			
			//create + evaluate an annulus using the given parameters
			
			state->WasSuccessful = A.QuickAnalysis();
		
			//std::cout << state->WasSuccessful << std::endl;
		
			if (state->WasSuccessful)
			{
				state->MeetsValueLimits = A.ValueAnalysis(false);
				
				if (state->MeetsValueLimits)
				{
					++miniGrid[0][i][j];
					++state->nSuccess;
					state->SuccessfulFracs.push_back(i);
					state->SuccessfulTaus.push_back(j);
					A.SaveDerivedParams();
				}
				
			}
			

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
	
	ActiveGalaxies.resize(pp.NThreads);
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
						if (ActiveGalaxies[j].nSuccess > 0)
						{
							SolvedGalaxies.push_back(ActiveGalaxies[j]);
							if (SolvedGalaxies.size() > copy.IterationSaveInterval)
							{
								SaveGrid(copy,false,k);
								SolvedGalaxies.resize(0);
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
		ActiveGalaxies[currentThread] = copy;
		persei[currentThread] = std::thread(LaunchProcess,&miniGrids[currentThread],k,currentThread);
		
		noThreadAssigned = true;	
		
		const int nPrints = pp.NThreads;
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
			if (ActiveGalaxies[j].nSuccess > 0)
			{
				SolvedGalaxies.push_back(ActiveGalaxies[j]);
			}
		}
	}
	std::cout << "All threads closed" << std::endl;
	
}

void VerificationMode(ParameterPack pp)
{
	int N = 100;
	srand (time(NULL));

	ParameterPack currentState;
	for (int i = 0; i < N; ++i)
	{
		currentState = RandomiseGalaxy(pp);
		int s1 = rand() % pp.NGrid;
		int s2 = rand() % pp.NGrid;
		
		currentState.tauColls.IterateValue(s1);
		currentState.collFrac.IterateValue(s2);
		
		std::string name = "verificationRun_" + std::to_string(i);
		Annulus A = Annulus(&currentState);
	
	
		A.SaveAnnulus(name);
		
		bool success = A.ValueAnalysis(true);
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
		
		Annulus A = Annulus(&pp);
		A.SaveAnnulus("SingleEvaluation");
		
	}
	
	if (pp.Mode == 1)
	{
		InitialiseSaveGrids(pp);
		IterationMode(pp);
		SaveGrid(pp,true,pp.NRandomGalaxies);
	}

	if (pp.Mode == 2)
	{
		VerificationMode(pp);
	}


}


