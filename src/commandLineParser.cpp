#include "commandLineParser.h"
#include "ParameterPack.h"
#include <string.h>
#include <iostream>
#include <vector>
#include <sys/types.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h> 

ParameterPack pp = ParameterPack();

std::vector<string> integerGlobalTriggers = {"-steps","-mode","-threads","-random","-save","-grid"};
std::vector<int *> integerGlobalPointers = {&pp.IntegrationSteps,&pp.Mode,&pp.NThreads,&pp.NRandomGalaxies,&pp.SaveValue,&pp.NGrid};


std::vector<string> doubleGlobalTriggers =  {};
std::vector<double *> doubleGlobalPointers = { };


//fractions are doubles that are constrained to be between zero and 1
std::vector<string> fractionGlobalTriggers = {}; 
std::vector<double *> fractionGlobalPointers= {};


//fraction pairs are pairs of doubles between 0 and 1, which sum to 1. Therefore two variables need to be passed, formatted as a std::vector.
std::vector<string> fractionPairGlobalTriggers = {};
std::vector<std::vector<double *>> fractionPairGlobalPointers= {};


std::vector<string> booleanGlobalTriggers = {};
std::vector<bool *> booleanGlobalPointers = {};


//special triggers are those which need a dedicated function to do their job, such as -dir, which needs to create directories etc. 
//special functions are forward-declared in the .h file. 
std::vector<string> specialGlobalTriggers= {"-h", "-dir","grid"};
parseFunctions specialFuncs[] = {help, changeFileRoot};




bool changeIntegerParameter(char * newValue, int * host, char * callName)
{
	try
	{
		double convertedVal = stoi(newValue);
		host[0] = convertedVal;
		std::cout << "Changing the " << callName << " value. Value is now: " << convertedVal << std::endl;
		return true;
	}
	catch (const std::exception& e)
	{
		std::cout << "\n\n\nERROR: A problem was encountered trying to parse" << callName << ". Error message is as follows" << std::endl;
		std::cout << e.what() << std::endl;
		return false;
	}
}

bool changeDoubleParameter(char * newValue, double * host, char * callName)
{
	try
	{
		double convertedVal = stod(newValue);
		host[0] = convertedVal;
		std::cout << "Changing the " << callName << " value. Value is now: " << convertedVal << std::endl;
		return true;
	}
	catch (const std::exception& e)
	{
		std::cout << "\n\n\nERROR: A problem was encountered trying to parse" << callName << ". Error message is as follows" << std::endl;
		std::cout << e.what() << std::endl;
		return false;
	}
}

bool changeFractionParameter(char * newValue, double * host, char * callName)
{
	try
	{
		double convertedVal = stod(newValue);
		
		if (convertedVal < 0 || convertedVal > 1)
		{
			std::cout << "\n\nERROR: The value passed to " << callName << " was not a valid fraction." << std::endl;
			return false;
		}
		
		host[0] = convertedVal;
		std::cout << "Changing the " << callName << " value. Value is now: " << convertedVal << std::endl;
		return true;
	}
	catch (const std::exception& e)
	{
		std::cout << "\n\n\nERROR: A problem was encountered trying to parse" << callName << ". Error message is as follows" << std::endl;
		std::cout << e.what() << std::endl;
		return false;
	}
}

bool changeFractionPairParameter(char * newValue, std::vector<double *> host, char * callName)
{
	try
	{
		double convertedVal = stod(newValue);
		if (convertedVal < 0 || convertedVal > 1)
		{
			std::cout << "\n\nERROR: The value passed to " << callName << " was not a valid fraction." << std::endl;
			return false;
		}
		
		host[0][0] = convertedVal;
		host[1][0] = 1.0 - convertedVal;
		std::cout << "Changing the " << callName << " value. Value is now: " << convertedVal << std::endl;
		return true;
	}
	catch (const std::exception& e)
	{
		std::cout << "\n\n\nERROR: A problem was encountered trying to parse" << callName << ". Error message is as follows" << std::endl;
		std::cout << e.what() << std::endl;
		return false;
	}
}

bool changeBooleanParameter(char * newValue, bool * host, char * callName)
{
	try
	{
		bool convertedVal = (bool)stoi(newValue);
		host[0] = convertedVal;
		std::cout << "Changing the " << callName << " value. Value is now: " << convertedVal << std::endl;
		return true;
	}
	catch (const std::exception& e)
	{
		std::cout << "\n\n\nERROR: A problem was encountered trying to parse" << callName << ". Error message is as follows" << std::endl;
		std::cout << e.what() << std::endl;
		return false;
	}
}

bool help(char* arg)
{
	ifstream helpFile("DataFiles/helpList.dat");
	if (!helpFile.is_open())
	{
		std::cout << "\n\nERROR: Could not find the help files. Something terrible has occured.\n\n" << std::endl;
		return false;
	}
	
	string rawLine;
	while (getline(helpFile, rawLine) )
	{
		std::vector<string> line = split(rawLine,'\t');
		for (int i =0; i < line.size(); ++i)
		{
			std::cout << "\t" << setw(15) << left << line[i];
		}
		std::cout << "\n";
	}
	return true;
}

bool changeGridSize(char* arg)
{
	int NGrid = stoi(arg);
	
	std::vector<IterableParameter<double> *> iters = {&pp.tauColls, &pp.collFrac};
	
	for (int i = 0; i < iters.size(); ++i)
	{
		iters[i][0] = IterableParameter<double>(iters[i]->Value, iters[i]->MinValue, iters[i]->MaxValue, NGrid);
	}
	
}

bool changeFileRoot(char* arg)
{

	string root = (string)arg;
	pp.FILEROOT = root;
	int n = root.size();
	std::string lastCharacter = root.substr(n -1);
	bool needsSlash = (lastCharacter.compare("/"));
	std::cout << "Saving Files to directory '" << root << "'. ";
	const char *fileChar = root.c_str();
	DIR *dir = opendir(fileChar);
	if (dir)
	{
		std::cout << "Dir exists" << std::endl;
	}
	else
	{
		std::cout << "Dir does not exist. Making Dir" << std::endl;
		try
		{
			string command = "mkdir -p ";
			command.append(root);
			const char *commandChar = command.c_str(); 
			const int dir_err = system(commandChar);
			
			command = "mkdir -p " + root+ "/IterationChecker/";
			const char *commandChar2 = command.c_str(); 
			const int dir_err2 = system(commandChar2);
		}
		catch (const std::exception& e)
		{
			std::cout << "\n\n\nERROR: A problem was encountered trying to create directory. Error message is as follows" << std::endl;
			std::cout << e.what() << std::endl;
			return false;
		}
	}
	if (needsSlash)
	{
		pp.FILEROOT.append("/");
	}	
	return true;
}


bool helpNeeded(int argc, char** args)
{
	for (int i = 1; i < argc;++i)
	{
		string temp = args[i];
		if (temp.compare("-h")==0)
		{
			return true;
		}
	}
	return false;
}
ParameterPack parseCommandLine(int argc, char** argv)
{
	//parse command line arguments
	bool linesParsed = true;
	
	
	if (helpNeeded(argc,argv))
	{
		help(argv[0]);
		exit(0);
	}
	
	
	if (argc > 1)
	{
		std::cout << "Processing " << argc - 1 << " command line arguments" << std::endl;
		for (int i = 1; i < argc-1; i +=2)
		{
			string temp = argv[i];
			bool found = false;
			std::vector<std::vector<string>> triggers = {specialGlobalTriggers, integerGlobalTriggers, doubleGlobalTriggers, booleanGlobalTriggers, fractionGlobalTriggers, fractionPairGlobalTriggers};
			
			//loop through the standard assignment std::vectors
			for (int k = 0; k < triggers.size(); ++k)
			{
				int N = triggers[k].size();
				for (int p = 0; p < N; ++p)
				{
					if (temp.compare(triggers[k][p]) == 0)
					{
						//although the function calls here look identical (so could be compacted into a looped-array over a std::vector?), note that
						//the vector-of-pointer-arguments contain different data types, so you'd need to construct some kind of wrapper to manage it. Too much effort. Switch/case does the job.
						found = true;
						switch(k)
						{
							case 0:
								linesParsed = specialFuncs[p](argv[i+1]);	
								break;
							case 1:
								linesParsed = changeIntegerParameter(argv[i+1], integerGlobalPointers[p], argv[i]);
								break;
							case 2:
								linesParsed = changeDoubleParameter(argv[i+1], doubleGlobalPointers[p],argv[i]);
								break;
							case 3:
								linesParsed = changeBooleanParameter(argv[i+1], booleanGlobalPointers[p],argv[i]);
								break;
							case 4:
								linesParsed = changeFractionParameter(argv[i+1], fractionGlobalPointers[p], argv[i]);
								break;
							case 5:
								linesParsed = changeFractionPairParameter(argv[i+1], fractionPairGlobalPointers[p], argv[i]);
								break;
							default:
								linesParsed = false;
								found = false;
								std::cout << "Something went wrong in the default parameter assignment routine" << std::endl;
								i = argc;
								break;
						}							
						p = triggers[k].size(); // might as well jump to end of loop. Not performance critical, however. 
						k = triggers.size();
					}
				}				
			}
			
			if (!found)
			{
				std::cout << "\n\n ERROR: An unknown command (" << temp << ") was passed to the CLP. Quitting.";
				i = argc;
				linesParsed = false;
			}
			
		}
	}
	
	if (argc % 2 == 0 && linesParsed)
	{
		std::cout << "\nWARNING: An uneven number of parameters was passed, so the final command (" << argv[argc-1] << ") was not examined. \nThis is a non-critical error, so code will continue.\n\n" << std::endl; 
	}
	
	pp.InitialisedCorrectly = linesParsed;
	return pp;
}

