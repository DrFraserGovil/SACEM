#include "commandLineParser.h"
#include "ParameterPack.h"
#include <string.h>
#include <iostream>
#include <vector>
#include <sys/types.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h> 
using namespace std;


ParameterPack pp = ParameterPack();

vector<string> integerGlobalTriggers = {"-steps","-mode"};
vector<int *> integerGlobalPointers = {&pp.IntegrationSteps,&pp.Mode};


vector<string> doubleGlobalTriggers =  {"-FeH","-MgFeSat","-MgFePlat","-EuMg","-M0","-M1", "-M2","-b1","-b2","-Rd","-nuSFR","-nuColl","-nuSN","-tauColl","-tauSN","-tauInf","-width"};
vector<double *> doubleGlobalPointers = { &pp.FeH_SN,&pp.MgFe_Sat,&pp.MgFe_SN,&pp.EuMg_SN,&pp.galaxyM0,&pp.galaxyM1,&pp.galaxyM2,&pp.galaxyB1,&pp.galaxyB2,&pp.galaxyScaleLength,&pp.nuSFR,&pp.nuColls,&pp.nuSNIa,&pp.tauColls,&pp.tauSNIa,&pp.tauInf,&pp.collWidth};


//fractions are doubles that are constrained to be between zero and 1
vector<string> fractionGlobalTriggers = {"-sFrac","-collFrac","-hotFrac"}; 
vector<double *> fractionGlobalPointers= {&pp.sProcFrac,&pp.collFrac, &pp.hotFrac};


//fraction pairs are pairs of doubles between 0 and 1, which sum to 1. Therefore two variables need to be passed, formatted as a vector.
vector<string> fractionPairGlobalTriggers = {};
vector<vector<double *>> fractionPairGlobalPointers= {};


vector<string> booleanGlobalTriggers = {};
vector<bool *> booleanGlobalPointers = {};


//special triggers are those which need a dedicated function to do their job, such as -dir, which needs to create directories etc. 
//special functions are forward-declared in the .h file. 
vector<string> specialGlobalTriggers= {"-h", "-dir"};
parseFunctions specialFuncs[] = {help, changeFileRoot};


bool changeIntegerParameter(char * newValue, int * host, char * callName)
{
	try
	{
		double convertedVal = stoi(newValue);
		host[0] = convertedVal;
		cout << "Changing the " << callName << " value. Value is now: " << convertedVal << endl;
		return true;
	}
	catch (const std::exception& e)
	{
		cout << "\n\n\nERROR: A problem was encountered trying to parse" << callName << ". Error message is as follows" << endl;
		cout << e.what() << endl;
		return false;
	}
}

bool changeDoubleParameter(char * newValue, double * host, char * callName)
{
	try
	{
		double convertedVal = stod(newValue);
		host[0] = convertedVal;
		cout << "Changing the " << callName << " value. Value is now: " << convertedVal << endl;
		return true;
	}
	catch (const std::exception& e)
	{
		cout << "\n\n\nERROR: A problem was encountered trying to parse" << callName << ". Error message is as follows" << endl;
		cout << e.what() << endl;
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
			cout << "\n\nERROR: The value passed to " << callName << " was not a valid fraction." << endl;
			return false;
		}
		
		host[0] = convertedVal;
		cout << "Changing the " << callName << " value. Value is now: " << convertedVal << endl;
		return true;
	}
	catch (const std::exception& e)
	{
		cout << "\n\n\nERROR: A problem was encountered trying to parse" << callName << ". Error message is as follows" << endl;
		cout << e.what() << endl;
		return false;
	}
}

bool changeFractionPairParameter(char * newValue, vector<double *> host, char * callName)
{
	try
	{
		double convertedVal = stod(newValue);
		if (convertedVal < 0 || convertedVal > 1)
		{
			cout << "\n\nERROR: The value passed to " << callName << " was not a valid fraction." << endl;
			return false;
		}
		
		host[0][0] = convertedVal;
		host[1][0] = 1.0 - convertedVal;
		cout << "Changing the " << callName << " value. Value is now: " << convertedVal << endl;
		return true;
	}
	catch (const std::exception& e)
	{
		cout << "\n\n\nERROR: A problem was encountered trying to parse" << callName << ". Error message is as follows" << endl;
		cout << e.what() << endl;
		return false;
	}
}

bool changeBooleanParameter(char * newValue, bool * host, char * callName)
{
	try
	{
		bool convertedVal = (bool)stoi(newValue);
		host[0] = convertedVal;
		cout << "Changing the " << callName << " value. Value is now: " << convertedVal << endl;
		return true;
	}
	catch (const std::exception& e)
	{
		cout << "\n\n\nERROR: A problem was encountered trying to parse" << callName << ". Error message is as follows" << endl;
		cout << e.what() << endl;
		return false;
	}
}

bool help(char* arg)
{
	ifstream helpFile("DataFiles/helpList.dat");
	if (!helpFile.is_open())
	{
		cout << "\n\nERROR: Could not find the help files. Something terrible has occured.\n\n" << endl;
		return false;
	}
	
	string rawLine;
	while (getline(helpFile, rawLine) )
	{
		std::vector<string> line = split(rawLine,'\t');
		for (int i =0; i < line.size(); ++i)
		{
			cout << "\t" << setw(15) << left << line[i];
		}
		cout << "\n";
	}
	return true;
}

bool changeFileRoot(char* arg)
{

	string FILEROOT = (string)arg;
	pp.FILEROOT = FILEROOT;
	int n = FILEROOT.size();
	std::string lastCharacter =FILEROOT.substr(n -1);
	bool needsSlash = (lastCharacter.compare("/"));
	cout << "Saving Files to directory '" << FILEROOT << "'. ";
	const char *fileChar = FILEROOT.c_str();
	DIR *dir = opendir(fileChar);
	if (dir)
	{
		cout << "Dir exists" << endl;
		if (needsSlash)
		{
			pp.FILEROOT.append("/");
		}	
	}
	else
	{
		cout << "Dir does not exist. Making Dir" << endl;
		try
		{
			string command = "mkdir -p ";
			command.append(FILEROOT);
			const char *commandChar = command.c_str(); 
			//~const int dir_err = system(commandChar);
			//~if (needsSlash)
			//~{
				//~FILEROOT.append("/");
			//~}	
		}
		catch (const std::exception& e)
		{
			cout << "\n\n\nERROR: A problem was encountered trying to create directory. Error message is as follows" << endl;
			cout << e.what() << endl;
			return false;
		}
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
		cout << "Processing " << argc - 1 << " command line arguments" << endl;
		for (int i = 1; i < argc-1; i +=2)
		{
			string temp = argv[i];
			bool found = false;
			vector<vector<string>> triggers = {specialGlobalTriggers, integerGlobalTriggers, doubleGlobalTriggers, booleanGlobalTriggers, fractionGlobalTriggers, fractionPairGlobalTriggers};
			
			//loop through the standard assignment vectors
			for (int k = 0; k < triggers.size(); ++k)
			{
				int N = triggers[k].size();
				for (int p = 0; p < N; ++p)
				{
					if (temp.compare(triggers[k][p]) == 0)
					{
						//although the function calls here look identical (so could be compacted into a looped-array over a vector?), note that
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
								cout << "Something went wrong in the default parameter assignment routine" << endl;
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
				cout << "\n\n ERROR: An unknown command (" << temp << ") was passed to the CLP. Quitting.";
				i = argc;
				linesParsed = false;
			}
			
		}
	}
	
	if (argc % 2 == 0 && linesParsed)
	{
		cout << "\nWARNING: An uneven number of parameters was passed, so the final command (" << argv[argc-1] << ") was not examined. \nThis is a non-critical error, so code will continue.\n\n" << endl; 
	}
	
	pp.InitialisedCorrectly = linesParsed;
	return pp;
}

