#pragma once
#include <vector>
#include <string>
#include <iomanip>
#include "ParameterPack.h"
#include "stringManipulator.h"
//See .cpp for initialisation values
typedef bool (*parseFunctions) (char* arg);


//designated "special functions" (see .cpp for more detail)
bool help(char* arg);
bool changeFileRoot(char* arg);
bool changeGridSize(char* arg);
bool changeConstraints(char * arg);
bool changeGradientBounds(char * arg);
//Caller
ParameterPack parseCommandLine(int argc, char** argv);



