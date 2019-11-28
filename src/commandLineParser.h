#pragma once
#include <vector>
#include <string>
#include <iomanip>
#include "defaultValues.h"
#include "stringManipulator.h"
//See .cpp for initialisation values
typedef bool (*parseFunctions) (char* arg);


//designated "special functions" (see .cpp for more detail)
bool help(char* arg);
bool changeFileRoot(char* arg);


//Caller
bool parseCommandLine(int argc, char** argv);



