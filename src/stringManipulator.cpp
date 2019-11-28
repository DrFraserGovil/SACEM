#include "stringManipulator.h"
std::vector<string> split(const string& s, char delimiter)
{
	//dumb brute-force string splitter based on a delimiter
	std::vector<std::string> tokens;
	string token;
	istringstream tokenStream(s);
	while (getline(tokenStream, token, delimiter))
	{
		if( token.length() > 0) // empty rows not allowed
		{
			tokens.push_back(token);
		}
	}
	return tokens;
}
