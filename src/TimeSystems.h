#pragma once
#include <chrono>
#include <string.h>
#include <stdio.h>

std::string convertTime(double secs)
{
	std::ostringstream time;
	int hours = floor(secs/3600.0);
	int mins = floor((secs - 3600.0*hours)/60.0);
	int cumulative = hours*3600 + mins*60;
	double roundingAccuracy = 4;
	secs = round((secs - cumulative)*(pow(10,roundingAccuracy)))/pow(10,roundingAccuracy);
	if (hours > 0)
	{
		time << hours << " hour";
		if (hours > 1)
		{
			time  << "s";
		}
		time << " ";
	}
	if (mins > 0)
	{
		time << mins << " minute";
		if (mins > 1)
		{
			time << "s";
		}
		time << " ";
	}
	if ((int)secs > 0)
	{
		time << (int)secs << " second";
		if ( (int)secs > 1)
		{
			time << "s";
		}
	}
	else
	{
		time << round(secs*100)/100 << " seconds ";
	}
	return time.str();
}

void printTimeSince(std::chrono::time_point<std::chrono::high_resolution_clock> start,int len,int total)
{
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	double elapsedSecs = elapsed.count();
	
	double perc = round( 1000 * (float)len/total)/10;
	std::cout << "\033[A\33[2K\r[";
	std::cout << std::setw(1);
	for (int j = 0 ; j < floor(perc/5); ++j)
	{
		std::cout << "|";
	}
	for (int j = floor(perc/5); j < 20; ++j)
	{
		std::cout << " ";
	}
	std::cout << "] ";
	std::cout << "Elapsed: " << convertTime(elapsedSecs) << ", remaining: ";
	double secs = elapsedSecs*((float)total/len - 1);
	
	
	std::cout << convertTime(secs) << std::endl;
}
