#pragma once

#include <string>

// global variable store for command-line modification

extern std::string FILEROOT;
extern int IntegrationSteps;
extern int Mode;

//calibration data

extern double FeH_SN;
extern double MgFe_SN;
extern double MgFe_Sat;
extern double EuMg_SN;
extern double sProcFrac;
extern double collFrac;


//galaxy parameters
extern double galaxyM0;
extern double galaxyM1;
extern double galaxyM2;
extern double galaxyB1;
extern double galaxyB2;
extern double galaxyScaleLength;
extern double nuSFR;


//uncalibrated stuff
extern double nuColls;
extern double nuSNIa;
extern double tauColls;
extern double tauSNIa;
extern double tauNSM;
extern double nuNSM;
extern double tauInf;
extern double collWidth;
