#include "IteratorClasses.h"
#include "MassReservoir.h"


inline double Cutoff(double t, double tau, double w)
{	
	if (t < tau - w/2)
	{
		return 1.0;
	}
	if (t > tau + 2/w)
	{
		return 0.0;
	}
	
	return (tau + w/2 - t)/w;
}

ISMIterator::ISMIterator(ParameterPack pp,std::vector<double> & TimeVector)
{
	PP = pp;
	
	
	Iterate(TimeVector);
}


void ISMIterator::Iterate(std::vector<double> & TimeVector)
{

	ISM = MassReservoir(TimeVector,PP,false);

	EFHIterator(PP,this);
}


EFHIterator::EFHIterator(ParameterPack pp, ISMIterator * parent)
{
	PP = pp;
	Parent = parent;
	
	LoadVectors();
	
	Iterate();
}

double EFHIterator::tauSN(int n)
{
	return PP.tauSNIa_Min + (PP.tauSNIa_Max - PP.tauSNIa_Min)*(float)n/(PP.tauSNIa_N - 1);
}

void EFHIterator::LoadVectors()
{
	//Load the basic vector, note that the E0 value changes as tauSNia changes the evaluation time
	EInf = 0;
	E0 = std::vector(PP.tauSNIa_N,0.0);

	
	std::cout << PP.tauSNIa_N << std::endl;
	std::cout << PP.tauColls_N << std::endl;
	
	F0_List = std::vector(PP.tauSNIa_N,std::vector(PP.tauColls_N,std::vector(PP.collWidth_N,0.0)));
	FInf_List= std::vector(PP.tauColls_N,std::vector(PP.collWidth_N,0.0));
	
	H0_NSM_List = std::vector(PP.tauSNIa_N, std::vector(PP.tauNSM_N, std::vector(PP.nuNSM_N,0.0)));
	HInf_NSM_List = std::vector(PP.tauNSM_N, std::vector(PP.nuNSM_N,0.0));
	
	HInf_SNIa_List = std::vector(PP.tauSNIa_N, std::vector(PP.nuSNIa_N,0.0));
	
	//loop over tauSNIa
	double prevT = 0;
	for (int i = 0; i < PP.tauSNIa_N; ++i)
	{
		
		double t = tauSN(i);
		double prev = 0;
		if (i > 0)
		{
			prev = E0[i - 1];
		}
		E0[i] = Evalues(t,prev,prevT);		
		
		
		//loop over tau and w
		for (int j = 0; j < PP.tauColls_N; ++j)
		{
			double tau = PP.tauColls_Min + (PP.tauColls_Max - PP.tauColls_Min)/(PP.tauColls_N - 1) * j;
			for (int k = 0; k < PP.collWidth_N; ++k)
			{
				double width = PP.collWidth_Min + (PP.collWidth_Max - PP.collWidth_Min)/(PP.collWidth_N - 1) * k;
				
				double val;
				if (t < tau - width/2)
				{
					val = E0[i];
				}
				else
				{
					prev = 0;
					if (i > 0)
					{
						prev = F0_List[i-1][j][k];
					}
					val =Fvalues(t,tau,width,prev,prevT);
				}
				F0_List[i][j][k] = val;
			}
		}
		
		
		//loop over 
		
		
		prevT = t;
	}
	
	if (PP.tMax >= PP.tauSNIa_Max)
	{
		EInf = Evalues(PP.tMax,E0[PP.tauSNIa_N-1],PP.tauSNIa_Max);
	}
	else
	{
		EInf = Evalues(PP.tMax,0,0);
	}

	for (int j = 0; j < PP.tauColls_N; ++j)
	{
		double tau = PP.tauColls_Min + (PP.tauColls_Max - PP.tauColls_Min)/(PP.tauColls_N - 1) * j;
		for (int k = 0; k < PP.collWidth_N; ++k)
		{
			double width = PP.collWidth_Min + (PP.collWidth_Max - PP.collWidth_Min)/(PP.collWidth_N - 1) * k;
			
			if (PP.tMax >= PP.tauSNIa_Max)
			{
				FInf_List[j][k] = Fvalues(PP.tMax,tau,width,F0_List[PP.tauSNIa_N -1][j][k],PP.tauSNIa_Max);   
			}
			else
			{
				FInf_List[j][k] = Fvalues(PP.tMax,tau,width,0,0);
			}
		}
	}

	
	//now for the difficult loops
	prevT = 0;
	for (int i = 0; i < PP.tauSNIa_N; ++i)
	{	
		double tSN = tauSN(i);
		for (int j = 0; j < PP.nuSNIa_N; ++j)
		{
			double nu = PP.nuSNIa_Min + (PP.nuSNIa_Max - PP.nuSNIa_Min)*j/(PP.nuSNIa_N - 1);
			
				
			HInf_SNIa_List[i][j] = Hvalues(PP.tMax,tSN,nu,0,tSN);
		}
		
		
		
		for (int k = 0; k < PP.tauNSM_N; ++k)
		{
			double tNSM = PP.tauNSM_Min + (PP.tauNSM_Max - PP.tauNSM_Min)* k/(PP.tauNSM_N - 1);
			for (int m = 0; m < PP.nuNSM_N; ++m)
			{
				double nu = PP.nuNSM_Min + (PP.nuNSM_Max - PP.nuNSM_Min)*m/(PP.nuNSM_N - 1);
				
				double prev = 0;
				if (i > 0)
				{
					prev = H0_NSM_List[i-1][k][m];
				}
				
				H0_NSM_List[i][k][m] = Hvalues(tSN,tNSM,nu,prev,prevT);
			}
		}
		
		prevT = tSN;
	}
	
	for (int k = 0; k < PP.tauNSM_N; ++k)
	{
		double tNSM = PP.tauNSM_Min + (PP.tauNSM_Max - PP.tauNSM_Min)* k/(PP.tauNSM_N - 1);
		for (int m = 0; m < PP.nuNSM_N; ++m)
		{
			double nu = PP.nuNSM_Min + (PP.nuNSM_Max - PP.nuNSM_Min)*m/(PP.nuNSM_N - 1);
			
			
			if (PP.tMax >= PP.tauSNIa_Max)
			{
				HInf_NSM_List[k][m] = Hvalues(PP.tMax,tNSM,nu,H0_NSM_List[PP.tauSNIa_N -1][k][m],PP.tauSNIa_Max);
			}
			else
			{
				HInf_NSM_List[k][m] = Hvalues(PP.tMax,tNSM,nu,0,0);
			}
			
		}
	}

	
	
	
}


double EFHIterator::Gvalues(double t, double tau, double nu)
{
	double sum = 0;
	
	if (t <= tau)
	{
		return 0.0;
	}
	
	int N = (t - tau)/(PP.timeStep*3)+1;
	double dt = (t - tau)/(N);
	double mult = PP.nuSFR * dt;
	sum = mult/2* ( Parent->ISM.ColdGas(t-tau) + Parent->ISM.ColdGas(0) * exp(-PP.nuSFR*(t - tau)));
	
	for (int i = 1; i < N; ++i)
	{
		double tpp = tau + i*dt;
		double sfr = mult * Parent->ISM.ColdGas(tpp);
		double ex = exp(-PP.nuSFR*(tpp - tau));
		
		sum += sfr * ex;
	}
	return sum;
}

double EFHIterator::Evalues(double t, double prevVal, double prevT)
{
	double sum = prevVal;
	double dt = PP.timeStep;
	
	int N = (t - prevT)/dt + 1;
	
	double dtPrime = (t-prevT)/N;
	double mult = PP.nuSFR * dtPrime;
	
	sum += (Parent->ISM.ColdGas(prevT) + Parent->ISM.ColdGas(t))*mult/2;
	for (int i = 1; i < N; ++i)
	{
		double tp = prevT + i*dtPrime;
		sum += Parent->ISM.ColdGas(tp) * mult;
	}
	
	return sum;
}

double EFHIterator::Fvalues(double t, double tau, double w, double prevVal, double prevT)
{
	double sum = prevVal;
	
	double dt = PP.timeStep;	
	int N = (t - prevT)/dt + 1;
	double dtPrime = (t-prevT)/N;
	double mult = PP.nuSFR * dtPrime;
	
	sum += (Parent->ISM.ColdGas(prevT)*Cutoff(prevT,tau,w) + Parent->ISM.ColdGas(t)*Cutoff(t,tau,w))*mult/2;
	for (int i = 1; i < N; ++i)
	{
		double tp = prevT + i*dtPrime;
		sum += Parent->ISM.ColdGas(tp) * mult * Cutoff(tp,tau,w);
	}
	
	return sum;
}

double EFHIterator::Hvalues(double t, double tau, double nu, double prevVal, double prevT)
{
	if (t <= tau)
	{
		return 0.0;
	}
	double sum = prevVal;
	
	
	double dt = PP.timeStep;	
	int N = (t - prevT)/dt + 1;
	double dtPrime = (t-prevT)/N;
	
	sum += dtPrime/2 * (Gvalues(prevT,tau,nu) + Gvalues(t,tau,nu));
	for (int i = 0; i < N; ++i)
	{
		double tp = prevT + i * dtPrime;
		sum += Gvalues(tp,tau,nu) * dtPrime;
	}
	
	return sum;
}

void EFHIterator::Iterate()
{
	
}
