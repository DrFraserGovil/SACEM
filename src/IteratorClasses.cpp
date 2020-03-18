#include "IteratorClasses.h"
#include "MassReservoir.h"


void print(std::vector<double> &v)
{
	int width = 10;
	for (int i = 0; i < v.size(); ++i)
	{
		std::cout << std::setw(width) << std::left << v[i] << "\t";
	}
	std::cout << std::endl;
}


void printi(std::vector<int> &v)
{
	for (int i = 0; i < v.size(); ++i)
	{
		std::cout << v[i] << "\t";
	}
	std::cout << std::endl;
}
inline double Cutoff(double t, double tau, double w)
{	
	double v = 0;
	if (t < tau - w/2)
	{
		v = 1.0;
	}
	else
	{
		if (t > tau + w/2)
		{
			v =  0.0;
		}
		else
		{
			v = (tau + w/2 - t)/w;
		}
	}
		
	return v;
}

ISMIterator::ISMIterator(ParameterPack pp,std::vector<double> & TimeVector, std::vector<std::vector<int>> * grid)
{
	PP = pp;
	
	SuccessCount = grid;
	
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
	
	//LoadVectors();
	
	Iterate();
}


void EFHIterator::LoadVectors()
{
	//Load the basic vector, note that the E0 value changes as tauSNia changes the evaluation time
	int nTauSN = PP.tauSNIa.NSteps;
	int nNuSN = PP.nuSNIa.NSteps;
	int nTauColls = PP.tauColls.NSteps;
	int nCollWidths = PP.collWidth.NSteps;
	int nTauNSM = PP.tauNSM.NSteps;
	int nnuNSM = PP.nuNSM.NSteps;
	
	EInf = 0;
	E0_List= std::vector(nTauSN,0.0);

	
	F0_List = std::vector(nTauSN,std::vector(nTauColls,std::vector(nCollWidths,0.0)));
	FInf_List= std::vector(nTauColls,std::vector(nCollWidths,0.0));
	
	H0_NSM_List = std::vector(nTauSN, std::vector(nTauNSM, std::vector(nnuNSM,0.0)));
	HInf_NSM_List = std::vector(nTauNSM, std::vector(nnuNSM,0.0));
	
	HInf_SNIa_List = std::vector(nTauSN, std::vector(nNuSN,0.0));
	
	//loop over tauSNIa
	double prevT = 0;
	for (int i = 0; i < nTauSN; ++i)
	{
		double t = PP.tauSNIa.IntermediateValue(i);
		double prev = 0;
		if (i > 0)
		{
			prev = E0_List[i - 1];
		}
		E0_List[i] = Evalues(t,prev,prevT);		
				
		//loop over tau and w
		for (int j = 0; j < nTauColls; ++j)
		{
			double tau = PP.tauColls.IntermediateValue(j);
			std::cout << tau << std::endl;
			for (int k = 0; k < nCollWidths; ++k)
			{
				
				double width = PP.collWidth.IntermediateValue(k);
				
				double val;
				
				prev = 0;
				if (i > 0)
				{
					prev = F0_List[i-1][j][k];
				}
				val =Fvalues(t,tau,width,prev,prevT);
			
				std::cout << "\t" << val <<std::endl;
				F0_List[i][j][k] = val;
			}
		}
		
		
		//loop over 
		
		
		prevT = t;
	}
	

	
	if (PP.tMax >= PP.tauSNIa.MaxValue)
	{
		EInf = Evalues(PP.tMax,E0_List[nTauSN-1],PP.tauSNIa.MaxValue);
	}
	else
	{
		EInf = Evalues(PP.tMax,0,0);
	}


	std::cout << "FInfs" <<std::endl;
	for (int j = 0; j < nTauColls; ++j)
	{
		double tau = PP.tauColls.IntermediateValue(j);
		std::cout << tau <<std::endl;
		for (int k = 0; k < nCollWidths; ++k)
		{
			double width = PP.collWidth.IntermediateValue(k);
			
			if (PP.tMax >= PP.tauSNIa.MaxValue)
			{
				FInf_List[j][k] = Fvalues(PP.tMax,tau,width,F0_List[nTauSN -1][j][k],PP.tauSNIa.MaxValue);   
			}
			else
			{
				FInf_List[j][k] = Fvalues(PP.tMax,tau,width,0,0);
			}
			
			std::cout << "\t" << FInf_List[j][k] << std::endl;
		}
	}

	
	//now for the difficult loops
	prevT = 0;
	for (int i = 0; i < nTauSN; ++i)
	{	
		double tSN = PP.tauSNIa.IntermediateValue(i);
		for (int j = 0; j < nNuSN; ++j)
		{
			double nu = PP.nuSNIa.IntermediateValue(j);
			
			HInf_SNIa_List[i][j] = Hvalues(PP.tMax,tSN,nu,0,tSN);
		}
		
		
		
		for (int k = 0; k < nTauNSM; ++k)
		{
			double tNSM = PP.tauNSM.IntermediateValue(k);
			for (int m = 0; m < nnuNSM; ++m)
			{
				double nu = PP.nuNSM.IntermediateValue(m);
				
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
	
	for (int k = 0; k < nTauNSM; ++k)
	{
		double tNSM = PP.tauNSM.IntermediateValue(k);
		for (int m = 0; m < nnuNSM; ++m)
		{
			double nu = PP.nuNSM.IntermediateValue(m);
			
			
			if (PP.tMax >= PP.tauSNIa.MaxValue)
			{
				HInf_NSM_List[k][m] = Hvalues(PP.tMax,tNSM,nu,H0_NSM_List[nTauSN -1][k][m],PP.tauSNIa.MaxValue);
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
	double mult = PP.nuSFR.Value * dt;
	sum = mult/2* ( Parent->ISM.ColdGas(t-tau) + Parent->ISM.ColdGas(0) * exp(-PP.nuSFR.Value*(t - tau)));
	
	for (int i = 1; i < N; ++i)
	{
		double tpp = tau + i*dt;
		double sfr = mult * Parent->ISM.ColdGas(tpp);
		double ex = exp(-PP.nuSFR.Value*(tpp - tau));
		
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
	double mult = PP.nuSFR.Value * dtPrime;
	
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
	double mult = PP.nuSFR.Value * dtPrime;
	
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
	//iteration order:
	//1: tauSN, 2: nuSN, 3: tauNSM, 4: nuNSM, 5: tauColl, 6: collWidth
	
	ParameterPack copy = PP;
	
	double E0;
	double F0;
	double FInf;
	double H0_NSM;
	double HInf_NSM;
	double HInf_SN;
	
	//~ std::cout << "E List: " << std::endl;
	//~ print(E0_List);
	//~ std::cout << "\n\n";
	//~ std::cout << "F List:" << std::endl;
	//~ for (int m = 0; m < PP.tauColls.NSteps; ++m)
	//~ {
		//~ std::cout << std::setw(10) << std::left << PP.tauColls.IntermediateValue(m) << ":\t"; 
		//~ print(FInf_List[m]);
	//~ }
	
	EInf = Evalues(PP.tauInf,0,0);
	for (int i = 0; i < PP.tauSNIa.NSteps; ++i)
	{
		copy.tauSNIa.UpdateValue(i);
		E0 = Evalues(copy.tauSNIa.Value,0,0);
		for (int j = 0; j < PP.nuSNIa.NSteps; ++j)
		{
			copy.nuSNIa.UpdateValue(j);
			HInf_SN = Hvalues(copy.tauInf,copy.tauSNIa.Value, copy.nuSNIa.Value,0,0);
			for (int k = 0; k < PP.tauNSM.NSteps; ++k)
			{
				copy.tauNSM.UpdateValue(k);
				for (int l = 0; l < PP.nuNSM.NSteps; ++l)
				{
					copy.nuNSM.UpdateValue(l);
					H0_NSM = Hvalues(copy.tauSNIa.Value, copy.tauNSM.Value,copy.nuNSM.Value,0,0);
					HInf_NSM = Hvalues(copy.tauInf, copy.tauNSM.Value,copy.nuNSM.Value,0,0);
					for (int m = 0; m < PP.tauColls.NSteps; ++m)
					{
						copy.tauColls.UpdateValue(m);
						for (int n = 0; n < PP.collWidth.NSteps; ++n)
						{
							copy.collWidth.UpdateValue(n);
							F0 = Fvalues(copy.tauSNIa.Value, copy.tauColls.Value, copy.collWidth.Value,0,0);
							FInf =  Fvalues(copy.tauInf, copy.tauColls.Value, copy.collWidth.Value,0,0);
							
							//~ std::vector<double> v = {copy.tauColls.Value, copy.collWidth.Value, F0, E0, EInf, FInf};
							//~ print(v);
							
							CalibrationIterator ci = CalibrationIterator(copy,E0,EInf,F0,FInf,H0_NSM, HInf_NSM, HInf_SN, this);
							ci.Iterate();
						}
					}
				}
			}
			
		}
	}
}

CalibrationIterator::CalibrationIterator(ParameterPack pp, double e0, double eInf, double f0, double fInf, double h0_nsm, double hInf_nsm, double hInf_SN, EFHIterator* parent)
{
	PP = pp;
	
	E0 = e0;
	EInf = eInf;
	F0 = fInf;
	FInf = fInf;
	H0_NSM = h0_nsm;
	HInf_NSM = hInf_nsm;
	HInf_SN = hInf_SN;
	
	Parent = parent;
	
	alpha = 0;
	beta = 0;
	gamma = 0;
	delta = 0;
	epsilon = 0;
	eta = 0;
}

void CalibrationIterator::Calibrate(double mass)
{
	double Fcal = PP.FeH_SN.Value;
	double M0 = PP.MgFe_SN.Value;
	double MInf = PP.MgFe_Sat.Value;
	double Eps = PP.EuMg_SN.Value;
	double zeta = PP.sProcFrac.Value;
	double omega = PP.collFrac.Value;
	double nsmFrac = 1.0 - zeta - omega;
	
	
	alpha = mass/E0 * pow(10.0,Fcal);
	beta = alpha * EInf/HInf_SN*(pow(10.0,M0 - MInf) - 1.0);
	eta = alpha * pow(10.0,M0);
	
	
	double gammaFrac = zeta*E0/EInf + omega * F0/FInf + nsmFrac*H0_NSM/HInf_NSM;
	gamma = eta * pow(10.0,Eps) * zeta * E0/EInf / gammaFrac;
	delta = omega/zeta * EInf/FInf * gamma;
	epsilon = nsmFrac/zeta * gamma * EInf;
	
	//std::cout << alpha << " " << beta << " " << gamma << " " << delta << " " << epsilon << " " << eta << std::endl;
}

void CalibrationIterator::Iterate()
{
	double mass = Parent->Parent->ISM.ColdGas(PP.tauSNIa.Value) + Parent->Parent->ISM.HotGas(PP.tauSNIa.Value) + Parent->Parent->ISM.Stars(PP.tauSNIa.Value);
	for (int i = 0; i < PP.FeH_SN.NSteps; ++i)
	{
		
		PP.FeH_SN.UpdateValue(i);
		
		for(int j = 0; j < PP.MgFe_SN.NSteps; ++j)
		{
			PP.MgFe_SN.UpdateValue(j);
			
			for (int k = 0; k < PP.MgFe_Sat.NSteps; ++k)
			{
				PP.MgFe_Sat.UpdateValue(k);
				for (int l = 0; l < PP.EuMg_SN.NSteps; ++l)
				{
					PP.MgFe_SN.UpdateValue(l);
					for (int m = 0; m < PP.sProcFrac.NSteps; ++m)
					{
						
						PP.sProcFrac.UpdateValue(m);
						for (int n = 0; n < PP.collFrac.NSteps; ++n)
						{
							PP.collFrac.UpdateValue(n);
							Calibrate(mass);
							Evaluate();
						}
					}
				}
			}
			
		}
		
	}
}


void CalibrationIterator::Evaluate()
{
	double finalFe = alpha * EInf + beta * HInf_SN;
	double finalEu = gamma * EInf + delta * FInf + epsilon * HInf_NSM;
	
	double finalVal = log10(finalEu/finalFe);
	
	//~ std::vector<double> v1 = {PP.FeH_SN.Value, PP.MgFe_SN.Value, PP.MgFe_Sat.Value, PP.EuMg_SN.Value, PP.sProcFrac.Value, PP.collFrac.Value};
	//~ std::vector<double> v15 = {E0,EInf,F0,FInf, H0_NSM, HInf_NSM, HInf_SN};
	//~ std::vector<double> v2 = {alpha, beta, gamma, delta, epsilon, eta};

	
	//~ std::vector<double> vals = {PP.tauSNIa.Value, PP.tauColls.Value, PP.collWidth.Value, PP.collFrac.Value, finalFe, finalEu, finalVal};
	//~ print(v1);
	//~ print(v15);
	//~ print(v2);
	//~ print(vals);
	//~ std::cout << "\n\n";
	
	if (finalVal >= PP.finalEuFe_Min && finalVal <= PP.finalEuFe_Max)
	{
		int index1 = PP.tauColls.Index();
		int index2 = PP.collFrac.Index();
		//std::cout << index1 << " " << index2 << std::endl;
		
		
		
		
		++Parent->Parent->SuccessCount[0][index1][index2];
	}
	
}
