#include "ProcessClasses.h"


double Accretion::Mass(double t)
{
	double M = MBar;
	for (int i = 0; i < PP.galaxyMs.size(); ++i)
	{
		M-=PP.galaxyMs[i] * exp(-PP.Betas[i] * t);
	}
	return M;
}

double StarFormation::Rho(double t)
{
	return RealDot(SFRVector, E(t));
}


void CCSN::GeneratePrefactors()
{
	CCSNVector = CCSNStyleVector();

	TotalVector = Hadamard(TotalIntegrator(),SFRVector);
	
}

double CCSN::Total(double t)
{
	return RealDot(TotalVector, F(t) - F(0));
}

double CCSN::Count(double t)
{
	return RealDot(CCSNVector, E(t) );
}

void Collapsar::GeneratePrefactors()
{
	UntouchedVector = CCSNStyleVector();
	
	ComplexVector Hlambda = Integrator(CoolingID);
	ComplexVector Hnu = Integrator(FormationID);
	//ccsn hot vector
	ComplexVector Q = Hadamard(SFRVector,Hlambda);
	Q.Real[CoolingID] -= RealDot(SFRVector, Hlambda);
	Q = hotFrac * Q;
	

	//Collapsar hot
	
	ComplexVector constTerm = SFRVector * T;
	ComplexVector linearTerm = SFRVector;
	ComplexVector Hsq = Hadamard(Hlambda,Hlambda);
	
	ComplexVector partialS = hotFrac * (Hadamard(constTerm,Hlambda) + Hadamard(linearTerm, Hsq));
	ComplexVector TVec = hotFrac * ( Hadamard(linearTerm, Hlambda));
	
	double C;
	if (T - Delta > 0)
	{
		ComplexVector vBit = Delta * Q + (T - Delta) * TVec - partialS;
		
		C = RealDot(vBit, E(T - Delta)) * exp(Powers.Real[CoolingID]*(T - Delta) );
	}
	else
	{
		C = -1 * RealDot(partialS, E(0) );
	}
	
	ComplexVector S = partialS;
	S.Real[CoolingID] = C;


	//collapsar mid cold
	
	ComplexVector U = (1.0 - hotFrac) * T * SFRVector + Powers.Real[CoolingID] * S;
	ComplexVector V = (1.0 - hotFrac) * SFRVector + Powers.Real[CoolingID] * TVec;
	
	
	ComplexVector Z = 1.0/Delta * Hadamard(V,Hnu);
	ComplexVector  a = 1.0/Delta * Hadamard(Hnu, U + Hadamard(V,Hnu) );

	double coldConst;
	if (T - Delta > 0)
	{
		ComplexVector vBit = UntouchedVector - a + (T-Delta)*Z;
		coldConst = exp(Powers.Real[FormationID]*(T-Delta) ) * RealDot(vBit,E(T - Delta) );
	}
	else
	{
		coldConst = - RealDot(a,E(0) );
	}
	
	IntermediaryVector_Linear = Z;
	IntermediaryVector_Constant = a;
	IntermediaryVector_Constant.Real[FormationID] = coldConst;
	
	
	
	/// final domain constraints
	
	double HotLimit = 1.0/Delta * RealDot( S - T*TVec, E(T) );
	double ColdLimit = Middle(T);
	
	 FinalVector = ComplexVector(Powers.size());
	 double lambdaFac = Powers.Real[CoolingID]/(Powers.Real[FormationID] - Powers.Real[CoolingID] ) * HotLimit ;
	
	 FinalVector.Real[CoolingID] = lambdaFac * exp(Powers.Real[CoolingID] * T);
	
	FinalVector.Real[FormationID] = ( ColdLimit - lambdaFac) *exp(Powers.Real[FormationID] * T);
	
	
	InitialTotal = Hadamard(TotalIntegrator(),SFRVector);
	
	ComplexVector s1 = Hadamard(TotalIntegrator(),TotalIntegrator() );
	s1.Real[0] = 0;
	MidConstantTotal = Hadamard(s1 + T * TotalIntegrator(), SFRVector);
	
	ComplexVector s2 = TotalIntegrator();
	s2.Real[0] = 0.5;
	MidLinearTotal = Hadamard(s2, SFRVector);
	
	
	TotalConst = 0; 
	
	if (T - Delta > 0)
	{
		TotalConst = RealDot(InitialTotal, F(T-Delta) - F(0) );
		TotalConst -= TotalMid(T- Delta);
	}
	else
	{
		TotalConst = -TotalMid(0);
	}
}



double Collapsar::Total(double t)
{
	if (T<=0)
	{
		return 0;
	}
	if (t < T - Delta)
	{
		return RealDot(InitialTotal, F(t) - F(0) );
	}
	
	if (t > T)
	{
		return TotalMid(T);
	}
	return TotalMid(t);
	
}

double Collapsar::TotalMid(double t)
{
	return TotalConst + 1.0/Delta * RealDot( MidConstantTotal - t* MidLinearTotal, F(t) );
}

double Collapsar::Count(double t)
{
	if (T <= 0)
	{
		return 0;
	}
	if (t < T - Delta)
	{
		return Initial(t);
	}
	if (t > T)
	{
		return End(t);
	}
	
	
	return Middle(t);
}

double Collapsar::Initial(double t)
{
	return RealDot(UntouchedVector, E(t) );
}

double Collapsar::End(double t)
{
	return RealDot(FinalVector,E(t) );
}

double Collapsar::Middle(double t)
{
	ComplexVector vectorBit = IntermediaryVector_Constant - t * IntermediaryVector_Linear;

	return RealDot( vectorBit, E(t) );	
}


void Decayer::GeneratePrefactors()
{
	ComplexVector W = nu * GenericIntegrator(SFRVector,DelayID);
	
	TotalVector = Hadamard(TotalIntegrator(),W);
	
	
	ComplexVector beta = hotFrac * GenericIntegrator(W,CoolingID);
	
	ComplexVector gamma = (1.0 - hotFrac) * W + Powers.Real[CoolingID] * beta;
	
	DecayVector = GenericIntegrator(gamma, FormationID);
	
	
};

double Decayer::Count(double t)
{
	if (t < Tau)
	{
		return 0;
	}
	
	return RealDot(DecayVector, E(t - Tau));
}

double Decayer::Total(double t)
{
	if (t < Tau)
	{
		return 0;
	}
	
	return RealDot(TotalVector, F(t - Tau) - F(0) );
}

double Decayer::Yield(double t)
{
	if (t < Tau)
	{
			return 0;
	}
	
	ComplexVector H = Integrator(DelayID);
	ComplexVector prefac = Hadamard(H,SFRVector);
	prefac.Real[DelayID] = - RealDot(H,SFRVector);
	
	prefac = nu*prefac;
	
	return RealDot(prefac,E(t-Tau));
	
	
}
