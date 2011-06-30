/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "physconst.h"
#include "thirdparty.h"
#include "continuum.h"

STATIC double RealF2_1( double alpha, double beta, double gamma, double chi );
STATIC complex<double> Hypergeometric2F1( complex<double> a, complex<double> b, complex<double> c,
		double chi, long *NumRenorms, long *NumTerms );
STATIC complex<double> F2_1( complex<double> alpha, complex<double> beta, complex<double> gamma,
		double chi, long *NumRenormalizations, long *NumTerms );
STATIC complex<double> HyperGeoInt( double v );
STATIC complex<double> qg32complex(	double xl, double xu, complex<double> (*fct)(double) );
STATIC double GauntIntegrand( double y );
STATIC double FreeFreeGaunt( double x );
STATIC double DoBeckert_etal( double etai, double etaf, double chi );
STATIC double DoSutherland( double etai, double etaf, double chi );

/* used to keep intermediate results from over- or underflowing.	*/
static const complex<double> Normalization(1e100, 1e100);
static complex<double> CMinusBMinus1, BMinus1, MinusA;
static double GlobalCHI;
static double Zglobal, HNUglobal, TEglobal;

double cont_gaunt_calc( double temp, double z, double photon )
{
	double gaunt, u, gamma2;

	Zglobal = z;
	HNUglobal = photon;
	TEglobal = temp;

	u = TE1RYD*photon/temp;
	gamma2 = TE1RYD*z*z/temp;

	if( log10(u)<-5. )
	{
		/* this cutoff is where these two approximations are equal. */
		if( log10( gamma2 ) < -0.75187 )
		{
			/* given as eqn 3.2 in Hummer 88, original reference is
			 * >>refer gaunt	ff	Elwert, G. 1954, Zs. Naturforshung, 9A, 637 */
			gaunt = 0.551329 * ( 0.80888 - log(u) );
		}
		else
		{
			/* given as eqn 3.1 in Hummer 88, original reference is
			 * >>refer gaunt	ff	Scheuer, P.A.G. 1960, MNRAS, 120, 231 */
			gaunt = -0.551329 * (0.5*log(gamma2) + log(u) + 0.056745);
		}
	}
	else
	{
		/* Perform integration.	*/
		gaunt =  qg32( 0.01, 1., GauntIntegrand );
		gaunt += qg32( 1., 5., GauntIntegrand );
	}
	ASSERT( gaunt>0. && gaunt<100. );

	return gaunt;
}

STATIC double GauntIntegrand( double y )
{
	double value;
	value =	FreeFreeGaunt( y ) * exp(-y);
	return value;
}

STATIC double FreeFreeGaunt( double x )
{
	double Csum, zeta, etaf, etai, chi, gaunt, z, InitialElectronEnergy, FinalElectronEnergy, photon;
	bool lgSutherlandOn = false;
	long i;

	z = Zglobal;
	photon = HNUglobal;
	ASSERT( z > 0. );
	ASSERT( photon > 0. );

	/* The final electron energy should be xkT + hv and the initial, just xkT	*/
	InitialElectronEnergy = sqrt(x) * TEglobal/TE1RYD;
	FinalElectronEnergy = photon + InitialElectronEnergy;
	ASSERT( InitialElectronEnergy > 0. );

	/* These are the free electron analog to a bound state principal quantum number.*/
	etai = z/sqrt(InitialElectronEnergy);	
	etaf = z/sqrt(FinalElectronEnergy);
	ASSERT( etai > 0. );
	ASSERT( etaf > 0. );
	chi = -4. * etai * etaf / POW2( etai - etaf );
	zeta = etai-etaf;

	if( etai>=130.) 
	{
		if( etaf < 1.7 )
		{
			/* Brussard and van de Hulst (1962), as given in Hummer 88, eqn 2.23b	*/
			gaunt = 1.1027 * (1.-exp(-2.*PI*etaf));
		}
		else if( etaf < 0.1*etai )
		{
			/* Hummer 88, eqn 2.23a	*/
			gaunt = 1. + 0.17282604*pow(etaf,-0.67) - 0.04959570*pow(etaf,-1.33) 
				- 0.01714286*pow(etaf,-2.) + 0.00204498*pow(etaf,-2.67) 
				- 0.00243945*pow(etaf,-3.33) - 0.00120387*pow(etaf,-4.) 
				+ 0.00071814*pow(etaf,-4.67) + 0.00026971*pow(etaf,-5.33);
		}
		else if( zeta > 0.5 )
		{
			/* Grant 1958, as given in Hummer 88, eqn 2.25a	*/
			gaunt = 1. + 0.21775*pow(zeta,-0.67) - 0.01312*pow(zeta, -1.33);
		}
		else 
		{
			double a[10] = {1.47864486, -1.72329012, 0.14420320, 0.05744888, 0.01668957,
				0.01580779,  0.00464268, 0.00385156, 0.00116196, 0.00101967};

			Csum = 0.;
			for( i = 0; i <=9; i++ )
			{
				/* The Chebyshev of the first kind is just a special kind of hypergeometric.	*/
				Csum += a[i]*RealF2_1( (double)(-i), (double)i, 0.5, 0.5*(1.-zeta) );
			}
			gaunt = fabs(0.551329*(0.57721 + log(zeta/2.))*exp(PI*zeta)*Csum);
			ASSERT( gaunt < 10. );
		}
	}
	else if( lgSutherlandOn )
		gaunt = DoSutherland( etai, etaf, chi ); 		
	else
		gaunt = DoBeckert_etal( etai, etaf, chi );

	/*if( gaunt*exp(-x) > 2. && TEglobal < 1000. )
		fprintf( ioQQQ,"ni %.3e nf %.3e chi %.3e u %.3e gam2 %.3e gaunt %.3e x %.3e\n", 
				etai, etaf, chi, TE1RYD*HNUglobal/TEglobal,
				TE1RYD*z*z/TEglobal, gaunt, x);*/

	/** \todo	2	- These are liberal bounds, in final product, this
	 * ASSERT should be much more demanding.	*/
	ASSERT( gaunt > 0. && gaunt<BIGFLOAT );

	if( gaunt == 0. )
	{
		fprintf( ioQQQ, "Uh-Oh! Gaunt is zero!  Is this okay?\n");
		/* assign some small value */
		gaunt = 1e-5;
	}

	return gaunt;
}


#if defined(__ICC) && defined(__i386) && __INTEL_COMPILER < 910
#pragma optimize("", off)
#endif
/************************************
 * This part calculates the gaunt factor as per Beckert et al 2000	*/
/** \todo	2	- insert reference	*/
STATIC double DoBeckert_etal( double etai, double etaf, double chi )
{
	double Delta, BeckertGaunt, MaxFReal, LnBeckertGaunt;
	long NumRenorms[2]={0,0}, NumTerms[2]={0,0};
	int IndexMinNumRenorms, IndexMaxNumRenorms;
	complex<double> a,b,c,F[2];

	a = complex<double>( 1., -etai );
	b = complex<double>( 0., -etaf );
	c = 1.;

	/* evaluate first hypergeometric function.	*/	
	F[0] = Hypergeometric2F1( a, b, c, chi, &NumRenorms[0], &NumTerms[0] );

	a = complex<double>( 1., -etaf );
	b = complex<double>( 0., -etai );

	/* evaluate second hypergeometric function.	*/	
	F[1] = Hypergeometric2F1( a, b, c, chi, &NumRenorms[1], &NumTerms[1] );

	/* If there is a significant difference in the number of terms used, 
	 * they should be recalculated with the max	number of terms in initial calculations */
	/* If NumTerms[i]=-1, the hypergeometric was calculated by the use of an integral instead
	 * of series summation...hence NumTerms has no meaning, and no need to recalculate.	*/
	if( ( MAX2(NumTerms[1],NumTerms[0]) - MIN2(NumTerms[1],NumTerms[0]) >= 2  )
		&& NumTerms[1]!=-1 && NumTerms[0]!=-1)
	{
		a = complex<double>( 1., -etai );
		b = complex<double>( 0., -etaf );
		c = 1.;

		NumTerms[0] = MAX2(NumTerms[1],NumTerms[0])+1;
		NumTerms[1] = NumTerms[0];
		NumRenorms[0] = 0;
		NumRenorms[1] = 0;

		/* evaluate first hypergeometric function.	*/	
		F[0] = Hypergeometric2F1( a, b, c, chi, &NumRenorms[0], &NumTerms[0] );

		a = complex<double>( 1., -etaf );
		b = complex<double>( 0., -etai );

		/* evaluate second hypergeometric function.	*/	
		F[1] = Hypergeometric2F1( a, b, c, chi, &NumRenorms[1], &NumTerms[1] );

		ASSERT( NumTerms[0] == NumTerms[1] );
	}

	/* if magnitude of unNormalized F's are vastly different, zero out the lesser	*/
	if( log10(abs(F[0])/abs(F[1])) + (NumRenorms[0]-NumRenorms[1])*log10(abs(Normalization)) > 10. )
	{
		F[1] = 0.;
		/*  no longer need to keep track of differences in NumRenorms	*/
		NumRenorms[1] = NumRenorms[0];
	}
	else if( log10(abs(F[1])/abs(F[0])) + (NumRenorms[1]-NumRenorms[0])*log10(abs(Normalization)) > 10. )
	{
		F[0] = 0.;
		/*  no longer need to keep track of differences in NumRenorms	*/
		NumRenorms[0] = NumRenorms[1];
	}

	/* now must fix if NumRenorms[0] != NumRenorms[1], because the next calculation is the
	 * difference of squares...information is lost and cannot be recovered if this calculation
	 * is done with NumRenorms[0] != NumRenorms[1]	*/
	MaxFReal = (fabs(F[1].real())>fabs(F[0].real())) ? fabs(F[1].real()):fabs(F[0].real());
	while( NumRenorms[0] != NumRenorms[1] )
	{
		/* but must be very careful to prevent both overflow and underflow.	*/
		if( MaxFReal > 1e50 )
		{
			IndexMinNumRenorms = ( NumRenorms[0] > NumRenorms[1] ) ? 1:0;
			F[IndexMinNumRenorms] /= Normalization;
			++NumRenorms[IndexMinNumRenorms];
		}
		else
		{
			IndexMaxNumRenorms = ( NumRenorms[0] > NumRenorms[1] ) ? 0:1;
			F[IndexMaxNumRenorms] = F[IndexMaxNumRenorms]*Normalization;
			--NumRenorms[IndexMaxNumRenorms];
		}
	}

	ASSERT( NumRenorms[0] == NumRenorms[1] );

	/* Okay, now we are guaranteed (?) a small dynamic range, but may still have to renormalize	*/

	/* Are we gonna have an overflow or underflow problem?	*/
	ASSERT( (fabs(F[0].real())<1e+150) && (fabs(F[1].real())<1e+150) &&
		(fabs(F[0].imag())<1e+150) && (fabs(F[1].real())<1e+150) );
	ASSERT( (fabs(F[0].real())>1e-150) && ((fabs(F[0].imag())>1e-150) || (abs(F[0])==0.)) );
	ASSERT( (fabs(F[1].real())>1e-150) && ((fabs(F[1].real())>1e-150) || (abs(F[1])==0.)) );

	/* guard against spurious underflow/overflow by braindead implementations of the complex class */
	complex<double> CDelta = F[0]*F[0] - F[1]*F[1];
	double renorm = MAX2(fabs(CDelta.real()),fabs(CDelta.imag()));
	ASSERT( renorm > 0. );
	/* carefully avoid complex division here... */
	complex<double> NCDelta( CDelta.real()/renorm, CDelta.imag()/renorm );

	Delta = renorm * abs( NCDelta );

	ASSERT( Delta > 0. );

	/* Now multiply by the coefficient in Beckert 2000, eqn 7	*/
	if( etaf > 100. )
	{
		/* must compute logarithmically	if etaf too big for linear computation.	*/
		LnBeckertGaunt = 1.6940360 + log(Delta) + log(etaf) + log(etai) - log(fabs(etai-etaf)) - 6.2831853*etaf;
		LnBeckertGaunt += 2. * NumRenorms[0] * log(abs(Normalization));
		BeckertGaunt = exp( LnBeckertGaunt );
		NumRenorms[0] = 0;
	}
	else
	{
		BeckertGaunt = Delta*5.4413981*etaf*etai/fabs(etai - etaf)
			/(1.-exp(-6.2831853*etai) )/( exp(6.2831853*etaf) - 1.);

		while( NumRenorms[0] > 0 )
		{
			BeckertGaunt *= abs(Normalization);
			BeckertGaunt *= abs(Normalization);
			ASSERT( BeckertGaunt < BIGDOUBLE );
			--NumRenorms[0];
		} 
	}

	ASSERT( NumRenorms[0] == 0 );

	/*fprintf( ioQQQ,"etai %.3e etaf %.3e u %.3e B %.3e \n", 
		etai, etaf, TE1RYD * HNUglobal / TEglobal, BeckertGaunt );	*/

	return BeckertGaunt;
}
#if defined(__ICC) && defined(__i386) && __INTEL_COMPILER < 910
#pragma optimize("", on)
#endif


/************************************
 * This part calculates the gaunt factor as per Sutherland 98	*/
/** \todo	2	- insert reference	*/
STATIC double DoSutherland( double etai, double etaf, double chi )
{
	double Sgaunt, ICoef, weightI1, weightI0;
	long i, NumRenorms[2]={0,0}, NumTerms[2]={0,0};
	complex<double> a,b,c,GCoef,kfac,etasum,G[2],I[2],ComplexFactors,GammaProduct;

	kfac = complex<double>( fabs((etaf-etai)/(etaf+etai)), 0. );
	etasum = complex<double>( 0., etai + etaf );

	GCoef = pow(kfac, etasum);
	/* GCoef is a complex vector that should be contained within the unit circle.
	 * and have a non-zero magnitude.  Or is it ON the unit circle?	*/
	ASSERT( fabs(GCoef.real())<1.0 && fabs(GCoef.imag())<1.0 && ( GCoef.real()!=0. || GCoef.imag()!=0. ) );

	for( i = 0; i <= 1; i++ )
	{
		a = complex<double>( i + 1., -etaf );
		b = complex<double>( i + 1., -etai );
		c = complex<double>( 2.*i + 2., 0. );

		/* First evaluate hypergeometric function.	*/	
		G[i] = Hypergeometric2F1( a, b, c, chi, &NumRenorms[i], &NumTerms[i] );
	}

	/* If there is a significant difference in the number of terms used, 
	 * they should be recalculated with the max	number of terms in initial calculations */
	/* If NumTerms[i]=-1, the hypergeometric was calculated by the use of an integral instead
	 * of series summation...hence NumTerms has no meaning, and no need to recalculate.	*/
	if( MAX2(NumTerms[1],NumTerms[0]) - MIN2(NumTerms[1],NumTerms[0]) > 2 
		&& NumTerms[1]!=-1 && NumTerms[0]!=-1  )
	{
		NumTerms[0] = MAX2(NumTerms[1],NumTerms[0]);
		NumTerms[1] = NumTerms[0];
		NumRenorms[0] = 0;
		NumRenorms[1] = 0;

		for( i = 0; i <= 1; i++ )
		{
			a = complex<double>( i + 1., -etaf );
			b = complex<double>( i + 1., -etai );
			c = complex<double>( 2.*i + 2., 0. );

			G[i] = Hypergeometric2F1( a, b, c, chi, &NumRenorms[i], &NumTerms[i] );
		}

		ASSERT( NumTerms[0] == NumTerms[1] );
	}

	for( i = 0; i <= 1; i++ )
	{
		/** \todo	2	- this check may also too liberal.  */
		ASSERT( fabs(G[i].real())>0. && fabs(G[i].real())<1e100 &&
			fabs(G[i].imag())>0. && fabs(G[i].imag())<1e100 );

		/* Now multiply by the coefficient in Sutherland 98, eqn 9	*/
		G[i] *= GCoef;

		/* This is the coefficient in equation 8 in Sutherland	*/
		/* Karzas and Latter give tgamma(2.*i+2.), Sutherland gives tgamma(2.*i+1.)	*/
		ICoef = 0.25*pow(-chi, (double)i+1.)*exp( 1.5708*fabs(etai-etaf) )/factorial(2*i);
		GammaProduct = cdgamma(complex<double>(i+1.,etai))*cdgamma(complex<double>(i+1.,etaf));
		ICoef *= abs(GammaProduct);

		ASSERT( ICoef > 0. );

		I[i] = ICoef*G[i];

		while( NumRenorms[i] > 0 )
		{
			I[i] *= Normalization;
			ASSERT( fabs(I[i].real()) < BIGDOUBLE && fabs(I[i].imag()) < BIGDOUBLE );
			--NumRenorms[i];
		} 

		ASSERT( NumRenorms[i] == 0 );
	}

	weightI0 = POW2(etaf+etai);
	weightI1 = 2.*etaf*etai*sqrt(1. + etai*etai)*sqrt(1. + etaf*etaf);

	ComplexFactors = I[0] * ( weightI0*I[0] - weightI1*I[1] );

	/* This is Sutherland equation 13	*/
	Sgaunt  = 1.10266 / etai / etaf * abs( ComplexFactors );

	return Sgaunt;
}

#if defined(__ICC) && defined(__i386) && __INTEL_COMPILER < 910
#pragma optimize("", off)
#endif
/* This routine is a wrapper for F2_1	*/
STATIC complex<double> Hypergeometric2F1( complex<double> a, complex<double> b, complex<double> c,
					  double chi, long *NumRenorms, long *NumTerms )
{
	complex<double> a1, b1, c1, a2, b2, c2, Result, Part[2], F[2];
	complex<double> chifac, GammaProduct, Coef, FIntegral;
	/** \todo	2	- pick these interface values and stick with it...best results have been 0.4, 1.5  */
	double Interface1 = 0.4, Interface2 = 10.;
	long N_Renorms[2], N_Terms[2], IndexMaxNumRenorms, lgDoIntegral = false;

	N_Renorms[0] = *NumRenorms;
	N_Renorms[1] = *NumRenorms;
	N_Terms[0] = *NumTerms;
	N_Terms[1] = *NumTerms;

	/* positive and zero chi are not possible.	*/
	ASSERT( chi < 0. );

	/* We want to be careful about evaluating the hypergeometric 
	 * in the vicinity of chi=1.  So we employ three different methods...*/

	/* for small chi, we pass the parameters to the hypergeometric function as is.	*/
	if( fabs(chi) < Interface1 )
	{
		Result = F2_1( a, b, c, chi, &*NumRenorms, &*NumTerms );
	}
	/* for large chi, we use a relation given as eqn 5 in Nicholas 89.	*/
	else if( fabs(chi) > Interface2 )
	{
		a1 = a;
		b1 = 1.-c+a;
		c1 = 1.-b+a;

		a2 = b;
		b2 = 1.-c+b;
		c2 = 1.-a+b;

		chifac = -chi;

		F[0] = F2_1(a1,b1,c1,1./chi,&N_Renorms[0], &N_Terms[0]);
		F[1] = F2_1(a2,b2,c2,1./chi,&N_Renorms[1], &N_Terms[1]);

		/* do it again if significant difference in number of terms.	*/
		if( MAX2(N_Terms[1],N_Terms[0]) - MIN2(N_Terms[1],N_Terms[0]) >= 2 )
		{
			N_Terms[0] = MAX2(N_Terms[1],N_Terms[0]);
			N_Terms[1] = N_Terms[0];
			N_Renorms[0] = *NumRenorms;
			N_Renorms[1] = *NumRenorms;

			F[0] = F2_1(a1,b1,c1,1./chi,&N_Renorms[0], &N_Terms[0]);
			F[1] = F2_1(a2,b2,c2,1./chi,&N_Renorms[1], &N_Terms[1]);
			ASSERT( N_Terms[0] == N_Terms[1] );
		}

		*NumTerms = MAX2(N_Terms[1],N_Terms[0]);

		/************************************************************************/
		/* Do the first part	*/
		GammaProduct = (cdgamma(b-a)/cdgamma(b))*(cdgamma(c)/cdgamma(c-a));

		/* divide the hypergeometric by (-chi)^a and multiply by GammaProduct	*/
		Part[0] = F[0]/pow(chifac,a)*GammaProduct;

		/************************************************************************/
		/* Do the second part	*/
		GammaProduct = (cdgamma(a-b)/cdgamma(a))*(cdgamma(c)/cdgamma(c-b));

		/* divide the hypergeometric by (-chi)^b and multiply by GammaProduct	*/
		Part[1] = F[1]/pow(chifac,b)*GammaProduct;

		/************************************************************************/
		/* Add the two parts to get the result.	*/

		/* First must fix it if N_Renorms[0] != N_Renorms[1]	*/
		if( N_Renorms[0] != N_Renorms[1] )
		{
			IndexMaxNumRenorms = ( N_Renorms[0] > N_Renorms[1] ) ? 0:1;
			Part[IndexMaxNumRenorms] *= Normalization;
			--N_Renorms[IndexMaxNumRenorms];
			/* Only allow at most a difference of one in number of renormalizations...
			 * otherwise something is really screwed up	*/
			ASSERT( N_Renorms[0] == N_Renorms[1] );
		}

		*NumRenorms = N_Renorms[0];

		Result = Part[0] + Part[1];
	}
	/* And for chi of order 1, we use Nicholas 89, eqn 27.	*/
	else
	{
		/* the hypergeometric integral does not seem to work well.	*/
		if( lgDoIntegral /* && fabs(chi+1.)>0.1 */)
		{
			/* a and b are always interchangeable, assign the lesser to b to 
			 * prevent Coef from blowing up	*/
			if( abs(b) > abs(a) )
			{
				complex<double> btemp = b;
				b = a;
				a = btemp;
			}
			Coef = cdgamma(c)/(cdgamma(b)*cdgamma(c-b));
			CMinusBMinus1 = c-b-1.;
			BMinus1 = b-1.;
			MinusA = -a;
			GlobalCHI = chi;
			FIntegral = qg32complex( 0., 0.5, HyperGeoInt );
			FIntegral += qg32complex( 0.5, 1., HyperGeoInt );

			Result = Coef + FIntegral;
			*NumTerms = -1;
			*NumRenorms = 0;
		}
		else
		{
			/*	Near chi=1 solution	*/
			a1 = a;
			b1 = c-b;
			c1 = c;
			chifac = 1.-chi;

			Result = F2_1(a1,b1,c1,chi/(chi-1.),&*NumRenorms,&*NumTerms)/pow(chifac,a);
		}
	}

	/* Limit the size of the returned value	*/
	while( fabs(Result.real()) >= 1e50 )
	{
		Result /= Normalization;
		++*NumRenorms;
	}

	return Result;
}

/* This routine calculates hypergeometric functions */
STATIC complex<double> F2_1(
	complex<double> alpha, complex<double> beta, complex<double> gamma,
	double chi, long *NumRenormalizations, long *NumTerms )
{
	long  i = 3, MinTerms;
	bool lgNotConverged = true;
	complex<double> LastTerm, Term, Sum;

	MinTerms = MAX2( 3, *NumTerms );

	/* This is the first term of the hypergeometric series.	*/
	Sum = 1./Normalization;
	++*NumRenormalizations;

	/* This is the second term	*/
	LastTerm = Sum*alpha*beta/gamma*chi;

	Sum += LastTerm;

	/* Every successive term is easily found by multiplying the last term
	 * by (alpha + i - 2)*(beta + i - 2)*chi/(gamma + i - 2)/(i-1.)	*/
	do{
		alpha += 1.;
		beta += 1.;
		gamma += 1.;

		/* multiply old term by incremented alpha++*beta++/gamma++.  Also multiply by chi/(i-1.)	*/
		Term = LastTerm*alpha*beta/gamma*chi/(i-1.); 

		Sum += Term;

		/* Renormalize if too big	*/
		if( Sum.real() > 1e100 )
		{
			Sum /= Normalization;
			LastTerm = Term/Normalization;
			++*NumRenormalizations;
			/* notify of renormalization, and print the number of the term	*/
			fprintf( ioQQQ,"Hypergeometric: Renormalized at term %li.  Sum = %.3e %.3e\n",
				i, Sum.real(), Sum.imag());
		}
		else
			LastTerm = Term;

		/* Declare converged if this term does not affect Sum by much	*/
		/* Must do this with abs because terms alternate sign.	*/
		if( fabs(LastTerm.real()/Sum.real())<0.001 && fabs(LastTerm.imag()/Sum.imag())<0.001 )
			lgNotConverged = false;

		if( *NumRenormalizations >= 5 )
		{
			fprintf( ioQQQ, "We've got too many (%li) renorms!\n",*NumRenormalizations );
		}

		++i;

	}while ( lgNotConverged || i<MinTerms );

	*NumTerms = i;

	return Sum;
}
#if defined(__ICC) && defined(__i386) && __INTEL_COMPILER < 910
#pragma optimize("", on)
#endif

/* This routine calculates hypergeometric functions */
STATIC double RealF2_1( double alpha, double beta, double gamma, double chi )
{
	long  i = 3;
	bool lgNotConverged = true;
	double LastTerm, Sum;

	/* This is the first term of the hypergeometric series.	*/
	Sum = 1.;

	/* This is the second term	*/
	LastTerm = alpha*beta*chi/gamma;

	Sum += LastTerm;

	/* Every successive term is easily found by multiplying the last term
	 * by (alpha + i - 2)*(beta + i - 2)*chi/(gamma + i - 2)/(i-1.)	*/
	do{
		alpha++;
		beta++;
		gamma++;

		/* multiply old term by incremented alpha++*beta++/gamma++.  Also multiply by chi/(i-1.)	*/
		LastTerm *= alpha*beta*chi/gamma/(i-1.); 

		Sum += LastTerm;

		/* Declare converged if this term does not affect Sum by much	*/
		/* Must do this with abs because terms alternate sign.	*/
		if( fabs(LastTerm/Sum)<0.001 )
			lgNotConverged = false;

		++i;

	}while ( lgNotConverged );

	return Sum;
}

STATIC complex<double> HyperGeoInt( double v )
{
	return pow(v,BMinus1)*pow(1.-v,CMinusBMinus1)*pow(1.-v*GlobalCHI,MinusA);
}

/*complex 32 point Gaussian quadrature, originally given to Gary F by Jim Lattimer */
/* modified to handle complex numbers by Ryan Porter.	*/
STATIC complex<double> qg32complex(
	double xl, /*lower limit to integration range*/
	double xu, /*upper limit to integration range*/
	/*following is the pointer to the function that will be evaulated*/
	complex<double> (*fct)(double) )
{
	double a, 
	  b, 
	  c;
	complex<double> y;


	/********************************************************************************
	 *                                                                              *
	 *  32-point Gaussian quadrature                                                *
	 *  xl  : the lower limit of integration                                        *
	 *  xu  : the upper limit                                                       *
	 *  fct : the (external) function                                               *
	 *  returns the value of the integral                                           *
	 *                                                                              *
	 * simple call to integrate sine from 0 to pi                                   *
	 * double agn = qg32( 0., 3.141592654 ,  sin );                                 *
	 *                                                                              *
	 *******************************************************************************/

	a = .5*(xu + xl);
	b = xu - xl;
	c = .498631930924740780*b;
	y = .35093050047350483e-2 * ( (*fct)(a+c) + (*fct)(a-c) );
	c = b*.49280575577263417;
	y += .8137197365452835e-2 * ( (*fct)(a+c) + (*fct)(a-c) );
	c = b*.48238112779375322;
	y += .1269603265463103e-1 * ( (*fct)(a+c) + (*fct)(a-c) );
	c = b*.46745303796886984;
	y += .17136931456510717e-1* ( (*fct)(a+c) + (*fct)(a-c) );
	c = b*.44816057788302606;
	y += .21417949011113340e-1* ( (*fct)(a+c) + (*fct)(a-c) );
	c = b*.42468380686628499;
	y += .25499029631188088e-1* ( (*fct)(a+c) + (*fct)(a-c) );
	c = b*.3972418979839712;
	y += .29342046739267774e-1* ( (*fct)(a+c) + (*fct)(a-c) );
	c = b*.36609105937014484;
	y += .32911111388180923e-1* ( (*fct)(a+c) + (*fct)(a-c) );
	c = b*.3315221334651076;
	y += .36172897054424253e-1* ( (*fct)(a+c) + (*fct)(a-c) );
	c = b*.29385787862038116;
	y += .39096947893535153e-1* ( (*fct)(a+c) + (*fct)(a-c) );
	c = b*.2534499544661147;
	y += .41655962113473378e-1* ( (*fct)(a+c) + (*fct)(a-c) );
	c = b*.21067563806531767;
	y += .43826046502201906e-1* ( (*fct)(a+c) + (*fct)(a-c) );
	c = b*.16593430114106382;
	y += .45586939347881942e-1* ( (*fct)(a+c) + (*fct)(a-c) );
	c = b*.11964368112606854;
	y += .46922199540402283e-1* ( (*fct)(a+c) + (*fct)(a-c) );
	c = b*.7223598079139825e-1;
	y += .47819360039637430e-1* ( (*fct)(a+c) + (*fct)(a-c) );
	c = b*.24153832843869158e-1;
	y += .4827004425736390e-1 * ( (*fct)(a+c) + (*fct)(a-c) );
	y *= b;

	/* the answer */

	return( y );
}
