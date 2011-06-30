/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* CS_VS80 compute thermally averaged collision strength for collisional deexcitation 
 * of hydrogenic atoms, from
 * >>refer	H1	collision	Vriens, L., & Smeets, A.H.M. 1980, Phys Rev A 22, 940
 *hydro_vs_ioniz generate hydrogenic collisional ionization rate coefficients */
/*Hion_coll_ioniz_ratecoef calculate hydrogenic ionization rates for ions 
 * with all n, and Z*/
#include "cddefines.h"
#include "dense.h"
#include "phycon.h"
#include "physconst.h"
#include "iso.h"
#include "hydro_vs_rates.h"
#include "lines_service.h"
#include "taulines.h"

STATIC double hydro_vs_coll_str( double energy, long ipISO, long nelem, long ipHi, long ipLo, long Collider, double Aul );

namespace {
	class my_Integrand : public std::unary_function<double, double> 
	{
	public:
		long ipISO, nelem, ipHi, ipLo, Collider;
		double Aul;
		double temp;
		
		double operator() (double EOverKT)
			{
				double col_str = hydro_vs_coll_str( EOverKT * EVRYD * temp/TE1RYD, ipISO, nelem, ipHi, ipLo, Collider, Aul );
				return exp( -1.*EOverKT ) * col_str;
			}
	};
}

/* These are masses relative to the proton mass of the electron, proton, and alpha particle. */
static const double ColliderMass[4] = {ELECTRON_MASS/PROTON_MASS, 1.0, 4.0, 4.0};

/*
 Neither of these rates can be modified to account for impact by non-electrons because they 
 are fits to thermally-averaged rates for electron impact...it can not be unravelled to 
 expose a projectile energy that can then be scaled to account for other projectiles.
 Instead, we have included the original cross section formula (eq 14) from 
 Vriens & Smeets (1980) below.
*/

/* VS80 stands for Vriens and Smeets 1980 */
/* This routine calculates thermally-averaged collision strengths. */
double CS_VS80( long int ipISO, long int nelem, long int ipHi, long int ipLo, double Aul, double temp, long int Collider )
{
	double coll_str;

	if( Collider == ipELECTRON )
	{
		coll_str = hydro_vs_deexcit( ipISO, nelem, ipHi, ipLo, Aul);
	}
	else
	{
		/* evaluate collision strength, two options,
		 * do thermal average (very slow) if set with
		 * SET COLLISION STRENGTH AVERAGE command,
		 * default just do point at kT
		 * tests show no impact on test suite, much faster */
		if( iso.lgCollStrenThermAver )
		{
			my_Integrand func;

			func.ipISO = ipISO;
			func.nelem = nelem;
			func.ipHi = ipHi;
			func.ipLo = ipLo;
			func.temp = temp;
			func.Collider = Collider;
			func.Aul = Aul;

			Integrator<my_Integrand,Gaussian32> VS80;
			/* do average over Maxwellian */
			coll_str = VS80.sum( 0., 1., func );
			coll_str += VS80.sum( 1., 10., func );
		}
		else
		{
			/* the argument to the function is equivalent to evaluating the function at
			 * hnu = kT */
			coll_str = hydro_vs_coll_str( EVRYD*temp/TE1RYD, ipISO, nelem, ipHi, ipLo, Collider, Aul );
		}
	}

	ASSERT( coll_str >= 0. );
	return coll_str;
}

/*hydro_vs_deexcit compute collision strength for collisional deexcitation for hydrogen atom, 
 * from Vriens and Smeets */
STATIC double hydro_vs_coll_str( double energy, long ipISO, long nelem, long ipHi, long ipLo, long Collider, double Aul )
{
	double Apn, bp, Bpn, delta;
	double Epi, Epn;
	double gamma, p, n;
	double ryd, s;
	double cross_section, coll_str, gLo, gHi, abs_osc_str, reduced_mass;

	DEBUG_ENTRY( "hydro_vs_coll_str()" );

	// number of colliders is 4 in def, should be macro
	ASSERT( Collider >= 0.&& Collider <4 );
	reduced_mass = dense.AtomicWeight[nelem]*ColliderMass[Collider]/
		(dense.AtomicWeight[nelem]+ColliderMass[Collider])*ATOMIC_MASS_UNIT;

	gLo = StatesElemNEW[nelem][nelem-ipISO][ipLo].g;
	gHi = StatesElemNEW[nelem][nelem-ipISO][ipHi].g;

	/* This comes from equations 14,15, and 16 of Vriens and Smeets. */
	/* >>refer he-like cs Vriens, L. \& Smeets, A. H. M. Phys. Rev. A, 22, 940 */ 
	/* Computes the Vriens and Smeets
	 * rate coeff. for collisional dexcitation between any two levels of H.
	 * valid for all nhigh, nlow
	 * at end converts this excitation rate to collision strength	*/

	n = (double)StatesElemNEW[nelem][nelem-ipISO][ipHi].n;
	p = (double)StatesElemNEW[nelem][nelem-ipISO][ipLo].n;

	ryd = EVRYD;
	s = fabs(n-p);
	ASSERT( s > 0. );

	Epn = EVRYD*(iso.xIsoLevNIonRyd[ipISO][nelem][ipLo] - iso.xIsoLevNIonRyd[ipISO][nelem][ipHi]);
	Epi = EVRYD*iso.xIsoLevNIonRyd[ipISO][nelem][ipLo];

	/* This is an absorption oscillator strength. */
	abs_osc_str = GetGF( Aul, Epn*RYD_INF/EVRYD, gHi)/gLo;
	Apn = 2.*ryd/Epn*abs_osc_str;

	bp = 1.4*log(p)/p - .7/p - .51/p/p + 1.16/p/p/p - 0.55/p/p/p/p;

	Bpn = 4.*ryd*ryd/n/n/n*(1./Epn/Epn + 4./3.*Epi/POW3(Epn) + bp*Epi*Epi/powi(Epn,4));

	delta = exp(-Bpn/Apn) - 0.4*Epn/ryd;

	gamma = ryd*(8. + 23.*POW2(s/p))/
		( 8. + 1.1*n*s + 0.8/pow2(s) + 0.4*sqrt(pow3(n))/sqrt(s)*fabs(s-1.0) );

	/* Scale the energy to get an equivalent electron energy. */
	energy *= ColliderMass[ipELECTRON]/ColliderMass[Collider];

	/* cross section in units of PI*a_o^2 */
	if( energy/2./ryd+delta <= 0 /*|| energy < Epn*/ )
	{
		cross_section = 0.;
	}
	else
	{
		cross_section = 2.*ryd/(energy + gamma)*(Apn*log(energy/2./ryd+delta) + Bpn);
		cross_section = MAX2( cross_section, 0. );
	}

	/* convert to collision strength */
	coll_str = ConvCrossSect2CollStr( cross_section * PI*BOHR_RADIUS_CM*BOHR_RADIUS_CM, gLo, energy/EVRYD, reduced_mass );

	ASSERT( coll_str >= 0. );

	return( coll_str );
}

/*hydro_vs_coll_recomb generate hydrogenic collisional recombination rate coefficients */
double hydro_vs_coll_recomb( double ionization_energy_Ryd, double Te, double stat_level, double stat_ion )
{
	double coef, 
	  denom, 
	  epi, 
	  t_eV;

	DEBUG_ENTRY( "hydro_vs_coll_recomb()" );

	/* This is equation 9 of
	 * >>refer	H1	collision recomb	Vriens, L., & Smeets, A.H.M. 1980, Phys Rev A 22, 940 */

	/* this is kT in eV */
	t_eV = Te/EVDEGK;

	/* this is the ionization energy relative to kT, dimensionless */
	epi = ionization_energy_Ryd * EVRYD / t_eV;

	/* this is the denominator of equation 8 of VS80. */
	denom = pow(epi,2.33) + 4.38*pow(epi,1.72) + 1.32*epi;

	/* this is equation 9 of VS80 */
	coef = 3.17e-27 / pow3(t_eV) * stat_level / stat_ion / denom;

	ASSERT( coef >= 0. );
	return( coef );
}


/*hydro_vs_ioniz generate hydrogenic collisional ionization rate coefficients */
double hydro_vs_ioniz( double ionization_energy_Ryd, double Te )
{
	double coef, 
	  denom, 
	  epi, 
	  t_eV;

	DEBUG_ENTRY( "hydro_vs_ioniz()" );

	/* a function written to calculate the rate coefficients 
	 * for hydrogen collisional ionizations from
	 * Jason Ferguson, summer 94
	 * valid for all n
	 * >>refer	H1	collision	Vriens, L., & Smeets, A.H.M. 1980, Phys Rev A 22, 940
	 * */

	/* this is kT in eV */
	t_eV = Te/EVDEGK;

	/* this is the ionization energy relative to kT, dimensionless */
	epi = ionization_energy_Ryd * EVRYD / t_eV;

	/* this is the denominator of equation 8 of VS80. */
	denom = pow(epi,2.33) + 4.38*pow(epi,1.72) + 1.32*epi;

	/* this is equation 8 of VS80 */
	coef = 9.56e-6 / sqrt(pow3(t_eV)) * dsexp( epi ) / denom;

	ASSERT( coef >= 0. );
	return( coef );
}

/*Hion_coll_ioniz_ratecoef calculate hydrogenic ionization rates for all n, and Z*/
double Hion_coll_ioniz_ratecoef(
		/* the isoelectronic sequence */
		long int ipISO ,
		/* element, >=1 since only used for ions 
		 * nelem = 1 is helium the least possible charge */
		long int nelem,
		/* principal quantum number, > 1
		 * since only used for excited states */
		long int n,
		double ionization_energy_Ryd,
		double Te )
{
	double H, 
	  HydColIon_v, 
	  Rnp, 
	  charge,
	  chim, 
	  eone, 
	  etwo, 
	  ethree, 
	  g, 
	  rate, 
	  rate2, 
	  boltz,
	  t1, 
	  t2, 
	  t3, 
	  t4, 
	  tev, 
	  xn, 
	  y;
	static const double arrH[4] = {1.48,3.64,5.93,8.32};
	static const double arrRnp[8] = {2.20,1.90,1.73,1.65,1.60,1.56,1.54,1.52};
	static const double arrg[10] = {0.8675,0.932,0.952,0.960,0.965,0.969,0.972,0.975,0.978,0.981};

	static double small = 0.;

	DEBUG_ENTRY( "Hion_coll_ioniz_ratecoef()" );
	/*calculate hydrogenic ionization rates for all n, and Z
	 * >>refer	HI	cs	Allen 1973, Astro. Quan. for low Te.
	 * >>refer	HI	cs	Sampson and Zhang 1988, ApJ, 335, 516 for High Te.
	 * */

	charge = nelem - ipISO;
	/* this routine only for ions, nelem=0 is H, nelem=1 he, etc */
	ASSERT( charge > 0);
	ASSERT( n>1 );

	if( n > 4 )
	{
		H = 2.15*n;
	}
	else
	{
		H = arrH[n-1];
	}

	if( n > 8 )
	{
		Rnp = 1.52;
	}
	else
	{
		Rnp = arrRnp[n-1];
	}

	if( n > 10 )
	{
		g = arrg[9];
	}
	else
	{
		g = arrg[n-1];
	}

	tev = EVRYD/TE1RYD*Te;
	xn = (double)n;
	chim = EVRYD * ionization_energy_Ryd;
	y = chim/tev;
	boltz = dsexp( chim/tev );

	eone = ee1(y);
	etwo = boltz - y*eone;
	ethree = (boltz - y*etwo)/2.;

	t1 = 1/xn*eone;
	t2 = 1./3./xn*(boltz - y*ethree);
	t3 = 3.*H/xn/(3. - Rnp)*(y*etwo - 2.*y*eone + boltz);
	t4 = 3.36*y*(eone - etwo);
	rate = 7.69415e-9*sqrt(Te)*9.28278e-3*powi(xn/(charge+1),4)*g*y;
	rate *= t1 - t2 + t3 + t4;
	rate2 = 2.1e-8*sqrt(Te)/chim/chim*dsexp(2.302585*5040.*chim/Te);

	/* don't let the rates go negative */
	rate = MAX2(rate,small);
	rate2 = MAX2(rate2,small);

	/* Take the lowest of the two, they fit nicely together... */
	/* >>chng 10 Sept 02, sometimes one of these is zero and the other is positive.
	 * in that case take the bigger one.	*/
	if( rate==0. || rate2==0. )
		HydColIon_v = MAX2(rate,rate2);
	else
		HydColIon_v = MIN2(rate,rate2);

	ASSERT( HydColIon_v >= 0. );
	return( HydColIon_v );
}

/*hydro_vs_deexcit compute collisional deexcitation rates for hydrogen atom, 
 * from Vriens and Smeets 1980 */
double hydro_vs_deexcit( long ipISO, long nelem, long ipHi, long ipLo, double Aul )
{
	double Anp, bn, Bnp, delta_np;
	double Eni, Enp;
	double Gamma_np, p, n, g_p, g_n;
	double ryd, s, kT_eV, rate, col_str, abs_osc_str;

	DEBUG_ENTRY( "hydro_vs_deexcit()" );

	kT_eV = EVRYD * phycon.te / TE1RYD;

	/* This comes from equations 24 of Vriens and Smeets. */
	/* >>refer he-like cs Vriens, L. \& Smeets, A. H. M. Phys. Rev. A, 22, 940 */ 
	/* Computes the Vriens and Smeets
	 * rate coeff. for collisional dexcitation between any two levels of H.
	 * valid for all nhigh, nlow
	 * at end converts this excitation rate to collision strength	*/

	n = (double)StatesElemNEW[nelem][nelem-ipISO][ipLo].n;
	p = (double)StatesElemNEW[nelem][nelem-ipISO][ipHi].n;

	ASSERT( n!=p );

	g_n = StatesElemNEW[nelem][nelem-ipISO][ipLo].g;
	g_p = StatesElemNEW[nelem][nelem-ipISO][ipHi].g;

	ryd = EVRYD;
	s = fabs(n-p);

	Enp = EVRYD*(iso.xIsoLevNIonRyd[ipISO][nelem][ipLo] - iso.xIsoLevNIonRyd[ipISO][nelem][ipHi]);
	Eni = EVRYD*iso.xIsoLevNIonRyd[ipISO][nelem][ipHi];

	ASSERT( Enp > 0. );

	/* This is an absorption oscillator strength. */
	abs_osc_str = GetGF( Aul, Enp*RYD_INF/EVRYD, g_p)/g_n;
	Anp = 2.*ryd/Enp*abs_osc_str;

	bn = 1.4*log(n)/n - .7/n - .51/n/n + 1.16/n/n/n - 0.55/n/n/n/n;

	Bnp = 4.*ryd*ryd/p/p/p*(1./Enp/Enp + 4./3.*Eni/POW3(Enp) + bn*Eni*Eni/powi(Enp,4));

	delta_np = exp(-Bnp/Anp) + 0.1*Enp/ryd;

	Gamma_np = ryd*log(1. + n*n*n*kT_eV/ryd) * (3. + 11.*s*s/n/n) /
		( 6. + 1.6*p*s + 0.3/pow2(s) + 0.8*sqrt(pow3(p))/sqrt(s)*fabs(s-0.6) );

	if( 0.3*kT_eV/ryd+delta_np <= 0 )
	{
		rate = 0.;
	}
	else
	{
		rate = 1.6E-7 * sqrt(kT_eV) * g_n / g_p / ( kT_eV + Gamma_np ) *
			( Anp * log(0.3*kT_eV/ryd + delta_np) + Bnp );
	}

	col_str = rate / COLL_CONST * phycon.sqrte * StatesElemNEW[nelem][nelem-ipISO][ipHi].g;

	return( col_str );
}

