/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*atom_pop5 do five level atom population and cooling */
#include "cddefines.h"
#include "physconst.h"
#include "phycon.h"
#include "thermal.h"
#include "dense.h"
#include "thirdparty.h"
#include "atoms.h"

/*atom_pop5 do five level atom population and cooling */
void atom_pop5(
	/* vector giving statistical weights on the five levels */
	const double g[], 
	/* vector giving the excitation energy differences of the 5 levels.  The energies
	 * are the energy in wavenumbers between adjacent levels.  So EnerWN[0] is the energy
	 * 1-2, EnerWN[1] is the energy between 2-3, etc */
	const double EnerWN[], 
	/* the collision strengths for the levels */
	double cs12, 
	double cs13, 
	double cs14, 
	double cs15, 
	double cs23, 
	double cs24, 
	double cs25, 
	double cs34, 
	double cs35, 
	double cs45, 
	/* the transition probabilities between the various levels */
	double a21, 
	double a31, 
	double a41, 
	double a51, 
	double a32, 
	double a42, 
	double a52, 
	double a43, 
	double a53, 
	double a54, 
	/* the destroyed level populations (units cm-3) for the five levels */
	double p[], 
	/* the total density of this ion */
	realnum abund,
	double *Cooling, 
	double *CoolingDeriv,
	double pump01,
	double pump02,
	double pump03,
	double pump04
	)
{

	DEBUG_ENTRY( "atom_pop5()" );

	/* quit if no species present */
	ASSERT( abund>=0. );
	if( abund == 0. )
	{
		p[0] = 0.;
		p[1] = 0.;
		p[2] = 0.;
		p[3] = 0.;
		p[4] = 0.;
		*Cooling = 0.;
		*CoolingDeriv = 0.;
		return;
	}

	// tf = 1.438 / te, converts energy in wavenumbers into Boltzmann factor
	double tf = T1CM/phycon.te;

	// define Boltzmann factors
	double BoltzFac[5][5];

	BoltzFac[0][1] = sexp(EnerWN[0]*tf);
	BoltzFac[1][2] = sexp(EnerWN[1]*tf);
	BoltzFac[2][3] = sexp(EnerWN[2]*tf);
	BoltzFac[3][4] = sexp(EnerWN[3]*tf);
	BoltzFac[0][2] = BoltzFac[0][1]*BoltzFac[1][2];
	BoltzFac[0][3] = BoltzFac[0][2]*BoltzFac[2][3];
	BoltzFac[0][4] = BoltzFac[0][3]*BoltzFac[3][4];
	BoltzFac[1][3] = BoltzFac[1][2]*BoltzFac[2][3];
	BoltzFac[1][4] = BoltzFac[1][3]*BoltzFac[3][4];
	BoltzFac[2][4] = BoltzFac[2][3]*BoltzFac[3][4];

	/* quit it highest level Boltzmann factor too large */
	if( (BoltzFac[0][4]+pump04) == 0. )
	{
		p[0] = 0.;
		p[1] = 0.;
		p[2] = 0.;
		p[3] = 0.;
		p[4] = 0.;
		*Cooling = 0.;
		*CoolingDeriv = 0.;
		return;
	}

	// get collision rates, dense.cdsqte is 8.629e-6 / sqrte * eden */
	// rates units are s^-1
	double col[5][5];
	col[1][0] = dense.cdsqte*cs12/g[1];
	col[0][1] = col[1][0]*g[1]/g[0]*BoltzFac[0][1];

	col[2][0] = dense.cdsqte*cs13/g[2];
	col[0][2] = col[2][0]*g[2]/g[0]*BoltzFac[0][2];

	col[3][0] = dense.cdsqte*cs14/g[3];
	col[0][3] = col[3][0]*g[3]/g[0]*BoltzFac[0][3];

	col[4][0] = dense.cdsqte*cs15/g[4];
	col[0][4] = col[4][0]*g[4]/g[0]*BoltzFac[0][4];

	col[2][1] = dense.cdsqte*cs23/g[2];
	col[1][2] = col[2][1]*g[2]/g[1]*BoltzFac[1][2];

	col[3][1] = dense.cdsqte*cs24/g[3];
	col[1][3] = col[3][1]*g[3]/g[1]*BoltzFac[1][3];

	col[4][1] = dense.cdsqte*cs25/g[4];
	col[1][4] = col[4][1]*g[4]/g[1]*BoltzFac[1][4];

	col[3][2] = dense.cdsqte*cs34/g[3];
	col[2][3] = col[3][2]*g[3]/g[2]*BoltzFac[2][3];

	col[4][2] = dense.cdsqte*cs35/g[4];
	col[2][4] = col[4][2]*g[4]/g[2]*BoltzFac[2][4];

	col[4][3] = dense.cdsqte*cs45/g[4];
	col[3][4] = col[4][3]*g[4]/g[3]*BoltzFac[3][4];

	double amat[5][5], bvec[5];
	// homogeneous matrix - no source or sink - use conservation at 5th equation
	for( long i=0; i<5; ++i )
	{
		amat[i][4] = 1.;
		bvec[i] = 0.;
	}
	bvec[4] = abund;

	/* level one balance equation */
	amat[0][0] = col[0][1] + col[0][2] + col[0][3] + col[0][4] +pump01+pump02+pump03+pump04;
	amat[1][0] = -a21 - col[1][0];
	amat[2][0] = -a31 - col[2][0];
	amat[3][0] = -a41 - col[3][0];
	amat[4][0] = -a51 - col[4][0];

	/* level two balance equation */
	amat[0][1] = -col[0][1] - pump01;
	amat[1][1] = col[1][0] + a21 + col[1][2] + col[1][3] + col[1][4];
	amat[2][1] = -col[2][1] - a32;
	amat[3][1] = -col[3][1] - a42;
	amat[4][1] = -col[4][1] - a52;

	/* level three balance equation */
	amat[0][2] = -col[0][2] - pump02;
	amat[1][2] = -col[1][2];
	amat[2][2] = a31 + a32 + col[2][0] + col[2][1] + col[2][3] + col[2][4];
	amat[3][2] = -col[3][2] - a43;
	amat[4][2] = -col[4][2] - a53;

	/* level four balance equation */
	amat[0][3] = -col[0][3] - pump03;
	amat[1][3] = -col[1][3];
	amat[2][3] = -col[2][3];
	amat[3][3] = a41 + col[3][0] + a42 + col[3][1] + a43 + col[3][2] + col[3][4];
	amat[4][3] = -col[4][3] - a54;

#	if 0
	// is it necessary to precondition the vars?
	/* divide both sides of equation by largest number to stop overflow */
	double dmax = -1e0;
	for( i=0; i < 6; i++ )
		for( j=0; j < 5; j++ )
			dmax = MAX2(zz[i][j],dmax);

	for( i=0; i < 6; i++ )
		for( j=0; j < 5; j++ )
			zz[i][j] /= dmax;
#	endif

	int32 ipiv[5], ner=0;
	/* solve matrix */
	getrf_wrapper(5,5,(double*)amat,5,ipiv,&ner);
	getrs_wrapper('N',5,1,(double*)amat,5,ipiv,bvec,5,&ner);

	if( ner != 0 )
	{
		fprintf( ioQQQ, "DISASTER PROBLEM atom_pop5: dgetrs finds singular or ill-conditioned matrix\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* p(5) was very slightly negative (1e-40) for SII in dqher.in - highest
	 * level has smallest excitation rates may be closest to ill conditioned*/
	p[1] = MAX2(0.e0,bvec[1]);
	p[2] = MAX2(0.e0,bvec[2]);
	p[3] = MAX2(0.e0,bvec[3]);
	p[4] = MAX2(0.e0,bvec[4]);
	// this should be majority and so numerically more
	p[0] = abund - p[1] - p[2] - p[3] - p[4];

	// level energies in ergs, needed for energy exchange
	double Erg[5] , EnergyKelvin[5];
	Erg[0] = 0.;
	EnergyKelvin[0] = 0.;
	for( long i=1; i<5; ++i )
	{
		Erg[i] = Erg[i-1] + EnerWN[i-1]*ERG1CM;
		EnergyKelvin[i] = EnergyKelvin[i-1] + EnerWN[i-1]*T1CM;
	}

	*Cooling = 0.;
	*CoolingDeriv = 0.;
	for( long ihi=1; ihi<5; ++ihi )
	{
		for( long ilo=0; ilo<ihi; ++ilo )
		{
			double CoolOne = (p[ilo]*col[ilo][ihi] - p[ihi]*col[ihi][ilo])* 
				(Erg[ihi]-Erg[ilo]);

			*Cooling += CoolOne;
			*CoolingDeriv += CoolOne*( EnergyKelvin[ihi]*thermal.tsq1 - thermal.halfte );
		}
	}
	
	return;
}
