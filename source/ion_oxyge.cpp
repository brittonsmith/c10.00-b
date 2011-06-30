/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*IonOxyge derive ionization balance for oxygen */
#include "cddefines.h"
#include "opacity.h"
#include "oxy.h"
#include "thermal.h"
#include "dense.h"
#include "iso.h"
#include "trace.h"
#include "rfield.h"
#include "atmdat.h"
#include "atoms.h"
#include "gammas.h"
#include "ionbal.h"

void IonOxyge(void)
{
	const int NDIM = ipOXYGEN+1;

	static const double dicoef[2][NDIM] = {
		{1.11e-3,5.07e-3,1.48e-2,1.84e-2,4.13e-3,1.06e-1,6.23e-2,0.},
		{.0925,.181,.305,.1,.162,.34,.304,0.}
	};
	static const double dite[2][NDIM] = {
		{1.75e5,1.98e5,2.41e5,2.12e5,1.25e5,6.25e6,7.01e6,0.},
		{1.45e5,3.35e5,2.83e5,2.83e5,2.27e5,1.12e6,1.47e6,0.}
	};
	static const double ditcrt[NDIM] = {2.7e4,2.2e4,2.4e4,2.5e4,1.6e4,1.0e6,1.5e6,1e20};
	static const double aa[NDIM] = {0.,-0.0036,0.,0.0061,-2.8425,0.,0.,0.};
	static const double bb[NDIM] = {0.0238,0.7519,21.8790,0.2269,0.2283,0.,0.,0.};
	static const double cc[NDIM] = {0.0659,1.5252,16.2730,32.1419,40.4072,0.,0.,0.};
	static const double dd[NDIM] = {0.0349,-0.0838,-0.7020,1.9939,-3.4956,0.,0.,0.};
	static const double ff[NDIM] = {0.5334,0.2769,1.1899,-0.0646,1.7558,0.,0.,0.};

	bool lgDebug = false;
	long int iup;
	double aeff;

	DEBUG_ENTRY( "IonOxyge()" );

	/* oxygen, atomic number 8 */
	if( !dense.lgElmtOn[ipOXYGEN] )
	{
		oxy.poiii2Max = 0.;
		oxy.poiii3Max = 0.;
		oxy.r4363Max = 0.;
		oxy.r5007Max = 0.;
		oxy.poiii2 = 0.;
		oxy.p1666 = 0.;
		oxy.AugerO3 = 0.;
		oxy.p1401 = 0.;
		oxy.s3727 = 0.;
		oxy.s7325 = 0.;
		thermal.heating[7][9] = 0.;
		oxy.poimax = 0.;
		return;
	}

	ion_zero(ipOXYGEN);

	ion_photo(ipOXYGEN,false);

	/* find collisional ionization rates */
	ion_collis(ipOXYGEN);

	/* get recombination coefficients */
	/*lint -e64 type mismatch */
	ion_recomb(false,(const double*)dicoef,(const double*)dite,ditcrt,aa,bb,cc,dd,ff,ipOXYGEN);
	/*lint +e64 type mismatch */

	/* photoexcitation of O III 1666 and O IV 1401 */
	/** \todo	2	this will be zero in current form of atmdat_phfit
	 * set 2s**2 rate to rate for O V */
	oxy.p1666 = ionbal.PhotoRate_Shell[ipOXYGEN][3][1][0];

	oxy.p1401 = ionbal.PhotoRate_Shell[ipOXYGEN][2][1][0];

	t_phoHeat dummy;

	/* photoionization from O++ 1D
	 *
	 * estimate gamma function by assuming no frequency dependence
	 * betwen 1D and O++3P edge */
	/* destroy upper level of OIII 5007*/
	oxy.d5007r = (realnum)(GammaK(opac.ipo3exc[0],opac.ipo3exc[1],
	  opac.ipo3exc[2] , 1., &dummy ));

	/* destroy upper level of OIII 4363*/
	oxy.d4363 = (realnum)(GammaK(opac.ipo3exc3[0],opac.ipo3exc3[1],
	  opac.ipo3exc3[2] , 1., &dummy ));

	/* destroy upper level of OI 6300*/
	oxy.d6300 = (realnum)(GammaK(opac.ipo1exc[0],opac.ipo1exc[1],
	  opac.ipo1exc[2] , 1., &dummy ));

	/* A21 = 0.0263 */
	aeff = 0.0263 + oxy.d5007r;

	/* 1. as last arg makes this the relative population */
	oxy.poiii2 = (realnum)(atom_pop2(2.5,9.,5.,aeff,2.88e4,1.)/aeff);
	{
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC )
		{
			fprintf(ioQQQ,"pop rel  %.1e rate %.1e  grnd rate %.1e\n", 
				oxy.poiii2 , oxy.d5007r ,ionbal.PhotoRate_Shell[ipOXYGEN][2][2][0] );
		}
	}

	/* photoionization from excited states */
	if( nzone > 0 ) 
	{
		/* neutral oxygen destruction */
		ionbal.PhotoRate_Shell[ipOXYGEN][0][2][0] = ionbal.PhotoRate_Shell[ipOXYGEN][0][2][0]*
		  (1. - oxy.poiexc) + oxy.d6300*oxy.poiexc;

		/* doubly ionized oxygen destruction */
		ionbal.PhotoRate_Shell[ipOXYGEN][2][2][0] = ionbal.PhotoRate_Shell[ipOXYGEN][2][2][0]*
		  (1. - oxy.poiii2 - oxy.poiii3) + oxy.d5007r*oxy.poiii2 + 
		  oxy.d4363*oxy.poiii3;

		if( ionbal.PhotoRate_Shell[ipOXYGEN][2][2][0] > 1e-30 && dense.IonLow[ipOXYGEN] <= 2 )
		{
			if( (oxy.d5007r*oxy.poiii2 + oxy.d4363*oxy.poiii3)/
			  ionbal.PhotoRate_Shell[ipOXYGEN][2][2][0] > (oxy.r4363Max + 
			  oxy.r5007Max) )
			{
				oxy.poiii2Max = (realnum)(oxy.d5007r*oxy.poiii2/ionbal.PhotoRate_Shell[ipOXYGEN][2][2][0]);
				oxy.poiii3Max = (realnum)(oxy.d4363*oxy.poiii3/ionbal.PhotoRate_Shell[ipOXYGEN][2][2][0]);
			}
			oxy.r4363Max = (realnum)(MAX2(oxy.r4363Max,oxy.d4363));
			oxy.r5007Max = (realnum)(MAX2(oxy.r5007Max,oxy.d5007r));
		}

		/* ct into excited states */
		if( dense.IonLow[ipOXYGEN] <= 0 && (ionbal.PhotoRate_Shell[ipOXYGEN][0][2][0] + 
		  atmdat.HCharExcIonOf[ipOXYGEN][0]*dense.xIonDense[ipHYDROGEN][1]) > 1e-30 )
		{
			oxy.poimax = (realnum)(MAX2(oxy.poimax,oxy.d6300*oxy.poiexc/
			  (ionbal.PhotoRate_Shell[ipOXYGEN][0][2][0]+
			  atmdat.HCharExcIonOf[ipOXYGEN][0]* dense.xIonDense[ipHYDROGEN][1])));
		}
	}
	else
	{
		oxy.poiii2Max = 0.;
		oxy.poiii3Max = 0.;
		oxy.r4363Max = 0.;
		oxy.r5007Max = 0.;
		oxy.poimax = 0.;
	}

	/* save atomic oxygen photodistruction rate for 3727 creation */
	if( dense.IonLow[ipOXYGEN] == 0 && oxy.i2d < rfield.nflux )
	{
		oxy.s3727 = (realnum)(GammaK(oxy.i2d,oxy.i2p,opac.iopo2d , 1., &dummy ));

		iup = MIN2(iso.ipIsoLevNIonCon[ipH_LIKE][1][0],rfield.nflux);
		oxy.s7325 = (realnum)(GammaK(oxy.i2d,iup,opac.iopo2d , 1., &dummy ));

		oxy.s7325 -= oxy.s3727;
		oxy.s3727 = oxy.s3727 + oxy.s7325;

		/* ratio of cross sections */
		oxy.s7325 *= 0.66f;
	}
	else
	{
		oxy.s3727 = 0.;
		oxy.s7325 = 0.;
	}

	oxy.AugerO3 = (realnum)ionbal.PhotoRate_Shell[ipOXYGEN][0][0][0];

	/* solve for ionization balance */
	/**/if(0 &&  nzone > 100 )
		lgDebug = true;
	else
		lgDebug = false;
	ion_solver(ipOXYGEN,lgDebug);
	if( lgDebug )
		fprintf(ioQQQ,"DEBUG O\t%.3e\t%.3e\tH\t%.3e\t%.3e\n",
		dense.xIonDense[ipOXYGEN][0],
		dense.xIonDense[ipOXYGEN][1],
		dense.xIonDense[ipHYDROGEN][0],
		dense.xIonDense[ipHYDROGEN][1]);

	/* 1666 ratio corrected for phot crs at 50ev */
	oxy.p1666 *= dense.xIonDense[ipOXYGEN][1]*0.3;
	oxy.p1401 *= dense.xIonDense[ipOXYGEN][2]*0.43;
	oxy.s3727 *= dense.xIonDense[ipOXYGEN][0];
	oxy.s7325 *= dense.xIonDense[ipOXYGEN][0];
	oxy.AugerO3 *= dense.xIonDense[ipOXYGEN][0];

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "     IonOxyge returns; frac=" );
		for( int i=1; i <= 9; i++ )
		{
			fprintf( ioQQQ, " %10.3e", dense.xIonDense[ipOXYGEN][i-1]/
			  dense.gas_phase[ipOXYGEN] );
		}
		fprintf( ioQQQ, "\n" );
	}
	return;
}
