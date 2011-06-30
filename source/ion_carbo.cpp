/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*IonCarbo compute ionization balance for carbon */
#include "cddefines.h"
#include "thermal.h"
#include "carb.h"
#include "trace.h"
#include "dense.h"
#include "phycon.h"
#include "hmi.h"
#include "mole.h"
#include "rfield.h"
#include "save.h"
#include "ionbal.h"

void IonCarbo(void)
{
	const int NDIM = ipCARBON+1;

	static const double dicoef[2][NDIM] = {
		{2.54e-3,6.15e-3,1.62e-3,4.78e-2,3.22e-2,0.}, {4.42e-2,5.88e-2,0.343,0.362,0.315,0.}
	};
	static const double dite[2][NDIM] = {
		{1.57e5,1.41e5,8.19e4,3.44e6,4.06e6,0.}, {3.74e5,1.41e5,1.59e5,5.87e5,8.31e5,0.}
	};
	static const double ditcrt[NDIM] = {1.2e4,1.2e4,1.1e4,4.4e5,7.0e5,1e20};
	/* >>refer	C	DR	Nussbaumer H. & Storey, P 1983, A&A, 126, 75 */
	static const double aa[NDIM] = {.0108,1.8267,2.3196,0.,0.,0.};
	static const double bb[NDIM] = {-0.1075,4.1012,10.7328,0.,0.,0.};
	static const double cc[NDIM] = {.2810,4.8443,6.8830,0.,0.,0.};
	static const double dd[NDIM] = {-0.0193,.2261,-0.1824,0.,0.,0.};
	/* NB C+ - C0 recom numerically unstable below 2000K, use est instead */
	static const double ff[NDIM] = {0.000,0.5960,0.4101,0.1,0.1,0.};

	DEBUG_ENTRY( "IonCarbo()" );

	if( trace.lgTrace && trace.lgCarBug )
	{
		fprintf( ioQQQ, "     IonCarbo called.\n" );
	}

	if( !dense.lgElmtOn[ipCARBON] )
	{
		carb.p1909 = 0.;
		carb.p2326 = 0.;
		thermal.heating[ipCARBON][9] = 0.;
		return;
	}

	/* zero out ionization balance arrays */
	ion_zero(ipCARBON);

	ion_photo(ipCARBON,false);

	/* >>chng 05 aug 10, add Leiden hack here to get same C0 photo rate
	 * as UMIST - negates difference in grain opacities */
	if(!co.lgUMISTrates)
	{
		int nelem=ipCARBON , ion=0 , ns=2;
		ionbal.PhotoRate_Shell[nelem][ion][ns][0] = 
			(HMRATE((1e-10)*3.0,0,0)*(hmi.UV_Cont_rel2_Habing_TH85_face*
			exp(-(3.0*rfield.extin_mag_V_point))/1.66));
		/* heating rates */
		ionbal.PhotoRate_Shell[nelem][ion][ns][1] = 0.;
		ionbal.PhotoRate_Shell[nelem][ion][ns][2] = 0.;
		/*fprintf(ioQQQ,"DEBUG hack %li\tav\t%.2e\tgo\t%.2e\trate\t%.2e\tabun\t%.2e\n",
			nzone ,
			rfield.extin_mag_V_point,
			hmi.UV_Cont_rel2_Habing_TH85_face,
			ionbal.PhotoRate_Shell[nelem][ion][ns][0] ,
			dense.xIonDense[ipCARBON][1]/SDIV(dense.xIonDense[ipCARBON][0]) );*/
	}

	/* find collisional ionization rates */
	ion_collis(ipCARBON);

	/* get recombination coefficients */
	/*lint -e64 double = pointer */
	ion_recomb(false,(const double*)dicoef,(const double*)dite,ditcrt,aa,bb,cc,dd,ff,ipCARBON);
	/*lint +e64 double = pointer */

	/* Photoproduction of 3P */
	carb.p1909 = ionbal.PhotoRate_Shell[ipCARBON][1][1][0];

	/* photoproduction of C II] 2326 */
	carb.p2326 = ionbal.PhotoRate_Shell[ipCARBON][0][1][0];


	bool lgDEBUG=false;
	/*if( iteration>1 )*/
	/*if( nzone > 258  )
		lgDEBUG = true;
	else
		lgDEBUG = false;*/
	ion_solver(ipCARBON,lgDEBUG);

	/* rate will be cm-3 s-1, into triplets only
	 * >>chng 96 nov 22, rid of carb() */
	carb.p1909 *= dense.xIonDense[ipCARBON][1]*0.62;
	carb.p2326 *= dense.xIonDense[ipCARBON][0]*0.1;

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "     IonCarbo returns;  fracs=" );
		for( int i=0; i < 7; i++ )
		{
			fprintf( ioQQQ, " %10.3e", dense.xIonDense[ipCARBON][i]/
			  dense.gas_phase[ipCARBON] );
		}
		fprintf( ioQQQ, "\n" );
	}

	enum {AGN=false};
	/* if set true there will be two problems at very low temps */
	if( AGN )
	{
		phycon.te=10.;
		/* this tells ion_recomb to save data */
		save.lgioRecom = true;
		/* this is where it will be punched */
		save.ioRecom = ioQQQ;
		while( phycon.te<1e7 )
		{
			/* get recombination coefficients */
			/*lint -e64 double = pointer */
			ion_recomb(false,(double*)dicoef,(double*)dite,ditcrt,aa,bb,cc,
				   dd,ff,ipCARBON);
			/*lint +e64 double = pointer */
			TempChange(phycon.te *2.f , true);
		}
		/* bail out */
		cdEXIT(EXIT_SUCCESS);
	}
	return;
}
