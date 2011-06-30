/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolNick compute nickel cooling */
#include "cddefines.h"
#include "taulines.h"
#include "phycon.h"
#include "lines_service.h"
#include "dense.h"
#include "atoms.h"
#include "cooling.h"

void CoolNick(void)
{
	double cs,
		tused;
	realnum rate;

	DEBUG_ENTRY( "CoolNick()" );

	/*>>refer	Ni1	cs	Hollenbach, D. & McKee, C.F. 1989, ApJ, 342, 306 */
	/* rates are said to be ok over range 30 - 3000K */
	tused = MAX2( 30. , phycon.te );
	tused = MIN2( 3000. , phycon.te  );
	tused /= 100.;

	/* the 7.5 micron line */
	/* >>chng 03 nov 15, add these lines */
	rate = (realnum)(1.2e-7 * dense.eden + 
		/*8.0e-10*pow(tused, 0.17 )*dense.xIonDense[ipHYDROGEN][0]) / dense.eden);*/
		/* >>chng 05 jul 05, eden to cdsqte */
		8.0e-10*pow(tused, 0.17 )*dense.xIonDense[ipHYDROGEN][0] );
	LineConvRate2CS( &TauLines[ipNi1_7m]  , rate );

	/* the 11.3 micron line */
	rate = (realnum)(9.3e-8 * dense.eden + 
		/* >>chng 05 jul 05, eden to cdsqte */
		/*5.3e-10*pow(tused, 0.17 )*dense.xIonDense[ipHYDROGEN][0]) / dense.eden);*/
		5.3e-10*pow(tused, 0.17 )*dense.xIonDense[ipHYDROGEN][0] );
	LineConvRate2CS( &TauLines[ipNi1_11m]  , rate );

	rate = (realnum)(1.2e-7 * dense.eden + 
		/* >>chng 05 jul 05, eden to cdsqte */
		/*6.9e-10*pow(tused, 0.17 )*dense.xIonDense[ipHYDROGEN][0]) / dense.eden);*/
		6.9e-10*pow(tused, 0.17 )*dense.xIonDense[ipHYDROGEN][0] );
	TauDummy.Hi->g = TauLines[ipNi1_11m].Hi->g;
	LineConvRate2CS( &TauDummy  , rate );
	/* this says that line is a dummy, not real one */
	TauDummy.Hi->g = 0.;

	atom_level3(&TauLines[ipNi1_7m],&TauLines[ipNi1_11m],&TauDummy);

	/* [Ni XII] 4230.8A
	 * Y(ik) from 
	 * >>refer	ni12	cs	Pelan, J., & Berrington, K.A. 1995, A&A Suppl, 110, 209
	 * call PutCS( 1.68 , xNi1242 )
	 * >>chng 98 jan 09, update collision strength to 
	 * >>refer	ni12	cs	Mathews, A., Ramsbottom, C.A., Bell, K.L., & Keenan, F.P., 1998,
	 * >>refercon ApJ 492, 415 */
	if( phycon.alogte < 5.4 )
	{
		cs = MAX2(0.6,0.6+0.8182*(phycon.alogte-3.2));
	}
	else
	{
		cs = MAX2(0.6,2.4-1.5*(phycon.alogte-5.4));
	}
	PutCS(cs,&TauLines[ipxNi1242]);
	atom_level2(&TauLines[ipxNi1242]);
	return;
}
