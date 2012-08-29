/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolSave save coolants */
#include "cddefines.h"
#include "thermal.h"
#include "dynamics.h"
#include "radius.h"
#include "conv.h"
#include "phycon.h"
#include "save.h"

/* this is limit to number of coolants to print out */
static const int IPRINT = 100;

/*CoolSave save coolants */
void CoolSave(FILE * io)
{
	long int i, 
	  ip, 
	  is;

	int nFail;

	double cset, 
		cool_total,
		heat_total;

	/*
	  8-29-2012 - Britton changed precision to double.
	realnum
	*/
	double
		*csav,
		*sgnsav;
	long int *index;

	DEBUG_ENTRY( "CoolSave()" );

	/* cannot do one-time init since thermal.ncltot can change */
	index = (long int *)CALLOC((size_t)thermal.ncltot,sizeof(long int));
	/*
	  8-29-2012 - Britton changed precision to double.
	csav = (realnum *)CALLOC((size_t)thermal.ncltot,sizeof(realnum));
	sgnsav = (realnum *)CALLOC((size_t)thermal.ncltot,sizeof(realnum));
	*/
	csav = (double *)CALLOC((size_t)thermal.ncltot,sizeof(double));
	sgnsav = (double *)CALLOC((size_t)thermal.ncltot,sizeof(double));

	cool_total = thermal.ctot;
	heat_total = thermal.htot;

	/* >>chng 06 mar 17, comment out following block and replace with this 
	 * removing dynamics heating & cooling and report only physical
	 * heating and cooling 
	 * NB the heating and cooling as punched no longer need be
	 * equal for a converged model */
	cool_total -= dynamics.Cool();
	heat_total -= dynamics.Heat();
#	if 0
	if(dynamics.Cool > dynamics.Heat()) 
	{
		cool_total -= dynamics.Heat();
		heat_total -= dynamics.Heat();
	} 
	else
	{
		cool_total -= dynamics.Cool;
		heat_total -= dynamics.Cool;
	}
#	endif

	/* cset will be weakest cooling to consider
	 * WeakHeatCool set with 'set weakheatcool' command
	 * default is 0.05 */
	cset = cool_total*save.WeakHeatCool;

	/* first find all strong lines, both + and - sign */
	ip = thermal.ncltot;

	for( i=0; i < ip; i++ )
	{
	  /* 
	     8-29-2012 - Britton changed this so cooling rates are printed 
	     instead of fractions of the total cooling rate.
		csav[i] = (realnum)(MAX2(thermal.cooling[i],thermal.heatnt[i])/
			SDIV(cool_total));
	  */
	        csav[i] = (double)(MAX2(thermal.cooling[i],thermal.heatnt[i]));

		/* save sign to remember if heating or cooling line */
		if( thermal.heatnt[i] == 0. )
		{
			sgnsav[i] = 1.;
		}
		else
		{
			sgnsav[i] = -1.;
		}
	}

	/* order strongest to weakest */
	/* now sort by decreasing importance */
	/*spsort netlib routine to sort array returning sorted indices */
	spsort_double(
		  /* input array to be sorted */
		  csav, 
		  /* number of values in x */
		  ip, 
		  /* permutation output array */
		  index, 
		  /* flag saying what to do - 1 sorts into increasing order, not changing
		   * the original routine */
		  -1, 
		  /* error condition, should be 0 */
		  &nFail);

	/* warn if tcovergence failure occurred */
	if( !conv.lgConvTemp )
	{
		fprintf( io, "#>>>>  Temperature not converged.\n" );
	}
	else if( !conv.lgConvEden )
	{
		fprintf( io, "#>>>>  Electron density not converged.\n" );
	}
	else if( !conv.lgConvIoniz )
	{
		fprintf( io, "#>>>>  Ionization not converged.\n" );
	}
	else if( !conv.lgConvPres )
	{
		fprintf( io, "#>>>>  Pressure not converged.\n" );
	}

	/*>>chng 06 jun 06, change start of save to give same info as heating 
	 * as per comment by Yumihiko Tsuzuki */
	/* begin the print out with zone number, total heating and cooling */
	/* 
	   8-29-2012 - Britton changed the heating and cooling output 
	   precision to six digits.
	fprintf( io, "%.5e\t%.4e\t%.4e\t%.4e", 
	*/
	fprintf( io, "%.5e\t%.4e\t%.6e\t%.6e", 
		radius.depth_mid_zone, 
		phycon.te, 
		heat_total, 
		cool_total );

	/* print only up to IPRINT, which is defined above */
	ip = MIN2( ip , IPRINT );

	/* now print the coolants 
	 * keep sign of coolant, for strong negative cooling 
	 * order is ion, wavelength, fraction of total */
	for( is=0; is < ip; is++ )
	{
		if(is > 4 && (thermal.cooling[index[is]] < cset && thermal.heatnt[index[is]] < cset))
			break;
		/*
		  8-29-2012 - Britton changed this so cooling rates are output
		  instead of fractions of the total cooling rate.
		fprintf( io, "\t%s %.1f\t%.7f", 
		*/
		fprintf( io, "\t%s %.1f\t%.6e", 
			thermal.chClntLab[index[is]], 
			thermal.collam[index[is]], 
			sign(csav[index[is]],sgnsav[index[is]]) );
	}
	fprintf( io, " \n" );

	/* finished, now free space */
	free(sgnsav);
	free(csav);
	free(index);
	return;
}
