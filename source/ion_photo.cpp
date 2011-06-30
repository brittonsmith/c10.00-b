/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ion_photo fill array PhotoRate with photoionization rates for heavy elements */
#include "cddefines.h"
#include "yield.h"
#include "heavy.h"
#include "opacity.h"
#include "dense.h"
#include "thermal.h"
#include "conv.h"
#include "grainvar.h"
#include "elementnames.h"
#include "gammas.h"
#include "ionbal.h"

void ion_photo(
	/* nlem is atomic number on C scale, 0 for H */
	long int nelem , 
	/* debugging flag to turn on print */
	bool lgPrintIt )
{
	long int ion, 
	  iphi, 
	  iplow, 
	  ipop, 
	  limit_hi, 
	  limit_lo,
	  ns;

	DEBUG_ENTRY( "ion_photo()" );

	/* IonLow(nelem) and IonHigh(nelem) are bounds for evaluation*/

	/*begin sanity checks */
	ASSERT( nelem < LIMELM );
	ASSERT( dense.IonLow[nelem] >= 0 );
	ASSERT( dense.IonLow[nelem] <= nelem);
	ASSERT( dense.IonHigh[nelem] <= nelem + 1);
	/*end sanity checks */

	/* NB - in following, nelem is on c scale, so is 29 for Zn */

	/* min since iso-like atom rates produced in iso_photo.
	 * IonHigh and IonLow range from 0 to nelem+1, but for photo rates
	 * we want 0 to nelem since cannot destroy fully stripped ion.  
	 * since iso-seq done elsewhere, want to actually do IonHigh-xxx.*/
	/* >>chng 00 dec 07, logic on limit_hi now precisely identical to ion_solver */
	/* >>chng 02 mar 31, change limit_hi to < in loop (had been <=) and
	 * also coded for arbitrary number of iso sequences */
	/*limit_hi = MIN2(nelem,dense.IonHigh[nelem]);
	limit_hi = MIN2(nelem-2,dense.IonHigh[nelem]-1);*/
	limit_hi = MIN2( dense.IonHigh[nelem] , nelem+1-NISO );

	/* >>chng 03 sep 26, always do atom itself since may be needed for molecules */
	limit_hi = MAX2( 1 , limit_hi );

	/* when grains are present want to use atoms as lower bound to number of stages of ionization,
	 * since atomic rates needed for species within grains */
	if( !conv.nPres2Ioniz && gv.lgDustOn() )
	{
		limit_lo = 0;
	}
	else
	{
		limit_lo = dense.IonLow[nelem];
	}

	/* >>chng 01 dec 11, lower bound now limit_lo */
	/* loop over all ions for this element */
	/* >>chng 02 mar 31, now ion < limit_hi not <= */
	for( ion=limit_lo; ion < limit_hi; ion++ )
	{
		/* loop over all shells for this ion */
		for( ns=0; ns < Heavy.nsShells[nelem][ion]; ns++ )
		{
			/* always reevaluate the outer shell, and all shells if lgRedoStatic is set */
			if( (ns==(Heavy.nsShells[nelem][ion]-1) || opac.lgRedoStatic) )
			{
				/* option to redo the rates only on occasion */
				iplow = opac.ipElement[nelem][ion][ns][0];
				iphi = opac.ipElement[nelem][ion][ns][1];
				ipop = opac.ipElement[nelem][ion][ns][2];

				t_phoHeat photoHeat;

				/* compute the photoionization rate, ionbal.lgPhotoIoniz_On is 1, set 0
				 * with "no photoionization" command */
				ionbal.PhotoRate_Shell[nelem][ion][ns][0] =
					GammaK(iplow,iphi,
					ipop,t_yield::Inst().elec_eject_frac(nelem,ion,ns,0),
					&photoHeat )*ionbal.lgPhotoIoniz_On;

				/* these three lines must be kept parallel with the lines
				 * in GammaK ion*/

				/* the heating rate */
				ionbal.PhotoRate_Shell[nelem][ion][ns][1] = photoHeat.HeatLowEnr*ionbal.lgPhotoIoniz_On;
				ionbal.PhotoRate_Shell[nelem][ion][ns][2] = photoHeat.HeatHiEnr*ionbal.lgPhotoIoniz_On;
			}
		}

		/* add on compton recoil ionization for atoms to outer shell */
		/* >>chng 02 mar 24, moved here from ion_solver */
		/* this is the outer shell */
		ns = (Heavy.nsShells[nelem][ion]-1);
		/* this must be moved to photoionize and have code parallel to iso_photo code */
		ionbal.PhotoRate_Shell[nelem][ion][ns][0] += ionbal.CompRecoilIonRate[nelem][ion];
		/* add the heat as secondary-ionization capable heating */
		ionbal.PhotoRate_Shell[nelem][ion][ns][2] += ionbal.CompRecoilHeatRate[nelem][ion];
	}

	/* option to print information about these rates for this element */
	if( lgPrintIt )
	{
		/* option to print rates for particular shell */
		ns = 5;
		ion = 1;
		GammaPrt(
		  opac.ipElement[nelem][ion][ns][0], 
		  opac.ipElement[nelem][ion][ns][1], 
		  opac.ipElement[nelem][ion][ns][2], 
		  ioQQQ, /* io unit we will write to */
		  ionbal.PhotoRate_Shell[nelem][ion][ns][0], 
		  0.05);

		/* outer loop is from K to most number of shells present in atom */
		for( ns=0; ns < Heavy.nsShells[nelem][0]; ns++ )
		{
			fprintf( ioQQQ, "\n %s", elementnames.chElementNameShort[nelem] );
			fprintf( ioQQQ, " %s" , Heavy.chShell[ns]);
			/* MB hydrogenic photo rate may not be included in beow */
			for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
			{
				if( Heavy.nsShells[nelem][ion] > ns )
				{
					fprintf( ioQQQ, " %8.1e", ionbal.PhotoRate_Shell[nelem][ion][ns][0] );
				}
				else
				{
					break;
				}
			}
		}
		fprintf(ioQQQ,"\n");
	}
	return;
}
