/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*lines main routine to put emission line intensities into line stack,
 * calls lineset1, 2, 3, 4 */
/*Drive_cdLine do the drive cdLine command */
/*FindStrongestLineLabels find strongest lines contributing to point in continuum energy mesh, output in some save commands */
#include "cddefines.h"
#include "physconst.h"
#include "taulines.h"
#include "thermal.h"
#include "yield.h"
#include "ipoint.h"
#include "ionbal.h"
#include "cddrive.h"
#include "iterations.h"
#include "trace.h"
#include "dense.h"
#include "prt.h"
#include "rt.h"
#include "coolheavy.h"
#include "rfield.h"
#include "phycon.h"
#include "elementnames.h"
#include "iso.h"
#include "hyperfine.h"
#include "hydrogenic.h"
#include "lines_service.h"
#include "atmdat.h"
#include "lines.h"
#include "radius.h"

STATIC void Drive_cdLine( void );

void lines(void)
{
	char chLabel[5];
	long int i, 
	  ipnt,
	  nelem;
	double BigstExtra, 
	  ExtraCool,  
	  f2, sum; 

	DEBUG_ENTRY( "lines()" );

	/* LineSave.ipass
	 * -1 - space not yet allocated - just count number of lines entered into stack
	 *  0 - first call with space allocated - must create labels and add in wavelengths
	 * +1 - later calls - just add intensity 
	 */

	/* major routines used here:
	 *
	 * PutLine( tarray )
	 * this uses information in tarray to give various
	 * contributions to lines, and their intensities
	 *
	 * PutExtra( extra )
	 * extra is some extra intensity to add to the line
	 * it will go into the totl contribution put out by PutLine,
	 * and this contribution should be indicated by independent
	 * call to linadd
	 * */

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " lines called\n" );
	}

	/* the drive cdline command, which checks that all lines can be pulled out by cdLine */
	if( trace.lgDrv_cdLine  && LineSave.ipass > 0 )
		Drive_cdLine();

	/* total luminosity radiated by this model - will be compared with energy in incident
	 * continuum when calculation is complete */
	thermal.power += thermal.htot*radius.dVeffAper;

	/* remember the total free-free heating */
	thermal.FreeFreeTotHeat += CoolHeavy.brems_heat_total*radius.dVeffAper;

	/* total Compton heating - cooling */
	rfield.comtot += rfield.cmheat*radius.dVeffAper;
	thermal.totcol += thermal.ctot*radius.dVeffAper;

	/* up up induced recombination cooling */
	for( nelem=0; nelem<LIMELM; ++nelem )
	{
		hydro.cintot += iso.RecomInducCool_Rate[ipH_LIKE][nelem]*radius.dVeffAper;
	}

	/* nsum is line pointer within large stack of line intensities */
	LineSave.nsum = 0;
	LineSave.nComment = 0;

	/* this is used by lindst to proportion inward and outward.  should be 50% for
	 * optically thin line.  putline sets this to actual value for particular line
	 * and calls lindst then rests to 50% */
	rt.fracin = 0.5;

	/* last arg in call to lindst and linadd is info on line
	 * info is char variable indicating type of line this is
	 * 'c' cooling
	 * 'h' heating
	 * 'i' information only
	 * 'r' recombination line
	 *
	 * all components of lines are entered into the main line stack here
	 * when printing, filters exist to not print Inwd component */

	/* initialize routine that can generate pointers for forbidden lines,
	 * these are lines that are not transferred otherwise,
	 * in following routines there will be pairs of calls, first to
	 * PntForLine to get pointer, then lindst to add to stack */
	PntForLine(0.,"FILL",&i);

	/* evaluate rec coefficient for rec lines of C, N, O
	 * some will be used in LineSet2 and then zeroed out,
	 * others left alone and used below */
	t_ADfA::Inst().rec_lines(phycon.te,LineSave.RecCoefCNO);

	/* initialize ExtraInten with a zero */
	PutExtra(0.);

	/* put in something impossible in element 0 of line stack */
	linadd(0.f,0,"zero",'i' , "null placeholder");

	/* this is compared with true volume in final.  The number can't
	 * actually be unity since this would overflow on a 32 bit machine */
	/* integrate the volume as a sanity check */
	linadd( 1.e-10 , 1 , "Unit" , 'i' , "unit integration placeholder");
	static long int ipOneAng=-1;
	if( LineSave.ipass<0 )
		ipOneAng = ipoint( RYDLAM );
	lindst( 1.e-10 , 1. , "UntD" , ipOneAng , 'i' , false,"unit integration placeholder");

	/* initial set of general properties */
	lines_general();

	/* do all continua */
	lines_continuum();

	/* information about grains */
	lines_grains();

	/* update all satellite lines */
	iso_satellite_update();

	/* do all hydrogenic ions */
	lines_hydro();

	/* enter He-iso sequence lines */
	lines_helium();

#if	0
	/* This is Ryan's code for dumping lots of Helium lines according to
	 * quantum number rather than wavelength, principally for comparison with Rob
	 * Bauman's code. */
	if( iteration > 1 )
	{
		fprintf( ioQQQ,"ipHi\tipLo\tnu\tlu\tsu\tnl\tll\tsl\tWL\tintens\n" );
		for( long ipHi=5; ipHi<= iso.numLevels_local[ipHE_LIKE][ipHELIUM] - iso.nCollapsed_local[ipHE_LIKE][ipHELIUM]; ipHi++ )
		{
			for( long ipLo=0; ipLo<ipHi; ipLo++ )
			{
				if( Transitions[ipHE_LIKE][ipHELIUM][ipHi][ipLo].ipCont > 0 )
				{
					double relint, absint;

					if( cdLine("He 1", 
						Transitions[ipHE_LIKE][ipHELIUM][ipHi][ipLo].WLAng,
						&relint, &absint ) )
					{
						//Transitions[ipHE_LIKE][ipHELIUM][ipHi][ipLo].Hi->chLabel

						if( Transitions[ipHE_LIKE][ipHELIUM][ipHi][ipLo].WLAng < 1.1E4 &&
							Transitions[ipHE_LIKE][ipHELIUM][ipHi][ipLo].WLAng > 3.59E3 &&
							ipLo!=3 && ipLo!=4 && relint >= 0.0009 )
						{
							fprintf( ioQQQ,"%li\t%li\t%li\t%li\t%li\t%li\t%li\t%li\t%e\t%e\n",
								ipHi,
								ipLo,
								StatesElemNEW[ipHELIUM][ipHELIUM-ipHE_LIKE][ipHi].n,
								StatesElemNEW[ipHELIUM][ipHELIUM-ipHE_LIKE][ipHi].l,
								StatesElemNEW[ipHELIUM][ipHELIUM-ipHE_LIKE][ipHi].S,
								StatesElemNEW[ipHELIUM][ipHELIUM-ipHE_LIKE][ipLo].n,
								StatesElemNEW[ipHELIUM][ipHELIUM-ipHE_LIKE][ipLo].l,
								StatesElemNEW[ipHELIUM][ipHELIUM-ipHE_LIKE][ipLo].S,
								Transitions[ipHE_LIKE][ipHELIUM][ipHi][ipLo].WLAng,
								relint );
						}
					}
				}
			}
		}
	}
#endif

	/* do heavies, lithium through neon */
	lines_lv1_li_ne();

	/* do heavies, sodium through argon */
	lines_lv1_na_ar();

	/* do heavies, potassium through zinc */
	lines_lv1_k_zn();

	/* add up line intensities for certain set of lines */
	sum = PrtLineSum();
	/* zero out the location that will receive this information, 
	 * remember that memory not allocated until ipass >= 0 */
	if( LineSave.ipass > 0 )
	{
		LineSv[LineSave.nsum].SumLine[0] = 0.;
		LineSv[LineSave.nsum].SumLine[1] = 0.;
	}
	/* optional sum of certain emission lines, set with "print sum" */
	linadd(sum/radius.dVeffAper,0,"Stoy",'i' , 
		"Stoy method energy sum ");

	/* next come some recombination lines */
	i = StuffComment( "recombination" );
	linadd( 0., (realnum)i , "####", 'i' ,
		"recombination lines");

	/***********************************************************************
	 * large number of C, N, and O recombination lines                     *
	 *************************************************************************/

	for( i=0; i < 471; i++ )
	{
		/* generate label for the line if ipass is -1 or 0 - saved in arrays
		 * so no need to do this during production */
		if( LineSave.ipass <= 0 )
		{
			/* generate label for the line */
			strcpy( chLabel, elementnames.chElementSym[(long)(LineSave.RecCoefCNO[0][i])-1] );
			strcat( chLabel, elementnames.chIonStage[(long)(LineSave.RecCoefCNO[0][i]-
				LineSave.RecCoefCNO[1][i]+1.01)-1] );
		}
		else
			chLabel[0] = ' ';

		/* number of rec per unit vol
		 * do not predict these simple reccombination intensities at high densities
		 * since lines neglect collisional deexciation and line optical depth effects.
		 * They were not intended for high densities or column densities.
		 * As a result they become unphysically bright at high densities and
		 * violate the black body limit.  There would be major
		 * energy conservation problems if they were added in the outward beam in
		 * dense simulations.
		 * */
		if( dense.eden < 1e8  )
		{
			nelem = (long)LineSave.RecCoefCNO[0][i]-1;
			long int ion = (long)(LineSave.RecCoefCNO[0][i]-LineSave.RecCoefCNO[1][i]+2)-1;
			f2 = LineSave.RecCoefCNO[3][i]*dense.eden*
				dense.xIonDense[nelem][ion];

			/* convert to intensity */
			f2 = f2*1.99e-8/LineSave.RecCoefCNO[2][i];
		}
		else
		{
			f2 = 0.;
		}
		/* stuff it into the stack */
		PntForLine(LineSave.RecCoefCNO[2][i],chLabel,&ipnt);
		lindst(f2,LineSave.RecCoefCNO[2][i],chLabel,ipnt, 'r',true ,
			"recombination line");
	}

	/* next come the atom_level2 lines */
	i = StuffComment( "level2 lines" );
	linadd( 0., (realnum)i , "####", 'i' ,
		"level2 lines");

	/* add in all the other level 2 wind lines
	 * Dima's 6k lines */
	ExtraCool = 0.;
	BigstExtra = 0.;
	for( i=0; i < nWindLine; i++ )
	{
		if( TauLine2[i].Hi->IonStg < TauLine2[i].Hi->nelem+1-NISO )
		{
			PutLine(&TauLine2[i],"level 2 line");
			if( TauLine2[i].Coll.cool > BigstExtra )
			{
				BigstExtra = TauLine2[i].Coll.cool;
				thermal.ipMaxExtra = i+1;
			}
			ExtraCool += TauLine2[i].Coll.cool;
		}
	}
	/* keep track of how important this is */
	thermal.GBarMax = MAX2(thermal.GBarMax,(realnum)(ExtraCool/thermal.ctot));

	/* next come the hyperfine structure lines */
	i = StuffComment( "hyperfine structure" );
	linadd( 0., (realnum)i , "####", 'i' ,
		"hyperfine structure lines ");

	/* this is total cooling due to all HF lines */
	linadd( hyperfine.cooling_total, 0., "hfin", 'i' ,
		"total cooling all hyperfine structure lines ");

	/* remember largest local cooling for possible printout in comments */
	hyperfine.cooling_max = (realnum)MAX2(hyperfine.cooling_max,hyperfine.cooling_total/thermal.ctot);

	/* the hyperfine lines */
	for( i=0; i < nHFLines; i++ )
	{
		PutLine(&HFLines[i],
			"hyperfine structure line");
	}

	/* external database lines */
	i = StuffComment( "database lines" );
	linadd( 0., (realnum)i , "####", 'i' ,
		"database lines ");

	/* external database lines */
	for( i=0; i < linesAdded2; i++ )
	{
		/* \todo 2 say which database in the comment */
		PutLine(dBaseLines[i].tran, "lines from third-party databases", dBaseLines[i].tran->Hi->chLabel);
	}

	/* next come the inner shell fluorescent lines */
	i = StuffComment( "inner shell" );
	linadd( 0., (realnum)i , "####", 'i' ,
		"inner shell lines");

	/* the group of inner shell fluorescent lines */
	for( i=0; i < t_yield::Inst().nlines(); ++i )
	{
		double xInten = 
			/* density of parent ion, cm-3 */
			dense.xIonDense[t_yield::Inst().nelem(i)][t_yield::Inst().ion(i)] *
			/* photo rate per atom per second, s-1 */
			ionbal.PhotoRate_Shell[t_yield::Inst().nelem(i)][t_yield::Inst().ion(i)][t_yield::Inst().nshell(i)][0]*
			/* fluor yield - dimensionless */
			t_yield::Inst().yield(i) *
			/* photon energy - ryd, converted into ergs */
			t_yield::Inst().energy(i) * EN1RYD;

		/* create label if initializing line stack */
		if( LineSave.ipass == 0 )
		{
			/* only generate the line label if it is going to be used */
			strcpy( chLabel , elementnames.chElementSym[t_yield::Inst().nelem(i)] );
			strcat( chLabel , elementnames.chIonStage[t_yield::Inst().ion_emit(i)] );
#			if 0
			/* only print yields for atoms */
			if( t_yield::Inst().ion(i) == 0 && t_yield::Inst().nelem(i) == ipIRON )
			fprintf(ioQQQ,"DEBUGyeild\t%s\t%.3f\t%.3e\n",
				/* line designation, energy in eV, yield */
				chLabel , t_yield::Inst().energy(i)*EVRYD, t_yield::Inst().yield(i) );
#			endif
		}

		/* the group of inner shell fluorescent lines */
		lindst(
			/* intensity of line */
			xInten,
			/* wavelength of line in Angstroms */
			(realnum)RYDLAM / t_yield::Inst().energy(i),
			/* label */
			chLabel ,
			/* continuum array offset for line as set in ipoint */
			t_yield::Inst().ipoint(i), 
			/* type of line - count as a recombination line */
			'r',
			/* include line in continuum? */
			true ,
			"inner shell line");
	}

	/* >>chng 06 jan 03, confirm that number of lines never changes once we have
	 * created the labels */
	{
		static long nLineSave=-1 , ndLineSave=-1;
		if( LineSave.ipass == 0 )
		{
			nLineSave = LineSave.nsum;
			ndLineSave = LineSave.nsum;
		}
		else if( LineSave.ipass > 0 )
		{
			/* this can't happen */
			if( nLineSave<= 0 || ndLineSave < 0 )
				TotalInsanity();

			/* now make sure that we have the same number of lines as we had previously
			 * created labels.  This would not pass if we did not add in exactly the same
			 * number of lines on each pass */
			if( nLineSave != LineSave.nsum )
			{
				fprintf( ioQQQ, "DISASTER number of lines in LineSave.nsum changed between pass 0 and 1 - this is impossible\n" );
				fprintf( ioQQQ, "DISASTER LineSave.nsum is %li and nLineSave is %li\n",
					LineSave.nsum , 
					nLineSave);
				ShowMe();
				cdEXIT(EXIT_FAILURE);
			}
			if( ndLineSave != LineSave.nsum )
			{
				fprintf( ioQQQ, "DISASTER number of lines in LineSave.nsum changed between pass 0 and 1 - this is impossible\n" );
				fprintf( ioQQQ, "DISASTER LineSave.nsum is %li and ndLineSave is %li\n",
					LineSave.nsum , 
					ndLineSave);
				ShowMe();
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	/* now do all molecules - do last since so many H2 lines */
	lines_molecules();

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " lines returns\n" );
	}
	return;
}

/*Drive_cdLine do the drive cdLine command */
STATIC void Drive_cdLine( void )
{
	long int j;
	bool lgMustPrintHeader = true;
	double absval , rel;

	DEBUG_ENTRY( "Drive_cdLine()" );

	for( j=1; j < LineSave.nsum; j++ )
	{
		if( cdLine( LineSv[j].chALab , LineSv[j].wavelength , &absval , &rel ) <= 0 )
		{
			/* print header if first miss */
			if( lgMustPrintHeader )
			{
				fprintf(ioQQQ,"n\tlab\twl\n");
				lgMustPrintHeader = false;
			}

			fprintf(ioQQQ,"%li\t%s\t%f\n", j, LineSv[j].chALab , LineSv[j].wavelength );
		}
	}
	fprintf( ioQQQ, " Thanks for checking on the cdLine routine!\n" );
	cdEXIT(EXIT_FAILURE);
}

