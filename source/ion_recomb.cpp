/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ion_recomb generate recombination coefficients for any species */
/*ion_recombAGN generate recombination coefficients for AGN table */
#include "cddefines.h"
#include "phycon.h"
#include "heavy.h"
#include "hmi.h"
#include "grainvar.h"
#include "dense.h"
#include "conv.h"
#include "thermal.h"
#include "iso.h"
#include "abund.h"
#include "save.h"
#include "elementnames.h"
#include "atmdat.h"
#include "ionbal.h"

/*ion_recomb generate recombination coefficients for any species */
void ion_recomb(
  /* this is debug flag */
  bool lgPrintIt,
  const double *dicoef, 
  const double *dite, 
  const double ditcrt[], 
  const double aa[], 
  const double bb[], 
  const double cc[], 
  const double dd[], 
  const double ff[], 
  /* nelem is the atomic number on the C scale, 0 for H */
  long int nelem/*, 
  double tlow[]*/)
{
#define DICOEF(I_,J_)	(*(dicoef+(I_)*(nelem+1)+(J_)))
#define DITE(I_,J_)	(*(dite+(I_)*(nelem+1)+(J_)))
	long int i, 
	  ion, 
	  limit;
	double 
	  fac2, 
	  factor,
	  DielRecomRateCoef_HiT[LIMELM],
	  DielRecomRateCoef_LowT[LIMELM],
	  ChargeTransfer[LIMELM];

	/* these are used for adding noise to rec coef */
	static double RecNoise[LIMELM][LIMELM];
	static bool lgNoiseNeedEval=true;

	DEBUG_ENTRY( "ion_recomb(false,)" );

	/* set up ionization balance matrix, C(I,1)=destruction, 2=creation
	 * heating rates saved in array B(I) in same scratch block
	 * factor is for Aldrovandi+Pequignot fit, FAC2 is for Nuss+Storey
	 * fit for dielectronic recombination
	 * GrnIonRec is rate ions recombine on grain surface, normally zero;
	 * set in hmole, already has factor of hydrogen density
	 * */

	/* routine only used for Li on up */
	ASSERT( nelem < LIMELM);
	ASSERT( nelem > 1 );

	/* check that range of ionization is correct */
	ASSERT( dense.IonLow[nelem] >= 0 );
	ASSERT( dense.IonLow[nelem] <= nelem+1 );

	atmdat.nsbig = MAX2(dense.IonHigh[nelem]+1,atmdat.nsbig);
	fac2 = 1e-14*phycon.sqrte;

	/* option to put noise into rec coefficient -
	 * one time initialization of noise  - this is set with command
	 * set dielectronic recombination kludge noise */
	if( lgNoiseNeedEval )
	{
		int n;
		if( ionbal.guess_noise !=0. )
		{
			for( n=ipHYDROGEN; n<LIMELM; ++n )
			{
				for( ion=0; ion<=n; ++ion )
				{
					/* log normal noise with dispersion entered on command line */
					/* NB the seed for rand was set when the command was parsed */
					RecNoise[n][ion] = pow(10., RandGauss( 0. , ionbal.guess_noise ) );
				}
			}
		}
		else
		{
			for( n=ipHYDROGEN; n<LIMELM; ++n )
			{
				for( ion=0; ion<=n; ++ion )
				{
					RecNoise[n][ion] = 1.;
				}
			}
		}
		lgNoiseNeedEval = false;
	}

	/* this routine only does simple two-level species, 
	 * loop over ions will be <= limit, IonHigh is -1 since very
	 * highest stage of ionization is not recombined into.  
	 * for Li, will do only atom, since ions are H and He like,
	 * so limit is zero */
	limit = MIN2(nelem-NISO,dense.IonHigh[nelem]-1);
	ASSERT( limit >= -1 );

	/* zero-out loop comes before main loop since there are off-diagonal
	 * elements in the main ionization loop, due to multi-electron processes */
	/* >>chng 00 dec 07, limit changed to identical to ion_solver */
	for( ion=0; ion <= limit; ion++ )
	{
		ionbal.RateRecomTot[nelem][ion] = 0.;
		ChargeTransfer[ion] = 0.;
		DielRecomRateCoef_LowT[ion] = 0.;
		DielRecomRateCoef_HiT[ion] = 0.;
		ionbal.RR_rate_coef_used[nelem][ion] = 0.;
		ionbal.DR_rate_coef_used[nelem][ion] = 0.;
	}
	for( ion=limit+1; ion < LIMELM; ion++ )
	{
		/* >>chng 01 dec 18, do not set this to FLT_MAX since it clobbers what
		 * had been set in h-like and he-like routines - that would only affect
		 * the printout */
		ChargeTransfer[ion] = -FLT_MAX;
		DielRecomRateCoef_LowT[ion] = -FLT_MAX;
		DielRecomRateCoef_HiT[ion] = -FLT_MAX;
	}

	DielRecomRateCoef_HiT[nelem] = 0.;
	DielRecomRateCoef_HiT[nelem-1] = 0.;

	DielRecomRateCoef_LowT[nelem] = 0.;
	DielRecomRateCoef_LowT[nelem-1] = 0.;

	/* these are counted elsewhere and must not be added here */
	Heavy.xLyaHeavy[nelem][nelem] = 0.;
	Heavy.xLyaHeavy[nelem][nelem-1] = 0.;

	/* IonLow is 0 for the atom, limit chosen to NOT include iso sequences  */
	for( ion=dense.IonLow[nelem]; ion <= limit; ion++ )
	{
		/* number of bound electrons of the ion after recombination, 
		 * for an atom (ion=0) this is
		 * equal to nelem+1, the element on the physical scale, since nelem is 
		 * on the C scale, being 5 for carbon */
		long n_bnd_elec_after_recomb = nelem+1 - ion;

#if	0
		/* >>chng 02 nov 06 add charge transfer here rather than in ion_solver */
		/* charge transfer recombination of this species by ionizing hydrogen and helium */
		ChargeTransfer[ion] = 
			/* He0 + ion charge transfer recombination */
			atmdat.HeCharExcRecTo[nelem][ion]*
			/* following is density [cm-3] of ground state of He0 */
			StatesElemNEW[ipHELIUM][ipHELIUM-ipHE_LIKE][ipHe1s1S].Pop + 
			/* H0 + ion charge transfer recombination */
			atmdat.HCharExcRecTo[nelem][ion]*
			/* following is density [cm-3] of ground state of H0 */
			StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop;
#else
		ChargeTransfer[ion] = 
			/* He0 + ion charge transfer recombination */
			atmdat.HeCharExcRecTo[nelem][ion]*
			/* following is density [cm-3] of He0 */
			dense.xIonDense[ipHELIUM][0] + 
			/* H0 + ion charge transfer recombination */
			atmdat.HCharExcRecTo[nelem][ion]*
			/* following is density [cm-3] of H0 */
			dense.xIonDense[ipHYDROGEN][0];
#endif
		/*>>chng 04 feb 20, add this, had always been in for destruction for H- */
		/* charge transfer recombination of first ion to neutral, by reaction with H- 
		 * the ion==0 is right, the first array element is the */
		if( ion==0 && nelem>ipHELIUM && atmdat.lgCTOn )
			ChargeTransfer[ion] += hmi.hmin_ct_firstions * hmi.Hmolec[ipMHm];

		/* >>chng 06 feb 01, add option to use Badnell RR data rather than Verner */
		if( ionbal.lgRR_recom_Badnell_use && ionbal.lgRR_Badnell_rate_coef_exist[nelem][ion] )
		{
			ionbal.RR_rate_coef_used[nelem][ion] = ionbal.RR_Badnell_rate_coef[nelem][ion];
		}
		else
		{
			ionbal.RR_rate_coef_used[nelem][ion] = ionbal.RR_Verner_rate_coef[nelem][ion];
		}

		/* Burgess or high-T dielectronic recombination */
		DielRecomRateCoef_HiT[ion] = 0.;
		/* >>chng 03 oct 29, do fe here rather than after loop */
		if( nelem==ipIRON )
		{
			/* implements dn > 0 DR from Arnaud & Raymond 1992 */
			DielRecomRateCoef_HiT[ion] = atmdat_dielrec_fe(ion+1,phycon.te);
		}
		else if( phycon.te > (ditcrt[ion]*0.1) )
		{
			DielRecomRateCoef_HiT[ion] = ionbal.DielSupprs[0][ion]/phycon.te32*
			  DICOEF(0,ion)*exp(-DITE(0,ion)/phycon.te)*
			  (1. + DICOEF(1,ion)*
			  sexp(DITE(1,ion)/phycon.te));
		}

		/* begin dn = 0 dielectronic recombination
		 * do not include it for rec from
		 * a closed shell, n_bnd_elec_after_recomb-1 is number of bound electrons in parent ion */
		DielRecomRateCoef_LowT[ion] = 0.;
		if( ((n_bnd_elec_after_recomb-1) !=  2) && 
		    ((n_bnd_elec_after_recomb-1) != 10) && 
		    ((n_bnd_elec_after_recomb-1) != 18) )
		{
			/* do not do iron here since all dn = 0 DR either Badnell or a guess */
			if( ff[ion] != 0. && nelem != ipIRON )
			{
				double t4m1 = 1e4/phycon.te;
				double tefac = ff[ion]*t4m1;
				/* >>chng 06 feb 14, O+3 ff[ion] is very negative, as a result exp goes to +inf
				 * at very low temperatures.  This is a error in the Nussbaumer & Storey fits to DR.
				 * do not use them is tefac = ff[ion] / t4 is very negative */
				if( tefac > -5. )
				{
					factor = (((aa[ion]*t4m1+bb[ion])*t4m1+cc[ion])*t4m1+dd[ion])* sexp(tefac);
					DielRecomRateCoef_LowT[ion] = ionbal.DielSupprs[1][ion]*fac2*
						MAX2(0.,factor );
				}
				else
				{
					DielRecomRateCoef_LowT[ion] = 0.;
				}
			}
			else if( ionbal.lg_guess_coef )
			{
				/* use mean of Badnell DR rates */
				DielRecomRateCoef_LowT[ion] = ionbal.DR_Badnell_rate_coef_mean_ion[ion];
				/* include optional noise here 
				 * >>chng 06 feb 07, move noise down to here so that use for both
				 * guesses of DR rates */
				DielRecomRateCoef_LowT[ion] *= RecNoise[nelem][ion];
			}
		}
		/* >>chng 05 dec 19, add option to use Badnell numbers */
		/* this is total old DR rates - may not use it */
		ionbal.DR_old_rate_coef[nelem][ion] = DielRecomRateCoef_HiT[ion] + DielRecomRateCoef_LowT[ion];

		/* set total DR rate - either Badnell if it exists or on with set dielectronic recombination badnell */
		if( ionbal.lgDR_recom_Badnell_use && ionbal.lgDR_Badnell_rate_coef_exist[nelem][ion] )
		{
			ionbal.DR_rate_coef_used[nelem][ion] = ionbal.DR_Badnell_rate_coef[nelem][ion];
		}
		else
		{
			ionbal.DR_rate_coef_used[nelem][ion] = ionbal.DR_old_rate_coef[nelem][ion];
		}

		/* sum of recombination rates [units s-1] for radiative, three body, charge transfer */
		ionbal.RateRecomTot[nelem][ion] = 
			dense.eden* (
			ionbal.RR_rate_coef_used[nelem][ion] + 
			ionbal.DR_rate_coef_used[nelem][ion] +
			ionbal.CotaRate[ion] ) + 
			ChargeTransfer[ion];

		/* >>chng 01 jun 30, FRAC_LINE was 0.1, not 1, did not include anything except
		 * radiative recombination, the radrec term */
#		define FRAC_LINE 1.
		/* was 0.1 */
		/*Heavy.xLyaHeavy[nelem][ion] = (realnum)(dense.eden*radrec*FRAC_LINE );*/
		Heavy.xLyaHeavy[nelem][ion] = (realnum)(dense.eden*
			(ionbal.RR_rate_coef_used[nelem][ion]+ionbal.DR_rate_coef_used[nelem][ion])*FRAC_LINE );
	}

	/* option to save rec coefficients */
	if( save.lgioRecom || lgPrintIt )
	{
		/* >>chng 04 feb 22, make option to print ions for single element */
		FILE *ioOut;
		if( lgPrintIt )
			ioOut = ioQQQ;
		else
			ioOut = save.ioRecom;

		/* print name of element */
		fprintf( ioOut, 
			" %s recombination coefficients fnzone:%.2f \tte\t%.4e\tne\t%.4e\n", 
			elementnames.chElementName[nelem] , fnzone , phycon.te , dense.eden );

		/*limit = MIN2(11,dense.IonHigh[nelem]);*/
		/* >>chng 05 sep 24, just print one long line - need info */
		limit = dense.IonHigh[nelem];
		// give ion stage
		for( i=0; i<limit; ++i )
			fprintf( ioOut, "%10ld",i+1);
		fprintf( ioOut, "\n");
		
		for( i=0; i < limit; i++ )
		{
			fprintf( ioOut, "%10.2e", ionbal.RR_rate_coef_used[nelem][i] );
		}
		fprintf( ioOut, " radiative used vs Z\n" );

		for( i=0; i < limit; i++ )
		{
			fprintf( ioOut, "%10.2e", ionbal.RR_Verner_rate_coef[nelem][i] );
		}
		fprintf( ioOut, " old Verner vs Z\n" );

		for( i=0; i < limit; i++ )
		{
			fprintf( ioOut, "%10.2e", ionbal.RR_Badnell_rate_coef[nelem][i] );
		}
		fprintf( ioOut, " new Badnell vs Z\n" );

		for( i=0; i < limit; i++ )
		{
			/* >>chng 06 jan 19, from div by eden to div by H0 - want units of cm3 s-1 but
			 * no single collider does this so not possible to get rate coefficient easily
			 * H0 is more appropriate than electron density */
			fprintf( ioOut, "%10.2e", ChargeTransfer[i]/SDIV(dense.xIonDense[ipHYDROGEN][0]) );
		}
		fprintf( ioOut, " CT/n(H0)\n" );

		for( i=0; i < limit; i++ )
		{
			fprintf( ioOut, "%10.2e", ionbal.CotaRate[ion] );
		}
		fprintf( ioOut, " 3body vs Z /ne\n" );

		/* note different upper limit - this routine does grain rec for all ions */
		for( i=0; i < dense.IonHigh[nelem]; i++ )
		{
			fprintf( ioOut, "%10.2e", gv.GrainChTrRate[nelem][i+1][i]/dense.eden );
		}
		fprintf( ioOut, " Grain vs Z /ne\n" );

		for( i=0; i < limit; i++ )
		{
			fprintf( ioOut, "%10.2e", DielRecomRateCoef_HiT[i] );
		}
		fprintf( ioOut, " Burgess vs Z\n" );

		for( i=0; i < limit; i++ )
		{
			fprintf( ioOut, "%10.2e", ionbal.DR_old_rate_coef[nelem][i] );
		}
		fprintf( ioOut, " old Nussbaumer Storey DR vs Z\n" );

		for( i=0; i < limit; i++ )
		{
			fprintf( ioOut, "%10.2e", ionbal.DR_Badnell_rate_coef[nelem][i] );
		}
		fprintf( ioOut, " new Badnell DR vs Z\n" );

		for( i=0; i < limit; i++ )
		{
			fprintf( ioOut, "%10.2e", DielRecomRateCoef_LowT[i] );
		}
		fprintf( ioOut, " low T DR used vs Z\n" );

		/* total recombination rate, with density included - this goes into the matrix */
		for( i=0; i < limit; i++ )
		{
			fprintf( ioOut, "%10.2e", ionbal.RateRecomTot[nelem][i] );
		}
		fprintf( ioOut, 
			" total rec rate (with density) for %s\n", 
			elementnames.chElementSym[nelem] );
		for( i=0; i < limit; i++ )
		{
			fprintf( ioOut, "%10.2e", ionbal.RateRecomTot[nelem][i]/dense.eden );
		}
		fprintf( ioOut, 
			" total rec rate / ne for %s\n\n", 
			elementnames.chElementSym[nelem] );

		/* spill over to next line for many stages of ionization */
		if( dense.IonHigh[nelem] > 11 )
		{
			limit = MIN2(29,dense.IonHigh[nelem]);
			fprintf( ioOut, " R " );
			for( i=11; i < limit; i++ )
			{
				fprintf( ioOut, "%10.2e", dense.eden*ionbal.CotaRate[ion] );
			}
			fprintf( ioOut, "\n" );

			fprintf( ioOut, " B " );
			for( i=11; i < limit; i++ )
			{
				fprintf( ioOut, "%10.2e", DielRecomRateCoef_HiT[i] );
			}
			fprintf( ioOut, "\n" );

			fprintf( ioOut, " NS" );
			for( i=11; i < limit; i++ )
			{
				fprintf( ioOut, "%10.2e", DielRecomRateCoef_LowT[i] );
			}
			fprintf( ioOut, "\n" );

			fprintf( ioOut, "   " );
			for( i=11; i < limit; i++ )
			{
				fprintf( ioOut, "%10.2e", ionbal.RateRecomTot[nelem][i] );
			}
			fprintf( ioOut, "\n\n" );
		}
	}

	/* >>chng 02 nov 09, from -2 to -NISO */
	/*limit = MIN2(nelem-2,dense.IonHigh[nelem]-1);*/
	limit = MIN2(nelem-NISO,dense.IonHigh[nelem]-1);
	for( i=dense.IonLow[nelem]; i <= limit; i++ )
	{
		ASSERT( Heavy.xLyaHeavy[nelem][i] > 0. );
		ASSERT( ionbal.RateRecomTot[nelem][i] > 0. );
	}
	return;
#	undef	DITE
#	undef	DICOEF
}

/*ion_recombAGN generate recombination coefficients for AGN table */
void ion_recombAGN( FILE * io )
{
#	define N1LIM 3
#	define N2LIM 4
	double te1[N1LIM]={ 5000., 10000., 20000.};
	double te2[N2LIM]={ 20000.,50000.,100000.,1e6};
	/* this is boundary between two tables */
	double BreakEnergy = 100./13.0;
	long int nelem, ion , i;
	/* this will hold element symbol + ionization */
	char chString[100],
		chOutput[100];
	/* save temp here	*/
	double TempSave = phycon.te;
	/* save ne here	*/
	double EdenSave = dense.eden;

	DEBUG_ENTRY( "ion_recomb(false,)" );

	dense.eden = 1.;
	/*atmdat_readin();*/

	/* first put header on file */
	fprintf(io,"X+i\\Te");
	for( i=0; i<N1LIM; ++i )
	{
		phycon.te = te1[i];
		fprintf(io,"\t%.0f K",phycon.te);
	}
	fprintf(io,"\n");

	/* now do loop over temp, but add elements */
	for( nelem=ipLITHIUM; nelem<LIMELM; ++nelem )
	{
		/* this list of elements included in the AGN tables is defined in zeroabun.c */
		if( abund.lgAGN[nelem] )
		{
			for( ion=0; ion<=nelem; ++ion )
			{
				ASSERT( Heavy.Valence_IP_Ryd[nelem][ion] > 0.05 );

				if( Heavy.Valence_IP_Ryd[nelem][ion] > BreakEnergy )
					break;

				/* print chemical symbol */
				sprintf(chOutput,"%s", 
					elementnames.chElementSym[nelem]);
				/* some elements have only one letter - this avoids leaving a space */
				if( chOutput[1]==' ' )
					chOutput[1] = chOutput[2];
				/* now ionization stage */
				if( ion==0 )
				{
					sprintf(chString,"0 ");
				}
				else if( ion==1 )
				{
					sprintf(chString,"+ ");
				}
				else
				{
					sprintf(chString,"+%li ",ion);
				}
				strcat( chOutput , chString );
				fprintf(io,"%5s",chOutput );

				for( i=0; i<N1LIM; ++i )
				{
					TempChange(te1[i] , false);
					dense.IonLow[nelem] = 0;
					dense.IonHigh[nelem] = nelem+1;
					if( ConvBase(0) )
						fprintf(ioQQQ,"PROBLEM ConvBase returned error.\n");
					fprintf(io,"\t%.2e",ionbal.RateRecomTot[nelem][ion]);
				}
				fprintf(io,"\n");
			}
			fprintf(io,"\n");
		}
	}

	/* second put header on file */
	fprintf(io,"X+i\\Te");
	for( i=0; i<N2LIM; ++i )
	{
		TempChange(te2[i] , false);
		fprintf(io,"\t%.0f K",phycon.te);
	}
	fprintf(io,"\n");

	/* now do same loop over temp, but add elements */
	for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		/* this list of elements included in the AGN tables is defined in zeroabun.c */
		if( abund.lgAGN[nelem] )
		{
			for( ion=0; ion<=nelem; ++ion )
			{
				ASSERT( Heavy.Valence_IP_Ryd[nelem][ion] > 0.05 );

				if( Heavy.Valence_IP_Ryd[nelem][ion] <= BreakEnergy )
					continue;

				/* print chemical symbol */
				fprintf(io,"%s", 
					elementnames.chElementSym[nelem]);
				/* now ionization stage */
				if( ion==0 )
				{
					fprintf(io,"0 ");
				}
				else if( ion==1 )
				{
					fprintf(io,"+ ");
				}
				else
				{
					fprintf(io,"+%li",ion);
				}

				for( i=0; i<N2LIM; ++i )
				{
					TempChange(te2[i] , false);
					dense.IonLow[nelem] = 0;
					dense.IonHigh[nelem] = nelem+1;
					if( ConvBase(0) )
						fprintf(ioQQQ,"PROBLEM ConvBase returned error.\n");
					fprintf(io,"\t%.2e",ionbal.RateRecomTot[nelem][ion]);
				}
				fprintf(io,"\n");
			}
			fprintf(io,"\n");
		}
	}

	TempChange(TempSave , true);
	dense.eden = EdenSave;
	return;
}
