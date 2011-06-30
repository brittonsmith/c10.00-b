/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ConvIterCheck check whether model has converged or whether more iterations
 * are needed - implements the iter to converg comnd */
#include "cddefines.h"
#include "taulines.h"
#include "iso.h"
#include "phycon.h"
#include "cddrive.h"
#include "mole.h"
#include "elementnames.h"
#include "dynamics.h"
#include "stopcalc.h"
#include "dense.h"
#include "iterations.h"
#include "colden.h"
#include "save.h"
#include "rt.h"
#include "conv.h"

/*ConvIterCheck check whether model has converged or whether more iterations
 * are needed - implements the iterate to convergence command */
void ConvIterCheck( void )
{
	bool lgConverged;
	long int nelem, 
		i,
		ipISO,
		ipHi, ipLo;

	DEBUG_ENTRY( "ConvIterCheck()" );

	/* =======================================================================*/
	/* this is an option to keep iterating until it converges
	 * iterate to convergence option
	 * autocv is percentage difference in optical depths allowed,
	 * =0.20 in block data
	 * checking on Ly and Balmer lines */
	/*>>chng 04 oct 19, promote loop to do all iso-electronic series */
	lgConverged = true;
	strcpy( conv.chNotConverged, "Converged!" );

	// set up intensities used to converge outward intensity of Hb
	static double HbFracOutOld=-1. , HbFracOutNew=-1.;
	HbFracOutOld = HbFracOutNew;
	double a, total, BeamedIn;
	long int ipTotal = cdLine( "TOTL" , 4861.36f , &a , &total );
	long int ipInwd = cdLine( "INWD" , 4861.36f , &a , &BeamedIn );
	HbFracOutNew = 1. - pow(10. , (BeamedIn-total));
	ASSERT( iteration == 1 || (HbFracOutNew>=0 && HbFracOutNew<=1.) );
	// this disables the test on the outward Hb
	ipInwd = -1;

	bool lgReasonGiven = false;
	if( iteration > 1 && conv.lgAutoIt )
	{
		if( nzone>3 && ipInwd>=0 && ipTotal>=0 )
		{
			// check whether outward intensity of Hb has converged
			if( fabs(HbFracOutNew-HbFracOutOld)/HbFracOutNew> conv.autocv )
			{
				lgConverged = false;
				sprintf( conv.chNotConverged, "change in outward Hb");
				if( save.lgPunConv )
				{
					lgReasonGiven = true;
					fprintf( save.ipPunConv, " Change in outward Hb, "
						"old=%.3e new=%.3e \n",
						HbFracOutOld , HbFracOutNew);
				}
			}
		}
		for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			for( nelem=ipISO; nelem < LIMELM; nelem++ )
			{
				if( dense.lgElmtOn[nelem] )
				{
					/* now check if major subordinate line is converged - for H-like this will
					 * be Ha, and for He-like, the 23P - 23S transition - this will not work for
					 * NISO > 2 so must check against this */
					if(ipISO==ipH_LIKE )
					{
						ipHi = ipH3p;
						ipLo = ipH2s;
					}
					else if( ipISO==ipHE_LIKE )
					{
						ipHi = ipHe2p3P2;
						ipLo = ipHe2s3S;
					}
					else
						/* fails when NISO increased, add more sequences */
						TotalInsanity();

					/* check both H-alpha and Ly-alpha for all species - 
					 * only if Balmer lines thick 
					 * so check if Ha optical depth significant */
					if( Transitions[ipISO][nelem][ipHi][ipLo].Emis->TauIn > 0.5 )
					{
						/* test if Lya converged, nLyaLevel is upper level of Lya for iso seq */
						if( fabs(Transitions[ipISO][nelem][iso.nLyaLevel[ipISO]][0].Emis->TauTot -
							 Transitions[ipISO][nelem][iso.nLyaLevel[ipISO]][0].Emis->TauIn*rt.DoubleTau) > 
						    conv.autocv*fabs(Transitions[ipISO][nelem][iso.nLyaLevel[ipISO]][0].Emis->TauIn*rt.DoubleTau) )
						{
							/* not converged to within AUTOCV, normally 15 percent */
							lgConverged = false;

							/* for iterate to convergence, print reason why it was not converged 
							* on 3rd and higher iterations */
							sprintf( conv.chNotConverged, "%s-like Lya",elementnames.chElementSym[ipISO] );

							if( save.lgPunConv )
							{
								lgReasonGiven = true;
								fprintf( save.ipPunConv, " %s-like Lya thick, "
									"nelem= %s iteration %li old %.3e new %.3e \n",
									elementnames.chElementSym[ipISO] ,
									elementnames.chElementSym[nelem], 
									iteration,
									Transitions[ipISO][nelem][iso.nLyaLevel[ipISO]][0].Emis->TauTot ,
									Transitions[ipISO][nelem][iso.nLyaLevel[ipISO]][0].Emis->TauIn);
							}
						}

						if( fabs(Transitions[ipISO][nelem][ipHi][ipLo].Emis->TauTot -
							 Transitions[ipISO][nelem][ipHi][ipLo].Emis->TauIn*rt.DoubleTau) >
						    conv.autocv*fabs(Transitions[ipISO][nelem][ipHi][ipLo].Emis->TauIn*rt.DoubleTau) )
						{
							/* not converged to within AUTOCV, normally 15 percent */
							lgConverged = false;

							/* for iterate to convergence, print reason why it was not converged 
							* on 3rd and higher iterations */
							sprintf( conv.chNotConverged, "%s-like subord",elementnames.chElementSym[ipISO] );

							if( save.lgPunConv )
							{
								lgReasonGiven = true;
								fprintf( save.ipPunConv, " %s-like subord, nelem= %s iteration %li old %.3e new %.3e \n" ,
									elementnames.chElementSym[ipISO],
									elementnames.chElementSym[nelem], 
									iteration,
									Transitions[ipISO][nelem][ipHi][ipLo].Emis->TauTot ,
									Transitions[ipISO][nelem][ipHi][ipLo].Emis->TauIn	);
							}
						}
					}
				}
			}
		}

		/* >>chng 03 sep 07, add this test */
		/* check on changes in major column densities */
		for( i=0; i<NCOLD; ++i )
		{
			/* was the species column density significant relative to
			 * the total H column density, and was its abundance changing? */
			if( colden.colden[i]/colden.colden[ipCOL_HTOT] > 1e-5 &&
			    fabs(colden.colden_old[i]-colden.colden[i]) > conv.autocv*colden.colden[i] )
			{
				/* not converged to within conv.autocv, normally 20 percent */
				lgConverged = false;

				/* for iterate to convergence, print reason why it was not converged 
				 * on 3rd and higher iterations */
				strcpy( conv.chNotConverged, "H mole col" );

				if( save.lgPunConv )
				{
					lgReasonGiven = true;
					fprintf( save.ipPunConv, " H mole col species %li iteration %li old %.2e new %.2e H col den %.2e\n",
						i,iteration,
						colden.colden_old[i],
						colden.colden[i],
						colden.colden[ipCOL_HTOT] );
				}
			}
		}

		/* >>chng 03 sep 07, add this test */
		/* check on changes in major column densities */
		for( i=0; i<mole.num_comole_calc; ++i )
		{
 			if(COmole[i]->n_nuclei == 1)
 				continue;

 			/* was the species abundance and changing? */
 			if( COmole[i]->hevcol/colden.colden[ipCOL_HTOT] > 1e-5 &&
			    fabs(COmole[i]->hevcol_old-COmole[i]->hevcol) > conv.autocv*COmole[i]->hevcol )
			{
				/* not converged to within conv.autocv, normally 20 percent */
				lgConverged = false;

				/* for iterate to convergence, print reason why it was not converged 
				 * on 3rd and higher iterations */
				strcpy( conv.chNotConverged, "CO mol col" );
				/*fprintf(ioQQQ,"debugggreset\t CO mole %li %li %.2e %.2e\n",
					i,iteration,COmole[i]->hevcol_old,COmole[i]->hevcol);*/

				if( save.lgPunConv )
				{
					lgReasonGiven = true;
					fprintf( save.ipPunConv, "CO mol col, old:%.3e new:%.3e\n" ,
						COmole[i]->hevcol_old ,
						COmole[i]->hevcol );
				}
			}
		}

		/* check on dynamical convergence in wind model with negative velocity */
		if( dynamics.lgAdvection )
		{
			/* >>chng 02 nov 29, as per Will Henney email */
			if( dynamics.convergence_error > conv.autocv*dynamics.error_scale2*dynamics.convergence_tolerance ||
			    dynamics.discretization_error > conv.autocv*dynamics.error_scale2 )
			{
				lgConverged = false;
				/* for iterate to convergence, print reason why it was not converged 
				 * on 3rd and higher iterations */
				strcpy( conv.chNotConverged, "Dynamics  " );
				if( save.lgPunConv )
				{
					lgReasonGiven = true;
					fprintf( save.ipPunConv, " Dynamics\n" );
				}
			}
		}

		if( save.lgPunConv && lgConverged )
		{
			lgReasonGiven = true;
			fprintf( save.ipPunConv, " converged\n" );
		}

		/* lower limit to number of iterations if converged */
		if( lgConverged )
			iterations.itermx = MIN2(iterations.itermx,iteration);

		/* >>chng 96 dec 20, moved following to within if on lgAutoIt
		 * this is test for stopping on first zone */
		if( phycon.te < StopCalc.TempLoStopZone && nzone == 1 )
		{
			lgConverged = true;
			strcpy( conv.chNotConverged, "          " );
			iterations.itermx = MIN2(iterations.itermx,iteration);
		}
		// fails if we have not fully implemented save convergence reason
		ASSERT( !save.lgPunConv || lgReasonGiven );
	}
	return;
}
