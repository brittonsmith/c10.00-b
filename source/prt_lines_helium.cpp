/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*lines_helium put He-like iso sequence into line intensity stack */
/*TempInterp interpolates on a grid of values to produce predicted value at current Te.*/
#include "cddefines.h"
#include "dense.h"
#include "prt.h"
#include "helike.h"
#include "iso.h"
#include "atmdat.h"
#include "lines.h"
#include "lines_service.h"
#include "phycon.h"
#include "physconst.h"
#include "taulines.h"
#include "thirdparty.h"
#include "trace.h"

#define NUMTEMPS	22

typedef struct 
{
	/* index for upper and lower levels of line */
	long int ipHi;
	long int ipLo;

	char label[5];

} stdLines;

STATIC void GetStandardHeLines(void);
STATIC double TempInterp2( double* TempArray , double* ValueArray, long NumElements );
STATIC void DoSatelliteLines( long nelem );

static bool lgFirstRun = true;
static double CaABTemps[NUMTEMPS];
static long NumLines;
static double ***CaABIntensity;
static stdLines **CaABLines;

void lines_helium(void)
{
	long ipISO = ipHE_LIKE;
	long int i, nelem, ipHi, ipLo;
	char chLabel[5]="    ";

	long int j;

	double 
	  sum,
	  Pop2_3S,
	  photons_3889_plus_7065 = 0.;

	DEBUG_ENTRY( "lines_helium()" );

	if( trace.lgTrace )
		fprintf( ioQQQ, "   prt_lines_helium called\n" );

	// this can be changed with the atom levels command but must be at
	// least 3.
	ASSERT( iso.n_HighestResolved_max[ipHE_LIKE][ipHELIUM] >= 3 );

	i = StuffComment( "He-like iso-sequence" );
	linadd( 0., (realnum)i , "####", 'i',
		" start He-like iso sequence");

	linadd(MAX2(0.,iso.xLineTotCool[ipHE_LIKE][ipHELIUM]),506,"Clin",'c',
		"  total collisional cooling due to all HeI lines ");

	linadd(MAX2(0.,-iso.xLineTotCool[ipHE_LIKE][ipHELIUM]),506,"Hlin",'h'	,
		"  total collisional heating due to all HeI lines ");

	/* read in Case A and B lines from data file	*/
	if( lgFirstRun )
	{
		GetStandardHeLines();
		lgFirstRun = false;
	}

	/* store labels for all case b HeI lines in case we assert case b 
	 * ipass == -1 only counting number of lines, =0, malloc then set wl */
	static bool lgMustMalloc=true;
	if( LineSave.ipass == 0 && atmdat.nCaseBHeI>0 && lgMustMalloc )
	{
		/* second time through - on ipass=-1 we counted number of lines
		 * atmdat.nCaseBHeI, now create space but only if there are He I lines 
		 * this is not done if He is turned off */
		atmdat.CaseBWlHeI = (realnum*)MALLOC( sizeof(realnum)*atmdat.nCaseBHeI);
		lgMustMalloc=false;
	}
	atmdat.nCaseBHeI = 0;

	/* this is the main printout, where line intensities are entered into the stack */
	for( nelem=ipISO; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem] )
		{
			ASSERT( iso.n_HighestResolved_max[ipHE_LIKE][nelem] >= 3 );

			if( nelem == ipHELIUM )
			{
				double *qTotEff;

				/* >>chng 06 aug 17, all of these from _max to _local */
				/* >>chng 06 dec 21, mistake - change back to _max */
				qTotEff = (double*)MALLOC(sizeof(double)*(unsigned)(iso.numLevels_max[ipHE_LIKE][nelem]) );

				qTotEff[0] = 0.;
				qTotEff[1] = 0.;

				for( i=ipHe2s3S+1; i<iso.numLevels_max[ipHE_LIKE][nelem]-iso.nCollapsed_max[ipHE_LIKE][nelem]; i++ )
				{
					qTotEff[i] = 0.;
					for( j = i; j<iso.numLevels_max[ipHE_LIKE][nelem]-iso.nCollapsed_max[ipHE_LIKE][nelem]; j++ )
					{
						/*if( StatesElemNEW[nelem][nelem-ipHE_LIKE][i].S == 3 )
						{*/
							qTotEff[i] += 
								Transitions[ipHE_LIKE][nelem][j][ipHe2s3S].Coll.ColUL*dense.EdenHCorr*
								iso.Boltzmann[ipHE_LIKE][nelem][j][ipHe2s3S] *
								(double)Transitions[ipHE_LIKE][nelem][j][ipHe2s3S].Hi->g / 
								(double)Transitions[ipHE_LIKE][nelem][j][ipHe2s3S].Lo->g*
								iso.CascadeProb[ipISO][nelem][j][i];
						/*}*/
					}
				}
				
				/* get simple 2^3S pop, assume recombinations in are just 0.75 * case B */
				Pop2_3S = dense.eden*(0.75*iso.RadRec_caseB[ipHE_LIKE][nelem])/
					( Transitions[ipHE_LIKE][nelem][ipHe2s3S][ipHe1s1S].Emis->Aul+ dense.eden*iso.qTot2S[ipISO][nelem]);

				for( i=0; i< NumLines; i++ )
				{
					ipHi = CaABLines[nelem][i].ipHi;
					ipLo = CaABLines[nelem][i].ipLo;

					/* >>chng 06 aug 17, from _max to _local */
					/* >>chng 06 dec 21, mistake - change back to _max */
					if( ipHi <= iso.n_HighestResolved_max[ipHE_LIKE][nelem]*(iso.n_HighestResolved_max[ipHE_LIKE][nelem]+1))
					{
						double intens = TempInterp2( CaABTemps , CaABIntensity[nelem][i], NUMTEMPS );
						intens = pow( 10., intens ) * dense.xIonDense[nelem][nelem+1-ipISO]*dense.eden;
						ASSERT( intens >= 0. );

						linadd( intens,
							Transitions[ipHE_LIKE][nelem][ipHi][ipLo].WLAng,
							CaABLines[nelem][i].label,'i',
							"Case B intensity ");

						if( nMatch("Ca B",CaABLines[nelem][i].label) )
						{
							/* all lines to/from 2^3Pj are stored as lines to/from 2^3P1, so make sure this loop never tries to 
							 * explicitly consider 2^3P0 or 2^3P2 */
							ASSERT( ipLo!=ipHe2p3P0 && ipLo!=ipHe2p3P2 );
							ASSERT( ipHi!=ipHe2p3P0 && ipHi!=ipHe2p3P2 );

							double totBranch = iso.BranchRatio[ipISO][nelem][ipHi][ipLo];
							if( ipLo==4 )
								totBranch += iso.BranchRatio[ipISO][nelem][ipHi][3] + iso.BranchRatio[ipISO][nelem][ipHi][5];

							if( LineSave.ipass < 0 )
								++atmdat.nCaseBHeI;
							else if( LineSave.ipass == 0 )
							{
								/* save wavelengths */
								atmdat.CaseBWlHeI[atmdat.nCaseBHeI] = 
									Transitions[ipHE_LIKE][nelem][ipHi][ipLo].WLAng;
								++atmdat.nCaseBHeI;
							}

							if( ipHi==4 )
							{
								linadd( intens + 
									Pop2_3S*dense.xIonDense[nelem][nelem+1-ipISO]*
									(
									qTotEff[ipHe2p3P0]*iso.BranchRatio[ipISO][nelem][ipHe2p3P0][ipLo]+
									qTotEff[ipHe2p3P1]*iso.BranchRatio[ipISO][nelem][ipHe2p3P1][ipLo]+
									qTotEff[ipHe2p3P2]*iso.BranchRatio[ipISO][nelem][ipHe2p3P2][ipLo]
									)*
									Transitions[ipHE_LIKE][nelem][ipHi][ipLo].EnergyErg,
									Transitions[ipHE_LIKE][nelem][ipHi][ipLo].WLAng,
									"+Col",'i',
									"Case B intensity with collisions included");

							}
							else
							{
								/* chng 05 dec 14, branching ratio was missing here! 
								 * not a big effect because lines with biggest collision
								 * enhancements tend to be dominant decay route from upper level. */
								linadd( intens + 
									Pop2_3S*qTotEff[ipHi]*dense.xIonDense[nelem][nelem+1-ipISO]*totBranch*
									Transitions[ipHE_LIKE][nelem][ipHi][ipLo].EnergyErg,
									Transitions[ipHE_LIKE][nelem][ipHi][ipLo].WLAng,
									"+Col",'i',
									"Case B intensity with collisions included");
							}
						}
					}
					/* Make sure to at least do 4471 */
					else if( ipLo==ipHe2p3P1 && ipHi==ipHe4d3D && nMatch("Ca B",CaABLines[nelem][i].label) )
					{
						double intens = TempInterp2( CaABTemps , CaABIntensity[nelem][i], NUMTEMPS );
						intens = pow( 10., intens ) * dense.xIonDense[nelem][nelem+1-ipISO]*dense.eden;
						ASSERT( intens >= 0. );

						linadd( intens, 4471, CaABLines[nelem][i].label, 'i',
							"Case B intensity ");
					}

				}
				free( qTotEff );
			}

			/* NB NB - low and high must be in this order so that all balmer, paschen,
			 * etc series line up correctly in final printout */
			/* >>chng 01 jun 13, bring 23P lines back together */
			/* two photon is special, not a line and no ipCont array index, add here */
			Transitions[ipHE_LIKE][nelem][ipHe2s1S][ipHe1s1S].Emis->phots = 
				Transitions[ipHE_LIKE][nelem][ipHe2s1S][ipHe1s1S].Emis->Aul*
				StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHe2s1S].Pop*
				Transitions[ipHE_LIKE][nelem][ipHe2s1S][ipHe1s1S].Emis->Pesc;

			Transitions[ipHE_LIKE][nelem][ipHe2s1S][ipHe1s1S].Emis->xIntensity = 
				Transitions[ipHE_LIKE][nelem][ipHe2s1S][ipHe1s1S].Emis->phots*
				Transitions[ipHE_LIKE][nelem][ipHe2s1S][ipHe1s1S].EnergyErg;

			if( LineSave.ipass == 0 )
			{
				/* chIonLbl is function that generates a null terminated 4 char string, of form "C  2" 
				 * the result, chLable, is only used when ipass == 0, can be undefined otherwise */
				/* total two photon emission */
				chIonLbl(chLabel, &Transitions[ipHE_LIKE][nelem][ipHe2s1S][ipHe1s1S]);
			}

			linadd( Transitions[ipHE_LIKE][nelem][ipHe2s1S][ipHe1s1S].Emis->xIntensity , 0,chLabel,'r',
				" two photon continuum ");

			linadd(
				StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHe2s1S].Pop*
				iso.TwoNu_induc_dn[ipHE_LIKE][nelem]*
				Transitions[ipHE_LIKE][nelem][ipHe2s1S][ipHe1s1S].EnergyErg,
				22, chLabel ,'i',
				" induced two photon emission ");

			/* here we will create an entry for the three lines 
			 * coming from 2 3P to 1 1S - the entry called TOTL will
			 * appear before the lines of the multiplet */
			sum = 0.;
			for( i=ipHe2p3P0; i <= ipHe2p3P2; i++ )
			{
				if( Transitions[ipHE_LIKE][nelem][i][ipHe1s1S].ipCont <= 0 )
					continue;

				sum += 
					Transitions[ipHE_LIKE][nelem][i][ipHe1s1S].Emis->Aul*
					StatesElemNEW[nelem][nelem-ipHE_LIKE][i].Pop*
					Transitions[ipHE_LIKE][nelem][i][ipHe1s1S].Emis->Pesc*
					Transitions[ipHE_LIKE][nelem][i][ipHe1s1S].EnergyErg;
			}

			linadd(sum,Transitions[ipHE_LIKE][nelem][ipHe2p3P1][ipHe1s1S].WLAng,"TOTL",'i' ,
				" total emission in He-like intercombination lines from 2p3P to ground ");

			/* set number of levels we want to print, first is default,
			 * only print real levels, second is set with "print line
			 * iso collapsed" command */
			long int nLoop  = iso.numLevels_max[ipHE_LIKE][nelem] - iso.nCollapsed_max[ipHE_LIKE][nelem];
			if( prt.lgPrnIsoCollapsed )
				nLoop  = iso.numLevels_max[ipHE_LIKE][nelem];

			/* now do real permitted lines */
			for( ipLo=0; ipLo < ipHe2p3P0; ipLo++ )
			{
				for( ipHi=ipLo+1; ipHi < nLoop; ipHi++ )
				{
					/* >>chng 01 may 30, do not add fake he-like lines (majority) to line stack */
					/* >>chng 01 dec 11, use variable for smallest A */
					if( Transitions[ipHE_LIKE][nelem][ipHi][ipLo].ipCont < 1 ) 
						continue;

					/* recombine fine-structure lines since the energies are
					 * not resolved anyway.	*/
					if( iso.lgFSM[ipISO] && ( abs(StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].l -
						StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].l)==1 )
						&& (StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].l>1) 
						&& (StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].l>1) 
						&& ( StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].n ==
						StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].n ) )
					{
						/* skip if both singlets. */
						if( (StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].S==1) 
							&& (StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].S==1) )
						{
							continue;
						}
						else if( (StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].S==3) 
							&& (StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].S==3) )
						{

							/* singlet to singlet*/
							Transitions[ipHE_LIKE][nelem][ipHi+1][ipLo+1].Emis->phots = 
								Transitions[ipHE_LIKE][nelem][ipHi][ipLo+1].Emis->Aul*
								StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].Pop*
								Transitions[ipHE_LIKE][nelem][ipHi][ipLo+1].Emis->Pesc +
								Transitions[ipHE_LIKE][nelem][ipHi+1][ipLo+1].Emis->Aul*
								StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi+1].Pop*
								Transitions[ipHE_LIKE][nelem][ipHi+1][ipLo+1].Emis->Pesc;

							Transitions[ipHE_LIKE][nelem][ipHi+1][ipLo+1].Emis->xIntensity = 
								Transitions[ipHE_LIKE][nelem][ipHi+1][ipLo+1].Emis->phots *
								Transitions[ipHE_LIKE][nelem][ipHi+1][ipLo+1].EnergyErg;

							PutLine(&Transitions[ipHE_LIKE][nelem][ipHi+1][ipLo+1],
								" ");

							/* triplet to triplet */
							Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->phots = 
								Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->Aul*
								StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].Pop*
								Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->Pesc +
								Transitions[ipHE_LIKE][nelem][ipHi+1][ipLo].Emis->Aul*
								StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi+1].Pop*
								Transitions[ipHE_LIKE][nelem][ipHi+1][ipLo].Emis->Pesc;

							Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->xIntensity = 
								Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->phots *
								Transitions[ipHE_LIKE][nelem][ipHi][ipLo].EnergyErg;

							PutLine(&Transitions[ipHE_LIKE][nelem][ipHi][ipLo],
								" ");
						}
					}

					else if( ipLo==ipHe2s3S && ipHi == ipHe2p3P0 )
					{
						/* here we will create an entry for the three lines 
						 * coming from 2 3P to 2 3S - the entry called TOTL will
						 * appear before the lines of the multiplet 
						 * for He I this is 10830 */

						realnum av_wl = 0.;
						sum = 0.;
						for( i=ipHe2p3P0; i <= ipHe2p3P2; i++ )
						{
							sum += 
								Transitions[ipHE_LIKE][nelem][i][ipLo].Emis->Aul*
								StatesElemNEW[nelem][nelem-ipHE_LIKE][i].Pop*
								Transitions[ipHE_LIKE][nelem][i][ipLo].Emis->Pesc*
								Transitions[ipHE_LIKE][nelem][i][ipLo].EnergyErg;
							av_wl += Transitions[ipHE_LIKE][nelem][i][ipLo].WLAng;
						}
						av_wl /= 3.;
#						if 0
						{
#						include "elementnames.h"
#						include "prt.h"
						fprintf(ioQQQ,"DEBUG 2P - 2S multiplet wl %s ",
							elementnames.chElementSym[nelem] );
						prt_wl( ioQQQ , av_wl );
						fprintf(ioQQQ,"\n" );
						}
#						endif

						linadd(sum,av_wl,"TOTL",'i',
							"total emission in He-like lines, use average of three line wavelengths " );

						/* also add this with the regular label, so it is correctly picked up by assert case-b command */
						linadd(sum,av_wl,chLabel,'i',
							"total emission in He-like lines, use average of three line wavelengths " );

						/*>>chng 05 sep 8, added the following so that the component
						 * from ipHe2p3P0 is printed, in addition to the total. */

						/* find number of photons escaping cloud */
						Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->phots = 
							Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->Aul*
							StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].Pop*
							Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->Pesc;

						/* now find line intensity */
						Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->xIntensity = 
							Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->phots*
							Transitions[ipHE_LIKE][nelem][ipHi][ipLo].EnergyErg;

						if( iso.lgRandErrGen[ipISO] )
						{
							Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->phots *=
								iso.ErrorFactor[ipISO][nelem][ipHi][ipLo][IPRAD];
							Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->xIntensity *= 
								iso.ErrorFactor[ipISO][nelem][ipHi][ipLo][IPRAD];
						}
						PutLine(&Transitions[ipHE_LIKE][nelem][ipHi][ipLo],
							" ");
					}

					else
					{

						/* find number of photons escaping cloud */
						Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->phots = 
							Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->Aul*
							StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].Pop*
							Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->Pesc;

						/* now find line intensity */
						/* >>chng 01 jan 15, put cast double to force double evaluation */
						Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->xIntensity = 
							Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->phots*
							Transitions[ipHE_LIKE][nelem][ipHi][ipLo].EnergyErg;

						if( iso.lgRandErrGen[ipISO] )
						{
							Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->phots *=
								iso.ErrorFactor[ipISO][nelem][ipHi][ipLo][IPRAD];
							Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->xIntensity *= 
								iso.ErrorFactor[ipISO][nelem][ipHi][ipLo][IPRAD];
						}

						/* 
						fprintf(ioQQQ,"1 loop %li %li %.1f\n", ipLo,ipHi, 
							Transitions[ipHE_LIKE][nelem][ipHi][ipLo].WLAng ); */
						PutLine(&Transitions[ipHE_LIKE][nelem][ipHi][ipLo],
							"total intensity of He-like line");
						{
							/* option to print particulars of some line when called
							 * a prettier print statement is near where chSpin is defined below*/
							enum {DEBUG_LOC=false};
							if( DEBUG_LOC )
							{
								if( nelem==1 && ipLo==0 && ipHi==1 )
								{
									fprintf(ioQQQ,"he1 626 %.2e %.2e \n", 
										Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->TauIn,
										Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->TauTot
										);
								}
							}
						}
					}
				}
			}

			/* this sum will bring together the three lines going to J levels within 23P */
			for( ipHi=ipHe2p3P2+1; ipHi < nLoop; ipHi++ )
			{
				double sumcool , sumheat ,
					save , savecool , saveheat;

				sum = 0;
				sumcool = 0.;
				sumheat = 0.;
				for( ipLo=ipHe2p3P0; ipLo <= ipHe2p3P2; ++ipLo )
				{
					if( Transitions[ipHE_LIKE][nelem][ipHi][ipLo].ipCont <= 0 )
						continue;

					/* find number of photons escaping cloud */
					Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->phots = 
						Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->Aul*
						StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].Pop*
						Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->Pesc;

					/* now find line intensity */
					/* >>chng 01 jan 15, put cast double to force double evaluation */
					Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->xIntensity = 
						Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->phots*
						Transitions[ipHE_LIKE][nelem][ipHi][ipLo].EnergyErg;

					if( iso.lgRandErrGen[ipISO] )
					{
						Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->phots *=
							iso.ErrorFactor[ipISO][nelem][ipHi][ipLo][IPRAD];
						Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->xIntensity *= 
							iso.ErrorFactor[ipISO][nelem][ipHi][ipLo][IPRAD];
					}

					sumcool += Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Coll.cool;
					sumheat += Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Coll.heat;
					sum += Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->xIntensity;
				}

				/* skip non-radiative lines */
				if( Transitions[ipHE_LIKE][nelem][ipHi][ipHe2p3P2].ipCont < 1 ) 
					continue;

				/* this will enter .xIntensity into the line stack */
				save = Transitions[ipHE_LIKE][nelem][ipHi][ipHe2p3P2].Emis->xIntensity;
				savecool = Transitions[ipHE_LIKE][nelem][ipHi][ipHe2p3P2].Coll.cool;
				saveheat = Transitions[ipHE_LIKE][nelem][ipHi][ipHe2p3P2].Coll.heat;

				Transitions[ipHE_LIKE][nelem][ipHi][ipHe2p3P2].Emis->xIntensity = sum;
				Transitions[ipHE_LIKE][nelem][ipHi][ipHe2p3P2].Coll.cool = sumcool;
				Transitions[ipHE_LIKE][nelem][ipHi][ipHe2p3P2].Coll.heat = sumheat;

				/*fprintf(ioQQQ,"2 loop %li %li %.1f\n", ipHe2p3P2,ipHi, 
					Transitions[ipHE_LIKE][nelem][ipHi][ipHe2p3P2].WLAng );*/
				PutLine(&Transitions[ipHE_LIKE][nelem][ipHi][ipHe2p3P2],
					"predicted line, all processes included");

				Transitions[ipHE_LIKE][nelem][ipHi][ipHe2p3P2].Emis->xIntensity = save;
				Transitions[ipHE_LIKE][nelem][ipHi][ipHe2p3P2].Coll.cool = savecool;
				Transitions[ipHE_LIKE][nelem][ipHi][ipHe2p3P2].Coll.heat = saveheat;
			}
			for( ipLo=ipHe2p3P2+1; ipLo < nLoop-1; ipLo++ )
			{
				for( ipHi=ipLo+1; ipHi < nLoop; ipHi++ )
				{
					/* skip non-radiative lines */
					if( Transitions[ipHE_LIKE][nelem][ipHi][ipLo].ipCont < 1 ) 
						continue;

					/* find number of photons escaping cloud */
					Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->phots = 
						Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->Aul*
						StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].Pop*
						Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->Pesc;

					/* now find line intensity */
					/* >>chng 01 jan 15, put cast double to force double evaluation */
					Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->xIntensity = 
						Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->phots*
						Transitions[ipHE_LIKE][nelem][ipHi][ipLo].EnergyErg;

					if( iso.lgRandErrGen[ipISO] )
					{
						Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->phots *=
							iso.ErrorFactor[ipISO][nelem][ipHi][ipLo][IPRAD];
						Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->xIntensity *= 
							iso.ErrorFactor[ipISO][nelem][ipHi][ipLo][IPRAD];
					}

					/* this will enter .xIntensity into the line stack */
					PutLine(&Transitions[ipHE_LIKE][nelem][ipHi][ipLo],
						"predicted line, all processes included");
				}
			}

			/* Now put the satellite lines in */
			if( iso.lgDielRecom[ipISO] )
				DoSatelliteLines(nelem);
		}
	}

	if( iso.n_HighestResolved_max[ipHE_LIKE][ipHELIUM] >= 4 &&
		( iso.n_HighestResolved_max[ipH_LIKE][ipHYDROGEN] + iso.nCollapsed_max[ipH_LIKE][ipHYDROGEN] ) >=8 )
	{
		/* For info only, add the total photon flux in 3889 and 7065,
		* and 3188, 4713, and 5876. */
		photons_3889_plus_7065 =
			/* to 2p3P2 */
			Transitions[ipHE_LIKE][ipHELIUM][ipHe3s3S][ipHe2p3P2].Emis->xIntensity/
			Transitions[ipHE_LIKE][ipHELIUM][ipHe3s3S][ipHe2p3P2].EnergyErg +
			Transitions[ipHE_LIKE][ipHELIUM][ipHe3d3D][ipHe2p3P2].Emis->xIntensity/
			Transitions[ipHE_LIKE][ipHELIUM][ipHe3d3D][ipHe2p3P2].EnergyErg +
			Transitions[ipHE_LIKE][ipHELIUM][ipHe4s3S][ipHe2p3P2].Emis->xIntensity/
			Transitions[ipHE_LIKE][ipHELIUM][ipHe4s3S][ipHe2p3P2].EnergyErg +
			/* to 2p3P1 */
			Transitions[ipHE_LIKE][ipHELIUM][ipHe3s3S][ipHe2p3P1].Emis->xIntensity/
			Transitions[ipHE_LIKE][ipHELIUM][ipHe3s3S][ipHe2p3P1].EnergyErg +
			Transitions[ipHE_LIKE][ipHELIUM][ipHe3d3D][ipHe2p3P1].Emis->xIntensity/
			Transitions[ipHE_LIKE][ipHELIUM][ipHe3d3D][ipHe2p3P1].EnergyErg +
			Transitions[ipHE_LIKE][ipHELIUM][ipHe4s3S][ipHe2p3P1].Emis->xIntensity/
			Transitions[ipHE_LIKE][ipHELIUM][ipHe4s3S][ipHe2p3P1].EnergyErg +
			/* to 2p3P0 */
			Transitions[ipHE_LIKE][ipHELIUM][ipHe3s3S][ipHe2p3P0].Emis->xIntensity/
			Transitions[ipHE_LIKE][ipHELIUM][ipHe3s3S][ipHe2p3P0].EnergyErg +
			Transitions[ipHE_LIKE][ipHELIUM][ipHe3d3D][ipHe2p3P0].Emis->xIntensity/
			Transitions[ipHE_LIKE][ipHELIUM][ipHe3d3D][ipHe2p3P0].EnergyErg +
			Transitions[ipHE_LIKE][ipHELIUM][ipHe4s3S][ipHe2p3P0].Emis->xIntensity/
			Transitions[ipHE_LIKE][ipHELIUM][ipHe4s3S][ipHe2p3P0].EnergyErg +
			/* to 2s3S */
			Transitions[ipHE_LIKE][ipHELIUM][ipHe3p3P][ipHe2s3S].Emis->xIntensity/
			Transitions[ipHE_LIKE][ipHELIUM][ipHe3p3P][ipHe2s3S].EnergyErg +
			Transitions[ipHE_LIKE][ipHELIUM][ipHe4p3P][ipHe2s3S].Emis->xIntensity/
			Transitions[ipHE_LIKE][ipHELIUM][ipHe4p3P][ipHe2s3S].EnergyErg;

		long upperIndexofH8 = iso.QuantumNumbers2Index[ipH_LIKE][ipHYDROGEN][8][1][2];

		/* Add in photon flux of H8 3889 */
		photons_3889_plus_7065 += 
			Transitions[ipH_LIKE][ipHYDROGEN][upperIndexofH8][1].Emis->xIntensity/
			Transitions[ipH_LIKE][ipHYDROGEN][upperIndexofH8][1].EnergyErg;

		/* now multiply by ergs of normalization line, so that relative flux of
		* this line will be ratio of photon fluxes. */
		photons_3889_plus_7065 *= (ERG1CM*1.e8)/LineSave.WavLNorm;
		linadd( photons_3889_plus_7065, 3889., "Pho+", 'i',
			"photon sum given in Porter et al. 2007 (astro-ph/0611579)");
	}

	/* ====================================================
	 * end helium
	 * ====================================================*/

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "   lines_helium returns\n" );
	}
	return;
}


STATIC void GetStandardHeLines(void)
{
	FILE *ioDATA;
	bool lgEOL, lgHIT;
	long i, i1, i2, j, nelem;

#	define chLine_LENGTH 1000
	char chLine[chLine_LENGTH];

	DEBUG_ENTRY( "GetStandardHeLines()" );

	if( trace.lgTrace )
		fprintf( ioQQQ," lines_helium opening he1_case_ab.dat\n");

	ioDATA = open_data( "he1_case_ab.dat", "r" );

	/* check that magic number is ok */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " lines_helium could not read first line of he1_case_ab.dat.\n");
		cdEXIT(EXIT_FAILURE);
	}
	i = 1;
	/* first number is magic number, second is number of lines in file	*/
	i1 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
	i2 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
	NumLines = i2;

	/* the following is to check the numbers that appear at the start of he1_case_ab.dat */
	if( i1 !=CASEABMAGIC )
	{
		fprintf( ioQQQ, 
			" lines_helium: the version of he1_case_ab.dat is not the current version.\n" );
		fprintf( ioQQQ, 
			" lines_helium: I expected to find the number %i and got %li instead.\n" ,
			CASEABMAGIC, i1 );
		fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
		cdEXIT(EXIT_FAILURE);
	}

	/* get the array of temperatures */
	lgHIT = false;
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
	{
		/* only look at lines without '#' in first col */
		if( chLine[0] != '#')
		{
			lgHIT = true;
			j = 0;
			lgEOL = false;
			i = 1;
			while( !lgEOL && j < NUMTEMPS)
			{
				CaABTemps[j] = FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
				++j;
			}
			break;
		}
	}

	if( !lgHIT )
	{
		fprintf( ioQQQ, " lines_helium could not find line of temperatures in he1_case_ab.dat.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* create space for array of values, if not already done */
	CaABIntensity = (double ***)MALLOC(sizeof(double **)*(unsigned)LIMELM );
	/* create space for array of values, if not already done */
	CaABLines = (stdLines **)MALLOC(sizeof(stdLines *)*(unsigned)LIMELM );

	for( nelem=ipHELIUM; nelem<LIMELM; ++nelem )
	{
		/** \todo	2	- this structure is currently only used for helium itself...
		 * stuff numbers in for other elements, or drop the [nelem] dimension off
		 * of CaABLines	*/
		if( nelem != ipHELIUM )
			continue;

		/* only grab core for elements that are turned on */
		if( nelem == ipHELIUM || dense.lgElmtOn[nelem])
		{
			/* create space for array of values, if not already done */
			CaABIntensity[nelem] = (double **)MALLOC(sizeof(double *)*(unsigned)(i2) );
			CaABLines[nelem] = (stdLines *)MALLOC(sizeof(stdLines )*(unsigned)(i2) );

			/* avoid allocating 0 bytes, some OS return NULL pointer, PvH 
			CaABIntensity[nelem][0] = NULL;*/
			for( j = 0; j < i2; ++j )
			{
				CaABIntensity[nelem][j] = (double *)MALLOC(sizeof(double)*(unsigned)NUMTEMPS );

				CaABLines[nelem][j].ipHi = -1;
				CaABLines[nelem][j].ipLo = -1;
				strcpy( CaABLines[nelem][j].label , "    " );

				for( i=0; i<NUMTEMPS; ++i )
				{
					CaABIntensity[nelem][j][i] = 0.;
				}
			}
		}
	}

	/* now read in the case A and B line data */
	lgHIT = false;
	nelem = ipHELIUM;
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
	{
		static long line = 0;
		char *chTemp;

		/* only look at lines without '#' in first col */
		if( (chLine[0] == ' ') || (chLine[0]=='\n') )
			break;
		if( chLine[0] != '#')
		{
			lgHIT = true;

			/* get lower and upper level index */
			j = 1;
			/* the first number is the wavelength, which is not used
			 * for anything, but must skip over it. */
			i1 = (long)FFmtRead(chLine,&j,INPUT_LINE_LENGTH,&lgEOL);
			CaABLines[nelem][line].ipLo = (long)FFmtRead(chLine,&j,INPUT_LINE_LENGTH,&lgEOL);
			CaABLines[nelem][line].ipHi = (long)FFmtRead(chLine,&j,INPUT_LINE_LENGTH,&lgEOL);

			ASSERT( CaABLines[nelem][line].ipLo < CaABLines[nelem][line].ipHi );
			/*if( CaABLines[nelem][line].ipHi >= iso.numLevels_max[ipHE_LIKE][nelem] )
				continue;*/

			chTemp = chLine;
			/* skip over 4 tabs to start of data */
			for( i=0; i<3; ++i )
			{
				if( (chTemp = strchr_s( chTemp, '\t' )) == NULL )
				{
					fprintf( ioQQQ, " lines_helium no init case A and B\n" );
					cdEXIT(EXIT_FAILURE);
				}
				++chTemp;
			}

			strncpy( CaABLines[nelem][line].label, chTemp , 4 );
			CaABLines[nelem][line].label[4] = 0;

			for( i=0; i<NUMTEMPS; ++i )
			{
				double b;
				if( (chTemp = strchr_s( chTemp, '\t' )) == NULL )
				{
					fprintf( ioQQQ, " lines_helium could not scan case A and B lines, current indices: %li %li\n",
						CaABLines[nelem][line].ipHi,
						CaABLines[nelem][line].ipLo );
					cdEXIT(EXIT_FAILURE);
				}
				++chTemp;
				sscanf( chTemp , "%le" , &b );
				CaABIntensity[nelem][line][i] = b;
			}
			line++;
		}
	}

	/* close the data file */
	fclose( ioDATA );
	return;
}

/** \todo	there is a virtually identical routine in helike_recom.cpp -> combine */
STATIC double TempInterp2( double* TempArray , double* ValueArray, long NumElements )
{
	long int ipTe=-1;
	double rate = 0.;
	long i0;

	DEBUG_ENTRY( "TempInterp2()" );

	ASSERT( fabs( 1. - (double)phycon.alogte/log10(phycon.te) ) < 0.0001 );

	/* te totally unknown */
	if( phycon.alogte <= TempArray[0] )
	{
		return ValueArray[0];
	}
	else if( phycon.alogte >= TempArray[NumElements-1] )
	{
		return ValueArray[NumElements-1];
	}

	/* now search for temperature */
	ipTe = hunt_bisect( TempArray, NumElements, phycon.alogte );			

	ASSERT( (ipTe >=0) && (ipTe < NumElements-1)  );

	/* Do a four-point interpolation */
	const int ORDER = 3; /* order of the fitting polynomial */
	i0 = max(min(ipTe-ORDER/2,NumElements-ORDER-1),0);
	rate = lagrange( &TempArray[i0], &ValueArray[i0], ORDER+1, phycon.alogte );

	return rate;
}

/** \todo	2	say where these come from	*/	
/* For double-ionization discussions, see Lindsay, Rejoub, & Stebbings 2002	*/
/* Also read Itza-Ortiz, Godunov, Wang, and McGuire 2001.	*/
STATIC void DoSatelliteLines( long nelem )
{
	long ipISO = ipHE_LIKE;
	
	DEBUG_ENTRY( "DoSatelliteLines()" );

	ASSERT( iso.lgDielRecom[ipISO] && dense.lgElmtOn[nelem] );

	for( long i=0; i < iso.numLevels_max[ipISO][nelem] - iso.nCollapsed_max[ipISO][nelem]; i++ )
	{
		double dr_rate = iso.DielecRecomb[ipISO][nelem][i];

		SatelliteLines[ipISO][nelem][i].Emis->phots = 
			dr_rate * dense.eden * dense.xIonDense[nelem][nelem+1-ipISO];
		
		SatelliteLines[ipISO][nelem][i].Emis->xIntensity = 
			SatelliteLines[ipISO][nelem][i].Emis->phots * ERG1CM * SatelliteLines[ipISO][nelem][i].EnergyWN;
		SatelliteLines[ipISO][nelem][i].Emis->pump = 0.;

		PutLine( &SatelliteLines[ipISO][nelem][i], "iso satellite line" );
	}

	return;
}
