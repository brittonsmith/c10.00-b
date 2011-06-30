/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
* others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "phycon.h"
#include "taulines.h"
#include "mole.h"
#include "atoms.h"
#include "string.h"
#include "thirdparty.h"
#include "dense.h"
#include "conv.h"
#include "h2.h"
#include "physconst.h"
#include "secondaries.h"
#include "thermal.h"
#include "cooling.h"
#include "lines_service.h"

void states_popfill( void);
STATIC double LeidenCollRate(long, long, long, long,double);
STATIC double CHIANTI_Upsilon(long, long, long, long,double);

#define DEBUGSTATE false

static double *g, *ex, *pops, *depart, *source, *sink;
static double **AulEscp, **col_str, **AulDest, **AulPump, **CollRate;

/*Solving for the level populations*/

void dBase_solve(void )
{
	realnum abund;
	DEBUG_ENTRY( "dBase_solve()" );
	static bool lgFirstPass = true;
	static long maxNumLevels = 1;
	double totalHeating = 0.;

	if( nSpecies==0 )
		return;

	if( lgFirstPass )
	{
		for( long ipSpecies=0; ipSpecies<nSpecies; ipSpecies++ )
			maxNumLevels = MAX2( maxNumLevels, Species[ipSpecies].numLevels_max );

		g = (double *)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double));
		ex = (double *)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double));
		pops = (double *)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double));
		depart = (double *)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double));
		source = (double *)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double));
		sink = (double *)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double));

		AulEscp = (double **)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double *)); 
		col_str = (double **)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double *)); 
		AulDest = (double **)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double *)); 
		AulPump = (double **)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double *));  
		CollRate = (double **)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double *)); 

		for( long j=0; j< maxNumLevels; j++ )
		{
			AulEscp[j] = (double *)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double)); 
			col_str[j] = (double *)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double)); 
			AulDest[j] = (double *)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double)); 
			AulPump[j] = (double *)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double));  
			CollRate[j] = (double *)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double)); 
		}

		lgFirstPass = false;
	}

	// zero all of these values
	memset( g, 0, (unsigned long)(maxNumLevels)*sizeof(double) );
	memset( ex, 0, (unsigned long)(maxNumLevels)*sizeof(double) );
	memset( pops, 0, (unsigned long)(maxNumLevels)*sizeof(double) );
	memset( depart, 0, (unsigned long)(maxNumLevels)*sizeof(double) );
	memset( source, 0, (unsigned long)(maxNumLevels)*sizeof(double) );
	memset( sink, 0, (unsigned long)(maxNumLevels)*sizeof(double) );
	for( long j=0; j< maxNumLevels; j++ )
	{
		memset( AulEscp[j], 0, (unsigned long)(maxNumLevels)*sizeof(double) );
		memset( col_str[j], 0, (unsigned long)(maxNumLevels)*sizeof(double) );
		memset( AulDest[j], 0, (unsigned long)(maxNumLevels)*sizeof(double) );
		memset( AulPump[j], 0, (unsigned long)(maxNumLevels)*sizeof(double) );
		memset( CollRate[j], 0, (unsigned long)(maxNumLevels)*sizeof(double) );
	}


	for( long ipSpecies=0; ipSpecies<nSpecies; ipSpecies++ )
	{
		const char *spName = Species[ipSpecies].chLabel;
		double fcolldensity[9];
		double cooltl, coolder;
		double fupsilon;
		long intCollNo;
		int nNegPop;
		bool lgZeroPop, lgDeBug = false;
		double totaldCdT = 0.;

		Species[ipSpecies].CoolTotal = 0.;

#if 0
		//limit for now to small number of levels
		Species[ipSpecies].numLevels_local = MIN2( Species[ipSpecies].numLevels_local, 10 );
#endif

		/* first find current density (cm-3) of species */
		if( Species[ipSpecies].lgMolecular )
		{
			struct molecule *SpeciesCurrent;
			/** \todo	0	this pointer should be cached one time, and the species
			 * removed from the list if it is not computed */
			if( (SpeciesCurrent = findspecies(Species[ipSpecies].chLabel)) == &null_mole )
			{
				/* did not find the species - print warning for now */
				if( !conv.nTotalIoniz )
					fprintf(ioQQQ," PROBLEM dBase_solve did not find molecular species %li\n",ipSpecies);
			}
			abund = SpeciesCurrent->hevmol;
		}
		else
		{
			/* an atom or ion */
			ASSERT( dBaseStates[ipSpecies][0].nelem<=LIMELM && dBaseStates[ipSpecies][0].IonStg<=dBaseStates[ipSpecies][0].nelem+1 );
			abund = dense.xIonDense[ dBaseStates[ipSpecies][0].nelem-1 ][ dBaseStates[ipSpecies][0].IonStg-1 ];
		}

		abund *= Species[ipSpecies].fracType * Species[ipSpecies].fracIsotopologue;

		// initialization at start of each iteration
		if( conv.nTotalIoniz == 0)
			Species[ipSpecies].lgActive = true;

		bool lgMakeInActive = (abund <= 1e-20 * dense.xNucleiTotal);
		if( lgMakeInActive && Species[ipSpecies].lgActive )
		{
			// zero out populations and intensities, if previously not set
			dBaseStates[ipSpecies][0].Pop = 0.;
			for(long ipHi = 1; ipHi < Species[ipSpecies].numLevels_max; ipHi++ )
			{	
				dBaseStates[ipSpecies][ipHi].Pop = 0.;
				transition *t = &dBaseTrans[ipSpecies][ipHi][0];
				for(long ipLo = 0; ipLo < ipHi; ipLo++ )
				{
					t->Emis->phots = 0.;
					t->Emis->xIntensity = 0.;
					t->Coll.col_str = 0.;
					t->Coll.cool = 0.;
					t->Coll.heat = 0.;
					++t;
				}
			}
			Species[ipSpecies].lgActive = false;
		}

		if( !lgMakeInActive )
			Species[ipSpecies].lgActive = true;

		if( !Species[ipSpecies].lgActive )
			continue;

		// we always hit phase first, reset number of levels
		if( conv.lgSearch )
			Species[ipSpecies].numLevels_local = Species[ipSpecies].numLevels_max;
		
		for( long ipLo = 0; ipLo < Species[ipSpecies].numLevels_local; ipLo++ )
		{
			/* statistical weights & Excitation Energies*/
			g[ipLo] = dBaseStates[ipSpecies][ipLo].g ;
			// parts of the code assert that ground is at zero energy - this is
			// not true for the stored molecular data - so rescale to zero
			ex[ipLo] = dBaseStates[ipSpecies][ipLo].energy.WN() - 
					  dBaseStates[ipSpecies][0].energy.WN();
			/* zero some rates */	
			source[ipLo] = 0.;
			sink[ipLo] = 0.;
			AulEscp[ipLo][ipLo] = 0.;
			AulDest[ipLo][ipLo] = 0.;
			AulPump[ipLo][ipLo] = 0.;
		}

		// non-zero was due to roundoff errors on 32-bit
		if( ex[0] <= dBaseStates[ipSpecies][0].energy.WN()* 10. *DBL_EPSILON )
			ex[0] = 0.;
		else
			TotalInsanity();

		for( long ipHi=1; ipHi<Species[ipSpecies].numLevels_local; ipHi++ )
		{
			for( long ipLo=0; ipLo<ipHi; ipLo++ )
			{
				if( dBaseTrans[ipSpecies][ipHi][ipLo].ipCont > 0 )
				{	
					AulEscp[ipHi][ipLo] = dBaseTrans[ipSpecies][ipHi][ipLo].Emis->Aul*
						(dBaseTrans[ipSpecies][ipHi][ipLo].Emis->Pesc +
						dBaseTrans[ipSpecies][ipHi][ipLo].Emis->Pelec_esc);
					AulDest[ipHi][ipLo] = dBaseTrans[ipSpecies][ipHi][ipLo].Emis->Aul*
						dBaseTrans[ipSpecies][ipHi][ipLo].Emis->Pdest;
					AulPump[ipLo][ipHi] = dBaseTrans[ipSpecies][ipHi][ipLo].Emis->pump;
				}
				else
				{
					AulEscp[ipHi][ipLo] = SMALLFLOAT;
					AulDest[ipHi][ipLo] = SMALLFLOAT;
					AulPump[ipLo][ipHi] = SMALLFLOAT;
				}

				AulEscp[ipLo][ipHi] = 0.;
				AulDest[ipLo][ipHi] = 0.;
				AulPump[ipHi][ipLo] = 0.;
			}
		}

		/*Setting all the collision strengths and collision rate to zero*/
		for( long ipHi= 0; ipHi<Species[ipSpecies].numLevels_local; ipHi++)
		{
			for( long ipLo= 0; ipLo<Species[ipSpecies].numLevels_local; ipLo++)
			{
				col_str[ipHi][ipLo] = 0.;
				CollRate[ipHi][ipLo] = 0.;
			}
		}

		/*Updating the CollRatesArray*/
		/*You need to do this for a species indexed by ipSpecies*/
		for( intCollNo=0; intCollNo<NUM_COLLIDERS; intCollNo++)
		{
			for( long ipHi=1; ipHi<Species[ipSpecies].numLevels_local; ipHi++)
			{
				for( long ipLo=0; ipLo<ipHi; ipLo++)
				{
					/* molecule */
					if( Species[ipSpecies].lgMolecular )
					{
						if(CollRatesArray[ipSpecies][intCollNo]!=NULL)
						{
							/*using the collision rate coefficients directly*/
							CollRatesArray[ipSpecies][intCollNo][ipHi][ipLo] = 
								LeidenCollRate(ipSpecies, intCollNo, ipHi, ipLo, phycon.te);
						}
					}
					/* atom or ion */
					else
					{
						if(CollRatesArray[ipSpecies][intCollNo]!=NULL)
						{
							fupsilon = CHIANTI_Upsilon(ipSpecies, intCollNo, ipHi, ipLo, phycon.te);

							/* NB NB - if proton colliders, the upsilons returned here are actually already rate coefficients. */
							/* these are designated by a collider index and a transition type */
							if( intCollNo==1 && AtmolCollSplines[ipSpecies][ipHi][ipLo][intCollNo].intTranType==6 )
							{
								CollRatesArray[ipSpecies][intCollNo][ipHi][ipLo] = fupsilon;
							}
							else
							{
								/* convert the collision strength to a collision rate coefficient */
								/*This formula converting collision strength to collision rate coefficent works fine for the electrons*/
								/*For any other collider the mass would be different*/
								ASSERT( intCollNo == 0 );
								CollRatesArray[ipSpecies][intCollNo][ipHi][ipLo] = (COLL_CONST*fupsilon)/(g[ipHi]*phycon.sqrte);
							}
						}
					}

					/* now get the corresponding excitation rate */
					if(CollRatesArray[ipSpecies][intCollNo]!=NULL)
					{
						CollRatesArray[ipSpecies][intCollNo][ipLo][ipHi] = 
							CollRatesArray[ipSpecies][intCollNo][ipHi][ipLo] * g[ipHi] / g[ipLo] *
							sexp( dBaseTrans[ipSpecies][ipHi][ipLo].EnergyK / phycon.te );
					}
				}
			}
		}
		/*Get the densities separately*/
		/*Electron*/
		fcolldensity[0] = dense.eden;
		/*Proton*/
		fcolldensity[1] = dense.xIonDense[ipHYDROGEN][1];
		/*Atomic Hydrogen*/
		fcolldensity[2] = dense.xIonDense[ipHYDROGEN][0];
		/*Helium*/
		fcolldensity[3] = dense.xIonDense[ipHELIUM][0];
		/*Helium+*/
		fcolldensity[4] = dense.xIonDense[ipHELIUM][1];
		/*Helium ++*/
		fcolldensity[5] = dense.xIonDense[ipHELIUM][2];
		/*Molecular Hydrogen*/
		fcolldensity[8] = findspecies("H2")->hevmol;
		/*Ortho(6) and Para(7) Mol Hydrogen*/
		if(h2.lgH2ON)
		{
			fcolldensity[6] = h2.ortho_density;
			fcolldensity[7] = h2.para_density;
		}
		else
		{		
			fcolldensity[6] = (fcolldensity[8])*(0.75);
			fcolldensity[7] = (fcolldensity[8])*(0.25);
		}

		/*Updating the CollRate*/
		for( long ipHi=0; ipHi<Species[ipSpecies].numLevels_local; ipHi++ )
		{	
			for( long ipLo=0; ipLo<Species[ipSpecies].numLevels_local; ipLo++ )
			{
				if( ipHi == ipLo )
				{
					CollRate[ipHi][ipLo] = 0.;
					continue;
				}

				for( intCollNo=0; intCollNo<NUM_COLLIDERS; intCollNo++ )
				{
					if(CollRatesArray[ipSpecies][intCollNo]!=NULL)
					{
						CollRate[ipHi][ipLo] += CollRatesArray[ipSpecies][intCollNo][ipHi][ipLo]*fcolldensity[intCollNo];
					}
				}
				/*Correction for Helium*/
				/*To include the effects of collision with Helium,the user must multiply the density by 1.14*/
				/*pg 5,Schoier et al 2004*/
				if( Species[ipSpecies].lgMolecular )
				{
					/*The collision rate coefficients for helium should not be present and that for molecular hydrogen should be present*/
					if(CollRatesArray[ipSpecies][3]==NULL && CollRatesArray[ipSpecies][8]!=NULL )
					{
						CollRate[ipHi][ipLo] += 0.7 *CollRatesArray[ipSpecies][8][ipHi][ipLo]*fcolldensity[3];
					}

					/* Put in something for hydrogen collisions if not in database */
					if(CollRatesArray[ipSpecies][2]==NULL )
					{
						double colliderDensity = fcolldensity[2];
						if( CollRatesArray[ipSpecies][3]!=NULL ) //He0
						{
							CollRate[ipHi][ipLo] += 2.0 *CollRatesArray[ipSpecies][3][ipHi][ipLo]*colliderDensity;
						}
						else if( CollRatesArray[ipSpecies][6]!=NULL ) //ortho-H2
						{
							CollRate[ipHi][ipLo] += 1.4 *CollRatesArray[ipSpecies][6][ipHi][ipLo]*colliderDensity;
						}
						else if( CollRatesArray[ipSpecies][7]!=NULL ) //para-H2
						{
							CollRate[ipHi][ipLo] += 1.4 *CollRatesArray[ipSpecies][7][ipHi][ipLo]*colliderDensity;
						}
						else if( CollRatesArray[ipSpecies][8]!=NULL ) // total H2
						{
							CollRate[ipHi][ipLo] += 1.4 *CollRatesArray[ipSpecies][8][ipHi][ipLo]*colliderDensity;
						}
						else
						{
							/* make something up for hydrogen collisions, but not too large */
							CollRate[ipHi][ipLo] += 1e-13 * colliderDensity;
						}
					}

					/* Put in something for proton collisions if not in database */
					if(CollRatesArray[ipSpecies][1]==NULL )
					{
						double colliderDensity = fcolldensity[1];
						if( CollRatesArray[ipSpecies][4]!=NULL ) //He+
						{
							CollRate[ipHi][ipLo] += 2.0 *CollRatesArray[ipSpecies][3][ipHi][ipLo]*colliderDensity;
						}
						else
						{
							/* make something up for proton collisions, but not too large */
							CollRate[ipHi][ipLo] += 1e-13 * colliderDensity;
						}
					}
				}

				/* now add in excitations resulting from cosmic ray secondaries */
				if( ipHi>ipLo )
				{
					fixit();// this g-bar only works for permitted lines
					// these lines are mostly forbidden - need to add branch to
					// do forbidden transitions
					// 2010 jul 17, robin points out that the scaling below evaluates
					// to unity - comment out maths and have posted ticket #160
					if( dBaseTrans[ipSpecies][ipHi][ipLo].ipCont > 0 )
					{
						/* get secondaries for all iso lines by scaling LyA 
						 * excitation by ratio of cross section (oscillator strength/energy) 
						 * Born approximation or plane-wave approximation */
						double RateUp = secondaries.x12tot /*  *
							(dBaseTrans[ipSpecies][ipHi][ipLo].Emis->gf/
							dBaseTrans[ipSpecies][ipHi][ipLo].EnergyWN) /
							(dBaseTrans[ipSpecies][ipHi][ipLo].Emis->gf/
							dBaseTrans[ipSpecies][ipHi][ipLo].EnergyWN)*/;

						double RateDown = RateUp * dBaseTrans[ipSpecies][ipHi][ipLo].Lo->g /
							dBaseTrans[ipSpecies][ipHi][ipLo].Hi->g;

						CollRate[ipLo][ipHi] += RateUp;
						CollRate[ipHi][ipLo] += RateDown;
					}
				}
			}
		}

		/* solve the n-level atom */
		atom_levelN(
			/* Species[ipSpecies].numLevels_local is the number of levels to compute*/ 
			Species[ipSpecies].numLevels_local, 
			/* ABUND is total abundance of species, used for nth equation
			 * if balance equations are homogeneous */
			abund, 
			/* g(Species[ipSpecies].numLevels_local) is statistical weight of levels */
			g, 
			/* EX(Species[ipSpecies].numLevels_local) is excitation potential of levels, deg K or wavenumbers
			 * 0 for lowest level, all are energy rel to ground NOT d(ENER) */
			ex, 
			/* this is 'K' for ex[] as Kelvin deg, is 'w' for wavenumbers */
			'w',
			/* populations [cm-3] of each level as deduced here */
			pops, 
			/* departure coefficient, derived below */
			depart,
			/* net transition rate, A * esc prob, s-1 */
			&AulEscp, 
			/* col str from up to low */
			&col_str, 
			/* AulDest[ihi][ilo] is destruction rate, trans from ihi to ilo, A * dest prob,
			 * asserts confirm that [ihi][ilo] is zero */
			&AulDest, 
			/* AulPump[ilo][ihi] is pumping rate from lower to upper level (s^-1), (hi,lo) must be zero  */
			&AulPump, 
			/* collision rates (s^-1), evaluated here and returned for cooling by calling function,
			 * unless following flag is true.  If true then calling function has already filled
			 * in these rates.  CollRate[ipSpecies][j] is rate from ipSpecies to j */
			&CollRate,
			/* this is an additional creation rate from continuum, normally zero, units cm-3 s-1 */
			source,
			/* this is an additional destruction rate to continuum, normally zero, units s-1 */
			sink,
			/* flag saying whether CollRate already done (true), or we need to do it here (false),
			 * this is stored in data)[ihi][ilo] as either downward rate or collis strength*/
			true,
			/* total cooling and its derivative, set here but nothing done with it*/
			&cooltl, 
			&coolder, 
			/* string used to identify calling program in case of error */
			spName, 
			/* nNegPop flag indicating what we have done
			 * positive if negative populations occurred
			 * zero if normal calculation done
			 * negative if too cold (for some atoms other routine will be called in this case) */
			&nNegPop,
			/* true if populations are zero, either due to zero abundance of very low temperature */
			&lgZeroPop,
			/* option to print debug information */
			lgDeBug );

		if( nNegPop > 0 )
		{
			/* negative populations occurred */
			fprintf(ioQQQ," PROBLEM in dBase_solve, atom_levelN returned negative population .\n");
			cdEXIT( EXIT_FAILURE );
		}

		// highest levels may have no population
		while( (pops[Species[ipSpecies].numLevels_local-1]<=0 ) &&
			(Species[ipSpecies].numLevels_local > 1) )
				--Species[ipSpecies].numLevels_local;

		for( long j=0;j< Species[ipSpecies].numLevels_local; j++ )
			dBaseStates[ipSpecies][j].Pop = MAX2(pops[j],SMALLFLOAT);

		for( long j=Species[ipSpecies].numLevels_local;
			j< Species[ipSpecies].numLevels_max; j++ )
				dBaseStates[ipSpecies][j].Pop = 0.;

		/*Atmol  line*/
		for(long ipHi = 1; ipHi < Species[ipSpecies].numLevels_local; ipHi++ )
		{	
			for(long ipLo = 0; ipLo < ipHi; ipLo++ )
			{
				double qLU = dBaseTrans[ipSpecies][ipHi][ipLo].Lo->Pop * CollRate[ipLo][ipHi];
				double qUL = dBaseTrans[ipSpecies][ipHi][ipLo].Hi->Pop * CollRate[ipHi][ipLo];
				double netCooling = (qLU - qUL) * dBaseTrans[ipSpecies][ipHi][ipLo].EnergyErg;

				if( netCooling > 0. )
				{
					dBaseTrans[ipSpecies][ipHi][ipLo].Coll.cool = netCooling;
					dBaseTrans[ipSpecies][ipHi][ipLo].Coll.heat = 0.;
				}
				else
				{
					dBaseTrans[ipSpecies][ipHi][ipLo].Coll.cool = 0;
					dBaseTrans[ipSpecies][ipHi][ipLo].Coll.heat = -netCooling;
				}

				Species[ipSpecies].CoolTotal += netCooling;
				totaldCdT += netCooling *
					(dBaseTrans[ipSpecies][ipHi][ipLo].EnergyK * thermal.tsq1 - thermal.halfte);
	
				if( dBaseTrans[ipSpecies][ipHi][ipLo].ipCont > 0 )
				{
					dBaseTrans[ipSpecies][ipHi][ipLo].Emis->phots =
						dBaseTrans[ipSpecies][ipHi][ipLo].Emis->Aul * 
						dBaseTrans[ipSpecies][ipHi][ipLo].Hi->Pop *
						(dBaseTrans[ipSpecies][ipHi][ipLo].Emis->Pesc + 
						dBaseTrans[ipSpecies][ipHi][ipLo].Emis->Pelec_esc);

					dBaseTrans[ipSpecies][ipHi][ipLo].Emis->xIntensity = 
						dBaseTrans[ipSpecies][ipHi][ipLo].Emis->phots * 
						dBaseTrans[ipSpecies][ipHi][ipLo].EnergyErg;

					/* population of lower level rel to ion, corrected for stim em */
					dBaseTrans[ipSpecies][ipHi][ipLo].Emis->PopOpc =
						dBaseTrans[ipSpecies][ipHi][ipLo].Lo->Pop - dBaseTrans[ipSpecies][ipHi][ipLo].Hi->Pop*
						dBaseTrans[ipSpecies][ipHi][ipLo].Lo->g/dBaseTrans[ipSpecies][ipHi][ipLo].Hi->g;

					/* it's possible for this sum to be zero.  Set ratio to zero in that case. */
					if( CollRate[ipLo][ipHi]+AulPump[ipLo][ipHi] > 0. )
					{
						dBaseTrans[ipSpecies][ipHi][ipLo].Emis->ColOvTot = CollRate[ipLo][ipHi]/
							(CollRate[ipLo][ipHi]+AulPump[ipLo][ipHi]);
					}
					else
						dBaseTrans[ipSpecies][ipHi][ipLo].Emis->ColOvTot = 0.;

					// this is only used for save line printout.  Maybe colliders may be involved, but
					// this simple approximation of a "collision strength" should be good enough for
					// the purposes of that printout.
					dBaseTrans[ipSpecies][ipHi][ipLo].Coll.col_str = (realnum)(CollRate[ipHi][ipLo] *
						(g[ipHi]*phycon.sqrte)/COLL_CONST);
				}
				else
					dBaseTrans[ipSpecies][ipHi][ipLo].Coll.col_str = 0.;
			}
		}

		fixit();
		//ASSERT( fp_equal( Species[ipSpecies].CoolTotal, cooltl ) );
		
		for(long ipHi = Species[ipSpecies].numLevels_local; ipHi< Species[ipSpecies].numLevels_max; ipHi++ )
		{	
			for(long ipLo = 0; ipLo < ipHi; ipLo++ )
			{	
				EmLineZero( dBaseTrans[ipSpecies][ipHi][ipLo].Emis );		
			}
		}

		thermal.dCooldT += totaldCdT;
		CoolAdd( Species[ipSpecies].chLabel, 0. , MAX2(0.,Species[ipSpecies].CoolTotal) );
		totalHeating += MAX2(0., -Species[ipSpecies].CoolTotal);

		/* option to print departure coefficients */
		{
			enum {DEBUG_LOC=false};

			if( DEBUG_LOC )
			{
				fprintf( ioQQQ, " Departure coefficients for species %li\n", ipSpecies );
				for( long j=0; j< Species[ipSpecies].numLevels_local; j++ )
				{
					fprintf( ioQQQ, " level %li \t Depar Coef %e\n", j, depart[j] );
				}
			}
		}
	}

	// total heating for all dBase Species	
	thermal.heating[0][27] = totalHeating;

	return;
}

/*This function fills the population of the states to a dummy value of 0*/
void states_popfill( void)
{
	DEBUG_ENTRY( "states_popfill()" );

	for( long i=0; i<nSpecies; i++)
	{
		for( long j=0; j<Species[i].numLevels_max; j++)
		{
			dBaseStates[i][j].Pop = 0.;
		}
	}
	return;
}

/*Leiden*/
STATIC double LeidenCollRate(long ipSpecies, long ipCollider, long ipHi, long ipLo, double ftemp)
{
	double ret_collrate;
	int inttemps;
	DEBUG_ENTRY("LeidenCollRate()");
	inttemps = AtmolCollRateCoeff[ipSpecies][ipCollider].ntemps;

	if( AtmolCollRateCoeff[ipSpecies][ipCollider].temps == NULL )
	{
		return 0.;
	}

	if(ftemp < AtmolCollRateCoeff[ipSpecies][ipCollider].temps[0])
	{
		ret_collrate = AtmolCollRateCoeff[ipSpecies][ipCollider].collrates[ipHi][ipLo][0];
	}
	else if(ftemp > AtmolCollRateCoeff[ipSpecies][ipCollider].temps[inttemps-1])
	{
		ret_collrate = AtmolCollRateCoeff[ipSpecies][ipCollider].collrates[ipHi][ipLo][inttemps-1];
	}
	else if( AtmolCollRateCoeff[ipSpecies][ipCollider].ntemps==1 )
	{
		// lamda data files can have only one temperture (see cn.dat as of Feb. 10, 2009)
		ret_collrate = AtmolCollRateCoeff[ipSpecies][ipCollider].collrates[ipHi][ipLo][0];
	}
	else
	{
		ret_collrate = linint(AtmolCollRateCoeff[ipSpecies][ipCollider].temps,
			AtmolCollRateCoeff[ipSpecies][ipCollider].collrates[ipHi][ipLo],
			AtmolCollRateCoeff[ipSpecies][ipCollider].ntemps,
			ftemp);
	}

	if(DEBUGSTATE)
		printf("\nThe collision rate at temperature %f is %e",ftemp,ret_collrate);

	ASSERT( !isnan( ret_collrate ) );
	return(ret_collrate);

}

/*CHIANTI*/
STATIC double CHIANTI_Upsilon(long ipSpecies, long ipCollider, long ipHi, long ipLo, double ftemp)
{
	double fdeltae,fscalingparam,fkte,fxt,fsups,fups;
	int intxs,inttype,intsplinepts;

	DEBUG_ENTRY("CHIANTI_Upsilon()");

	if( AtmolCollSplines[ipSpecies][ipHi][ipLo][ipCollider].collspline == NULL )
	{
		return 0.;
	}

	intsplinepts = AtmolCollSplines[ipSpecies][ipHi][ipLo][ipCollider].nSplinePts;
	inttype = AtmolCollSplines[ipSpecies][ipHi][ipLo][ipCollider].intTranType;
	fdeltae = AtmolCollSplines[ipSpecies][ipHi][ipLo][ipCollider].EnergyDiff;
	fscalingparam = AtmolCollSplines[ipSpecies][ipHi][ipLo][ipCollider].ScalingParam;

	fkte = ftemp/fdeltae/1.57888e5;

	/*Way the temperature is scaled*/
	/*Burgess&Tully 1992:Paper gives only types 1 to 4*/
	/*Found that the expressions were the same for 5 & 6 from the associated routine DESCALE_ALL*/
	/*What about 7,8&9?*/
	if( inttype ==1 || inttype==4 )
	{
		fxt = 1-(log(fscalingparam)/(log(fkte+fscalingparam)));
	}
	else if(inttype  == 2 || inttype == 3||inttype == 5 || inttype == 6)
	{
		fxt = fkte/(fkte+fscalingparam);
	}
	else
		TotalInsanity();

	double xs[9],*spl,*spl2;
	/*Creating spline points array*/
	if(intsplinepts == 5)
	{
		spl = AtmolCollSplines[ipSpecies][ipHi][ipLo][ipCollider].collspline;
		for(intxs=0;intxs<5;intxs++)
		{
			xs[intxs] = 0.25*intxs;
			if(DEBUGSTATE)
			{
				printf("The xs and spl values are %f and %f \n",xs[intxs],spl[intxs]);
				getchar();
			}
		}
	}
	else if(intsplinepts == 9)
	{
		spl = AtmolCollSplines[ipSpecies][ipHi][ipLo][ipCollider].collspline;	
		for( intxs=0; intxs<9; intxs++ )
		{
			xs[intxs] = 0.125*intxs;
			if(DEBUGSTATE)
			{
				printf("The xs and spl values are %f and %f \n",xs[intxs],spl[intxs]);
				getchar();
			}
		}
	}
	else
	{
		TotalInsanity();
	}

	/*Finding the second derivative*/
	spl2 = AtmolCollSplines[ipSpecies][ipHi][ipLo][ipCollider].SplineSecDer;

	if(DEBUGSTATE)
	{
		printf("\n");
		for(intxs=0;intxs<intsplinepts;intxs++)
		{
			printf("The %d value of 2nd derivative is %f \n",intxs+1,spl2[intxs]);
		}
	}

	/*Extracting out the value*/
	//splint(xs,spl,spl2,intsplinepts,fxt,&fsups);

	fsups = linint( xs, spl, intsplinepts, fxt);

	/*Finding upsilon*/
	if(inttype == 1)
	{
		fups = fsups*log(fkte+exp(1.0));
	}
	else if(inttype == 2)
	{
		fups = fsups;
	}
	else if(inttype == 3)
	{
		fups = fsups/(fkte+1.0) ;
	}
	else if(inttype == 4)
	{
		fups = fsups*log(fkte+fscalingparam) ;
	}
	else if(inttype == 5)
	{
		fups = fsups/fkte ;
	}
	else if(inttype == 6)
	{
		fups = pow(10.0,fsups) ;
	}
	else
	{
		TotalInsanity();
	}

	if( fups < 0. ) 
	{
		fprintf( ioQQQ," WARNING: Negative upsilon in species %s, collider %li, indices %4li %4li, Te = %e.\n",
				Species[ipSpecies].chLabel, ipCollider, ipHi, ipLo, ftemp );
		fups = 0.;
	}
	ASSERT(fups>=0);
	return(fups);
}

