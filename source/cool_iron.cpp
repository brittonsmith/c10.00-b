/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolIron compute iron cooling */
/*Fe_10_11_13_cs evaluate collision strength for Fe */
/*fe14cs compute collision strengths for forbidden transitions */
/*Fe3Lev14 compute populations and cooling due to 14 level Fe III ion */
/*Fe4Lev12 compute populations and cooling due to 12 level Fe IV ion */
/*Fe7Lev8 compute populations and cooling due to 8 level Fe VII ion */
/*Fe2_cooling compute cooling due to FeII emission */
/*Fe11Lev5 compute populations and cooling due to 5 level Fe 11 ion */
/*Fe13Lev5 compute populations and cooling due to 5 level Fe 13 ion */
#include "cddefines.h"
#include "physconst.h"
#include "dense.h"
#include "coolheavy.h"
#include "taulines.h"
#include "phycon.h"
#include "iso.h"
#include "conv.h"
#include "trace.h"
#include "hydrogenic.h"
#include "ligbar.h"
#include "cooling.h"
#include "thermal.h"
#include "lines_service.h"
#include "atoms.h"
#include "atomfeii.h"
#include "fe.h"

/*Fe11Lev5 compute populations and cooling due to 5 level Fe 11 ion */
STATIC void Fe11Lev5(void);

/*Fe13Lev5 compute populations and cooling due to 5 level Fe 13 ion */
STATIC void Fe13Lev5(void);

/*fe14cs compute collision strengths for forbidden transitions */
STATIC void fe14cs(double te1, 
	  double *csfe14);

/*Fe7Lev8 compute populations and cooling due to 8 level Fe VII ion */
STATIC void Fe7Lev8(void);

/*Fe3Lev14 compute populations and cooling due to 14 level Fe III ion */
STATIC void Fe3Lev14(void);

/*Fe4Lev12 compute populations and cooling due to 12 level Fe IV ion */
STATIC void Fe4Lev12(void);

/*Fe_10_11_13_cs evaluate collision strength for Fe */
STATIC double Fe_10_11_13_cs(
	/* ion, one of 10, 11, or 13 on physics scale */
	int ion, 
	/* initial and final index of levels, with lowest energy 1, highest of 5 
	 * on physics scale */
	int init, 
	int final )
{
	const int N = 10;
	static double Fe10cs[6][6][2];
	static double Fe11cs[6][6][2];
	static double Fe13cs[6][6][2];
	int i, j;
	double cs;
	int index = 0;
	double temp_max, temp_min = 4;
	double temp_log = phycon.alogte;
	static bool lgFirstTime = true;

	DEBUG_ENTRY( "Fe_10_11_13_cs()" );

	if( lgFirstTime )
	{
		/* fit coefficients for collision strengths */
		double aFe10[N] = {10.859,-1.1541,11.593,22.333,-0.4283,7.5663,3.087,1.0937,0.8261,59.678};
		double bFe10[N] = {-1.4804,0.4956,-2.1096,-4.1288,0.1929,-1.3525,-0.5531,-0.1748,-0.1286,-11.081};
		double aFe11[N] = {5.7269,1.2885,4.0877,0.4571,1.2911,2.2339,0.3621,0.7972,0.2225,1.1021};
		double bFe11[N] = {-0.7559,-0.1671,-0.5678,-0.0653,-0.1589,-0.2924,-0.0506,-0.1038,-0.0302,-0.1062};
		double aFe13[N] = {2.9102,1.8682,-0.353,0.0622,14.229,-4.3845,0.0375,-6.9222,0.688,-0.0609};
		double bFe13[N] = {-0.4158,-0.242,0.1417,0.0023,-2.0643,1.2573,0.0286,2.0919,-0.083,0.1487};
		/* do one-time initialization 
		 * assigning array initially to zero */ 
		for( i=0; i < 6; i++ )
		{
			for( j=0; j < 6; j++ )
			{
				set_NaN( Fe10cs[i][j], 2L );
				set_NaN( Fe11cs[i][j], 2L );
				set_NaN( Fe13cs[i][j], 2L );
			}
		}

		/* reading coefficients into 3D array */
		for( i=1; i < 6; i++ )
		{
			for( j=i+1; j < 6; j++ )
			{
				Fe10cs[i][j][0] = aFe10[index];
				Fe10cs[i][j][1] = bFe10[index];
				Fe11cs[i][j][0] = aFe11[index];
				Fe11cs[i][j][1] = bFe11[index];
				Fe13cs[i][j][0] = aFe13[index];
				Fe13cs[i][j][1] = bFe13[index];
				index++;
			}
		}
		lgFirstTime = false;
	}

	/* Invalid entries returns '-1':the initial indices are smaller than 
	 * the final indices */
	if( init >= final )
	{
		cs = -1;
	}
	/* Invalid returns '-1': the indices are greater than 5 or smaller than 0 */
	else if( init < 1 || init > 4 || final < 2 || final > 5 )
	{
		cs = -1;
	}
	else
	{
		/* cs = a + b*log(T) 
		 * if temp is out of range, return boundary values */
		if( ion == 10 )
		{
			temp_max = 5;
			temp_log = MAX2(temp_log, temp_min);
			temp_log = MIN2(temp_log, temp_max);
			cs = Fe10cs[init][final][0] + Fe10cs[init][final][1]*temp_log;
		}
		else if( ion == 11 )
		{
			temp_max = 6.7;
			temp_log = MAX2(temp_log, temp_min);
			temp_log = MIN2(temp_log, temp_max);
			cs = Fe11cs[init][final][0] + Fe11cs[init][final][1]*temp_log;
		}
		else if( ion ==13 )
		{
			temp_max = 5;
			temp_log = MAX2(temp_log, temp_min);
			temp_log = MIN2(temp_log, temp_max);
			cs = Fe13cs[init][final][0] + Fe13cs[init][final][1]*temp_log;
		}
		else
			/* this can't happen */
			TotalInsanity();
	}

	return cs;
}

/*Fe2_cooling compute cooling due to FeII emission */
STATIC void Fe2_cooling( void )
{
	const int NLFE2 = 6;
	long int i, j;
	int nNegPop;

	static double **AulPump,
		**CollRate,
		**AulEscp,
		**col_str ,
		**AulDest, 
		*depart,
		*pops,
		*destroy,
		*create;

	static bool lgFirst=true;
	bool lgZeroPop;

	/* stat weights for Fred's 6 level model FeII atom */
	static double gFe2[NLFE2]={1.,1.,1.,1.,1.,1.};
	/* excitation energies (Kelvin) for Fred's 6 level model FeII atom */
	static double ex[NLFE2]={0.,3.32e4,5.68e4,6.95e4,1.15e5,1.31e5};

	/* these are used to only evaluated FeII one time per temperature, zone, 
	 * and abundance */
	static double TUsed = 0.; 
	static realnum AbunUsed = 0.;
	/* remember which zone last evaluated with */
	static long int nZUsed=-1,
		/* make sure at least two calls per zone */
		nCall=0;

	DEBUG_ENTRY( "Fe2_cooling()" );

	/* return if nothing to do */
	if( dense.xIonDense[ipIRON][1] == 0. )
	{
		/* zero abundance, do nothing */
		/* heating cooling and derivative from large atom */
		FeII.Fe2_large_cool = 0.;
		FeII.Fe2_large_heat = 0.;
		FeII.ddT_Fe2_large_cool = 0.;

		/* cooling and derivative from simple UV atom */
		FeII.Fe2_UVsimp_cool = 0.;
		FeII.ddT_Fe2_UVsimp_cool = 0.;

		/* now zero out intensities of all FeII lines and level populations */
		FeIIIntenZero();
		return;
	}

	/* evaluate FeII model atom, by default only the lowest 16 levels
	 * but number of levels can be increased with atom feii command */

	/* totally reevaluate FeII atom if new zone, or cooling is significant
	 * and temperature changed, we are in search phase,
	 * lgSlow option set true with atom FeII slow, forces constant
	 * evaluation of atom */
	if( FeII.lgSlow || 
		conv.lgFirstSweepThisZone || conv.lgLastSweepThisZone ||
	    /* check whether things have changed on later calls */
	    ( !fp_equal( phycon.te, TUsed ) && fabs(FeII.Fe2_large_cool/thermal.ctot) > 0.002 &&  
	      fabs(dense.xIonDense[ipIRON][1]-AbunUsed)/SDIV(AbunUsed) > 0.002 ) ||
	    ( !fp_equal( phycon.te, TUsed ) && fabs(FeII.Fe2_large_cool/thermal.ctot) > 0.01) )
	{

		if( conv.lgFirstSweepThisZone )
			/* first call this zone set nCall to zero*/
			nCall = 0;
		else
			/* not first call, increment, check above to make sure at least
			 * two evaluations */
			++nCall;

		/* option to trace convergence and FeII calls */
		if( trace.nTrConvg >= 5 )
		{
			fprintf( ioQQQ, "        CoolIron5 calling FeIILevelPops since ");
			if( ! fp_equal( phycon.te, TUsed ) )
			{
				fprintf( ioQQQ, 
					"temperature changed, old new are %g %g, nCall %li ", 
					TUsed, phycon.te , nCall);
			}
			else if( nzone != nZUsed )
			{
				fprintf( ioQQQ, 
					"new zone, nCall %li ", nCall );
			}
			else if( FeII.lgSlow )
			{
				fprintf( ioQQQ, 
					"FeII.lgSlow set  %li", nCall );
			}
			else if( conv.lgSearch )
			{
				fprintf( ioQQQ, 
					" in search phase  %li", nCall );
			}
			else if( nCall < 2 )
			{
				fprintf( ioQQQ, 
					"not second nCall %li " , nCall );
			}
			else if( ! fp_equal( phycon.te, TUsed ) && FeII.Fe2_large_cool/thermal.ctot > 0.001 )
			{
				fprintf( ioQQQ, 
					"temp or cooling changed, new are %g %g nCall %li ", 
					phycon.te, FeII.Fe2_large_cool, nCall );
			}
			else
			{
				fprintf(ioQQQ, "????");
			}
			fprintf(ioQQQ, "\n");
		}

		/* remember parameters for current conditions */
		TUsed = phycon.te;
		AbunUsed = dense.xIonDense[ipIRON][1];
		nZUsed = nzone;

		/* this print turned on with atom FeII print command */
		if( FeII.lgPrint )
		{
			fprintf(ioQQQ,
				" FeIILevelPops called zone %4li te %5f abun %10e c(fe/tot):%6f nCall %li\n", 
				nzone,phycon.te,AbunUsed,FeII.Fe2_large_cool/thermal.ctot,nCall);
		}

		/* this solves the multi-level problem, 
		 * sets FeII.Fe2_large_cool, 
				FeII.Fe2_large_heat, & 
				FeII.ddT_Fe2_large_cool 
				but does nothing with them */

		// line RT update, followed by level population solver
		FeII_RT_Make();
		FeIILevelPops();
		{
			enum{DEBUG_LOC=false};
			if( DEBUG_LOC && iteration > 1 && nzone >=4 )
			{
				fprintf(ioQQQ,"DEBUG1\t%li\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n", 
					nzone,
					phycon.te,
					dense.gas_phase[ipHYDROGEN],
					dense.eden,
					FeII.Fe2_large_cool , 
					FeII.ddT_Fe2_large_cool ,
					FeII.Fe2_large_cool/dense.eden/dense.gas_phase[ipHYDROGEN] , 
					thermal.ctot );
			}
		}

		if( trace.nTrConvg >= 5 || FeII.lgPrint)
		{
			/* spacing needed to get proper trace convergence output */
			fprintf( ioQQQ, "           FeIILevelPops5 returned cool=%.2e heat=%.2e derivative=%.2e\n",
					FeII.Fe2_large_cool,FeII.Fe2_large_heat ,FeII.ddT_Fe2_large_cool);
		}

	}
	else if( ! fp_equal( dense.xIonDense[ipIRON][1], AbunUsed ) )
	{
		realnum ratio;
		/* this branch, same zone and temperature, but small change in abundance, so just
		 * rescale cooling and derivative by this change.  assumption is that very small changes
		 * in abundance occurs as ots rates damp out */
		if( trace.nTrConvg >= 5 )
		{
			fprintf( ioQQQ, 
				"       CoolIron rescaling FeIILevelPops since small change, CFe2=%.2e CTOT=%.2e\n",
				FeII.Fe2_large_cool,thermal.ctot);
		}
		ratio = dense.xIonDense[ipIRON][1]/AbunUsed;
		FeII.Fe2_large_cool *= ratio;
		FeII.ddT_Fe2_large_cool *= ratio;
		FeII.Fe2_large_heat *= ratio;
		AbunUsed = dense.xIonDense[ipIRON][1];
	}
	else
	{
		/* this is case where temp is unchanged, so heating and cooling same too */
		if( trace.nTrConvg >= 5 )
		{
			fprintf( ioQQQ, "       CoolIron NOT calling FeIILevelPops\n");
		}
	}

	/* now update total cooling and its derivative 
	 * all paths flow through here */
	/* FeII.Fecool = FeII.Fe2_large_cool; */
	CoolAdd("Fe 2",0,MAX2(0.,FeII.Fe2_large_cool));

	/* add negative cooling to heating stack */
	thermal.heating[25][27] = MAX2(0.,FeII.Fe2_large_heat);

	/* counts as heating derivative if negative cooling */
	if( FeII.Fe2_large_cool > 0. )
	{
		/* >>chng 01 mar 16, add factor of 3 due to conv problems after changing damper */
		thermal.dCooldT += 3.*FeII.ddT_Fe2_large_cool;
	}

	if( trace.lgTrace && trace.lgCoolTr )
	{
		fprintf( ioQQQ, " Large FeII returns te, cooling, dc=%11.3e%11.3e%11.3e\n", 
		  phycon.te, FeII.Fe2_large_cool, FeII.ddT_Fe2_large_cool );
	}

	/* >>chng 05 nov 29, still do simple UV atom if only ground term is done */
	if( !FeII.lgFeIILargeOn )
	{

		/* following treatment of Fe II follows
		 * >>refer	fe2	model	Wills, B.J., Wills, D., Netzer, H. 1985, ApJ, 288, 143
		 * all elements are used, and must be set to zero if zero */

		/* set up space for simple  model of UV FeII emission */
		if( lgFirst )
		{
			/* will never do this again in this core load */
			lgFirst = false;
			/* allocate the 1D arrays*/
			pops = (double *)MALLOC( sizeof(double)*(NLFE2) );
			create = (double *)MALLOC( sizeof(double)*(NLFE2) );
			destroy = (double *)MALLOC( sizeof(double)*(NLFE2) );
			depart = (double *)MALLOC( sizeof(double)*(NLFE2) );
			/* create space for the 2D arrays */
			AulPump = ((double **)MALLOC((NLFE2)*sizeof(double *)));
			CollRate = ((double **)MALLOC((NLFE2)*sizeof(double *)));
			AulDest = ((double **)MALLOC((NLFE2)*sizeof(double *)));
			AulEscp = ((double **)MALLOC((NLFE2)*sizeof(double *)));
			col_str = ((double **)MALLOC((NLFE2)*sizeof(double *)));

			for( i=0; i < NLFE2; ++i )
			{
				AulPump[i] = ((double *)MALLOC((NLFE2)*sizeof(double )));
				CollRate[i] = ((double *)MALLOC((NLFE2)*sizeof(double )));
				AulDest[i] = ((double *)MALLOC((NLFE2)*sizeof(double )));
				AulEscp[i] = ((double *)MALLOC((NLFE2)*sizeof(double )));
				col_str[i] = ((double *)MALLOC((NLFE2)*sizeof(double )));
			}
		}

		/*zero out all arrays, then check that upper diagonal remains zero below */
		for( i=0; i < NLFE2; i++ )
		{
			create[i] = 0.;
			destroy[i] = 0.;
			for( j=0; j < NLFE2; j++ )
			{
				/*data[j][i] = 0.;*/
				col_str[j][i] = 0.;
				AulEscp[j][i] = 0.;
				AulDest[j][i] = 0.;
				AulPump[j][i] = 0.;
			}
		}

		/* now put in real data for lines */
		AulEscp[1][0] = 1.;
		AulEscp[2][0] = ( TauLines[ipTuv3].Emis->Pesc + TauLines[ipTuv3].Emis->Pelec_esc)*TauLines[ipTuv3].Emis->Aul;
		AulDest[2][0] = TauLines[ipTuv3].Emis->Pdest*TauLines[ipTuv3].Emis->Aul;
		AulPump[0][2] = TauLines[ipTuv3].Emis->pump;

		AulEscp[5][0] = (TauLines[ipTFe16].Emis->Pesc + TauLines[ipTFe16].Emis->Pelec_esc)*TauLines[ipTFe16].Emis->Aul;
		AulDest[5][0] = TauLines[ipTFe16].Emis->Pdest*TauLines[ipTFe16].Emis->Aul;
		/* continuum pumping of n=6 */
		AulPump[0][5] = TauLines[ipTFe16].Emis->pump;
		/* Ly-alpha pumping */

		double PumpLyaFeII = StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH2p].Pop*
			Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->Aul*
			hydro.dstfe2lya/SDIV(dense.xIonDense[ipIRON][1]);

		AulPump[0][5] += PumpLyaFeII;

		AulEscp[2][1] = (TauLines[ipTr48].Emis->Pesc + TauLines[ipTr48].Emis->Pelec_esc)*TauLines[ipTr48].Emis->Aul;
		AulDest[2][1] = TauLines[ipTr48].Emis->Pdest*TauLines[ipTr48].Emis->Aul;
		AulPump[1][2] = TauLines[ipTr48].Emis->pump;

		AulEscp[5][1] = (TauLines[ipTFe26].Emis->Pesc + TauLines[ipTFe26].Emis->Pelec_esc)*TauLines[ipTFe26].Emis->Aul;
		AulDest[5][1] = TauLines[ipTFe26].Emis->Pdest*TauLines[ipTFe26].Emis->Aul;
		AulPump[1][5] = TauLines[ipTFe26].Emis->pump;

		AulEscp[3][2] = (TauLines[ipTFe34].Emis->Pesc + TauLines[ipTFe34].Emis->Pelec_esc)*TauLines[ipTFe34].Emis->Aul;
		AulDest[3][2] = TauLines[ipTFe34].Emis->Pdest*TauLines[ipTFe34].Emis->Aul;
		AulPump[2][3] = TauLines[ipTFe34].Emis->pump;

		AulEscp[4][2] = (TauLines[ipTFe35].Emis->Pesc + TauLines[ipTFe35].Emis->Pelec_esc)*TauLines[ipTFe35].Emis->Aul;
		AulDest[4][2] = TauLines[ipTFe35].Emis->Pdest*TauLines[ipTFe35].Emis->Aul;
		AulPump[2][4] = TauLines[ipTFe35].Emis->pump;

		AulEscp[5][3] = (TauLines[ipTFe46].Emis->Pesc + TauLines[ipTFe46].Emis->Pelec_esc)*TauLines[ipTFe46].Emis->Aul;
		AulDest[5][3] = TauLines[ipTFe46].Emis->Pdest*TauLines[ipTFe46].Emis->Aul;
		AulPump[3][5] = TauLines[ipTFe46].Emis->pump;

		AulEscp[5][4] = (TauLines[ipTFe56].Emis->Pesc + TauLines[ipTFe56].Emis->Pelec_esc)*TauLines[ipTFe56].Emis->Aul;
		AulDest[5][4] = TauLines[ipTFe56].Emis->Pdest*TauLines[ipTFe56].Emis->Aul;
		AulPump[4][5] = TauLines[ipTFe56].Emis->pump;

		/* these are collision strengths */
		col_str[1][0] = 1.;
		col_str[2][0] = 12.;
		col_str[3][0] = 1.;
		col_str[4][0] = 1.;
		col_str[5][0] = 12.;
		col_str[2][1] = 6.;
		col_str[3][1] = 1.;
		col_str[4][1] = 1.;
		col_str[5][1] = 12.;
		col_str[3][2] = 6.;
		col_str[4][2] = 12.;
		col_str[5][2] = 1.;
		col_str[4][3] = 1.;
		col_str[5][3] = 12.;
		col_str[5][4] = 6.;

		/*void atom_levelN(long,long,realnum,double[],double[],double[],double*,
		double*,double*,long*,realnum*,realnum*,STRING,int*);*/
		atom_levelN(NLFE2,
			dense.xIonDense[ipIRON][1],
			gFe2,
			ex,
			'K',
			pops,
			depart,
			&AulEscp ,
			&col_str,
			&AulDest,
			&AulPump,
			&CollRate,
			create,
			destroy,
			false,/* say atom_levelN should evaluate coll rates from cs */
			/*&&ipdest,*/
			&FeII.Fe2_UVsimp_cool,
			&FeII.ddT_Fe2_UVsimp_cool,
			"FeII",
			&nNegPop,
			&lgZeroPop,
			false );

		/* nNegPop positive if negative pops occurred, negative if too cold */
		if( nNegPop > 0 /*negative if too cold - that is not negative and is OK */ )
		{
			fprintf(ioQQQ," PROBLEM, atom_levelN returned negative population for simple UV FeII.\n");
		}

		/* add heating - cooling by this process */;
		CoolAdd("Fe 2",0,MAX2(0.,FeII.Fe2_UVsimp_cool));
		thermal.heating[25][27] = MAX2(0.,-FeII.Fe2_UVsimp_cool);
		thermal.dCooldT += FeII.ddT_Fe2_UVsimp_cool;

		/* LIMLEVELN is the dim of the PopLevels vector */
		ASSERT( NLFE2 <= LIMLEVELN );
		for( i=0; i < NLFE2; ++i )
		{
			atoms.PopLevels[i] = pops[i];
			atoms.DepLTELevels[i] = depart[i];
		}

		TauLines[ipTuv3].Lo->Pop = pops[0];
		TauLines[ipTuv3].Hi->Pop = pops[2];
		TauLines[ipTuv3].Emis->PopOpc = (pops[0] - pops[2]);
		TauLines[ipTuv3].Emis->phots = pops[2]*AulEscp[2][0];
		TauLines[ipTuv3].Emis->xIntensity = 
			TauLines[ipTuv3].Emis->phots*TauLines[ipTuv3].EnergyErg;

		TauLines[ipTr48].Lo->Pop = pops[1];
		TauLines[ipTr48].Hi->Pop = pops[2];
		TauLines[ipTr48].Emis->PopOpc = (pops[1] - pops[2]);
		TauLines[ipTr48].Emis->phots = pops[2]*AulEscp[2][1];
		TauLines[ipTr48].Emis->xIntensity = 
			TauLines[ipTr48].Emis->phots*TauLines[ipTr48].EnergyErg;

		FeII.for7 = pops[1]*AulEscp[1][0]*4.65e-12;

		TauLines[ipTFe16].Lo->Pop = pops[0];
		TauLines[ipTFe16].Hi->Pop = pops[5];
		/* Lyman alpha optical depths are not known on first iteration,
		 * inward optical depths used, so line trapping overestimated,
		 * this can cause artificial maser in FeII - prevent by not
		 * including stimulated emission correction on first iteration */
		TauLines[ipTFe16].Emis->PopOpc = (pops[0] - pops[5]);
		TauLines[ipTFe16].Emis->phots = pops[5]*AulEscp[5][0];
		TauLines[ipTFe16].Emis->xIntensity = 
			TauLines[ipTFe16].Emis->phots*TauLines[ipTFe16].EnergyErg;

		TauLines[ipTFe26].Lo->Pop = pops[1];
		TauLines[ipTFe26].Hi->Pop = pops[5];
		TauLines[ipTFe26].Emis->PopOpc = (pops[1] - pops[5]);
		TauLines[ipTFe26].Emis->phots = pops[5]*AulEscp[5][1];
		TauLines[ipTFe26].Emis->xIntensity = 
			TauLines[ipTFe26].Emis->phots*TauLines[ipTFe26].EnergyErg;

		TauLines[ipTFe34].Lo->Pop = pops[2];
		TauLines[ipTFe34].Hi->Pop = pops[3];
		TauLines[ipTFe34].Emis->PopOpc = (pops[2] - pops[3]);
		TauLines[ipTFe34].Emis->phots = pops[3]*AulEscp[3][2];
		TauLines[ipTFe34].Emis->xIntensity = 
			TauLines[ipTFe34].Emis->phots*TauLines[ipTFe34].EnergyErg;

		TauLines[ipTFe35].Lo->Pop = pops[2];
		TauLines[ipTFe35].Hi->Pop = pops[4];
		TauLines[ipTFe35].Emis->PopOpc = (pops[2] - pops[4]);
		TauLines[ipTFe35].Emis->phots = pops[4]*AulEscp[4][2];
		TauLines[ipTFe35].Emis->xIntensity = 
			TauLines[ipTFe35].Emis->phots*TauLines[ipTFe35].EnergyErg;

		TauLines[ipTFe46].Lo->Pop = pops[3];
		TauLines[ipTFe46].Hi->Pop = pops[5];
		TauLines[ipTFe46].Emis->PopOpc = (pops[3] - pops[5]);
		TauLines[ipTFe46].Emis->phots = pops[5]*AulEscp[5][3];
		TauLines[ipTFe46].Emis->xIntensity = 
			TauLines[ipTFe46].Emis->phots*TauLines[ipTFe46].EnergyErg;

		TauLines[ipTFe56].Lo->Pop = pops[4];
		TauLines[ipTFe56].Hi->Pop = pops[5];
		TauLines[ipTFe56].Emis->PopOpc = (pops[4] - pops[5]);
		TauLines[ipTFe56].Emis->phots = pops[5]*AulEscp[5][4];
		TauLines[ipTFe56].Emis->xIntensity = 
			TauLines[ipTFe56].Emis->phots*TauLines[ipTFe56].EnergyErg;

		/* Jack's funny FeII lines, data from 
		 * >>refer	fe2	energy	Johansson, S., Brage, T., Leckrone, D.S., Nave, G. &
		 * >>refercon Wahlgren, G.M. 1995, ApJ 446, 361 */
		PutCS(10.,&TauLines[ipT191]);
		atom_level2(&TauLines[ipT191]);
	}

	{
		/*@-redef@*/
		enum{DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC && iteration > 1 && nzone >=4 )
		{
			fprintf(ioQQQ,"DEBUG2\t%.2e\t%.2e\t%.2e\n",
				phycon.te,
				FeII.Fe2_large_cool , 
				FeII.Fe2_UVsimp_cool );
		}
	}

	return;

}

/*CoolIron - calculation total cooling due to Fe */
void CoolIron(void)
{
	long int i;
	double cs ,
	  cs12, cs13, cs23,
	  cs2s2p, 
	  cs2s3p;
	realnum p2, 
		rate;

	DEBUG_ENTRY( "CoolIron()" );

	/* cooling by FeI 24m, 34.2 m */
	/* >>chng 03 nov 15, add these lines */
	/** \todo	2	- ground term is actually a fix level system, the vectors are
	 * created, with pointers ipFe1_54m , ipFe1_111m, must add collision date, use
	 * larger model atom */
	/*>>refer	Fe1	cs	Hollenbach & McKee 89 */
	/* the 24.0 micron line */
	rate = (realnum)(1.2e-7 * dense.eden + 
		/* >>chng 05 jul 05, eden to cdsqte */
		/*8.0e-10*pow((phycon.te/100.), 0.17 )*dense.xIonDense[ipHYDROGEN][0]) / dense.eden);*/
		8.0e-10*pow((phycon.te/100.), 0.17 )*dense.xIonDense[ipHYDROGEN][0]);
	LineConvRate2CS( &TauLines[ipFe1_24m]  , rate );

	rate = (realnum)(9.3e-8 * dense.eden + 
		/* >>chng 05 jul 05, eden to cdsqte */
		/*5.3e-10*pow((phycon.te/100.), 0.17 )*dense.xIonDense[ipHYDROGEN][0]) / dense.eden);*/
		5.3e-10*pow((phycon.te/100.), 0.17 )*dense.xIonDense[ipHYDROGEN][0]);
	LineConvRate2CS( &TauLines[ipFe1_35m]  , rate );

	rate = (realnum)(1.2e-7 * dense.eden + 
		/* >>chng 05 jul 05, eden to cdsqte */
		/*6.9e-10*pow((phycon.te/100.), 0.17 )*dense.xIonDense[ipHYDROGEN][0]) / dense.eden);*/
		6.9e-10*pow((phycon.te/100.), 0.17 )*dense.xIonDense[ipHYDROGEN][0]);
	TauDummy.Hi->g = TauLines[ipFe1_35m].Hi->g;
	LineConvRate2CS( &TauDummy  , rate );
	/* this says that line is a dummy, not real one */
	TauDummy.Hi->g = 0.;

	atom_level3(&TauLines[ipFe1_24m],&TauLines[ipFe1_35m],&TauDummy);

	/* series of FeI lines from Dima Verner's list, each 2-lev atom
	 *
	 * Fe I 3884 */
	MakeCS(&TauLines[ipFeI3884]);
	atom_level2(&TauLines[ipFeI3884]);

	/* Fe I 3729 */
	MakeCS(&TauLines[ipFeI3729]);
	atom_level2(&TauLines[ipFeI3729]);

	/* Fe I 3457 */
	MakeCS(&TauLines[ipFeI3457]);
	atom_level2(&TauLines[ipFeI3457]);

	/* Fe I 3021 */
	MakeCS(&TauLines[ipFeI3021]);
	atom_level2(&TauLines[ipFeI3021]);

	/* Fe I 2966 */
	MakeCS(&TauLines[ipFeI2966]);
	atom_level2(&TauLines[ipFeI2966]);

	/* >>chng 05 dec 03, move Fe2 FeII Fe II cooling into separate routine */
	Fe2_cooling();

	/* lump 3p and 3f together; cs=
	 * >>refer	fe3	as	Garstang, R.H., Robb, W.D., Rountree, S.P. 1978, ApJ, 222, 384
	 * A from
	 * >>refer	fe3	as	Garstang, R.H., 1957, Vistas in Astronomy, 1, 268
	 * FE III 5270, is 20.9% of total 
	 * >>chng 05 feb 18, Kevin Blagrave email
	 * average wavelength is 4823 with statistical weight averaging of upper energy level,
	 * as per , change 5th number from 2.948 to 2.984, also photon energy
	 * from 3.78 to 4.12 */

	/* >>chng 05 dec 16, FeIII update by Kevin Blagrave */
	/* FeIII 1122 entire multiplet - atomic data=A's Dima, CS = guess */
	PutCS(25.,&TauLines[ipT1122]);
	atom_level2(&TauLines[ipT1122]);

	/* call 14 level atom */
	Fe3Lev14();

	/* call 12 level atom */
	Fe4Lev12();

	/* FE V 3892 + 3839, data from Shields */
	CoolHeavy.c3892 = atom_pop2(7.4,25.,5.,0.6,3.7e4,dense.xIonDense[ipIRON][4])*
	  5.11e-12;
	CoolAdd("Fe 5",3892,CoolHeavy.c3892);

	/* FE VI 5177+5146+4972+4967
	 * data from 
	 * >>refer	fe6	as	Garstang, R.H., Robb, W.D., Rountree, S.P. 1978, ApJ, 222, 384 */
	CoolHeavy.c5177 = atom_pop2(1.9,28.,18.,0.52,2.78e4,dense.xIonDense[ipIRON][5])*
	  3.84e-12;
	CoolAdd("Fe 6",5177,CoolHeavy.c5177);

	/* >>chng 04 nov 04, add multi-level fe7 */
	Fe7Lev8();

	/* Fe IX 245,242
	 * all atomic data from 
	 * >>refer	fe9	all	Flower, D.R. 1976, A&A, 56, 451 */
	/* >>chng 01 sep 09, AtomSeqBeryllium will reset this to 1/3 so critical density correct */
	PutCS(0.123,&TauLines[ipT245]);
	/* AtomSeqBeryllium(cs23,cs24,cs34,tarray,a41)
	 * C245 = AtomSeqBeryllium( .087 ,.038 ,.188 , t245 ,71. )  * 8.14E-11 */
	AtomSeqBeryllium(.087,.038,.188,&TauLines[ipT245],71.);

	CoolHeavy.c242 = atoms.PopLevels[3]*8.22e-11*71.;

	/* Fe X Fe 10 Fe10 6374, A from 
	 * >>referold	fe10	as	Mason, H. 1975, MNRAS 170, 651
	 * >>referold	fe10	cs	Mohan, M., Hibbert, A., & Kingston, A.E. 1994, ApJ, 434, 389
	 * >>referold	fe10	as	Pelan, J., & Berrington, K.A. 1995, A&A Suppl, 110, 209
	 * does not agree with Mohan et al. by factor of 4
	 * 20x larger than Mason numbers used pre 1994 */
	/*cs = 13.3-2.*MIN2(7.0,phycon.alogte); */
	/*cs = MAX2(0.1,cs); */
	/** \todo	2	following to stop insane FeX strengths
	 * >>chng 96 jul 11, removed 1 / 10 factor, so using real data, 90.01
	 * cs = cs * 0.1
	 * >>chng 96 jun 03, transferred following
	 * >>chng 97 jan 31, I give up on this mess, use cs of unity */
	/*cs = 1.;*/
	/* >>chng 00 dec 05, update Fe10 collisions strength to Tayal 2000 */
	/* >>refer	fe10	cs	Tayal, S.S., 2000, ApJ, 544, 575-580 */
	/* >>chng 04 nov 10, 172.9 was mult rather than div by temp law,
	 * had no effect due to min on cs that lie below it */
	/*cs = 172.9 /( phycon.te30 * phycon.te05 * phycon.te02 *  phycon.te005 );*/
	/* tabulated cs ends at 30,000K, do not allow cs to grow larger than largest
	 * tabulated value */
	/*cs = MIN2( cs, 3.251 );*/
	/* >>chng 05 dec 19, update As to 
	 * >>refer	fe10	As	Aggarwal, K.M. & Keenan, F.P. 2004, A&A, 427, 763 
	 * value changed from 54.9 to 68.9 */
	/* >>chng 05 dec 19, update cs to
	 *>>refer	fe10	cs	Aggarwal, K.M. & Keenan, F.P. 2005, A&A, 439, 1215 */
	cs = Fe_10_11_13_cs(
	/* ion, one of 10, 11, or 13 on physics scale */
		10, 
	/* initial and final index of levels, with lowest energy 1, highest of 5 
	 * on physics scale */
		1, 
		2 );

	PutCS(cs,&TauLines[ipFe106375]);
	atom_level2(&TauLines[ipFe106375]);
	/* c6374 = atom_pop2(cs ,4.,2.,69.5,2.26e4,dense.xIonDense(26,10))*3.12e-12
	 * dCooldT = dCooldT + c6374 * 2.26e4 * tsq1
	 * call CoolAdd( 'Fe10' , 6374 , c6374 )
	 *
	 * this is the E1 line that can pump [FeX] */
	cs = 0.85*sexp(0.045*1.259e6/phycon.te);
	cs = MAX2(0.05,cs);
	PutCS(cs,&TauLines[ipT352]);
	atom_level2(&TauLines[ipT352]);

	/* Fe XI fe11 fe 11, 7892, 6.08 mic, CS from 
	 * >>refer	fe11	as	Mason, H. 1975, MNRAS 170, 651
	 * A from 
	 * >>refer	fe11	as	Mendoza, C., & Zeippen, C.J. 1983, MNRAS, 202, 981 
	 * following reference has very extensive energy levels and As
	 * >>refer	fe11	as	Fritzsche, S., Dong, C.Z., & Traebert, E., 2000, MNRAS, 318, 263 */
	/*cs = 0.27;*/
	/* >>chng 97 jan 31, set cs to unity, above ignored resonances */
	/*cs = 1.;*/
	/* >>chng 00 dec 05, update Fe11 CS to Tayal 2000 */
	/* >>referold	fe11	cs	Tayal, S.S., 2000, ApJ, 544, 575-580 */
	/* this is the low to mid transition */
	/*cs = 0.0282 * phycon.te30*phycon.te05*phycon.te01;*/
	/* CS is about 2.0 across broad temp range in following reference:
	 * >>refer	Fe11	cs	Aggarwal, K.M., & Keenan, F.P. 2003, MNRAS, 338, 412 */
	/*cs = 2.;
	PutCS(cs,&TauLines[ipTFe07]);*/
	/* Tayal value is 10x larger than previous values */
	/* Aggarwal & Keenan is about same as Tayal */
	/*cs = 0.48; */
	/*cs = 0.50;
	PutCS(cs,&TauLines[ipTFe61]);*/
	/* Tayal value is 3x larger than previous values */
	/*cs = 0.35; tayal number */
	/*cs = 0.55;
	PutCS(cs,&TauDummy);
	atom_level3(&TauLines[ipTFe07],&TauLines[ipTFe61],&TauDummy);*/

	/* >>refer	fe11	cs	Kafatos, M., & Lynch, J.P. 1980, ApJS, 42, 611 */
	/*CoolHeavy.c1467 = atom_pop3(9.,5.,1.,.24,.028,.342,100.,1012.,9.3,5.43e4,
	  6.19e4,&p2,dense.xIonDense[ipIRON][11-1],0.,0.,0.)*1012.*1.36e-11;
	CoolHeavy.c2649 = p2*91.0*7.512e-12;*/
	/*CoolAdd("Fe11",1467,CoolHeavy.c1467);
	CoolAdd("Fe11",2469,CoolHeavy.c2649);*/

	/* >>chng 05 dec 18, use give level Fe 11 model */
	Fe11Lev5();

	/* A's for 2-1 and 3-1 from 
	 * >>refer	fe12	as	Tayal, S.S., & Henry, R.J.W. 1986, ApJ, 302, 200
	 * CS fro 
	 * >>refer	fe12	cs	Tayal, S.S., Henry, R.J.W., Pradhan, A.K. 1987, ApJ, 319, 951
	 * a(3-2) scaled from both ref */
	CoolHeavy.c2568 = atom_pop3(4.,10.,6.,0.72,0.69,2.18,8.1e4,1.84e6,1.33e6,
	  6.37e4,4.91e4,&p2,dense.xIonDense[ipIRON][12-1],0.,0.,0.)*1.33e6*6.79e-12;
	CoolAdd("Fe12",2568,CoolHeavy.c2568);
	CoolHeavy.c1242 = CoolHeavy.c2568*2.30*1.38;
	CoolAdd("Fe12",1242,CoolHeavy.c1242);
	CoolHeavy.c2170 = p2*8.09e4*8.82e-12;
	CoolAdd("Fe12",2170,CoolHeavy.c2170);

	/* [Fe XIII] fe13 fe 13 1.07, 1.08 microns */
	/* >>chng 00 dec 05, update Fe13 CS to Tayal 2000 */
	/* >>refer	fe13	cs	Tayal, S.S., 2000, ApJ, 544, 575-580
	cs = 5.08e-3 * phycon.te30* phycon.te10;
	PutCS(cs,&TauLines[ipFe1310]);
	PutCS(2.1,&TauLines[ipFe1311]);
	PutCS(.46,&TauDummy);
	atom_level3(&TauLines[ipFe1310],&TauLines[ipFe1311],&TauDummy); */

	/* Fe13 lines from 1D and 1S -
	 >>chng 01 aug 07, added these lines */
	/* >>refer	fe11	cs	Tayal, S.S., 2000, ApJ, 544, 575-580 */
	/* >>refer	fe13	as	Shirai, T., Sugar, J., Musgrove, A., & Wiese, W.L., 2000., 
	 * >>refercon	J Phys Chem Ref Data Monograph 8 */
	/* POP3(G1,G2,G3,O12,O13,O23,A21,A31,A32,E12,E23,P2,ABUND,GAM2) 
	CoolHeavy.fe13_1216 = atom_pop3(  9.,5.,1.,  5.08,0.447,0.678,  103.2 ,1010.,7.99,
	  48068.,43440.,&p2,dense.xIonDense[ipIRON][13-1],0.,0.,0.)*1010.*
	  1.64e-11;
	CoolHeavy.fe13_2302 = CoolHeavy.fe13_1216*0.528*0.00791;
	CoolHeavy.fe13_3000 = p2*103.2*7.72e-12;*/
	/*CoolAdd("Fe13",1216,CoolHeavy.fe13_1216);
	CoolAdd("Fe13",3000,CoolHeavy.fe13_3000);
	CoolAdd("Fe13",2302,CoolHeavy.fe13_2302);*/

	/* >>chng 05 dec 18, use give level Fe 13 model */
	Fe13Lev5();

	/* Fe XIV 5303, cs from 
	 * >>refer	Fe14	cs	Storey, P.J., Mason, H.E., Saraph, H.E., 1996, A&A, 309, 677 */
	fe14cs(phycon.alogte,&cs);
	/* >>chng 97 jan 31, set cs to unity, above is VERY large, >10 */
	/* >>chng 01 may 30, back to their value, as per discussion at Lexington 2000 */
	/* >>chng 05 aug 04, following line commented out, had set 5303 to constant value */
	/*cs = 3.1;*/
	/* >>chng 05 dec 22, A from 60.3 to 59.7, experimental value 59.7 +/- 0.4 in
	 * >>refer	Fe14	as	Trabert, E. 2004, A&A, 415, L39 */
	CoolHeavy.c5303 = atom_pop2(cs,2.,4.,59.7,2.71e4,dense.xIonDense[ipIRON][14-1])*
	  3.75e-12;
	thermal.dCooldT += CoolHeavy.c5303*2.71e4*thermal.tsq1;
	CoolAdd("Fe14",5303,CoolHeavy.c5303);

	// [Fe XVII] 17.1 A
	atom_level2(&TauLines[ipFe17_17]);

	/* Fe 18 974.8A;, cs from 
	 * >>referold	fe18	cs	Saraph, H.E. & Tully, J.A. 1994, A&AS, 107, 29 */
	/* >>refer	fe18	cs	Berrington,K.A.,Saraph, H.E. & Tully, J.A. 1998, A&AS, 129, 161 */
	/*>>chng 06 jul 19 Changes made-Humeshkar Nemala*/
	/*cs = MIN2(0.143,0.804/(phycon.te10*phycon.te05));*/
	if( phycon.te < 1.29E6 )
	{
		cs = (realnum)(0.3465/((phycon.te10/phycon.te01)*phycon.te001*phycon.te0003));
	}
	else if( phycon.te < 5.135E6 )
	{
		cs = (realnum)((1.1062E-02)*phycon.te10*phycon.te05*phycon.te003*phycon.te0005);
	}
	else
	{
		cs = (realnum)((60.5728)/(phycon.te40*phycon.te003*phycon.te0001*phycon.te0005));
	}

	PutCS(cs,&TauLines[ipFe18975]);
	atom_level2(&TauLines[ipFe18975]);

	/* O-like Fe19, 3P ground term, 7046.72A vac wl, 1328.90A */
	/* >>refer	fe19	cs	Butler, K., & Zeippen, C.J., 2001, A&A, 372, 1083 */
	cs12 = 0.0627 / phycon.te03;
	cs13 = 0.692 /(phycon.te10*phycon.te01);
	cs23 = 0.04;
	/*CoolHeavy.c7082 = atom_pop3(5.,1.,3.,0.0132,0.0404,0.0137,0.505,1.46e4,*/
	/* POP3(G1,G2,G3,O12,O13,O23,A21,A31,A32,E12,E23,P2,ABUND,GAM2) */
	CoolHeavy.c7082 = atom_pop3(5.,1.,3.,cs12,cs13,cs23,0.505,1.46e4,
	  41.2,1.083e5,2.03e4,&p2,dense.xIonDense[ipIRON][19-1],0.,0.,0.)*41.2*
	  2.81e-12;
	CoolHeavy.c1118 = CoolHeavy.c7082*354.4*6.335;
	CoolHeavy.c1328 = p2*0.505*1.50e-11;
	CoolAdd("Fe19",7082,CoolHeavy.c7082);
	CoolAdd("Fe19",1118,CoolHeavy.c1118);
	CoolAdd("Fe19",1328,CoolHeavy.c1328);

	CoolHeavy.c592 = atom_pop2(0.0913,9.,5.,1.64e4,2.428e5,dense.xIonDense[ipIRON][19-1])*
	  3.36e-11;
	CoolAdd("Fe19",592,CoolHeavy.c592);

	/* >>chng 01 aug 10, add this line */
	/* Fe20 721.40A, 578*/
	/* >>refer	fe20	cs	Butler, K., & Zeippen, C.J., 2001, A&A, 372, 1078 */
	/*>>refer	fe20	as	Merkelis, G., Martinson, I., Kisielius, R., & Vilkas, M.J., 1999, 
	  >>refercon	Physica Scripta, 59, 122 */
	cs = 1.17 /(phycon.te20/phycon.te01);
	PutCS(cs , &TauLines[ipTFe20_721]);
	cs = 0.248 /(phycon.te10/phycon.te01);
	PutCS(cs , &TauLines[ipTFe20_578]);
	cs = 0.301 /(phycon.te10/phycon.te02);
	PutCS(cs , &TauDummy);
	atom_level3(&TauLines[ipTFe20_721],&TauLines[ipTFe20_578],&TauDummy);

	/* >>chng 00 oct 26, much larger CS for following transition */
	/* Fe 21 1354, 2299 A, cs eval at 1e6K
	 * >>refer	Fe21	cs	Aggarwall, K.M., 1991, ApJS 77, 677 */
	PutCS(0.072,&TauLines[ipTFe13]);
	PutCS(0.269,&TauLines[ipTFe23]);
	PutCS(0.055,&TauDummy);
	atom_level3(&TauLines[ipTFe13],&TauLines[ipTFe23],&TauDummy);

	/*>>refer	Fe22	energy	Feldman, U., Curdt, W., Landi, E., Wilhelm, K., 2000, ApJ, 544, 508 */
	/*>>chng 00 oct 26, added these fe 21 lines, removed old CoolHeavy.fs846, the lowest */

	static vector< pair<transition*,double> > Fe22Pump;
	Fe22Pump.reserve(96);

	/* one time initialization if first call */
	if( Fe22Pump.empty() )
	{
		// set up level 2 pumping lines
		for( i=0; i < nWindLine; ++i )
		{
			/* don't test on nelem==ipIRON since lines on physics, not C, scale */
			if( TauLine2[i].Hi->nelem == 26 && TauLine2[i].Hi->IonStg == 22 )
			{
#				if	0
				DumpLine( &TauLine2[i] );
#				endif
				double branch_ratio;
				// the branching ratios used here ignore cascades via intermediate levels
				// usually the latter are much slower, so this should be reasonable
				if( fp_equal( TauLine2[i].Hi->g, realnum(2.) ) )
					branch_ratio = 2./3.; // 2S upper level
				else if( fp_equal( TauLine2[i].Hi->g, realnum(6.) ) )
					branch_ratio = 1./2.; // 2P upper level
				else if( fp_equal( TauLine2[i].Hi->g, realnum(10.) ) )
					branch_ratio = 1./6.; // 2D upper level
				else
					TotalInsanity();
				pair<transition*,double> pp2( &TauLine2[i], branch_ratio ); 
				Fe22Pump.push_back( pp2 );
			}
		}
	}

	/* now sum pump rates */
	double Fe22_pump_rate = 0.;
	vector< pair<transition*,double> >::const_iterator fe22p;
	for( fe22p=Fe22Pump.begin(); fe22p != Fe22Pump.end(); ++fe22p )
	{
		const transition* t = fe22p->first;
		double branch_ratio = fe22p->second;
		Fe22_pump_rate += t->Emis->pump*branch_ratio;
#		if	0
		dprintf( ioQQQ, "Fe XXII %.3e %.3e\n",
			 t->WLAng , t->Emis->pump*branch_ratio );
#		endif
	}

	/*AtomSeqBoron compute cooling from 5-level boron sequence model atom */
	/*>>refer	fe22	cs	Zhang, H.L., Graziani, M., & Pradhan, A.K., 1994, A&A 283, 319*/
	/*>>refer	fe22	as	Dankwort, W., & Trefftz, E., 1978, A&A 65, 93-98 */
	AtomSeqBoron(&TauLines[ipFe22_846], 
	  &TauLines[ipFe22_247], 
	  &TauLines[ipFe22_217], 
	  &TauLines[ipFe22_348], 
	  &TauLines[ipFe22_292], 
	  &TauLines[ipFe22_253], 
	  /*double cs40,
		double cs32,
		double cs42,
		double cs43, */
	  0.00670 , 0.0614 , 0.0438 , 0.122 , Fe22_pump_rate ,"Fe22");

	/* Fe 22 845.6, C.S.=guess, A from
	 * >>refer	fe22	as	Froese Fischer, C. 1983, J.Phys. B, 16, 157 
	CoolHeavy.fs846 = atom_pop2(0.2,2.,4.,1.5e4,1.701e5,dense.xIonDense[ipIRON][22-1])*
	  2.35e-11;
	CoolAdd("Fe22",846,CoolHeavy.fs846);*/

	/** \todo	2	update atomic data to Chidichimo et al 1999 AASup 137, 175*/
	/* FE 23 262.6, 287A, 1909-LIKE, 
	 * cs from 
	 * >>refer	fe23	cs	Bhatia, A.K., & Mason, H.E. 1986, A&A, 155, 413 */
	CoolHeavy.c263 = atom_pop2(0.04,1.,9.,1.6e7,5.484e5,dense.xIonDense[ipIRON][23-1])*
	  7.58e-11;
	CoolAdd("Fe23",262,CoolHeavy.c263);

	/* FE 24, 255, 192 Li seq doublet
	 * >>refer	fe24	cs	Cochrane, D.M., & McWhirter, R.W.P. 1983, PhyS, 28, 25 */
	ligbar(26,&TauLines[ipT192],&TauLines[ipT11],&cs2s2p,&cs2s3p);

	PutCS(cs2s2p,&TauLines[ipT192]);
	atom_level2(&TauLines[ipT192]);

	/* funny factor (should have been 0.5) due to energy change */
	PutCS(cs2s2p*0.376,&TauLines[ipT255]);
	atom_level2(&TauLines[ipT255]);

	PutCS(cs2s3p,&TauLines[ipT11]);
	atom_level2(&TauLines[ipT11]);

	/** \todo	2		 * following not in cooling function */
	TauLines[ipT353].Emis->PopOpc = dense.xIonDense[ipIRON][11-1];
	TauLines[ipT353].Lo->Pop = dense.xIonDense[ipIRON][11-1];
	TauLines[ipT353].Hi->Pop = 0.;
	TauLines[ipT347].Emis->PopOpc = dense.xIonDense[ipIRON][14-1];
	TauLines[ipT347].Lo->Pop = dense.xIonDense[ipIRON][14-1];
	TauLines[ipT347].Hi->Pop = 0.;

	return;
}

/*fe14cs compute collision strength for forbidden transition 
 * w/in Fe XIV ground term. From 
 * >>refer	fe14	cs	Storey, P.J., Mason, H.E., Saraph, H.E., 1996, A&A, 309, 677
 * */
STATIC void fe14cs(double te1, 
	  double *csfe14)
{
	double a, 
	  b, 
	  c, 
	  d, 
	  telog1, 
	  telog2, 
	  telog3;

	DEBUG_ENTRY( "fe14cs()" );

	/* limit range in log T: */
	telog1 = te1;
	telog1 = MIN2(9.0,telog1);
	telog1 = MAX2(4.0,telog1);

	/* compute square and cube */
	telog2 = telog1*telog1;
	telog3 = telog2*telog1;

	/* compute cs:
	 * */
	if( telog1 <= 5.0 )
	{
		a = 557.05536;
		b = -324.56109;
		c = 63.437974;
		d = -4.1365147;
		*csfe14 = a + b*telog1 + c*telog2 + d*telog3;
	}
	else
	{
		a = 0.19515493;
		b = 2.9404407;
		c = 4.9578944;
		d = 0.79887506;
		*csfe14 = a + b*exp(-0.5*((telog1-c)*(telog1-c)/d));
	}
	return;
}

/*Fe4Lev12 compute populations and cooling due to 12 level Fe IV ion */
STATIC void Fe4Lev12(void)
{
	const int NLFE4 = 12;
	bool lgZeroPop;
	int nNegPop;
	long int i, 
	  j;
	static bool lgFirst=true;

	double dfe4dt;

	/*static long int **ipdest; */
	static double 
		**AulEscp,
		**col_str,
		**AulDest, 
		depart[NLFE4],
		pops[NLFE4], 
		destroy[NLFE4], 
		create[NLFE4],
		**CollRate,
		**AulPump;

	static const double Fe4A[NLFE4][NLFE4] = {
		{0.,0.,0.,1.e-5,0.,1.368,.89,0.,1.3e-3,1.8e-4,.056,.028},
		{0.,0.,2.6e-8,0.,0.,0.,0.,0.,1.7e-7,0.,0.,0.},
		{0.,0.,0.,0.,3.5e-7,6.4e-10,0.,0.,6.315e-4,0.,6.7e-7,0.},
		{0.,0.,0.,0.,1.1e-6,6.8e-5,8.6e-6,3.4e-10,7.6e-5,1.e-7,5.8e-4,2.8e-4},
		{0.,0.,0.,0.,0.,1.5e-5,1.3e-9,0.,7.6e-4,0.,1.1e-6,6.0e-7},
		{0.,0.,0.,0.,0.,0.,1.1e-5,1.2e-13,.038,9.9e-7,.022,.018},
		{0.,0.,0.,0.,0.,0.,0.,3.7e-5,2.9e-6,.034,3.5e-3,.039},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,.058,3.1e-6,1.4e-3},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.3e-4,3.1e-14},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.9e-19,1.0e-5},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.3e-7},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}
	};
	static const double Fe4CS[NLFE4][NLFE4] = {
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
		{0.98,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
		{0.8167,3.72,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
		{0.49,0.0475,0.330,0.,0.,0.,0.,0.,0.,0.,0.,0.},
		{0.6533,0.473,2.26,1.64,0.,0.,0.,0.,0.,0.,0.,0.},
		{0.45,0.686,0.446,0.106,0.254,0.,0.,0.,0.,0.,0.,0.},
		{0.30,0.392,0.152,0.269,0.199,0.605,0.,0.,0.,0.,0.,0.},
		{0.15,0.0207,0.190,0.0857,0.166,0.195,0.327,0.,0.,0.,0.,0.},
		{0.512,1.23,0.733,0.174,0.398,0.623,0.335,0.102,0.,0.,0.,0.},
		{0.128,0.0583,0.185,0.200,0.188,0.0835,0.127,0.0498,0.0787,0.,0.,0.},
		{0.384,0.578,0.534,0.363,0.417,0.396,0.210,0.171,0.810,0.101,0.,0.},
		{0.256,0.234,0.306,0.318,0.403,0.209,0.195,0.112,0.195,0.458,0.727,0.}
	};

	static const double gfe4[NLFE4]={6.,12.,10.,6.,8.,6.,4.,2.,8.,2.,6.,4.};

	/* excitation energies in Kelvin 
	static double ex[NLFE4]={0.,46395.8,46464.,46475.6,46483.,50725.,
	  50838.,50945.,55796.,55966.,56021.,56025.};*/
	/*>>refer	Fe3	energies	version 3 NIST Atomic Spectra Database */
	/*>>chng 05 dec 17, from Kelvin above to excitation energies in wn */
	static const double excit_wn[NLFE4]={0.,32245.5,32292.8,32301.2,32305.7,35253.8,
	  35333.3,35406.6,38779.4,38896.7,38935.1,38938.2};

	DEBUG_ENTRY( "Fe4Lev12()" );

	if( lgFirst )
	{
		/* will never do this again */
		lgFirst = false;
		/* allocate the 1D arrays*/
		/* create space for the 2D arrays */
		AulPump = ((double **)MALLOC((NLFE4)*sizeof(double *)));
		CollRate = ((double **)MALLOC((NLFE4)*sizeof(double *)));
		AulDest = ((double **)MALLOC((NLFE4)*sizeof(double *)));
		AulEscp = ((double **)MALLOC((NLFE4)*sizeof(double *)));
		col_str = ((double **)MALLOC((NLFE4)*sizeof(double *)));
		for( i=0; i < NLFE4; ++i )
		{
			AulPump[i] = ((double *)MALLOC((NLFE4)*sizeof(double )));
			CollRate[i] = ((double *)MALLOC((NLFE4)*sizeof(double )));
			AulDest[i] = ((double *)MALLOC((NLFE4)*sizeof(double )));
			AulEscp[i] = ((double *)MALLOC((NLFE4)*sizeof(double )));
			col_str[i] = ((double *)MALLOC((NLFE4)*sizeof(double )));
		}
	}

	/* bail if no Fe+3 */
	if( dense.xIonDense[ipIRON][3] <= 0. )
	{
		fe.Fe4CoolTot = 0.;
		fe.fe40401 = 0.;
		fe.fe42836 = 0.;
		fe.fe42829 = 0.;
		fe.fe42567 = 0.;
		fe.fe41207 = 0.;
		fe.fe41206 = 0.;
		fe.fe41106 = 0.;
		fe.fe41007 = 0.;
		fe.fe41008 = 0.;
		fe.fe40906 = 0.;
		CoolAdd("Fe 4",0,0.);

		/* level populations */
		/* LIMLEVELN is the dimension of the atoms vectors */
		ASSERT( NLFE4 <= LIMLEVELN);
		for( i=0; i < NLFE4; i++ )
		{
			atoms.PopLevels[i] = 0.;
			atoms.DepLTELevels[i] = 1.;
		}
		return;
	}
	/* number of levels in model ion */

	/* these are in wavenumbers
	 * data excit_wn/ 0., 32245.5, 32293., 32301.2, 
	 *  1  32306., 35254., 35333., 35407., 38779., 38897., 38935.,
	 *  2  38938./ 
	 * excitation energies in Kelvin */

	/* A's are from Garstang, R.H., MNRAS 118, 572 (1958).
	 * each set is for a lower level indicated by second element in array,
	 * index runs over upper level
	 * A's are saved into arrays as data(up,lo) */

	/* collision strengths from Berrington and Pelan  Ast Ap S 114, 367.
	 * order is cs(low,up) */

	/* all elements are used, and must be set to zero if zero */
	for( i=0; i < NLFE4; i++ )
	{
		create[i] = 0.;
		destroy[i] = 0.;
		for( j=0; j < NLFE4; j++ )
		{
			/*data[j][i] = 1e33;*/
			col_str[j][i] = 0.;
			AulEscp[j][i] = 0.;
			AulDest[j][i] = 0.;
			AulPump[j][i] = 0.;
		}
	}

	/* fill in einstein As and collision strengths */
	for( i=0; i < NLFE4; i++ )
	{
		for( j=i+1; j < NLFE4; j++ )
		{
			/*data[i][j] = Fe4A[i][j];*/
			AulEscp[j][i] = Fe4A[i][j];
			/*data[j][i] = Fe4CS[j][i];*/
			col_str[j][i] = Fe4CS[j][i];
		}
	}

	/* leveln itself is well-protected against zero abundances,
	 * low temperatures */

	atom_levelN(NLFE4,
		dense.xIonDense[ipIRON][3],
		gfe4,
		excit_wn,
		'w',
		pops,
		depart,
		&AulEscp ,
		&col_str ,
		&AulDest,
		&AulPump,
		&CollRate,
		create,
		destroy,
		/* say atom_levelN should evaluate coll rates from cs */
		false,
		&fe.Fe4CoolTot,
		&dfe4dt,
		"FeIV",
		/* nNegPop positive if negative pops occured, negative if too cold */
		&nNegPop,
		&lgZeroPop,
		false );

	/* LIMLEVELN is the dim of the PopLevels vector */
	ASSERT( NLFE4 <= LIMLEVELN );
	for( i=0; i < NLFE4; ++i )
	{
		atoms.PopLevels[i] = pops[i];
		atoms.DepLTELevels[i] = depart[i];
	}

	if( nNegPop > 0 )
	{
		fprintf( ioQQQ, " fe4levl2 found negative populations\n" );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	CoolAdd("Fe 4",0,fe.Fe4CoolTot);

	thermal.dCooldT += dfe4dt;

	/* three transitions hst can observe
	 * 4 -1 3095.9A and 5 -1 3095.9A */
	fe.fe40401 = (pops[3]*Fe4A[0][3]*(excit_wn[3] - excit_wn[0]) + 
		pops[4]*Fe4A[0][4]*(excit_wn[4] - excit_wn[0]) )*T1CM*BOLTZMANN;

	fe.fe42836 = pops[5]*Fe4A[0][5]*(excit_wn[5] - excit_wn[0])*T1CM*BOLTZMANN;

	fe.fe42829 = pops[6]*Fe4A[0][6]*(excit_wn[5] - excit_wn[0])*T1CM*BOLTZMANN;

	fe.fe42567 = (pops[10]*Fe4A[0][10]*(excit_wn[10] - excit_wn[0]) + 
		pops[11]*Fe4A[0][11]*(excit_wn[10] - excit_wn[0]))*T1CM*BOLTZMANN;

	fe.fe41207 = pops[11]*Fe4A[6][11]*(excit_wn[11] - excit_wn[6])*T1CM*BOLTZMANN;
	fe.fe41206 = pops[11]*Fe4A[5][11]*(excit_wn[11] - excit_wn[5])*T1CM*BOLTZMANN;
	fe.fe41106 = pops[10]*Fe4A[5][10]*(excit_wn[10] - excit_wn[5])*T1CM*BOLTZMANN;
	fe.fe41007 = pops[9]*Fe4A[6][9]*(excit_wn[9] - excit_wn[6])*T1CM*BOLTZMANN;
	fe.fe41008 = pops[9]*Fe4A[7][9]*(excit_wn[9] - excit_wn[7])*T1CM*BOLTZMANN;
	fe.fe40906 = pops[8]*Fe4A[5][8]*(excit_wn[8] - excit_wn[5])*T1CM*BOLTZMANN;
	return;
}

/*Fe7Lev8 compute populations and cooling due to 8 level Fe VII ion */
STATIC void Fe7Lev8(void)
{
	bool lgZeroPop;
	int nNegPop;
	double scale;
	long int i, 
	  j;
	static bool lgFirst=true;
	static long int ipPump=-1;

	double dfe7dt,
		FUV_pump;

	long int ihi , ilo;
	static double 
	  *depart,
	  *pops, 
	  *destroy, 
	  *create ,
	  **AulDest, 
	  **CollRate,
	  **AulPump,
	  **AulNet, 
	  **col_str;
	/*
	following are FeVII lines predicted in limit_laser_2000
	Fe 7  5721A -24.399   0.1028 
	Fe 7  4989A -24.607   0.0637 
	Fe 7  4894A -24.838   0.0374
	Fe 7  4699A -25.693   0.0052
	Fe 7  6087A -24.216   0.1566
	Fe 7  5159A -24.680   0.0538
	Fe 7  4943A -25.048   0.0231
	Fe 7  3587A -24.285   0.1336
	Fe 7  5277A -25.053   0.0228
	Fe 7  3759A -24.139   0.1870
	Fe 7 3.384m -25.621   0.0062
	Fe 7 2.629m -25.357   0.0113
	Fe 7 9.510m -24.467   0.0878
	Fe 7 7.810m -24.944   0.0293
	*/

	/* statistical weights for the n levels */
	static double gfe7[NLFE7]={5.,7.,9.,5.,1.,3.,5.,9.};

	/*refer	Fe7	energies	Ekbert, J.O. 1981, Phys. Scr 23, 7 */
	/*static double ex[NLFE7]={0.,1047.,2327.,17475.,20037.,20428.,21275.,28915.};*/
	/* excitation energies in wavenumbers 
	 * >>chng 05 dec 15, rest of code had assumed that there were energies in Kelvin
	 * rather than wavenumbers.  Bug caught by Kevin Blagrave 
	 * modified atom_levelN to accept either Kelvin or wavenumbers */
	static double excit_wn[NLFE7]={0. , 1051.5 , 2331.5 , 17475.5 , 20040.3 , 20430.1 , 21278.6 , 28927.3 };

	DEBUG_ENTRY( "Fe7Lev8()" );

	if( lgFirst )
	{
		/* will never do this again */
		lgFirst = false;
		/* allocate the 1D arrays*/
		/* create space for the 2D arrays */
		depart = ((double *)MALLOC((NLFE7)*sizeof(double)));
		pops = ((double *)MALLOC((NLFE7)*sizeof(double)));
		destroy = ((double *)MALLOC((NLFE7)*sizeof(double)));
		create = ((double *)MALLOC((NLFE7)*sizeof(double)));
		/* now the 2-d arrays */
		fe.Fe7_wl = ((double **)MALLOC((NLFE7)*sizeof(double *)));
		fe.Fe7_emiss = ((double **)MALLOC((NLFE7)*sizeof(double *)));
		AulNet = ((double **)MALLOC((NLFE7)*sizeof(double *)));
		col_str = ((double **)MALLOC((NLFE7)*sizeof(double *)));
		AulPump = ((double **)MALLOC((NLFE7)*sizeof(double *)));
		CollRate = ((double **)MALLOC((NLFE7)*sizeof(double *)));
		AulDest = ((double **)MALLOC((NLFE7)*sizeof(double *)));
		for( i=0; i < NLFE7; ++i )
		{
			fe.Fe7_wl[i] = ((double *)MALLOC((NLFE7)*sizeof(double )));
			fe.Fe7_emiss[i] = ((double *)MALLOC((NLFE7)*sizeof(double )));
			AulNet[i] = ((double *)MALLOC((NLFE7)*sizeof(double )));
			col_str[i] = ((double *)MALLOC((NLFE7)*sizeof(double )));
			AulPump[i] = ((double *)MALLOC((NLFE7)*sizeof(double )));
			CollRate[i] = ((double *)MALLOC((NLFE7)*sizeof(double )));
			AulDest[i] = ((double *)MALLOC((NLFE7)*sizeof(double )));
		}

		/* set some to constant values after zeroing out */
		for( i=0; i < NLFE7; ++i )
		{
			create[i] = 0.;
			destroy[i] = 0.;
			for( j=0; j < NLFE7; ++j )
			{
				AulNet[i][j] = 0.;
				col_str[i][j] = 0.;
				CollRate[i][j] = 0.;
				AulDest[i][j] = 0.;
				AulPump[i][j] = 0.;
				fe.Fe7_wl[i][j] = 0.;
				fe.Fe7_emiss[i][j] = 0.;
			}
		}
		set_NaN( fe.Fe7_wl[2][1] );
		set_NaN( fe.Fe7_emiss[2][1] );
		set_NaN( fe.Fe7_wl[1][0] );
		set_NaN( fe.Fe7_emiss[1][0] );
		/* two levels within ground term must never be addressed in this array -
		 * they are fully transferred */
		for( ilo=0; ilo < NLFE7-1; ++ilo )
		{
			/* must not do 1-0 or 2-1, which are transferred lines */
			for( ihi=MAX2(3,ilo+1); ihi < NLFE7; ++ihi )
			{
				fe.Fe7_wl[ihi][ilo] = 1e8/(excit_wn[ihi]-excit_wn[ilo]) / RefIndex( excit_wn[ihi]-excit_wn[ilo] );
			}
		}

		/* assume FeVII is optically thin - just use As as net escape */
		/*>>refer	Fe7	As	Berrington, K.A., Nakazaki, S., & Norrington, P.H. 2000, A&AS, 142, 313 */
		AulNet[1][0] = 0.0325;
		AulNet[2][0] = 0.167e-8;
		/* following corrected from 3.25 to 0.372 as per Keith Berrington email */
		AulNet[3][0] = 0.372;
		AulNet[4][0] = 0.135;
		AulNet[5][0] = 0.502e-1;
		/* following corrected from 0.174 to 0.0150 as per Keith Berrington email */
		AulNet[6][0] = 0.0150;
		AulNet[7][0] = 0.959e-3;

		AulNet[2][1] = 0.466e-1;
		AulNet[3][1] = 0.603;
		AulNet[5][1] = 0.762e-1;
		AulNet[6][1] = 0.697e-1;
		AulNet[7][1] = 0.343;

		AulNet[3][2] = 0.139e-2;
		AulNet[6][2] = 0.735e-1;
		AulNet[7][2] = 0.503;

		AulNet[4][3] = 0.472e-6;
		AulNet[5][3] = 0.572e-1;
		/* following corrected from 0.191 to 0.182 as per Keith Berrington email */
		AulNet[6][3] = 0.182;
		AulNet[7][3] = 0.414e-2;

		AulNet[5][4] = 0.115e-2;
		AulNet[6][4] = 0.139e-7;

		AulNet[6][5] = 0.743e-2;

		AulNet[7][6] = 0.454e-4;

		/* collision strengths at log T 4.5 */
		/** \todo	2	put in temperature dependence */
		/*>>refer	Fe7	cs	Berrington, K.A., Nakazaki, S., & Norrington, P.H. 2000, A&AS, 142, 313 */
#		if 0
		if( fudge(-1) )
		{
			fprintf(ioQQQ,"DEBUG fudge call cool_iron\n");
			scale = fudge(0);
		}
		else
#		endif
		scale = 1.;

		col_str[1][0] = 3.35;
		col_str[2][0] = 1.17;
		col_str[3][0] = 0.959;
		col_str[4][0] = 0.299;
		col_str[5][0] = 0.633;
		col_str[6][0] = 0.549;
		col_str[7][0] = 1.24*scale;

		col_str[2][1] = 4.11;
		col_str[3][1] = 1.29;
		col_str[4][1] = 0.235;
		col_str[5][1] = 0.833;
		col_str[6][1] = 1.06;
		col_str[7][1] = 1.74*scale;

		col_str[3][2] = 1.60;
		col_str[4][2] = 0.187;
		col_str[5][2] = 0.690;
		col_str[6][2] = 1.94;
		col_str[7][2] = 2.25*scale;

		col_str[4][3] = 0.172;
		col_str[5][3] = 0.531;
		col_str[6][3] = 1.06;
		col_str[7][3] = 2.02;

		col_str[5][4] = 0.370;
		col_str[6][4] = 0.324;
		col_str[7][4] = 0.164;

		col_str[6][5] = 1.17;
		col_str[7][5] = 0.495;

		col_str[7][6] = 0.903;

		/* check whether level 2 lines are on, and if so, find the driving lines that
		 * can pump the upper levels of this atom */
		if( nWindLine > 0 )
		{
			ipPump = -1;
			for( i=0; i < nWindLine; ++i )
			{
				/* don't test on nelem==ipIRON since lines on physics, not C, scale */
				if( TauLine2[i].Hi->nelem == 26 && TauLine2[i].Hi->IonStg == 7 &&
					/* >>chng 05 jul 10, move line to wavelength that agrees with nist tables
					 * here and in level2 data file 
					 * NB must add a few wn to second number to get hit -
					 * logic is that this is lowest E1 transition */
					(TauLine2[i].EnergyWN-4.27360E+05) < 0. )
				{
					ipPump = i;
					break;
				}
			}
			if( ipPump < 0 )
			{
				fprintf(ioQQQ,"PROBLEM Fe7Lev8 cannot identify the FUV driving line.\n");
				TotalInsanity();
			}
		}
		else
		{
			ipPump = 0;
		}
	}

	/* bail if no ions */
	if( dense.xIonDense[ipIRON][6] <= 0. )
	{
		CoolAdd("Fe 7",0,0.);

		fe.Fe7CoolTot = 0.;
		for( ilo=0; ilo < NLFE7-1; ++ilo )
		{
			/* must not do 1-0 or 2-1, which are transferred lines */
			for( ihi=MAX2(3,ilo+1); ihi < NLFE7; ++ihi )
			{
				fe.Fe7_emiss[ihi][ilo] = 0.;
			}
		}
		TauLines[ipFe0795].Lo->Pop = 0.;
		TauLines[ipFe0795].Emis->PopOpc = 0.;
		TauLines[ipFe0795].Hi->Pop = 0.;
		TauLines[ipFe0795].Emis->xIntensity = 0.;
		TauLines[ipFe0795].Coll.cool = 0.;
		TauLines[ipFe0795].Emis->phots = 0.;
		TauLines[ipFe0795].Emis->ColOvTot = 0.;
		TauLines[ipFe0795].Coll.heat = 0.;
		CoolAdd( "Fe 7", TauLines[ipFe0795].WLAng , 0.);
		TauLines[ipFe0778].Lo->Pop = 0.;
		TauLines[ipFe0778].Emis->PopOpc = 0.;
		TauLines[ipFe0778].Hi->Pop = 0.;
		TauLines[ipFe0778].Emis->xIntensity = 0.;
		TauLines[ipFe0778].Coll.cool = 0.;
		TauLines[ipFe0778].Emis->phots = 0.;
		TauLines[ipFe0778].Emis->ColOvTot = 0.;
		TauLines[ipFe0778].Coll.heat = 0.;
		CoolAdd( "Fe 7", TauLines[ipFe0778].WLAng , 0.);
		/* level populations */
		/* LIMLEVELN is the dimension of the atoms vectors */
		ASSERT( NLFE7 <= LIMLEVELN);
		for( i=0; i < NLFE7; i++ )
		{
			atoms.PopLevels[i] = 0.;
			atoms.DepLTELevels[i] = 1.;
		}
		return;
	}

	/* do pump rate for FUV excitation of 3P (levels 5-7 on physics scale, not C scale) */
	if( ipPump )
	{
		FUV_pump = TauLine2[ipPump].Emis->pump * 0.3 /(0.3+TauLine2[ipPump].Emis->Pesc);
	}
	else
	{
		FUV_pump = 0.;
	}

	/* this represents photoexcitation of 3P from ground level
	 * >>chng 04 nov 22, upper array elements were on physics not c scale, off by one, TE */
	AulPump[0][4] = FUV_pump;
	AulPump[1][4] = FUV_pump;
	AulPump[2][4] = FUV_pump;
	AulPump[0][5] = FUV_pump;
	AulPump[1][5] = FUV_pump;
	AulPump[2][5] = FUV_pump;
	AulPump[0][6] = FUV_pump;
	AulPump[1][6] = FUV_pump;
	AulPump[2][6] = FUV_pump;
	/* these were set in the previous call to atom_levelN as the previous pump times
	 * ratio of statistical weights, so this is the downward pump rate */
	AulPump[4][0] = 0;
	AulPump[4][1] = 0;
	AulPump[4][2] = 0;
	AulPump[5][0] = 0;
	AulPump[5][1] = 0;
	AulPump[5][2] = 0;
	AulPump[6][0] = 0;
	AulPump[6][1] = 0;
	AulPump[6][2] = 0;

	/* within ground term update to rt results */
	AulNet[1][0] = TauLines[ipFe0795].Emis->Aul*(TauLines[ipFe0795].Emis->Pesc + TauLines[ipFe0795].Emis->Pelec_esc);
	AulDest[1][0] = TauLines[ipFe0795].Emis->Aul*TauLines[ipFe0795].Emis->Pdest;
	AulPump[0][1] = TauLines[ipFe0795].Emis->pump;
	AulPump[1][0] = 0.;

	AulNet[2][1] = TauLines[ipFe0778].Emis->Aul*(TauLines[ipFe0778].Emis->Pesc + TauLines[ipFe0778].Emis->Pelec_esc);
	AulDest[2][1] = TauLines[ipFe0778].Emis->Aul*TauLines[ipFe0778].Emis->Pdest;
	AulPump[1][2] = TauLines[ipFe0778].Emis->pump;
	AulPump[2][1] = 0.;

	/* nNegPop positive if negative pops occurred, negative if too cold */
	atom_levelN(
		/* number of levels */
		NLFE7,
		/* the abundance of the ion, cm-3 */
		dense.xIonDense[ipIRON][6],
		/* the statistical weights */
		gfe7,
		/* the excitation energies in wavenumbers */
		excit_wn,
		/* units are wavenumbers */
		'w',
		/* the derived populations - cm-3 */
		pops,
		/* the derived departure coefficients */
		depart,
		/* the net emission rate, Aul * escp prob */
		&AulNet ,
		/* the collision strengths */
		&col_str ,
		/* A * destruction prob */
		&AulDest,
		/* pumping rate */
		&AulPump,
		/* collision rate, s-1, must defined if no collision strengths */
		&CollRate,
		/* creation vector */
		create,
		/* destruction vector */
		destroy,
		/* say atom_levelN should evaluate coll rates from cs */
		false,
		&fe.Fe7CoolTot,
		&dfe7dt,
		"Fe 7",
		&nNegPop,
		&lgZeroPop,
		false );

	/* LIMLEVELN is the dim of the PopLevels vector */
	ASSERT( NLFE7 <= LIMLEVELN );
	for( i=0; i < NLFE7; ++i )
	{
		atoms.PopLevels[i] = pops[i];
		atoms.DepLTELevels[i] = depart[i];
	}

	if( lgZeroPop )
	{
		/* this case, too cool to excite upper levels of atom, but will want 
		 * to do IR lines in ground term */
		PutCS(col_str[1][0],&TauLines[ipFe0795]);
		PutCS(col_str[2][1],&TauLines[ipFe0778]);
		PutCS(col_str[2][0],&TauDummy);
		atom_level3(&TauLines[ipFe0795],&TauLines[ipFe0778],&TauDummy);
		atoms.PopLevels[0] = TauLines[ipFe0795].Lo->Pop;
		atoms.PopLevels[1] = TauLines[ipFe0795].Hi->Pop;
		atoms.PopLevels[2] = TauLines[ipFe0778].Hi->Pop;
		for( ilo=0; ilo < NLFE7-1; ++ilo )
		{
			/* must not do 1-0 or 2-1, which are transferred lines */
			for( ihi=MAX2(3,ilo+1); ihi < NLFE7; ++ihi )
			{
				fe.Fe7_emiss[ihi][ilo] = 0.;
			}
		}
	}
	else
	{
		/* this case non-zero pops, but must set up vars within transition structure */
		TauLines[ipFe0795].Lo->Pop = atoms.PopLevels[0];
		TauLines[ipFe0795].Hi->Pop = atoms.PopLevels[1];
		TauLines[ipFe0795].Emis->PopOpc = (pops[0] - pops[1]*gfe7[0]/gfe7[1]);
		TauLines[ipFe0795].Emis->xIntensity = TauLines[ipFe0795].Emis->phots*TauLines[ipFe0795].EnergyErg;
		TauLines[ipFe0795].Emis->phots = TauLines[ipFe0795].Emis->Aul*(TauLines[ipFe0795].Emis->Pesc + TauLines[ipFe0795].Emis->Pelec_esc)*pops[1];
		TauLines[ipFe0795].Emis->ColOvTot = CollRate[0][1]/(CollRate[0][1]+TauLines[ipFe0795].Emis->pump);
		TauLines[ipFe0795].Coll.cool = 0.;
		TauLines[ipFe0795].Coll.heat = 0.;

		TauLines[ipFe0778].Lo->Pop = atoms.PopLevels[1];
		TauLines[ipFe0778].Hi->Pop = atoms.PopLevels[2];
		TauLines[ipFe0778].Emis->PopOpc = (pops[1] - pops[2]*gfe7[1]/gfe7[2]);
		TauLines[ipFe0778].Emis->xIntensity = TauLines[ipFe0778].Emis->phots*TauLines[ipFe0778].EnergyErg;
		TauLines[ipFe0778].Emis->phots = TauLines[ipFe0778].Emis->Aul*(TauLines[ipFe0778].Emis->Pesc + TauLines[ipFe0778].Emis->Pelec_esc)*pops[2];
		TauLines[ipFe0778].Emis->ColOvTot = CollRate[1][2]/(CollRate[1][2]+TauLines[ipFe0778].Emis->pump);
		TauLines[ipFe0778].Coll.heat = 0.;
		TauLines[ipFe0778].Coll.cool = 0.;
	}

	if( nNegPop > 0 )
	{
		fprintf( ioQQQ, "PROBLEM Fe7Lev8 found negative populations\n" );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	/* add cooling then its derivative */
	CoolAdd("Fe 7",0,fe.Fe7CoolTot);
	/* derivative of cooling */
	thermal.dCooldT += dfe7dt;

	/* find emission in each line */
	for( ilo=0; ilo < NLFE7-1; ++ilo )
	{
		/* must not do 1-0 or 2-1, which are transferred lines */
		for( ihi=MAX2(3,ilo+1); ihi < NLFE7; ++ihi )
		{
			/* emission in these lines */
			fe.Fe7_emiss[ihi][ilo] = pops[ihi]*AulNet[ihi][ilo]*(excit_wn[ihi] - excit_wn[ilo])*ERG1CM;
		}
	}
	return;
}

/*Fe3Lev14 compute populations and cooling due to 14 level Fe III ion 
 * >>chng 05 dec 17, code provided by Kevin Blagrave and described in
 *>>refer	Fe3	model	Blagrave, K.P.M., Martin, P.G. & Baldwin, J.A. 
 *>>refercon 2006, ApJ, in press (astro-ph 0603167) */
STATIC void Fe3Lev14(void)
{
	bool lgZeroPop;
	int nNegPop;
	long int i,
		j;
	static bool lgFirst=true;

	double dfe3dt;

	long int ihi , ilo;
	static double
		*depart,
		*pops,
		*destroy,
		*create ,
		**AulDest,
		**CollRate,
		**AulPump,
		**AulNet,
		**col_str;

	/* statistical weights for the n levels */
	static double gfe3[NLFE3]={9.,7.,5.,3.,1.,5.,13.,11.,9.,3.,1.,9.,7.,5.};

	/*refer Fe3	energies	NIST version 3 Atomic Spectra Database */
	/* from smallest to largest energy (not in multiplet groupings) */

	/* energy in wavenumbers */
	static double excit_wn[NLFE3]={
		0.0    ,   436.2,   738.9,   932.4,  1027.3,
		19404.8, 20051.1, 20300.8, 20481.9, 20688.4,
		21208.5, 21462.2, 21699.9, 21857.2 };

	DEBUG_ENTRY( "Fe3Lev14()" );

	if( lgFirst )
	{
		/* will never do this again */
		lgFirst = false;
		/* allocate the 1D arrays*/
		/* create space for the 2D arrays */
		depart = ((double *)MALLOC((NLFE3)*sizeof(double)));
		pops = ((double *)MALLOC((NLFE3)*sizeof(double)));
		destroy = ((double *)MALLOC((NLFE3)*sizeof(double)));
		create = ((double *)MALLOC((NLFE3)*sizeof(double)));
		/* now the 2-d arrays */
		fe.Fe3_wl = ((double **)MALLOC((NLFE3)*sizeof(double *)));
		fe.Fe3_emiss = ((double **)MALLOC((NLFE3)*sizeof(double *)));
		AulNet = ((double **)MALLOC((NLFE3)*sizeof(double *)));
		col_str = ((double **)MALLOC((NLFE3)*sizeof(double *)));
		AulPump = ((double **)MALLOC((NLFE3)*sizeof(double *)));
		CollRate = ((double **)MALLOC((NLFE3)*sizeof(double *)));
		AulDest = ((double **)MALLOC((NLFE3)*sizeof(double *)));
		for( i=0; i < NLFE3; ++i )
		{
			fe.Fe3_wl[i] = ((double *)MALLOC((NLFE3)*sizeof(double )));
			fe.Fe3_emiss[i] = ((double *)MALLOC((NLFE3)*sizeof(double )));
			AulNet[i] = ((double *)MALLOC((NLFE3)*sizeof(double )));
			col_str[i] = ((double *)MALLOC((NLFE3)*sizeof(double )));
			AulPump[i] = ((double *)MALLOC((NLFE3)*sizeof(double )));
			CollRate[i] = ((double *)MALLOC((NLFE3)*sizeof(double )));
			AulDest[i] = ((double *)MALLOC((NLFE3)*sizeof(double )));
		}

		/* set some to constant values after zeroing out */
		for( i=0; i < NLFE3; ++i )
		{
			create[i] = 0.;
			destroy[i] = 0.;
			for( j=0; j < NLFE3; ++j )
			{
				AulNet[i][j] = 0.;
				col_str[i][j] = 0.;
				CollRate[i][j] = 0.;
				AulDest[i][j] = 0.;
				AulPump[i][j] = 0.;
				fe.Fe3_wl[i][j] = 0.;
				fe.Fe3_emiss[i][j] = 0.;
			}
		}
		/* calculates wavelengths of transitions */
		/* dividing by RefIndex converts the vacuum wavelength to air wavelength */
		for( ihi=1; ihi < NLFE3; ++ihi )
		{
			for( ilo=0; ilo < ihi; ++ilo )
			{
				fe.Fe3_wl[ihi][ilo] = 1e8/(excit_wn[ihi]-excit_wn[ilo]) / 
					RefIndex( (excit_wn[ihi]-excit_wn[ilo]) );
			}
		}

		/* assume FeIII is optically thin - just use As as net escape */
		/*>>refer	Fe3	as	Quinet, P., 1996, A&AS, 116, 573      */
		AulNet[1][0] = 2.8e-3;
		AulNet[7][0] = 4.9e-6;
		AulNet[8][0] = 5.7e-3;
		AulNet[11][0] = 4.5e-1;
		AulNet[12][0] = 4.2e-2;

		AulNet[2][1] = 1.8e-3;
		AulNet[5][1] = 4.2e-1;
		AulNet[8][1] = 1.0e-3;
		AulNet[11][1] = 8.4e-2;
		AulNet[12][1] = 2.5e-1;
		AulNet[13][1] = 2.7e-2;

		AulNet[3][2] = 7.0e-4;
		AulNet[5][2] = 5.1e-5;
		AulNet[9][2] = 5.4e-1;
		AulNet[12][2] = 8.5e-2;
		AulNet[13][2] = 9.8e-2;

		AulNet[4][3] = 1.4e-4;
		AulNet[5][3] = 3.9e-2;
		AulNet[9][3] = 4.1e-5;
		AulNet[10][3] = 7.0e-1;
		AulNet[13][3] = 4.7e-2;

		AulNet[9][4] = 9.3e-2;

		AulNet[9][5] = 4.7e-2;
		AulNet[12][5] = 2.5e-6;
		AulNet[13][5] = 1.7e-5;

		AulNet[7][6] = 2.7e-4;

		AulNet[8][7] = 1.2e-4;
		AulNet[11][7] = 6.6e-4;

		AulNet[11][8] = 1.6e-3;
		AulNet[12][8] = 7.8e-4;

		AulNet[10][9] = 8.4e-3;
		AulNet[13][9] = 2.8e-7;

		AulNet[12][11] = 3.0e-4;

		AulNet[13][12] = 1.4e-4;

		/* collision strengths at log T 4 */
		/** \todo	2	put in temperature dependence */
		/*>>refer	Fe3	cs	Zhang, H.  1996, 119, 523 */
		col_str[1][0] = 2.92;
		col_str[2][0] = 1.24;
		col_str[3][0] = 0.595;
		col_str[4][0] = 0.180;
		col_str[5][0] = 0.580;
		col_str[6][0] = 1.34; 
		col_str[7][0] = 0.489;
		col_str[8][0] = 0.0926;
		col_str[9][0] = 0.165;
		col_str[10][0] = 0.0213;
		col_str[11][0] = 1.07;
		col_str[12][0] = 0.435;
		col_str[13][0] = 0.157;

		col_str[2][1] = 2.06;
		col_str[3][1] = 0.799;
		col_str[4][1] = 0.225;
		col_str[5][1] = 0.335;
		col_str[6][1] = 0.555;
		col_str[7][1] = 0.609; 
		col_str[8][1] = 0.367;
		col_str[9][1] = 0.195;
		col_str[10][1] = 0.0698;
		col_str[11][1] = 0.538;
		col_str[12][1] = 0.484;
		col_str[13][1] = 0.285;

		col_str[3][2] = 1.29;
		col_str[4][2] = 0.312;
		col_str[5][2] = 0.173;
		col_str[6][2] = 0.178;
		col_str[7][2] = 0.430;
		col_str[8][2] = 0.486;
		col_str[9][2] = 0.179;
		col_str[10][2] = 0.0741;
		col_str[11][2] = 0.249;
		col_str[12][2] = 0.362;
		col_str[13][2] = 0.324;

		col_str[4][3] = 0.493;
		col_str[5][3] = 0.0767;
		col_str[6][3] = 0.0348;
		col_str[7][3] = 0.223;
		col_str[8][3] = 0.401;
		col_str[9][3] = 0.126;
		col_str[10][3] = 0.0528;
		col_str[11][3] = 0.101;
		col_str[12][3] = 0.207;
		col_str[13][3] = 0.253;

		col_str[5][4] = 0.0211;
		col_str[6][4] = 0.00122; 
		col_str[7][4] = 0.0653;
		col_str[8][4] = 0.154;
		col_str[9][4] = 0.0453;
		col_str[10][4] = 0.0189;
		col_str[11][4] = 0.0265;
		col_str[12][4] = 0.0654;
		col_str[13][4] = 0.0950;

		col_str[6][5] = 0.403;
		col_str[7][5] = 0.213;
		col_str[8][5] = 0.0939;
		col_str[9][5] = 1.10;
		col_str[10][5] = 0.282;
		col_str[11][5] = 0.942;
		col_str[12][5] = 0.768;
		col_str[13][5] = 0.579;

		col_str[7][6] = 2.84; /* 10-9 */
		col_str[8][6] = 0.379; /* 11-9 */
		col_str[9][6] = 0.0876;  /* 7-9 */
		col_str[10][6] = 0.00807; /* 8-9 */
		col_str[11][6] = 1.85; /* 12-9 */
		col_str[12][6] = 0.667; /* 13-9 */
		col_str[13][6] = 0.0905; /* 14-9 */

		col_str[8][7] = 3.07; /* 11-10 */
		col_str[9][7] = 0.167;   /* 7-10 */
		col_str[10][7] = 0.0526;  /* 8-10 */
		col_str[11][7] = 0.814; /* 12-10 */
		col_str[12][7] = 0.837; /* 13-10 */
		col_str[13][7] = 0.626; /* 14-10 */

		col_str[9][8] = 0.181; /* 7-11 */
		col_str[10][8] = 0.0854; /* 8-11 */
		col_str[11][8] = 0.180; /* 12-11 */
		col_str[12][8] = 0.778; /* 13-11 */
		col_str[13][8] = 0.941; /* 14-11 */

		col_str[10][9] = 0.377; /* 8-7 */
		col_str[11][9] = 0.603; /* 12-7 */
		col_str[12][9] = 0.472; /* 13-7 */
		col_str[13][9] = 0.302; /* 14-7 */

		col_str[11][10] = 0.216; /* 12-8 */
		col_str[12][10] = 0.137; /* 13-8 */
		col_str[13][10] = 0.106; /* 14-8 */

		col_str[12][11] = 1.25;
		col_str[13][11] = 0.292;

		col_str[13][12] = 1.10;
	}

	/* bail if no ions */
	if( dense.xIonDense[ipIRON][2] <= 0. )
	{
		CoolAdd("Fe 3",0,0.);

		fe.Fe3CoolTot = 0.;   
		for( ihi=1; ihi < NLFE3; ++ihi )
		{
			for( ilo=0; ilo < ihi; ++ilo )
			{
				fe.Fe3_emiss[ihi][ilo] = 0.;
			}
		}
		/* level populations */
		/* LIMLEVELN is the dimension of the atoms vectors */
		ASSERT( NLFE3 <= LIMLEVELN);
		for( i=0; i < NLFE3; i++ )
		{
			atoms.PopLevels[i] = 0.;
			atoms.DepLTELevels[i] = 1.;
		}
		return;
	}

	/* nNegPop positive if negative pops occurred, negative if too cold */
	atom_levelN(
		/* number of levels */
		NLFE3,
		/* the abundance of the ion, cm-3 */
		dense.xIonDense[ipIRON][2],
		/* the statistical weights */
		gfe3,
		/* the excitation energies */
		excit_wn,
		'w',
		/* the derived populations - cm-3 */
		pops,
		/* the derived departure coefficients */
		depart,
		/* the net emission rate, Aul * escape prob */
		&AulNet ,
		/* the collision strengths */
		&col_str ,
		/* A * destruction prob */
		&AulDest,
		/* pumping rate */
		&AulPump,
		/* collision rate, s-1, must defined if no collision strengths */
		&CollRate,
		/* creation vector */
		create,
		/* destruction vector */
		destroy,
		/* say atom_levelN should evaluate coll rates from cs */
		false,   
		&fe.Fe3CoolTot,
		&dfe3dt,
		"Fe 3",
		&nNegPop,
		&lgZeroPop,
		false );

	/* LIMLEVELN is the dim of the PopLevels vector */
	ASSERT( NLFE3 <= LIMLEVELN );
	for( i=0; i < NLFE3; ++i )
	{
		atoms.PopLevels[i] = pops[i];
		atoms.DepLTELevels[i] = depart[i];
	}

	if( nNegPop > 0 )
	{
		fprintf( ioQQQ, " Fe3Lev14 found negative populations\n" );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	/* add cooling then its derivative */
	CoolAdd("Fe 3",0,fe.Fe3CoolTot);
	/* derivative of cooling */
	thermal.dCooldT += dfe3dt;

	/* find emission in each line */
	for( ihi=1; ihi < NLFE3; ++ihi )
	{
		for( ilo=0; ilo < ihi; ++ilo )
		{
			/* emission in these lines */
			fe.Fe3_emiss[ihi][ilo] = pops[ihi]*AulNet[ihi][ilo]*(excit_wn[ihi] - excit_wn[ilo])*T1CM*BOLTZMANN;
		}
	}
	return;
}

/*Fe11Lev5 compute populations and cooling due to 5 level Fe 11 ion */
STATIC void Fe11Lev5(void)
{
	bool lgZeroPop;
	int nNegPop;
	long int i,
		j;
	static bool lgFirst=true;

	double dCool_dT;

	long int ihi , ilo;
	static double
		*depart,
		*pops,
		*destroy,
		*create ,
		**AulDest,
		**CollRate,
		**AulPump,
		**AulNet,
		**col_str ,
		TeUsed=-1.;

	/* statistical weights for the n levels */
	static double stat_wght[NLFE11]={5.,3.,1.,5.,1.};

	/*refer Fe11	energies	NIST version 3 Atomic Spectra Database */
	/* from smallest to largest energy (not in multiplet groupings) */

	/* energy in wavenumbers */
	static double excit_wn[NLFE11]={
		0.0 , 12667.9 , 14312. , 37743.6 , 80814.7 };

	DEBUG_ENTRY( "Fe11Lev5()" );

	if( lgFirst )
	{
		/* will never do this again */
		lgFirst = false;
		/* allocate the 1D arrays*/
		/* create space for the 2D arrays */
		depart = ((double *)MALLOC((NLFE11)*sizeof(double)));
		pops = ((double *)MALLOC((NLFE11)*sizeof(double)));
		destroy = ((double *)MALLOC((NLFE11)*sizeof(double)));
		create = ((double *)MALLOC((NLFE11)*sizeof(double)));
		/* now the 2-d arrays */
		fe.Fe11_wl = ((double **)MALLOC((NLFE11)*sizeof(double *)));
		fe.Fe11_emiss = ((double **)MALLOC((NLFE11)*sizeof(double *)));
		AulNet = ((double **)MALLOC((NLFE11)*sizeof(double *)));
		col_str = ((double **)MALLOC((NLFE11)*sizeof(double *)));
		AulPump = ((double **)MALLOC((NLFE11)*sizeof(double *)));
		CollRate = ((double **)MALLOC((NLFE11)*sizeof(double *)));
		AulDest = ((double **)MALLOC((NLFE11)*sizeof(double *)));
		for( i=0; i < NLFE11; ++i )
		{
			fe.Fe11_wl[i] = ((double *)MALLOC((NLFE11)*sizeof(double )));
			fe.Fe11_emiss[i] = ((double *)MALLOC((NLFE11)*sizeof(double )));
			AulNet[i] = ((double *)MALLOC((NLFE11)*sizeof(double )));
			col_str[i] = ((double *)MALLOC((NLFE11)*sizeof(double )));
			AulPump[i] = ((double *)MALLOC((NLFE11)*sizeof(double )));
			CollRate[i] = ((double *)MALLOC((NLFE11)*sizeof(double )));
			AulDest[i] = ((double *)MALLOC((NLFE11)*sizeof(double )));
		}

		/* set some to constant values after zeroing out */
		for( i=0; i < NLFE11; ++i )
		{
			create[i] = 0.;
			destroy[i] = 0.;
			for( j=0; j < NLFE11; ++j )
			{
				AulNet[i][j] = 0.;
				col_str[i][j] = 0.;
				CollRate[i][j] = 0.;
				AulDest[i][j] = 0.;
				AulPump[i][j] = 0.;
				fe.Fe11_wl[i][j] = 0.;
				fe.Fe11_emiss[i][j] = 0.;
			}
		}
		/* calculates wavelengths of transitions */
		/* dividing by RefIndex converts the vacuum wavelength to air wavelength */
		for( ihi=1; ihi < NLFE11; ++ihi )
		{
			for( ilo=0; ilo < ihi; ++ilo )
			{
				fe.Fe11_wl[ihi][ilo] = 1e8/(excit_wn[ihi]-excit_wn[ilo]) / 
					RefIndex( (excit_wn[ihi]-excit_wn[ilo]) );
			}
		}

		/*>>refer	Fe11	as	Mendoza C. & Zeippen, C. J. 1983, MNRAS, 202, 981 */
		AulNet[4][0] = 1.7;
		AulNet[4][1] = 990.;
		AulNet[4][3] = 8.3;
		AulNet[3][0] = 92.3;
		AulNet[3][1] = 9.44;
		AulNet[3][2] = 1.4e-3;
		AulNet[2][0] = 9.9e-3;
		AulNet[1][0] = 43.7;
		AulNet[2][1] = 0.226;

	}

	/* bail if no ions */
	if( dense.xIonDense[ipIRON][10] <= 0. )
	{
		CoolAdd("Fe11",0,0.);

		fe.Fe11CoolTot = 0.;   
		for( ihi=1; ihi < NLFE11; ++ihi )
		{
			for( ilo=0; ilo < ihi; ++ilo )
			{
				fe.Fe11_emiss[ihi][ilo] = 0.;
			}
		}
		/* level populations */
		/* LIMLEVELN is the dimension of the atoms vectors */
		ASSERT( NLFE11 <= LIMLEVELN);
		for( i=0; i < NLFE11; i++ )
		{
			atoms.PopLevels[i] = 0.;
			atoms.DepLTELevels[i] = 1.;
		}
		return;
	}

	/* evaluate collision strengths if Te changed */
	if( !fp_equal(phycon.te,TeUsed) )
	{
		TeUsed = phycon.te;

		/* collision strengths at current temperature */
		for( ihi=1; ihi < NLFE11; ++ihi )
		{
			for( ilo=0; ilo < ihi; ++ilo )
			{
				col_str[ihi][ilo] = Fe_10_11_13_cs( 11 , ilo+1 , ihi+1 );
				ASSERT( col_str[ihi][ilo] > 0. );
			}
		}
	}

	/* nNegPop positive if negative pops occurred, negative if too cold */
	atom_levelN(
		/* number of levels */
		NLFE11,
		/* the abundance of the ion, cm-3 */
		dense.xIonDense[ipIRON][10],
		/* the statistical weights */
		stat_wght,
		/* the excitation energies */
		excit_wn,
		'w',
		/* the derived populations - cm-3 */
		pops,
		/* the derived departure coefficients */
		depart,
		/* the net emission rate, Aul * escp prob */
		&AulNet ,
		/* the collision strengths */
		&col_str ,
		/* A * destruction prob */
		&AulDest,
		/* pumping rate */
		&AulPump,
		/* collision rate, s-1, must defined if no collision strengths */
		&CollRate,
		/* creation vector */
		create,
		/* destruction vector */
		destroy,
		/* say atom_levelN should evaluate coll rates from cs */
		false,   
		&fe.Fe11CoolTot,
		&dCool_dT,
		"Fe11",
		&nNegPop,
		&lgZeroPop,
		false );

	/* LIMLEVELN is the dim of the PopLevels vector */
	ASSERT( NLFE11 <= LIMLEVELN );
	for( i=0; i < NLFE11; ++i )  
	{
		atoms.PopLevels[i] = pops[i];
		atoms.DepLTELevels[i] = depart[i];
	}

	if( nNegPop > 0 )  
	{
		fprintf( ioQQQ, " Fe11Lev5 found negative populations\n" );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	/* add cooling then its derivative */
	CoolAdd("Fe11",0,fe.Fe11CoolTot);
	/* derivative of cooling */
	thermal.dCooldT += dCool_dT;

	/* find emission in each line */
	for( ihi=1; ihi < NLFE11; ++ihi )
	{
		for( ilo=0; ilo < ihi; ++ilo )
		{
			/* emission in these lines */
			fe.Fe11_emiss[ihi][ilo] = pops[ihi]*AulNet[ihi][ilo]*
				(excit_wn[ihi] - excit_wn[ilo])*T1CM*BOLTZMANN;
		}
	}
	return;
}

/*Fe13Lev5 compute populations and cooling due to 5 level Fe 13 ion */
STATIC void Fe13Lev5(void)
{
	bool lgZeroPop;
	int nNegPop;
	long int i,
		j;
	static bool lgFirst=true;

	double dCool_dT;

	long int ihi , ilo;
	static double
		*depart,
		*pops,
		*destroy,
		*create ,
		**AulDest,
		**CollRate,
		**AulPump,
		**AulNet,
		**col_str ,
		TeUsed=-1.;

	/* statistical weights for the n levels */
	static double stat_wght[NLFE13]={1.,3.,5.,5.,1.};

	/*refer Fe13	energies	NIST version 3 Atomic Spectra Database */
	/* energy in wavenumbers */
	static double excit_wn[NLFE13]={
		0.0 , 9302.5 , 18561.0 , 48068. , 91508. };

	DEBUG_ENTRY( "Fe13Lev5()" );

	if( lgFirst )
	{
		/* will never do this again */
		lgFirst = false;
		/* allocate the 1D arrays*/
		/* create space for the 2D arrays */
		depart = ((double *)MALLOC((NLFE13)*sizeof(double)));
		pops = ((double *)MALLOC((NLFE13)*sizeof(double)));
		destroy = ((double *)MALLOC((NLFE13)*sizeof(double)));
		create = ((double *)MALLOC((NLFE13)*sizeof(double)));
		/* now the 2-d arrays */
		fe.Fe13_wl = ((double **)MALLOC((NLFE13)*sizeof(double *)));
		fe.Fe13_emiss = ((double **)MALLOC((NLFE13)*sizeof(double *)));
		AulNet = ((double **)MALLOC((NLFE13)*sizeof(double *)));
		col_str = ((double **)MALLOC((NLFE13)*sizeof(double *)));
		AulPump = ((double **)MALLOC((NLFE13)*sizeof(double *)));
		CollRate = ((double **)MALLOC((NLFE13)*sizeof(double *)));
		AulDest = ((double **)MALLOC((NLFE13)*sizeof(double *)));
		for( i=0; i < NLFE13; ++i )
		{
			fe.Fe13_wl[i] = ((double *)MALLOC((NLFE13)*sizeof(double )));
			fe.Fe13_emiss[i] = ((double *)MALLOC((NLFE13)*sizeof(double )));
			AulNet[i] = ((double *)MALLOC((NLFE13)*sizeof(double )));
			col_str[i] = ((double *)MALLOC((NLFE13)*sizeof(double )));
			AulPump[i] = ((double *)MALLOC((NLFE13)*sizeof(double )));
			CollRate[i] = ((double *)MALLOC((NLFE13)*sizeof(double )));
			AulDest[i] = ((double *)MALLOC((NLFE13)*sizeof(double )));
		}

		/* set some to constant values after zeroing out */
		for( i=0; i < NLFE13; ++i )
		{
			create[i] = 0.;
			destroy[i] = 0.;
			for( j=0; j < NLFE13; ++j )
			{
				AulNet[i][j] = 0.;
				col_str[i][j] = 0.;
				CollRate[i][j] = 0.;
				AulDest[i][j] = 0.;
				AulPump[i][j] = 0.;
				fe.Fe13_wl[i][j] = 0.;
				fe.Fe13_emiss[i][j] = 0.;
			}
		}
		/* calculates wavelengths of transitions */
		/* dividing by RefIndex converts the vacuum wavelength to air wavelength */
		for( ihi=1; ihi < NLFE13; ++ihi )
		{
			for( ilo=0; ilo < ihi; ++ilo )
			{
				fe.Fe13_wl[ihi][ilo] = 1e8/(excit_wn[ihi]-excit_wn[ilo]) / 
					RefIndex( (excit_wn[ihi]-excit_wn[ilo]) );
			}
		}

		/*>>refer	Fe13	as	Huang, K.-N. 1985, At. Data Nucl. Data Tables 32, 503 */
		AulNet[4][1] = 1010.;
		AulNet[4][2] = 3.8;
		AulNet[4][3] = 8.1;
		AulNet[3][1] = 45.7;
		AulNet[3][2] = 57.6;
		AulNet[2][0] = 0.0063;
		AulNet[1][0] = 14.0;
		AulNet[2][1] = 9.87;

	}

	/* bail if no ions */
	if( dense.xIonDense[ipIRON][12] <= 0. )
	{
		CoolAdd("Fe13",0,0.);

		fe.Fe13CoolTot = 0.;   
		for( ihi=1; ihi < NLFE13; ++ihi )
		{
			for( ilo=0; ilo < ihi; ++ilo )
			{
				fe.Fe13_emiss[ihi][ilo] = 0.;
			}
		}
		/* level populations */
		/* LIMLEVELN is the dimension of the atoms vectors */
		ASSERT( NLFE13 <= LIMLEVELN);
		for( i=0; i < NLFE13; i++ )
		{
			atoms.PopLevels[i] = 0.;
			atoms.DepLTELevels[i] = 1.;
		}
		return;
	}

	/* evaluate collision strengths if Te changed */
	if( !fp_equal(phycon.te,TeUsed) )
	{
		TeUsed = phycon.te;

		/* collision strengths at current temperature */
		for( ihi=1; ihi < NLFE13; ++ihi )
		{
			for( ilo=0; ilo < ihi; ++ilo )
			{
				col_str[ihi][ilo] = Fe_10_11_13_cs( 13 , ilo+1 , ihi+1 );
				ASSERT( col_str[ihi][ilo] > 0. );
			}
		}
	}

	/*fprintf(ioQQQ,"DEBUG Fe13 N %.2e %.2e %.2e %.2e %.2e %.2e\n",
		col_str[1][0] , col_str[2][1] , col_str[2][0] ,
		AulNet[1][0] , AulNet[2][1] , AulNet[2][0]);*/

	/* nNegPop positive if negative pops occurred, negative if too cold */
	atom_levelN(
		/* number of levels */
		NLFE13,
		/* the abundance of the ion, cm-3 */
		dense.xIonDense[ipIRON][12],
		/* the statistical weights */
		stat_wght,
		/* the excitation energies */
		excit_wn,
		'w',
		/* the derived populations - cm-3 */
		pops,
		/* the derived departure coefficients */
		depart,
		/* the net emission rate, Aul * escape prob */
		&AulNet ,
		/* the collision strengths */
		&col_str ,
		/* A * destruction prob */
		&AulDest,
		/* pumping rate */
		&AulPump,
		/* collision rate, s-1, must defined if no collision strengths */
		&CollRate,
		/* creation vector */
		create,
		/* destruction vector */
		destroy,
		/* say atom_levelN should evaluate coll rates from cs */
		false,   
		&fe.Fe13CoolTot,
		&dCool_dT,
		"Fe13",
		&nNegPop,
		&lgZeroPop,
		false );

	/* LIMLEVELN is the dim of the PopLevels vector */
	ASSERT( NLFE13 <= LIMLEVELN );
	for( i=0; i < NLFE13; ++i )  
	{
		atoms.PopLevels[i] = pops[i];
		atoms.DepLTELevels[i] = depart[i];
	}

	if( nNegPop > 0 )  
	{
		fprintf( ioQQQ, " Fe13Lev5 found negative populations\n" );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	/* add cooling then its derivative */
	CoolAdd("Fe13",0,fe.Fe13CoolTot);
	/* derivative of cooling */
	thermal.dCooldT += dCool_dT;

	/* find emission in each line */
	for( ihi=1; ihi < NLFE13; ++ihi )
	{
		for( ilo=0; ilo < ihi; ++ilo )
		{
			/* emission in these lines */
			fe.Fe13_emiss[ihi][ilo] = pops[ihi]*AulNet[ihi][ilo]*(excit_wn[ihi] - excit_wn[ilo])*T1CM*BOLTZMANN;
		}
	}
	return;
}
