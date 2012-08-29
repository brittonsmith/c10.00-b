/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*SaveHeat save contributors to local heating, with save heat command, called by punch_do */
#include "cddefines.h"
#include "thermal.h"
#include "radius.h"
#include "conv.h"
#include "lines_service.h"
#include "dense.h"
#include "taulines.h"
#include "phycon.h"
#include "elementnames.h"
#include "dynamics.h"
#include "save.h"

/* limit for number of heat agents that are saved */
/* limit to number to print */
static const int IPRINT= 9;

/*SaveHeat save contributors to local heating, with save heat command, 
 * called by punch_do */
void SaveHeat(FILE* io)
{
	char **chLabel, 
	  chLbl[11];
	bool lgHeatLine;
	int nFail;
	long int i, 
	  ipStrong, 
	  ipnt, 
	  *ipOrdered,
	  *ipsave, 
	  j, 
	  *jpsave, 
	  k, 
	  level;
	double CS, 
	  ColHeat, 
	  EscP, 
	  Pump, 
	  Strong, 
	  TauIn,
		cool_total,
		heat_total;
	realnum *SaveVal;

	DEBUG_ENTRY( "SaveHeat()" );

	SaveVal = (realnum *) CALLOC(LIMELM*LIMELM, sizeof(realnum));
	ipsave = (long int *) CALLOC(LIMELM*LIMELM, sizeof(long int));
	jpsave = (long int *) CALLOC(LIMELM*LIMELM, sizeof(long int));
	ipOrdered = (long int *) CALLOC(LIMELM*LIMELM, sizeof(long int));
	chLabel = (char **) CALLOC(LIMELM*LIMELM, sizeof(char *));

	for( i=0; i<LIMELM*LIMELM; ++i )
	{
		ipsave[i] = INT_MIN;
		jpsave[i] = INT_MIN;
		SaveVal[i] = -FLT_MAX;
		chLabel[i] = (char *) CALLOC(10, sizeof(char));
	}

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
	if(dynamics.Cool() > dynamics.Heat()) 
	{
		cool_total -= dynamics.Heat();
		heat_total -= dynamics.Heat();
	} 
	else
	{
		cool_total -= dynamics.Cool();
		heat_total -= dynamics.Cool();
	}
#	endif

	ipnt = 0;

	/* heat sources are saved in a 2d square array */
	/* WeakHeatCool set with 'set weakheatcool' command
	 * default is 0.05 */
	for( i=0; i < LIMELM; i++ )
	{
		for( j=0; j < LIMELM; j++ )
		{
			if( thermal.heating[i][j]/SDIV(heat_total) > SMALLFLOAT )
			{
				ipsave[ipnt] = i;
				jpsave[ipnt] = j;
				/*
				  8-29-2012 - Britton changed this so heating rates
				  are saved instead of fractions of the total.
				SaveVal[ipnt] = (realnum)(thermal.heating[i][j]/SDIV(heat_total));
				*/
				SaveVal[ipnt] = (realnum)(thermal.heating[i][j]);
				ipnt++;
			}
		}
	}

	/* now check for possible line heating not counted in 1,23
	 * there should not be any significant heating source here
	 * since they would not be counted in derivative correctly */
	for( i=0; i < thermal.ncltot; i++ )
	{
		if( thermal.heatnt[i]/SDIV(heat_total) > save.WeakHeatCool )
		{
			realnum awl;
			awl = thermal.collam[i];
			/* force to save wavelength convention as printout */
			if( awl > 100000 )
				awl /= 10000;
			fprintf( io, " Negative coolant was %s %.2f %.2e\n", 
			  thermal.chClntLab[i], awl, thermal.heatnt[i]/SDIV(heat_total) );
		}
	}

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

	/* this is mainly to keep the compiler from flagging possible paths 
	 * with j not being set */
	i = INT_MIN;
	j = INT_MIN;
	/* following loop tries to identify strongest agents and turn to labels */
	for( k=0; k < ipnt; k++ )
	{
		/* generate labels that make sense in printout 
		 * if not identified with a specific agent, will print indices as [i][j],
		 * heating is thermal.heating[i][j] */
		i = ipsave[k];
		j = jpsave[k];
		/* i >= j indicates agent is one of the first LIMELM elements */
		if( i >= j )
		{
			if( dense.xIonDense[i][j] == 0. && thermal.heating[i][j]>SMALLFLOAT )
				fprintf(ioQQQ,"DISASTER assert about to be thrown - search for hit it\n");
			/*fprintf(ioQQQ,"DEBUG %li %li %.2e %.2e\n", i , j , 
				dense.xIonDense[i][j],
				thermal.heating[i][j]);*/
			ASSERT( dense.xIonDense[i][j] > 0. || thermal.heating[i][j]<SMALLFLOAT );
			/* this is case of photoionization of atom or ion */
			strcpy( chLabel[k], elementnames.chElementSym[i] );
			strcat( chLabel[k], elementnames.chIonStage[j]  );
		}
		/* notice that in test i and j are swapped from order in heating array */
		else if( i == 0 && j == 1 )
		{
			/* photoionization from all excited states of Hydrogenic species,
			 * heating[0][1] */
			strcpy( chLabel[k], "Hn=2" );
		}
		else if( i == 0 && j == 3 )
		{
			/* collisional ionization of all iso-seq from all levels, 
			 * heating[0][3] */
			strcpy( chLabel[k], "Hion" );
		}
		else if( i == 0 && j == 7 )
		{
			/* UTA ionization heating[0][7] */
			strcpy( chLabel[k], " UTA" );
		}
		else if( i == 0 && j == 8 )
		{
			/* thermal.heating[0][8] is heating due to collisions within 
			 * X of H2, code var hmi.HeatH2Dexc_used, hmi.HeatH2Dexc_BigH2,
			 * hmi.HeatH2Dexc_TH85 */
			strcpy( chLabel[k], "H2vH" );
		}
		else if( i == 0 && j == 17 )
		{
			/* thermal.heating[0][17] is heating due to photodissociation 
			 * heating of X within H2,
			 * code var hmi.HeatH2Dish_used */
			strcpy( chLabel[k], "H2dH" );
		}
		else if( i == 0 && j == 9 )
		{
			/* CO dissociation, co.CODissHeat, heating[0][9] */
			strcpy( chLabel[k], "COds" );
		}
		else if( i == 0 && j == 20 )
		{
			/* extra heat thermal.heating[0][20]*/
			strcpy( chLabel[k], "extH" );
		}
		else if( i == 0 && j == 21 )
		{
			/* pair heating thermal.heating[0][21]*/
			strcpy( chLabel[k], "pair" );
		}
		else if( i == 0 && j == 11 )
		{
			/* free free heating, heating[0][11] */
			strcpy( chLabel[k], "H FF" );
		}
		else if( i == 0 && j == 12 )
		{
			/* heating coolant (not line), physical cooling process, often a bug, heating[0][12] */
			strcpy( chLabel[k], "Hcol" );
		}
		else if( i == 0 && j == 13 )
		{
			/* grain photoionization, heating[0][13] */
			strcpy( chLabel[k], "GrnP" );
		}
		else if( i == 0 && j == 14 )
		{
			/* grain collisions, heating[0][14] */
			strcpy( chLabel[k], "GrnC" );
		}
		else if( i == 0 && j == 15 )
		{
			/* H- heating, heating[0][15] */
			strcpy( chLabel[k], "H-  " );
		}
		else if( i == 0 && j == 16 )
		{
			/* H2+ heating, heating[0][16] */
			strcpy( chLabel[k], "H2+ " );
		}
		else if( i == 0 && j == 18 )
		{
			/* H2 photoionization heating, heating[0][18] */
			strcpy( chLabel[k], "H2ph" );
		}
		else if( i == 0 && j == 19 )
		{
			/* Compton heating, heating[0][19] */
			strcpy( chLabel[k], "Comp" );
		}
		else if( i == 0 && j == 22 )
		{
			/* line heating, heating[0][22] */
			strcpy( chLabel[k], "line" );
		}
		else if( i == 0 && j == 23 )
		{
			/* iso-sequence line heating - all elements together, 
			 * heating[0][23] */
			strcpy( chLabel[k], "Hlin" );
		}
		else if( i == 0 && j == 24 )
		{
			/* charge transfer heating, heating[0][24] */
			strcpy( chLabel[k], "ChaT" );
		}
		else if( i == 1 && j == 3 )
		{
			/* helium triplet line heating, heating[1][3] */
			strcpy( chLabel[k], "He3l" );
		}
		else if( i == 1 && j == 5 )
		{
			/* advective heating, heating[1][5] */
			strcpy( chLabel[k], "adve" );
		}
		else if( i == 1 && j == 6 )
		{
			/* cosmic ray heating thermal.heating[1][6]*/
			strcpy( chLabel[k], "CR H" );
		}
		else if( i == 25 && j == 27 )
		{
			/* Fe 2 line heating, heating[25][27] */
			strcpy( chLabel[k], "Fe 2" );
		}
		else
		{
			sprintf( chLabel[k], "[%ld][%ld]" , i , j );
		}
	}

	/* now sort by decreasing importance */
	/*spsort netlib routine to sort array returning sorted indices */
	spsort(
		  /* input array to be sorted */
		  SaveVal, 
		  /* number of values in x */
		  ipnt, 
		  /* permutation output array */
		  ipOrdered, 
		  /* flag saying what to do - 1 sorts into increasing order, not changing
		   * the original routine */
		  -1, 
		  /* error condition, should be 0 */
		  &nFail);

	/*>>chng 06 jun 06, change start of save to give same info as cooling 
	 * as per comment by Yumihiko Tsuzuki */
	/* begin the print out with zone number, total heating and cooling */
	/*
	  8-29-2012 - Britton changed the output precision to six digits.
	fprintf( io, "%.5e\t%.4e\t%.4e\t%.4e", 
	 */
	fprintf( io, "%.5e\t%.4e\t%.6e\t%.6e", 
		radius.depth_mid_zone, 
		phycon.te, 
		heat_total, 
		cool_total );

	/* we only want to print the IPRINT strongest of the group */
	ipnt = MIN2( ipnt , IPRINT );

	for( k=0; k < ipnt; k++ )
	{
		int ip = ipOrdered[k];
		i = ipsave[ip];
		j = jpsave[ip];
		ASSERT( i<LIMELM && j<LIMELM );
		if(k > 4 && thermal.heating[i][j]/SDIV(heat_total) < save.WeakHeatCool )
			break;
		/*
		  8-29-2012 - Britton changed this so heating rates are written out 
		  instead of fractions of the total heating rate.
		fprintf( io, "\t%s\t%.7f ", 
		*/
		fprintf( io, "\t%s\t%.6e ", 
			chLabel[ip], SaveVal[ip] );
	}
	fprintf( io, " \n" );

	/* a negative pointer in the heating array is probably a problem,
	 * indicating that some line has become a heat source */
	lgHeatLine = false;

	/* check if any lines were major heat sources */
	for( i=0; i < ipnt; i++ )
	{
		/* heating[22][0] is line heating - identify line if important */
		if( ipsave[i] == 0 && jpsave[i] == 22 )
			lgHeatLine = true;
	}

	if( lgHeatLine )
	{
		/* a line was a major heat source - identify it */
		FndLineHt(&level,&ipStrong,&Strong);
		if( Strong/SDIV(heat_total) > 0.005 )
		{
			if( level == 1 )
			{
				strcpy( chLbl, chLineLbl(&TauLines[ipStrong] ) );
				TauIn = TauLines[ipStrong].Emis->TauIn;
				Pump = TauLines[ipStrong].Emis->pump;
				EscP = TauLines[ipStrong].Emis->Pesc;
				CS = TauLines[ipStrong].Coll.col_str;
				/* ratio of line to total heating */
				ColHeat = TauLines[ipStrong].Coll.heat/
				  SDIV(heat_total);
			}
			else if( level == 2 )
			{
				strcpy( chLbl, chLineLbl(&TauLine2[ipStrong]) );
				TauIn = TauLine2[ipStrong].Emis->TauIn;
				Pump = TauLine2[ipStrong].Emis->pump;
				EscP = TauLine2[ipStrong].Emis->Pesc;
				CS = TauLine2[ipStrong].Coll.col_str;
				ColHeat = TauLine2[ipStrong].Coll.heat/
				  SDIV(heat_total);
			}
			else if( level == 3 )
				/* hyperfine levels */
			{
				strcpy( chLbl, chLineLbl(&HFLines[ipStrong]) );
				TauIn = HFLines[ipStrong].Emis->TauIn;
				Pump = HFLines[ipStrong].Emis->pump;
				EscP = HFLines[ipStrong].Emis->Pesc;
				CS = HFLines[ipStrong].Coll.col_str;
				ColHeat = HFLines[ipStrong].Coll.heat/
				  SDIV(heat_total);
			}
			else if( level == 4 )
				/* 3rd party database levels */
			{
				strcpy( chLbl, chLineLbl(dBaseLines[ipStrong].tran) );
				TauIn = dBaseLines[ipStrong].TauIn;
				Pump = dBaseLines[ipStrong].pump;
				EscP = dBaseLines[ipStrong].Pesc;
				CS = dBaseLines[ipStrong].tran->Coll.col_str;
				ColHeat = dBaseLines[ipStrong].tran->Coll.heat/
				  SDIV(heat_total);
			}
			else
			{
				fprintf( ioQQQ, " SaveHeat insane level%4ld\n", 
				  level );
				cdEXIT(EXIT_FAILURE);
			}
			fprintf( io, "  LHeat lv%2ld %10.10s TIn%10.2e Pmp%9.1e EscP%9.1e CS%9.1e Hlin/tot%10.2e\n", 
			  level, chLbl, TauIn, Pump, EscP, CS, ColHeat );
		}
	}
	for( i=0; i<LIMELM*LIMELM; ++i )
	{
		free(chLabel[i]);
	}

	free(chLabel);
	free(ipOrdered);
	free(jpsave);
	free(ipsave);
	free(SaveVal);
	return;
}
