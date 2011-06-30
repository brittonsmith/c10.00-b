/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*H2_ParseSave parse the save h2 command */
/*H2_PunchDo save some properties of the large H2 molecule */
/*chMolBranch returns a char with the spectroscopic branch of a transition */
/*H2_Prt_line_tau print line optical depths, called from premet in response to print line optical depths command*/
/*H2_PunchLineStuff include H2 lines in punched optical depths, etc, called from SaveLineStuff */
/*H2_Punch_line_data save line data for H2 molecule */
/*H2_Read_hminus_distribution read distribution function for H2 population following formation from H minus */
/*H2_ReadDissprob read dissociation probabilities and kinetic energies for all electronic levels */
/*H2_ReadEnergies read energies for all electronic levels */
/*H2_ReadTransprob read transition probabilities */
/*H2_Prt_Zone print H2 info into zone results, called from prtzone for each printed zone */
/*H2_ParseSave parse the save h2 command */
/*H2_Prt_column_density print H2 info into zone results, called from prtzone for each printed zone */
/*H2_LinesAdd add in explicit lines from the large H2 molecule, called by lines_molecules */
 /*cdH2_Line returns 1 if we found the line, 
  * or false==0 if we did not find the line because ohoto-para transition
  * or upper level has lower energy than lower level */
#include "cddefines.h" 
#include "physconst.h" 
#include "save.h" 
#include "hmi.h"
#include "prt.h"
#include "secondaries.h"
#include "grainvar.h"
#include "input.h"
#include "phycon.h"
#include "rfield.h"
#include "hyperfine.h"
#include "thermal.h"
#include "lines.h"
#include "lines_service.h"
#include "dense.h"
#include "radius.h"
#include "colden.h"
#include "taulines.h"
#include "h2.h"
#include "h2_priv.h"
#include "cddrive.h"
#include "mole.h"
#include "doppvel.h"
#include "parser.h"

/* this will say whether ortho or para,
 * H2_lgOrtho is 0 or 1 depending on whether or not ortho, 
 * so chlgPara[H2_lgOrtho] gives P or O for printing */
static char chlgPara[2]={'P','O'};

/* intensity, relative to normalization line, for faintest line to save */
static realnum thresh_punline_h2;

/*H2_LinesAdd add in explicit lines from the large H2 molecule, called by lines_molecules */
void H2_LinesAdd(void)
{
	/* these are the quantum designations of the lines we will output */
	int iRotHi, iVibHi, iElecHi ,iRotLo, iVibLo, iElecLo;

	/* H2 not on, so space not allocated */
	if( !h2.lgH2ON )
		return;

	DEBUG_ENTRY( "H2_LinesAdd()" );

	/* >>chng 05 nov 04, make info copies of these lines up here
	 * these are among the strongest lines in the 2 micron window and some have nearly the same
	 * wavelength as far weaker lines that may come before them in the line stack.  in that case
	 * cdLine would find the much weaker line with the same wavelength.
	 * put strong H2 lines first so that line search will find these, and not far weaker
	 * lines with nearly the same wavelength - these will be duplicated in the output but 
	 * these are here for into (the 'i) so does no harm
	 *
	 * the array indices in the following structures give upper and lower electronic, vib, rot quantum
	 * indices. 
	 * H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo]
	  * >>chng 05 dec 22, had hand entered wavelength in A as second parameter.  This gave
	  * rounded off result when set line precision 5 was used.  now uses same logic that
	  * PutLine will eventually use - simply enter same wl in Ang 
	 * 1-0 S(4) - 18910 */
	lindst( &H2Lines[0][1][6][0][0][4], "H2  ", 'i', false, "H2 line");
	/* 1-0 S(3) - 19570 */
	lindst( &H2Lines[0][1][5][0][0][3], "H2  ", 'i', false, "H2 line");
	/* 1-0 S(2) - 20330 */
	lindst( &H2Lines[0][1][4][0][0][2], "H2  ", 'i', false, "H2 line");
	/* 1-0 S(1) - 21210 */
	lindst( &H2Lines[0][1][3][0][0][1], "H2  ", 'i', false, "H2 line");
	/* 1-0 S(0) - 22230 */
	lindst( &H2Lines[0][1][2][0][0][0], "H2  ", 'i', false, "H2 line");
	/* start Q branch - selection rule requires that J be non-zero, so no Q(0) */
	/* 1-0 Q(2) - 24130 */
	lindst( &H2Lines[0][1][2][0][0][2], "H2  ", 'i', false, "H2 line");
	/* 1-0 Q(1) - 24060 */
	lindst( &H2Lines[0][1][1][0][0][1], "H2  ", 'i', false, "H2 line");

	/* print all lines from lowest n levels within X */
	/* loop over all possible lines and set H2_populations, 
	 * and quantities that depend on them */
	for( iElecHi=0; iElecHi<h2.nElecLevelOutput; ++iElecHi )
	{
		for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
		{
			for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
			{
				long int lim_elec_lo = 0;
				/* now the lower levels */
				/* NB - X is the only lower level considered here, since we are only 
				* concerned with excited electronic levels as a photodissociation process
				* code exists to relax this assumption - simply change following to iElecHi */
				for( iElecLo=0; iElecLo<=lim_elec_lo; ++iElecLo )
				{
					/* want to include all vib states in lower level if different electronic level,
					* but only lower vib levels if same electronic level */
					long int nv = h2.nVib_hi[iElecLo];
					if( iElecLo==iElecHi )
						nv = iVibHi;
					for( iVibLo=0; iVibLo<=nv; ++iVibLo )
					{
						long nr = h2.nRot_hi[iElecLo][iVibLo];
						if( iElecLo==iElecHi && iVibHi==iVibLo )
							nr = iRotHi-1;

						for( iRotLo=h2.Jlowest[iElecLo]; iRotLo<=nr; ++iRotLo )
						{
							if( lgH2_line_exists[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo] )
							{
								/* all ground vib state rotation lines - first is J to J-2 */
								PutLine(&H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo],
									"H2 lines");
								if( LineSave.ipass == 0 )
								{
									H2_SaveLine[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo] =  0.;
								}
								else if( LineSave.ipass == 1 )
								{
									H2_SaveLine[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo] += (realnum)(
										radius.dVeffAper*H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->xIntensity);
								}
							}
						}
					}
				}
			}
		}
	}
	return;
}

/*H2_ParseSave parse the save h2 command */
void H2_ParseSave( Parser &p ,
				   char *chHeader)
{
	DEBUG_ENTRY( "H2_ParseSave()" );

	/* this provides info on the large H2 molecule */
	if( p.nMatch("COLU") )
	{
		/* save column density */
		strcpy( save.chSave[save.nsave], "H2cl" );

		/* this is an option to scan off highest vib and rot states 
		 * to save pops - first is limit to vibration, then rotation 
		 * if no number is entered then 0 is set and all levels punched */
		/* now get vib limit */
		save.punarg[save.nsave][0] = (realnum)p.getNumberDefault(
			"H2 vibration state",0.0);

		/* highest rotation */
		save.punarg[save.nsave][1] = (realnum)p.getNumberDefault(
			"H2 rotation state",0.0);
		/* this says whether to save triplets or a matrix for output -
		 * default is triplets, so only check for matrix */
		if( p.nMatch( "MATR"  ) )
		{
			/* matrix */
			save.punarg[save.nsave][2] = 1;
			sprintf( chHeader, "#vib\trot\tcolumn density\n" );
		}
		else
		{
			/* triplets */
			save.punarg[save.nsave][2] = -1;
			sprintf( chHeader, "#vib\trot\tEner(K)\tcolden\tcolden/stat wght\tLTE colden\tLTE colden/stat wght\n" );
		}
	}
	else if( p.nMatch("COOL") )
	{
		/* heating and cooling rates */
		strcpy( save.chSave[save.nsave], "H2co" );
		sprintf( chHeader, 
			"#H2 depth\ttot cool\tTH Sol\tBig Sol\tTH pht dis\tpht dis\tTH Xcool\tXcool \n" );
	}

	else if( p.nMatch("CREA") )
	{
		/* H2 creation rates */
		strcpy( save.chSave[save.nsave], "H2cr" );
		sprintf( chHeader, 
			"#H2 depth\tH2_rate_create\tH2_rate_destroy\trate_h2_form_grains_used_total\tassoc_detach"
			"\tbh2dis\tbh2h2p\tradasc\th3ph2p\th2phmh2h\tbh2h22hh2\th3phmh2hh\th3phm2h2\th32h2\teh3_h2h\th3ph2hp"
			"\tH_CH_C_H2\tH_CHP_CP_H2\tH_CH2_CH_H2\tH_CH3P_CH2P_H2\tH_OH_O_H2\tHminus_HCOP_CO_H2\tHminus_H3OP_H2O_H2\tHminus_H3OP_OH_H2_H"
			"\tHP_CH2_CHP_H2\tHP_SiH_SiP_H2\tH2P_CH_CHP_H2\tH2P_CH2_CH2P_H2\tH2P_CO_COP_H2\tH2P_H2O_H2OP_H2\tH2P_O2_O2P_H2"
			"\tH2P_OH_OHP_H2\tH3P_C_CHP_H2\tH3P_CH_CH2P_H2\tH3P_CH2_CH3P_H2\tH3P_OH_H2OP_H2\tH3P_H2O_H3OP_H2\tH3P_CO_HCOP_H2"
			"\tH3P_O_OHP_H2\tH3P_SiH_SiH2P_H2\tH3P_SiO_SiOHP_H2\tH_CH3_CH2_H2\tH_CH4P_CH3P_H2\tH_CH5P_CH4P_H2\tH2P_CH4_CH3P_H2"
			"\tH2P_CH4_CH4P_H2\tH3P_CH3_CH4P_H2\tH3P_CH4_CH5P_H2\tHP_CH4_CH3P_H2\tHP_HNO_NOP_H2\tHP_HS_SP_H2\tH_HSP_SP_H2"
			"\tH3P_NH_NH2P_H2\tH3P_NH2_NH3P_H2\tH3P_NH3_NH4P_H2\tH3P_CN_HCNP_H2\tH3P_NO_HNOP_H2\tH3P_S_HSP_H2\tH3P_CS_HCSP_H2"
			"\tH3P_NO2_NOP_OH_H2\tH2P_NH_NHP_H2\tH2P_NH2_NH2P_H2\tH2P_NH3_NH3P_H2\tH2P_CN_CNP_H2\tH2P_HCN_HCNP_H2\tH2P_NO_NOP_H2"
			"\tH3P_Cl_HClP_H2\tH3P_HCl_H2ClP_H2\tH2P_C2_C2P_H2\tHminus_NH4P_NH3_H2\tH3P_HCN_HCNHP_H2"
			"\tdestruction/creation\tav Einstein A\n");
	}
	else if( p.nMatch("DEST") )
	{
		/* save H2 destruction - output destruction rates */
		strcpy( save.chSave[save.nsave], "H2ds" );
		sprintf( chHeader, 
			"#depth\ttot H2 rate create\ttot H2 rate destroy\ttot H- backwards\tSolomon H2g\tSolomon H2s\tphotodissoc H2s\tphotodissoc H2g"
			"\te- dissoc\trh2h2p\th2hph3p\tH0 dissoc\tCR\trheph2hpheh\theph2heh2p\thehph2h3phe\th3petc\tH2Ph3p"
			"\th2sh2g\th2h22hh2\th2sh2sh2g2h\th2sh2sh2s2h\tH2_CHP_CH2P_H\tH2_CH2P_CH3P_H\tH2_OHP_H2OP_H\tH2_H2OP_H3OP_H\tH2_COP_HCOP_H"
			"\tH2_OP_OHP_H\tH2_SiOP_SiOHP_H\tH2_C_CH_H\tH2_CP_CHP_H\tH2_CH_CH2_H\tH2_OH_H2O_H\tH2_O_OH_H"
			"\th2s_ch_ch2_h\th2s_o_oh_h\th2s_oh_h2o_h\th2s_c_ch_h\th2s_cp_chp_h\tH2_CH2_CH3_H\tH2_CH3_CH4_H"
			"\tH2_CH4P_CH5P_H\tH2s_CH2_CH3_H\tH2s_CH3_CH4_H\th2s_op_ohp_h\tH2_N_NH_H\tH2_NH_NH2_H\tH2_NH2_NH3_H\tH2_CN_HCN_H\tH2_NP_NHP_H"
			"\tH2_NHP_N_H3P\tH2_NHP_NH2P_H\tH2_NH2P_NH3P_H\tH2_NH3P_NH4P_H\tH2_CNP_HCNP_H\tH2_SP_HSP_H\tH2_CSP_HCSP_H"
			"\tH2_ClP_HClP_H\tH2_HClP_H2ClP_H\tH2_HCNP_HCNHP_H"
			"\tfrac H2g\tfrac H2s\n");
	}

	else if( p.nMatch("HEAT") )
	{
		/* heating and cooling rates */
		strcpy( save.chSave[save.nsave], "H2he" );
		sprintf( chHeader, 
			"#H2 depth\ttot Heat\tHeat(big)\tHeat(TH85)\tDissoc(Big)\tDissoc(TH85) \n" );
	}

	else if( p.nMatch("LEVE") )
	{
		/* save H2 level energies */
		strcpy( save.chSave[save.nsave], "H2le" );
		sprintf( chHeader, 
			"#H2 v\tJ\tenergy(wn)\tstat wght\tSum As" );
		char chHoldit[chN_X_COLLIDER+12];
		for( int nColl=0; nColl<N_X_COLLIDER; ++nColl )
		{
			/* labels for all colliders */
			sprintf(chHoldit,"\tCritDen %s",chH2ColliderLabels[nColl]);
			strcat( chHeader , chHoldit );
		}
		strcat( chHeader , "\n" );
	}

	else if( p.nMatch("LINE") )
	{
		/* save H2 lines - all in X */
		strcpy( save.chSave[save.nsave], "H2ln" );
		sprintf( chHeader, 
			"#H2 line\tEhi\tVhi\tJhi\tElo\tVlo\tJlo\twl(mic)\twl(lab)\tlog L or I\tI/Inorm\tExcit(hi, K)\tg_u h nu * Aul\n" );
		/* first optional number changes the threshold of weakest line to print*/
		/* fe2thresh is intensity relative to normalization line,
		 * normally Hbeta, and is set to zero in zero.c */

		/* threshold for faintest line to save, default is 1e-4 of norm line */
		thresh_punline_h2 = (realnum)p.getNumberDefaultNegImplLog(
			"faintest line to save",1e-4);

		/* lines from how many electronic states?  default is one, just X, and is
		 * obtained with GROUND keyword.  ALL will produce all lines from all levels.
		 * else, if a number is present, will be the number.  if no number, no keyword,
		 * appear then just ground */
		if( p.nMatch( "ELEC"  ) ) 
		{
			if( p.nMatch(" ALL") )
			{
				/* all electronic levels - when done, will set upper limit, the
				 * number of electronic levels actually computed, don't know this yet,
				 * so signify with negative number */
				h2.nElecLevelOutput = -1;
			}
			else if( p.nMatch("GROU") )
			{
				/* just the ground electronic state */
				h2.nElecLevelOutput = 1;
			}
			else
			{
				h2.nElecLevelOutput = (int)p.getNumberDefault(
					"electronic levels for output",1.0);
			}
		}
	}

	else if( p.nMatch(" PDR") )
	{
		/* creation and destruction processes */
		strcpy( save.chSave[save.nsave], "H2pd" );
		sprintf( chHeader, "#H2 creation, destruction. \n" );
	}
	else if( p.nMatch("POPU") )
	{
		/* save H2_populations */
		strcpy( save.chSave[save.nsave], "H2po" );

		/* this is an option to scan off highest vib and rot states 
		 * to save pops - first is limit to vibration, then rotation 
		 * if no number is entered then 0 is set and all levels punched */
		/* now get vib lim */
		save.punarg[save.nsave][0] = (realnum)p.getNumberDefault(
			"highest H2 save vibration state",0.0);

		/* this is limit to rotation quantum index */
		save.punarg[save.nsave][1] = (realnum)p.getNumberDefault(
			"highest H2 save rotation state",0.0);

		if( p.nMatch( "ZONE"  ) )
		{
			/* save v=0 pops for each zone, all along one line */
			save.punarg[save.nsave][2] = 0;
			sprintf( chHeader, "#depth\torth\tpar\te=1 rel pop\te=2 rel pop\tv,J rel pops\n" );
		}
		else
		{
			/* will not do zone output, only output at the end of the calculation
			 * now check whether to save triplets or a matrix for output -
			 * default is triplets, so only check for matrix */
			if( p.nMatch( "MATR"  ) )
			{
				/* matrix */
				save.punarg[save.nsave][2] = 1;
				sprintf( chHeader, "#vib\trot\tpops\n" );
			}
			else
			{
				/* triplets */
				save.punarg[save.nsave][2] = -1;
				sprintf( chHeader, "#vib\trot\ts\tenergy(wn)\tpops/H2\told/H2\tpops/g/H2\tdep coef\tFin(Col)\tFout(col)\tRCout\tRRout\tRCin\tRRin\n" );
			}
		}
	}

	else if( p.nMatch("RATE") )
	{
		/* save h2 rates - creation and destruction rates */
		strcpy( save.chSave[save.nsave], "H2ra" );
		sprintf( chHeader, 
			"#depth\tN(H2)\tN(H2)/u(H2)\tA_V(star)\tn(Eval)"
			"\tH2/Htot\trenorm\tfrm grn\tfrmH-\tdstTH85\tBD96\tELWERT\tBigH2\telec->H2g\telec->H2s"
			"\tG(TH85)\tG(DB96)\tCR\tEleclife\tShield(BD96)\tShield(H2)\tBigh2/G0(spc)\ttot dest"
			"\tHeatH2Dish_TH85\tHeatH2Dexc_TH85\tHeatH2Dish_BigH2\tHeatH2Dexc_BigH2\thtot\n" );
	}
	else if( p.nMatch("SOLO") )
	{
		/* rate of Solomon process then fracs of exits from each v, J level */
		strcpy( save.chSave[save.nsave], "H2so" );
		sprintf( chHeader, 
			"#depth\tSol tot\tpump/dissoc\tpump/dissoc BigH2\tavH2g\tavH2s\tH2g chem/big H2\tH2s chem/big H2\tfrac H2g BigH2\tfrac H2s BigH2\teHi\tvHi\tJHi\tvLo\tJLo\tfrac\twl(A)\n" );
	}
	else if( p.nMatch("SPEC") )
	{
		/* special save command*/
		strcpy( save.chSave[save.nsave], "H2sp" );
		sprintf( chHeader, 
			"#depth\tspecial\n" );
	}
	else if( p.nMatch("TEMP") )
	{
		/* various temperatures for neutral/molecular gas */
		strcpy( save.chSave[save.nsave], "H2te" );
		sprintf( chHeader, 
			"#depth\tH2/H\tn(1/0)\tn(ortho/para)\tT(1/0)\tT(2/0)\tT(3/0)\tT(3/1)\tT(4/0)\tT(kin)\tT(21cm)\tT_sum(1/0)\tT_sum(2/0)\tT_sum(3/0)\tT_sum(3/1)\tT_sum(4/0) \n");
	}
	else if( p.nMatch("THER") )
	{
		/* thermal heating cooling processes involving H2 */
		strcpy( save.chSave[save.nsave], "H2th" );
		sprintf( chHeader, 
			"#depth\tH2/H\tn(1/0)\tn(ortho/para)\tT(1/0)\tT(2/0)\tT(3/0)\tT(3/1)\tT(4/0)\tT(kin)\tT(21cm)\tT_sum(1/0)\tT_sum(2/0)\tT_sum(3/0)\tT_sum(3/1)\tT_sum(4/0) \n");
	}
	else
	{
		fprintf( ioQQQ, 
			" There must be a second key; they are  RATE, LINE, COOL, COLUMN, _PDR, SOLOmon, TEMP, and POPUlations\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return;
}


/*H2_Prt_Zone print H2 info into zone results, called from prtzone for each printed zone */
void H2_Prt_Zone(void)
{
	int iElecHi , iVibHi;

	/* no print if H2 not turned on, or not computed for these conditions */
	if( !h2.lgH2ON || !h2.nCallH2_this_zone )
		return;

	DEBUG_ENTRY( "H2_Prt_Zone()" );

	fprintf( ioQQQ, " H2 density   ");
	fprintf(ioQQQ,PrintEfmt("%9.2e", hmi.H2_total));

	fprintf( ioQQQ, " orth/par");
	fprintf(ioQQQ,PrintEfmt("%9.2e", h2.ortho_density / SDIV( h2.para_density )));

	iElecHi = 0;
	iVibHi = 0;
	fprintf( ioQQQ, " v0 J=0,3");
	fprintf(ioQQQ,PrintEfmt("%9.2e", H2_populations[iElecHi][iVibHi][0] / hmi.H2_total));
	fprintf(ioQQQ,PrintEfmt("%9.2e", H2_populations[iElecHi][iVibHi][1] / hmi.H2_total));
	fprintf(ioQQQ,PrintEfmt("%9.2e", H2_populations[iElecHi][iVibHi][2] / hmi.H2_total));
	fprintf(ioQQQ,PrintEfmt("%9.2e", H2_populations[iElecHi][iVibHi][3] / hmi.H2_total));

	fprintf( ioQQQ, " TOTv=0,3");
	fprintf(ioQQQ,PrintEfmt("%9.2e", pops_per_vib[iElecHi][0] / hmi.H2_total));
	fprintf(ioQQQ,PrintEfmt("%9.2e", pops_per_vib[iElecHi][1] / hmi.H2_total));
	fprintf(ioQQQ,PrintEfmt("%9.2e", pops_per_vib[iElecHi][2] / hmi.H2_total));
	fprintf(ioQQQ,PrintEfmt("%9.2e", pops_per_vib[iElecHi][3] / hmi.H2_total));
	fprintf( ioQQQ, "\n");
	return;
}

/*H2_Prt_column_density print H2 info into zone results, called from prtzone for each printed zone */
void H2_Prt_column_density(	
	/* this is stream used for io, is stdout when called by final,
	 * is save unit when save output generated */
	 FILE *ioMEAN )

{
	int iVibHi;

	/* no print if H2 not turned on, or not computed for these conditions */
	if( !h2.lgH2ON || !h2.nCallH2_this_zone )
		return;

	DEBUG_ENTRY( "H2_Prt_column_density()" );

	fprintf( ioMEAN, " H2 total   ");
	fprintf(ioMEAN,"%7.3f", log10(SDIV(colden.colden[ipCOL_H2g]+colden.colden[ipCOL_H2s])));

	fprintf( ioMEAN, " H2 ortho   ");
	fprintf(ioMEAN,"%7.3f", log10(SDIV(h2.ortho_colden)));

	fprintf( ioMEAN, " para");
	fprintf(ioMEAN,"%7.3f", log10(SDIV(h2.para_colden)));

	iVibHi = 0;
	fprintf( ioMEAN, " v0 J=0,3");
	fprintf(ioMEAN,"%7.3f", log10(SDIV(H2_X_colden[iVibHi][0])));
	fprintf(ioMEAN,"%7.3f", log10(SDIV(H2_X_colden[iVibHi][1])));
	fprintf(ioMEAN,"%7.3f", log10(SDIV(H2_X_colden[iVibHi][2])));
	fprintf(ioMEAN,"%7.3f", log10(SDIV(H2_X_colden[iVibHi][3])));

#	if 0
	fprintf( ioMEAN, "    v=0,3");
	fprintf(ioMEAN,PrintEfmt("%9.2e", pops_per_vib[iElecHi][0] / hmi.H2_total));
	fprintf(ioMEAN,PrintEfmt("%9.2e", pops_per_vib[iElecHi][1] / hmi.H2_total));
	fprintf(ioMEAN,PrintEfmt("%9.2e", pops_per_vib[iElecHi][2] / hmi.H2_total));
	fprintf(ioMEAN,PrintEfmt("%9.2e", pops_per_vib[iElecHi][3] / hmi.H2_total));
	fprintf( ioMEAN, "\n");
#	endif
	return;
}


/*H2_ReadTransprob read transition probabilities */
void H2_ReadTransprob( long int nelec )
{
	const char* cdDATAFILE[N_H2_ELEC] = 
	{
		"H2_transprob_X.dat",
		"H2_transprob_B.dat", 
		"H2_transprob_C_plus.dat",
		"H2_transprob_C_minus.dat", 
		"H2_transprob_B_primed.dat", 
		"H2_transprob_D_plus.dat",
		"H2_transprob_D_minus.dat" 
	};
	FILE *ioDATA;
	char chLine[FILENAME_PATH_LENGTH_2];
	long int i, n1, n2, n3;
	long int iVibHi , iVibLo , iRotHi , iRotLo , iElecHi , iElecLo;
	bool lgEOL;

	DEBUG_ENTRY( "H2_ReadTransprob()" );

	/* now open the data file */
	char chPath[FILENAME_PATH_LENGTH_2];
	strcpy( chPath, "h2" );
	strcat( chPath, input.chDelimiter );
	strcat( chPath, cdDATAFILE[nelec] );
	ioDATA = open_data( chPath , "r" );

	/* read the first line and check that magic number is ok */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " H2_ReadTransprob could not read first line of %s\n", cdDATAFILE[nelec]);
		cdEXIT(EXIT_FAILURE);
	}
	i = 1;
	/* magic number */
	n1 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
	n2 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
	n3 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);

	/* magic number
	 * the following is the set of numbers that appear at the start of level1.dat 01 08 10 */
	if( ( n1 != 2 ) || ( n2 != 4 ) || ( n3 != 29 ) )
	{
		fprintf( ioQQQ, 
			" H2_ReadTransprob: the version of %s is not the current version.\n", cdDATAFILE[nelec] );
		fprintf( ioQQQ, 
			" I expected to find the number 2 4 29 and got %li %li %li instead.\n" ,
			n1 , n2 , n3 );
		fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
		cdEXIT(EXIT_FAILURE);
	}

	/* read until not a comment */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
		BadRead();

	while( chLine[0]=='#' )
	{
		if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
			BadRead();
	}
	iVibHi = 1;
	while( iVibHi >= 0 )
	{
		double Aul;
		sscanf(chLine,"%li\t%li\t%li\t%li\t%li\t%li\t%le", 
			&iElecHi , &iVibHi ,&iRotHi , &iElecLo , &iVibLo , &iRotLo , &Aul );
		ASSERT( iElecHi == nelec );
		/* negative iVibHi says end of data */
		if( iVibHi < 0 )
			continue;

		ASSERT( iElecHi < N_H2_ELEC );
		ASSERT( iElecLo < N_H2_ELEC );
		/** \todo 2 the "50" here and in h2.h should be made a macro.  */
		ASSERT( iVibHi < 50 );
		ASSERT( iVibLo < 50 );

		/* check that we actually included the levels in the model representation */
		if( iVibHi <= h2.nVib_hi[iElecHi] && 
		    iVibLo <= h2.nVib_hi[iElecLo] && 
			iRotHi <= h2.nRot_hi[iElecHi][iVibHi] && 
			iRotLo <= h2.nRot_hi[iElecLo][iVibLo])
		{
			double ener = energy_wn[iElecHi][iVibHi][iRotHi] - energy_wn[iElecLo][iVibLo][iRotLo];

			/* only lines that have real Aul are added to stack.  */
			H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis = AddLine2Stack( true );

			H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->Aul = (realnum)Aul;
			/* say that this line exists */
			lgH2_line_exists[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo] = true;
			/* prints transitions with negative energies  -  should not happen */
			if( ener <= 0. )
			{
				fprintf(ioQQQ,"negative energy H2 transition\t%li\t%li\t%li\t%li\t%.2e\t%.2e\n", 
					iVibHi,iVibLo,iRotHi,iRotLo,Aul,
					H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].EnergyWN);
				ShowMe();
				cdEXIT(EXIT_FAILURE);
			}
		}
#		if 0
		/* this prints all levels with As but without energies */
		else
		{
			fprintf(ioQQQ,"no\t%li\t%li\t%li\t%li\t%.2e\n", 
				iVibHi,iVibLo,iRotHi,iRotLo,Aul);
		}
#		endif

		if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
			BadRead();
		while( chLine[0]=='#' )
		{
			if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
				BadRead();
		}
	}
	fclose( ioDATA );
	return;
}

#if 0
/*H2_Read_Cosmicray_distribution read distribution function for H2 population following cosmic ray collisional excitation */
void H2_Read_Cosmicray_distribution(void)
{
	/*>>refer	H2	cr excit	Tine, S., Lepp, S., Gredel, R., & Dalgarno, A. 1997, ApJ, 481, 282 */
	FILE *ioDATA;
	char chLine[FILENAME_PATH_LENGTH_2];
	long int i, n1, n2, n3, iVib , iRot;
	long neut_frac;
	bool lgEOL;

	DEBUG_ENTRY( "H2_Read_Cosmicray_distribution()" );

	/* now open the data file */
	char chPath[FILENAME_PATH_LENGTH_2];
	strcpy( chPath, "h2" );
	strcat( chPath, input.chDelimiter );
	strcat( chPath, "H2_CosmicRay_collision.dat" );
	ioDATA = open_data( chPath, "r" );

	/* read the first line and check that magic number is ok */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " H2_Read_Cosmicray_distribution could not read first line of %s\n", "H2_Cosmic_collision.dat");
		cdEXIT(EXIT_FAILURE);
	}

	i = 1;
	/* magic number */
	n1 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
	n2 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
	n3 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);

	/* magic number
	 * the following is the set of numbers that appear at the start of H2_Cosmic_collision.dat 01 21 03 */
	if( ( n1 != 1 ) || ( n2 != 21 ) || ( n3 != 3 ) )
	{
		fprintf( ioQQQ, 
			" H2_Read_Cosmicray_distribution: the version of %s is not the current version.\n", "H2_Cosmic_collision.dat" );
		fprintf( ioQQQ, 
			" I expected to find the number 1 21 3 and got %li %li %li instead.\n" ,
			n1 , n2 , n3 );
		fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
		cdEXIT(EXIT_FAILURE);
	}

	/* read until not a comment */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
		BadRead();

	while( chLine[0]=='#' )
	{
		if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
			BadRead();
	}

	iRot = 1;
	iVib = 1;
	neut_frac = 0;
	while( iVib >= 0 )
	{
		long int j_minus_ji;
		double a[10];

		sscanf(chLine,"%li\t%li\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", 
			&iVib ,&j_minus_ji , &a[0],&a[1],&a[2],&a[3],&a[4],&a[5],&a[6],&a[7],&a[8],&a[9] 
			);
		/* negative iVib says end of data */
		if( iVib < 0 )
			continue;

		/* cr_rate[CR_X][CR_VIB][CR_J][CR_EXIT];*/
		/* check that we actually included the levels in the model representation */
		ASSERT( iVib < CR_VIB );
		ASSERT( j_minus_ji == -2 || j_minus_ji == +2 || j_minus_ji == 0 );
		ASSERT( neut_frac < CR_X );

		/* now make i_minus_ji an array index */
		j_minus_ji = 1 + j_minus_ji/2;
		ASSERT( j_minus_ji>=0 && j_minus_ji<=2 );

		/* option to add Gaussian random mole */
		for( iRot=0; iRot<CR_J; ++iRot )
		{
			cr_rate[neut_frac][iVib][iRot][j_minus_ji] = (realnum)a[iRot];
		}
		if( mole.lgH2_NOISECOSMIC )
		{
			realnum r;
			r = (realnum)RandGauss( mole.xMeanNoise , mole.xSTDNoise );

			for( iRot=0; iRot<CR_J; ++iRot )
			{
				cr_rate[neut_frac][iVib][iRot][j_minus_ji] *= (realnum)pow(10.,(double)r);
			}
		}

		if( CR_PRINT )
		{
			fprintf(ioQQQ,"cr rate\t%li\t%li", iVib , j_minus_ji ); 
			for( iRot=0; iRot<CR_J; ++iRot )
			{ 
				fprintf(ioQQQ,"\t%.3e", cr_rate[neut_frac][iVib][iRot][j_minus_ji] );
			} 
			fprintf(ioQQQ,"\n" );
		}

		/* now get next line */
		if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
			BadRead();
		while( chLine[0]=='#' )
		{
			if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
				BadRead();
		}
	}
	fclose( ioDATA );

	return;
}
#endif

/*H2_ReadEnergies read energies for all electronic levels */
void H2_ReadEnergies( long int nelec )
{
	const char* cdDATAFILE[N_H2_ELEC] = 
	{
		"H2_energy_X.dat",
		"H2_energy_B.dat", 
		"H2_energy_C_plus.dat",
		"H2_energy_C_minus.dat", 
		"H2_energy_B_primed.dat", 
		"H2_energy_D_plus.dat",
		"H2_energy_D_minus.dat"
	};
	FILE *ioDATA;
	char chLine[FILENAME_PATH_LENGTH_2];
	long int i, n1, n2, n3, iVib , iRot;
	bool lgEOL;

	DEBUG_ENTRY( "H2_ReadEnergies()" );

	/* now open the data file */
	char chPath[FILENAME_PATH_LENGTH_2];
	strcpy( chPath, "h2" );
	strcat( chPath, input.chDelimiter );
	strcat( chPath, cdDATAFILE[nelec] );
	ioDATA = open_data( chPath, "r" );

	/* read the first line and check that magic number is ok */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " H2_ReadEnergies could not read first line of %s\n", cdDATAFILE[nelec]);
		cdEXIT(EXIT_FAILURE);
	}
	i = 1;
	/* magic number */
	n1 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
	n2 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
	n3 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);

	/* magic number
	 * the following is the set of numbers that appear at the start of level1.dat 01 08 10 */
	if( ( n1 != 2 ) || ( n2 != 4 ) || ( n3 != 29 ) )
	{
		fprintf( ioQQQ, 
			" H2_ReadEnergies: the version of %s is not the current version.\n", cdDATAFILE[nelec] );
		fprintf( ioQQQ, 
			" I expected to find the number 2 4 29 and got %li %li %li instead.\n" ,
			n1 , n2 , n3 );
		fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
		cdEXIT(EXIT_FAILURE);
	}

	/* read until not a comment */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
		BadRead();

	while( chLine[0]=='#' )
	{
		if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
			BadRead();
	}

	/* this will count the number of levels within each electronic state */
	nLevels_per_elec[nelec] = 0;

	for( iVib=0; iVib<=h2.nVib_hi[nelec]; ++iVib )
	{
		for( iRot=h2.Jlowest[nelec]; iRot<=h2.nRot_hi[nelec][iVib]; ++iRot )
		{
			i = 1;
			sscanf(chLine,"%li\t%li\t%le", &n1 , &n2 , &energy_wn[nelec][iVib][iRot] );
			ASSERT( n1 == iVib );
			ASSERT( n2 == iRot );
#			if 0
			/* in atomic units, or 1 Hartree, or two rydbergs */
			if( nelec == 0 )
			{
				/* only do this for Phillip Stancil's file */
				/* corrections are to get lowest rotation level to have energy of zero */
				energy_wn[0][iVib][iRot] = -( energy_wn[0][iVib][iRot]- 3.6118114E+04 );
			}
#			endif
			ASSERT( energy_wn[nelec][iVib][iRot]> 0. || (nelec==0 && iVib==0 && iRot==0 ) );
			/* increment number of levels within this electronic state */
			++nLevels_per_elec[nelec];

			/* now start reading next line */
			if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
				BadRead();
			while( chLine[0]=='#' )
			{
				if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
					BadRead();
			}
		}
	}
	fclose( ioDATA );
	return;
}

/*H2_ReadDissprob read dissociation probabilities and kinetic energies for all electronic levels */
void H2_ReadDissprob( long int nelec )
{
	const char* cdDATAFILE[N_H2_ELEC] = 
	{
		"H2_dissprob_X.dat",/* this does not exist and nelec == 0 is not valid */
		"H2_dissprob_B.dat", 
		"H2_dissprob_C_plus.dat",
		"H2_dissprob_C_minus.dat", 
		"H2_dissprob_B_primed.dat", 
		"H2_dissprob_D_plus.dat",
		"H2_dissprob_D_minus.dat"
	};
	FILE *ioDATA;
	char chLine[FILENAME_PATH_LENGTH_2];
	long int i, n1, n2, n3, iVib , iRot;
	bool lgEOL;

	DEBUG_ENTRY( "H2_ReadDissprob()" );

	ASSERT( nelec > 0 );

	/* now open the data file */
	char chPath[FILENAME_PATH_LENGTH_2];
	strcpy( chPath, "h2" );
	strcat( chPath, input.chDelimiter );
	strcat( chPath, cdDATAFILE[nelec] );
	ioDATA = open_data( chPath, "r" );

	/* read the first line and check that magic number is ok */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " H2_ReadDissprob could not read first line of %s\n", cdDATAFILE[nelec]);
		cdEXIT(EXIT_FAILURE);
	}
	i = 1;
	/* magic number */
	n1 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
	n2 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
	n3 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);

	/* magic number
	 * the following is the set of numbers that appear at the start of level1.dat 01 08 10 */
	if( ( n1 != 3 ) || ( n2 != 2 ) || ( n3 != 11 ) )
	{
		fprintf( ioQQQ, 
			" H2_ReadDissprob: the version of %s is not the current version.\n", cdDATAFILE[nelec] );
		fprintf( ioQQQ, 
			" I expected to find the number 3 2 11 and got %li %li %li instead.\n" ,
			n1 , n2 , n3 );
		fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
		cdEXIT(EXIT_FAILURE);
	}

	/* read until not a comment */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
		BadRead();

	while( chLine[0]=='#' )
	{
		if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
			BadRead();
	}

	for( iVib=0; iVib<=h2.nVib_hi[nelec]; ++iVib )
	{
		for( iRot=h2.Jlowest[nelec]; iRot<=h2.nRot_hi[nelec][iVib]; ++iRot )
		{
			double a, b;
			i = 1;
			sscanf(chLine,"%li\t%li\t%le\t%le", 
				&n1 , &n2 , 
				/* dissociation probability */
				&a ,
				/* dissociation kinetic energy - eV not ergs */
				&b);

			/* these have to agree if data file is valid */
			ASSERT( n1 == iVib );
			ASSERT( n2 == iRot );

			/* dissociation probability */
			H2_dissprob[nelec][iVib][iRot] = (realnum)a;
			/* dissociation kinetic energy - eV not ergs */
			H2_disske[nelec][iVib][iRot] = (realnum)b;

			/* now get next line */
			if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
				BadRead();
			while( chLine[0]=='#' )
			{
				if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
					BadRead();
			}
		}
	}
	fclose( ioDATA );
	return;
}


/*H2_Read_hminus_distribution read distribution function for H2 population following formation from H minus */
void H2_Read_hminus_distribution(void)
{
	FILE *ioDATA;
	char chLine[FILENAME_PATH_LENGTH_2];
	long int i, n1, n2, n3, iVib , iRot;
	bool lgEOL;
	double sumrate[nTE_HMINUS];
	/* set true for lots of printout */
#	define H2HMINUS_PRT	false

	DEBUG_ENTRY( "H2_Read_hminus_distribution()" );

	/* now open the data file */
	char chPath[FILENAME_PATH_LENGTH_2];
	strcpy( chPath, "h2" );
	strcat( chPath, input.chDelimiter );
	strcat( chPath, "H2_hminus_deposit.dat" );
	ioDATA = open_data( chPath, "r" );

	/* read the first line and check that magic number is ok */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " H2_Read_hminus_distribution could not read first line of %s\n", "H2_hminus_deposit.dat");
		cdEXIT(EXIT_FAILURE);
	}

	i = 1;
	/* magic number */
	n1 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
	n2 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
	n3 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);

	/* magic number
	 * the following is the set of numbers that appear at the start of H2_hminus_deposit.dat 01 08 10 */
	if( ( n1 != 2 ) || ( n2 != 10 ) || ( n3 != 17 ) )
	{
		fprintf( ioQQQ, 
			" H2_Read_hminus_distribution: the version of %s is not the current version.\n", "H2_hminus_deposit.dat" );
		fprintf( ioQQQ, 
			" I expected to find the number 2 10 17 and got %li %li %li instead.\n" ,
			n1 , n2 , n3 );
		fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
		cdEXIT(EXIT_FAILURE);
	}

	/* read until not a comment */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
		BadRead();

	while( chLine[0]=='#' )
	{
		if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
			BadRead();
	}

	/* convert temps to log */
	for(i=0; i<nTE_HMINUS; ++i )
	{
		H2_te_hminus[i] = (realnum)log10(H2_te_hminus[i]);
		sumrate[i] = 0.;
	}

	iRot = 1;
	iVib = 1;
	while( iVib >= 0 )
	{
		/* set true to print rates */

		double a[nTE_HMINUS] , ener;
		sscanf(chLine,"%li\t%li\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", 
			&iVib ,&iRot , &ener, &a[0],&a[1],&a[2] , &a[3],&a[4],&a[5] ,&a[6] 
			);
		/* negative iVib says end of data */
		if( iVib < 0 )
			continue;

		/* check that we actually included the levels in the model representation */
		ASSERT( iVib <= h2.nVib_hi[0] && 
			iRot <= h2.nRot_hi[0][iVib] );

		if( H2HMINUS_PRT )
			fprintf(ioQQQ,"hminusss\t%li\t%li", iVib , iRot );
		for( i=0; i<nTE_HMINUS; ++i )
		{
			H2_X_hminus_formation_distribution[i][iVib][iRot] = (realnum)pow(10.,-a[i]);
			sumrate[i] += H2_X_hminus_formation_distribution[i][iVib][iRot];
			if( H2HMINUS_PRT )
				fprintf(ioQQQ,"\t%.3e", H2_X_hminus_formation_distribution[i][iVib][iRot] );
		}
		if( H2HMINUS_PRT )
			fprintf(ioQQQ,"\n" );

		if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
			BadRead();
		while( chLine[0]=='#' )
		{
			if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
				BadRead();
		}
	}
	fclose( ioDATA );

	if( H2HMINUS_PRT )
	{
		/* print total rate */
		fprintf(ioQQQ," total H- formation rate ");
		/* convert temps to log */
		for(i=0; i<nTE_HMINUS; ++i )
		{
			fprintf(ioQQQ,"\t%.3e" , sumrate[i]);
		}
		fprintf(ioQQQ,"\n" );
	}

	/* convert to dimensionless factors that add to unity */
	for( iVib=0; iVib<=h2.nVib_hi[0]; ++iVib )
	{
		for( iRot=h2.Jlowest[0]; iRot<=h2.nRot_hi[0][iVib]; ++iRot )
		{
			for(i=0; i<nTE_HMINUS; ++i )
			{
				H2_X_hminus_formation_distribution[i][iVib][iRot] /= (realnum)sumrate[i];
			}
		}
	}

	if( H2HMINUS_PRT )
	{
		/* print total rate */
		fprintf(ioQQQ,"  H- distribution function ");
		for( iVib=0; iVib<=h2.nVib_hi[0]; ++iVib )
		{
			for( iRot=h2.Jlowest[0]; iRot<=h2.nRot_hi[0][iVib]; ++iRot )
			{
				fprintf(ioQQQ,"%li\t%li", iVib , iRot );
				for(i=0; i<nTE_HMINUS; ++i )
				{
					fprintf(ioQQQ,"\t%.3e", H2_X_hminus_formation_distribution[i][iVib][iRot] );
				}
				fprintf(ioQQQ,"\n" );
			}
		}
	}
	return;
}

/* ===================================================================== */
/*H2_Punch_line_data save line data for H2 molecule */
void H2_Punch_line_data(
	/* io unit for save */
	FILE* ioPUN ,
	/* save all levels if true, only subset if false */
	bool lgDoAll )
{
	long int iElecHi , iElecLo , iVibHi , iVibLo , iRotHi , iRotLo;

	if( !h2.lgH2ON )
		return;

	DEBUG_ENTRY( "H2_Punch_line_data()" );

	if( lgDoAll )
	{
		fprintf( ioQQQ, 
			" H2_Punch_line_data ALL option not implemented in H2_Punch_line_data yet 1\n" );
		cdEXIT(EXIT_FAILURE);
	}
	else
	{

		/* save line date, looping over all possible lines */
		for( iElecHi=0; iElecHi<mole.n_h2_elec_states; ++iElecHi )
		{
			for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
			{
				for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
				{
					/* now the lower levels */
					/* NB - X is the only lower level considered here, since we are only 
					* concerned with excited electronic levels as a photodissociation process
					* code exists to relax this assumption - simply change following to iElecHi */
					long int lim_elec_lo = 0;
					for( iElecLo=0; iElecLo<=lim_elec_lo; ++iElecLo )
					{
						/* want to include all vib states in lower level if different electronic level,
						 * but only lower vib levels if same electronic level */
						long int nv = h2.nVib_hi[iElecLo];
						if( iElecLo==iElecHi )
							nv = iVibHi;
						for( iVibLo=0; iVibLo<=nv; ++iVibLo )
						{
							long nr = h2.nRot_hi[iElecLo][iVibLo];
							if( iElecLo==iElecHi && iVibHi==iVibLo )
								nr = iRotHi-1;

							for( iRotLo=h2.Jlowest[iElecLo]; iRotLo<=nr; ++iRotLo )
							{
								/* only save if radiative transition exists */
								if( H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].ipCont > 0 )
								{
									/** \todo	1	add logic to deduce cs */
									H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Coll.col_str = 0.;
									/* print quantum indices */
									fprintf(ioPUN,"%2li %2li %2li %2li %2li %2li ",
										iElecHi,iVibHi,iRotHi,iElecLo,iVibLo,iRotLo );
									Save1LineData( &H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo] , ioPUN , false);
								}
							}
						}
					}
				}
			}
		}
		fprintf( ioPUN , "\n");
	}
	return;
}

/*H2_PunchLineStuff include H2 lines in punched optical depths, etc, called from SaveLineStuff */
void H2_PunchLineStuff( FILE * io , realnum xLimit  , long index)
{
	long int iElecHi , iElecLo , iVibHi , iVibLo , iRotHi , iRotLo;

	if( !h2.lgH2ON )
		return;

	DEBUG_ENTRY( "H2_PunchLineStuff()" );

	/* loop over all possible lines */
	for( iElecHi=0; iElecHi<mole.n_h2_elec_states; ++iElecHi )
	{
		for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
		{
			for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
			{
				/* now the lower levels */
				/* NB - X is the only lower level considered here, since we are only 
				* concerned with excited electronic levels as a photodissociation process
				* code exists to relax this assumption - simply change following to iElecHi */
				long int lim_elec_lo = 0;
				for( iElecLo=0; iElecLo<=lim_elec_lo; ++iElecLo )
				{
					/* want to include all vib states in lower level if different electronic level,
					* but only lower vib levels if same electronic level */
					long int nv = h2.nVib_hi[iElecLo];
					if( iElecLo==iElecHi )
						nv = iVibHi;
					for( iVibLo=0; iVibLo<=nv; ++iVibLo )
					{
						long nr = h2.nRot_hi[iElecLo][iVibLo];
						if( iElecLo==iElecHi && iVibHi==iVibLo )
							nr = iRotHi-1;

						for( iRotLo=h2.Jlowest[iElecLo]; iRotLo<=nr; ++iRotLo )
						{
							/* only save if radiative transition exists */
							if( H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].ipCont > 0 )
							{
								Save1Line( &H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo] , io , xLimit  , index,
									GetDopplerWidth(2.f*dense.AtomicWeight[ipHYDROGEN]));
							}
						}
					}
				}
			}
		}
	}

	return;
}


/*H2_Prt_line_tau print line optical depths, called from premet in response to print line optical depths command*/
void H2_Prt_line_tau(void)
{
	long int iElecHi , iElecLo , iVibHi , iVibLo , iRotHi , iRotLo;

	if( !h2.lgH2ON )
		return;

	DEBUG_ENTRY( "H2_Prt_line_tau()" );

	/* loop over all possible lines */
	for( iElecHi=0; iElecHi<mole.n_h2_elec_states; ++iElecHi )
	{
		for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
		{
			for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
			{
				/* now the lower levels */
				/* NB - X is the only lower level considered here, since we are only 
				* concerned with excited electronic levels as a photodissociation process
				* code exists to relax this assumption - simply change following to iElecHi */
				long int lim_elec_lo = 0;
				for( iElecLo=0; iElecLo<=lim_elec_lo; ++iElecLo )
				{
					/* want to include all vib states in lower level if different electronic level,
					* but only lower vib levels if same electronic level */
					long int nv = h2.nVib_hi[iElecLo];
					if( iElecLo==iElecHi )
						nv = iVibHi;
					for( iVibLo=0; iVibLo<=nv; ++iVibLo )
					{
						long nr = h2.nRot_hi[iElecLo][iVibLo];
						if( iElecLo==iElecHi && iVibHi==iVibLo )
							nr = iRotHi-1;

						for( iRotLo=h2.Jlowest[iElecLo]; iRotLo<=nr; ++iRotLo )
						{
							/* only print if radiative transition exists */
							if( H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].ipCont > 0 )
							{
								prme(" c",&H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo] );
							}
						}
					}
				}
			}
		}
	}

	return;
}


/*chMolBranch returns a char with the spectroscopic branch of a transition */
STATIC char chMolBranch( long iRotHi , long int iRotLo )
{
	/* these are the spectroscopic branches */
	char chBranch[5] = {'O','P','Q','R','S'};
	/* this is the index within the chBranch array */
	int ip = 2 + (iRotHi - iRotLo);
	if( ip<0 || ip>=5 )
	{
		fprintf(ioQQQ," chMolBranch called with insane iRotHi=%li iRotLo=%li ip=%i\n",
			iRotHi , iRotLo , ip );
		ip = 0;
	}

	return( chBranch[ip] );
}

/*H2_PunchDo save some properties of the large H2 molecule */
void H2_PunchDo( FILE* io ,  char chJOB[] , const char chTime[] , long int ipPun )
{
	long int iVibHi , iElecHi , iRotHi , iVibLo , iElecLo , iRotLo,
		ip;
	long int iRot , iVib;
	long int LimVib , LimRot;

	DEBUG_ENTRY( "H2_PunchDo()" );

	/* which job are we supposed to do? This routine is active even when H2 is not turned on
	 * so do not test on h2.lgH2ON initially */

	/* H2 populations computed in last zone - 
	 * give all of molecule in either matrix or triplet format */
	if( (strcmp( chJOB , "H2po" ) == 0) && (strcmp(chTime,"LAST") == 0) &&
		(save.punarg[ipPun][2] != 0) )
	{
		/* >>chng 04 feb 19, do not save if H2 not yet evaluated */
		if( h2.lgH2ON  && hmi.lgBigH2_evaluated )
		{
			iVibHi= 0;
			iRotHi = 0;
			iElecHi=0;
			/* the limit to the number of vibration levels punched -
			* default is all, but first two numbers on save h2 pops command
			* reset limit */
			/* this is limit to vibration */
			if( save.punarg[ipPun][0] > 0 )
			{
				LimVib = (long)save.punarg[ipPun][0];
			}
			else
			{
				LimVib = h2.nVib_hi[iElecHi];
			}

			/* first save the current ortho, para, and total H2 density */
			fprintf(io,"%i\t%i\t%.3e\tortho\n", 
				103 , 
				103 ,
				h2.ortho_density );
			fprintf(io,"%i\t%i\t%.3e\tpara\n", 
				101 , 
				101 ,
				h2.para_density );
			fprintf(io,"%i\t%i\t%.3e\ttotal\n", 
				0 , 
				0 ,
				hmi.H2_total );

			/* now save the actual H2_populations, first part both matrix and triplets */
			for( iVibHi=0; iVibHi<=LimVib; ++iVibHi )
			{
				/* this is limit to rotation quantum index */
				if( save.punarg[ipPun][1] > 0 )
				{
					LimRot = (long)MIN2(
						save.punarg[ipPun][1] , (realnum)h2.nRot_hi[iElecHi][iVibHi]);
				}
				else
				{
					LimRot = h2.nRot_hi[iElecHi][iVibHi];
				}
				if( save.punarg[ipPun][2] > 0 )
				{
					long int i;
					/* this option save matrix */
					if( iVibHi == 0 )
					{
						fprintf(io,"vib\\rot");
						/* this is first vib, so make row of rot numbs */
						for( i=0; i<=LimRot; ++i )
						{
							fprintf(io,"\t%li",i);
						}
						fprintf(io,"\n");
					}
					fprintf(io,"%li",iVibHi );
					for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=LimRot; ++iRotHi )
					{
						fprintf(io,"\t%.3e", 
							H2_populations[iElecHi][iVibHi][iRotHi]/hmi.H2_total );
					}
					fprintf(io,"\n" );
				}
				else if( save.punarg[ipPun][2] < 0 )
				{
					/* this option save triplets - the default */
					for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=LimRot; ++iRotHi )
					{
						fprintf(io,"%li\t%li\t%c\t%.1f\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n", 
							/* upper vibration and rotation quantum numbers */
							iVibHi , iRotHi ,
							/* an 'O' or 'P' for ortho or para */
							chlgPara[H2_lgOrtho[iElecHi][iVibHi][iRotHi]],
							/* the level excitation energy in wavenumbers */
							energy_wn[iElecHi][iVibHi][iRotHi],
							/* actual population relative to total H2 */
							H2_populations[iElecHi][iVibHi][iRotHi]/hmi.H2_total ,
							/* old level H2_populations for comparison */
							H2_old_populations[iElecHi][iVibHi][iRotHi]/hmi.H2_total ,
							/* H2_populations per h2 and per statistical weight */
							H2_populations[iElecHi][iVibHi][iRotHi]/hmi.H2_total/H2_stat[iElecHi][iVibHi][iRotHi] ,
							/* LTE departure coefficient */
							/* >>chng 05 jan 26, missing factor of H2 abundance LTE is norm to unity, not tot abund */
							H2_populations[iElecHi][iVibHi][iRotHi]/SDIV(H2_populations_LTE[iElecHi][iVibHi][iRotHi]*hmi.H2_total ) ,
							/* fraction of exits that were collisional */
							H2_col_rate_out[iVibHi][iRotHi]/SDIV(H2_col_rate_out[iVibHi][iRotHi]+H2_rad_rate_out[0][iVibHi][iRotHi]) ,
							/* fraction of entries that were collisional */
							H2_col_rate_in[iVibHi][iRotHi]/SDIV(H2_col_rate_in[iVibHi][iRotHi]+H2_rad_rate_in[iVibHi][iRotHi]),
							/* collisions out */
							H2_col_rate_out[iVibHi][iRotHi],
							/* radiation out */
							H2_rad_rate_out[0][iVibHi][iRotHi] ,
							/* radiation out */
							H2_col_rate_in[iVibHi][iRotHi],
							/* radiation in */
							H2_rad_rate_in[iVibHi][iRotHi]
							);
					}
				}
			}
		}
	}
	/* save H2 populations for each zone 
	 * H2_populations of v=0 for each zone */
	else if( (strcmp( chJOB , "H2po" ) == 0) && (strcmp(chTime,"LAST") != 0) &&
		(save.punarg[ipPun][2] == 0) )
	{
		/* >>chng 04 feb 19, do not save if h2 not yet evaluated */
		if( h2.lgH2ON  && hmi.lgBigH2_evaluated )
		{
			fprintf(io,"%.5e\t%.3e\t%.3e", radius.depth_mid_zone , 
				h2.ortho_density , h2.para_density);
			/* rel pops of first two excited electronic states */
			fprintf(io,"\t%.3e\t%.3e", 
				pops_per_elec[1] , pops_per_elec[2]);
			iElecHi = 0;
			iVibHi = 0;
			/* this is limit to vibration quantum index */
			if( save.punarg[ipPun][0] > 0 )
			{
				LimVib = (long)save.punarg[ipPun][1];
			}
			else
			{
				LimVib = h2.nRot_hi[iElecHi][iVibHi];
			}
			LimVib = MIN2( LimVib , h2.nVib_hi[iElecHi] );
			/* this is limit to rotation quantum index */
			if( save.punarg[ipPun][1] > 0 )
			{
				LimRot = (long)save.punarg[ipPun][1];
			}
			else
			{
				LimRot = h2.nRot_hi[iElecHi][iVibHi];
			}
			for( iVibHi = 0; iVibHi<=LimVib; ++iVibHi )
			{
				fprintf(io,"\tv=%li",iVibHi);
				long int LimRotVib = MIN2( LimRot , h2.nRot_hi[iElecHi][iVibHi] );
				for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=LimRotVib; ++iRotHi )
				{
					fprintf(io,"\t%.3e", 
						H2_populations[iElecHi][iVibHi][iRotHi]/hmi.H2_total );
				}
			}
			fprintf(io,"\n");
		}
	}

	/* save column densities */
	else if( (strcmp( chJOB , "H2cl" ) == 0) && (strcmp(chTime,"LAST") == 0) )
	{
		iVibHi= 0;
		iRotHi = 0;
		iElecHi=0;
		/* the limit to the number of vibration levels punched -
		 * default is all, but first two numbers on save h2 pops command
		 * reset limit */
		/* this is limit to vibration */
		if( save.punarg[ipPun][0] > 0 )
		{
			LimVib = (long)save.punarg[ipPun][0];
		}
		else
		{
			LimVib = h2.nVib_hi[iElecHi];
		}

		/* first save ortho and para H2_populations */
		fprintf(io,"%i\t%i\t%.3e\tortho\n", 
			103 , 
			103 ,
			h2.ortho_colden );
		fprintf(io,"%i\t%i\t%.3e\tpara\n", 
			101 , 
			101 ,
			h2.para_colden );
		/* total H2 column density */
		fprintf(io,"%i\t%i\t%.3e\ttotal\n", 
			0 , 
			0 ,
			colden.colden[ipCOL_H2g]+colden.colden[ipCOL_H2s]);

		/* save level column densities */
		for( iVibHi=0; iVibHi<=LimVib; ++iVibHi )
		{
		if( h2.lgH2ON )
		{
			/* this is limit to rotation quantum index */
			if( save.punarg[ipPun][1] > 0 )
			{
				LimRot = (long)save.punarg[ipPun][1];
			}
			else
			{
				LimRot = h2.nRot_hi[iElecHi][iVibHi];
			}
			if( save.punarg[ipPun][2] > 0 )
			{
				long int i;
				/* save matrix */
				if( iVibHi == 0 )
				{
					fprintf(io,"vib\\rot");
					/* this is first vib, so make row of rot numbs */
					for( i=0; i<=LimRot; ++i )
					{
						fprintf(io,"\t%li",i);
					}
					fprintf(io,"\n");
				}
				fprintf(io,"%li",iVibHi );
				for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=LimRot; ++iRotHi )
				{
					fprintf(io,"\t%.3e", 
						H2_X_colden[iVibHi][iRotHi]/hmi.H2_total );
				}
				fprintf(io,"\n" );
			}
			else
			{
				/* save triplets - the default */
				for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=LimRot; ++iRotHi )
				{
					fprintf(io,"%li\t%li\t%.1f\t%.3e\t%.3e\t%.3e\t%.3e\n", 
						iVibHi , 
						iRotHi ,
						/* energy relative to 0,0, T1CM converts wavenumber to K */
						energy_wn[iElecHi][iVibHi][iRotHi]*T1CM,
						/* these are column densities for actual molecule */
						H2_X_colden[iVibHi][iRotHi] ,
						H2_X_colden[iVibHi][iRotHi]/H2_stat[iElecHi][iVibHi][iRotHi] ,
						/* these are same column densities but for LTE populations */
						H2_X_colden_LTE[iVibHi][iRotHi] ,
						H2_X_colden_LTE[iVibHi][iRotHi]/H2_stat[iElecHi][iVibHi][iRotHi]);
				}
			}
		}
		}
	}
	else if( (strcmp(chJOB , "H2pd" ) == 0) && (strcmp(chTime,"LAST") != 0) )
	{
		/* save PDR 
		 * output some PDR information (densities, rates) for each zone */
		fprintf(io,"%.5e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n", 
			/* depth in cm */
			radius.depth_mid_zone ,
			/* the computed ortho and para densities */
			h2.ortho_density , 
			h2.para_density ,
			/* the Lyman Werner band dissociation, Tielens & Hollenbach */
			hmi.H2_Solomon_dissoc_rate_TH85_H2g , 
			/* the Lyman Werner band dissociation, Bertoldi & Draine */
			hmi.H2_Solomon_dissoc_rate_BD96_H2g,
			/* the Lyman Werner band dissociation, big H2 mole */
			hmi.H2_Solomon_dissoc_rate_BigH2_H2g);
	}
	else if( (strcmp(chJOB , "H2co" ) == 0) && (strcmp(chTime,"LAST") != 0) )
	{
		/* save H2 cooling - do heating cooling for each zone old new H2 */
		fprintf(io,"%.5e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n", 
			/* depth in cm */
			radius.depth_mid_zone ,
			/* total cooling, equal to total heating */
			thermal.ctot , 
			/* H2 destruction by Solomon process, TH85 rate */
			hmi.H2_Solomon_dissoc_rate_TH85_H2g,
			/* H2 destruction by Solomon process, big H2 model rate */
			hmi.H2_Solomon_dissoc_rate_BigH2_H2g +
				hmi.H2_Solomon_dissoc_rate_BigH2_H2s,
			/* H2 photodissociation heating, eqn A9 of Tielens & Hollenbach 1985a */
			hmi.HeatH2Dish_TH85,
			/* heating due to dissociation of electronic excited states */
			hmi.HeatH2Dish_BigH2 , 
			/* cooling (usually neg and so heating) due to collisions within X */
			hmi.HeatH2Dexc_TH85,
			hmi.HeatH2Dexc_BigH2 
			);

	}
	else if( (strcmp(chJOB , "H2cr" ) == 0) && (strcmp(chTime,"LAST") != 0) )
	{
		/* PUNCH H2 CREATION - show H2 creation processes for each zone */
		/* >>chng 05 jul 15, TE, save all H2 creation processes, unit cm-3 s-1 */
		fprintf(io,"%.5e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e", 
			/* depth in cm */
			radius.depth_mid_zone ,
			/* creation cm-3 s-1, destruction rate, s-1 */
			hmi.H2_rate_create, 
			hmi.H2_rate_destroy, 
			gv.rate_h2_form_grains_used_total * dense.xIonDense[ipHYDROGEN][0] / hmi.H2_rate_create, 
			hmi.assoc_detach * hmi.Hmolec[ipMH]*hmi.Hmolec[ipMHm] / hmi.H2_rate_create, 
			hmi.bh2dis * dense.xIonDense[ipHYDROGEN][0] / hmi.H2_rate_create, 
			hmi.bh2h2p * dense.xIonDense[ipHYDROGEN][0] * hmi.Hmolec[ipMH2p] / hmi.H2_rate_create, 
			hmi.radasc * dense.xIonDense[ipHYDROGEN][0] / hmi.H2_rate_create, 
			hmi.h3ph2p * dense.xIonDense[ipHYDROGEN][0] * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create, 
			hmi.h2phmh2h * hmi.Hmolec[ipMH2p] * hmi.Hmolec[ipMHm] / hmi.H2_rate_create,
			hmi.bh2h22hh2 * 2 * dense.xIonDense[ipHYDROGEN][0] * hmi.Hmolec[ipMH2g] / hmi.H2_rate_create,
			hmi.h3phmh2hh * hmi.Hmolec[ipMH3p] * hmi.Hmolec[ipMHm] / hmi.H2_rate_create,
			hmi.h3phm2h2 * hmi.Hmolec[ipMH3p] * hmi.Hmolec[ipMHm] / hmi.H2_rate_create,
			hmi.h32h2 * hmi.Hmolec[ipMH2p] * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create,
			hmi.eh3_h2h * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create,
			hmi.h3ph2hp * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create
			);
		fprintf(io,"\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e", 
		/*chemical network*/
			/*light elements*/
			co.H_CH_C_H2 * dense.xIonDense[ipHYDROGEN][0] / hmi.H2_rate_create,
			co.H_CHP_CP_H2 * dense.xIonDense[ipHYDROGEN][0] / hmi.H2_rate_create,
			co.H_CH2_CH_H2 * dense.xIonDense[ipHYDROGEN][0] / hmi.H2_rate_create,
			co.H_CH3P_CH2P_H2 * dense.xIonDense[ipHYDROGEN][0] / hmi.H2_rate_create,
			co.H_OH_O_H2 * dense.xIonDense[ipHYDROGEN][0] / hmi.H2_rate_create,
			co.Hminus_HCOP_CO_H2 * hmi.Hmolec[ipMHm] / hmi.H2_rate_create,
			co.Hminus_H3OP_H2O_H2 * hmi.Hmolec[ipMHm] / hmi.H2_rate_create,
			co.Hminus_H3OP_OH_H2_H * hmi.Hmolec[ipMHm] / hmi.H2_rate_create,
			co.HP_CH2_CHP_H2 * hmi.Hmolec[ipMHp] / hmi.H2_rate_create,
			co.HP_SiH_SiP_H2* hmi.Hmolec[ipMHp] / hmi.H2_rate_create,
			co.H2P_CH_CHP_H2 * hmi.Hmolec[ipMH2p] / hmi.H2_rate_create,
			co.H2P_CH2_CH2P_H2 * hmi.Hmolec[ipMH2p] / hmi.H2_rate_create,
			co.H2P_CO_COP_H2 * hmi.Hmolec[ipMH2p] / hmi.H2_rate_create,
			co.H2P_H2O_H2OP_H2 * hmi.Hmolec[ipMH2p] / hmi.H2_rate_create,
			co.H2P_O2_O2P_H2 * hmi.Hmolec[ipMH2p] / hmi.H2_rate_create,
			co.H2P_OH_OHP_H2 * hmi.Hmolec[ipMH2p] / hmi.H2_rate_create,
			co.H3P_C_CHP_H2 * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create,
			co.H3P_CH_CH2P_H2 * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create,
			co.H3P_CH2_CH3P_H2 * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create,
			co.H3P_OH_H2OP_H2 * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create,
			co.H3P_H2O_H3OP_H2 * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create,
			co.H3P_CO_HCOP_H2 * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create,
			co.H3P_O_OHP_H2 * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create,
			co.H3P_SiH_SiH2P_H2 * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create,
			co.H3P_SiO_SiOHP_H2	 * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create,
			co.H_CH3_CH2_H2 * dense.xIonDense[ipHYDROGEN][0] / hmi.H2_rate_create,
			co.H_CH4P_CH3P_H2 * dense.xIonDense[ipHYDROGEN][0] / hmi.H2_rate_create,
			co.H_CH5P_CH4P_H2 * dense.xIonDense[ipHYDROGEN][0] / hmi.H2_rate_create,
			co.H2P_CH4_CH3P_H2 * hmi.Hmolec[ipMH2p] / hmi.H2_rate_create,
			co.H2P_CH4_CH4P_H2 * hmi.Hmolec[ipMH2p] / hmi.H2_rate_create,
			co.H3P_CH3_CH4P_H2 * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create,
			co.H3P_CH4_CH5P_H2 * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create,
			co.HP_CH4_CH3P_H2 * hmi.Hmolec[ipMHp] / hmi.H2_rate_create,
			/* heavy elements */
			co.HP_HNO_NOP_H2  * hmi.Hmolec[ipMHp] / hmi.H2_rate_create,
			co.HP_HS_SP_H2 * hmi.Hmolec[ipMHp] / hmi.H2_rate_create,
			co.H_HSP_SP_H2 * dense.xIonDense[ipHYDROGEN][0] / hmi.H2_rate_create,
			co.H3P_NH_NH2P_H2 * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create,
			co.H3P_NH2_NH3P_H2 * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create,
			co.H3P_NH3_NH4P_H2 * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create,
			co.H3P_CN_HCNP_H2 * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create,
			co.H3P_NO_HNOP_H2 * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create,
			co.H3P_S_HSP_H2 * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create,
			co.H3P_CS_HCSP_H2 * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create,
			co.H3P_NO2_NOP_OH_H2 * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create,
			co.H2P_NH_NHP_H2 * hmi.Hmolec[ipMH2p] / hmi.H2_rate_create,
			co.H2P_NH2_NH2P_H2 * hmi.Hmolec[ipMH2p] / hmi.H2_rate_create,
			co.H2P_NH3_NH3P_H2 * hmi.Hmolec[ipMH2p] / hmi.H2_rate_create,
			co.H2P_CN_CNP_H2 * hmi.Hmolec[ipMH2p] / hmi.H2_rate_create,
			co.H2P_HCN_HCNP_H2 * hmi.Hmolec[ipMH2p] / hmi.H2_rate_create,
			co.H2P_NO_NOP_H2 * hmi.Hmolec[ipMH2p] / hmi.H2_rate_create, 
			co.H3P_Cl_HClP_H2 * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create,
			co.H3P_HCl_H2ClP_H2 * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create,
			co.H2P_C2_C2P_H2 * hmi.Hmolec[ipMH2p] / hmi.H2_rate_create,
			co.Hminus_NH4P_NH3_H2 * hmi.Hmolec[ipMHm] / hmi.H2_rate_create,
			co.H3P_HCN_HCNHP_H2 * hmi.Hmolec[ipMH3p] / hmi.H2_rate_create
			);
		fprintf(io,"\t%.3e\t%.3e\n",
			hmi.H2_rate_destroy  * hmi.H2_total / hmi.H2_rate_create,
			hmi.h2s_sp_decay
			);
	}
	else if( (strcmp(chJOB , "H2ds" ) == 0) && (strcmp(chTime,"LAST") != 0) )
	{
		/* save H2 destruction - show H2 destruction processes for each zone 
		 * >>chng 05 nov 17, TE, added the new reaction H2s + O+ -> OH+ + H 
		 * >>chng 05 oct 04, TE, remove eh2hhm(was double) and include dissociation by electrons, H2g/H2s + e -> 2H  
		 * >>chng 05 jul 15, TE, save all H2 destruction rates, weighted by fractional abundance of H2g, H2s, unit s-1 */
		fprintf(io,"%.5e\t%.2e\t%.2e", 
			/* depth in cm */
			radius.depth_mid_zone ,
			/* total H2 creation rate, cm-3 s-1 */
			hmi.H2_rate_create, 
			/* destruction rate, s-1 */
			hmi.H2_rate_destroy );
		/* this can be zero following a high-Te init temperature search
		 * abort */
		if( hmi.H2_total > 0. )
		{
			fprintf(io,"\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e", 
				/* H2 + e -> H- + H0 */
				(hmi.assoc_detach_backwards_grnd * hmi.Hmolec[ipMH2g] + hmi.assoc_detach_backwards_exct * hmi.Hmolec[ipMH2s]) / hmi.H2_rate_destroy  /	hmi.H2_total, 
				/*photons*/
				hmi.H2_Solomon_dissoc_rate_used_H2g / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total, 
				hmi.H2_Solomon_dissoc_rate_used_H2s / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2s] / hmi.H2_total,
				hmi.H2_photodissoc_used_H2s / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2s] / hmi.H2_total, 
				hmi.H2_photodissoc_used_H2g / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total, 
				/*electrons*/
				(hmi.h2ge2h * hmi.Hmolec[ipMH2g] + hmi.h2se2h * hmi.Hmolec[ipMH2s]) / hmi.H2_rate_destroy  / hmi.H2_total,
				/*H+*/
				hmi.rh2h2p*dense.xIonDense[ipHYDROGEN][1] / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total, 
				hmi.h2hph3p*dense.xIonDense[ipHYDROGEN][1] / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				/*H0*/
				(hmi.rh2dis * hmi.Hmolec[ipMH2g] + hmi.h2sh * hmi.Hmolec[ipMH2s]) * dense.xIonDense[ipHYDROGEN][0] / hmi.H2_rate_destroy  / hmi.H2_total,
				/*CR*/
				(hmi.CR_reac_H2g * hmi.Hmolec[ipMH2g] + hmi.CR_reac_H2s * hmi.Hmolec[ipMH2s]) / hmi.H2_rate_destroy / hmi.H2_total,
				/*He+,HeH+*/
				hmi.rheph2hpheh*dense.xIonDense[ipHELIUM][1] / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				hmi.heph2heh2p*dense.xIonDense[ipHELIUM][1] / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				hmi.hehph2h3phe*hmi.Hmolec[ipMHeHp] / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				/*H3+*/
				hmi.h3petc*hmi.Hmolec[ipMH3p] / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				/*H2+*/
				hmi.h2ph3p*hmi.Hmolec[ipMH2p] / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				/*H2s+H2g -> H2g + 2H*/
				hmi.h2sh2g * hmi.Hmolec[ipMH2g] / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2s] / hmi.H2_total,
				/*H2g+H2g -> H2g + 2H*/
				hmi.h2h22hh2*2*hmi.Hmolec[ipMH2g] / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				/*H2s+H2s -> H2g + 2H*/
				hmi.h2sh2sh2g2h*2*hmi.Hmolec[ipMH2s] / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2s] / hmi.H2_total,
				/*H2s+H2s -> H2s + 2H*/
				hmi.h2sh2sh2s2h*2*hmi.Hmolec[ipMH2s] / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2s] / hmi.H2_total,
				/*chemical network*/
					/*light elements*/
				co.H2_CHP_CH2P_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_CH2P_CH3P_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_OHP_H2OP_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_H2OP_H3OP_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_COP_HCOP_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_OP_OHP_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_SiOP_SiOHP_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_C_CH_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_CP_CHP_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_CH_CH2_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_OH_H2O_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_O_OH_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2s_CH_CH2_H /  hmi.H2_rate_destroy * hmi.Hmolec[ipMH2s] / hmi.H2_total,
				co.H2s_O_OH_H /  hmi.H2_rate_destroy * hmi.Hmolec[ipMH2s] / hmi.H2_total,
				co.H2s_OH_H2O_H /  hmi.H2_rate_destroy * hmi.Hmolec[ipMH2s] / hmi.H2_total,
				co.H2s_C_CH_H /  hmi.H2_rate_destroy * hmi.Hmolec[ipMH2s] / hmi.H2_total,	
				co.H2s_CP_CHP_H /  hmi.H2_rate_destroy * hmi.Hmolec[ipMH2s] / hmi.H2_total,
				co.H2_CH2_CH3_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_CH3_CH4_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_CH4P_CH5P_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2s_CH2_CH3_H /  hmi.H2_rate_destroy * hmi.Hmolec[ipMH2s] / hmi.H2_total,
				co.H2s_CH3_CH4_H /  hmi.H2_rate_destroy * hmi.Hmolec[ipMH2s] / hmi.H2_total,
				co.H2s_OP_OHP_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2s] / hmi.H2_total,
				/*heavy elements*/
				co.H2_N_NH_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_NH_NH2_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_NH2_NH3_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_CN_HCN_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_NP_NHP_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_NHP_N_H3P / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_NHP_NH2P_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_NH2P_NH3P_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_NH3P_NH4P_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_CNP_HCNP_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_SP_HSP_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_CSP_HCSP_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_ClP_HClP_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_HClP_H2ClP_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total,
				co.H2_HCNP_HCNHP_H / hmi.H2_rate_destroy * hmi.Hmolec[ipMH2g] / hmi.H2_total
				);
			fprintf(io,"\t%.4e\t%.4e",
				/*H2g/Htot, H2s/Htot chemical network and from big molecule model*/
				hmi.Hmolec[ipMH2g] / hmi.H2_total,
				hmi.Hmolec[ipMH2s] / hmi.H2_total
				);
		}
		fprintf(io,"\n");
	}

	else if( (strcmp(chJOB , "H2le" ) == 0) && (strcmp(chTime,"LAST") == 0) )
	{
		/* save H2 levels */
		for( long int ipHi=0; ipHi < nLevels_per_elec[0]; ipHi++ )
		{
			long int ipVJ_Hi = H2_ipX_ener_sort[ipHi];
			iRotHi = ipRot_H2_energy_sort[ipVJ_Hi];
			iVibHi = ipVib_H2_energy_sort[ipVJ_Hi];
			double Asum , Csum[N_X_COLLIDER];
			long int nColl;
			Asum = 0;
			for( nColl=0; nColl<N_X_COLLIDER; ++nColl )
				Csum[nColl] = 0.;
			for( long int ipLo=0; ipLo<ipHi; ++ipLo )
			{
				/* all lower levels */
				long int ipVJ_Lo = H2_ipX_ener_sort[ipLo];
				iRotLo = ipRot_H2_energy_sort[ipVJ_Lo];
				iVibLo = ipVib_H2_energy_sort[ipVJ_Lo];

				/* radiative decays down */
				if( ( abs(iRotHi-iRotLo) == 2 || (iRotHi-iRotLo) == 0 ) && iVibLo <= iVibHi &&
				    lgH2_line_exists[0][iVibHi][iRotHi][0][iVibLo][iRotLo] )
				{
					Asum += H2Lines[0][iVibHi][iRotHi][0][iVibLo][iRotLo].Emis->Aul*(
						H2Lines[0][iVibHi][iRotHi][0][iVibLo][iRotLo].Emis->Pesc + 
						H2Lines[0][iVibHi][iRotHi][0][iVibLo][iRotLo].Emis->Pdest +
						H2Lines[0][iVibHi][iRotHi][0][iVibLo][iRotLo].Emis->Pelec_esc);
				}
				/* all collisions down */
				mr5ci H2cr = H2_CollRate.begin(iVibHi,iRotHi,iVibLo,iRotLo);
				for( nColl=0; nColl<N_X_COLLIDER; ++nColl )
					Csum[nColl] += H2cr[nColl];
			}

			/* save H2 level energies */
			fprintf(io,"%li\t%li\t%.2f\t%li\t%.3e", 
				iVibHi , iRotHi,
				energy_wn[0][iVibHi][iRotHi],
				(long)H2_stat[0][iVibHi][iRotHi],
				Asum );
			for( nColl=0; nColl<N_X_COLLIDER; ++nColl )
				/* sum over all lower levels */
				fprintf(io,"\t%.3e",Csum[nColl]);
			fprintf(io,"\n");
		}
	}

	else if( (strcmp(chJOB , "H2ra" ) == 0) && (strcmp(chTime,"LAST") != 0) )
	{
		/* save h2 rates - some rates and lifetimes */
		double sumpop = 0. , sumlife = 0.;

		/* this block, find lifetime against photo excitation into excited electronic states */
		iElecLo = 0;
		iVibLo = 0;
		if( h2.lgH2ON && hmi.lgBigH2_evaluated )
		{
			for( iElecHi=1; iElecHi<mole.n_h2_elec_states; ++iElecHi )
			{
				for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
				{
					for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
					{
						for( iRotLo=h2.Jlowest[iElecLo]; iRotLo<=h2.nRot_hi[iElecLo][iVibLo]; ++iRotLo )
						{
							/* only do if radiative transition exists */
							if( H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].ipCont > 0 )
							{
								sumlife +=
									H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->pump *
									H2_populations[iElecLo][iVibLo][iRotLo];
								sumpop +=
									H2_populations[iElecLo][iVibLo][iRotLo];
							}
						}
					}
				}
			}
		}

		/* continue output from save h2 rates command */
		/* find photoexcitation rates from v=0 */
		/* PDR information for each zone */
		fprintf(io,
			"%.5e\t%.3e\t%.3e\t%.3e\t%li", 
			/* depth in cm */
			radius.depth_mid_zone ,
			/* the column density (cm^-2) in H2 */
			colden.colden[ipCOL_H2g]+colden.colden[ipCOL_H2s],
			/* this is a special form of column density - should be proportional
			 * to total shielding */
			colden.coldenH2_ov_vel ,
			/* visual extinction due to dust alone, of point source (star)*/
			rfield.extin_mag_V_point,
			/* number of large molecule evaluations in this zone */
			h2.nCallH2_this_zone );
		fprintf(io,
			"\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e",
			/* total H2 fraction */
			hmi.H2_total/dense.gas_phase[ipHYDROGEN] ,
			/* chemistry renorm factor */
			H2_renorm_chemistry,
			/* rate H2 forms on grains */
			gv.rate_h2_form_grains_used_total , 
			/* rate H2 forms by H minus route */
			hmi.Hmolec[ipMHm]*1.35e-9,
			/* H2 destruction by Solomon process, TH85 rate */
			hmi.H2_Solomon_dissoc_rate_TH85_H2g + hmi.H2_Solomon_dissoc_rate_TH85_H2s,
			/* H2 destruction by Solomon process, Bertoldi & Draine rate */
			hmi.H2_Solomon_dissoc_rate_BD96_H2g + hmi.H2_Solomon_dissoc_rate_BD96_H2s,
			/* H2 destruction by Solomon process, Elwert et al. in preparation */
			hmi.H2_Solomon_dissoc_rate_ELWERT_H2g + hmi.H2_Solomon_dissoc_rate_ELWERT_H2g,
			/* H2 destruction by Solomon process, big H2 model rate */
			hmi.H2_Solomon_dissoc_rate_BigH2_H2g + hmi.H2_Solomon_dissoc_rate_BigH2_H2s,
			/* rate s-1 H2 electronic excit states decay into H2g */
			hmi.H2_Solomon_elec_decay_H2g ,
			/* rate s-1 H2 electronic excit states decay into H2s */
			hmi.H2_Solomon_elec_decay_H2s 
			);
		fprintf(io,
			"\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e", 
			/* The TH85 estimate of the radiation field relative to the Habing value */
			hmi.UV_Cont_rel2_Habing_TH85_depth,
			/* The DB96 estimate of the radiation field relative to the Habing value */
			hmi.UV_Cont_rel2_Draine_DB96_depth,
			/* cosmic ray ionization rate */
			secondaries.csupra[ipHYDROGEN][0]*0.93,
			sumlife/SDIV( sumpop ) ,
			hmi.H2_Solomon_dissoc_rate_BD96_H2g/SDIV(hmi.UV_Cont_rel2_Habing_TH85_depth) ,
			hmi.H2_Solomon_dissoc_rate_BigH2_H2g/SDIV(hmi.UV_Cont_rel2_Habing_TH85_depth),
			hmi.H2_Solomon_dissoc_rate_BigH2_H2g/SDIV(hmi.UV_Cont_rel2_Habing_spec_depth),
			hmi.H2_rate_destroy);
		fprintf(io,
			"\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
			hmi.HeatH2Dish_TH85,
			hmi.HeatH2Dexc_TH85,
			hmi.HeatH2Dish_BigH2, 
			hmi.HeatH2Dexc_BigH2,
			thermal.htot);
	}
	/* save h2 solomon */
	else if( (strcmp(chJOB , "H2so" ) == 0) && (strcmp(chTime,"LAST") != 0) )
	{
		/* remember as many as NSOL lines contributing to total Solomon process */
#		define NSOL 100
		double sum, one;
		long int jlosave[NSOL] , ivlosave[NSOL],
			iehisave[NSOL] ,jhisave[NSOL] , ivhisave[NSOL],
			nsave,
			ipOrdered[NSOL];
		int nFail;
		int i;
		realnum fsave[NSOL], wlsave[NSOL];
		/* Solomon process, and where it came from */
		fprintf(io,"%.5e\t%.3e", 
			/* depth in cm */
			radius.depth_mid_zone ,
			/* H2 destruction by Solomon process, big H2 model rate */
			hmi.H2_Solomon_dissoc_rate_BigH2_H2g +
				hmi.H2_Solomon_dissoc_rate_BigH2_H2s);
		sum = 0.;
		iElecLo = 0;
		/* find sum of all radiative exits from X into excited electronic states */
		if( h2.lgH2ON && hmi.lgBigH2_evaluated )
		{
			for( iElecHi=1; iElecHi<mole.n_h2_elec_states; ++iElecHi )
			{
				for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
				{
					for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
					{
						long int nv = h2.nVib_hi[iElecLo];
						for( iVibLo=0; iVibLo<=nv; ++iVibLo )
						{
							long nr = h2.nRot_hi[iElecLo][iVibLo];
							for( iRotLo=h2.Jlowest[iElecLo]; iRotLo<=nr; ++iRotLo )
							{
								/* only do if radiative transition exists */
								if( H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].ipCont > 0 )
								{
									one = H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Lo->Pop *
										H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->pump;
									sum += one;
								}
							}
						}
					}
				}
			}

			/* make sure it is safe to div by sum */
			sum = SDIV( sum );
			nsave = 0;
			/* now loop back over X and print all those which contribute more than FRAC of the total */
#			define FRAC	0.01
			for( iElecHi=1; iElecHi<mole.n_h2_elec_states; ++iElecHi )
			{
				for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
				{
					for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
					{
						long int nv = h2.nVib_hi[iElecLo];
						for( iVibLo=0; iVibLo<=nv; ++iVibLo )
						{
							long nr = h2.nRot_hi[iElecLo][iVibLo];
							for( iRotLo=h2.Jlowest[iElecLo]; iRotLo<=nr; ++iRotLo )
							{
								/* only do if radiative transition exists */
								if( H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].ipCont > 0 )
								{
									one = H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Lo->Pop *
										H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->pump;
									if( one/sum > FRAC && nsave<NSOL)
									{
										fsave[nsave] = (realnum)(one/sum);
										jlosave[nsave] = iRotLo;
										ivlosave[nsave] = iVibLo;
										jhisave[nsave] = iRotHi;
										ivhisave[nsave] = iVibHi;
										iehisave[nsave] = iElecHi;
										wlsave[nsave] = H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].WLAng;
										++nsave;
										/*fprintf(io,"\t%li\t%li\t%li\t%li\t%li\t%.3f", 
											iElecHi,iVibHi,iRotHi,iVibLo , iRotLo , one/sum );*/
									}
								}
							}
						}
					}
				}
			}/* iElecHi */
			/* now sort these into decreasing order */

			/* now sort by decreasing importance */
			/*spsort netlib routine to sort array returning sorted indices */
			spsort(
				/* input array to be sorted */
				fsave, 
				/* number of values in x */
				nsave, 
				/* permutation output array */
				ipOrdered, 
				/* flag saying what to do - 1 sorts into increasing order, not changing
				* the original routine */
				-1, 
				/* error condition, should be 0 */
				&nFail);

			/* print ratio of pumps to dissociations - this is 9:1 in TH85 */
			/*>>chng 05 jul 20, TE, save average energy in H2s and renormalization factors for H2g and H2s */
			/* >>chng 05 sep 16, TE, chng denominator to do g and s with proper dissoc rates */
			fprintf(io,"\t%.3f\t%.3f\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e",
				/* this is sum of photons and CRs */
				(sum + secondaries.csupra[ipHYDROGEN][0]*2.02f)/SDIV((hmi.H2_Solomon_dissoc_rate_BigH2_H2g * hmi.Hmolec[ipMH2g] +
					hmi.H2_Solomon_dissoc_rate_BigH2_H2s * hmi.Hmolec[ipMH2s]) ), 
				/* this is sum of photons and CRs */
				(sum + secondaries.csupra[ipHYDROGEN][0]*2.02f) /SDIV((hmi.H2_Solomon_dissoc_rate_BigH2_H2g *  hmi.H2g_BigH2 +
					hmi.H2_Solomon_dissoc_rate_BigH2_H2s * hmi.H2s_BigH2) ),
				hmi.H2_BigH2_H2g_av, hmi.H2_BigH2_H2s_av, 
				hmi.H2_chem_BigH2_H2g, hmi.H2_chem_BigH2_H2s,
				hmi.H2g_BigH2/SDIV(hmi.H2_total_BigH2), hmi.H2s_BigH2/SDIV(hmi.H2_total_BigH2)
				);
			for( i=0; i<nsave; ++i )
			{
				ip = ipOrdered[i];
				/*lint -e644 not init */
				fprintf(io,"\t%li\t%li\t%li\t%li\t%li\t%.3f\t%.3f", 
					iehisave[ip],ivhisave[ip],jhisave[ip],ivlosave[ip] , jlosave[ip] , fsave[ip] , wlsave[ip] );
				/*lint +e644 not init */
			}
			fprintf(io,"\n"); 
		}
		/*fprintf(io,"DEBUG tau\t%.3e\t%.3f\n",
			H2Lines[1][0][1][0][0][0].Emis->TauIn, H2Lines[1][0][1][0][0][0].WLAng); */
#		undef NSOL
	}

	else if( (strcmp(chJOB , "H2te" ) == 0) && (strcmp(chTime,"LAST") != 0) )
	{
		/* save h2 temperatures */
		double pop_ratio10,pop_ratio20,pop_ratio30,pop_ratio31,pop_ratio40;
		double T10,T20,T30,T31,T40;
		/* subscript"sum" denotes integrated quantities */
		double T10_sum,T20_sum,T30_sum,T31_sum,T40_sum;
		double pop_ratio10_sum,pop_ratio20_sum,pop_ratio30_sum,pop_ratio31_sum,pop_ratio40_sum;
		if( h2.lgH2ON && h2.nCallH2_this_zone )
		{
			double energyK = T1CM*(energy_wn[0][0][1] - energy_wn[0][0][0]);
			/* the ratio of H2_populations of J=1 to 0 */
			pop_ratio10 = H2_populations[0][0][1]/SDIV(H2_populations[0][0][0]);
			pop_ratio10_sum = H2_X_colden[0][1]/SDIV(H2_X_colden[0][0]);
			/* the corresponding temperature */
			T10 = -170.5/log(SDIV(pop_ratio10) * H2_stat[0][0][0]/H2_stat[0][0][1]);
			T10_sum = -170.5/log(SDIV(pop_ratio10_sum) * H2_stat[0][0][0]/H2_stat[0][0][1]);

			energyK = T1CM*(energy_wn[0][0][2] - energy_wn[0][0][0]);
			pop_ratio20 = H2_populations[0][0][2]/SDIV(H2_populations[0][0][0]);
			T20 = -energyK/log(SDIV(pop_ratio20) * H2_stat[0][0][0]/H2_stat[0][0][2]);

			pop_ratio20_sum = H2_X_colden[0][2]/SDIV(H2_X_colden[0][0]);
			T20_sum = -energyK/log(SDIV(pop_ratio20_sum) * H2_stat[0][0][0]/H2_stat[0][0][2]);

			energyK = T1CM*(energy_wn[0][0][3] - energy_wn[0][0][0]);
			pop_ratio30 = H2_populations[0][0][3]/SDIV(H2_populations[0][0][0]);
			T30 = -energyK/log(SDIV(pop_ratio30) * H2_stat[0][0][0]/H2_stat[0][0][3]);

			pop_ratio30_sum = H2_X_colden[0][3]/SDIV(H2_X_colden[0][0]);
			T30_sum = -energyK/log(SDIV(pop_ratio30_sum) * H2_stat[0][0][0]/H2_stat[0][0][3]);

			energyK = T1CM*(energy_wn[0][0][3] - energy_wn[0][0][1]);
			pop_ratio31 = H2_populations[0][0][3]/SDIV(H2_populations[0][0][1]);
			T31 = -energyK/log(SDIV(pop_ratio31) * H2_stat[0][0][1]/H2_stat[0][0][3]);

			pop_ratio31_sum = H2_X_colden[0][3]/SDIV(H2_X_colden[0][1]);
			T31_sum = -energyK/log(SDIV(pop_ratio31_sum) * H2_stat[0][0][1]/H2_stat[0][0][3]);

			energyK = T1CM*(energy_wn[0][0][4] - energy_wn[0][0][0]);
			pop_ratio40 = H2_populations[0][0][4]/SDIV(H2_populations[0][0][0]);
			T40 = -energyK/log(SDIV(pop_ratio40) * H2_stat[0][0][0]/H2_stat[0][0][4]);

			pop_ratio40_sum = H2_X_colden[0][4]/SDIV(H2_X_colden[0][0]);
			T40_sum = -energyK/log(SDIV(pop_ratio40_sum) * H2_stat[0][0][0]/H2_stat[0][0][4]);
		}
		else
		{
			pop_ratio10 = 0.;
			pop_ratio10_sum = 0.;
			T10 = 0.;
			T20 = 0.;
			T30 = 0.;
			T31 = 0.;
			T40 = 0.;
			T10_sum = 0.;
			T20_sum = 0.;
			T30_sum = 0.;
			T31_sum = 0.;
			T40_sum = 0.;
		}

		/* various temperatures for neutral/molecular gas */
		fprintf( io, 
			"%.5e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n" ,
			/* depth in cm */
			radius.depth_mid_zone ,
			/* total H2 fraction */
			hmi.H2_total/dense.gas_phase[ipHYDROGEN] ,
			/* ratio of H2_populations of 1 to 0 only */
			pop_ratio10 ,
			/* sum of all ortho and para */
			h2.ortho_density / SDIV(h2.para_density),
			T10,T20,T30,T31,T40,
			phycon.te ,
			hyperfine.Tspin21cm,T10_sum,T20_sum,T30_sum,T31_sum,T40_sum  );
	}
	else if( (strcmp(chJOB , "H2ln" ) == 0) && (strcmp(chTime,"LAST") == 0) )
	{
		/* save H2 lines - output the full emission-line spectrum */
		double thresh;
		double renorm;
		/* first test, is H2 turned on?  Second test, have lines arrays
		 * been set up - nsum is negative if abort occurs before lines
		 * are set up */
		if( h2.lgH2ON && LineSave.nsum > 0)
		{
			ASSERT( LineSave.ipNormWavL >= 0 );
			/* get the normalization line */
			if( LineSv[LineSave.ipNormWavL].SumLine[0] > SMALLFLOAT )
				renorm = LineSave.ScaleNormLine/
				LineSv[LineSave.ipNormWavL].SumLine[0];
			else
				renorm = 1.;

			if( renorm > SMALLFLOAT )
			{
				/* this is threshold for faintest line, normally 0, set with 
				 * number on save H2 command */
				thresh = thresh_punline_h2/(realnum)renorm;
			}
			else
				thresh = 0.f;

			/* save H2 line intensities at end of iteration 
			 * h2.nElecLevelOutput is electronic level with 1 for ground, so this loop is < h2.nElecLevelOutput */
			for( iElecHi=0; iElecHi < h2.nElecLevelOutput; ++iElecHi )
			{
				for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
				{
					for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
					{
						/* now the lower levels */
						/* NB - X is the only lower level considered here, since we are only 
						* concerned with excited electronic levels as a photodissociation process
						* code exists to relax this assumption - simply change following to iElecHi */
						long int lim_elec_lo = 0;
						for( iElecLo=0; iElecLo<=lim_elec_lo; ++iElecLo )
						{
							long int nv = h2.nVib_hi[iElecLo];
							if( iElecLo==iElecHi )
								nv = iVibHi;
							for( iVibLo=0; iVibLo<=nv; ++iVibLo )
							{
								long nr = h2.nRot_hi[iElecLo][iVibLo];
								if( iElecLo==iElecHi && iVibHi==iVibLo )
									nr = iRotHi-1;
								for( iRotLo=h2.Jlowest[iElecLo]; iRotLo<nr; ++iRotLo )
								{
									if( H2_SaveLine[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo] > thresh )
									{
										/* air wavelength in microns */
										/* WLAng contains correction for index of refraction of air */
										double wl = H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].WLAng/1e4;
										/*ASSERT( abs(iRotHi-iRotLo)<=2 );*/

										fprintf(io, "%li-%li %c(%li)", 
											iVibHi , 
											iVibLo ,
											chMolBranch( iRotHi , iRotLo ) ,
											iRotLo );
										fprintf( io, "\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld", 
											iElecHi , iVibHi , iRotHi , iElecLo , iVibLo , iRotLo);
										/* WLAng contains correction for index of refraction of air */
										fprintf( io, "\t%.7f\t", wl );
										/*prt_wl print floating wavelength in Angstroms, in output format */
										prt_wl( io , H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].WLAng );
										/* the log of the line intensity or luminosity */
										fprintf( io, "\t%.3f\t%.3e", 
											log10(MAX2(1e-37,H2_SaveLine[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo])) + radius.Conv2PrtInten, 
											H2_SaveLine[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo]*renorm );
										/* excitation energy of upper level in K */
										fprintf( io, "\t%.3f", energy_wn[iElecHi][iVibHi][iRotHi]*T1CM );
										/* the product g_hi h nu * Aul */
										ASSERT( H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].ipCont > 0 );
										fprintf( io, "\t%.3e", H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->Aul* 
											H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].EnergyErg * 
											H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Hi->g);
										fprintf( io, "\n");
									}
								}
							}
						}
					}
				}
			}
		}
	}
	else if( (strcmp(chJOB , "H2sp" ) == 0)  )
	{
		iVib = 0;
		iRot = 0;
#		if 0
		/* save h2 special */
		fprintf(io,"%.4e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
			radius.depth_mid_zone , 
			H2_populations[0][iVib][iRot] ,
			radius.depth_mid_zone * H2_populations[0][iVib][iRot] ,
			H2_rad_rate_out[0][iVib][iRot] ,
			H2_rad_rate_in[iVib][iRot] ,
			H2_col_rate_out_ old[iVib][iRot] ,
			H2_col_rate_in_ old[iVib][iRot] );
#		endif
		fprintf(io,"%.4e\t%.2e\t%.2e\t%.2e\t%.2e\n",
			radius.depth_mid_zone , 
			H2_populations[0][iVib][iRot] ,
			H2Lines[1][1][1][0][iVib][iRot].Emis->pump,
			H2Lines[1][1][1][0][iVib][iRot].Emis->TauIn,
			H2Lines[1][1][1][0][iVib][iRot].Emis->TauCon);
	}
	return;
}
 /*cdH2_Line determines intensity and luminosity of and H2 line.  The first
  * six arguments give the upper and lower quantum designation of the levels.
  * The function returns 1 if we found the line, 
  * and false==0 if we did not find the line because ohoto-para transition
  * or upper level has lower energy than lower level  */
long int cdH2_Line(
	  /* indices for the upper level */
	  long int iElecHi, 
	  long int iVibHi ,
	  long int iRotHi ,
	  /* indices for lower level */
	  long int iElecLo, 
	  long int iVibLo ,
	  long int iRotLo ,
	  /* linear intensity relative to normalization line*/
	  double *relint, 
	  /* log of luminosity or intensity of line */
	  double *absint )
{

	DEBUG_ENTRY( "cdH2_Line()" );

	/* these will be return values if we can't find the line */
	*relint = 0.;
	*absint = 0.;

	/* for now both electronic levels must be zero */
	if( iElecHi!=0 || iElecLo!=0 )
	{
		return 0;
	}

	/* check that energy of first level is higher than energy of second level */
	if( energy_wn[iElecHi][iVibHi][iRotHi] < energy_wn[iElecLo][iVibHi][iRotHi] )
	{
		return 0;
	}

	/* check that ortho-para does not change */
	if( H2_lgOrtho[iElecHi][iVibHi][iRotHi] - H2_lgOrtho[iElecLo][iVibLo][iRotLo] != 0 )
	{
		return 0;
	}

	/* exit if lines does not exist */
	if( !lgH2_line_exists[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo] )
	{
		return 0;
	}

	ASSERT( LineSave.ipNormWavL >= 0 );
	/* does the normalization line have a positive intensity*/
	if( LineSv[LineSave.ipNormWavL].SumLine[0] > 0. )
	{
		*relint = H2_SaveLine[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo]/
			LineSv[LineSave.ipNormWavL].SumLine[0] * LineSave.ScaleNormLine;
	}
	else
	{
		*relint = 0.;
	}

	/* return log of line intensity if it is positive */
	if( H2_SaveLine[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo] > 0. )
	{
		*absint = log10(H2_SaveLine[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo]) + 
			radius.Conv2PrtInten;
	}
	else
	{
		/* line intensity is actually zero, return small number */
		*absint = -37.;
	}
	/* this indicates success */
	return 1;
}
