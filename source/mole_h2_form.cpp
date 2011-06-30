/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*mole_H2_form find state specific rates grains and H- form H2 */
#include "cddefines.h" 
#include "grainvar.h" 
#include "phycon.h" 
#include "hmi.h" 
#include "dense.h" 
#include "h2.h" 
#include "h2_priv.h" 

/*mole_H2_form find state specific rates grains and H- form H2 */
void mole_H2_form( void )
{

	long int iVib, iRot;
	long int ipT;
	double rate, oldrate ,
		frac_lo , frac_hi; 

	DEBUG_ENTRY( "mole_H2_form()" );

	/* rate of entry into X from H- and formation on grain surfaces 
	 * will one of several distribution functions derived elsewhere
	 * first zero out formation rates and rates others collide into particular level */
	for( iVib = 0; iVib <= h2.nVib_hi[0]; ++iVib )
	{
		for( iRot=h2.Jlowest[0]; iRot<=h2.nRot_hi[0][iVib]; ++iRot )
		{
			/* this will be the rate formation (s-1) of H2 due to
			 * both formation on grain surfaces and the H minus route,
			 * also H2+ + H => H2 + H+ into one vJ level */
			/* units cm-3 s-1 */
			H2_X_formation[iVib][iRot] = 0.;
			H2_X_Hmin_back[iVib][iRot] = 0.;
		}
	}

	/* loop over all grain types, finding total formation rate into each ro-vib level,
	 * also keeps trace of total formation into H2 ground and star, as defined by Tielens & Hollenbach,
	 * these are used in the H molecular network */
	hmi.H2_forms_grains = 0.;
	hmi.H2star_forms_grains = 0.;
	/* >>chng 02 oct 08, resolved grain types */
	for( size_t nd=0; nd < gv.bin.size(); ++nd )
	{
		int ipH2 = (int)gv.which_H2distr[gv.bin[nd]->matType];
		for( iVib = 0; iVib <= h2.nVib_hi[0]; ++iVib )
		{
			for( iRot=h2.Jlowest[0]; iRot<=h2.nRot_hi[0][iVib]; ++iRot )
			{
				/* >>chng 02 nov 14, changed indexing into H2_X_grain_formation_distribution and gv.bin, PvH */
				realnum one = 
					/* H2_X_grain_formation_distribution is normalized to a summed total of unity */
					H2_X_grain_formation_distribution[ipH2][iVib][iRot] * 
					/* units of following are s-1 */
					(realnum)gv.bin[nd]->rate_h2_form_grains_used;
				/* final units are s-1 */
				/* units cm-3 s-1 */
				/* >>chng 04 may 05, added atomic hydrogen density, units cm-3 s-1 */
				H2_X_formation[iVib][iRot] += one*dense.xIonDense[ipHYDROGEN][0];

				/* save rates for formation into "H2" and "H2*" in the chemical
				 * network - it resolves the H2 into two species, as in 
				 * Hollenbach / Tielens work - these rates will be used in the
				 * chemistry solver to get H2 and H2* densities */
				if( energy_wn[0][iVib][iRot] < ENERGY_H2_STAR )
				{   
					/*  unit s-1, h2 means h2g*/ 
					hmi.H2_forms_grains += one;
				}
				else
				{
					hmi.H2star_forms_grains += one;
				}
			}
		}
	}

	/* formation of H2 in excited states from H- H minus */
	/* >>chng 02 oct 17, temp dependent fit to rate, updated reference,
	 * about 40% larger than before */
	/* rate in has units cm-3 s-1 */
	rate = hmi.Hmolec[ipMHm] * dense.xIonDense[ipHYDROGEN][0] * hmi.assoc_detach;
	/*rate = hmi.hminus*1.35e-9f;*/
	/* convert to dimensionless factors that add to unity */
	/* >>chng 02 oct 17, use proper distribution function */
	hmi.H2star_forms_hminus = 0.;
	hmi.H2_forms_hminus = 0.;
	oldrate = 0.;
	/* which temperature point to use? */
	if( phycon.alogte<=1. )
	{
		ipT = 0;
		frac_lo = 1.;
		frac_hi = 0.;
	}
	else if( phycon.alogte>=4. )
	{
		ipT = nTE_HMINUS-2;
		frac_lo = 0.;
		frac_hi = 1.;
	}
	else
	{
		/* find the temp */
		for( ipT=0; ipT<nTE_HMINUS-1; ++ipT )
		{
			if( H2_te_hminus[ipT+1]>phycon.alogte )
				break;
		}
		frac_hi = (phycon.alogte-H2_te_hminus[ipT])/(H2_te_hminus[ipT+1]-H2_te_hminus[ipT]);
		frac_lo = 1.-frac_hi;
	}

	/* we now know how to interpolate, now fill in H- formation sites */
	for( iVib=0; iVib<=h2.nVib_hi[0]; ++iVib )
	{
		for( iRot=h2.Jlowest[0]; iRot<=h2.nRot_hi[0][iVib]; ++iRot )
		{
			/* the temperature-interpolated distribution function, adds up to unity, 
			 * dimensionless */
			double rate_interp =
				frac_lo*H2_X_hminus_formation_distribution[ipT][iVib][iRot] +
				frac_hi*H2_X_hminus_formation_distribution[ipT+1][iVib][iRot];

			/* above rate was set, had dimensions cm-3 s-1 
			 * rate is product of parent densities and total formation rate */
			realnum one = (realnum)(rate * rate_interp);

			/* save this rate [cm3 s-1] for back reaction in levels solver for v,J */
			H2_X_Hmin_back[iVib][iRot] = (realnum)(rate_interp * hmi.assoc_detach);

			/* units cm-3 s-1 */
			H2_X_formation[iVib][iRot] += one;

			oldrate += rate_interp;

			/* save rates to pass back into molecule network */
			if( energy_wn[0][iVib][iRot] < ENERGY_H2_STAR )
			{	
				/*  unit cm-3s-1, h2 means h2g*/
				hmi.H2_forms_hminus += one;
			}
			else
			{
				hmi.H2star_forms_hminus += one;
			}
		}
	}
	/* confirm that shape function is normalized correctly */
	ASSERT( fabs(1.-oldrate)<1e-4 );

	/* >>chng 03 feb 10, add this population process */
	/* H2+ + H => H2 + H+,
	 * >>refer	H2	population	Krstic, P.S., preprint 
	 * all goes into v=4 but no J information, assume into J = 0 */
	/* >>chng 04 may 05, add density at end */
	rate = hmi.bh2h2p * hmi.Hmolec[ipMH2p] * dense.xIonDense[ipHYDROGEN][0];
	iVib = 4;
	iRot = 0;
	/* units cm-3 s-1 */
	H2_X_formation[iVib][iRot] += (realnum)rate;

	return;
}
