/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*H2_ContPoint set the ipCont struc element for the H2 molecule, called by ContCreatePointers */
/*H2_Accel radiative acceleration due to H2 */
/*H2_RadPress rad pressure due to h2 lines called in PresTotCurrent */
/*H2_InterEnergy internal energy of H2 called in PresTotCurrent */
/*H2_RT_diffuse do emission from H2 - called from RT_diffuse */
/*H2_itrzn - average number of H2 pop evaluations per zone */
/*H2_RTMake do RT for H2 - called from RT_line_all */
/*H2_RT_tau_inc increment optical depth for the H2 molecule, called from RT_tau_inc */
/*H2_LineZero initialize optical depths in H2, called from RT_tau_init */
/*H2_RT_tau_reset the large H2 molecule, called from RT_tau_reset */
/*H2_Colden maintain H2 column densities within X */
/*H2_LevelPops do level H2_populations for H2, called by Hydrogenic */
/*H2_Level_low_matrix evaluate CO rotation cooling */
/*H2_cooling evaluate cooling and heating due to H2 molecule */
/*H2_X_coll_rate_evaluate find collisional rates within X */
/*cdH2_colden return column density in H2, negative -1 if cannot find state,
 * header is cddrive */
/*H2_DR choose next zone thickness based on H2 big molecule */
/* turn this flag on to do minimal debug print of pops */
#define	PRT_POPS	false
/* this is limit to number of loops over H2 pops before giving up */
#define	LIM_H2_POP_LOOP	100
/* this is a typical dissociation cross section (cm2) for H2 + Hnu -> 2H + ke */
/* >>chng 05 may 11, had been 2.5e-19 */
#define	H2_DISS_ALLISON_DALGARNO	6e-19f
#include "cddefines.h" 
#include "cddrive.h" 
#include "physconst.h" 
#include "taulines.h" 
#include "atoms.h" 
#include "conv.h" 
#include "secondaries.h" 
#include "pressure.h" 
#include "trace.h" 
#include "hmi.h" 
#include "hextra.h" 
#include "rt.h" 
#include "radius.h" 
#include "ipoint.h" 
#include "phycon.h" 
#include "thermal.h" 
#include "dense.h" 
#include "rfield.h" 
#include "lines_service.h" 
#include "mole.h"
#include "h2.h"
#include "h2_priv.h"

/* this counts how many times we go through the H2 level H2_populations loop */
static long int loop_h2_pops;

realnum H2_te_hminus[nTE_HMINUS] = {10.,30.,100.,300.,1000.,3000.,10000.};

/* this will contain a vector for collisions within the X ground electronic state,
 * CollRateFit[vib_up][rot_up][vib_lo][rot_lo][coll_type][3] */
static realnum collider_density[N_X_COLLIDER];
static realnum collider_density_total_not_H2;

/* this integer is added to rotation quantum number J for the test of whether
 * a particular J state is ortho or para - the state is ortho if J+below is odd,
 * and para if J+below is even */
int H2_nRot_add_ortho_para[N_H2_ELEC] = {0 , 1 , 1 , 0, 1, 1 , 0};

/* dissociation energies (cm-1) for each electronic state, from 
 * >>refer	H2	energies	Sharp, T. E., 1971, Atomic Data, 2, 119 
 * energy for H_2 + hnu -> 2H,
 * note that energies for excited states are in the Lyman continuum
 * (1 Ryd = 109670 cm-1) */
/* >>chng 02 oct 08, improved energies */
double H2_DissocEnergies[N_H2_ELEC] = 
{ 36118.11, 118375.6, 118375.6, 118375.6, 118375.6,133608.6,133608.6 };
/* original values
 { 36113., 118372., 118372., 118372., 118372.,0.,0. };*/

double exp_disoc; 

/*H2_X_coll_rate_evaluate find collisional rates within X - 
 * this is one time upon entry into H2_LevelPops */
STATIC void H2_X_coll_rate_evaluate( void )
{
	long int nColl,
		ipHi ,
		ipLo,
		ip,
		iVibHi,
		iRotHi,
		iVibLo,
		iRotLo;
	realnum colldown;
	double exph2_big, 
		exph2_big_0,
		rel_pop_LTE_H2_big;

	DEBUG_ENTRY( "H2_X_coll_rate_evaluate()" );

	/* set collider density 
	 * the colliders are:
	 * [0] = H
	 * [1], [5] = He (old and new cs data)
	 * [2] = H2 ortho
	 * [3] = H2 para
	 * [4] = H+ + H3+ */
	/* atomic hydrogen */
	collider_density[0] = dense.xIonDense[ipHYDROGEN][0];
	/* all ortho h2 */
	/* He - H2 */
	collider_density[1] = dense.xIonDense[ipHELIUM][0];
	/* H2 - H2(ortho) */
	collider_density[2] = (realnum)h2.ortho_density;
	/* all para H2 */
	collider_density[3] = (realnum)h2.para_density;
	/* protons - ionized hydrogen */
	collider_density[4] = dense.xIonDense[ipHYDROGEN][1];
	/* H3+ - assume that H3+ has same rates as proton */
	collider_density[4] += hmi.Hmolec[ipMH3p];

	ASSERT( fp_equal(hmi.H2_total ,collider_density[2]+collider_density[3]) );

	/* this is total density of all colliders, is only used for collisional dissociation
	 * rates for H2 are not included here, will be added separately*/
	collider_density_total_not_H2 = collider_density[0] + 
		collider_density[1] + collider_density[4] + 
		(realnum)dense.eden;

	if( mole.nH2_TRACE >= mole.nH2_trace_full )
	{
		fprintf(ioQQQ," Collider densities are:");
		for( nColl=0; nColl<N_X_COLLIDER; ++nColl )
		{
			fprintf(ioQQQ,"\t%.3e", collider_density[nColl]);
		}
		fprintf(ioQQQ,"\n");
	}

	for( ipHi=0; ipHi<nLevels_per_elec[0]; ++ipHi )
	{
		H2_X_source[ipHi] = 0.;
		H2_X_sink[ipHi] = 0.;
	}
	H2_X_coll_rate.zero();

	/*>>chng 05 sep 18, GS, determine LTE populations for H2+ e = H- + H*/
	exp_disoc = sexp(H2_DissocEnergies[0]/phycon.te_wn);
	exph2_big_0 = exp_disoc/SDIV(H2_Boltzmann[0][0][0]);
	/*>>chng 05 oct 17, GS,  (2*m_e/m_p)^3/2 = 3.634e-5 */
	rel_pop_LTE_H2_big =SAHA/SDIV((phycon.te32*exph2_big_0))*(H2_stat[0][0][0]/(2.*2.))*3.634e-5;
	/* now find rates for all collisions within X */
	/* this is special J=1 to J=0 collision, which is only fast at
	 * very low grain temperatures 
	H2_X_coll_rate[1][0] = 
		(realnum)(hmi.rate_grain_h2_J1_to_J0);*/

	/* count formation from grains and H- as a collisional formation process */
	/* cm-3 s-1, evaluated in mole_H2_form */
	H2_X_source[0] += H2_X_formation[0][0];
	/*>>chng 05 sep 18, GS, H2+ e = H- + H, unit s-1
	 * H2_X_Hmin_back[iVib][iRot] is resolved formation rate, units cm3 s-1 */
	H2_X_sink[0] +=  (realnum)(H2_X_Hmin_back[0][0]*hmi.rel_pop_LTE_Hmin/
		SDIV(rel_pop_LTE_H2_big)*dense.eden);
	/* this represents collisional dissociation into continuum of X,
	 * rates are just guesses */
	H2_X_sink[0] += collider_density_total_not_H2 *
		H2_coll_dissoc_rate_coef[0][0] * mole.lgColl_deexec_Calc;
	/*>>chng 05 jul 20, GS, collisional dissociation with H2g and H2s are added here*/
	H2_X_sink[0] +=  hmi.H2_total*
		H2_coll_dissoc_rate_coef_H2[0][0] * mole.lgColl_deexec_Calc;

	/* this is cosmic ray or secondary photodissociation into mainly H2+ */
	H2_X_sink[0] += secondaries.csupra[ipHYDROGEN][0]*2.02f * mole.lgColl_deexec_Calc;

	/*>>chng 05 sep 26, GS, H2 + CRP -> H+ + H-  */
	H2_X_sink[0] += (realnum)(3.9e-21 * hextra.cryden_ov_background * co.lgUMISTrates);

	/*>>chng 05 sep 26, GS, H2 + CRP -> H+ + H + e */
	H2_X_sink[0] += (realnum)(2.2e-19 * hextra.cryden_ov_background * co.lgUMISTrates);

	/*>>chng 05 sep 26, GS, H2 + CR -> H+ + H + e */  
	H2_X_sink[0] += (realnum)(secondaries.csupra[ipHYDROGEN][0]*0.0478);

	/* collisional dissociation by non-thermal electrons 
	 * and cosmic rays to the triplet state of H2 (a,b,c ) from 
	 *>>>refer Dalgarno, Yan, & Liu 1999 ApJs
	 * rates depend on energy of secondary electrons, this is an estimate 
	 * from fig 4b around incident electron energy 15 to 20 eV. 
	 * The cross section of state b is divided by the cross-section of Lya 
	 * (the ratio is ~ 10)and then multiplied by the cosmic ray excitation  
	 * rate of Lya (secondaries.x12tot) the triplet state of H2 (b ) 
	 *>>chng 07 apr 08, GS from 3 to 10 after study of Dalgarno et al */
	H2_X_sink[0] += 10.f*secondaries.x12tot; 

	/*>>chng 05 jul 07, GS, ionization by hard photons are added to be consistent with chemistry-network */
	H2_X_sink[0] += (realnum)hmi.H2_photoionize_rate; 

	/* rate (s-1) out of this level, this is dissociation into the continuum of singlet B and C states*/
	H2_X_sink[0] += rfield.flux_accum[H2_ipPhoto[0][0]-1]*H2_DISS_ALLISON_DALGARNO;

	/*>>chng 05 sept 28, GS, H2 + H+ = H2+ + H; note that H2 is destroyed from v=0, J=0,
	 * forward and backward reactions are not in detailed balance, astro-ph/0404288, final unit s-1
	 * ihmi.rh2h2p is a rate coefficient, cm3 s-1, and needs a density to become
	 * the sink inverse lifetime.  The process is H2(v=0, J=0) + H+ -> H0 + H2+ */
	H2_X_sink[0] += (realnum)(hmi.rh2h2p*dense.xIonDense[ipHYDROGEN][1]);

	/* now find total coll rates - this loop is over the energy-sorted levels */
	ipHi = nLevels_per_elec[0];
	while( (--ipHi) > 0 )
	{
		/* array of energy sorted indices within X */
		ip = H2_ipX_ener_sort[ipHi];
		iVibHi = ipVib_H2_energy_sort[ip];
		iRotHi = ipRot_H2_energy_sort[ip];

		/*>>chng 05 sep 18, GS, determine LTE populations for H2+ e = H- + H, unit s-1*/
		exph2_big = exp_disoc/SDIV(H2_Boltzmann[0][iVibHi][iRotHi]);
		/* >>chng 05 oct 13, replace 1 with stat wght of level */
		/*>>chng 05 oct 17, GS,  (2*m_e/m_p)^3/2 = 3.634e-5 */
		rel_pop_LTE_H2_big = SAHA/SDIV((phycon.te32*exph2_big))*(H2_stat[0][iVibHi][iRotHi]/(2.*2.))*3.634e-5;
		/* formation */
		/* count formation from grains and H- as a collisional formation process */
		/* cm-3 s-1, evaluated in mole_H2_form */
		H2_X_source[ipHi] += H2_X_formation[iVibHi][iRotHi];
		/*>>chng 05 sep 18, GS, H2+ e = H- + H*, H2_X_Hmin_back has units cm3 s-1 so final units s-1 */
		H2_X_sink[ipHi] += (realnum)(H2_X_Hmin_back[iVibHi][iRotHi]*hmi.rel_pop_LTE_Hmin/
			SDIV(rel_pop_LTE_H2_big)*dense.eden);
		/* this represents collisional dissociation into continuum of X,
		 * rates are just guesses */
		H2_X_sink[ipHi] += collider_density_total_not_H2 *
			H2_coll_dissoc_rate_coef[iVibHi][iRotHi] * mole.lgColl_deexec_Calc;

		/*>>chng 05 jul 20, GS, collisional dissociation with H2g and H2s are added here*/
		H2_X_sink[ipHi] +=  hmi.H2_total*
			H2_coll_dissoc_rate_coef_H2[iVibHi][iRotHi] * mole.lgColl_deexec_Calc;

		/* this is cosmic ray or secondary photodissociation into mainly H2+ */
		H2_X_sink[ipHi] += secondaries.csupra[ipHYDROGEN][0]*2.02f * mole.lgColl_deexec_Calc;

		/*>>chng 05 sep 26, GS, H2 + CRP -> H+ + H-  */
		H2_X_sink[ipHi] += (realnum)(3.9e-21 * hextra.cryden_ov_background * co.lgUMISTrates);

		/*>>chng 05 sep 26, GS, H2 + CRP -> H+ + H + e */
		H2_X_sink[ipHi] += (realnum)(2.2e-19 * hextra.cryden_ov_background * co.lgUMISTrates);

		/*>>chng 05 sep 26, GS, H2 + CR -> H+ + H + e */  
		H2_X_sink[ipHi] += (realnum)(secondaries.csupra[ipHYDROGEN][0]*0.0478);

		/* collisional dissociation by non-thermal electrons 
		* and cosmic rays to the triplet state of H2 (a,b,c ) from 
		*>>>refer Dalgarno, Yan, & Liu 1999 ApJs
		* rates depend on energy of secondary electrons, this is an estimate 
		*from fig 4b around incident electron energy 15 to 20 eV. 
		* The cross section of state b is divided by the cross-section of Lya 
		* (the ratio is ~ 10)and then multiplied by the cosmic ray excitation  
		* rate of Lya (secondaries.x12tot) the triplet state of H2 (b ) 
		*>>chng 07 apr 08, GS from 3 to 10 after study of Dalgarno et al */
		H2_X_sink[ipHi] += 10.f*secondaries.x12tot;

		/*>>chng 05 jul 07, GS, ionization by hard photons are added to be consistent with chemistry-network */
		H2_X_sink[ipHi] += (realnum)hmi.H2_photoionize_rate;

		/* rate (s-1) out of this level */
		H2_X_sink[ipHi] += rfield.flux_accum[H2_ipPhoto[iVibHi][iRotHi]-1]*H2_DISS_ALLISON_DALGARNO;

		if( mole.lgColl_deexec_Calc )
		{
			/* excitation within X due to thermal particles */
			for( ipLo=0; ipLo<ipHi; ++ipLo )
			{
				/* array of energy sorted indices within X */
				ip = H2_ipX_ener_sort[ipLo];
				iVibLo = ipVib_H2_energy_sort[ip];
				iRotLo = ipRot_H2_energy_sort[ip];

				/* collisional interactions with upper levels within X */
				colldown = 0.;
				mr5ci H2CollRate = H2_CollRate.begin(iVibHi,iRotHi,iVibLo,iRotLo);
				for( nColl=0; nColl<N_X_COLLIDER; ++nColl )
				{
					/* downward collision rate, units s-1 */
					colldown += H2CollRate[nColl]*collider_density[nColl];
					ASSERT( H2CollRate[nColl]*collider_density[nColl] >= 0. );
				}
				/* rate in from upper level, units cm-3 s-1 */
				H2_X_coll_rate[ipHi][ipLo] += colldown;
			}/* end loop over ipLo */
		}
	}/* end loop over ipHi */
	return;
}

/*H2_itrzn - average number of H2 pop evaluations per zone */
double H2_itrzn( void )
{
	if( h2.lgH2ON && nH2_zone>0 )
	{
		return( (double)nH2_pops / (double)nH2_zone );
	}
	else
	{
		return 0.;
	}
}

/* set the ipCont struc element for the H2 molecule, called by ContCreatePointers */
void H2_ContPoint( void )
{
	long int iElecHi , iElecLo , iVibHi , iVibLo , iRotHi , iRotLo;

	if( !h2.lgH2ON )
		return;

	DEBUG_ENTRY( "H2_ContPoint()" );

	/* set array index for line energy within continuum array */
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
					/* want to include all vibration states in lower level if different electronic level,
					 * but only lower vibration levels if same electronic level */
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
							/* >>chng 05 feb 07, use lgH2_line_exists */
							if( lgH2_line_exists[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo] )
							{
								ASSERT( H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->Aul > 0. );
								H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].ipCont = 
									ipLineEnergy(
									H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].EnergyWN * WAVNRYD , 
									"H2  " , 0 );
								H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->ipFine = 
									ipFineCont(
									H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].EnergyWN * WAVNRYD );
							}
						}
					}
				}
			}
		}
	}
	return;
}

/* ===================================================================== */
/* radiative acceleration due to H2 called in rt_line_driving */
double H2_Accel(void)
{
	long int iElecHi , iElecLo , iVibHi , iVibLo , iRotHi , iRotLo;
	double h2_drive;

	/* >>chng 05 jan 26, pops now set to LTE for small abundance case, so do this */
	if( !h2.lgH2ON /*|| !h2.nCallH2_this_zone*/ )
		return 0.;

	DEBUG_ENTRY( "H2_Accel()" );

	/* this routine computes the line driven radiative acceleration
	 * due to H2 molecule*/

	h2_drive = 0.;
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
					/* want to include all vibration states in lower level if different electronic level,
					 * but only lower vibration levels if same electronic level */
					long int nv = h2.nVib_hi[iElecLo];
					if( iElecLo==iElecHi )
						nv = iVibHi;
					for( iVibLo=0; iVibLo<=nv; ++iVibLo )
					{
						long nr = h2.nRot_hi[iElecLo][iVibLo];
						if( iElecLo==iElecHi && iVibHi==iVibLo )
							nr = iRotHi-1;

						mb6ci lgH2le = lgH2_line_exists.ptr(iElecHi,iVibHi,iRotHi,iElecLo,iVibLo,h2.Jlowest[iElecLo]);
						mt6ci H2L = H2Lines.ptr(iElecHi,iVibHi,iRotHi,iElecLo,iVibLo,h2.Jlowest[iElecLo]);
						for( iRotLo=h2.Jlowest[iElecLo]; iRotLo<=nr; ++iRotLo )
						{
							/* >>chng 03 feb 14, from !=0 to >0 */
							/* >>chng 05 feb 07, use lgH2_line_exists */
							if( *lgH2le++ )
							{
								ASSERT( H2L->ipCont > 0 );
								h2_drive += H2L->Emis->pump*H2L->EnergyErg*H2L->Emis->PopOpc;
							}
							++H2L;
						}
					}
				}
			}
		}
	}
	return h2_drive;
}

/* ===================================================================== */
/* rad pressure due to H2 lines called in PresTotCurrent */
double H2_RadPress(void)
{
	long int iElecHi , iElecLo , iVibHi , iVibLo , iRotHi , iRotLo;
	double press;

	/* will be used to check on size of opacity, was capped at this value */
	realnum smallfloat=SMALLFLOAT*10.f;

	/* radiation pressure sum is expensive - do not evaluate if we did not
	 * bother evaluating large molecule */
	if( !h2.lgH2ON || !h2.nCallH2_this_zone )
		return 0.;

	DEBUG_ENTRY( "H2_RadPress()" );

	press = 0.;
	/* loop over all possible lines */
	for( iElecHi=0; iElecHi<mole.n_h2_elec_states; ++iElecHi )
	{
		for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
		{
			for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
			{
				/* now the lower levels - X is the only lower level considered 
				 * here, we only do excited electronic levels as a 
				 * photodissociation process - code exists to relax this 
				 * assumption - simply change following to iElecHi */
				long int lim_elec_lo = 0;
				for( iElecLo=0; iElecLo<=lim_elec_lo; ++iElecLo )
				{
					/* want to include all vibration states in lower level if different electronic level,
					 * but only lower vibration levels if same electronic level */
					long int nv = h2.nVib_hi[iElecLo];
					if( iElecLo==iElecHi )
						nv = iVibHi;
					for( iVibLo=0; iVibLo<=nv; ++iVibLo )
					{
						long nr = h2.nRot_hi[iElecLo][iVibLo];
						if( iElecLo==iElecHi && iVibHi==iVibLo )
							nr = iRotHi-1;

						mb6ci lgH2le = lgH2_line_exists.ptr(iElecHi,iVibHi,iRotHi,iElecLo,iVibLo,h2.Jlowest[iElecLo]);
						mt6ci H2L = H2Lines.ptr(iElecHi,iVibHi,iRotHi,iElecLo,iVibLo,h2.Jlowest[iElecLo]);
						for( iRotLo=h2.Jlowest[iElecLo]; iRotLo<=nr; ++iRotLo )
						{
							if( *lgH2le++ )
							{
								ASSERT( H2L->ipCont > 0 );

								if( H2L->Hi->Pop > smallfloat && H2L->Emis->PopOpc > smallfloat )
								{
									double RadPres1 =  PressureRadiationLine( &(*H2L), GetDopplerWidth(2.f*dense.AtomicWeight[ipHYDROGEN])  );
									press += RadPres1;
								}
							}
							++H2L;
						}
					}
				}
			}
		}
	}

	if(mole.nH2_TRACE >= mole.nH2_trace_full) 
		fprintf(ioQQQ,
		"  H2_RadPress returns, radiation pressure is %.2e\n", 
		press );
	return press;
}

#if 0
/* ===================== */
/* internal energy of H2 */
double H2_InterEnergy(void)
{
	double energy;

	/* >>chng 05 jan 26, pops now set to LTE for small abundance case, so do this */
	if( !h2.lgH2ON /*|| !h2.nCallH2_this_zone*/ )
		return 0.;

	DEBUG_ENTRY( "H2_InterEnergy()" );

	/* There is no stack of levels, so we access all levels 
	 * uniquely by via the line stack with H2Lines[iElecHi][iVibHi][iRotHi][0][0][0].Hi */

	energy = 0.;
	/* loop over all possible levels */
	for( long iElecHi=0; iElecHi<mole.n_h2_elec_states; ++iElecHi )
	{
		for( long iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
		{
			for( long iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
			{
				fixit(); // how is ground defined here?  Must skip that level. 
				energy += H2Lines[iElecHi][iVibHi][iRotHi][0][0][0].Hi->Pop *
					H2Lines[iElecHi][iVibHi][iRotHi][0][0][0].EnergyErg;
			}
		}
	}
	return energy;
}
#endif

/*H2_RT_diffuse do emission from H2 - called from RT_diffuse */
void H2_RT_diffuse(void)
{
	long int iElecHi , iElecLo , iVibHi , iVibLo , iRotHi , iRotLo;

	if( !h2.lgH2ON || !h2.nCallH2_this_zone )
		return;

	DEBUG_ENTRY( "H2_RT_diffuse()" );

	/* loop over all possible lines */
	/* NB - this loop does not include the electronic lines */
	for( iElecHi=0; iElecHi<1; ++iElecHi )
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
					/* want to include all vibration states in lower level if different electronic level,
					* but only lower vibration levels if same electronic level */
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
							/* >>chng 03 feb 14, from !=0 to > 0 */
							/* >>chng 05 feb 07, use lgH2_line_exists */
							if( lgH2_line_exists[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo] )
							{
								H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].outline_resonance();
							}
						}
					}
				}
			}
		}
	}
	return;
}

/* RT for H2 lines */
void H2_RTMake( void )
{
	long int iElecHi , iElecLo , iVibHi , iVibLo , iRotHi , iRotLo;

	if( !h2.lgH2ON )
		return;

	DEBUG_ENTRY( "H2_RTMake()" );

	/* this routine drives calls to make RT relations for H2 molecule */
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
					/* want to include all vibration states in lower level if different electronic level,
					* but only lower vibration levels if same electronic level */
					long int nv = h2.nVib_hi[iElecLo];
					if( iElecLo==iElecHi )
						nv = iVibHi;
					for( iVibLo=0; iVibLo<=nv; ++iVibLo )
					{
						long nr = h2.nRot_hi[iElecLo][iVibLo];
						if( iElecLo==iElecHi && iVibHi==iVibLo )
							nr = iRotHi-1;

						mb6ci lgH2le = lgH2_line_exists.ptr(iElecHi,iVibHi,iRotHi,iElecLo,iVibLo,h2.Jlowest[iElecLo]);
						for( iRotLo=h2.Jlowest[iElecLo]; iRotLo<=nr; ++iRotLo )
						{
							/* >>chng 03 feb 14, change test from !=0 to >0 */
							/* >>chng 05 feb 07, use lgH2_line_exists */
							if( *lgH2le++ )
							{
								/* >>chng 03 jun 18, added 4th parameter in call to this routine - says to not
								 * include self-shielding of line across this zone.  This introduces a dr dependent
								 * variation in the line pumping rate, which made H2 abundance fluctuate due to
								 * Solomon process having slight dr-caused mole. */
								RT_line_one( &H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo], false, 0.f, 
									GetDopplerWidth(2.f*dense.AtomicWeight[ipHYDROGEN]) );
							}
						}
					}
				}
			}
		}
	}

	/* debug print population weighted transition probability */
	{
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			double sumpop = 0.;
			double sumpopA = 0.;
			for( iElecHi=1; iElecHi<mole.n_h2_elec_states; ++iElecHi )
			{
				for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
				{
					for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
					{
						/* now the lower levels */
						/* NB - X is the only lower level considered here, since we are only 
						 * concerned with excited electronic levels as a photodissociation process
						 * code exists to relax this assumption - simply change following to iElecHi */
						iElecLo = 0;
						/* want to include all vibration states in lower level if different electronic level,
						 * but only lower vibration levels if same electronic level */
						for( iVibLo=0; iVibLo<=h2.nVib_hi[iElecLo]; ++iVibLo )
						{
							long nr = h2.nRot_hi[iElecLo][iVibLo];
							for( iRotLo=h2.Jlowest[iElecLo]; iRotLo<=nr; ++iRotLo )
							{
								/* >>chng 03 feb 14, change test from !=0 to >0 */
								/* >>chng 05 feb 07, use lgH2_line_exists */
								/*if( H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Aul > 0. )*/
								if( lgH2_line_exists[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo] )
								{
									ASSERT( H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->Aul > 0. );

									sumpop += H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Lo->Pop;
									sumpopA += H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Lo->Pop*
										H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->Aul;
								}
							}
						}
					}
				}
			}
			fprintf(ioQQQ,"DEBUG sumpop = %.3e sumpopA= %.3e A=%.3e\n", sumpop, sumpopA, 
				sumpopA/SDIV(sumpop) );
		}
	}
	return;
}

/* increment optical depth for the H2 molecule, called from RT_tau_inc which is called  by cloudy,
 * one time per zone */
void H2_RT_tau_inc(void)
{
	long int iElecHi , iElecLo , iVibHi , iVibLo , iRotHi , iRotLo;

	/* >>chng 05 jan 26, now use LTE populations for small H2 abundance case, since electronic
	 * lines become self-shielding surprisingly quickly */
	if( !h2.lgH2ON /*|| !h2.nCallH2_this_zone*/ )
		return;

	DEBUG_ENTRY( "H2_RT_tau_inc()" );

	/* remember largest and smallest chemistry renormalization factor -
	 * if both networks are parallel will be unity,
	 * but only do this after we have stable solution */
	if( nzone > 0 && nCallH2_this_iteration>2 )
	{
		h2.renorm_max = MAX2( H2_renorm_chemistry , h2.renorm_max );
		h2.renorm_min = MIN2( H2_renorm_chemistry , h2.renorm_min );
	}

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
					/* want to include all vibration states in lower level if different electronic level,
					 * but only lower vibration levels if same electronic level */
					long int nv = h2.nVib_hi[iElecLo];
					if( iElecLo==iElecHi )
						nv = iVibHi;
					for( iVibLo=0; iVibLo<=nv; ++iVibLo )
					{
						long nr = h2.nRot_hi[iElecLo][iVibLo];
						if( iElecLo==iElecHi && iVibHi==iVibLo )
							nr = iRotHi-1;

						mb6ci lgH2le = lgH2_line_exists.ptr(iElecHi,iVibHi,iRotHi,iElecLo,iVibLo,h2.Jlowest[iElecLo]);
						for( iRotLo=h2.Jlowest[iElecLo]; iRotLo<=nr; ++iRotLo )
						{
							if( *lgH2le++ )
							{
								ASSERT( H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].ipCont > 0 );
								RT_line_one_tauinc( &H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo],
									-9, iRotHi, iVibLo, iRotLo, GetDopplerWidth(2.f*dense.AtomicWeight[ipHYDROGEN]) );
							}
						}/* iRotLo loop */
					}/* iVibLo loop */
				}/* iElecLo loop */
			}/* iRotHi loop */
		}/* iVibHi loop */
	}/* iElecHi loop */
	/*fprintf(ioQQQ,"\t%.3e\n",H2Lines[1][0][0][0][0][1].TauCon);*/
	return;
}


/* initialize optical depths in H2, called from RT_tau_init */
void H2_LineZero( void )
{
	long int iElecHi , iElecLo , iVibHi , iVibLo , iRotHi , iRotLo;

	if( !h2.lgH2ON )
		return;

	DEBUG_ENTRY( "H2_LineZero()" );

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
					/* want to include all vibration states in lower level if different electronic level,
					* but only lower vibration levels if same electronic level */
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
							/*if( H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Aul != 0. )*/
							/* >>chng 03 feb 14, from !=0 to > 0 */
							/*if( H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Aul > 0. )*/
							/* >>chng 05 feb 07, use lgH2_line_exists */
							if( lgH2_line_exists[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo] )
							{
								/* >>chng 03 feb 14, use TransitionZero rather than explicit sets */
								H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Zero();
							}
						}
					}
				}
			}
		}
	}
	return;
}

/* the large H2 molecule, called from RT_tau_reset */
void H2_RT_tau_reset( void )
{
	long int iElecHi , iElecLo , iVibHi , iVibLo , iRotHi , iRotLo;

	if( !h2.lgH2ON )
		return;

	DEBUG_ENTRY( "H2_RT_tau_reset()" );

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
					/* want to include all vibration states in lower level if different electronic level,
					* but only lower vibration levels if same electronic level */
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
							/* >>chng 03 feb 14, change test from !=0 to >0 */
							/*if( H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Aul > 0. )*/
							/* >>chng 05 feb 07, use lgH2_line_exists */
							if( lgH2_line_exists[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo] )
							{
								/* inward optical depth */
								RT_line_one_tau_reset( &H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo] );
							}
						}
					}
				}
			}
		}
	}
	return;
}

/* this is fraction of population that is within levels done with matrix */
static double frac_matrix;

/*H2_Level_low_matrix evaluate lower populations within X */
STATIC void H2_Level_low_matrix(
	/* total abundance within matrix */
	realnum abundance )
{

	/* will need to MALLOC space for these but only on first call */
	static double **AulEscp ,
		**col_str ,
		**AulDest, 
		/* AulPump[low][high] is rate (s^-1) from lower to upper level */
		**AulPump,
		**CollRate_levn,
		*pops,
		*create,
		*destroy,
		*depart,
		/* statistical weight */
		*stat_levn ,
		/* excitation energies in kelvin */
		*excit;
	bool lgDoAs;
	static long int levelAsEval=-1;
	long int ip;
	static bool lgFirst=true;
	long int i,
		j,
		ilo , 
		ihi,
		iElec,
		iElecHi,
		iVib,
		iRot,
		iVibHi,
		iRotHi;
	int nNegPop;
	bool lgDeBug,
		lgZeroPop;
	double rot_cooling , dCoolDT;
	static long int ndimMalloced = 0;
	double rateout , ratein;

	DEBUG_ENTRY( "H2_Level_low_matrix()" );

	/* option to not use the matrix */
	if( nXLevelsMatrix <= 1 )
	{
		return;
	}

	if( lgFirst )
	{
		/* check that not more levels than there are in X */
		if( nXLevelsMatrix > nLevels_per_elec[0] )
		{
			/* number is greater than number of levels within X */
			fprintf( ioQQQ, 
				" The total number of levels used in the matrix solver must be <= %li, the number of levels within X.\n Sorry.\n",
				nLevels_per_elec[0]);
			cdEXIT(EXIT_FAILURE);
		}
		/* will never do this again */
		lgFirst = false;
		/* remember how much space we malloced in case ever called with more needed */
		/* >>chng 05 jan 19, allocate max number of levels
		ndimMalloced = nXLevelsMatrix;*/
		ndimMalloced = nLevels_per_elec[0];
		/* allocate the 1D arrays*/
		excit = (double *)MALLOC( sizeof(double)*(size_t)(ndimMalloced) );
		stat_levn = (double *)MALLOC( sizeof(double)*(size_t)(ndimMalloced) );
		pops = (double *)MALLOC( sizeof(double)*(size_t)(ndimMalloced) );
		create = (double *)MALLOC( sizeof(double)*(size_t)(ndimMalloced) );
		destroy = (double *)MALLOC( sizeof(double)*(size_t)(ndimMalloced) );
		depart = (double *)MALLOC( sizeof(double)*(size_t)(ndimMalloced) );
		/* create space for the 2D arrays */
		AulPump = ((double **)MALLOC((size_t)(ndimMalloced)*sizeof(double *)));
		CollRate_levn = ((double **)MALLOC((size_t)(ndimMalloced)*sizeof(double *)));
		AulDest = ((double **)MALLOC((size_t)(ndimMalloced)*sizeof(double *)));
		AulEscp = ((double **)MALLOC((size_t)(ndimMalloced)*sizeof(double *)));
		col_str = ((double **)MALLOC((size_t)(ndimMalloced)*sizeof(double *)));
		for( i=0; i<(ndimMalloced); ++i )
		{
			AulPump[i] = ((double *)MALLOC((size_t)(ndimMalloced)*sizeof(double )));
			CollRate_levn[i] = ((double *)MALLOC((size_t)(ndimMalloced)*sizeof(double )));
			AulDest[i] = ((double *)MALLOC((size_t)(ndimMalloced)*sizeof(double )));
			AulEscp[i] = ((double *)MALLOC((size_t)(ndimMalloced)*sizeof(double )));
			col_str[i] = ((double *)MALLOC((size_t)(ndimMalloced)*sizeof(double )));
		}

		for( j=0; j < ndimMalloced; j++ )
		{
			stat_levn[j]=0;
			excit[j] =0;
		}
		/* the statistical weights of the levels
		 * and excitation potentials of each level relative to ground */
		for( j=0; j < ndimMalloced; j++ )
		{
			/* obtain the proper indices for the upper level */
			ip = H2_ipX_ener_sort[j];
			iVib = ipVib_H2_energy_sort[ip];
			iRot = ipRot_H2_energy_sort[ip];

			/* statistical weights for each level */
			stat_levn[j] = H2_stat[0][iVib][iRot];
			/* excitation energy of each level relative to ground, in K */
			excit[j] = energy_wn[0][iVib][iRot]*T1CM;
		}

		for( j=0; j < ndimMalloced-1; j++ )
		{
			/* make sure that the energies are ok */
			ASSERT( excit[j+1] > excit[j] );
		}
	}
	/* end malloc space and creating constant terms */

	/* this is test for call with too many rotation levels to handle - 
	 * logic needs for largest model atom to be called first */
	if( nXLevelsMatrix > ndimMalloced )
	{
		fprintf(ioQQQ," H2_Level_low_matrix has been called with the number of rotor levels greater than space allocated.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* all elements are used, and must be set to zero */
	for( i=0; i < nXLevelsMatrix; i++ )
	{
		pops[i] = 0.;
		depart[i] = 0;
		for( j=0; j < nXLevelsMatrix; j++ )
		{
			col_str[j][i] = 0.;
		}
	}

	/* do we need to reevaluate radiative quantities?  only do this one time per zone */
	if( nzone!=nzoneAsEval || iteration!=iterationAsEval || nXLevelsMatrix!=levelAsEval)
	{
		lgDoAs = true;
		nzoneAsEval = nzone;
		iterationAsEval = iteration;
		levelAsEval = nXLevelsMatrix;
		ASSERT( levelAsEval <= ndimMalloced );
	}
	else
	{
		lgDoAs = false;
	}

	/* all elements are used, and must be set to zero */
	if( lgDoAs )
	{
		for( i=0; i < nXLevelsMatrix; i++ )
		{
			pops[i] = 0.;
			depart[i] = 0;
			for( j=0; j < nXLevelsMatrix; j++ )
			{
				AulEscp[j][i] = 0.;
				AulDest[j][i] = 0.;
				AulPump[j][i] = 0.;
				CollRate_levn[j][i] = 0.;
			}
		}
	}

	/* find all radiative interactions within matrix, and between
	 * matrix and upper X and excited electronic states */
	iElec = 0;
	for( ilo=0; ilo < nXLevelsMatrix; ilo++ )
	{
		ip = H2_ipX_ener_sort[ilo];
		iRot = ipRot_H2_energy_sort[ip];
		iVib = ipVib_H2_energy_sort[ip];

		/* H2_X_sink[ilo] includes all processes that destroy H2 in one step, 
		 * these include cosmic ray ionization and dissociation, photodissociation,
		 * BUT NOT THE SOLOMON process, which, directly, only goes to excited
		 * electronic states */
		destroy[ilo] = H2_X_sink[ilo];

		/* rates H2 is created from grains and H- units cm-3 s-1, evaluated in mole_H2_form */
		create[ilo] = H2_X_source[ilo];

		/* this loop does radiative decays from upper states inside matrix, 
		 * and upward pumps within matrix region into this lower level */
		if( lgDoAs )
		{
			for( ihi=ilo+1; ihi<nXLevelsMatrix; ++ihi )
			{
				ip = H2_ipX_ener_sort[ihi];
				iRotHi = ipRot_H2_energy_sort[ip];
				iVibHi = ipVib_H2_energy_sort[ip];
				ASSERT( H2_energies[ip] <= H2_energies[H2_ipX_ener_sort[nXLevelsMatrix-1]] );
				/* general case - but line may not actually exist */
				if( (abs(iRotHi-iRot)==2 || (iRotHi-iRot)==0 ) && (iVib<=iVibHi) )
				{
					/* >>chng 05 feb 07, use lgH2_line_exists */
					if( lgH2_line_exists[0][iVibHi][iRotHi][0][iVib][iRot] )
					{
						ASSERT( H2Lines[0][iVibHi][iRotHi][0][iVib][iRot].ipCont > 0 );

						/* NB - the destruction probability is included in 
						 * the total and the destruction is set to zero
						 * since we want to only count one ots rate, in 
						 * main calling routine, and do not want matrix 
						 * solver below to include it */
						AulEscp[ihi][ilo] = H2Lines[0][iVibHi][iRotHi][0][iVib][iRot].Emis->Aul*(
							H2Lines[0][iVibHi][iRotHi][0][iVib][iRot].Emis->Pesc + 
							H2Lines[0][iVibHi][iRotHi][0][iVib][iRot].Emis->Pdest +
							H2Lines[0][iVibHi][iRotHi][0][iVib][iRot].Emis->Pelec_esc);
						AulDest[ilo][ihi] = 0.;
						AulPump[ilo][ihi] = H2Lines[0][iVibHi][iRotHi][0][iVib][iRot].Emis->pump;
					}
				}
			}
		}

		iElecHi = 0;
		iElec = 0;
		rateout = 0.;
		ratein = 0.;
		/* now do all levels within X, which are above nXLevelsMatrix,
		 * the highest level inside the matrix */
		for( ihi=nXLevelsMatrix; ihi<nLevels_per_elec[0]; ++ihi )
		{
			ip = H2_ipX_ener_sort[ihi];
			iRotHi = ipRot_H2_energy_sort[ip];
			iVibHi = ipVib_H2_energy_sort[ip];
			if( (abs(iRotHi-iRot)==2 || (iRotHi-iRot)==0 ) && (iVib<=iVibHi) )
			{
				if( lgH2_line_exists[iElecHi][iVibHi][iRotHi][iElec][iVib][iRot] )
				{
					ASSERT( H2Lines[iElecHi][iVibHi][iRotHi][iElec][iVib][iRot].ipCont > 0 );

					/* these will enter as net creation terms in creation vector, with
					 * units cm-3 s-1
					 * radiative transitions from above the matrix within X */
					ratein +=
						H2_populations[iElecHi][iVibHi][iRotHi] *
						(H2Lines[iElecHi][iVibHi][iRotHi][iElec][iVib][iRot].Emis->Aul*
						(H2Lines[iElecHi][iVibHi][iRotHi][iElec][iVib][iRot].Emis->Pesc + 
						H2Lines[iElecHi][iVibHi][iRotHi][iElec][iVib][iRot].Emis->Pelec_esc + 
						H2Lines[iElecHi][iVibHi][iRotHi][iElec][iVib][iRot].Emis->Pdest)+H2Lines[iElecHi][iVibHi][iRotHi][iElec][iVib][iRot].Emis->pump *
							H2Lines[iElecHi][iVibHi][iRotHi][iElec][iVib][iRot].Lo->g/
							H2Lines[iElecHi][iVibHi][iRotHi][iElec][iVib][iRot].Hi->g);
					/* rate out has units s-1 - destroys current lower level */
					rateout +=
						H2Lines[iElecHi][iVibHi][iRotHi][iElec][iVib][iRot].Emis->pump;
				}
			}
		}

		/* all states above the matrix but within X */
		create[ilo] += ratein;

		/* rates out of matrix into levels in X but above matrix */
		destroy[ilo] += rateout;

		/* Solomon process, this sum dos all pump and decays from all electronic excited states */
		/* radiative rates [cm-3 s-1] from electronic excited states into X only vibration and rot */
		create[ilo] += H2_X_rate_from_elec_excited[iVib][iRot];

		/* radiative & cosmic ray rates [s-1] to electronic excited states from X only vibration and rot */
		destroy[ilo] += H2_X_rate_to_elec_excited[iVib][iRot];
	}

	/* this flag set with atom H2 trace matrix */
	if( mole.nH2_TRACE >= mole.nH2_trace_matrix )
		lgDeBug = true;
	else
		lgDeBug = false;

	/* now evaluate the rates for all transitions within matrix */
	for( ilo=0; ilo < nXLevelsMatrix; ilo++ )
	{
		ip = H2_ipX_ener_sort[ilo];
		iRot = ipRot_H2_energy_sort[ip];
		iVib = ipVib_H2_energy_sort[ip];
		if( lgDoAs )
		{
			if(lgDeBug)fprintf(ioQQQ,"DEBUG H2_Level_low_matrix, ilo=%li",ilo);
			for( ihi=ilo+1; ihi < nXLevelsMatrix; ihi++ )
			{
				ip = H2_ipX_ener_sort[ihi];
				iRotHi = ipRot_H2_energy_sort[ip];
				iVibHi = ipVib_H2_energy_sort[ip];
				/* >>chng 05 may 31, replace with simple expresion */
				CollRate_levn[ihi][ilo] = H2_X_coll_rate[ihi][ilo];

				/*create[ilo] +=CollRate_levn[ihi][ilo]*H2_populations[0][iVibHi][iRotHi];*/
				if(lgDeBug)fprintf(ioQQQ,"\t%.1e",CollRate_levn[ihi][ilo]);

				/* now get upward excitation rate - units s-1 */
				CollRate_levn[ilo][ihi] = CollRate_levn[ihi][ilo]*
					H2_Boltzmann[0][iVibHi][iRotHi]/SDIV(H2_Boltzmann[0][iVib][iRot])*
					H2_stat[0][iVibHi][iRotHi] / 
					H2_stat[0][iVib][iRot];
			}
		}

		if(lgDeBug)fprintf(ioQQQ,"\n");

		/* now do all collisions for levels within X, which are above nXLevelsMatrix,
		 * the highest level inside the matrix */
		iElecHi = 0;

		for( ihi=nXLevelsMatrix; ihi<nLevels_per_elec[0]; ++ihi )
		{
			ip = H2_ipX_ener_sort[ihi];
			iRotHi = ipRot_H2_energy_sort[ip];
			iVibHi = ipVib_H2_energy_sort[ip];
			rateout = 0;
			/* first do downward deexcitation rate */
			/* >>chng 04 sep 14, do all levels */
			/* >>chng 05 may 31, use summed rate */
			ratein = H2_X_coll_rate[ihi][ilo];
			if(lgDeBug)fprintf(ioQQQ,"\t%.1e",ratein);

			/* now get upward excitation rate */
			rateout = ratein *
				H2_Boltzmann[0][iVibHi][iRotHi]/SDIV(H2_Boltzmann[0][iVib][iRot])*
				H2_stat[0][iVibHi][iRotHi]/H2_stat[0][iVib][iRot];

			/* these are general entries and exits going into vector */
			create[ilo] += ratein*H2_populations[iElecHi][iVibHi][iRotHi];
			destroy[ilo] += rateout;
		}
	}

	/* H2 grain interactions
	 * >>chng 05 apr 30,GS, Instead of hmi.H2_total, the specific populations are used because high levels have much less
	 * populations than ground levels which consists most of the H2 population.*/
	if( lgDoAs )
	{
		for( ihi=2; ihi < nXLevelsMatrix; ihi++ )
		{

			ip = H2_ipX_ener_sort[ihi];
			iVibHi = ipVib_H2_energy_sort[ip];
			iRotHi = ipRot_H2_energy_sort[ip];

			/* collisions with grains goes to either J=1 or J=0 depending on 
			* spin of upper level - this conserves op ratio - following
			* var is 1 if ortho, 0 if para, so this conserves op ratio
			* units are s-1 */
			CollRate_levn[ihi][H2_lgOrtho[0][iVibHi][iRotHi]] += hmi.rate_grain_h2_op_conserve;
		}

		/* H2 ortho - para conversion on grain surface,
		 * rate (s-1) all v,J levels go to 0 or 1 */
		CollRate_levn[1][0] += 
			(realnum)(hmi.rate_grain_h2_J1_to_J0);
	}

	/* now all levels in X above the matrix */
	for( ihi=nXLevelsMatrix; ihi<nLevels_per_elec[0]; ++ihi )
	{
		ip = H2_ipX_ener_sort[ihi];
		iVibHi = ipVib_H2_energy_sort[ip];
		iRotHi = ipRot_H2_energy_sort[ip];

		/* these collisions all go into 0 or 1 depending on whether upper level was ortho or para 
		 * units are cm-3 s-1 - rate new molecules appear in matrix */
		create[H2_lgOrtho[0][iVibHi][iRotHi]] += H2_populations[0][iVibHi][iRotHi]*hmi.rate_grain_h2_op_conserve;
	}

	/* debug print individual contributors to matrix elements */
	{
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC || lgDeBug)
		{
			fprintf(ioQQQ,"DEBUG H2 matexcit");
			for(ilo=0; ilo<nXLevelsMatrix; ++ilo )
			{
				fprintf(ioQQQ,"\t%li",ilo );
			}
			fprintf(ioQQQ,"\n");
			for(ihi=0; ihi<nXLevelsMatrix;++ihi)
			{
				fprintf(ioQQQ,"\t%.2e",excit[ihi] );
			}
			fprintf(ioQQQ,"\n");
			for(ihi=0; ihi<nXLevelsMatrix;++ihi)
			{
				fprintf(ioQQQ,"\t%.2e",stat_levn[ihi] );
			}
			fprintf(ioQQQ,"\n");

			fprintf(ioQQQ,"AulEscp[n][]\\[][n] = Aul*Pesc\n");
			for(ilo=0; ilo<nXLevelsMatrix; ++ilo )
			{
				fprintf(ioQQQ,"\t%li",ilo );
			}
			fprintf(ioQQQ,"\n");
			for(ihi=0; ihi<nXLevelsMatrix;++ihi)
			{
				fprintf(ioQQQ,"%li", ihi);
				for(ilo=0; ilo<nXLevelsMatrix; ++ilo )
				{
					fprintf(ioQQQ,"\t%.2e",AulEscp[ilo][ihi] );
				}
				fprintf(ioQQQ,"\n");
			}

			fprintf(ioQQQ,"AulPump [n][]\\[][n]\n");
			for(ilo=0; ilo<nXLevelsMatrix; ++ilo )
			{
				fprintf(ioQQQ,"\t%li",ilo );
			}
			fprintf(ioQQQ,"\n");
			for(ihi=0; ihi<nXLevelsMatrix;++ihi)
			{
				fprintf(ioQQQ,"%li", ihi);
				for(ilo=0; ilo<nXLevelsMatrix; ++ilo )
				{
					fprintf(ioQQQ,"\t%.2e",AulPump[ihi][ilo] );
				}
				fprintf(ioQQQ,"\n");
			}

			fprintf(ioQQQ,"CollRate_levn [n][]\\[][n]\n");
			for(ilo=0; ilo<nXLevelsMatrix; ++ilo )
			{
				fprintf(ioQQQ,"\t%li",ilo );
			}
			fprintf(ioQQQ,"\n");
			for(ihi=0; ihi<nXLevelsMatrix;++ihi)
			{
				fprintf(ioQQQ,"%li", ihi);
				for(ilo=0; ilo<nXLevelsMatrix; ++ilo )
				{
					fprintf(ioQQQ,"\t%.2e",CollRate_levn[ihi][ilo] );
				}
				fprintf(ioQQQ,"\n");
			}
			fprintf(ioQQQ,"SOURCE");
			for(ihi=0; ihi<nXLevelsMatrix;++ihi)
			{
				fprintf(ioQQQ,"\t%.2e",create[ihi]);
			}
			fprintf(ioQQQ,"\nSINK");
			for(ihi=0; ihi<nXLevelsMatrix;++ihi)
			{
				fprintf(ioQQQ,"\t%.2e",destroy[ihi]);
			}
			fprintf(ioQQQ,"\n");
		}
	}

	atom_levelN(
		/* number of levels */
		nXLevelsMatrix,
		abundance,
		stat_levn,
		excit,
		'K',
		pops,
		depart,
		/* net transition rate, A * escape prob, s-1, indices are [upper][lower] */
		&AulEscp, 
		/* col str from high to low */
		&col_str, 
		&AulDest,
		&AulPump,
		&CollRate_levn,
		create,
		destroy,
		/* say that we have evaluated the collision rates already */
		true,
		&rot_cooling,
		&dCoolDT,
		" H2 ",
		/* nNegPop positive if negative pops occurred, negative if too cold */
		&nNegPop,
		&lgZeroPop,
		lgDeBug );/* option to print stuff - set to true for debug printout */

	for( i=0; i< nXLevelsMatrix; ++i )
	{
		ip = H2_ipX_ener_sort[i];
		iRot = ipRot_H2_energy_sort[ip];
		iVib = ipVib_H2_energy_sort[ip];
		/* >>chng 05 feb 08, do not update h2_old_populations here, since not done
		 * like this anywhere else - only update H2_populations here, and let single loop
		 * in main calling routine handle updating various forms of the population
		H2_old_populations[0][iVib][iRot] = H2_populations[0][iVib][iRot]; */
		H2_populations[0][iVib][iRot] = pops[i];
	}

	if( 0 && mole.nH2_TRACE >= mole.nH2_trace_full) 
	{
		/*static int nn=0; ++nn; if( nn>5)cdEXIT(1);*/
		/* print pops that came out of matrix */
		fprintf(ioQQQ,"\n DEBUG H2_Level_lowJ hmi.H2_total: %.3e matrix rel pops\n",hmi.H2_total);
		fprintf(ioQQQ,"v\tJ\tpop\n");
		for( i=0; i<nXLevelsMatrix; ++i )
		{
			ip = H2_ipX_ener_sort[i];
			iRot = ipRot_H2_energy_sort[ip];
			iVib = ipVib_H2_energy_sort[ip];
			fprintf(ioQQQ,"%3li\t%3li\t%.3e\t%.3e\t%.3e\n",
				iVib , iRot , H2_populations[0][iVib][iRot]/hmi.H2_total , create[i] , destroy[i]);
		}
	}

	/* nNegPop positive if negative pops occurred, negative if too cold */
	if( nNegPop > 0 )
	{
		fprintf(ioQQQ," H2_Level_low_matrix called atom_levelN which returned negative H2_populations.\n");
		ConvFail( "pops" , "H2" );
	}
	return;
}
/* do level H2_populations for H2, called by Hydrogenic after ionization and H chemistry
 * has been recomputed */
void H2_LevelPops( void )
{
	static double TeUsedColl=-1.f;
	double H2_renorm_conserve=0.,
		H2_renorm_conserve_init=0. ,
		sumold, 
		H2_BigH2_H2s,
		H2_BigH2_H2g;
	double old_solomon_rate=-1.;
	long int iElecHi , iElecLo , iVibHi , iVibLo , iRotHi , iRotLo;
	long int i;
	long int n_pop_oscil = 0;
	int kase=0;
	bool lgConv_h2_soln,
		lgPopsConv_total,
		lgPopsConv_relative,
		lgHeatConv,
		lgSolomonConv,
		lgOrthoParaRatioConv;
	double quant_old=-1.,
		quant_new=-1.;

	bool lgH2_pops_oscil=false,
		lgH2_pops_ever_oscil=false;
	long int nEner,
		ipHi, ipLo;
	long int iElec , iVib , iRot,ip;
	double sum_pops_matrix;
	realnum collup;
	/* old and older ortho - para ratios, used to determine whether soln is converged */
	static double ortho_para_old=0. , ortho_para_older=0. , ortho_para_current=0.;
	realnum frac_new_oscil=1.f;

	/* keep track of changes in population */
	double PopChgMax_relative=0. , PopChgMaxOld_relative=0., PopChgMax_total=0., PopChgMaxOld_total=0.;
	long int iRotMaxChng_relative , iVibMaxChng_relative,
		iRotMaxChng_total , iVibMaxChng_total,
		nXLevelsMatrix_save;
	double popold_relative , popnew_relative , popold_total , popnew_total;
	/* reason not converged */
	char chReason[100];

	double flux_accum_photodissoc_BigH2_H2g, flux_accum_photodissoc_BigH2_H2s;
	long int ip_H2_level;

	/* these are convergence criteria - will be increased during search phase */
	double converge_pops_relative=0.1 ,
		converge_pops_total=1e-3, 
		converge_ortho_para=1e-2;

	DEBUG_ENTRY( "H2_LevelPops()" );

	/* H2 not on, so space not allocated and return,
	 * also return if calculation has been declared a failure */
	if( !h2.lgH2ON || lgAbort )
		return;

	if(mole.nH2_TRACE >= mole.nH2_trace_full ) 
	{
		fprintf(ioQQQ,
			"\n***************H2_LevelPops call %li this iteration, zone is %.2f, H2/H:%.e Te:%e ne:%e\n", 
			nCallH2_this_iteration,
			fnzone,
			hmi.H2_total/dense.gas_phase[ipHYDROGEN],
			phycon.te,
			dense.eden
			);
	}
	else if( mole.nH2_TRACE >= mole.nH2_trace_final )
	{
		static long int nzone_prt=-1;
		if( nzone!=nzone_prt )
		{
			nzone_prt = nzone;
			fprintf(ioQQQ,"DEBUG zone %li H2/H:%.3e Te:%.3e *ne:%.3e n(H2):%.3e\n",
				nzone,
				hmi.H2_total / dense.gas_phase[ipHYDROGEN],
				phycon.te,
				dense.eden,
				hmi.H2_total );
		}
	}

	/* evaluate Boltzmann factors and LTE unit population - for trivial abundances
	 * LTE populations are used in place of full solution */
	mole_H2_LTE();

	/* zero out H2_populations and cooling, and return, if H2 fraction is small
	 * but, if H2 has ever been done, redo irregardless of abundance -
	 * if large H2 is ever evaluated then mole.H2_to_H_limit is ignored */
	if( (!hmi.lgBigH2_evaluated && hmi.H2_total/dense.gas_phase[ipHYDROGEN] < mole.H2_to_H_limit )
		|| hmi.H2_total < 1e-20 )
	{
		/* will not compute the molecule */
		if( mole.nH2_TRACE >= mole.nH2_trace_full ) 
			fprintf(ioQQQ,
			"  H2_LevelPops pops too small, not computing, set to LTE and return, H2/H is %.2e and mole.H2_to_H_limit is %.2e.",
			hmi.H2_total/dense.gas_phase[ipHYDROGEN] ,
			mole.H2_to_H_limit);
		H2_zero_pops_too_low();
		/* end of zero abundance branch */
		return;
	}

	/* check whether we need to update the H2_Boltzmann factors, LTE level H2_populations,
	 * and partition function.  LTE level pops normalized by partition function,
	 * so sum of pops is unity */

	/* say that H2 has been computed, ignore previous limit to abund
	 * in future - this is to prevent oscillations as model is engaged */
	hmi.lgBigH2_evaluated = true;
	/* end loop setting H2_Boltzmann factors, partition function, and LTE H2_populations */

	/* >>chng 05 jun 21,
	 * during search phase we want to use full matrix - save number of levels so that
	 * we can restore it */
	nXLevelsMatrix_save = nXLevelsMatrix;
	if( conv.lgSearch )
	{
		nXLevelsMatrix = nLevels_per_elec[0];
	}

	/* 05 oct 27, had only reevaluated collision rates when 5% change in temperature
	 * caused temp failures in large G0 sims - 
	 * do not check whether we need to update the collision rates but
	 * reevaluate them always  
	 * >>chng 05 nov 04, above caused a 25% increase in the exec time for constant-T sims
	 * in test suite- original code had reevaluated if > 0.05 change in T - was too much
	 * change to 10x smaller, change > 0.005 */
	if( !fp_equal(phycon.te,TeUsedColl) )
	{
		H2_CollidRateEvalAll();
		TeUsedColl = phycon.te;
	}

	/* set the H2_populations when this is the first call to this routine on 
	 * current iteration- will use LTE H2_populations - populations were set by
	 * call to 	mole_H2_LTE before above block */
	if( nCallH2_this_iteration==0 || mole.lgH2_LTE )
	{
		/* very first call so use LTE H2_populations */
		if(mole.nH2_TRACE >= mole.nH2_trace_full ) 
			fprintf(ioQQQ,"H2 1st call - using LTE level pops\n");

		for( iElec=0; iElec<mole.n_h2_elec_states; ++iElec )
		{
			pops_per_elec[iElec] = 0.;
			for( iVib=0; iVib<=h2.nVib_hi[iElec]; ++iVib )
			{
				pops_per_vib[iElec][iVib] = 0.;
				for( iRot=h2.Jlowest[iElec]; iRot<=h2.nRot_hi[iElec][iVib]; ++iRot )
				{
					/* LTE populations are for unit H2 density, so need to multiply
					 * by total H2 density */
					H2_old_populations[iElec][iVib][iRot] = 
						(realnum)H2_populations_LTE[iElec][iVib][iRot]*hmi.H2_total;
					H2_populations[iElec][iVib][iRot] = H2_old_populations[iElec][iVib][iRot];
				}
			}
		}
		/* first guess at ortho and para densities */
		h2.ortho_density = 0.75*hmi.H2_total;
		h2.para_density = 0.25*hmi.H2_total;
		ortho_para_current = h2.ortho_density / SDIV( h2.para_density );
		ortho_para_older = ortho_para_current;
		ortho_para_old = ortho_para_current;
		/* this is the fraction of the H2 pops that are within the levels done with a matrix */
		frac_matrix = 1.;
	}

	/* find sum of all H2_populations in X */
	iElec = 0;
	pops_per_elec[0] = 0.;
	for( iVib=0; iVib<=h2.nVib_hi[iElec]; ++iVib )
	{
		pops_per_vib[0][iVib] = 0.;
		for( iRot=h2.Jlowest[iElec]; iRot<=h2.nRot_hi[iElec][iVib]; ++iRot )
		{
			pops_per_elec[0] += H2_populations[iElec][iVib][iRot];
			pops_per_vib[0][iVib] += H2_populations[iElec][iVib][iRot];
		}
	}
	ASSERT( pops_per_elec[0]>SMALLFLOAT );
	/* now renorm the old populations to the correct current H2 density, 
	 * At this point pops_per_elec[0] (pop of X)
	 * is the result of the previous evaluation of the H2 population, 
	 * following is the ratio of the current chemistry solution H2 to the previous H2 */
	H2_renorm_chemistry = hmi.H2_total/ SDIV(pops_per_elec[0]);

	/* >>chng 05 jul 13, TE, 
	 * evaluate ratio of H2g and H2s from chemical network and big molecule model */
	iElec = 0;
	hmi.H2_chem_BigH2_H2g = 0.;
	hmi.H2_chem_BigH2_H2s = 0.;
	for( iVib=0; iVib<=h2.nVib_hi[iElec]; ++iVib )
	{
		for( iRot=h2.Jlowest[iElec]; iRot<=h2.nRot_hi[iElec][iVib]; ++iRot )
		{ 
			if( energy_wn[0][iVib][iRot] > ENERGY_H2_STAR )
			{
				hmi.H2_chem_BigH2_H2s += H2_populations[iElec][iVib][iRot];

			}
			else
			{
				hmi.H2_chem_BigH2_H2g += H2_populations[iElec][iVib][iRot];

			}
		}
	}

	hmi.H2_chem_BigH2_H2g = hmi.Hmolec[ipMH2g]/SDIV(hmi.H2_chem_BigH2_H2g);
	hmi.H2_chem_BigH2_H2s = hmi.Hmolec[ipMH2s]/SDIV(hmi.H2_chem_BigH2_H2s);


	if(mole.nH2_TRACE >= mole.nH2_trace_full) 
		fprintf(ioQQQ,
			"H2 H2_renorm_chemistry is %.4e, hmi.H2_total is %.4e pops_per_elec[0] is %.4e\n",
			H2_renorm_chemistry ,
			hmi.H2_total,
			pops_per_elec[0]);

	/* renormalize all level populations for the current chemical solution */
	iElec = 0;
	for( iVib=0; iVib<=h2.nVib_hi[iElec]; ++iVib )
	{
		for( iRot=h2.Jlowest[iElec]; iRot<=h2.nRot_hi[iElec][iVib]; ++iRot )
		{
			H2_populations[iElec][iVib][iRot] *= H2_renorm_chemistry;
			H2_old_populations[iElec][iVib][iRot] = H2_populations[iElec][iVib][iRot];
		}
	}
	ASSERT( fabs(h2.ortho_density+h2.para_density-hmi.H2_total)/hmi.H2_total < 0.001 );

	if(mole.nH2_TRACE >= mole.nH2_trace_full )
		fprintf(ioQQQ,
		" H2 entry, old pops sumed to %.3e, renorm to htwo den of %.3e\n",
		pops_per_elec[0],
		hmi.H2_total);

	/* >>chng 05 feb 10, reset convergence criteria if we are in search phase */
	if( conv.lgSearch )
	{
		converge_pops_relative *= 2.; /*def is 0.1 */
		converge_pops_total *= 3.;    /*def is 1e-3*/
		converge_ortho_para *= 3.;    /*def is 1e-2*/
	}

	/* update state specific rates in H2_X_formation (cm-3 s-1) that H2 forms from grains and H- */
	mole_H2_form();

	/* evaluate total collision rates */
	H2_X_coll_rate_evaluate();

	/* this flag will say whether H2 H2_populations have converged,
	 * by comparing old and new values */
	lgConv_h2_soln = false;
	/* this will count number of passes around following loop */
	loop_h2_pops = 0;
	{
		static long int nzoneEval=-1;
		if( nzone != nzoneEval )
		{
			nzoneEval = nzone;
			/* this is number of zones with H2 solution in this iteration */
			++nH2_zone;
		}
	}

	/* begin - start level population solution
	 * first do electronic excited states, Lyman, Werner, etc
	 * using old solution for X
	 * then do matrix if used, then solve for pops of rest of X 
	 * >>chng 04 apr 06, subtract number of oscillations from limit - don't waste loops 
	 * if solution is unstable */
	while( loop_h2_pops < LIM_H2_POP_LOOP-n_pop_oscil && !lgConv_h2_soln && !mole.lgH2_LTE )
	{
		double rate_in , rate_out;
		static double old_HeatH2Dexc_BigH2=0., HeatChangeOld=0. , HeatChange=0.;

		/* this is number of trips around loop this time */
		++loop_h2_pops;
		/* this is number of times through this loop in entire iteration */
		++nH2_pops;

		/* beginning solution for electronic excited states
		 * loop over all possible pumping routes to excited electronic states
		 * to get radiative excitation and dissociation rates */
		hmi.H2_H2g_to_H2s_rate_BigH2 = 0.;

		rate_out = 0.;

		/* these will store radiative rates between electronic excited states and X */
		iElec = 0;
		for( iVib=0; iVib<=h2.nVib_hi[iElec]; ++iVib )
		{
			for( iRot=h2.Jlowest[iElec]; iRot<=h2.nRot_hi[iElec][iVib]; ++iRot )
			{
				/* radiative rates [cm-3 s-1] from electronic excited states into X vibration and rot */
				H2_X_rate_from_elec_excited[iVib][iRot] = 0.;
				/* radiative & cosmic ray rates [s-1] to electronic excited states from X */
				H2_X_rate_to_elec_excited[iVib][iRot] = 0.;
			}
		}

		iElecHi = -INT32_MAX;
		/* solve for population of electronic excited states by back-substitution */
		for( iElecHi=1; iElecHi<mole.n_h2_elec_states; ++iElecHi )
		{
			/* for safety, these loop vars are set to insane value, will crash
			 * if used without being set */
			iVib = -INT32_MAX;
			iRot = -INT32_MAX;
			iVibLo = -INT32_MAX;
			iRotLo = -INT32_MAX;
			iVibHi = -INT32_MAX;
			iRotHi = -INT32_MAX;

			/* this will be total population in each electronic state */
			pops_per_elec[iElecHi] = 0.;

			if(mole.nH2_TRACE >= mole.nH2_trace_full) 
				fprintf(ioQQQ," Pop(e=%li):",iElecHi);

			double CosmicRayHILyaExcitationRate = ( hmi.lgLeidenCRHack ) ? secondaries.x12tot : 0.;

			/* electronic excited state, loop over all vibration */
			for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
			{
				pops_per_vib[iElecHi][iVibHi] = 0.;
				/* ======================= EXCITED ELEC STATE POPS LOOP =====================*/
				for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
				{
					/* Solomon process done here,
					 * sum of all rates into and out of these upper levels 
					 * all inward rates have units cm-3 s-1 */
					rate_in = 0.;
					/* this term is spontaneous dissociation of excited electronic states into 
					 * the X continuum 
					 * all outward rates have units s-1 */
					rate_out = H2_dissprob[iElecHi][iVibHi][iRotHi];

					realnum H2gHi = H2Lines[iElecHi][iVibHi][iRotHi][0][0][0].Hi->g;

					/* now loop over all levels within X find Solomon rate */
					iElecLo=0;

					for( iVibLo=0; iVibLo<=h2.nVib_hi[iElecLo]; ++iVibLo )
					{
						mb6ci lgH2le = lgH2_line_exists.ptr(iElecHi,iVibHi,iRotHi,iElecLo,iVibLo,h2.Jlowest[iElecLo]);
						mt6ci H2L = H2Lines.ptr(iElecHi,iVibHi,iRotHi,iElecLo,iVibLo,h2.Jlowest[iElecLo]);
						long nr = h2.nRot_hi[iElecLo][iVibLo];
						for( iRotLo=h2.Jlowest[iElecLo]; iRotLo<=nr; ++iRotLo )
						{
							/* consider all radiatively permitted transitions between X and excit states */
							if( *lgH2le++ )
							{
								ASSERT( H2L->ipCont > 0 );
								/* solve electronic excited state, 
								 * rate lower level in X goes to electronic excited state, s-1 
								 * first term is direct pump, second is cosmic ray excitation */
								double rate_one = H2L->Emis->pump 
								/* collisional excitation of singlets by non-thermal electrons 
								 * this is stored ratio of electronic transition relative 
								 * cross section relative to the HI Lya cross section  */
								+CosmicRayHILyaExcitationRate*H2L->Coll.col_str;

								/* this is a permitted electronic transition, must preserve nuclear spin */
								ASSERT( H2_lgOrtho[iElecHi][iVibHi][iRotHi] == H2_lgOrtho[iElecLo][iVibLo][iRotLo] );

								/* this is the rate [cm-3 s-1] electrons move into the upper level from X */
								rate_in += H2_old_populations[iElecLo][iVibLo][iRotLo]*rate_one;

								/* this is total X -> excited electronic state rate, cm-3 s-1 */
								/*rate_tot += H2_old_populations[iElecLo][iVibLo][iRotLo]*rate_one;*/

								/* rate [s-1] from levels within X to electronic excited states,
								 * includes photoexcitation and cosmic ray excitation */
								H2_X_rate_to_elec_excited[iVibLo][iRotLo] += rate_one;

								/* excitation rate for Solomon process - this currently has units
								 * cm-3 s-1 but will be divided by total pop of X and become s-1 */
								/* this has unit s-1 and will be used for H2g->H2s */
								rate_one = H2L->Emis->Aul*
									/* escape and destruction */
									(H2L->Emis->Pesc + 
									H2L->Emis->Pelec_esc + 
									H2L->Emis->Pdest) +
									/* induced emission down */
									H2L->Emis->pump *
									H2L->Lo->g/H2gHi;

								/* this is the rate [s-1] electrons leave the excited electronic upper level
								 * and decay into X - will be used to get pops of electronic excited states */
								rate_out += rate_one;

								ASSERT( rate_in >= 0. && rate_out >= 0. );
							}
							++H2L;
						}
					}

					/* update population [cm-3] of the electronic excited state this only includes 
					 * radiative processes between X and excited electronic states, and cosmic rays - 
					 * thermal collisions are neglected
					 * X is done below and includes all processes */
					H2_rad_rate_out[iElecHi][iVibHi][iRotHi] = rate_out;
					double H2popHi = rate_in / SDIV( rate_out );
					H2_populations[iElecHi][iVibHi][iRotHi] = H2popHi;
					if( H2_old_populations[iElecHi][iVibHi][iRotHi]==0. )
						H2_old_populations[iElecHi][iVibHi][iRotHi] = H2popHi;

					/* for this population, find rate electronic states decay into X */
					/* now loop over all levels within X find Solomon rate */
					iElecLo=0;
					for( iVibLo=0; iVibLo<=h2.nVib_hi[iElecLo]; ++iVibLo )
					{
						mb6ci lgH2le = lgH2_line_exists.ptr(iElecHi,iVibHi,iRotHi,iElecLo,iVibLo,h2.Jlowest[iElecLo]);
						mt6ci H2L = H2Lines.ptr(iElecHi,iVibHi,iRotHi,iElecLo,iVibLo,h2.Jlowest[iElecLo]);

						long nr = h2.nRot_hi[iElecLo][iVibLo];
						for( iRotLo=h2.Jlowest[iElecLo]; iRotLo<=nr; ++iRotLo )
						{
							if( *lgH2le++ )
							{
								ASSERT( H2L->ipCont > 0 );

								double rate_one =
									/* >>chng 05 may 31, from old to new as in matrix  */
									/*H2_old_populations[iElecHi][iVibHi][iRotHi]*  */
									H2popHi*
									(H2L->Emis->Aul*
									/* escape and destruction */
									 (H2L->Emis->Pesc + H2L->Emis->Pelec_esc + H2L->Emis->Pdest) +
									 /* induced emission down */
									 H2L->Emis->pump * H2L->Lo->g/H2gHi);

								/* radiative rates [cm-3 s-1] from electronic excited states to X  */
								H2_X_rate_from_elec_excited[iVibLo][iRotLo] += rate_one;
							}
							++H2L;
						}
					}

					ASSERT( H2popHi >= 0. && H2popHi <= hmi.H2_total );

					/* this is total pop in this vibration state */
					pops_per_vib[iElecHi][iVibHi] += H2popHi;

					/* ======================= POPS EXCITED ELEC CONVERGE LOOP =====================*/
				}/* end excit electronic state rot pops loop */

				if(mole.nH2_TRACE >= mole.nH2_trace_full) 
					fprintf(ioQQQ,"\t%.2e",pops_per_vib[iElecHi][iVibHi]/hmi.H2_total);

				/* total pop in each electronic state */
				pops_per_elec[iElecHi] += pops_per_vib[iElecHi][iVibHi];
			}/* end excit electronic state all vibration loop */

			/* end excited electronic pops loop */
			if(mole.nH2_TRACE >= mole.nH2_trace_full) 
				fprintf(ioQQQ,"\n");
		} /* end loop over all electronic excited states */
		/*fprintf(ioQQQ,"DEBUG\t%.3e\t%.3e\n",
			H2_X_rate_from_elec_excited[0][0],
			pops_per_elec[0] );*/
		/* above set pops of excited electronic levels and found rates between them and X - 
		 * now solve highly excited levels within the X state by back-substitution */
		/* these will do convergence check */
		PopChgMaxOld_relative = PopChgMax_relative;
		PopChgMaxOld_total = PopChgMax_total;
		PopChgMax_relative = 0.;
		PopChgMax_total = 0.;
		iElec = 0;
		iElecHi = 0;
		iRotMaxChng_relative =-1;
		iVibMaxChng_relative = -1;
		iRotMaxChng_total =-1;
		iVibMaxChng_total = -1;
		popold_relative = 0.;
		popnew_relative = 0.;
		popold_total = 0.;
		popnew_total = 0.;

		/* now evaluate total rates for all levels within X */
		for( nEner=0; nEner<nLevels_per_elec[0]; ++nEner )
		{
			/* array of energy sorted indices within X */
			ip = H2_ipX_ener_sort[nEner];
			iVib = ipVib_H2_energy_sort[ip];
			iRot = ipRot_H2_energy_sort[ip];

			realnum H2stat = H2_stat[0][iVib][iRot];
			double H2boltz = H2_Boltzmann[0][iVib][iRot];

			/* these will be total rates into and out of the level */
			double col_rate_in = 0.;
			double col_rate_out = 0.;
			for( ipLo=0; ipLo<nEner; ++ipLo )
			{
				ip = H2_ipX_ener_sort[ipLo];
				iVibLo = ipVib_H2_energy_sort[ip];
				iRotLo = ipRot_H2_energy_sort[ip];

				/* this is rate from this level down to lower level, units s-1 */
				col_rate_out += H2_X_coll_rate[nEner][ipLo];
				/*if( nEner==288 ) fprintf(ioQQQ,"DEBUG %3li %3li %.3e\n", nEner,ipLo,H2_X_coll_rate[nEner][ipLo]);*/

				/* inverse, rate up, cm-3 s-1 */
				collup = (realnum)(H2_old_populations[0][iVibLo][iRotLo] * H2_X_coll_rate[nEner][ipLo] *	
					H2stat / H2_stat[0][iVibLo][iRotLo] *
					H2boltz / SDIV( H2_Boltzmann[0][iVibLo][iRotLo] ) );

				col_rate_in += collup;
			}

			for( ipHi=nEner+1; ipHi<nLevels_per_elec[0]; ++ipHi )
			{
				double colldn;
				ip = H2_ipX_ener_sort[ipHi];
				iVibHi = ipVib_H2_energy_sort[ip];
				iRotHi = ipRot_H2_energy_sort[ip];

				/* this is rate from this level up to higher level, units s-1 */
				col_rate_out += H2_X_coll_rate[ipHi][nEner] * 
					H2_stat[0][iVibHi][iRotHi] / H2stat *
					(realnum)(H2_Boltzmann[0][iVibHi][iRotHi] / SDIV( H2boltz ) );

				/* rate down from higher level, cm-3 s-1 */
				colldn = H2_old_populations[0][iVibHi][iRotHi] * H2_X_coll_rate[ipHi][nEner];
				/*if( nEner==288 ) fprintf(ioQQQ,"DEBUG %3li %3li %.3e\n", nEner,ipHi,H2_X_coll_rate[ipHi][nEner] * 
					H2_stat[0][iVibHi][iRotHi] / H2_stat[0][iVib][iRot] *
					(realnum)(H2_Boltzmann[0][iVibHi][iRotHi] /
					SDIV( H2_Boltzmann[0][iVib][iRot] ) ));*/

				col_rate_in += colldn;
			}

			H2_col_rate_in[iVib][iRot] = col_rate_in;
			H2_col_rate_out[iVib][iRot] = col_rate_out;
		}

		/* =======================INSIDE X POPULATIONS CONVERGE LOOP =====================*/
		/* begin solving for X by back-substitution
		 * this is the main loop that determines H2_populations within X 
		 * units of all rates in are cm-3 s-1, all rates out are s-1  
		 * nLevels_per_elec is number of levels within electronic 0 - so nEner is one
		 * beyond end of array here - but will be decremented at start of loop 
		 * this starts at the highest energy wihtin X and moves down to lower energies */
		nEner = nLevels_per_elec[0];
		while( (--nEner) >= nXLevelsMatrix )
		{

			/* array of energy sorted indices within X - we are moving down
			 * starting from highest level within X */
			ip = H2_ipX_ener_sort[nEner];
			iVib = ipVib_H2_energy_sort[ip];
			iRot = ipRot_H2_energy_sort[ip];

			if( nEner+1 < nLevels_per_elec[0] )
				ASSERT( H2_energies[H2_ipX_ener_sort[nEner]] < H2_energies[H2_ipX_ener_sort[nEner+1]] );

			/* >>chng 05 apr 30,GS, Instead of hmi.H2_total, the specific populations are used because high levels have much less
			 * populations than ground levels which consists most of the H2 population.
			 * only do this if working level is not v=0, J=0, 1 */ 
			if( nEner >1 )
			{
				H2_col_rate_out[iVib][iRot] += 
					/* H2 grain interactions
					 * rate (s-1) all v,J levels go to 0 or 1 preserving spin */
					(realnum)(hmi.rate_grain_h2_op_conserve);

				/* this goes into v=0, and J=0 or 1 depending on whether initial
				 * state is ortho or para */
				H2_col_rate_in[0][H2_lgOrtho[0][iVib][iRot]] += 
					/* H2 grain interactions
					 * rate (cm-3 s-1) all v,J levels go to 0 or 1 preserving spin,
					 * in above lgOrtho says whether should go to 0 or 1 */
					(realnum)(hmi.rate_grain_h2_op_conserve*H2_old_populations[0][iVib][iRot]);
			}
			else if( nEner == 1 )
			{
				/* this is special J=1 to J=0 collision, which is only fast at
				 * very low grain temperatures */
				H2_col_rate_out[0][1] += 
					/* H2 grain interactions
					 * H2 ortho - para conversion on grain surface,
					 * rate (s-1) all v,J levels go to 0 or 1, preserving nuclear spin */
					(realnum)(hmi.rate_grain_h2_J1_to_J0);

				H2_col_rate_in[0][0] += 
					/* H2 grain interactions
					 * H2 ortho - para conversion on grain surface,
					 * rate (s-1) all v,J levels go to 0 or 1, preserving nuclear spin */
					(realnum)(hmi.rate_grain_h2_J1_to_J0 *H2_old_populations[0][0][1]);
			}

			/* will become rate (cm-3 s-1) other levels have radiative transitions to here */
			H2_rad_rate_in[iVib][iRot] = 0.;
			H2_rad_rate_out[0][iVib][iRot] = 0.;

			/* the next two account for the Solomon process, 
			 * the first is the sum of decays from electronic excited into X
			 * second is X going into all excited electronic states 
			 * units cm-3 s-1 */
			H2_rad_rate_in[iVib][iRot] += H2_X_rate_from_elec_excited[iVib][iRot];

			/* radiative & cosmic ray rates [s-1] to electronic excited states from X only vibration and rot */
			H2_rad_rate_out[0][iVib][iRot] += H2_X_rate_to_elec_excited[iVib][iRot];

			/* now sum over states within X which are higher than current state */
			iElecHi = 0;
			for( ipHi = nEner+1; ipHi<nLevels_per_elec[0]; ++ipHi )
			{
				ip = H2_ipX_ener_sort[ipHi];
				iVibHi = ipVib_H2_energy_sort[ip];
				iRotHi = ipRot_H2_energy_sort[ip];
				/* =======================INSIDE POPULATIONS CONVERGE LOOP =====================*/
				/* the rate we enter this state from more highly excited states within X
				 * by radiative decays, which have delta J = 0 or 2 */
				/* note test on vibration is needed - iVibHi<iVib, energy order ok and space not allocated */
				/* >>chng 05 feb 07, tried to use use lgH2_line_exists but cant */
				if( ( abs(iRotHi-iRot) ==2 || iRotHi==iRot ) && (iVib <= iVibHi) &&
					H2Lines[iElecHi][iVibHi][iRotHi][iElec][iVib][iRot].ipCont > 0 )
				{
					double rateone;
					rateone =
						H2_old_populations[iElecHi][iVibHi][iRotHi]*
						 H2Lines[iElecHi][iVibHi][iRotHi][iElec][iVib][iRot].Emis->Aul*
						(H2Lines[iElecHi][iVibHi][iRotHi][iElec][iVib][iRot].Emis->Pesc + 
						 H2Lines[iElecHi][iVibHi][iRotHi][iElec][iVib][iRot].Emis->Pelec_esc + 
						 H2Lines[iElecHi][iVibHi][iRotHi][iElec][iVib][iRot].Emis->Pdest);
					ASSERT( rateone >=0 );

					/* units cm-3 s-1 */
					H2_rad_rate_in[iVib][iRot] += rateone;
				}
			}

			/* =======================INSIDE POPULATIONS X CONVERGE LOOP =====================*/
			/* we now have total rate this state is populated from above, now get rate
			 * this state interacts with levels that are below */
			iElecLo = 0;
			for( ipLo = 0; ipLo<nEner; ++ipLo )
			{
				ip = H2_ipX_ener_sort[ipLo];
				iVibLo = ipVib_H2_energy_sort[ip];
				iRotLo = ipRot_H2_energy_sort[ip];
				/* radiative interactions between this level and lower levels */
				/* the test on vibration is needed - the energies are ok but the space does not exist */
				/* >>chng 05 feb 07, can't use lgH2_line_exists */
				if( ( abs(iRotLo-iRot) == 2 || iRotLo == iRot )  && (iVibLo <= iVib) && 
					H2Lines[iElec][iVib][iRot][iElecLo][iVibLo][iRotLo].ipCont > 0 ) 
				{
					H2_rad_rate_out[0][iVib][iRot] +=
						H2Lines[iElec][iVib][iRot][iElecLo][iVibLo][iRotLo].Emis->Aul*
						(H2Lines[iElec][iVib][iRot][iElecLo][iVibLo][iRotLo].Emis->Pesc + 
						H2Lines[iElec][iVib][iRot][iElecLo][iVibLo][iRotLo].Emis->Pelec_esc + 
						H2Lines[iElec][iVib][iRot][iElecLo][iVibLo][iRotLo].Emis->Pdest);
				}
			}
			/* =======================INSIDE X POPULATIONS CONVERGE LOOP =====================*/

			/* we now have the total rates into and out of this level, get its population 
			 * units cm-3 */
			H2_populations[iElec][iVib][iRot] = 
				(H2_col_rate_in[iVib][iRot]+ H2_rad_rate_in[iVib][iRot]+H2_X_source[nEner]) / 
				SDIV(H2_col_rate_out[iVib][iRot]+H2_rad_rate_out[0][iVib][iRot]+H2_X_sink[nEner]);

			ASSERT( H2_populations[iElec][iVib][iRot] >= 0.  );
		}
		/* >>chng 05 may 10, move to following back substitution part within X */
		/* =======================INSIDE POPULATIONS CONVERGE LOOP =====================*/
		/* now do lowest levels H2_populations with matrix, 
		 * these should be collisionally dominated */
		if( nXLevelsMatrix )
		{
			H2_Level_low_matrix(
				/* the total abundance - frac_matrix is fraction of pop that was in these
				 * levels the last time this was done */
				hmi.H2_total * (realnum)frac_matrix );
		}
		iElecHi = 0;
		if(mole.nH2_TRACE >= mole.nH2_trace_full) 
		{
			fprintf(ioQQQ," Rel pop(e=%li)" ,iElecHi);
		}

		/* find ortho and para densites, sum of pops in each vibration */
		/* this will become total pop is X, which will be renormed to equal hmi.H2_total */
		pops_per_elec[0] = 0.;
		for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
		{
			double sumv;
			sumv = 0.;
			pops_per_vib[0][iVibHi] = 0.;

			for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
			{
				pops_per_elec[0] += H2_populations[iElecHi][iVibHi][iRotHi];
				sumv += H2_populations[iElecHi][iVibHi][iRotHi];
				pops_per_vib[0][iVibHi] += H2_populations[iElecHi][iVibHi][iRotHi];
			}
			/* print sum of H2_populations in each vibration if trace on */
			if(mole.nH2_TRACE >= mole.nH2_trace_full) 
				fprintf(ioQQQ,"\t%.2e",sumv/hmi.H2_total);
		}
		ASSERT( pops_per_elec[0] > SMALLFLOAT );
		/* =======================INSIDE POPULATIONS CONVERGE LOOP =====================*/
		if(mole.nH2_TRACE >= mole.nH2_trace_full) 
		{
			fprintf(ioQQQ,"\n");
			/* print the ground vibration state */
			fprintf(ioQQQ," Rel pop(0,J)");
			iElecHi = 0;
			iVibHi = 0;
			for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
			{
				fprintf(ioQQQ,"\t%.2e",H2_populations[iElecHi][iVibHi][iRotHi]/hmi.H2_total);
			}
			fprintf(ioQQQ,"\n");
		}

		/* now find population in states done with matrix - this is only used to pass
		 * to matrix solver */
		iElec = 0;
		sum_pops_matrix = 0.;
		ip =0;
		for( i=0; i<nXLevelsMatrix; ++i )
		{
			ip = H2_ipX_ener_sort[i];
			iVib = ipVib_H2_energy_sort[ip];
			iRot = ipRot_H2_energy_sort[ip];
			sum_pops_matrix += H2_populations[iElec][iVib][iRot];
		}
		/* =======================INSIDE POPULATIONS CONVERGE LOOP =====================*/
		/* this is self consistent since pops_per_elec[0] came from current soln,
		* as did the matrix.  pops will be renormalized by results from the chemistry
		* a few lines down */
		frac_matrix = sum_pops_matrix / SDIV(pops_per_elec[0]);

		/* assuming that all H2 population is in X, this is the
		 * ratio of H2 that came out of the chemistry network to what we just obtained -
		 * we need to multiply the pops by renorm to agree with the chemistry,
		 * this routine does not alter hmi.H2_total, but does change pops_per_elec */
		H2_renorm_conserve = hmi.H2_total/ SDIV(pops_per_elec[0]);
		/*pops_per_elec[0] = hmi.H2_total;*/

		/* renormalize H2_populations  - the H2_populations were updated by renorm when routine entered, 
		 * before pops determined - but population determinations above do not have a sum rule on total
		 * population - this renorm is to preserve total population */
		H2_sum_excit_elec_den = 0.;
		for( iElecHi=0; iElecHi<mole.n_h2_elec_states; ++iElecHi )
		{
			pops_per_elec[iElecHi] *= H2_renorm_conserve;
			if( iElecHi > 0 )
				H2_sum_excit_elec_den += pops_per_elec[iElecHi];

			for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
			{
				for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
				{
					H2_populations[iElecHi][iVibHi][iRotHi] *= H2_renorm_conserve;
					/* =======================INSIDE POPULATIONS CONVERGE LOOP =====================*/
				}
			}
		}

		/* this loop first checks for largest changes in populations, to determine whether
		 * we have converged, then updates the population array with a new value,
		 * which may be a mean of old and new
		 * update populations check convergence converged */
		sumold = 0.;
		for( iElec=0; iElec<mole.n_h2_elec_states; ++iElec )
		{
			for( iVib=0; iVib<=h2.nVib_hi[iElec]; ++iVib )
			{
				for( iRot=h2.Jlowest[iElec]; iRot<=h2.nRot_hi[iElec][iVib]; ++iRot )
				{
					double rel_change;
					/* keep track of largest relative change in H2_populations to
					 * determines convergence */
					if( fabs(H2_populations[iElec][iVib][iRot] - 
						H2_old_populations[iElec][iVib][iRot])/
						/* on first call some very high J states can have zero pop ,
						 * hence the SDIV, will retain sign for checks on oscilations,
						 * hence the fabs */
						SDIV(H2_populations[iElec][iVib][iRot]) > fabs(PopChgMax_relative) &&
						/* >>chng 03 jul 19, this had simply been H2_populations > SMALLFLOAT,
						* change to relative pops > 1e-15, spent too much time converging
						* levels at pops = 1e-37 */
						/* >>chng 03 dec 27, from rel pop 1e-15 to 1e-6 since converging heating will
						* be main convergence criteria check convergence */
						/*H2_populations[iElecHi][iVibHi][iRotHi]/SDIV(hmi.H2_total)>1e-15 )*/
						H2_populations[iElec][iVib][iRot]/SDIV(hmi.H2_total)>1e-6 )
					{
						PopChgMax_relative = 
							(H2_populations[iElec][iVib][iRot] - 
							H2_old_populations[iElec][iVib][iRot])/
							SDIV(H2_populations[iElec][iVib][iRot]);
						iRotMaxChng_relative = iRot;
						iVibMaxChng_relative = iVib;
						popold_relative = H2_old_populations[iElec][iVib][iRot];
						popnew_relative = H2_populations[iElec][iVib][iRot];
					}
					/* >>chng 05 feb 08, add largest rel change in total, this will be converged
					 * down to higher accuracy than above 
					 * keep track of largest change in H2_populations relative to total H2 to
					 * determine convergence check convergence */
					rel_change = (H2_populations[iElec][iVib][iRot] - 
						      H2_old_populations[iElec][iVib][iRot])/SDIV(hmi.H2_total);
					/*  retain sign for checks on oscillations hence the fabs */
					if( fabs(rel_change) > fabs(PopChgMax_total) )
					{
						PopChgMax_total = rel_change;
						iRotMaxChng_total = iRot;
						iVibMaxChng_total = iVib;
						popold_total = H2_old_populations[iElec][iVib][iRot];
						popnew_total = H2_populations[iElec][iVib][iRot];
					}

					kase = -1;
					/* update populations - we used the old populations to update the
					 * current new populations - will do another iteration if they changed
					 * by much.  here old populations are updated for next sweep through molecule */
					/* pop oscillations have occurred - use small changes */
					/* >>chng 04 may 10, turn this back on - now with min on how small frac new
					 * can become */
					rel_change = fabs( H2_old_populations[iElec][iVib][iRot] - 
						H2_populations[iElec][iVib][iRot] )/
						SDIV( H2_populations[iElec][iVib][iRot] );

					/* this branch very large changes, use mean of logs but onlly if both are positive*/
					if( rel_change > 3. &&
						H2_old_populations[iElec][iVib][iRot]*H2_populations[iElec][iVib][iRot]>0  )
					{
						/* large changes or oscillations - take average in the log */
						H2_old_populations[iElec][iVib][iRot] = pow( 10. , 
							log10(H2_old_populations[iElec][iVib][iRot])/2. +
							log10(H2_populations[iElec][iVib][iRot])/2. );
						kase = 2;
					}

					/* modest change, use means of old and new */
					else if( rel_change> 0.1 )
					{
						realnum frac_old=0.25f;
						/* large changes or oscillations - take average */
						H2_old_populations[iElec][iVib][iRot] = 
							frac_old*H2_old_populations[iElec][iVib][iRot] +
							(1.f-frac_old)*H2_populations[iElec][iVib][iRot];
						kase = 3;
					}
					else
					{
						/* small changes, use new value */
						H2_old_populations[iElec][iVib][iRot] = 
							H2_populations[iElec][iVib][iRot];
						kase = 4;
					}
					sumold += H2_old_populations[iElec][iVib][iRot];
				}
			}
		}
		/* will renormalize so that total population is correct */
		H2_renorm_conserve_init = hmi.H2_total/sumold;

		/* renormalize H2_populations  - the H2_populations were updated by renorm when routine entered, 
		 * before pops determined - but population determinations above do not have a sum rule on total
		 * population - this renorm is to preserve total population */
		for( iElecHi=0; iElecHi<mole.n_h2_elec_states; ++iElecHi )
		{
			for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
			{
				for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
				{
					H2_old_populations[iElecHi][iVibHi][iRotHi] *= H2_renorm_conserve_init;
					/* =======================INSIDE POPULATIONS CONVERGE LOOP =====================*/
				}
			}
		}
		/* get current ortho-para ratio, will be used as test on convergence */
		iElecHi = 0;
		h2.ortho_density = 0.;
		h2.para_density = 0.;
		H2_den_s = 0.;
		H2_den_g = 0.;

		for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
		{
			for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
			{
				/* find current population in H2s and H2g */
				if( energy_wn[0][iVibHi][iRotHi]> ENERGY_H2_STAR )
				{
					H2_den_s += H2_populations[iElecHi][iVibHi][iRotHi];
				}
				else
				{
					H2_den_g += H2_populations[iElecHi][iVibHi][iRotHi];
				}
				if( H2_lgOrtho[iElecHi][iVibHi][iRotHi] )
				{
					h2.ortho_density += H2_populations[iElecHi][iVibHi][iRotHi];
				}
				else
				{
					h2.para_density += H2_populations[iElecHi][iVibHi][iRotHi];
				}
				/* =======================INSIDE POPULATIONS CONVERGE LOOP =====================*/
			}
		}
		/* these will be used to determine whether solution has converged */
		ortho_para_older = ortho_para_old;
		ortho_para_old = ortho_para_current;
		ortho_para_current = h2.ortho_density / SDIV( h2.para_density );

		/* this will be evaluated in call to routine that follows - will check
		 * whether this has converged */
		old_solomon_rate = hmi.H2_Solomon_dissoc_rate_BigH2_H2g;

		/* >>chng 05 jul 24, break code out into separate routine for clarify
		 * located in mole_h2_etc.c - true says to only do Solomon rate */
		H2_Solomon_rate(  );

		/* are changes too large? must decide whether population shave converged,
		 * will check whether H2_populations themselves have changed by much,
		 * but also change in heating by collisional deexcitation is stable */
		HeatChangeOld = HeatChange;
		HeatChange = old_HeatH2Dexc_BigH2 - hmi.HeatH2Dexc_BigH2;
		{
			static long int loop_h2_oscil=-1;
			/* check whether pops are oscillating, as evidenced by change in
			 * heating changing sign */
			if( loop_h2_pops>2 && (
				(HeatChangeOld*HeatChange<0. ) ||
				(PopChgMax_relative*PopChgMaxOld_relative<0. ) ) )
			{
				lgH2_pops_oscil = true;
				if( loop_h2_pops > 6 )
				{
					loop_h2_oscil = loop_h2_pops;
					lgH2_pops_ever_oscil = true;
					/* make this smaller in attempt to damp out oscillations,
					 * but don't let get too small*/
					frac_new_oscil *= 0.8f;
					frac_new_oscil = MAX2( frac_new_oscil , 0.1f);
					++n_pop_oscil;
				}
			}
			else
			{
				lgH2_pops_oscil = false;
				/* turn off flag if no oscillations for a while */
				if( loop_h2_pops -  loop_h2_oscil > 4 )
				{
					frac_new_oscil = 1.f;
					lgH2_pops_ever_oscil = false;
				}
			}
		}

		/* reevaluate heating - cooling if H2 molecule is significant source or either,
		 * since must have stable heating cooling rate */
		old_HeatH2Dexc_BigH2 = hmi.HeatH2Dexc_BigH2;
		if(fabs(hmi.HeatH2Dexc_BigH2)/thermal.ctot > conv.HeatCoolRelErrorAllowed/10. ||
			hmi.HeatH2Dexc_BigH2==0. )
			H2_Cooling("H2lup");

		/* begin check on whether solution is converged */
		lgConv_h2_soln = true;
		lgPopsConv_total = true;
		lgPopsConv_relative = true;
		lgHeatConv = true;
		lgSolomonConv = true;
		lgOrthoParaRatioConv = true;

		/* these are all the convergence tests 
		 * check convergence converged */
		if( fabs(PopChgMax_relative)>converge_pops_relative )
		{
			/*lgPopsConv = (fabs(PopChgMax_relative)<=0.1);*/
			lgConv_h2_soln = false;
			lgPopsConv_relative = false;
			/* >>chng 04 sep 08, set quant_new to new chng max gs */
			/*quant_old = PopChgMax_relative;*/
			quant_old = PopChgMaxOld_relative;
			/*quant_new = 0.;*/
			quant_new = PopChgMax_relative;

			strcpy( chReason , "rel pops changed" );
		}

		/* check largest change in a level population relative to total h2 
		 * population convergence converged check */
		else if( fabs(PopChgMax_total)>converge_pops_total)
		{
			lgConv_h2_soln = false;
			lgPopsConv_total = false;
			/* >>chng 04 sep 08, set quant_new to new chng max gs */
			/*quant_old = PopChgMax_relative;*/
			quant_old = PopChgMaxOld_total;
			/*quant_new = 0.;*/
			quant_new = PopChgMax_total;

			strcpy( chReason , "tot pops changed" );
		}

		/* >>chng 04 apr 30, look at change in ortho-para ratio, also that is not
		 * oscillating */
		/* >>chng 04 dec 15, only look at change, and don't make allowed change so tiny -
		 * these were attempts at fixing problems that were due to shielding not thin*/
		else if( fabs(ortho_para_current-ortho_para_old) / SDIV(ortho_para_current)> converge_ortho_para )
		/* else if( fabs(ortho_para_current-ortho_para_old) / SDIV(ortho_para_current)> 1e-3 
			&& (ortho_para_current-ortho_para_old)*(ortho_para_old-ortho_para_older)>0. )*/
		{
			lgConv_h2_soln = false;
			lgOrthoParaRatioConv = false;
			quant_old = ortho_para_old;
			quant_new = ortho_para_current;
			strcpy( chReason , "ortho/para ratio changed" );
		}
		/* >>chng 04 dec 16, reduce error allowed fm /5 to /2, to be similar to 
		 * logic in conv_base */
		else if( !thermal.lgTemperatureConstant &&
			fabs(hmi.HeatH2Dexc_BigH2-old_HeatH2Dexc_BigH2)/MAX2(thermal.ctot,thermal.htot) > 
			conv.HeatCoolRelErrorAllowed/2.
			/* >>chng 04 may 09, do not check on error in heating if constant temperature */
			/*&& !(thermal.lgTemperatureConstant || phycon.te <= phycon.TEMP_LIMIT_LOW  )*/ )
		{
			/* default on HeatCoolRelErrorAllowed is 0.02 */
			/*lgHeatConv = (fabs(hmi.HeatH2Dexc_BigH2-old_HeatH2Dexc_BigH2)/thermal.ctot <=
			 * conv.HeatCoolRelErrorAllowed/5.);*/
			lgConv_h2_soln = false;
			lgHeatConv = false;
			quant_old = old_HeatH2Dexc_BigH2/MAX2(thermal.ctot,thermal.htot);
			quant_new = hmi.HeatH2Dexc_BigH2/MAX2(thermal.ctot,thermal.htot);
			strcpy( chReason , "heating changed" );
			/*fprintf(ioQQQ,"DEBUG old new trip \t%.4e \t %.4e\n",
				old_HeatH2Dexc_BigH2,
				hmi.HeatH2Dexc_BigH2);*/
		}

		/* check on Solomon rate,
		 * >>chng 04 aug 28, do not do this check if induced processes are disabled,
		 * since Solomon process is then irrelevant */
		/* >>chng 04 sep 21, GS*/
		else if( rfield.lgInducProcess && 
			/* this is check that H2 abundance has not been set - if it has been
			 * then we don't care what the Solomon rate is doing */ 
			 hmi.H2_frac_abund_set==0 &&
			 /*>>chng 05 feb 10, rather than checking change in Solomon relative to Solomon,
			  * check it relative to total h2 destruction rate */
			fabs( hmi.H2_Solomon_dissoc_rate_BigH2_H2g - old_solomon_rate)/SDIV(hmi.H2_rate_destroy) > 
			conv.EdenErrorAllowed/5.)
		{
			lgConv_h2_soln = false;
			lgSolomonConv = false;
			quant_old = old_solomon_rate;
			quant_new = hmi.H2_Solomon_dissoc_rate_BigH2_H2g;
			strcpy( chReason , "Solomon rate changed" );
		}

		/* did we pass all the convergence test */
		if( !lgConv_h2_soln )
		{
			/* this branch H2 H2_populations within X are not converged,
			 * print diagnostic */

			if( PRT_POPS || mole.nH2_TRACE >=mole.nH2_trace_iterations )
			{
				/*fprintf(ioQQQ,"temppp\tnew\t%.4e\tnew\t%.4e\t%.4e\n",
					hmi.HeatH2Dexc_BigH2,
					old_HeatH2Dexc_BigH2,
					fabs(hmi.HeatH2Dexc_BigH2-old_HeatH2Dexc_BigH2)/thermal.ctot );*/
				fprintf(ioQQQ,"    loop %3li no conv oscl?%c why:%s ",
					loop_h2_pops,
					TorF(lgH2_pops_ever_oscil),
					chReason );
				if( !lgPopsConv_relative )
					fprintf(ioQQQ," PopChgMax_relative:%.4e v:%li J:%li old:%.4e new:%.4e",
					PopChgMax_relative,
					iVibMaxChng_relative,
					iRotMaxChng_relative ,
					popold_relative ,
					popnew_relative );
				else if( !lgPopsConv_total )
					fprintf(ioQQQ," PopChgMax_total:%.4e v:%li J:%li old:%.4e new:%.4e",
					PopChgMax_total,
					iVibMaxChng_total,
					iRotMaxChng_total ,
					popold_total ,
					popnew_total );
				else if( !lgHeatConv )
					fprintf(ioQQQ," heat:%.4e old:%.4e new:%.4e",
					(hmi.HeatH2Dexc_BigH2-old_HeatH2Dexc_BigH2)/MAX2(thermal.ctot,thermal.htot), 
					quant_old , 
					quant_new);
				/* Solomon rate changed */ 
				else if( !lgSolomonConv )
					fprintf(ioQQQ," d(sol rate)/tot dest\t%2e",(old_solomon_rate - hmi.H2_Solomon_dissoc_rate_BigH2_H2g)/SDIV(hmi.H2_rate_destroy));
				else if( !lgOrthoParaRatioConv )
					fprintf(ioQQQ," current, old, older ratios are %.4e %.4e %.4e",
					ortho_para_current , ortho_para_old, ortho_para_older );
				else
					TotalInsanity();
				fprintf(ioQQQ,"\n");
			}
		}
		/* end convergence criteria */

		/*fprintf(ioQQQ,"DEBUG h2 heat\t%3li\t%.2f\t%.4e\t%.4e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n",
			loop_h2_pops,
			fnzone,
			phycon.te,
			dense.eden,
			hmi.HeatH2Dexc_BigH2,
			hmi.HeatH2Dexc_BigH2/thermal.ctot ,
			hmi.H2_total,
			H2_renorm_chemistry ,
			H2_renorm_conserve,
			hmi.H2_H2g_to_H2s_rate_BigH2);*/
		if( trace.nTrConvg >= 5 )
		{
			fprintf( ioQQQ, 
				"     H2 5lev %li Conv?%c",
				loop_h2_pops ,
				TorF(lgConv_h2_soln) );

			if( fabs(PopChgMax_relative)>0.1 )
				fprintf(ioQQQ," pops, rel chng %.3e",PopChgMax_relative);
			else
				fprintf(ioQQQ," rel heat %.3e rel chng %.3e H2 heat/cool %.2e",
					hmi.HeatH2Dexc_BigH2/thermal.ctot ,
					fabs(hmi.HeatH2Dexc_BigH2-old_HeatH2Dexc_BigH2)/thermal.ctot ,
					hmi.HeatH2Dexc_BigH2/thermal.ctot);

			fprintf( ioQQQ, 
				" Oscil?%c Ever Oscil?%c",
				TorF(lgH2_pops_oscil) ,
				TorF(lgH2_pops_ever_oscil) );
			if( lgH2_pops_ever_oscil )
				fprintf(ioQQQ," frac_new_oscil %.4f",frac_new_oscil);
			fprintf(ioQQQ,"\n");
		}

		if( mole.nH2_TRACE >= mole.nH2_trace_full ) 
		{
			fprintf(ioQQQ,
			"H2 loop\t%li\tkase pop chng\t%i\tchem renorm fac\t%.4e\tortho/para ratio:\t%.3e\tfrac of pop in matrix: %.3f\n",
			loop_h2_pops,
			kase,
			H2_renorm_chemistry,
			h2.ortho_density / h2.para_density ,
			frac_matrix);

			/* =======================INSIDE POPULATIONS CONVERGE LOOP =====================*/
			if( iVibMaxChng_relative>=0 && iRotMaxChng_relative>=0 && PopChgMax_relative>1e-10 )
				fprintf(ioQQQ,
					"end loop %li H2 max rel chng=%.3e from %.3e to %.3e at v=%li J=%li\n\n",
					loop_h2_pops,
					PopChgMax_relative , 
					H2_old_populations[0][iVibMaxChng_relative][iRotMaxChng_relative],
					H2_populations[0][iVibMaxChng_relative][iRotMaxChng_relative],
					iVibMaxChng_relative , iRotMaxChng_relative
					);
		}
	}
	/* =======================END POPULATIONS CONVERGE LOOP =====================*/

	/* evaluate H2 rates over H2g and H2s for use in chemistry network */
	H2_gs_rates();

	/* >>chng 05 feb 08, do not print if we are in search phase */
	if( !lgConv_h2_soln && !conv.lgSearch )
	{
		conv.lgConvPops = false;
		strcpy( conv.chConvIoniz, "H2 pop cnv" );
		fprintf(ioQQQ,
			"  H2_LevelPops:  H2_populations not converged in %li tries; due to %s, old, new are %.4e %.4e, iteration %li zone %.2f.\n",
			loop_h2_pops, 
			chReason,
			quant_old,
			quant_new ,
			iteration , 
			fnzone );
		ConvFail("pops","H2");
	}

	/* loop over all possible lines and set H2_populations, 
	 * and quantities that depend on them */
	for( iElecHi=0; iElecHi<mole.n_h2_elec_states; ++iElecHi )
	{
		for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
		{
			for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
			{
				long int lim_elec_lo = 0;
				double H2popHi = H2_populations[iElecHi][iVibHi][iRotHi];
				// this will update Lo->Pop as well as Lo and Hi point to the same set of data structures
				// the loops below are OK since Lo->Pop is only accessed after Hi->Pop has been set.
				H2Lines[iElecHi][iVibHi][iRotHi][0][0][0].Hi->Pop = H2popHi;
				ASSERT( H2popHi >= 0. );
				realnum H2gHi = H2Lines[iElecHi][iVibHi][iRotHi][0][0][0].Hi->g;
				/* now the lower levels */
				/* NB - X is the only lower level considered here, since we are only 
				* concerned with excited electronic levels as a photodissociation process
				* code exists to relax this assumption - simply change following to iElecHi */
				for( iElecLo=0; iElecLo<=lim_elec_lo; ++iElecLo )
				{
					/* want to include all vibration states in lower level if different electronic level,
					 * but only lower vibration levels if same electronic level */
					long int nv = h2.nVib_hi[iElecLo];
					if( iElecLo==iElecHi )
						nv = iVibHi;
					for( iVibLo=0; iVibLo<=nv; ++iVibLo )
					{
						long nr = h2.nRot_hi[iElecLo][iVibLo];
						if( iElecLo==iElecHi && iVibHi==iVibLo )
							nr = iRotHi-1;

						mb6ci lgH2le = lgH2_line_exists.ptr(iElecHi,iVibHi,iRotHi,iElecLo,iVibLo,h2.Jlowest[iElecLo]);
						mt6i H2L = H2Lines.ptr(iElecHi,iVibHi,iRotHi,iElecLo,iVibLo,h2.Jlowest[iElecLo]);
						for( iRotLo=h2.Jlowest[iElecLo]; iRotLo<=nr; ++iRotLo )
						{
							if( *lgH2le++ )
							{
								/* following two heat exchange excitation, deexcitation */
								H2L->Coll.cool = 0.;
								H2L->Coll.heat = 0.;

								H2L->Emis->PopOpc = 
									H2L->Lo->Pop - H2popHi * H2L->Lo->g / H2gHi;

								/* number of photons in the line */
								H2L->Emis->phots = H2L->Emis->Aul * 
									(H2L->Emis->Pesc + H2L->Emis->Pelec_esc) * 
									 H2popHi; 

								/* intensity of line */
								H2L->Emis->xIntensity = 
								H2L->Emis->phots *
								H2L->EnergyErg;

								if( iElecHi==0 )
								{
									/** \todo	2	- put H2Lines in outward beams in RT_diffuse */
									/* the ground electronic state, most excitations are not direct pumping 
									 * (rather indirect, which does not count for ColOvTot) */
									H2L->Emis->ColOvTot = 1.;
								}
								else
								{
									/* these are excited electronic states, mostly pumped, except for supras */
									/** \todo	2	put supra thermal excitation into excitation of electronic bands */
									H2L->Emis->ColOvTot = 0.;
								}
							}
							++H2L;
						}
					}
				}
			}
		}
	}	

	/* add up H2 + hnu => 2H, continuum photodissociation,
	 * this is not the Solomon process, true continuum */
	/* >>chng 05 jun 16, GS, add dissociation to triplet states*/
	hmi.H2_photodissoc_BigH2_H2s = 0.;
	hmi.H2_photodissoc_BigH2_H2g = 0.;
	hmi.H2_tripletdissoc_H2s =0.;
	hmi.H2_tripletdissoc_H2g =0.;
	hmi.H2_BigH2_H2g_av = 0.;
	hmi.H2_BigH2_H2s_av = 0.;
	/* >>chng 05 jul 20, GS, add dissociation by H2 g and H2s*/
	hmi.Average_collH2s_dissoc = 0.;
	hmi.Average_collH2g_dissoc = 0.;

	iElec = 0;
	H2_BigH2_H2s = 0.;
	H2_BigH2_H2g = 0.;
	hmi.H2g_BigH2 =0;
	hmi.H2s_BigH2 = 0;
	hmi.H2_total_BigH2 =0;
	hmi.H2g_LTE_bigH2 =0.;
	hmi.H2s_LTE_bigH2 = 0.;
	/* >>chng 05 oct 20, no need to reset this var here */
	/*exp_disoc =  sexp(H2_DissocEnergies[0]/phycon.te_wn);*/

	/* >>chng 05 sep 12, TE, define a cutoff wavelength of 800 Angstrom 
	 * this is chosen as the cross sections given by 
	 *>>refer	H2	photo cs	Allison, A.C. & Dalgarno, A. 1969, Atomic Data, 1, 91 
	 * show a sharp decline in the cross section*/
	{
		static long ip_cut_off = -1;
		if( ip_cut_off < 0 )
		{
			/* one-time initialization of this pointer */
			ip_cut_off = ipoint( 1.14 );
		}

		/* >>chng 05 sep 12, TE, assume all H2s is at 2.5 eV
		 * the dissociation threshold is at 1.07896 Rydberg*/
		flux_accum_photodissoc_BigH2_H2s = 0;
		ip_H2_level = ipoint( 1.07896 - 2.5 / EVRYD);
		for( i= ip_H2_level; i < ip_cut_off; ++i )
		{
			flux_accum_photodissoc_BigH2_H2s += ( rfield.flux[0][i-1] + rfield.ConInterOut[i-1]+ 
				rfield.outlin[0][i-1]+ rfield.outlin_noplot[i-1]  );
		}

		/* sum over all levels to obtain s and g populations and dissociation rates */
		for( iVib=0; iVib<=h2.nVib_hi[iElec]; ++iVib )
		{
			for( iRot=h2.Jlowest[iElec]; iRot<=h2.nRot_hi[iElec][iVib]; ++iRot )
			{ 
				/* >>chng 05 mar 22, TE,  moved H2_photodissoc_BigH2_H2s in this statement and divide by the
				 *	density of H2s not total H2, we consider direct photodissociation only for H2s */
				/* >>chng 05 mar 22, TE, this should be for H2* rather than total */
				/* this is the total rate of direct photo-dissociation of excited electronic states into 
				 * the X continuum - this is continuum photodissociation, not the Solomon process */
				/* >>chng 03 sep 03, make sum of pops of excited states */
				if( energy_wn[0][iVib][iRot] > ENERGY_H2_STAR )
				{
					double arg_ratio;
					hmi.H2_photodissoc_BigH2_H2s += 
						H2_populations[iElec][iVib][iRot] * flux_accum_photodissoc_BigH2_H2s;

					/* cosmic ray & secondary electron excitation to triplets
					 * physics described where similar process is done for 
					 * big molecule 
					 *>>chng 07 apr 08, from 3 to 10 to better capture results 
					 * of Dalgarno et al 99 */
					hmi.H2_tripletdissoc_H2s += 
						H2_populations[iElec][iVib][iRot] * 10.f*secondaries.x12tot;

					/* sum of pops in levels in H2* for use in chemistry network */
					H2_BigH2_H2s += H2_populations[iElec][iVib][iRot];

					/* >>chng 05 jun 28, TE, determine average energy level in H2s */
					hmi.H2_BigH2_H2s_av += (H2_populations[iElec][iVib][iRot] * energy_wn[0][iVib][iRot]);

					/* >>chng 05 july 20, GS, collisional dissociation  by H2s, unit s-1*/
					hmi.Average_collH2s_dissoc += H2_populations[iElec][iVib][iRot] * H2_coll_dissoc_rate_coef_H2[iVib][iRot];

					/* >>chng 05 oct 17, GS, LTE populations of H2s*/
					arg_ratio = exp_disoc/SDIV(H2_Boltzmann[0][iVib][iRot]);
					if( arg_ratio > 0. )
					{
						/* >>chng 05 oct 21, GS, only add ratio if Boltzmann factor > 0 */
						hmi.H2s_LTE_bigH2 += H2_populations[0][iVib][iRot]*SAHA/SDIV(phycon.te32*arg_ratio)*
							(H2_stat[0][iVib][iRot]/(2.*2.))*3.634e-5;	
					}
				}
				else
				{
					double arg_ratio;
					/* >>chng 05 sep 12, TE, for H2g do the sum explicitly for every level*/
					flux_accum_photodissoc_BigH2_H2g = 0;
					/* this is the dissociation energy needed for the level*/
					ip_H2_level = ipoint( 1.07896 - energy_wn[0][iVib][iRot] * WAVNRYD);

					for( i= ip_H2_level; i < ip_cut_off; ++i )
					{
						flux_accum_photodissoc_BigH2_H2g += ( rfield.flux[0][i-1] + rfield.ConInterOut[i-1]+ 
							rfield.outlin[0][i-1]+ rfield.outlin_noplot[i-1] );
					}

					hmi.H2_photodissoc_BigH2_H2g += 
						H2_populations[iElec][iVib][iRot] * flux_accum_photodissoc_BigH2_H2g;


					/* cosmic ray & secondary electron excitation to triplets
					 * physics described where similar process is done for 
					 * big molecule 
					 *>>chng 07 apr 08, from 3 to 10 to better capture results 
					 * of Dalgarno et al 99 */
					hmi.H2_tripletdissoc_H2g += 
						H2_populations[iElec][iVib][iRot] * 10.f*secondaries.x12tot;

					/* sum of pops in levels in H2g for use in chemistry network */
					H2_BigH2_H2g += H2_populations[iElec][iVib][iRot];

					/* >>chng 05 jun 28, TE, determine average energy level in H2g */
					hmi.H2_BigH2_H2g_av += (H2_populations[iElec][iVib][iRot] * energy_wn[0][iVib][iRot]);

					/* >>chng 05 jul 20, GS, collisional dissociation  by H2s, unit s-1*/
					hmi.Average_collH2g_dissoc += H2_populations[iElec][iVib][iRot] * H2_coll_dissoc_rate_coef_H2[iVib][iRot];

					/* >>chng 05 oct 17, GS, LTE populations of H2g*/
					arg_ratio = exp_disoc/SDIV(H2_Boltzmann[0][iVib][iRot]);
					if( arg_ratio > 0. )
					{
						hmi.H2g_LTE_bigH2 += H2_populations[0][iVib][iRot]*(SAHA/SDIV(phycon.te32*arg_ratio)*
							(H2_stat[0][iVib][iRot]/(2.*2.))*3.634e-5);
					}
				}
			}
		}
	}
	hmi.H2g_BigH2 = (realnum)H2_BigH2_H2g;
	hmi.H2s_BigH2 = (realnum)H2_BigH2_H2s;
	hmi.H2_total_BigH2 =hmi.H2g_BigH2+hmi.H2s_BigH2;

	ASSERT( H2_BigH2_H2s > 0. );
	ASSERT( H2_BigH2_H2g > 0. );

	/* average energy in H2s */
	hmi.H2_BigH2_H2s_av = hmi.H2_BigH2_H2s_av / H2_BigH2_H2s;
	/* average energy in H2g */
	hmi.H2_BigH2_H2g_av = hmi.H2_BigH2_H2g_av / H2_BigH2_H2g;

	/* above sum was rate per unit vol since mult by H2 density, now div by H2* density to get rate s-1 */
	/* 0.25e-18 is wild guess of typical photodissociation cross section, from 
	 * >>refer	H2	dissoc	Allison, A.C. & Dalgarno, A. 1969, Atomic Data, 1, 91 
	 * this is based on an average of the highest v values they gave.  unfortunately, we want
	 * the highest J values - 
	 * final units are s-1*/ 
	hmi.H2_photodissoc_BigH2_H2s = hmi.H2_photodissoc_BigH2_H2s / SDIV(H2_BigH2_H2s) * H2_DISS_ALLISON_DALGARNO;
	hmi.H2_photodissoc_BigH2_H2g = hmi.H2_photodissoc_BigH2_H2g / SDIV(H2_BigH2_H2g) * H2_DISS_ALLISON_DALGARNO;
	hmi.H2_tripletdissoc_H2g = hmi.H2_tripletdissoc_H2g/SDIV(H2_BigH2_H2g);
	hmi.H2_tripletdissoc_H2s = hmi.H2_tripletdissoc_H2s/SDIV(H2_BigH2_H2s);
	hmi.Average_collH2g_dissoc = hmi.Average_collH2g_dissoc /SDIV(H2_BigH2_H2g);/* unit cm3s-1*/
	hmi.Average_collH2s_dissoc = hmi.Average_collH2s_dissoc /SDIV(H2_BigH2_H2s);/* unit cm3s-1*/
	hmi.H2s_LTE_bigH2 = hmi.H2s_LTE_bigH2/SDIV(H2_BigH2_H2s);
	hmi.H2g_LTE_bigH2 = hmi.H2g_LTE_bigH2/SDIV(H2_BigH2_H2g);

	/* >>chng 05 jul 09, GS*/ 
	/*  average Einstein value for H2* to H2g, GS*/
	double sumpop1 = 0.;
	double sumpopA1 = 0.;
	double sumpopcollH2O_deexcit = 0.;
	double sumpopcollH2p_deexcit = 0.;
	double sumpopcollH_deexcit = 0.;
	double popH2s = 0.;
	double sumpopcollH2O_excit = 0.;
	double sumpopcollH2p_excit = 0.;
	double sumpopcollH_excit = 0.;
	double popH2g = 0.;


	iElecLo = 0;
	for( iVibHi=0; iVibHi<=h2.nVib_hi[0]; ++iVibHi )
	{
		long nr1 = h2.nRot_hi[0][iVibHi];
		for( iRotHi=h2.Jlowest[0]; iRotHi<=nr1; ++iRotHi )
		{
			double ewnHi = energy_wn[0][iVibHi][iRotHi];
			for( iVibLo=0; iVibLo<=h2.nVib_hi[0]; ++iVibLo )
			{
				long nr2 = h2.nRot_hi[0][iVibLo];
				md3ci ewnLo = energy_wn.ptr(0,iVibLo,h2.Jlowest[0]);
				for( iRotLo=h2.Jlowest[0]; iRotLo<=nr2; ++iRotLo )
				{
					double ewnLo2 = energy_wn[0][iVibLo][iRotLo];
					/*if( ewnHi > ENERGY_H2_STAR && *ewnLo++ < ENERGY_H2_STAR )*/
					if( ewnHi > ENERGY_H2_STAR && ewnLo2 < ENERGY_H2_STAR )
					{
						/* >>chng 05 jul 10, GS*/ 
						/*  average collisional rate for H2* to H2g, GS*/
						if( H2_lgOrtho[0][iVibHi][iRotHi] == H2_lgOrtho[0][iVibLo][iRotLo] )
						{ 
							/* sums of populations */
							popH2s += H2_populations[0][iVibHi][iRotHi];
							popH2g += H2_populations[0][iVibLo][iRotLo];

							/* sums of deexcitation rates - H2* to H2g */
							sumpopcollH_deexcit += H2_populations[0][iVibHi][iRotHi]*H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][0];
							sumpopcollH2O_deexcit += H2_populations[0][iVibHi][iRotHi]*H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][2];
							sumpopcollH2p_deexcit += H2_populations[0][iVibHi][iRotHi]*H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][3];

							/* sums of excitation rates - H2g to H2* */
							sumpopcollH_excit += H2_populations[0][iVibLo][iRotLo]*H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][0]*H2_stat[0][iVibHi][iRotHi] / H2_stat[0][iVibLo][iRotLo] *
								H2_Boltzmann[0][iVibHi][iRotHi] /SDIV( H2_Boltzmann[0][iVibLo][iRotLo] );
							sumpopcollH2O_excit += H2_populations[0][iVibLo][iRotLo]*H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][2]*H2_stat[0][iVibHi][iRotHi] / H2_stat[0][iVibLo][iRotLo] *
								H2_Boltzmann[0][iVibHi][iRotHi] /SDIV( H2_Boltzmann[0][iVibLo][iRotLo] );
							sumpopcollH2p_excit += H2_populations[0][iVibLo][iRotLo]*H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][3]*H2_stat[0][iVibHi][iRotHi] / H2_stat[0][iVibLo][iRotLo] *
								H2_Boltzmann[0][iVibHi][iRotHi] /SDIV( H2_Boltzmann[0][iVibLo][iRotLo] );



							/* if( (abs((iRotHi-iRotLo) ))==2 || (iRotHi==iRotLo ) ) */
							if( lgH2_line_exists[0][iVibHi][iRotHi][0][iVibLo][iRotLo] )
							{
								sumpop1 += H2Lines[0][iVibHi][iRotHi][0][iVibLo][iRotLo].Hi->Pop;
								sumpopA1 += H2Lines[0][iVibHi][iRotHi][0][iVibLo][iRotLo].Hi->Pop*
									H2Lines[0][iVibHi][iRotHi][0][iVibLo][iRotLo].Emis->Aul;
							}
						}
					}
				}
			}
		}
	}
	hmi.Average_A = sumpopA1/SDIV(sumpop1);

	/* collisional excitation and deexcitation of H2g and H2s */
	hmi.Average_collH2_deexcit = (sumpopcollH2O_deexcit+sumpopcollH2p_deexcit)/SDIV(popH2s);
	hmi.Average_collH2_excit = (sumpopcollH2O_excit+sumpopcollH2p_excit)/SDIV(popH2g);
	hmi.Average_collH_excit = sumpopcollH_excit/SDIV(popH2g);
	hmi.Average_collH_deexcit = sumpopcollH_deexcit/SDIV(popH2s);
	/* populations are badly off during search phase */
	if( conv.lgSearch )
	{
		hmi.Average_collH_excit /=10.;
		hmi.Average_collH_deexcit/=10.;
	}

	/*fprintf(ioQQQ,
		"DEBUG Average_collH_excit sumpop = %.2e %.2e %.2e %.2e %.2e %.2e \n", 
		popH2g,popH2s,sumpopcollH_deexcit ,sumpopcollH_excit ,
		sumpopcollH_deexcit/SDIV(popH2s) ,sumpopcollH_excit/SDIV(popH2g));*/
	/*fprintf(ioQQQ,"sumpop = %le sumpopA = %le  Av= %le\n", 
	sumpop1,sumpopA1 , hmi.Average_A );*/

	if( mole.nH2_TRACE >= mole.nH2_trace_full|| (trace.lgTrace && trace.lgTr_H2_Mole) )
	{
		fprintf(ioQQQ,"  H2_LevelPops exit2 Sol dissoc %.2e (TH85 %.2e)",
			hmi.H2_Solomon_dissoc_rate_BigH2_H2g +
				hmi.H2_Solomon_dissoc_rate_BigH2_H2s , 
			hmi.H2_Solomon_dissoc_rate_TH85_H2g);

		/* Solomon process rate from X into the X continuum with units s-1
		 * rates are total rate, and rates from H2g and H2s */ 
		fprintf(ioQQQ," H2g Sol %.2e H2s Sol %.2e",
			hmi.H2_Solomon_dissoc_rate_used_H2g , 
			hmi.H2_Solomon_dissoc_rate_BigH2_H2s );

		/* photoexcitation from H2g to H2s */
		fprintf(ioQQQ," H2g->H2s %.2e (TH85 %.2e)",
			hmi.H2_H2g_to_H2s_rate_BigH2 , 
			hmi.H2_H2g_to_H2s_rate_TH85);

		/* add up H2s + hnu => 2H, continuum photodissociation,
		 * this is not the Solomon process, true continuum, units s-1 */
		fprintf(ioQQQ," H2 con diss %.2e (TH85 %.2e)\n",
			hmi.H2_photodissoc_BigH2_H2s , 
			hmi.H2_photodissoc_TH85);
	}
	else if( mole.nH2_TRACE )
	{
		fprintf(ioQQQ,"  H2_LevelPops exit1 %8.2f loops:%3li H2/H:%.3e Sol dis old %.3e new %.3e",
			fnzone ,
			loop_h2_pops ,
			hmi.H2_total / dense.gas_phase[ipHYDROGEN],
			old_solomon_rate,
			hmi.H2_Solomon_dissoc_rate_BigH2_H2g );
		fprintf(ioQQQ,"\n");
	}

	/* >>chng 03 sep 01, add this population - before had just used H2star from chem network */
	/* if big H2 molecule is turned on and used for this zone, use its
	 * value of H2* (pops of all states with v > 0 ) rather than simple network */

	/* update number of times we have been called */
	++nCallH2_this_iteration;

	/* this will say how many times the large H2 molecule has been called in this zone -
	 * if not called (due to low H2 abundance) then not need to update its line arrays */
	++h2.nCallH2_this_zone;

	/* >>chng 05 jun 21,
	 * during search phase we want to use full matrix - save number of levels so that
	 * we can restore it */
	nXLevelsMatrix = nXLevelsMatrix_save;

	/* >>chng 05 jan 19, check how many levels should be in the matrix if first call on
	 * new zone, and we have a solution */
	/* end loop setting very first LTE H2_populations */
	if( nCallH2_this_iteration && nzone != nzone_nlevel_set )
	{
		/* this is fraction of populations to include in matrix */
		const double FRAC = 0.99999;
		/* this loop is over increasing energy */
		double sum_pop = 0.;
		nEner = 0;
		iElec = 0;
		const bool PRT = false;
		if( PRT ) fprintf(ioQQQ,"DEBUG pops ");
		while( nEner < nLevels_per_elec[0] && sum_pop/hmi.H2_total < FRAC )
		{

			/* array of energy sorted indices within X */
			ip = H2_ipX_ener_sort[nEner];
			iVib = ipVib_H2_energy_sort[ip];
			iRot = ipRot_H2_energy_sort[ip];
			sum_pop += H2_old_populations[iElec][iVib][iRot];
			if( PRT ) fprintf(ioQQQ,"\t%.3e ", H2_old_populations[iElec][iVib][iRot]);
			++nEner;
		}
		if( PRT ) fprintf(ioQQQ,"\n");
		nzone_nlevel_set = nzone;
		/*fprintf(ioQQQ,"DEBUG zone %.2f old nmatrix %li proposed nmatrix %li sum_pop %.4e H2_total %.4e\n", 
			fnzone , nXLevelsMatrix ,nEner , sum_pop,hmi.H2_total);
		nXLevelsMatrix = nEner;*/
	}
	return;
}
/*lint -e802 possible bad pointer */

/*H2_cooling evaluate cooling and heating due to H2 molecule, called by 
 * H2_LevelPops in convergence loop when h2 heating is important, also 
 * called by CoolEvaluate to get final heating - argument is name of 
 * routine that called it */
#if defined(__ICC) && defined(__i386)
#pragma optimization_level 1
#endif
void H2_Cooling(
	/* string saying who called this routine, 
	 * "H2lup" call within H2 level populations solver
	 * "CoolEvaluate" call from main cooling routine */
	 const char *chRoutine)
{
	long int iElecHi , iElecLo , iVibHi , iVibLo , iRotHi , iRotLo;
	double heatone ,
		rate_dn_heat, 
		rate_up_cool;
	long int nColl,
		ipHi, ipLo;
	double Big1_heat , Big1_cool,
		H2_X_col_cool , H2_X_col_heat;
	long int ipVib_big_heat_hi,ipVib_big_heat_lo ,ipRot_big_heat_hi ,
		ipRot_big_heat_lo ,ipVib_big_cool_hi,ipVib_big_cool_lo ,
		ipRot_big_cool_hi ,	ipRot_big_cool_lo;
/*#	define DEBUG_DIS_DEAT*/
#	ifdef DEBUG_DIS_DEAT
	double  heatbig;
	long int iElecBig , iVibBig , iRotBig;
#	endif

	/* option to keep track of strongest single heating agent due to collisions
	 * within X */
	enum {DEBUG_H2_COLL_X_HEAT=false };

	DEBUG_ENTRY( "H2_Cooling()" );

	/* possible debug counters, counter itself, counter to turn on
	 * debug output, and counter for stopping */
	static long int nCount=-1;
	long int nCountDebug = 930,
		nCountStop = 940;
	++nCount;

	/* nCallH2_this_iteration is not incremented until after the level
	 * populations have converged the first time.  so for the first n calls
	 * this will return zero, a good idea since populations will be wildly
	 * incorrect during search for first valid pops */
	if( !h2.lgH2ON || !nCallH2_this_iteration )
	{
		hmi.HeatH2Dexc_BigH2 = 0.;
		hmi.HeatH2Dish_BigH2 = 0.;
		hmi.deriv_HeatH2Dexc_BigH2 = 0.;
		return;
	}

	hmi.HeatH2Dish_BigH2 = 0.;
	heatone = 0.;
#	ifdef DEBUG_DIS_DEAT
	heatbig = 0.;
	iElecBig = -1;
	iVibBig = -1;
	iRotBig = -1;
#	endif
	/* heating due to dissociation of electronic excited states */
	for( iElecHi=1; iElecHi<mole.n_h2_elec_states; ++iElecHi )
	{
		for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
		{
			for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
			{
					/* population, cm-3, in excited state */
				heatone = H2_populations[iElecHi][iVibHi][iRotHi] * 
					H2_dissprob[iElecHi][iVibHi][iRotHi] *
					H2_disske[iElecHi][iVibHi][iRotHi];
				hmi.HeatH2Dish_BigH2 += heatone;
#				ifdef DEBUG_DIS_DEAT
				if( heatone > heatbig )
				{
					heatbig = heatone;
					iElecBig = iElecHi;
					iVibBig = iVibHi;
					iRotBig = iRotHi;
				}
#				endif
			}
		}
	}
#	ifdef DEBUG_DIS_DEAT
	fprintf(ioQQQ,"DEBUG H2 dis heat\t%.2f\t%.3f\t%li\t\t%li\t%li\n",
		fnzone ,
		heatbig / SDIV( hmi.HeatH2Dish_BigH2 ) ,
		iElecBig ,
		iVibBig ,
		iRotBig );
#	endif
	/* dissociation heating HeatH2Dish_BigH2 was in eV - 
	 * convert to ergs */
	hmi.HeatH2Dish_BigH2 *= EN1EV;

	/*fprintf(ioQQQ,"DEBUG H2 heat/dissoc %.3e\n", 
		hmi.HeatH2Dish_BigH2/ SDIV(hmi.H2_rate_destroy*hmi.H2_total) );*/
	/* now work on collisional heating due to bound-bound
	 * collisional transitions within X */
	hmi.HeatH2Dexc_BigH2 = 0.;
	H2_X_col_cool = 0.;
	H2_X_col_heat = 0.;
	/* these are the colliders that will be considered as depopulating agents */
	/* the colliders are H, He, H2 ortho, H2 para, H+ */
	/* atomic hydrogen */
	long int nBug1 = 200 , nBug2 = 201; 
#	if 0
	/* fudge(-1) returns number of fudge factors entered on fudge command */
	if( fudge(-1) > 0 )
	{
		nBug1 = (long)fudge(0);
		nBug2 = (long)fudge(1);
	}
#	endif
	if( DEBUG_H2_COLL_X_HEAT && (nCount == nBug1 || nCount==nBug2) )
	{
		FILE *ioBAD=NULL;
		if( nCount==nBug1 )
		{
			ioBAD = fopen("firstpop.txt" , "w" );
		}
		else if( nCount==nBug2 )
		{
			ioBAD = fopen("secondpop.txt" , "w" );
		}
		for( ipHi=0; ipHi<nLevels_per_elec[0]; ++ipHi )
		{
			long int ip = H2_ipX_ener_sort[ipHi];
			iVibHi = ipVib_H2_energy_sort[ip];
			iRotHi = ipRot_H2_energy_sort[ip];
			fprintf(ioBAD , "%li\t%li\t%.2e\n",
				iVibHi , iRotHi , 
				H2_populations[0][iVibHi][iRotHi] );
		}
		fclose(ioBAD);
	}

	/* now make sum of all collisions within X itself */
	iElecHi = 0;
	iElecLo = 0;
	Big1_heat = 0.;
	Big1_cool = 0.;
	ipVib_big_heat_hi = -1;
	ipVib_big_heat_lo = -1;
	ipRot_big_heat_hi = -1;
	ipRot_big_heat_lo = -1;
	ipVib_big_cool_hi = -1;
	ipVib_big_cool_lo = -1;
	ipRot_big_cool_hi = -1;
	ipRot_big_cool_lo = -1;
	/* this will be derivative */
	hmi.deriv_HeatH2Dexc_BigH2 = 0.;
	for( ipHi=1; ipHi<nLevels_per_elec[iElecHi]; ++ipHi )
	{
		long int ip = H2_ipX_ener_sort[ipHi];
		iVibHi = ipVib_H2_energy_sort[ip];
		iRotHi = ipRot_H2_energy_sort[ip];
		if( iVibHi > VIB_COLLID )
			continue;

		realnum H2statHi = H2_stat[iElecHi][iVibHi][iRotHi];
		double H2boltzHi = H2_Boltzmann[iElecHi][iVibHi][iRotHi];
		double H2popHi = H2_populations[iElecHi][iVibHi][iRotHi];
		double ewnHi = energy_wn[iElecHi][iVibHi][iRotHi];

		for( ipLo=0; ipLo<ipHi; ++ipLo )
		{
			double coolone , oneline;
			ip = H2_ipX_ener_sort[ipLo];
			iVibLo = ipVib_H2_energy_sort[ip];
			iRotLo = ipRot_H2_energy_sort[ip];
			if( iVibLo > VIB_COLLID)
				continue;

			rate_dn_heat = 0.;

			/* this sum is total downward heating summed over all colliders */
			mr5ci H2cr = H2_CollRate.begin(iVibHi,iRotHi,iVibLo,iRotLo);
			for( nColl=0; nColl<N_X_COLLIDER; ++nColl )
				/* downward collision rate */
				rate_dn_heat += H2cr[nColl]*collider_density[nColl];

			/* now get upward collisional cooling by detailed balance */
			rate_up_cool = rate_dn_heat * H2_populations[iElecLo][iVibLo][iRotLo] *
				/* rest converts into upward collision rate */
				H2statHi / H2_stat[iElecLo][iVibLo][iRotLo] *
				H2boltzHi / SDIV( H2_Boltzmann[iElecLo][iVibLo][iRotLo] );

			rate_dn_heat *= H2popHi;

			/* net heating due to collisions within X - 
			 * positive if heating, negative is cooling
			 * this will usually be heating if X is photo pumped
			 * in printout and in save heating this is called "H2cX" */
			double conversion = (ewnHi - energy_wn[iElecLo][iVibLo][iRotLo]) * ERG1CM;
			heatone = rate_dn_heat * conversion;
			coolone = rate_up_cool * conversion;
			/* this is net heating, negative if cooling */
			oneline = heatone - coolone;
			hmi.HeatH2Dexc_BigH2 += oneline;

			/* keep track of heating and cooling separately */
			H2_X_col_cool += coolone;
			H2_X_col_heat += heatone;

			if( 0 && DEBUG_H2_COLL_X_HEAT && (nCount == 692 || nCount==693) )
			{
				static FILE *ioBAD=NULL;
				if(ipHi == 1 && ipLo == 0)
				{
					if( nCount==692 )
					{
						ioBAD = fopen("firstheat.txt" , "w" );
					}
					else if( nCount==693 )
					{
						ioBAD = fopen("secondheat.txt" , "w" );
					}
					fprintf(ioBAD,"DEBUG start \n");
				}
				fprintf(ioBAD,"DEBUG BAD DAY %li %li %li %li %.3e %.3e\n",
					iVibHi , iRotHi, iVibLo , iRotLo , heatone , coolone );/**/
				if( ipHi==nLevels_per_elec[iElecHi]-1 &&
					ipLo==ipHi-1 )
					fclose(ioBAD);
			}

			/* derivative wrt temperature - assume exp wrt ground - 
			 * this needs to be divided by square of temperature in wn - 
			 * done at end of loop */
			hmi.deriv_HeatH2Dexc_BigH2 +=  (realnum)(oneline * ewnHi);

			/* oneline is net heating, positive for heating, negative
			 * for cooling */
			if( DEBUG_H2_COLL_X_HEAT )
			{
				if( oneline < Big1_cool ) 
				{
					Big1_cool = oneline;
					ipVib_big_cool_hi = iVibHi;
					ipVib_big_cool_lo = iVibLo;
					ipRot_big_cool_hi = iRotHi;
					ipRot_big_cool_lo = iRotLo;
				}
				else if( oneline > Big1_heat )
				{
					Big1_heat = oneline;
					ipVib_big_heat_hi = iVibHi;
					ipVib_big_heat_lo = iVibLo;
					ipRot_big_heat_hi = iRotHi;
					ipRot_big_heat_lo = iRotLo;
				}
			}

			/* this would be a major logical error */
			ASSERT( 
				(rate_up_cool==0 && rate_dn_heat==0) || 
				(energy_wn[iElecHi][iVibHi][iRotHi] > energy_wn[iElecLo][iVibLo][iRotLo]) );
		}/* end loop over lower levels, all collisions within X */
	}/* end loop over upper levels, all collisions within X */

	/* this debug statement will identify the single strongest heating or
	 * cooling agent within X and give the lowest 7 level populations */
	if( DEBUG_H2_COLL_X_HEAT )
	{
		/* nCount counts number of calls through here, 
		 * nCountDebug is count to turn on debug prints, 
		 * nCountStop is count to abort code */
		if(nCount>nCountDebug )
		{
			fprintf(ioQQQ,
				"DEBUG H2_Cooling A %li %15s, Te %.3e net Heat(Xcol) %.2e "
				"heat %.2e cool %.2e H+/0 "
				"%.2e n(H2)%.3e Sol rat %.3e grn J1->0%.2e frac heat 1 line "
				"%.2e Hi(v,j)%li %li Lo(v,J)%li %li frac cool 1 line %.2e "
				"Hi(v,j)%li %li Lo(v,J)%li %li POP(J=1,13)",
				nCount , chRoutine,
				phycon.te, 
				hmi.HeatH2Dexc_BigH2 , 
				H2_X_col_cool ,
				H2_X_col_heat ,
				dense.xIonDense[ipHYDROGEN][1]/SDIV(dense.xIonDense[ipHYDROGEN][0]),
				hmi.H2_total,
				hmi.H2_Solomon_dissoc_rate_BigH2_H2g,
				hmi.rate_grain_h2_J1_to_J0 ,
				Big1_heat/hmi.HeatH2Dexc_BigH2 ,
				ipVib_big_heat_hi , ipRot_big_heat_hi , ipVib_big_heat_lo , ipRot_big_heat_lo ,
				Big1_cool/hmi.HeatH2Dexc_BigH2 ,
				ipVib_big_cool_hi , ipRot_big_cool_hi , ipVib_big_cool_lo , ipRot_big_cool_lo );

				for( iRotLo=0; iRotLo<14; ++iRotLo )
				{
					fprintf(ioQQQ,"\t%.2e" , 
						H2_populations[0][0][iRotLo]/hmi.H2_total );
				}
				fprintf(ioQQQ,"\t%li\n",nCount); 
				/* now give collision rates for strongest heat/cool level */
				fprintf(ioQQQ,"DEBUG H2_Cooling B heat Coll Rate (lrg col) dn,up" );
				double HeatNet = 0.;
				for( nColl=0; nColl<N_X_COLLIDER; ++nColl )
				{
					fprintf(ioQQQ,"\t%.2e" , 
						H2_CollRate[ipVib_big_heat_hi][ipRot_big_heat_hi][ipVib_big_heat_lo][ipRot_big_heat_lo][nColl]*
						collider_density[nColl] );
					fprintf(ioQQQ,"\t%.2e" , 
						/* downward collision rate */
						H2_CollRate[ipVib_big_heat_hi][ipRot_big_heat_hi][ipVib_big_heat_lo][ipRot_big_heat_lo][nColl]*
						collider_density[nColl]*
						/* rest converts into upward collision rate */
						H2_stat[0][ipVib_big_heat_hi][ipRot_big_heat_hi] / H2_stat[0][ipVib_big_heat_lo][ipRot_big_heat_lo] *
						H2_Boltzmann[0][ipVib_big_heat_hi][ipRot_big_heat_hi] /
						SDIV( H2_Boltzmann[0][ipVib_big_heat_lo][ipRot_big_heat_lo] ) );
					HeatNet += 
						H2_CollRate[ipVib_big_heat_hi][ipRot_big_heat_hi][ipVib_big_heat_lo][ipRot_big_heat_lo][nColl]*
						collider_density[nColl]*H2_populations[iElecLo][ipRot_big_heat_hi][ipVib_big_heat_lo]
					-
						H2_CollRate[ipVib_big_heat_hi][ipRot_big_heat_hi][ipVib_big_heat_lo][ipRot_big_heat_lo][nColl]*
						collider_density[nColl]*
						/* rest converts into upward collision rate */
						H2_stat[0][ipVib_big_heat_hi][ipRot_big_heat_hi] / H2_stat[0][ipVib_big_heat_lo][ipRot_big_heat_lo] *
						H2_Boltzmann[0][ipVib_big_heat_hi][ipRot_big_heat_hi] /
						SDIV( H2_Boltzmann[0][ipVib_big_heat_lo][ipRot_big_heat_lo] ) *
						H2_populations[iElecLo][ipVib_big_heat_lo][ipRot_big_heat_lo] ;
				}
				/* HeatNet includes the level populations, rates do not */
				fprintf( ioQQQ , " HeatNet %.2e",HeatNet);
				fprintf(ioQQQ,"\n"); 
				fprintf(ioQQQ,"DEBUG H2_Cooling C cool Coll Rate (lrg col) dn,up" );
				double CoolNet = 0.;
				for( nColl=0; nColl<N_X_COLLIDER; ++nColl )
				{
					fprintf(ioQQQ,"\t%.2e" , 
						H2_CollRate[ipVib_big_cool_hi][ipRot_big_cool_hi][ipVib_big_cool_lo][ipRot_big_cool_lo][nColl]*
						collider_density[nColl] );
					fprintf(ioQQQ,"\t%.2e" , 
						/* downward collision rate */
						H2_CollRate[ipVib_big_cool_hi][ipRot_big_cool_hi][ipVib_big_cool_lo][ipRot_big_cool_lo][nColl]*
						collider_density[nColl]*
						/* rest converts into upward collision rate */
						H2_stat[0][ipVib_big_cool_hi][ipRot_big_cool_hi] / H2_stat[0][ipVib_big_cool_lo][ipRot_big_cool_lo] *
						H2_Boltzmann[0][ipVib_big_cool_hi][ipRot_big_cool_hi] /
						SDIV( H2_Boltzmann[0][ipVib_big_cool_lo][ipRot_big_cool_lo] ) );
					CoolNet += 
						H2_CollRate[ipVib_big_cool_hi][ipRot_big_cool_hi][ipVib_big_cool_lo][ipRot_big_cool_lo][nColl]*
						collider_density[nColl]*H2_populations[iElecLo][ipRot_big_cool_hi][ipVib_big_cool_lo]
					-
						H2_CollRate[ipVib_big_cool_hi][ipRot_big_cool_hi][ipVib_big_cool_lo][ipRot_big_cool_lo][nColl]*
						collider_density[nColl]*
						/* rest converts into upward collision rate */
						H2_stat[0][ipVib_big_cool_hi][ipRot_big_cool_hi] / H2_stat[0][ipVib_big_cool_lo][ipRot_big_cool_lo] *
						H2_Boltzmann[0][ipVib_big_cool_hi][ipRot_big_cool_hi] /
						SDIV( H2_Boltzmann[0][ipVib_big_cool_lo][ipRot_big_cool_lo] ) *
						H2_populations[iElecLo][ipVib_big_cool_lo][ipRot_big_cool_lo] ;
				}
				/* CoolNet includes the level populations, rates do not */
				fprintf( ioQQQ , " CoolNet %.2e",CoolNet);
				fprintf(ioQQQ,"\n"); 
		}
	}

	/* this is inside h2 cooling, and is called extra times when H2 heating is important */
	if( PRT_POPS ) 
		fprintf(ioQQQ,
		"  DEBUG H2 heat fnzone\t%.2f\trenrom\t%.3e\tte\t%.4e\tdexc\t%.3e\theat/tot\t%.3e\n",
		fnzone , 
		H2_renorm_chemistry , 
		phycon.te , 
		hmi.HeatH2Dexc_BigH2,
		hmi.HeatH2Dexc_BigH2/thermal.ctot);

	/* this is derivative of collisional heating wrt temperature - needs 
	 * to be divided by square of temperature in wn */
	hmi.deriv_HeatH2Dexc_BigH2 /=  (realnum)POW2(phycon.te_wn);

	{
		enum {DEBUG_LOC=false };
		if( DEBUG_H2_COLL_X_HEAT && DEBUG_LOC && 
			(fabs(hmi.HeatH2Dexc_BigH2) > SMALLFLOAT) )
		{
			int iVib = 0;

			/*fprintf(ioQQQ," H2_cooling pops\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n",
				H2_populations[0][iVib][0]/hmi.H2_total,
				H2_populations[0][iVib][1]/hmi.H2_total,
				H2_populations[0][iVib][2]/hmi.H2_total,
				H2_populations[0][iVib][3]/hmi.H2_total,
				H2_populations[0][iVib][4]/hmi.H2_total,
				H2_populations[0][iVib][5]/hmi.H2_total);*/

			iElecHi = iElecLo = 0;
			iVibHi = iVibLo = 0;
			iRotHi = 7;
			iRotLo = 5;
			rate_dn_heat = rate_up_cool = 0.;
			/* this sum is total downward heating */
			for( nColl=0; nColl<N_X_COLLIDER; ++nColl )
			{
				rate_dn_heat +=
					H2_populations[iElecHi][iVibHi][iRotHi] * 
					H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][nColl]*
					collider_density[nColl];

				/* now get upward collisional cooling by detailed balance */
				rate_up_cool += 
					H2_populations[iElecLo][iVibLo][iRotLo] *
					/* downward collision rate */
					H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][nColl]*
					collider_density[nColl]*
					/* rest converts into upward collision rate */
					H2_stat[iElecHi][iVibHi][iRotHi] / H2_stat[iElecLo][iVibLo][iRotLo] *
					H2_Boltzmann[iElecHi][iVibHi][iRotHi] /
					SDIV( H2_Boltzmann[iElecLo][iVibLo][iRotLo] );
			}

			fprintf(ioQQQ,"DEBUG H2_cooling D pop %li ov %li\t%.3e\tdn up 31\t%.3e\t%.3e\n",
				iRotHi , iRotLo ,
				H2_populations[0][iVib][iRotHi]/H2_populations[0][iVib][iRotLo],
				rate_dn_heat,
				rate_up_cool);
			if( nCount>= nCountStop )
				cdEXIT(EXIT_FAILURE);
		}
	}
	{
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC  )
		{
			static long nzdone=-1 , nzincre;
			if( nzone!=nzdone )
			{
				nzdone = nzone;
				nzincre = -1;
			}
			++nzincre;
			fprintf(ioQQQ," H2 nz\t%.2f\tnzinc\t%li\tTe\t%.4e\tH2\t%.3e\tcXH\t%.2e\tdcXH/dt%.2e\tDish\t%.2e \n",
				fnzone, 
				nzincre,
				phycon.te,
				hmi.H2_total ,
				hmi.HeatH2Dexc_BigH2,
				hmi.deriv_HeatH2Dexc_BigH2 ,
				hmi.HeatH2Dish_BigH2);

		}
	}

#	if 0
	/* this can be noisy due to finite accuracy of solution, so take average with
	 * previous value */
	/*>>chng 04 mar 01, do not take average */
	if( 1 || nzone <1 || old_HeatH2Dexc==0. || nCallH2_this_iteration <2)
	{
		old_HeatH2Dexc = hmi.HeatH2Dexc_BigH2;
	}
	else
	{
		hmi.HeatH2Dexc_BigH2 = (hmi.HeatH2Dexc_BigH2+old_HeatH2Dexc)/2.f;
		old_HeatH2Dexc = hmi.HeatH2Dexc_BigH2;
	}
#	endif
	{
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC /*&& DEBUG_H2_COLL_X_HEAT*/ )
		{
			fprintf(ioQQQ,"DEBUG H2_cooling E %15s %c vib deex %li Te %.3e net heat %.3e cool %.3e heat %.3e\n", 
				chRoutine ,
				TorF(conv.lgSearch),
				nCount,
				phycon.te, hmi.HeatH2Dexc_BigH2,
				H2_X_col_cool ,
				H2_X_col_heat /*,
				H2_populations[0][0][7]/SDIV(H2_populations[0][0][5]) ,
				H2_populations[0][0][13]/SDIV(H2_populations[0][0][11])*/ );
			if( 0 && nCount > nCountStop )
			{
				cdEXIT( EXIT_FAILURE );
			}
		}
	}
	if( mole.nH2_TRACE >= mole.nH2_trace_full ) 
		fprintf(ioQQQ,
		" H2_Cooling Ctot\t%.4e\t HeatH2Dish_BigH2 \t%.4e\t HeatH2Dexc_BigH2 \t%.4e\n" ,
		thermal.ctot , 
		hmi.HeatH2Dish_BigH2 , 
		hmi.HeatH2Dexc_BigH2 );

	/* when we are very far from solution, during search phase, collisions within
	 * X can be overwhelmingly large heating and cooling terms, which nearly 
	 * cancel out.  Some dense cosmic ray heated clouds could not find correct
	 * initial solution due to noise introduced by large net heating which was
	 * the very noisy tiny difference between very large heating and cooling
	 * terms.  Do not include collisions with x as heat/cool during the
	 * initial search phase */
	if( conv.lgSearch )
		hmi.HeatH2Dexc_BigH2 = 0.;
	return;
}


/*cdH2_colden return column density in H2, negative -1 if cannot find state,
 * header is cdDrive */
double cdH2_colden( long iVib , long iRot )
{

	/*if iVib is negative, return
	 * total column density - iRot=0
	 * ortho column density - iRot 1
	 * para column density - iRot 2 
	 * else return column density in iVib, iRot */
	if( iVib < 0 )
	{
		if( iRot==0 )
		{
			/* return total H2 column density */
			return( h2.ortho_colden + h2.para_colden );
		}
		else if( iRot==1 )
		{
			/* return ortho H2 column density */
			return h2.ortho_colden;
		}
		else if( iRot==2 )
		{
			/* return para H2 column density */
			return h2.para_colden;
		}
		else
		{
			fprintf(ioQQQ," iRot must be 0 (total), 1 (ortho), or 2 (para), returning -1.\n");
			return -1.;
		}
	}
	else if( h2.lgH2ON )
	{
		/* this branch want state specific column density, which can only result from
		 * evaluation of big molecule */
		int iElec = 0;
		if( iRot <0 || iVib >h2.nVib_hi[iElec] || iRot > h2.nRot_hi[iElec][iVib])
		{
			fprintf(ioQQQ," iVib and iRot must lie within X, returning -2.\n");
			fprintf(ioQQQ," iVib must be <= %li and iRot must be <= %li.\n",
				h2.nVib_hi[iElec],h2.nRot_hi[iElec][iVib]);
			return -2.;
		}
		else
		{
			return H2_X_colden[iVib][iRot];
		}
	}
	/* error condition - no valid parameter */
	else
		return -1;
}

/*H2_Colden maintain H2 column densities within X */
void H2_Colden( const char *chLabel )
{
	long int iVib , iRot;

	/* >>chng 05 jan 26, pops now set to LTE for small abundance case, so do this */
	if( !h2.lgH2ON /*|| !h2.nCallH2_this_zone*/ )
		return;

	DEBUG_ENTRY( "H2_Colden()" );

	if( strcmp(chLabel,"ZERO") == 0 )
	{
		/* zero out formation rates and column densites */
		for( iVib = 0; iVib <= h2.nVib_hi[0]; ++iVib )
		{
			for( iRot=h2.Jlowest[0]; iRot<=h2.nRot_hi[0][iVib]; ++iRot )
			{
				/* space for the rotation quantum number */
				H2_X_colden[iVib][iRot] = 0.;
				H2_X_colden_LTE[iVib][iRot] = 0.;
			}
		}
	}

	else if( strcmp(chLabel,"ADD ") == 0 )
	{
		/*  add together column densities */
		for( iVib = 0; iVib <= h2.nVib_hi[0]; ++iVib )
		{
			for( iRot=h2.Jlowest[0]; iRot<=h2.nRot_hi[0][iVib]; ++iRot )
			{
				/* state specific H2 column density */
				H2_X_colden[iVib][iRot] += (realnum)(H2_populations[0][iVib][iRot]*radius.drad_x_fillfac);
				/* LTE state specific H2 column density - H2_populations_LTE is normed to unity
				 * so must be multiplied by total H2 density */
				H2_X_colden_LTE[iVib][iRot] += (realnum)(H2_populations_LTE[0][iVib][iRot]*
					hmi.H2_total*radius.drad_x_fillfac);
			}
		}
	}

	/* we will not print column densities so skip that - if not print then we have a problem */
	else if( strcmp(chLabel,"PRIN") != 0 )
	{
		fprintf( ioQQQ, " H2_Colden does not understand the label %s\n", 
		  chLabel );
		cdEXIT(EXIT_FAILURE);
	}

	return;
}

/*H2_DR choose next zone thickness based on H2 big molecule */
double H2_DR(void)
{
	return BIGFLOAT;
}

/*H2_RT_OTS - add H2 ots fields */
void H2_RT_OTS( void )
{

	long int iElecHi , iVibHi , iRotHi , iElecLo , iVibLo , iRotLo;

	/* do not compute if H2 not turned on, or not computed for these conditions */
	if( !h2.lgH2ON || !h2.nCallH2_this_zone )
		return;

	DEBUG_ENTRY( "H2_RT_OTS()" );

	/* loop over all possible lines and set H2_populations, and quantities that depend on escape prob, dest, etc */
	long int lim_elec_hi = 0;
	for( iElecHi=0; iElecHi<=lim_elec_hi; ++iElecHi )
	{
		for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
		{
			for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
			{
				double H2popHi = H2Lines[iElecHi][iVibHi][iRotHi][0][0][0].Hi->Pop;
				long int lim_elec_lo = 0;
				for( iElecLo=0; iElecLo<=lim_elec_lo; ++iElecLo )
				{
					/* want to include all vibration states in lower level if different electronic level,
					* but only lower vibration levels if same electronic level */
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
							/* >>chng 05 feb 07, use lgH2_line_exists */
							if( iElecHi==0 && lgH2_line_exists[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo] )
							{
								/* ots destruction rate */
								H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->ots = 
									H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->Aul * H2popHi *
									H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->Pdest;

								/* dump the ots rate into the stack - but only for ground electronic state*/
								RT_OTS_AddLine(
									H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->ots,
									H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].ipCont );
							}
						}
					}
				}
			}
		}
	}

	return;
}
/*lint +e802 possible bad pointer */

