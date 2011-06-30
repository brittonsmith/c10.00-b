/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* mole_H2_LTE sets Boltzmann factors and LTE unit population of large H2 molecular */
/* H2_init - called to initialize things from cdInit */
/*H2_init_coreload one time initialization */
/* H2_Zero zero out vars in the large H2 molecule, called from zero 
 * before any commands are parsed */
/* H2_zero_pops_too_low - zero out some H2 variables if we decide not to compute
 * the full sim, called by H2_LevelPops*/
/* H2_Solomon_rate find rates between H2s and H2g and other levels,
 * for eventual use in the chemistry */
/* H2_gs_rates evaluate rates between ground and star states of H2 for use in chemistry */
/* H2_He_coll_init initialize H2 - He collision data set
 * H2_He_coll interpolate on h2 - He collision data set to return rate at temp*/
#include "cddefines.h" 
#include "phycon.h" 
#include "mole.h" 
#include "hmi.h" 
#include "taulines.h" 
#include "h2.h" 
#include "h2_priv.h" 

/*H2_Solomon_rate find rates between H2s and H2g and other levels,
 * for eventual use in the chemistry */
void H2_Solomon_rate( void )
{

	long int iElecLo , iElecHi , iVibLo , iVibHi , iRotLo , iRotHi;

	DEBUG_ENTRY( "H2_Solomon_rate()" );

	/* iElecLo will always be X in this routine */
	iElecLo = 0;

	/* find rate (s-1) h2 dissociation into X continuum by Solomon process and
	 * assign to the TH85 g and s states 
	 * these will go back into the chemistry network */

	/* rates [s-1] for dissociation from s or g, into electronic excited states  
	 * followed by dissociation */
	hmi.H2_Solomon_dissoc_rate_BigH2_H2g = 0.;
	hmi.H2_Solomon_dissoc_rate_BigH2_H2s = 0.;

	/* these are used in a print statement - are they needed? */
	hmi.H2_Solomon_elec_decay_H2g = 0.;
	hmi.H2_Solomon_elec_decay_H2s = 0.;

	/* at this point we have already evaluated the sum of the radiative rates out
	 * of the electronic excited states - this is H2_rad_rate_out[electronic][vib][rot] 
	 * and this includes both decays into the continuum and bound states of X */

	/* sum over all electronic states, finding dissociation rate */
	for( iElecHi=1; iElecHi<mole.n_h2_elec_states; ++iElecHi )
	{
		for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
		{
			for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
			{
				double factor = (double)H2_dissprob[iElecHi][iVibHi][iRotHi]/
					H2_rad_rate_out[iElecHi][iVibHi][iRotHi];
				double H2popHi = H2_populations[iElecHi][iVibHi][iRotHi];

				/* loop over all levels within X to find 
				 * decay rates from H2g and H2s to continuum 
				 * distinction between H2g and H2s is determined 
				 * by ENERGY_H2_STAR */
				iElecLo = 0;

				for( iVibLo=0; iVibLo<=h2.nVib_hi[iElecLo]; ++iVibLo )
				{
					long nr = h2.nRot_hi[iElecLo][iVibLo];

					mb6ci lgH2le = lgH2_line_exists.ptr(iElecHi,iVibHi,iRotHi,iElecLo,iVibLo,h2.Jlowest[iElecLo]);
					mt6ci H2L = H2Lines.ptr(iElecHi,iVibHi,iRotHi,iElecLo,iVibLo,h2.Jlowest[iElecLo]);
					for( iRotLo=h2.Jlowest[iElecLo]; iRotLo<=nr; ++iRotLo )
					{
						if( *lgH2le++ )
						{
							/* this is the rate [cm-3 s-1] that mole goes from 
							 * lower level into electronic excited states then 
							 * into continuum */
							double rate_up_cont = 
								H2_populations[iElecLo][iVibLo][iRotLo]*
								H2L->Emis->pump*factor;

							/* rate electronic state decays into H2g */
							double elec_decay = 
								H2popHi*
								H2L->Emis->Aul*(H2L->Emis->Pesc+H2L->Emis->Pdest+H2L->Emis->Pelec_esc);

							if( energy_wn[0][iVibLo][iRotLo] > ENERGY_H2_STAR )
							{
								/* this is H2g up to excited then to continuum - 
								 * cm-3 s-1 at this point */
								hmi.H2_Solomon_dissoc_rate_BigH2_H2s += rate_up_cont;
								/* rate electronic state decays into H2g */
								hmi.H2_Solomon_elec_decay_H2s += elec_decay;
							}
							else
							{
								/* this is H2g up to excited then to continuum - 
								 * cm-3 s-1 at this point */
								hmi.H2_Solomon_dissoc_rate_BigH2_H2g += rate_up_cont;
								/* rate electronic state decays into H2g */
								hmi.H2_Solomon_elec_decay_H2g += elec_decay;
							}
						}
						++H2L;
					}
				}
			}/* end iRotHi */
		}/* end iVibHi */
	}/* end iElecHi */
	/* at this point units of hmi.H2_Solomon_elec_decay_H2g, H2s are cm-3 s-1 
	 * since H2_populations are included -
	 * div by pops to get actual dissocation rate, s-1 */
	if( hmi.H2_total > SMALLFLOAT )
	{
		hmi.H2_Solomon_elec_decay_H2g /= SDIV( H2_sum_excit_elec_den );
		hmi.H2_Solomon_elec_decay_H2s /= SDIV( H2_sum_excit_elec_den );

		/* will be used for H2s-> H + H */
		hmi.H2_Solomon_dissoc_rate_BigH2_H2s = hmi.H2_Solomon_dissoc_rate_BigH2_H2s / SDIV(H2_den_s);

		/* will be used for H2g-> H + H */
		hmi.H2_Solomon_dissoc_rate_BigH2_H2g = hmi.H2_Solomon_dissoc_rate_BigH2_H2g / SDIV(H2_den_g);

	}
	else
	{
		hmi.H2_Solomon_dissoc_rate_BigH2_H2s = 0; 
		hmi.H2_Solomon_dissoc_rate_BigH2_H2g = 0;	

	}
	/*fprintf(ioQQQ,"DEBUG H2 new %.2e %.2e %.2e %.2e \n",
		hmi.H2_Solomon_elec_decay_H2g ,
		hmi.H2_Solomon_elec_decay_H2s ,
		hmi.H2_Solomon_dissoc_rate_BigH2_H2s,
		hmi.H2_Solomon_dissoc_rate_BigH2_H2g );*/

	/* rate g goes to s */
	hmi.H2_H2g_to_H2s_rate_BigH2 = 0.;
	return;
}

/*H2_gs_rates evaluate rates between ground and star states of H2 for use in chemistry */
void H2_gs_rates( void )
{
	long int ipLoX , ip , iRotLoX , iVibLoX ,
		iElecHi , iVibHi , ipOther , iRotOther,
		iRotHi , iVibOther;

	DEBUG_ENTRY( "H2_gs_rates()" );

	/* rate g goes to s */
	hmi.H2_H2g_to_H2s_rate_BigH2 = 0.;

	/* loop over all levels in H2g */
	for( ipLoX=0; ipLoX < nEner_H2_ground; ++ipLoX )
	{
		ip = H2_ipX_ener_sort[ipLoX];
		iRotLoX = ipRot_H2_energy_sort[ip];
		iVibLoX = ipVib_H2_energy_sort[ip];
		/* now find all pumps up to electronic excited states */
		/* sum over all electronic states, finding dissociation rate */
		for( iElecHi=1; iElecHi<mole.n_h2_elec_states; ++iElecHi )
		{
			for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
			{
				for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
				{
					if( lgH2_line_exists[iElecHi][iVibHi][iRotHi][0][iVibLoX][iRotLoX] )
					{
						/* this is the rate {cm-3 s-1] that mole goes from 
						 * lower level into electronic excited states then 
						 * into continuum */
						double rate_up_cont = 
							H2_populations[0][iVibLoX][iRotLoX]*
							H2Lines[iElecHi][iVibHi][iRotHi][0][iVibLoX][iRotLoX].Emis->pump;

						double decay_star = H2_rad_rate_out[iElecHi][iVibHi][iRotHi] - H2_dissprob[iElecHi][iVibHi][iRotHi];
						/* loop over all other levels in H2g, subtracting 
						 * their rate - remainder is rate into star, this is 
						 * usually only a few levels */
						for( ipOther=0; ipOther < nEner_H2_ground; ++ipOther )
						{
							ip = H2_ipX_ener_sort[ipOther];
							iRotOther = ipRot_H2_energy_sort[ip];
							iVibOther = ipVib_H2_energy_sort[ip];
							if( lgH2_line_exists[iElecHi][iVibHi][iRotHi][0][iVibOther][iRotOther] )
							{
								decay_star -= 
									H2Lines[iElecHi][iVibHi][iRotHi][0][iVibOther][iRotOther].Emis->Aul*(
									H2Lines[iElecHi][iVibHi][iRotHi][0][iVibOther][iRotOther].Emis->Pesc+
									H2Lines[iElecHi][iVibHi][iRotHi][0][iVibOther][iRotOther].Emis->Pdest+
									H2Lines[iElecHi][iVibHi][iRotHi][0][iVibOther][iRotOther].Emis->Pelec_esc);
							}
						}
						/* MAX because may underflow to negative numbers is rates very large 
						 * this is fraction that returns to H2s */
						decay_star = MAX2(0., decay_star)/SDIV(H2_rad_rate_out[iElecHi][iVibHi][iRotHi]);
						hmi.H2_H2g_to_H2s_rate_BigH2 += rate_up_cont*decay_star;

					}/* end if line exists */
				}/* end loop rot electronic excited */
			}/* end loop vib electronic excited */
		}/* end loop electronic electronic excited */
	}

	/* at this point units are cm-3 s-1 - convert to rate s-1 */
	hmi.H2_H2g_to_H2s_rate_BigH2 /= SDIV( H2_den_g );
	return;
}

/* H2_zero_pops_too_low - zero out some H2 variables if we decide not to compute
 * the full sim, called by H2_LevelPops*/
void H2_zero_pops_too_low( void )
{

	long int iElec, iElecHi, iElecLo, iVib , iVibLo , iVibHi ,
		iRot , iRotHi , iRotLo;

	DEBUG_ENTRY( "H2_zero_pops_too_low()" );

	/* >>chng 05 jan 26, add this block to set populations to LTE value */
	for( iElec=0; iElec<mole.n_h2_elec_states; ++iElec )
	{
		for( iVib=0; iVib<=h2.nVib_hi[iElec]; ++iVib )
		{
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
	/* zero everything out - loop over all possible lines */
	/* >>chng 05 jan 26, set to LTE values, since we still need to accumulate Lyman line
		* optical depths to have correct self-shielding when large h2 does come on */
	for( iElecHi=0; iElecHi<mole.n_h2_elec_states; ++iElecHi )
	{
		pops_per_elec[iElecHi] = 0.;
		for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
		{
			pops_per_vib[iElecHi][iVibHi] = 0.;
			for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
			{
				long int lim_elec_lo = 0;
				// this will update Lo->Pop as well as Lo and Hi point to the same set of data structures
				// the loops below are OK since Lo->Pop is only accessed after Hi->Pop has been set.
				H2Lines[iElecHi][iVibHi][iRotHi][0][0][0].Hi->Pop = H2_populations[iElecHi][iVibHi][iRotHi];
				/* now the lower levels 
				 * NB - iElecLo the lower electronic level is only X 
				 * we don't consider excited electronic to excited electronic trans here, since we are only 
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
								/* population of lower level with correction for stim emission */
								H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->PopOpc = 
									H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Lo->Pop - 
									H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Hi->Pop*
									H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Lo->g / 
									H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Hi->g;

								/* following two heat exchange excitation, deexcitation */
								H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Coll.cool = 0.;
								H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Coll.heat = 0.;

								/* intensity of line */
								H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->xIntensity = 0.;

								H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->phots = 0.;
								H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->ots = 0.;
								H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->ColOvTot = 0.;
							}
						}
					}
				}
			}
		}
	}
	hmi.H2_photodissoc_BigH2_H2s = 0.;
	hmi.H2_photodissoc_BigH2_H2g = 0.;
	hmi.HeatH2Dish_BigH2 = 0.;
	hmi.HeatH2Dexc_BigH2 = 0.;
	hmi.deriv_HeatH2Dexc_BigH2 = 0.;
	hmi.H2_Solomon_dissoc_rate_BigH2_H2g = 0.;
	hmi.H2_Solomon_dissoc_rate_BigH2_H2s = 0.;
	hmi.H2_H2g_to_H2s_rate_BigH2 = 0.;
	return;
}

/*mole_H2_LTE sets Boltzmann factors and LTE unit population of large H2 molecular */
void mole_H2_LTE( void )
{
	/* used to recall the temperature used for last set of Boltzmann factors */
	static double TeUsedBoltz = -1.;
	double part_fun;
	long int iElec , iVib , iRot;

	DEBUG_ENTRY( "mole_H2_LTE()" );

	/* do we need to update the Boltzmann factors and unit LTE populations? */
	if( ! fp_equal( phycon.te, TeUsedBoltz ) )
	{
		part_fun = 0.;
		TeUsedBoltz = phycon.te;
		/* loop over all levels setting H2_Boltzmann and deriving partition function */
		for( iElec=0; iElec<mole.n_h2_elec_states; ++iElec )
		{
			for( iVib=0; iVib<=h2.nVib_hi[iElec]; ++iVib )
			{
				for( iRot=h2.Jlowest[iElec]; iRot<=h2.nRot_hi[iElec][iVib]; ++iRot )
				{
					H2_Boltzmann[iElec][iVib][iRot] = 
						/* energy is relative to lowest level in the molecule, v=0, J=0,
						 * so Boltzmann factor is relative to this level */
						sexp( energy_wn[iElec][iVib][iRot] / phycon.te_wn );
					/* sum the partition function - Boltzmann factor times statistical weight */
					part_fun += H2_Boltzmann[iElec][iVib][iRot] * H2_stat[iElec][iVib][iRot];
					ASSERT( part_fun > 0 );
				}
			}
		}
		/* have partition function, set H2_populations_LTE (populations for unit H2 density) */
		for( iElec=0; iElec<mole.n_h2_elec_states; ++iElec )
		{
			for( iVib=0; iVib<=h2.nVib_hi[iElec]; ++iVib )
			{
				for( iRot=h2.Jlowest[iElec]; iRot<=h2.nRot_hi[iElec][iVib]; ++iRot )
				{
					/* these are the H2 LTE populations for a unit H2 density -
					 * these populations will sum up to unity */
					H2_populations_LTE[iElec][iVib][iRot] = 
						H2_Boltzmann[iElec][iVib][iRot] * 
						H2_stat[iElec][iVib][iRot] / part_fun;
					/*if( iElec==0 && iVib < 2)
						fprintf(ioQQQ,"DEBUG LTE pop\t%i\t%i\t%e\n",
						iVib,iRot,H2_populations_LTE[iElec][iVib][iRot]*hmi.H2_total);*/
				}
			}
		}
		if( mole.nH2_TRACE >= mole.nH2_trace_full ) 
			fprintf(ioQQQ,
			"mole_H2_LTE set H2_Boltzmann factors, T=%.2f, partition function is %.2f\n",
			phycon.te,
			part_fun);
	}

	return;
}

/*H2_init_coreload one time initialization */
void H2_init_coreload( void )
{
	/* the order of the electronic states is
	 * X, B, C+, C-, B', D+, and D- */
	/* this will be the number of vibration levels within each electronic */
	/* number of vib states within electronic states from
	 * >>refer	H2	energies	Abgrall, */
	long int nVib_hi_init[N_H2_ELEC] = {14 , 37 , 13 , 13, 9, 2 , 2};

	/* this gives the first rotational state for each electronic state - J=0 does
	 * not exist when Lambda = 1 */
	long int Jlowest_init[N_H2_ELEC] = {0 , 0 , 1  , 1 , 0 , 1 , 1 };

	/* number of rotation levels within each electronic - vib */
	/*lint -e785 too few init for aggregate */
	long int nRot_hi_init[N_H2_ELEC][50]=
		/* ground, X */
		{ {31, 30, 28, 27, 25, 
			23, 22, 20, 18, 16, 
			14, 12, 10,  7,  3 } ,
		/* B */
		{25,25,25,25,25,25,25,25, 25,25,
		 25,25,25,25,25,25,25,25, 25,25,
		 25,25,25,25,25,25,25,25, 23,21,
		 19,17,15,15,11,9,7, 7},
		/* C plus */
		{ 25, 25, 25, 25, 24, 23, 21, 19, 17, 14, 12, 10, 6, 2 },
		/* C minus (the same) */
		{ 25, 25, 25, 25, 24, 23, 21, 19, 17, 15, 13, 10, 7, 2 },
		/* B primed */
		{19,17, 14, 12, 9, 8, 7, 7, 4, 1 },
		/* D plus */
		{13, 10, 5},
		/* D minus */
		{25 , 25 ,25 } }
		;
	/*lint +e785 too few init for aggregate */
	/* one time init */

	long int iElec;

	DEBUG_ENTRY( "H2_init_coreload()" );

	/* the order of the electronic states is
	 * X, B, C+, C-, B', D+, and D- */
	/* this will be the number of vibration levels within each electronic */
	/* number of vib states within electronic states from
	 * >>refer	H2	energies	Abgrall, */
	for( iElec=0; iElec<N_H2_ELEC; ++iElec )
	{
		int iVib;
		h2.nVib_hi[iElec] = nVib_hi_init[iElec];
		h2.Jlowest[iElec] = Jlowest_init[iElec];
		for( iVib=0; iVib<=h2.nVib_hi[iElec]; ++iVib )
		{
			h2.nRot_hi[iElec][iVib] = nRot_hi_init[iElec][iVib];
			/*fprintf(ioQQQ,"DEBUG h2 set %li %li %li \n", 
				iElec , 
				iVib , 
				h2.nRot_hi[iElec][iVib] );*/
		}
	}
	strcpy( chH2ColliderLabels[0] , "H0" );
	strcpy( chH2ColliderLabels[1] , "He" );
	strcpy( chH2ColliderLabels[2] , "H2 o" );
	strcpy( chH2ColliderLabels[3] , "H2 p" );
	strcpy( chH2ColliderLabels[4] , "H+" );

	return;
}

/*H2_init - called to initialize things from cdInit */
void H2_Init(void)
{

	DEBUG_ENTRY( "H2_Init()" );

	/* the number of electronic quantum states to include.
	 * To do both Lyman and Werner bands want nelec = 3,
	 * default is to do all bands included */
	mole.n_h2_elec_states = N_H2_ELEC;
	h2.nCallH2_this_zone = 0;
	return;
}


/*H2_Reset called by IterRestart to reset variables that are needed after an iteration */
void H2_Reset( void )
{

	DEBUG_ENTRY( "H2_Reset()" );

	if(mole.nH2_TRACE) 
		fprintf(ioQQQ,
		"\n***************H2_Reset called, resetting nCallH2_this_iteration, zone %.2f iteration %li\n", 
		fnzone,
		iteration );

	/* number of times large molecules evaluated in this iteration,
	 * is false if never evaluated, on next evaluation will start with LTE populations */
	nCallH2_this_iteration = 0;

	/* these remember the largest and smallest factors needed to
	 * renormalize the H2 chemistry */
	h2.renorm_max = 1.;
	h2.renorm_min = 1.;

	/* counters used by H2_itrzn to find number of calls of h2 per zone */
	nH2_pops  = 0;
	nH2_zone = 0;
	/* this is used to establish zone number for evaluation of number of levels in matrix */
	nzone_nlevel_set = 0;

	nzoneAsEval = -1;
	iterationAsEval = -1; 

	/* zero out array used to save emission line intensities */
	H2_SaveLine.zero();

	return;
}

/*H2_Zero zero out vars in the large H2 molecule, called from zero 
 * before any commands are parsed 
 * NB - this routine is called before space allocated - must not zero out
 * any allocated arrays */
void H2_Zero( void )
{

	DEBUG_ENTRY( "H2_Zero()" );

	/* this is the smallest ratio of H2/H where we will bother with the large H2 molecule
	 * this value was chosen so that large mole used at very start of TH85 standard PDR,
	 * NB - this appears in headinfo and must be updated there if changed here */
	/* >>chng 03 jun 02, from 1e-6 to 1e-8 - in orion veil large H2 turned on half way
	 * across, and Solomon process was very fast since all lines optically thin.  correct
	 * result has some shielding, so reset to lower value so that H2 comes on sooner. */
	mole.H2_to_H_limit = 1e-8;

	h2.lgH2ON = false;
	/* flag to force using LTE level populations */
	mole.lgH2_LTE = false;

	/* counters used by H2_itrzn to find number of calls of h2 per zone */
	nH2_pops  = 0;
	nH2_zone = 0;
	/* this is used to establish zone number for evaluation of number of levels in matrix */
	nzone_nlevel_set = 0;

	/* these remember the largest and smallest factors needed to
	 * renormalize the H2 chemistry */
	h2.renorm_max = 1.;
	h2.renorm_min = 1.;

	nCallH2_this_iteration = 0;
	h2.ortho_density = 0.;
	h2.para_density = 0.;

	hmi.H2_Solomon_dissoc_rate_BigH2_H2s = 0.;
	hmi.H2_Solomon_dissoc_rate_BigH2_H2g = 0.;

	hmi.H2_H2g_to_H2s_rate_BigH2 = 0.;
	hmi.H2_photodissoc_BigH2_H2s = 0.;
	hmi.H2_photodissoc_BigH2_H2g = 0.;

	hmi.HeatH2Dexc_BigH2 = 0.;

	/* say that H2 has never been computed */
	hmi.lgBigH2_evaluated = false;

	hmi.lgH2_Thermal_BigH2 = true;
	hmi.lgH2_Chemistry_BigH2 = true;

	if( !lgH2_READ_DATA )
	{
		/* the number of electronic levels in the H2 molecule,
		 * to just do the Lyman and Werner bands set to 3 -
		 * reset with atom h2 levels command,
		 * default is all levels with data */
		mole.n_h2_elec_states = N_H2_ELEC;
	}

	/* the number of levels used in the matrix solution
	 * of the level H2_populations - set with atom h2 matrix nlevel,
	 * >>chng 04 oct 05, make default 30 levels 
	 * >>chng 04 dec 23, make default 70 levels */
	nXLevelsMatrix = 70;

	/* this is used to establish zone number for evaluation of number of levels in matrix */
	nzone_nlevel_set = -1;

	nzoneAsEval = -1;
	iterationAsEval = -1; 
	return;
}
