/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef HMI_H_
#define HMI_H_

/**hmole determine populations of hydrogen molecules */
void hmole(void);

/**hmole_reactions - evaluates hydrogen chemistry reactions */
void hmole_reactions(void);

/** hmole_init - initialize some hmole vars */
void hmole_init(void);

/**hmirat computes radiative association rate for H- 
\param te
*/
double hmirat(double te);

/** Take one Newton step of the chemical network 
\param *nFixup
\param *error
*/
void hmole_step(int *nFixup, double *error);

/** hmi.h - parameters dealing with hydrogen molecules */
EXTERN struct t_hmi {

	/** densities (cm-3) of H2, H2+, H-, H3+ */
	realnum Hmolec[N_H_MOLEC];
	int nProton[N_H_MOLEC];
	int nElectron[N_H_MOLEC];

	/** average energy level of H2g and H2s */
	double H2_BigH2_H2g_av;
	double H2_BigH2_H2s_av;

	/** ratio of H2g and H2s from the chemical network and the big molecule model*/
	double H2_chem_BigH2_H2g;
	double H2_chem_BigH2_H2s;

	/** labels for the H-bearing molecules */
	char chLab[N_H_MOLEC][5];

	/** the associative detachment H- + H => H2 + e rate coef */
	double assoc_detach;

	/**rate for process H2+ + H => H2 + H+ */
	double bh2h2p;

	/** rate H2 goes from all X into either J=1 (ortho) 
	 * or (J=0) para on grain surfaces - units s-1*/
	double rate_grain_h2_J1_to_J0,
	  rate_grain_h2_op_conserve;

	/** the total H2 density [cm-3], NOT 2*n(H2), the sum of H2 and H2* */
	realnum H2_total;
	realnum H2_total_BigH2;
	realnum H2s_BigH2;
	realnum H2g_BigH2;
	

	/** rate ground of H2 is destroyed */
	double H2_rate_destroy;

	/** Average Einstein A for H2s to H2g transition*/
	double Average_A;
	double h2s_sp_decay;

	/** Average noreactive collisional rate for H2s to H2g transition*/
	double Average_collH2_deexcit;
	double Average_collH_deexcit;
	double Average_collH2_excit;
	double Average_collH_excit;

	/**Average collisional dissociation by H2g and H2s*/
	double Average_collH2s_dissoc;
	double Average_collH2g_dissoc;

	/** hminus heating, free bound */
	double hmihet, 
	  hmitot, 
	  hmicol;

	/** mean cross section (cm-2) for H2 Lyman absorption */
	realnum H2Opacity;

	/** these are departure coef for H-, H2, H2+, and HeH,
	 * defined in hmole */
	double hmidep, 
	  h2dep, 
	  h2pdep, 
	  h3pdep;

	/** heating due to photo dissoc of H2+ */
	double h2plus_heat;

	/** H- photo dissoc rate */
	double HMinus_photo_rate;

	realnum 
	  /** the largest fraction of total heat anywhere in model */
	  HeatH2DexcMax,
	  /** the largest fraction of total cooling anywhere in model */
	  CoolH2DexcMax,
	  h2dfrc, 
	  h2dtot,
	  /** fraqction of cooling carried by H2 lines */
	  h2line_cool_frac;

	double HMinus_induc_rec_cooling,
	  HMinus_induc_rec_rate,
	  HMinus_photo_heat,
	  hminus_rad_attach;

	long int iheh1, 
	  iheh2;

	/** rate hi dest H_2 */
	realnum rh2dis;

	/** Ha creation due to H- charge transfer */
	realnum HalphaHmin;

	/** max H2/H ratio that occured in the calculation, set in hmole */
	realnum BiggestH2;

	/** UV flux relative to Habing value, used for some simple molecular photodissociation rates,
	 * as defined by Tielens & Hollenbach 1985 */
	realnum UV_Cont_rel2_Habing_TH85_face,
	  UV_Cont_rel2_Habing_TH85_depth,
	  /** the special version of g0 with adjustable bounds */
	  UV_Cont_rel2_Habing_spec_depth;

	/** UV flux relative to Habing value, used for some simple molecular photodissociation rates,
	 * as defined by Draine & Bertoldi 1996 -0 we try to do this the way they describe,
	 * since they say that this will agree with their large H2 molecule, first
	 * define field at the illuminated face, then get value at depth using their
	 * form of the extinction and shielding, rather than our exact calculation */
	realnum UV_Cont_rel2_Draine_DB96_face , 
		UV_Cont_rel2_Draine_DB96_depth;

	/** the Solomon process excitation, H2g -> H2*, rate from Tielens & Hollenbach 85 */
	double H2_H2g_to_H2s_rate_TH85;

	/** the Solomon process excitation, H2g -> H2*, rate from Burton et al. 1990 */
	double H2_H2g_to_H2s_rate_BHT90;
	
	/** the Solomon process excitation, H2g -> H2*, rate for the Bertodi & Draine model */
	double H2_H2g_to_H2s_rate_BD96;
	
	/** the Solomon process excitation, H2g -> H2*, rate for Elwert et al. model in prep.*/
	double H2_H2g_to_H2s_rate_ELWERT;
	
	/** the Solomon process excitation, H2g -> H2*, rate (s-1) from large molecules */
	double H2_H2g_to_H2s_rate_BigH2;

	/** the Solomon process excitation, H2g -> H2*, - actually used */
	double H2_H2g_to_H2s_rate_used;

	/** the Solomon process dissociate rate from Tielens & Hollenbach 85 */
	double H2_Solomon_dissoc_rate_used_H2g;
	double H2_Solomon_dissoc_rate_BigH2_H2g;
	double H2_Solomon_dissoc_rate_TH85_H2g;
	double H2_Solomon_dissoc_rate_BHT90_H2g;
	double H2_Solomon_dissoc_rate_BD96_H2g;
	double H2_Solomon_dissoc_rate_ELWERT_H2g;

	double H2_Solomon_dissoc_rate_used_H2s;
	double H2_Solomon_dissoc_rate_BigH2_H2s;
	double H2_Solomon_dissoc_rate_TH85_H2s;
	double H2_Solomon_dissoc_rate_BHT90_H2s;
	double H2_Solomon_dissoc_rate_BD96_H2s;
	double H2_Solomon_dissoc_rate_ELWERT_H2s;

	/** the Solomon process rate H2 dissociates into X continuum - actually used */
	/**double H2_Solomon_dissoc_rate_used;*/
	/** H2 + hnu => 2H from TH85 */
	/** H2 + hnu => 2H actually used */
	double H2_photodissoc_used_H2g;
	double H2_photodissoc_used_H2s;
	double H2_photodissoc_BigH2_H2s;
	double H2_photodissoc_BigH2_H2g;
	double H2_photodissoc_ELWERT_H2g;
	double H2_photodissoc_ELWERT_H2s;
	double H2_photodissoc_TH85;
	double H2_photodissoc_BHT90;

	/** these are decay rates from electronic levels into H2g and H2s */
	double H2_Solomon_elec_decay_H2g ,
	  H2_Solomon_elec_decay_H2s;
	
	/** H2 + hnu(continuum) => 2H from big molecule */
	/** H2 dissociation to triplet state*/
	double H2_tripletdissoc_H2s,
	  H2_tripletdissoc_H2g;

	/** says whether big H2 has ever been evaluated in this run - if it has
	 * not been then use TH85 physics for mole balance and cooling */
	bool lgBigH2_evaluated;

	/** continuum array index for H minus threshold  */
	long int iphmin;

	/** largest local fraction heating due to dissoc of H2+ */
	realnum h2pmax;

	/** binding energy for change in H2 population while on grain surface,
	 * set with "set h2 Tad " command */
	realnum Tad;

	double 

	  /** HeatH2Dish_used is heating due to H2 dissociation actually used*/
	  HeatH2Dish_used, 
	  HeatH2Dish_BigH2, 
	  HeatH2Dish_TH85, 
	  HeatH2Dish_BD96 ,
	  HeatH2Dish_BHT90, 
	  HeatH2Dish_ELWERT ,

	  /** HeatH2Dexc_used is heating due to collisional deexcitation of vib-excited 
	   * H2 actually used */
	  HeatH2Dexc_used,
	  HeatH2Dexc_BigH2,
	  HeatH2Dexc_TH85,
	  HeatH2Dexc_BD96,
	  HeatH2Dexc_BHT90,
	  HeatH2Dexc_ELWERT;

	/** these are derivative wrt temp for collisional processes within X */
	realnum 
		deriv_HeatH2Dexc_used,
		deriv_HeatH2Dexc_BigH2 ,
		deriv_HeatH2Dexc_TH85 ,
		deriv_HeatH2Dexc_BD96 ,
		deriv_HeatH2Dexc_BHT90 ,
		deriv_HeatH2Dexc_ELWERT;

	/** these are the H- and grain formation rates, added above and below a
	 * certain energy (2.6 eV) for production of H2 or H2* in small network */
	double H2_forms_grains ,
		H2_forms_hminus,
		H2star_forms_grains,
		H2star_forms_hminus;

	/** say how to do thermal solution, if true (default) use results of large molecule,
	 * if false use TH85 approximations */
	bool lgH2_Thermal_BigH2,
	/** say how to do chemistry (formation and destruction), 
	 * if true (default) use results of large molecule,
	 * if false use TH85 approximations */
	  lgH2_Chemistry_BigH2;

	/** option to turn off H molecules */
	bool lgNoH2Mole;

	/** the set h2 small model command tells code says which of the small model H2
	 * to use.  Default is Elwert */
	char chH2_small_model_type;

	/** method used for grain formation pumping */
	char chGrainFormPump;

	/** the set h2 jura command tells code which treatment of H2 formation to use */
	char chJura;

	/** this is a scale factor to multiply the Jura rate, default is unity, changed
	 * with the set jura scale command */
	realnum ScaleJura;

	/** H2 formation rate as set with set h2 rate command units S^-1, actual depl */
	double rate_h2_form_grains_set;  

	/** this is set to zero, but to positive number with atom h2 fraction command
	 * this sets the H2 density by multiplying the hydrogen density to become the H2 density */
	double H2_frac_abund_set;

	/** scale the H2 formation rate.  */
	double H2_formation_scale;

	/** rate coefficient (cm3 s-1) for reaction He+ + H2 -> He + H+ + H,
	 * needed for both H2 and He solvers */
	/** chng 04 jun 30 -- He+ + H2 => He + H2+, also important for He solver */
	realnum rheph2hpheh,
	  heph2heh2p;	

	/** rate coefficient (cm3 s-1) for H- + A+ -> H + A */
	realnum hmin_ct_firstions;

	/** Boltzmann factor for hmi */
	double exphmi,
	/** related to the LTE populations of H-, H2, and H2+
	 * each is a constant with temperature dependence, and
	 * needs to be multiplied by the densities of the separated
	 * components to become the LTE density.  
	 * following is n(H-) / [ n(e) n(H) ], units cm3 */
	  rel_pop_LTE_Hmin,
	/** related to the LTE population of H2s, following is  
	 * n(H2s) / [n(H) n(H) ], units cm3 */
	  rel_pop_LTE_H2s;
	/** LTE population for H2+, following is
	 * n(H2+) / [n(H) n(p) ], units cm3 */
	double rel_pop_LTE_H2p,
	/** related to the LTE population of H2 in ground, following is  
	 * n(H2) / [n(H) n(H) ], units cm3 */
	  rel_pop_LTE_H2g,
	  /** related to population of H3+ */
	  rel_pop_LTE_H3p,
	  /** LTE pops of g and s used for H- back reactions */
	  H2g_LTE_bigH2,
	  H2s_LTE_bigH2;

	/** hack to kill effects of H2* in chemistry network "set leiden hack h2* off */
	bool lgLeiden_Keep_ipMH2s;
	bool lgLeidenCRHack;

	/** some H2 dest and creation rates set in hmole_step and output in save h2 creation or destruction */
	/** >>chng 05 mar 18, TE, add more rates to save in H2 destruction file*/
	  double assoc_detach_backwards_grnd, 
		  assoc_detach_backwards_exct, 
		  bh2h22hh2,
		  h3phmh2hh,
		  h3phm2h2,
		  h32h2,
		  eh3_h2h,
		  h3ph2hp,
		  h2sh,
		  CR_reac_H2g,
		  CR_reac_H2s,
		  h2phmh2h,
		  hehph2h3phe,
		  h2ph3p,
		  h2sh2g,
		  h2h22hh2, 
		  h2sh2sh2g2h, 
		  h2sh2sh2s2h, 
		  H2_photoionize_rate ,
		  H2_photo_heat_soft ,
		  H2_photo_heat_hard ,
		  rh2h2p ,
		  eh2hh ,
		  h2ge2h,
		  h2se2h,
		  h2hph3p,
		  bh2dis,
		  radasc, 
		  h3ph2p, 
		  h3petc,
		  /**< total H2 creation rate, cm-3 s-1 */
		  H2_rate_create;

	}	hmi;

/** labels for various H molecules */
enum {
	ipMH,    /**< 0 H0 */
	ipMHp,   /**< 1 H+ */
	ipMHm,   /**< 2 H- */
	ipMH2g,  /**< 3 H2g -ground - hmi.H2_total is total*/
	ipMH2p,  /**< 4 H2+ */
	ipMH3p,  /**< 5 H3+ */
	ipMH2s,  /**< 6 H2* -exct - s == "star" - hmi.H2_total is total */
	ipMHeHp  /**< 7 HeH+ */
};


#endif /* HMI_H_ */
