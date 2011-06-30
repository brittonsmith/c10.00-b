/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef MOLE_H_
#define MOLE_H_

/* mole.h */

/**CO_drive main driver for heavy molecular equilibrium routines */
extern void CO_drive(void);

/**CO_zero allocate + initialize workspace */
extern void CO_zero(void);

/**CO_create_react build reaction structures */
extern void CO_create_react(void);

/** called from cdInit to initialized co routines */
extern void CO_Init(void);
/**CO_update_rks update rate coefficients, only temp part */
extern void CO_update_rks( void );

extern void CO_update_species_cache(void);

extern double CO_sink_rate(const char chSpecies[]);

extern double CO_source_rate(const char chSpecies[]);

extern double CO_dissoc_rate(const char chSpecies[]);

struct COmole_rate_s *CO_findrate_s(const char buf[]);

extern double CO_findrk(const char buf[]);

extern double CO_findrate(const char buf[]);

/** \verbatim >>chng 03 feb 09, rm ipH3P_hev, since not used, and decrement NUM_HEAVY_MOLEC to 17 
 >>chng 03 aug 04, rm ipCTWO and ipC2P from hevmol since not included in balance,
 and always finds zero column density, so NUM_HEAVY_MOLEC from 17 to 15 
 >>chng 03 aug 05, rm ch2 and ch3, so n from 15 to 13 
 >>chng 03 nov 14  add Si chemistry & CH3+, so that now every
     reaction that is in the TH85 chemical network is also included
     in Cloudy.  Additionally, there are also reactions taken from other
     papers (mostly Hollenbach and McKee...see co.c).  In all 20 molecular
     species are calculated, along with the atomic and first ionization 
	 stages of C, O, and Si
 >>chng 04 May 13, Nick Abel.  Add CH3, CH4, CH4+, and CH5+ to network in order 
	to get the same chemical abundances vs. depth as other PDR codes in the Leiden
	meeting.  With changes we now can predict molecular abundances for 24 C, O, 
	and Si bearing molecules. 

 >>chng 04 jul 13, Nick Abel.  Add nitrogen and sulphur bearing molecules
    to the chemical network.  First added to generate a chemical model for
    eta carinae, but is applicable to all molecular clouds 

 >>chng 05 mar 11, Nick Abel.  Add C2 and C2+ to chemistry, reactions 
    involving these species affects the abundance of C

 >>chng 05 mar 23, Nick Abel.  Add Chlorine to chemistry 
 \endverbatim */

/** this includes the atomic and first ionized species
   of each element that can combine to form molecules.  
   This is the number of molecules, ions, and atoms that the co network uses
   This is used in comole
   to improve the calculation, as deep in molecular regions reactions with molecules
   can be important to the ionization balance */

extern struct molecule null_mole;

EXTERN struct t_co {

	/** CODissHeat is CO Photodissociation heating */
	realnum CODissHeat, 
	  /**< largest fraction of total heating */
	  codfrc, 
	  /**< total heating integrated over cloud */
	  codtot;

	/** flag to turn off CO molecules, set with no molecules command */
	bool lgNoCOMole;

	/** This lets the code know if the CO network is turned on */
	bool lgCODoCalc;

	/** molecules of the heavy elements.  order must exactly match that in hevcolumn
	 * number of molecules in heavy element network */

	/** the density (cm-3) of each species */
	realnum 
	  /** these are source and sink terms for the hydrogen molecular
	     network that appear in co.c */
	  hydro_sink[8],
	  hydro_source[8]; 

	/** Molecular mass and reduced mass between any two reactions.  
	 * This includes the pure hydrogen 
	 * molecules and He+, and so the dimensions of the array are increased 
	 * by the number of molecules in the hydrogen chemistry plus 1 */

	double hmole_mass[N_H_MOLEC+1];

	/** C12/C13 isotope ratio, sets the ratio of C12O16 to C13O16 and 
	* C13/C12 1909, initialized in zero.c  */
	realnum C12_C13_isotope_ratio;

	/** flag set true if H2O destruction rate went to zero */
	bool lgH2Ozer;

	/** set rates to that in UMIST */
	bool lgUMISTrates;

	bool lgFederman;

	/** option to use effective temperature as defined in
	 * >> refer Federman, S. R. & Zsargo, J. 2003, ApJ, 589, 319
	 * By default, this is false - changed with set chemistry command */
	bool lgNonEquilChem;

	/** option to set proton elimination rates to zero
	 * >>refer	CO	chemistry	Huntress, W. T., 1977, ApJS, 33, 495
	 * By default, this is false - changed with set chemistry command */
	bool lgProtElim;

	/** option to not include neutrals in the non-equilibrium scheme
	 * >> refer Federman, S. R. & Zsargo, J. 2003, ApJ, 589, 319
	 * By default, this is false - changed with set chemistry command */
	bool lgNeutrals;

	double h2lim;

	long co_nzone , iteration_co;

	/** total electron density in molecules */
	double comole_eden;

#	if 0
	/** this is element part of the CO chemistry network? */
	bool lgElem_in_CO_chem[LIMELM];
#	endif


	realnum nitro_dissoc_rate;

	/**  These reactions and included in hmole_step as formation and destruction processes */

	double H_CH_C_H_H;
	double H_OH_O_H_H	;
	double H_H2O_OH_H_H;
	double H_COP_CO_HP;	
	double H_CH_C_H2;
	double H_CHP_CP_H2;
	double H_CH2_CH_H2;
	double H_CH3P_CH2P_H2;
	double H_OH_O_H2;
	double H_H2O_OH_H2;
	double Hminus_HCOP_CO_H2;
	double Hminus_H3OP_H2O_H2;
	double Hminus_H3OP_OH_H2_H;
	double HP_CH_CHP_H;
	double HP_CH2_CH2P_H;
	double HP_H2O_H2OP_H;
	double HP_O2_O2P_H;
	double HP_OH_OHP_H;
	double HP_SiO_SiOP_H;
	double HP_CH2_CHP_H2;
	double HP_SiH_SiP_H2;
	double H2_CHP_CH2P_H;
	double H2_CH2P_CH3P_H;
	double H2_OHP_H2OP_H;
	double H2_H2OP_H3OP_H;
	double H2_COP_HCOP_H;
	double H2_OP_OHP_H;
	double H2_SiOP_SiOHP_H;
	double H2_C_CH_H;
	double H2_CP_CHP_H;
	double H2_CH_CH2_H;
	double H2_OH_H2O_H;
	double H2_O_OH_H;
	double H2_CH_C_H2_H;
	double H2_OH_O_H2_H;
	double H2_H2O_OH_H2_H;
	double H2_O2_O_O_H2;
	double H2_O2_OH_OH;
	double H2s_CH_C_H2_H;
	double H2s_OH_O_H2_H;
	double H2s_H2O_OH_H2_H;
	double H2s_O2_O_O_H2;
	double H2P_C_CHP_H;
	double H2P_CH_CH2P_H;	
	double H2P_CH2_CH3P_H;
	double H2P_OH_H2OP_H;
	double H2P_H2O_H3OP_H;
	double H2P_CO_HCOP_H;
	double H2P_O_OHP_H;	
	double H2P_CH_CHP_H2;
	double H2P_CH2_CH2P_H2;
	double H2P_CO_COP_H2;
	double H2P_H2O_H2OP_H2;
	double H2P_O2_O2P_H2;
	double H2P_OH_OHP_H2;
	double H3P_C_CHP_H2;
	double H3P_CH_CH2P_H2;
	double H3P_CH2_CH3P_H2;
	double H3P_OH_H2OP_H2;
	double H3P_H2O_H3OP_H2;
	double H3P_CO_HCOP_H2;
	double H3P_O_OHP_H2;
	double H3P_SiH_SiH2P_H2;
	double H3P_SiO_SiOHP_H2;
	double H2s_CH_CH2_H;
	double H2s_O_OH_H;
	double H2s_OH_H2O_H;
	double H2s_C_CH_H;
	double H2s_CP_CHP_H;
	double H_CH3_CH2_H2;
	double H_CH4P_CH3P_H2;
	double H_CH5P_CH4P_H2;
	double H2_CH2_CH3_H;
	double H2_CH3_CH4_H;
	double H2_CH4P_CH5P_H;
	double H2s_CH2_CH3_H;
	double H2s_CH3_CH4_H;
	double H2P_CH4_CH3P_H2;
	double H2P_CH4_CH4P_H2;
	double H2P_CH4_CH5P_H;
	double H3P_CH3_CH4P_H2;
	double H3P_CH4_CH5P_H2;
	double HP_CH3_CH3P_H;
	double HP_CH4_CH3P_H2;
	double HP_CH4_CH4P_H;
	double HP_HNC_HCN_HP;
	double H_HNC_HCN_H;
	double H2_HCNP_HCNHP_H;
	double H3P_HCN_HCNHP_H2;
	double H2s_OP_OHP_H;

	double C_H3OP_HCOP_H2_1,
	C_OH_CO_H_1,
	CP_OH_CO_HP_1,
	CP_H2O_HCOP_H_1,
	CP_OH_COP_H_1,
	O_CH_CO_H_1,
	O_CHP_COP_H_1,
	O_CH2_CO_H_H_1,
	O_CH2_CO_H2_1,
	O_CH2P_HCOP_H_1,
	O_CH3P_HCOP_H2_1,
	O_H2OP_O2P_H2_1,
	O_OH_O2_H_1,
	O_OHP_O2P_H_1,
	O_SiH_SiO_H_1,
	O_SiH2P_SiOHP_H_1,
	OP_CH_COP_H_1,
	OP_OH_O2P_H_1,
	Si_OH_SiO_H_1,
	SiP_H2O_SiOHP_H_1,
	SiP_OH_SiOP_H_1,
	CHP_H2O_HCOP_H2_1,
	CHP_OH_COP_H2_1,
	H_C_CH_nu,
	H_CP_CHP_nu,
	H_OH_H2O_nu,
	Hminus_CH_CH2_e,
	Hminus_C_CH_e,
	Hminus_OH_H2O_e,
	Hminus_O_OH_e,
	H2_C_CH2_nu,
	H2_CP_CH2P_nu,
	H2_SiP_SiH2P_nu,
	HeP_CH_CP_He_H,
	HeP_CH2_CHP_He_H,
	HeP_OH_OP_He_H,
	HeP_H2O_OHP_He_H,
	HeP_SiH_SiP_He_H,
	HeP_H2O_OH_He_HP,
	HeP_CH2_CP_He_H2,
	crnu_CH_C_H,
	crnu_CHP_CP_H,
	crnu_H2O_OH_H,
	crnu_OH_O_H,
	crnu_SiH_Si_H,
	nu_CH_C_H,
	nu_CHP_CP_H,
	nu_CH2_CH_H,
	nu_CH2P_CHP_H,
	nu_CH3P_CH2P_H,
	nu_CH3P_CHP_H2,
	nu_H2O_OH_H,
	nu_OH_O_H,
	nu_OHP_O_HP,
	nu_SiH_Si_H,
	e_CHP_C_H,
	e_CH2P_CH_H,
	e_CH2P_C_H_H,
	e_CH2P_C_H2,
	e_CH3P_C_H2_H,
	e_CH3P_CH2_H,
	e_CH3P_CH_H_H,
	e_CH3P_CH_H2,
	e_H2OP_OH_H,
	e_H2OP_O_H_H,
	e_H2OP_O_H2,
	e_H3OP_H2O_H,
	e_H3OP_OH_H_H,
	e_H3OP_OH_H2,
	e_H3OP_O_H2_H,
	e_HCOP_CO_H,
	e_OHP_O_H,
	e_SiH2P_SiH_H,
	e_SiH2P_Si_H_H,
	e_SiH2P_Si_H2,
	e_SiOHP_SiO_H,
	H2_CH_CH3_nu,
	H2_CH3P_CH5P_nu,
	H2s_CH_CH3_nu,
	Hminus_CH2_CH3_e,
	Hminus_CH3_CH4_e,
	nu_CH3_CH2_H,
	nu_CH3_CH_H2,
	nu_CH4_CH3_H,
	nu_CH4_CH2_H2,
	nu_CH4_CH_H2,
	crnu_CH3_CH2_H,
	crnu_CH3_CH_H2,
	crnu_CH4_CH2_H2,
	e_CH5P_CH3_H2,
	e_CH5P_CH4_H,
	e_CH4P_CH3_H,
	e_CH4P_CH2_H_H,
	H2_N_NH_H  ,      
	H2_NH_NH2_H ,        
	H2_NH2_NH3_H , 
	H2_CN_HCN_H   ,      
	HP_HNO_NOP_H2,
	HP_HS_SP_H2,
	H_HSP_SP_H2 ,
	H2P_N_NHP_H  ,     
	H2_NP_NHP_H   ,    
	H2_NHP_N_H3P   ,    
	H2P_NH_NH2P_H   ,  
	H2_NHP_NH2P_H    , 
	H2_NH2P_NH3P_H    , 
	H2_NH3P_NH4P_H     ,
	H2P_CN_HCNP_H     ,
	H2_CNP_HCNP_H     ,
	H2P_NO_HNOP_H      ,
	H2_SP_HSP_H        ,
	H2_CSP_HCSP_H      ,
	H3P_NH_NH2P_H2     ,
	H3P_NH2_NH3P_H2    ,
	H3P_NH3_NH4P_H2    ,
	H3P_CN_HCNP_H2     ,
	H3P_NO_HNOP_H2     ,
	H3P_S_HSP_H2       ,
	H3P_CS_HCSP_H2     ,
	H3P_NO2_NOP_OH_H2  ,
	HP_NH_NHP_H        ,
	HP_NH2_NH2P_H      ,
	HP_NH3_NH3P_H      ,
	H_CNP_CN_HP        ,
	HP_HCN_HCNP_H      ,
	H_HCNP_HCN_HP      ,
	H_N2P_N2_HP        ,
	HP_NO_NOP_H        ,
	HP_HS_HSP_H        ,
	HP_SiN_SiNP_H      ,
	HP_CS_CSP_H         ,
	HP_NS_NSP_H        ,
	HP_SO_SOP_H        ,
	HP_OCS_OCSP_H      ,
	HP_S2_S2P_H     ,
	H2P_NH_NHP_H2    ,
	H2P_NH2_NH2P_H2   ,
	H2P_NH3_NH3P_H2   ,
	H2P_CN_CNP_H2      ,
	H2P_HCN_HCNP_H2    ,
	H2P_NO_NOP_H2     ,
	H2_ClP_HClP_H	,
	H2_HClP_H2ClP_H	,
	H3P_Cl_HClP_H2,
	H3P_HCl_H2ClP_H2,
	HP_HCl_HClP_H,
	HP_C2_C2P_H,
	H2_S_HS_H,
	H2P_C2_C2P_H2,
	Hminus_NH4P_NH3_H2,
	Hminus_NP_N_H,
	HP_C2H2_C2H2P_H, 
	HP_C2H2_C2HP_H2 ,
	HP_C3H_C3HP_H ,
	HP_C3H_C3P_H2 ,
	H2P_C2H_C2H2P_H, 
	H2P_C2H2_C2H2P_H2, 
	H3P_C2H_C2H2P_H2 ,
	H3P_C3_C3HP_H2 ,
	H2_C2HP_C2H2P_H ,
	H2_C3P_C3HP_H,
	H_C2H3P_C2H2P_H2 ,
	H3P_C2H2_C2H3P_H2 ,
	H2P_C2H2_C2H3P_H ,
	HP_C3_C3P_H ,
	HP_C2H_C2HP_H, 
	H2P_C2_C2HP_H ,
	H2P_C2H_C2HP_H2, 
	H3P_C2_C2HP_H2 ,
	H2_C2P_C2HP_H ,
	HP_C2H_C2P_H2,
	N2_H3P_N2HP_H2;
}	co;

EXTERN struct t_mole {

		
	/** limit to the ratio H2/Htot - if ratio is below this, large atom is not called */
	double H2_to_H_limit;

	/** the number of electronic quantum states to include.
	* To do both Lyman and Werner bands want nelec = 3 */
	long int n_h2_elec_states;

	/** this is option to use estimates of the collision rates from g-bar approximations */
	/** turn mole.lgColl_gbar on/off with atom h2 gbar on off */
	bool lgColl_gbar;

	/** this is option to turn off the calculated collision rates */
	bool lgColl_deexec_Calc;

	/** this is option to turn off guesses of collisional dissociation rates */
	bool lgColl_dissoc_coll;

	/** include collision rates that come from real calculations,
	 * off with atom h2 collisions off command */
	bool lgH2_grain_deexcitation;

	/** flag to force LTE level populations, atom H2 LTE */
	bool lgH2_LTE;

	/** option to turn off ortho-para collisions, command ATOM H2 COLLISIONS ORTHO PARA OFF */
	bool lgH2_ortho_para_coll_on;

	/** which set of He - H2 collisions to use? default is ORNL, other
	 * is Le BOURlet */
	bool lgH2_He_ORNL;
	
	/*>>chng 08 feb 27, GS
	 * flag saying whether (true) or not to use ORNL H2 - H2 collisions*/
	bool lgH2_ORH2_ORNL;
	bool lgH2_PAH2_ORNL;

	/** turn on trace information */
	int nH2_TRACE;

	/** put noise into collision rates */
	bool lgH2_NOISE ,
		/** noise for the CR collisions */
		lgH2_NOISECOSMIC;

	/** this sets how fine a trace we want for atom  h2 trace */
	int nH2_trace_final , 
		nH2_trace_iterations , 
		nH2_trace_full,
		nH2_trace_matrix;
	
	 /** do we include capture of molecules onto grain surfaces?  default is true,
	  * turned off with NO GRAIN MOLECULES command */
	 bool lgGrain_mole_deplete;

	/** std and mean for the noise, log normal distribution */
	double xMeanNoise , xSTDNoise;

	/** flag saying whether an element is in the chemistry network */
	bool lgElem_in_chemistry[LIMELM];
	int num_comole_calc, num_comole_tot, num_elements;
	
	/** these are source and sink terms for the ionization ladder, for chemical
	 * processes that remove and add species */
	double **source , **sink;
	
	 /** rate s-1 for molecular charge transfer, nelem from to */
	realnum ***xMoleChTrRate;/***[LIMELM][LIMELM+1][LIMELM+1];*/
	double **amat, /* [NUM_COMOLE_CALC][NUM_COMOLE_CALC],  */
		*b, /* [NUM_COMOLE_CALC],  */
		**c; /* [NUM_COMOLE_TOT][NUM_COMOLE_CALC + 1]; */
}	mole;

enum {CHARS_SPECIES=7};
/* Structure containing molecule data, initially only CO */
EXTERN struct molecule {
	int nElem[LIMELM];	/** number of O, C, Si, N, S, and e- in each molecule */
	int nelem_hevmol;   /** this is the atomic number MINUS ONE of the main element within the molecule */
	char label[CHARS_SPECIES];      /** molecule name */
	int nElec;
	int Excit;
	bool lgGas_Phase;    /** is this in solid or gas phase? */
	int n_nuclei;       /** total number of nuclei in species */
	realnum hevmol;       /** the density (cm-3) of each species */
	realnum hev_reinit;   /* this save first valid solution from previous iteration */
	realnum *location;    /** Location of density in non-molecule code, if it exists */
	realnum hevcol;       /** total column density in this iteration */
	realnum hevcol_old;   /** column density in previous iteration */
	realnum pdr_mole_co;  /** default solution for initialization */
	/* realnum HevMolSav; ** the particle densities in first zone where CO computed,
										* in last iteration */
	realnum xMoleFracMax;
	realnum mole_mass;
	realnum co_save;	    /** previous solution to molecular network */
	realnum comole_save;
	realnum hevmol_save;
	int active;
	int index;
} **COmole;


extern struct molecule *findspecies(const char buf[]);

extern void CO_punch_mol(FILE *punit, const char chSpecies[], 
						 char header[], double depth);

#endif /* MOLE_H_ */
