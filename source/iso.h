/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef ISO_H_
#define ISO_H_

/**\file iso.h - information for isoelectronic sequences */

EXTERN long int max_num_levels;

/** temperatures at which CS for He1 are stored */
#define HE1CSARRAY 10

/** This macro is used to zero any radiative process with photon energy below
 * the plasma frequency.  The energy must be in Rydbergs!	*/
#define KILL_BELOW_PLASMA(E_)		( (rfield.lgPlasNu && ((E_)<rfield.plsfrq) ) ? 0.:1. )

#define NUM_DR_TEMPS	19

/** these macros are just an easy way to return the quantum numbers of a given level. */
#define N_(A_)	(StatesElemNEW[nelem][nelem-ipISO][A_].n)
#define L_(A_)	(StatesElemNEW[nelem][nelem-ipISO][A_].l)
#define S_(A_)	(StatesElemNEW[nelem][nelem-ipISO][A_].S)
#define J_(A_)	(StatesElemNEW[nelem][nelem-ipISO][A_].j)

#define IPRAD	0
#define IPCOLLIS	1
/*define IPENERGY	2*/

/* following two macros used to define recombination coef arrays */
/* Max n desired in RRCoef file.	*/
/** this is the number of levels used with the
 * atom xx-like levels large command */
/* Hydrogen and helium atoms will have precompiled recombination coefficients up to these maximum n.	*/
#define RREC_MAXN	40

/** Ions of the sequences will go up to this n, h-like He will get same as iso roots.	*/
#define LIKE_RREC_MAXN( A_ )		( A_ == ipHELIUM ? 40 : 20 )

#define N_ISO_TE_RECOMB		41

/** This is the n to go up to when calculating total recombination.	Any change 
 * here will not be reflected in total recomb until "compile xxlike" is run	*/
#define SumUpToThisN	1000
/** the magic number for the table of recombination coefficients, YYMMDD */
#define RECOMBMAGIC		(90521)

/**iso_cascade - calculate cascade probabilities, branching ratios, and associated errors
\param ipISO
\param nelem
*/
void iso_cascade( long ipISO, long nelem );

/**iso_charge_transfer_update - update rate coefficients for CT of H and He with everything else
 */
void iso_charge_transfer_update(void);

/**iso_collapsed_bnl_print - print departure coefficients for collapsed levels 
\param ipISO
\param nelem
*/
void iso_collapsed_bnl_print( long ipISO, long nelem );

/**iso_collapsed_bnl_set - set departure coefficients for collapsed levels 
\param ipISO
\param nelem
*/
void iso_collapsed_bnl_set( long ipISO, long nelem );

/**iso_collapsed_Aul_update - update decays from collapsed levels 
\param ipISO
\param nelem
*/
void iso_collapsed_Aul_update( long ipISO, long nelem );

/**iso_collapsed_lifetimes_update - update lifetimes of collapsed levels 
\param ipISO
\param nelem
*/
void iso_collapsed_lifetimes_update( long ipISO, long nelem );

/** iso_collide - calculate collision data for ipISO, nelem  
\param ipISO
\param nelem
*/
void iso_collide( long ipISO, long nelem );

/** iso_collisional_ionization - calculate collisional ionization rate for ipISO, nelem  
\param ipISO
\param nelem
*/
void iso_collisional_ionization( long ipISO, long nelem );

/** iso_continuum_lower - limit max prin. quan. no. due to continuum lowering processes 
\param ipISO
\param nelem
*/
void iso_continuum_lower( long ipISO , long nelem );

/**iso_cool compute net heating/cooling due to hydrogenc atom species 
\param ipISO the isoelectronic sequence, 0 for H 
\param nelem is element, so 0 for H itself
*/
void iso_cool( long ipISO , long nelem );

/**iso_create create storage space data for iso sequences, 1 one time per coreload 
*/
void iso_create( void );

/**iso_cross_section get cross section for a particular level of an iso sequence ion
\param ERyd 
\param EthRyd
\param n
\param l
\param S
\param Z 
\param ipISO
*/
double iso_cross_section( double ERyd , double EthRyd, long n, long l, long S, long globalZ, long globalISO );

/**iso_departure_coefficients - calculate departure coefficients
\param ipISO
\param nelem
*/
void iso_departure_coefficients( long ipISO, long nelem );

/**iso_dielec_recomb_rate - get state-specific dielectronic recombination rate 
\param ipISO
\param nelem
\param ipLo
*/
double iso_dielec_recomb_rate( long ipISO, long nelem, long ipLo );

/**iso_drive updates rates and solves level populations for all isoelectronic atoms and ions 
*/
void iso_drive( void );

/**iso_error_generation generate gaussian errors 
\param ipISO
\param nelem
*/
void iso_error_generation( long ipISO, long nelem );

/**iso_get_total_num_levels - get total number of levels with the given number of resolved and collapsed
\param ipISO
\param nmaxResolved
\param numCollapsed
*/
long iso_get_total_num_levels( long ipISO, long nmaxResolved, long numCollapsed );

/**IonHydro this controls hydrogen atomic and molecular crosstalk
*/
void IonHydro( void );

/**iso_ionize_recombine evaluate state specific creation and destruction processes 
\param ipISO
\param nelem
*/
void iso_ionize_recombine( long ipISO , long nelem );

/**iso_level solve for iso-sequence ionization balance 
\param ipISO
\param nelem
*/
void iso_level( const long ipISO, const long nelem);

/**iso_photo do photoionization rates for element nelem on the ipISO isoelectronic sequence 
\param ipISO
\param nelem
*/
void iso_photo( long ipISO , long nelem );

/**iso_prt_pops routine to print level pops or departure coefficients for iso sequences 
\param ipISO
\param nelem
\param lgPrtDeparCoef
*/
void iso_prt_pops( long ipISO, long nelem, bool lgPrtDeparCoef );

/**iso_put_error put an error bar on a piece of data, to be used with Gaussian random noise gen 
\param ipISO
\param nelem
\param ipHi
\param ipLo
\param whichData
\param errorOpt
\param errorPess
*/
void iso_put_error( long ipISO, long nelem, long ipHi, long ipLo, long whichData, realnum errorOpt, realnum errorPess);

/**iso_radiative_recomb - get rad recomb rate coefficients for iso sequences.
\param ipISO
\param nelem
*/
void iso_radiative_recomb( long ipISO, long nelem );


/**iso_radiative_recomb_effective - get effective recomb rate coefficients into each level (including indirect)
\param ipISO
\param nelem
*/
void iso_radiative_recomb_effective( long ipISO, long nelem );

/**iso_recomb_check - called by SanityCheck to confirm that recombination coef are ok,
 * return value is relative error between new calculation of recom, and interp value 
 \param ipISO
 \param nelem the chemical element, 1 for He
 \param level the level, 0 for ground
 \param temperature the temperature to be used
*/
double iso_recomb_check( long ipISO, long nelem, long level, double temperature );

/** iso_recomb_auxiliary_free - free up some auxiliary space associated with iso recombination tables.
*/
void iso_recomb_auxiliary_free( void );

/** iso_recomb_malloc - malloc space needed for iso recombination tables.
*/
void iso_recomb_malloc( void );

/** iso_recomb_setup - read in or compile iso recombination tables.
\param ipISO
*/
void iso_recomb_setup( long ipISO );

/** iso_RRCoef_Te - interpolate iso recomb coeff as function of temperature
\param ipISO
\param nelem
\param n
*/
double iso_RRCoef_Te( long ipISO, long nelem , long n );

/**iso_satellite_update - update iso satellite line information 
*/
void iso_satellite_update( void );

/* calculate radiative lifetime of an individual iso state 
\param ipISO
\param nelem
\param n
\param l
*/
double iso_state_lifetime( long ipISO, long nelem, long n, long l );

/**iso_solve - main routine to call iso_level and determine iso level balances
\param ipISO
*/
void iso_solve( long ipISO, long nelem );

/**iso_suprathermal - calculate secondary excitation by suprathermal electrons for iso sequences 
\param ipISO
\param nelem
*/
void iso_suprathermal( long ipISO, long nelem );

/**iso_update_num_levels - update level informations for iso sequences 
\param ipISO
\param nelem
*/
void iso_update_num_levels( long ipISO, long nelem );

/**iso_update_rates routine to set up iso rates, level balance is done elsewhere 
\param ipISO
\param nelem
*/
void iso_update_rates(long ipISO, long nelem );

EXTERN struct t_iso 
{
	bool lgPrintNumberOfLevels;

	const char *chISO[NISO];

	/** arrays for stark broadening in Puetter formalism */
	multi_arr<realnum,4> strkar;
	multi_arr<double,4> pestrk;

	/** Find index given quantum numbers
	 * Since separate j levels within a triplet term are only resolved in the case of 2tripP,
	 * allocating memory for a j dimension is unwarranted.  Instead 
	 * iso.QuantumNumbers2Index[ipISO][nelem][2][1][1] will point to 2^3P2, with 2^3P0 and 2^3P1
	 * easily accessed by subtracting 2 and 1 respectively from the returned index.	*/
	multi_arr<long,5> QuantumNumbers2Index;

	/** number of Lyman lines to include only as opacity sources, in each iso seq,
	 * all now set to 100 in zero.c */
	long int nLyman[NISO],
		/** number of levels actually malloced - probably greater than above */
		nLyman_malloc[NISO];

	/** a set of array indices for all atoms on the iso sequences,
	 * ipIsoLevNIonCon[ipISO][ipZ][n] */
	multi_arr<long int,3> ipIsoLevNIonCon;

	/** ionization potential of level N in Ryd */
	multi_arr<double,3> xIsoLevNIonRyd;

	/** the ratio of ion to atom for all iso species
	 * xIonSimple is simple estimate, should agree at low density */
	double xIonSimple[NISO][LIMELM];

	/** option to turn off l-mixing collisions */
	bool lgColl_l_mixing[NISO];

	/** option to turn off collisional excitation */
	bool lgColl_excite[NISO];

	/** option to turn off collisional ionization */
	bool lgColl_ionize[NISO];

	/** option to print departure coefficients */
	bool lgPrtDepartCoef[NISO][LIMELM];

	/** option to print level populations */
	bool lgPrtLevelPops[NISO][LIMELM];

	/** do thermal average of collision strenghts if true, false by default,
	 * set true with SET COLLISION STRENGTHS AVERAGE command */
	bool lgCollStrenThermAver;

	/** flag saying whether induced two photon is included 
	 * in the level pops for H- and He-like */
	bool lgInd2nu_On;

	/* option to disable continuum lowering due to stark broadening, particle packing, etc. */
	bool lgContinuumLoweringEnabled[NISO];

	/** true if the model is full size in the sense that all levels below
	 * continuum are considered.	*/
	bool lgLevelsLowered[NISO][LIMELM];

	/** This variable is set to true if the continuum was lowered at any point in the calculation.
	 * Necessary because some models will lowered continuum at intermediate points but not last zone. */
	bool lgLevelsEverLowered[NISO][LIMELM];

	/* flag that says we must reevaluate everything about this ion */
	bool lgMustReeval[NISO][LIMELM];

	/** the number of collapsed levels, these lie on top of resolved levels */
	long int nCollapsed_max[NISO][LIMELM];
	long int nCollapsed_local[NISO][LIMELM];

	/** total number of collapsed and resolve levels, 
	 * iso.numLevels_max[ipISO] is derived from total resolved and collapsed levels 
	 * it is the maximum number of levels ever to be used in this core load. */
	long int numLevels_max[NISO][LIMELM];

	/** total number of levels with continuum pressure lowering included 
	 * this varies from zone to zone, and from model to model, but cannot
	 * exceed numLevels_max  */
	long int numLevels_local[NISO][LIMELM];

	/** number of levels malloc'd in the core load, can't go over that later 
	 * in later sims can lower number of levels but not raise them  */
	long int numLevels_malloc[NISO][LIMELM];

	/** principal quantum number n of the highest resolved level */
	long int n_HighestResolved_max[NISO][LIMELM];
	/** the local (pressure lowered) version of the above */
	long int n_HighestResolved_local[NISO][LIMELM];

	/** statistical weight of the ground state of the parent ions for each 
	 * species, used for Milne relation and recombination */
	realnum stat_ion[NISO];

	/** the induced upward two-photon rate */
	double TwoNu_induc_up[NISO][LIMELM];

	/** the induced downward two-photon rate */
	double TwoNu_induc_dn[NISO][LIMELM];

	/** the largest induced downward two photon rate */
	double TwoNu_induc_dn_max[NISO][LIMELM];

	/** radiative recombination rate coefficient, RadRecomb[ipISO][ipZ][n][fcn]
	 * iso.RadRecomb[ipISO][ipZ][ipLo][ipRecEsc] escape prob
	 * iso.RadRecomb[ipISO][ipZ][n][ipRecNetEsc] net escape prob, accounting for absorption 
	 * iso.RadRecomb[ipISO][ipZ][ipLo][ipRecRad] rate coef, cm^3 s^-1 
	 * */
	multi_arr<double,4> RadRecomb;

	/** total radiative recombination continuum, RadRecCon[ipISO][nelem][n]
	 * units erg cm-3 s-1 */
	multi_arr<double,3> RadRecCon;

	/** tells whether dielectronic recombination is turned on	*/
	bool lgDielRecom[NISO];

	/** state specific dielectronic recombination rates, DielecRecomb[ipISO][ipZ][n]
	 * iso.DielecRecomb[ipISO][ipZ][ipLo] rate coef, cm^3 s^-1 */
	multi_arr<double,3> DielecRecomb;

	/** state specific dielectronic recombination rates as a function of temperature
	 * iso.DielecRecombVsTemp[ipISO][ipZ][ipLo][Temp] rate coef, cm^3 s^-1 */
	multi_arr<double,4> DielecRecombVsTemp;

	/** difference between actual case b photons in rtdiffuse, and correct case b */
	realnum CaseBCheck[NISO][LIMELM];

	/** case b recombination rate coefficient */
	double RadRec_caseB[NISO][LIMELM];

	/** the total effective radiative recombination rate coefficient (cm3 s-1), 
	 * radiative rate with correction for absorption and ionization */
	double RadRec_effec[NISO][LIMELM];

	/** all processes from level n to the continuum, units s-1 */
	multi_arr<double,3> RateLevel2Cont; 

	/** all processes from the continuum to level n, units s-1 */
	multi_arr<double,3> RateCont2Level; 

	/** ratio of collisional recombination rate to recom from all processes */
	double RecomCollisFrac[NISO][LIMELM];

	/** ipOpac pointers for photoionization cross sections of hydrogen
	 * iso.ipOpac[NISO][LIMELM][NHPLPHOT] */
	multi_arr<long int,3> ipOpac;

	/** continuum to total opacity factors for each level */
	multi_arr<double,3> ConOpacRatio;

	/** departure coefficient */
	multi_arr<double,3> DepartCoef;

	/** true is all lte populations positive for Hydrogenic atoms */
	bool lgPopLTE_OK[NISO][LIMELM];

	/** hlte is lte population of each level, cm^3 */
	multi_arr<double,3> PopLTE;

	/** collisional ionization rate coefficient from each level (cm3 s-1) */
	multi_arr<double,3> ColIoniz;

	/** net free bound cooling for this element */
	double FreeBnd_net_Cool_Rate[NISO][LIMELM];

	/** net cooling due to collisional ionization */
	double coll_ion[NISO][LIMELM];

	/** net cooling due to collisional excit of higher lines */
	double cRest_cool[NISO][LIMELM];

	/** net cooling due to total collisional excit of lines */
	double xLineTotCool[NISO][LIMELM];

	/** deriv of net cooling due to total collisional excit of lines */
	double dLTot[NISO][LIMELM];

	/** net cooling due to rad rec */
	double RadRecCool[NISO][LIMELM];

	/** net cooling due to collisional excit of balmer lines */
	double cBal_cool[NISO][LIMELM];

	/** net cooling due to collisional excit of higher lyman lines */
	double cLyrest_cool[NISO][LIMELM];

	/** net cooling due to collisional excit of Lya */
	double cLya_cool[NISO][LIMELM];

	/** photoionization rate, gammnc[iso][nelem][level] */
	multi_arr<double,3> gamnc;

	/** RecomInducRate will become induced recombination rate coefficient
	 * when multipled by lte population.  
	 * integral of photorate times exp(-hu/kt) 
	 * for ind rec, produced by gamma routine needs to be mult
	 * by lte pop to become real rate  */
	multi_arr<double,3> RecomInducRate;

	/** RecomInducCool_Coef becomes rate coef for incuded recombination cooling,
	 * when multipled by lte population. 
	 * this times hnu-hnuo0 to get cooling,
	 * evaluated in gamma routine and saved */
	multi_arr<double,3> RecomInducCool_Coef;

	/** the actual induced recom cooling rate, erg cm-3 s-1 */
	double RecomInducCool_Rate[NISO][LIMELM];

	/** Boltzmann factor from lower to upper level, [ISO][nelem][up][lo] */
	multi_arr<double,4> Boltzmann;

	/** photoelectric heating rate */
	multi_arr<double,3> PhotoHeat;

	/** this is the rate for the Aul given to bogus transitions,
	 * set to 1e-30 in zero */
	/** >>chng 04 may 17, esd 1e-20, changed to 1e-30 to allow
	 * rydberg levels to be treated with their small As */
	realnum SmallA;

	/** will become array of indices for induced two photon,
	 * series of symmetric indices 
	 * ipHy2nu[ipISO][ipZ][energy] */
	multi_arr<long int,3> ipSym2nu;

	/** will become array of two photon transition probabilities per energy bin	*/
	multi_arr<realnum,3> As2nu;

	/** These are pointers to the energies representing the two-photon gap,
	 * and half the gap, respectively.  */
	long ipTwoPhoE[NISO][LIMELM];
	long ipHalfTwoPhoE[NISO][LIMELM];

	/** types of redistribution functions for Lya, other resonances, and subordinate lines */
	int ipLyaRedist[NISO] , ipResoRedist[NISO] , ipSubRedist[NISO];

	/** this is the upper level for Lya */
	int nLyaLevel[NISO];

	/** flag to set which type of solution was used for level pops, "zero" or "popul" */
	char chTypeAtomUsed[NISO][LIMELM][10];

	/***********************************************/
	/***********************************************/

	/*   everything below here was originally part */
	/*   of the former helike structure, now       */
	/*   generalized for iso					   */
	/** all of these are initialized in zero */

	/***********************************************/
	/***********************************************/

	/** number of CS in the above array */
	long int nCS[NISO];

	/** flag set by compile he-like command, says to regenerate table of recombination coef */
	bool lgCompileRecomb[NISO];

	/** flag set by atom he-like no recomb interp command, 
	 * says to generate recombination coefficients
	 * on the fly */
	bool lgNoRecombInterp[NISO];

	/** parameters for changing gbar - set with set hegbar command */
	bool lgCS_Vriens[NISO] ,
		lgCS_None[NISO] ,
		lgCS_Vrinceanu[NISO],
		lgCS_therm_ave[NISO];
	int nCS_new[NISO];//vals are 0, 1, and 2

	/** used to print warning if density too low for first collapsed level to be l-mixed	*/
	bool lgCritDensLMix[NISO];

 	/** flag saying whether to include fine-structure mixing in spontaneous decays	
	 * set with ATOM HE-LIKE FSM command */
	bool lgFSM[NISO];

	/** error in sum of As. */
	multi_arr<double,3> SigmaAtot;
	
	/** array of collision strengths read from data file...this is interpolated upon.	*/
	multi_arr<realnum,5> HeCS;
	
	/** vector of temperatures corresponding to collision strengths stuffed into HeCS.	*/
	double CSTemp[NISO][HE1CSARRAY];

	/** This flag is set to true if the rates should be treated with a randomly generated error,
	 * on the range specifically set for each rate, before being entered into the rate matrix.	*/
	bool lgRandErrGen[NISO];

	/** this is flag saying that random gaussians have already been set...they should only
	 * be done once per model, and this must be reset to false at the beginning of each model.	*/
	bool lgErrGenDone[NISO][LIMELM];

	bool lgPessimisticErrors;

	bool lgTopoff[NISO];

	/** This is the used to set a unique seed in parallel gaussian runs */
	int modelRank[NISO];

	/** This is the array in which uncertainties are stored if iso.lgRandErrGen[NISO] is set. */
	/* first dimension is iso,
	 * second is nelem,
	 * third is upper level,
	 * fourth is lower level,
	 * fifth is for radiative, collisional, or energy errors.
	 * MACROS are used for the last dimension: IPRAD, IPCOLLIS, and IPENERGY. */
	multi_arr<realnum,5> Error;

	/** This is the array in which gaussian errors are generated, using the values in 
	 * the Error array above as the standard deviations */
	multi_arr<realnum,5> ErrorFactor;

	/** the effective collisional rate from 2S, for h-like and he-like sequences */
	double qTot2S[NISO][LIMELM];

	/** total brancing ratio and standard deviation in it */
	multi_arr<double,4> BranchRatio, CascadeProb, SigmaCascadeProb;

	/** effective recombination and standard deviation in it */
	multi_arr<double,3> RadEffec, SigmaRadEffec;

	multi_arr<realnum,5> CachedAs;
	/* the departure coefficients of collapsed levels */
	multi_arr<double,5> bnl_effective;

} iso;

#endif /* ISO_H_ */
