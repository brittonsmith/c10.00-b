/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef H2_PRIV_H_
#define H2_PRIV_H_

#include "transition.h"

/*H2_Read_Cosmicray_distribution read distribution function for H2 population following cosmic ray collisional excitation
void H2_Read_Cosmicray_distribution(void); */

/** this is the energy, above which levels are considered to be H2*, and below which
 * they are H2 - it is converted from eV into wavenumbers since energy scale
 * energy_wn[][][] is in these terms */
/* >> chng 05 jul 15, TE, H2g = sum (v=0, J=0,1) */
/* >>chng 05 jul 29, to 0.5 eV, this goes up to J=8 for v=0 */
/* >>chng 05 aug 03, slight upward change in energy to include the J=8 level,
 * also give energy in waveumbers for simplicity (save h2 levels give energy in ryd) */
/*#define	ENERGY_H2_STAR	(0.5/EVRYD/WAVNRYD)*/
/* energy of v=0, J=8 is 4051.73, J=9 is 5001.97
 * v=1, J=0 is 4161.14 */
const double ENERGY_H2_STAR = 4100.;

/**H2_He_coll Interpolate the rate coefficeints 
 * The range of the temperature is between 2K - 1e8K 
\param init
\param final
\param temp
*/
double H2_He_coll(int init, int final, double temp);

/**H2_He_coll_init receives the name of the file that contrains the fitting coefficeints 
 * of all transitions and read into 3d vectors. It outputs 'test.out' to test the arrays
 * return value is magic number 
 \param FILE_NAME_IN[]
*/
long int  H2_He_coll_init(const char FILE_NAME_IN[]);

/*>>chng 08 feb 27, GS - parameters for the ORNL H2 - H2 collision data*/
double H2_ORH2_coll(int init, int final, double temp);
double H2_PAH2_coll(int init, int final, double temp);
long int  H2_ORH2_coll_init(const char FILE_NAME_IN[] );
long int  H2_PAH2_coll_init(const char FILE_NAME_IN[] );

/** read energies for all electronic levels 
\param nelec
*/
void H2_ReadEnergies( long int nelec );

/** read dissociation probabilities and kinetic energies for all electronic levels 
\param nelec
*/
void H2_ReadDissprob( long int nelec );

/** H2_CollidRateEvalAll - set H2 collision rates */
void H2_CollidRateEvalAll( void );

/** read collision rates 
\param nColl
*/
void H2_CollidRateRead( long int nColl );

/** read transition probabilities 
\param nelec
*/
void H2_ReadTransprob( long int nelec );

/**H2_Read_hminus_distribution read distribution function for H2 population following formation from H minus */
void H2_Read_hminus_distribution(void);

/**mole_H2_form find state specific rates grains and H- form H2 */
void mole_H2_form( void );

/**mole_H2_LTE sets Boltzmann factors and LTE unit population of large H2 molecular */
void mole_H2_LTE( void );

/**H2_Solomon_rate find rates between H2s and H2g and other levels,
 * for use in the chemistry */
void H2_Solomon_rate( void );

/**H2_gs_rates evaluate rates between ground and star states of H2 for use in chemistry */
void H2_gs_rates( void );

/** H2_zero_pops_too_low - zero out some H2 variables if we decide not to compute
 * the full sim, called by H2_LevelPops*/
void H2_zero_pops_too_low( void );

const bool CR_PRINT = false;
const int CR_X = 1;
const int CR_VIB = 15;
const int CR_J = 10;
const int CR_EXIT = 3;

/** the number of different types of colliders 
 *  */
const int N_X_COLLIDER = 5;

/** labels for the colliders */
const int chN_X_COLLIDER = 10;
EXTERN char chH2ColliderLabels[N_X_COLLIDER][chN_X_COLLIDER];

/** this is the highest vib state that has collision data */
const int VIB_COLLID = 3;

/** the number of temperature points in the data file */
const int nTE_HMINUS = 7;

/* these vars are private for H2 but uses same style as all other header files -
 * the EXTERN is extern in all except cddefines */

/** number of levels in H2g */
EXTERN long int nEner_H2_ground;

EXTERN multi_arr<double,3> H2_populations;
EXTERN multi_arr<double,3> H2_rad_rate_out;

/** total population in each vib state */
EXTERN multi_arr<double,2> pops_per_vib;

/** the renorm factor for this H2 to the chemistry - should be unity */
EXTERN double H2_renorm_chemistry,
	H2_sum_excit_elec_den;

/** column density within X only vib and rot */
EXTERN multi_arr<realnum,2> H2_X_colden;

/** rates [cm-3 s-1] from elec excited states into X only vib and rot */
EXTERN multi_arr<double,2> H2_X_rate_from_elec_excited;

/** rates [s-1] to elec excited states from X only vib and rot */
EXTERN multi_arr<double,2> H2_X_rate_to_elec_excited;

/** rate [s-1} for collisions from ihi to ilo */
EXTERN multi_arr<realnum,2> H2_X_coll_rate;

/** LTE column density within X only vib and rot */
EXTERN multi_arr<realnum,2> H2_X_colden_LTE;

/** the number of ro-vib levels in each elec state */
EXTERN long int nLevels_per_elec[N_H2_ELEC];

/** the total population in each elec state */
EXTERN double pops_per_elec[N_H2_ELEC];

/** energy in wavenumbers */
EXTERN multi_arr<double,3> energy_wn;

/** this is the actual rate, cm^3 s^-1, for each collider
 * CollRate[coll_type][vib_up][rot_up][vib_lo][rot_lo] */
EXTERN multi_arr<realnum,6> CollRateFit;

/** these will mostly become xxx[elec][vib][rot] */
EXTERN multi_arr<realnum,3> H2_dissprob;
EXTERN multi_arr<realnum,3> H2_disske;
EXTERN multi_arr<realnum,5> H2_CollRate;

/** these will mostly become xxx[elec][vib][rot] */
EXTERN multi_arr<double,3> H2_old_populations;
EXTERN multi_arr<double,3> H2_Boltzmann;
EXTERN multi_arr<double,3> H2_populations_LTE;
/** this is total statistical weight, including nuclear spin */
EXTERN multi_arr<realnum,3> H2_stat;
/** this is true if state is para, false if ortho */
EXTERN multi_arr<bool,3> H2_lgOrtho;

EXTERN long int nzoneAsEval , iterationAsEval;

EXTERN multi_arr<int,2> H2_ipPhoto;
/*EXTERN realnum **H2_col_rate_in_old;
EXTERN realnum **H2_col_rate_out_old;*/
EXTERN multi_arr<double,2> H2_col_rate_in;
EXTERN multi_arr<double,2> H2_col_rate_out;
EXTERN multi_arr<double,2> H2_rad_rate_in;
EXTERN realnum *H2_X_source;
EXTERN realnum *H2_X_sink;

/** distribution function for formation on grain surfaces,
 * vib, rot, last dim is grain type */
EXTERN multi_arr<realnum,3> H2_X_grain_formation_distribution;

/** formation into specific states within X only vib and rot,
 * includes both H- and H2 routes */
EXTERN multi_arr<realnum,2> H2_X_formation;

/** backwards destruction of v,J levels due to the H- route */
EXTERN multi_arr<realnum,2> H2_X_Hmin_back;

/** save rate coef (cm3 s-1) for collisional dissociation */
EXTERN multi_arr<realnum,2> H2_coll_dissoc_rate_coef;

/** save rate coef (cm3 s-1) for collisional dissociation with H2g and H2s*/
EXTERN multi_arr<realnum,2> H2_coll_dissoc_rate_coef_H2;

/** density of H2s and H2g during current iteration */
EXTERN double H2_den_s , H2_den_g;

/** vib, rot, last dim is temperature */
EXTERN multi_arr<realnum,3> H2_X_hminus_formation_distribution;

/** these are energies and indices for levels within X ,
 * a vector of sorted energies*/
EXTERN realnum *H2_energies;
/** the total number of energies */
EXTERN long int nH2_energies;
EXTERN long int *H2_ipX_ener_sort;
EXTERN long int *ipVib_H2_energy_sort, *ipElec_H2_energy_sort;
EXTERN long int *ipRot_H2_energy_sort;
EXTERN multi_arr<long int,3> ipEnergySort;

/** number of levels within X which are done with matrix solver,
 * set with atom h2 matrix command */
EXTERN long int nXLevelsMatrix;

/** this is array of accumulated line intensities, used for save he lines command */
EXTERN multi_arr<realnum,6> H2_SaveLine;

/** fully defined array saying whether (true) or not (false) a radiative decay 
 * is defined by the standard emission line structure */
EXTERN multi_arr<bool,6> lgH2_line_exists;

/** counters used by H2_itrzn to find number of calls of h2 per zone */
EXTERN long int nH2_pops;
EXTERN long int nH2_zone;

/** this is used to establish zone number for evaluation of number of levels in matrix */
EXTERN long int nzone_nlevel_set;

/** the number of times the H2 molecules has been called in this iteration.  For the
 * very first call we will use lte for the level H2_populations, for later calls
 * use the last solution */
EXTERN long int nCallH2_this_iteration;

/** the remainder can't be EXTERN since init to values in h2.c */
extern int H2_nRot_add_ortho_para[N_H2_ELEC];

extern double H2_DissocEnergies[N_H2_ELEC];

/** temperature where H- distribution are set */
extern realnum H2_te_hminus[nTE_HMINUS];

typedef multi_arr<transition,6>::iterator mt6i;
typedef multi_arr<transition,6>::const_iterator mt6ci;

#endif /* H2_PRIV_H_ */
