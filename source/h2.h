/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef H2_H_
#define H2_H_


/** this is the number of electronic levels */
#define N_H2_ELEC	7

/** create H2 molecules, called by ContCreatePointers */
void H2_Create(void);

/** set the ipCont struc element for the H2 molecule, called by ContCreatePointers */
void H2_ContPoint( void );

/**H2_DR choose next zone thickness based on H2 big molecule */
double H2_DR(void);

/** H2_Init - called by cdInit to init H2 */
void H2_Init(void);

/** H2_init_coreload one time initialization */
void H2_init_coreload( void );

/** radiative acceleration due to H2 called in rt_line_driving */
double H2_Accel(void);

/**H2_RT_OTS - add H2 ots fields */
void H2_RT_OTS( void );

/** rad pre due to h2 lines called in PresTotCurrent*/
double H2_RadPress(void);

/** add in explicit lines from the large H2 molecule, called by lines_molecules */
void H2_LinesAdd(void);

/**H2_Reset called to reset variables that are needed after an iteration */
void H2_Reset( void );

/** internal energy of H2 called in PresTotCurrent */
double H2_InterEnergy(void);

/**H2_Colden maintain H2 column densities within X 
\param *chLabel
*/
void H2_Colden( const char *chLabel );

/**H2_cooling evaluate cooling and heating due to H2 molecule,
 * string is name of calling routine
 *\param *chString name of calling routine 
*/
void H2_Cooling(const char *chString);

/** save H2 line data 
\param ioPUN io unit for save
\param lgDoAll save all levels if true, only subset if false
*/
void H2_Punch_line_data(
	FILE* ioPUN ,
	bool lgDoAll );

/** include H2 lines in punched optical depths, etc, called from SaveLineStuff 
\param io
\param xLimit
\param index
*/
void H2_PunchLineStuff( FILE * io , realnum xLimit  , long index);

/** do emission from H2 - called from RT_diffuse */
void H2_RT_diffuse(void);

/** do RT for H2 lines 
*/
void H2_RTMake( void );

/** increment optical depth for the H2 molecule, called from RT_tau_inc */
void H2_RT_tau_inc(void);

/** zero out vars in the large H2 molecule, called from zero */
void H2_Zero( void );

/**H2_Prt_Zone print H2 info into zone results, called from prtzone for each printed zone */
void H2_Prt_Zone(void);

/** initialize optical depths in H2, called from RT_tau_init */
void H2_LineZero( void );

/** the large H2 molecule, called from RT_tau_reset */
void H2_RT_tau_reset( void );

/** do level populations for H2, called by iso_solve */
void H2_LevelPops( void );

/** save some properties of the large H2 molecule 
\param io
\param chJOB[]
\param chTime[]
\param ipPun
*/
void H2_PunchDo( FILE* io , char chJOB[] , const char chTime[] , long int ipPun );

/** print line optical depths, called from premet in response to print line optical depths command*/
void H2_Prt_line_tau(void);

class Parser;

/**H2_ParseSave parse the save h2 command */
void H2_ParseSave( Parser &p ,
				   char *chHeader);

/**H2_itrzn - average number of H2 pop evaluations per zone */
double H2_itrzn( void );

/**H2_Prt_column_density print H2 info into zone results, called from prtzone for each printed zone 
\param *ioMEAN this is stream used for io, is stdout when called by final,
	is save unit when save output generated
*/
void H2_Prt_column_density( FILE *ioMEAN );

/** flag saying whether molecular data have been read in yet */
extern bool lgH2_READ_DATA;

EXTERN struct t_h2 {

	/** the density (cm-3) of ortho H2 */
	double ortho_density,
	/** the density (cm-3) of para H2 */
		para_density;

	/** column density in ortho and para H2 */
	double ortho_colden ,
		para_colden;

	/** these remember the largest and smallest factors needed to
	 * renormalize the H2 chemistry */
	double renorm_max ,
		renorm_min;

	/** this will say how many times the large H2 molecule has been called in this zone -
	 * if not called (due to low H2 abundance) then not need to update its line arrays */
	long int nCallH2_this_zone;

	/** flag saying whether to bother with the large H2 molecule at all,
	 * default is false, set true with atom h2 on command */
	bool lgH2ON;

	/** this is the number of electronic levels to include in the output - default is 1,
	 * only X.  changed with PRINT LINES H2 ELECTRONIC and option on PUNCH H2 LINES commands */
	int nElecLevelOutput;

	/** number of vib states within electronic states from
	 * >>refer	H2	energies	Abgrall, */
	long int nVib_hi[N_H2_ELEC];

	/** number of rotation levels within each elec - vib */
	long int nRot_hi[N_H2_ELEC][50];

	/** this gives the first rotational state for each electronic state - J=0 does
	 * not exist when Lambda = 1 */
	long int Jlowest[N_H2_ELEC];

	/* true to use 2007 set of H2 - H collision rate, false use 1999 */
	bool lgH2_H_coll_07;

}	h2;

#endif /* H2_H_ */
