/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef IONBAL_H_
#define IONBAL_H_

#include "dense.h"

/**	ion_recom_calculate called by conv_base to calculate radiative and dielectronic
recombination rate coefficients */
void ion_recom_calculate( void );

/**ion_trim raise or lower most extreme stages of ionization considered 
\param nelem element number on the C scale, 5 for C
*/
void ion_trim(
	long int nelem );

/**ion_zero zero out heating and charge transfer save arrays */
void ion_zero(long int nelem);

/**ion_collis fill in collisional ionization rates, and resulting cooling 
\param nelem element number on C scale, H is 0
*/
void ion_collis(
	long nelem);

/**ion_solver solve the bi-diagonal matrix for ionization balance 
\param nelem - element number on C scale, He is 1
\param lgPrintIt - option to print details of matrix elements
*/
void ion_solver(long int nelem, bool lgPrintIt);

/**ion_solver same as above but solves two elements simultaneously
 \param nelem1
 \param nelem2
 \param lgPrintIt - option to print details of matrix elements
 */
void ion_solver(long int nelem1, long int nelem2, bool lgPrintIt);

/**ion_photo fill array PhotoRate with photoionization rates for heavy elements 
\param nelem is atomic number on C scale, 0 for H
\param lgPrintIt debugging flag to turn on print
*/
void ion_photo(
	long int nelem , 
	bool lgPrintIt );

/**ion_recomb generate recombination coefficients for any species */
void ion_recomb(bool,const double*,const double*,const double[],const double[],
		const double[],const double[],const double[],const double[],long);

/**ion_recombAGN generate recombination coefficients for AGN table */
void ion_recombAGN( FILE * io );

/**ion_wrapper a wrapper that redirects to IonHelium, IonCarbo, etc.. */
void ion_wrapper( long nelem );

/**Badnell_rec_init This code is written by Terry Yun, 2005 *
 * It reads rate coefficient fits into 3D arrays and output array.out for testing *
 * The testing can be commented out */
void Badnell_rec_init( void );

/* routines to do heavy element ionization balance */
void IonAlumi(void);
void IonArgon(void);
void IonBeryl(void);
void IonBoron(void);
void IonCalci(void);
void IonCarbo(void);
void IonChlor(void);
void IonChrom(void);
void IonCobal(void);
void IonCoppe(void);
void IonFluor(void);
void IonHelium(void);
void IonIron(void);
void IonLithi(void);
void IonMagne(void);
void IonManga(void);
void IonNeon(void);
void IonNicke(void);
void IonNitro(void);
void IonOxyge(void);
void IonPhosi(void);
void IonPotas(void);
void IonScand(void);
void IonSilic(void);
void IonSodiu(void);
void IonSulph(void);
void IonTitan(void);
void IonVanad(void);
void IonZinc(void);

/** max number of shells we ever have to deal with */
#define NSHELLS	7

/** class for vars dealing with ionization balance */
class t_ionbal 
{
public:
	/** limits for highest and lowest stages of ionization in ion_trim,
	 * these are set with command "set trim xx" where xx is log of upper
	 * and lower ionization fractions.  if only one number then both are
	 * set to it.  These variables are used in trimStages to adjust the
	 * range of ionization. <BR>
	 * limit to fractional abundance of high stage of ionization,
	 * set to 1e-6 in zero.c */
	double trimhi, 

	/** limit to fractional abundance of low stage of ionization,
	 * set to 1e-10 in zero.c */
	  trimlo;

	/** option to turn off upward ionization trimming, with set trim upper off  */
	bool lgTrimhiOn;

	/* ==============================================================
	 * all following deals with ionization processes */

	/** \verbatim store photoionization rates for all shells of all elements
	  first dim is nelem, the atomic number of element on the c scale, H is 0.
	  second dim is stage of ionization, on the c scale, atom is 0.
	  third dim is shell number, K shell is 0, valence shell depends on ion, up to 7
	  last dim: 0 is photo rate (s-1)
	            1 is low energy heating
	            2 is high energy (secondary-capable) total heating
	            both will be multiplied by ion abundance to get vol rates 
	  some special last pairs - 
	  [x][0][10][0] pair production in highen \endverbatim
	 */
	double ****PhotoRate_Shell/**[LIMELM][LIMELM][7][3]*/;

	/** set to 1 in zero, so have no effect, 
	 * set to 0 with 'no photoionization' command, 
	 * kills photoionization of everything */
	bool lgPhotoIoniz_On;

	/** should H - O charge transfer be done in ionization or chemistry? 
	 * default is chemistry, true */
	bool lgHO_ct_chem;

	/** collisional ionization rate for CollidRate[nelem][ion][0], s-1
	* cooling, erg/s in CollidRate[nelem][ion][1] */
	double ***CollIonRate_Ground/**[LIMELM][LIMELM][2]*/;

	/** cosmic ray ionization rate */
	double CosRayIonRate;

	/** cosmic ray heating rate - erg s-1 - must multiply by density of 
	 * absorbers - neutral hydrogen to get volume rate */
	double CosRayHeatNeutralParticles;

	/** cosmic ray heating of thermal electrons - must multiply by electron
	 * density to obtain erg cm-3 s-1 */
	double CosRayHeatThermalElectrons;

	/** local heating rate due to some "extra" process */
	double ExtraHeatRate;

	/** heating erg s-1 due to fast neutrons - energy flux times cross section
	 * but does not include density */
	double xNeutronHeatRate;

	/** ionization and heating due to pair production */
	double PairProducPhotoRate[3];

	/* ==============================================================
	 * following deal with Compton recoil ionization of bound electrons */
	
	/** flag saying that Compton recoil ionization of bound
	 * electrons is enabled,
	 * set false with no recoil ionization command */
	bool lgCompRecoil;

	/** the local heating due to Compton recoil ionization */
	double CompRecoilHeatLocal;

	/** array indices for continuum offset of Compton recoil ionization threshold */
	long int **ipCompRecoil;

	/** rate of bound electron ionization by Compton scattering */
	double **CompRecoilIonRate;

	/** save rate of bound electron ionization by Compton scattering */
	double **CompRecoilIonRateSave;

	/** heating rate due to bound electron ionization by Compton scattering */
	double **CompRecoilHeatRate;

	/** save heating rate due to bound electron ionization by Compton scattering */
	double **CompRecoilHeatRateSave;

	/** inner shell UTA ionization rate, includes autoionization probability */
	double **UTA_ionize_rate;
	/** inner shell UTA heating rate */
	double **UTA_heat_rate;

	/** this says whether to include inner shell absorption lines */
	bool lgInnerShellLine_on;
	/** says whether to include the new Romas data set */
	bool lgInnerShell_Kisielius;
	/** this says whether to replace the Behar 01 data with the Bu et al. 06
	 * data - default is true, to do so, set false with SET UTA BEHAR command */
	bool lgInnerShell_Gu06;

	/** stage-to-stage ionization rates (s-1), all processes
	 * dimensions [nelem][from_ion][to_ion] */
	double ***RateIoniz;

	/** number of valence electrons that can participate - multiplies since 
	 * electron rate */
	long int  nCompRecoilElec[LIMELM];

	double CompHeating_Max;
	/* ==============================================================
	 * end Compton recoil ionization of bound electrons */

	/* ==============================================================
	 * all following deals with recombination */

	/** total recombination rate (s-1) all processes */
	double **RateRecomTot;

	/** rate coefficients [cm3 s-1] for Badnell DR recombination */
	double **RR_Badnell_rate_coef ,
		   **DR_Badnell_rate_coef,
		   *DR_Badnell_rate_coef_mean_ion;

	/** do these rate coefficients exist? */
	int **lgDR_Badnell_rate_coef_exist ,
		**lgRR_Badnell_rate_coef_exist;

	/** do we use new Badnell rates? */
	bool lgDR_recom_Badnell_use,
		lgRR_recom_Badnell_use,
		/** option to print rates then exit */
		lgRecom_Badnell_print;

	/** rate coefficients [cm3 s-1] for older DR recombination */
	double **DR_old_rate_coef;

	/** radiative recombination rate coefficient (cm3 s-1) used by code */
	double **RR_rate_coef_used,
		   **DR_rate_coef_used;

	/** radiative recombination rate coefficient returned from Dima Verner's routine */
	double **RR_Verner_rate_coef;

	/** three cases for S DR - 
	 * 0, default larger of guess and Badnell
	 * 1, pure Badnell
	 * 3, scaled oxygen */
	int nDR_S_guess;
	realnum DR_S_scale[5];

	/** rate for recombination and ionization on grain surfaces */
	realnum **GrainCreat,
		**GrainDestr;

	int 
	 /**grecon usually true, set to 0 with no grain neutralization command	 */
	  lgGrainIonRecom;

	/** suppression factors for dielectronic recombination
	 * 1 is burgess and 2 is Storey */
	realnum DielSupprs[2][LIMELM];

	/** flag for guess of entire range of dr - false by default, true with kludge steve */
	realnum lg_guess_coef;
	/** log normal noise for guess, zero by default, turned on with noise option */
	realnum guess_noise;

	/** logical flag for suppression of dielectronic recombination
	 * 1 is burgess, 2 is Nussbaumer and Storey	 */
	bool lgSupDie[2];

	/** following all for 3-body recombination */
	/** lgNoCota flag set with no three body recombination */
	bool lgNoCota;

	/** the actual rates */
	realnum CotaRate[LIMELM];

	/** these are error flags for three-body recombination */
	long int ilt, 
	  iltln, 
	  ilthn, 
	  ihthn, 
	  ifail;

	// find the total ionization rate (including multiple-electron processes) 
	double RateIonizTot( long nelem, long ion )
	{
		double sum = 0.;
		
		for( long ion_to=ion+1; ion_to<=dense.IonHigh[nelem]; ion_to++ )
			sum += RateIoniz[nelem][ion][ion_to];

		return sum;
	}
};

#if 0
double t_ionbal::RateIonizTot(long nelem, long ion)
{
	double sum = 0.;
	
	for( long ion_to=ion+1; ion_to<=dense.IonHigh[nelem]; ion_to++ )
		sum += t_ionbal::RateIoniz[nelem][ion][ion_to];

	return sum;
}
#endif 

EXTERN t_ionbal ionbal;

#endif /* IONBAL_H_ */
