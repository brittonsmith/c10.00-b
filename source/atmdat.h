/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef ATMDAT_H_
#define ATMDAT_H_


 /**
  atmdat_2phot_shapefunction two photon emission function for all atomic and ionic species 
  \param  EbyE2nu 
  \param  ipISO
  \param  nelem 
 */ 
double atmdat_2phot_shapefunction( double EbyE2nu, long ipISO, long nelem );

/**
	atmdat_2phot_rate hydrogen two photon emission, including induced emission 
	\param ipZ   atomic number on c scale
	\param ipISO  the iso sequence, must be either 0 (hydrogen) or 1 (he-like)
 */ 
void atmdat_2phot_rate(
	long ipZ , 
	long ipISO );

 /**
  atmdat_readin read in some data files, but only if this is very first call 
 */ 
void atmdat_readin(void);

 /**
  atmdat_CHIANTI_readin read in data from CHIANTI database files
  \param intNS
  \param chFileName
 */ 
void atmdat_CHIANTI_readin( long intNS, char *chFileName );

 /**
  atmdat_LAMDA_readin read in data from LAMDA database files
  \param intNS
  \param chFileName
 */ 
void atmdat_LAMDA_readin( long intNS, char *chFileName );


 /**
  atmdat_outer_shell determine outer shell, and statistical weights of that and higher ion, for any ion
  written by Dima Verner
  \param [in] iz  atomic number from 1 to 30
  \param [in] in  number of electrons from 1 to iz
  \param [out] *imax  number of the outer shell
  \param [out] *ig0   statistical weight of (iz,in) ground state
  \param [out] *ig1 statistical weight of (iz,in-1) ground state
  \author Dima Verner
 */ 
void atmdat_outer_shell(
  long int iz, 
  long int in,
  long int *imax,
  long int *ig0, 
  long int *ig1);

 /**
  atmdat fill in the HCharExcIonOf and Rec arrays with Kingdon's fitted CT with H, 
  \param  O_HIonRate
  \return value is rate O ionizes H
 */ 
void ChargTranEval( double *O_HIonRate );

/**
 sum up the charge transfer heating
 \return 
*/ 
double ChargTranSumHeat(void);

/*ChargTranPun save charge transfer rate coefficients */
 /** 
  save charge transfer rate coefficients
  \param ipPnunit 
  \param chSave 
 */ 
void ChargTranPun( FILE* ipPnunit , char* chSave );

 /**
  atmdat_dielrec_fe Dielectronic recombination rates for Fe from Arnaud & Raymond 1992
  \param  ion
  \param  t
 */ 
double atmdat_dielrec_fe(long int ion, double t);

/** this initializes the arrays containing the fitting coefficients,
 * called by OpacityCreateAll, done once per coreload */
void atmdat_H_phot_cs(void);

/**atmdat_3body derive three-body recombination coefficients */
void atmdat_3body(void);

 /** 
 general utility to read in line emissivities from the Storey & Hummer tables
 of case B emissivities.  
\param iHi the principal quantum numbers, .	 
\param iLo upper and lower levels in any order
\param iZ charge of ion, only 1 and 2 for now	 
\param TempIn temperature, must lie within the range of the table, which depends on the ion charge, and is 500 - 30,000K for hydrogen
\param DenIn the density and must lie within the range of the table
\param chCase case - 'a' or 'b'
 */ 
double atmdat_HS_caseB( 
	long int iHi, long int iLo, long int iZ, double TempIn,
	double DenIn,
	char chCase 
	);

/**atmdat_DielSupres derive scale factors for suppression of Burgess dielectronic recombination */
void atmdat_DielSupres(void);

// arrays for Hummer & Storey 98 He1 cross sections and energies
EXTERN double ****HS_He1_Xsectn;
EXTERN double ****HS_He1_Energy;

// arrays for TOPbase Helike cross sections and energies
EXTERN double *****OP_Helike_Xsectn;
EXTERN double *****OP_Helike_Energy;
EXTERN long ****OP_Helike_NumPts;

/* these are the vectors that store the original Hummer and Storey case B
 * line data for H and He - the declaration for the interpolator follows */
#define NHSDIM 15 /**< used for following vectors*/
#define NLINEHS 300  /**< dimension of array with lines*/
#define HS_NZ 8 /**< number of elements that can be read in */
#define NHCSTE	8 /**< number of temperature points in h_coll_str arrays */
#define NUM_HS98_DATA_POINTS 811

EXTERN struct t_atmdat {
	/**
	 * ion, nelem
	 * these arrays save the charge transfer ionization and recombination
	 * rates for the heavy elements onto hydrogen.  ionization is
	 * of the heavy element, and so is a recombination for hydrogen
	 * 
	 * HCharExcIonOf( ion , nelem ), HCharExcRecTo( ion , nelem )
	 * charge transfer ionization of atomic oxygen = HCharExcIonOf[ipOXYGEN][0]*hii
	 * charge transfer recombination of ionized oxygen = HCharExcRecTo[ipOXYGEN][0]*hi
	 * HCharHeatMax, HCharCoolMax are largest fractions of local heating
	 * or cooling due to ct
	 * HCharHeatOn usually 1, set to 0 with no CTHeat command
	 */
	/** CharExcIon is ionization, */
	/** [0] is Atom^0 + H+ => Atom+1 + H0
	 * [n] is Atom^+n + H+ => Atom^+n-1 + H0 */

	/** CharExcRec is recombination */
	/** [0] is Atom^+1 + H0 => Atom^0 + H^+
	 * [n] is Atom^+n+1 + H0 => Atom^+n + H^+ */
	double HCharExcIonOf[LIMELM][LIMELM+1], //(cm3 s-1)
	  HCharExcRecTo[LIMELM][LIMELM+1],		//(cm3 s-1)
	  HCharHeatMax, 
	  HCharCoolMax, 
	  HCharHeatOn;

	/* rate coefficient (cm3 s-1) for N+(3P) + H+ -> N(2D) + H+ charge transfer*/
	double HCharExcRecTo_N0_2D;

	/** this is total rate (s-1) for ct ionization and recombination of H */
	double HCharExcIonTotal,
	  HCharExcRecTotal;

	/** this is total rate (s-1) for ct ionization and recombination of He */
	double HeCharExcIonTotal,
	  HeCharExcRecTotal;

	/** this is the current ratio of ct ionization of H, relative to total dest rate*/
	double HIonFrac;

	/** this is the largest ratio of ct ionization of H, relative to total dest rate*/
	double HIonFracMax;

	/**
	 *these are the charge ionization and recombination arrays
	 *that go into the ionization balance
	 *Ion[ion] is ionization, Rec[ion] is recombination
	 *units cm3 s-1
	 */
	double HeCharExcIonOf[LIMELM][LIMELM], 
	  HeCharExcRecTo[LIMELM][LIMELM];

	/** Dalgarno H charge transfer rate coefficient for high stages of ionization
	 * default is 1.92e-9 in zero, reset with 'set charge transfer' command */
	double HCTAlex;

	/** variable to turn on or off ct ionization-recombination of
	 * all elements - set off with no charge transfer command */
	bool lgCTOn;

	/** these are the density and temperature mesh points on the
	 * original Hummer & Storey data, for H[0] and He[1], */
	double Density[2][HS_NZ][NHSDIM], 
		ElecTemp[2][HS_NZ][NHSDIM],
		/**emiss[ipTemp][ipDens][ipLevel]*/
		Emiss[2][HS_NZ][NHSDIM][NHSDIM][NLINEHS];

	/** saves the number of density temperature mesh points for the two cases for 
	 * the HS_NZ elements */
	long int nDensity[2][HS_NZ] , ntemp[2][HS_NZ] , ncut[2][HS_NZ];

	/** following will be set false if we ever stop over bounds of HS table
	 * for any element.  first index is case A [0] or case B [1] -
	 * second is element number */
	bool lgHCaseBOK[2][HS_NZ];

	/** related to highest stage of ionization needed for Cota recom */
	long int nsbig;

	/** by default, include collisional ionization, option to not include it,
	 * with "no collisional ionization" command */
	bool lgCollIonOn;

	/** wavelengths of Hummer & Storey case B lines for H - O 
	 * first dimension is atomic number of C scale, H is 0
	 * next two are upper and lower configurations on physics 
	 * scale - Lya is 2-1, Lyb is 3-1, Ha is 3-2, etc */
	realnum WaveLengthCaseB[8][25][24];

	/** wavelengths for HeI case b */
	realnum *CaseBWlHeI;
	/** number of HeI case b lines included */
	long int nCaseBHeI;

	/** true if CHIANTI database is enabled **/
	bool lgChiantiOn;
	/** true if LAMDA database is enabled **/
	bool lgLamdaOn;

	/** Length +1 of the version read from /data/chianti/VERSION **/
	static const int iVersionLength = 10;
	/** Chianti version read from /data/chianti/VERSION **/
	char chVersion[iVersionLength];

	}	atmdat;

typedef enum { PHFIT_UNDEF, PHFIT95, PHFIT96 } phfit_version;

class t_ADfA : public Singleton<t_ADfA>
{
	friend class Singleton<t_ADfA>;
protected: 
	t_ADfA();
private:
	phfit_version version;
	/* phfit.dat */
	long int L[7];
	long int NINN[30];
	long int NTOT[30];
	realnum PH1[7][30][30][6];
	realnum PH2[30][30][7];
	/* hpfit.dat */
	realnum PHH[NHYDRO_MAX_LEVEL][5];
	/* rec_lines.dat */
	realnum P[8][110];
	realnum ST[9][405];
	/* rad_rec.dat */
	realnum rrec[30][30][2];
	realnum rnew[30][30][4];
	realnum fe[13][3];
	/* h_rad_rec */
	realnum HRF[NHYDRO_MAX_LEVEL][9];
	/* h_phot_cs.dat */
	/** array of cross sections for photoionization of hydrogen at threshold,
	 * 0 is 1s, 1 is 2s, 2 is 2p, up to 400 */
	realnum STH[NHYDRO_MAX_LEVEL];
	/* coll_ion.dat */
	double CF[30][30][5];
	/* h_coll_str.dat */
	/** array of EIE cross sections for hydrogen atom.  
	 * For all E1 transitions nl - n'l', with n' < n <= 5 */
	/* >>refer	H1	cs	Anderson, H., Ballance, C.P., Badnell, N.R., & Summers, H.P., 
	 * >>refercon	2000, J Phys B, 33, 1255; erratum, 2002 */
	double HCS[14][10][8];
public:
	/** set_version set version of phfit data to be used
	    \param val
	*/
	void set_version(phfit_version val) { version = val; }

	/** get_version which version of phfit data should be used? */
	phfit_version get_version() const { return version; }

	/** ph1 access elements of PH1 data block with parameters for photoionization cross section fits
	    \param i
	    \param j
	    \param k
	    \param l
	*/
	realnum ph1(int i, int j, int k, int l) const { return PH1[i][j][k][l]; }

	/** sth array of cross sections for photoionization of hydrogen at threshold,
	    0 is 1s, 1 is 2s, 2 is 2p, up to 400
	    \param i
	*/
	realnum sth(int i) const { return STH[i]; }

	/** phfit this subroutine calculates partial photoionization cross sections
	    for all ionization stages of all atoms from H to Zn (Z=30)
	    \param nz
	    \param ne
	    \param is
	    \param e
	    \author Dima Verner
	*/
	double phfit(long int nz, long int ne, long int is, double e);

	/** hpfit state specific photoionization cross sections for model hydrogen atom
	    \param iz 
	    \param n 
	    \param e 
	    \author Dima Verner
	*/ 
	double hpfit(long int iz, long int n, double e);

	/** rec_lines effective recombination coefficients for lines of C, N, O, by D. Verner
	    \param  t 
	    \param  r
	    \author Dima Verner
	*/ 
	void rec_lines(double t, realnum r[][471]);

	/** rad_rec calculates rates of radiative recombination for all ions
	    \param iz nuclear number on physics scale
	    \param in number of recombined electrons
	    \param t temperature K
	    \author Dima Verner
	*/ 
	double rad_rec(long int iz, long int in, double t);

	/** H_rad_rec calculates state-specific recombination rates for H and H-like ions
	    \param iz 
	    \param n 
	    \param t
	    \author Dima Verner
	*/ 
	double H_rad_rec(long int iz, long int n, double t);

	/** coll_ion D Verner's routine to compute collisional ionization rate coefficients 
	    \param iz 
	    \param in 
	    \param t
	    \author Dima Verner
	*/
	double coll_ion(long int iz, long int in, double t);

	/** h_coll_str routine to grab H cross sections from Anderson et al. 2002.
	    \param ipLo 
	    \param ipHi 
	    \param ipTe
	*/
	realnum h_coll_str( long ipLo, long ipHi, long ipTe );
};

#endif /* ATMDAT_H_ */
