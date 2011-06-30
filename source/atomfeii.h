/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef ATOMFEII_H_
#define ATOMFEII_H_

/**\file atomfeii.h
  these routines are in FeIILevelPops.c */

 /** 
called by assert feii depart coef command
\param   *pred
\param   *BigError 
\param   *StdDev 
 */ 
void AssertFeIIDep( double *pred , double *BigError , double *StdDev );

/** reads in feii data from disk, creates space for main arrays */ 
void FeIICreate(void);
/** FeIIPrint */
void FeIIPrint(void);

void FeIILevelPops(void);

 /**
    called in LineSet4, this sums over FeII bands and returns intensities args are lower and upper edges of bands as set in fe2bands.dat
    \param wl1 
    \param wl2
  * \param *SumBandInward
 */ 
double FeIISumBand(realnum wl1, 
				   realnum wl2,
				   double *SumBandInward );

/**FeII_RT_TauInc called once per zone in RT_tau_inc to increment large FeII atom line optical depths */
void FeII_RT_TauInc(void);

/**FeII_RT_tau_reset reset optical depths for large FeII atom, called by update after each iteration  */
void FeII_RT_tau_reset(void);

/**FeIIPoint called by ContCreatePointers to create pointers for lines in large FeII atom */
void FeIIPoint(void);


 /** 
  called by rt_line_driving to compute radiative acceleration due to FeII lines
  \param *fe2drive
 */ 
void FeIIAccel(double *fe2drive);

 /**
   do line RT for FeII model
 */ 
void FeII_RT_Make( void );

/** called by LineSet4, this adds feii line intensities together */ 
void FeIIAddLines( void );

/** called by parse_punch, save level energies and stat weights 
\param ioPUN
*/ 
void FeIIPunchLevels( FILE * ioPUN );

 /**
  FeIIPunchColden save level energies, stat weights, column density
 \param *ioPUN file we will save to
 */ 
void FeIIPunchColden(
  FILE *ioPUN );

 /** 
  FeII_Colden maintain H2 column densities within X
  \param *chLabel 
 */ 
void FeII_Colden( const char *chLabel );


 /** called by FeIIPunchLevels, this creates the save feii line intensity
 \param ioPUN 
 */ 
void FeIISaveLines( FILE * ioPUN );

 /** called by FeIIPunchLevels, this creates the save feii line optical depths
  \param ioPUN 
 */ 
void FeIIPunchOpticalDepth( FILE * ioPUN );

/** initialize optical depth arrays, called by TauOut */
void FeII_LineZero(void);

/**FeIIIntenZero zero out intensity of FeII atom */ 
void FeIIIntenZero(void);

/** rad pre due to FeII lines called in PresTotCurrent*/
double FeIIRadPress(void);

/** internal energy of FeII called in PresTotCurrent 
\return 
*/
double FeII_InterEnergy(void);

/**
save some departure coef for large atom, set with save feii departure command
\param ioPUN
\param lgDoAll option to save all dep coef if true
 */ 
void FeIIPunDepart(FILE* ioPUN , 
	bool lgDoAll );

void PunFeII( FILE * io );

 /** 
  send the departure coef for physical level nPUN to unit ioPUN
  \param ioPUN the io unit where the print should be directed
  \param nPUN the physical (not c) number of the level
 */ 
void FeIIPun1Depart(FILE * ioPUN , long int nPUN );
/**
 save line data for FeII atom
 \param ioPUN io unit for save
 \param lgDoAll save all levels if true, only subset if false
*/ 
void FeIIPunData(
	FILE* ioPUN ,
	bool lgDoAll );

/**
save some level pops for large atom, set with save feii level populations command
\param ioPUN 
\param lgPunchRange save range of levels if true, only selected subset if false
\param ipRangeLo if ipPunchRange is true, this is lower bound of range on C scale 
\param ipRangeHi if ipPunchRange is true, this is upper bound of range on C scale 
\param lgPunchDensity flag saying whether to save density (cm-3, true) or relative population (flase)
*/ 
void FeIIPunPop(FILE* ioPUN , 
	bool lgPunchRange ,
	long int ipRangeLo ,
	long int ipRangeHi ,
	bool lgPunchDensity );


 /**
include FeII lines in punched optical depths, etc, called from SaveLineStuff 
\param io 
\param xLimit
\param index
 */ 
void FeIIPunchLineStuff( FILE * io , realnum xLimit  , long index);

#if 0

/** 
send the level pops for physical level nPUN to unit ioPUN
\param ioPUN the io unit where the print should be directed 
\param nPUN the physical (not c) number of the level
*/ 
void FeIIPun1Pop( 
	FILE * ioPUN , 
	long int nPUN );
#endif

/** zero out variables that deal with FeII, called by zero */
void FeIIZero(void);

/** initialize some variables, called by zero */
void FeIIReset(void);

/** do OTS and outward parts of FeII lines, if large atom is turned on */
void FeII_OTS(void);

/** do outward rates for FeII, called by RT_diffuse */
void FeII_RT_Out(void);

/**ParseAtomFeII parse the atom feii command */
class Parser;

void ParseAtomFeII( Parser &p );

/** this is the number of levels for the large FeII atom */
#define	NFE2LEVN	371

/** this is set true when space is allocated for the FeII arrays,
 * once this happens the number of levels cannot be changed with the atom feii levels command 
 * set false in cddefines */
extern bool lgFeIIMalloc;

EXTERN struct t_FeII {

	/** number of levels for the large FeII atom, changed with the atom feii levels command 
	 * set to NFE2LEVN in cddefines */
	long int nFeIILevel_local;
	/** this remembers number of FeII levels allocated when MALLOC first called.  */
	long int nFeIILevel_malloc;

	/** this flag is true if the full FeII is being done and we do not want to use the
	 * simple FeII UV approximation.  Is false if only low levels of FeII are being done
	 * with full atom, and still want to use the simple atom */
	bool lgFeIILargeOn;

	/** option to always evaluate model atom, set with SLOW key on atom feii command */
	bool lgSlow;

	/** option to print calls to FeIILevelPops, set with print key on atom feii */
	bool lgPrint;

	/** option to only simulate calls to FeIILevelPops */
	bool lgSimulate;

	/** say which FeII atom this is, Verner or Netzer */
	char chFeIIAtom[7];

	/** save verner short for shorter save */
	bool lgShortFe2;

	/** says whether (true) or not (false) pumping of the FeiI model atom by HI Lya is included
	 * normally true, set false with the NO FEII command */
	bool lgLyaPumpOn;

	/** energy range for FeII lines output with save verner, 
	 * set with save verner range e ryd, e ryd */
	realnum fe2ener[2];

	/** energy range and threshold for FeII lines output with save verner */
	realnum fe2thresh;

	/** these are the lower and upper bounds to the FeII continuum, in Angstroms,*/
	realnum feconwlLo , feconwlHi;
	/** the number of intervals to break the FeII continuum into */
	long int nfe2con;

	/** redistribution function to use for resonance and subordinate lines */
	int ipRedisFcnResonance;
	int ipRedisFcnSubordinate;

	/** cooling computed by large FeII model atom
	 * Fe2_large_cool is total cooling (or heating if negative)
	 * and ddT_Fe2_large_cool is its derivative wrt temperature */
	double Fe2_large_cool, 
	  ddT_Fe2_large_cool, 
	  Fe2_large_heat;

	/** total cooling due to 16 level model FeII atom
	 * fe2cool is total cooling due to 16 level atom, */
	double 
		ddT_Fe2_UVsimp_cool ,
		Fe2_UVsimp_cool,
		for7;

	}	FeII;

	/* this info used to estimate destruction of Lya by overlap with FeII
 * in case where large atom is not turned on */

/**  number of FeII lines in Fred's atom */
const int NFEII = 373;
/** number of points in partition function table */
const int NFE2PR = 61;

class t_fe2ovr_la : public Singleton<t_fe2ovr_la>
{
	friend class Singleton<t_fe2ovr_la>;
protected:
	t_fe2ovr_la();
private:
	realnum fe2lam[NFEII];
	realnum fe2osc[NFEII];
	realnum fe2enr[NFEII];
	realnum fe2gs[NFEII];

	long int ipfe2[NFEII];

	/** opacity and optical depths for ground of Fred's FeII atom */
	realnum feopc[NFEII];
	realnum Fe2TauLte[NFEII];
	realnum Fe2PopLte[NFEII];

	double fe2pt[NFE2PR];
	double fe2pf[NFE2PR];

	/** fe2par evaluate FeII partition function */
	double fe2par(double te);
public:
	void zero_opacity();

	void init_pointers();

	/** tau_inc: update line opacities */
	void tau_inc();

	/** atoms_fe2ovr compute FeII overlap with Lya */
	void atoms_fe2ovr(void);
};

#endif /* ATOMFEII_H_ */
