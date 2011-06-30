/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef CONTINUUM_H_
#define CONTINUUM_H_

#include "thirdparty.h"


/**ContCreatePointers create pointers for lines and continua, one time per coreload */
void ContCreatePointers(void);

/**ContSetIntensity derive intensity of incident continuum, 
  \return returns 0 if continuum is ok, 1 if failure of some kind*/
int ContSetIntensity(void);

/**IncidentContinuumHere derive intensity of incident continuum*/
void IncidentContinuumHere( void );

/** set up continuum energy mesh if first call, otherwise reset to original mesh */
void ContCreateMesh(void);

/**ContNegative sanity check for negative continuum intensities */
void ContNegative(void);

/**ffun evaluate total flux for sum of all continuum sources 
 \param anu photon energy (Rydberg) where continuum is evaluated 
 \param frac_beam_time fraction of beamed continuum that is varies with time
 \param frac_beam_const fraction of beamed continuum that is constant
 \param frac_isotropic fraction of continuum that is isotropic
*/ 
double ffun(
			/* the energy in Rydbergs where the continuum will be evaluated */
			double anu , 
			/* fraction of beamed continuum that is varies with time */
			double *frac_beam_time,
			/* fraction of beamed continuum that is constant */
			double *frac_beam_const,
			/* fraction of continuum that is isotropic */
			double *frac_isotropic );

/**ffun version without fractions */
double ffun(double anu );

/**ffun1 derive flux at a specific energy, for one continuum 
\param anu photon energy (Rydberg) where continuum is evaluated 
 */
double ffun1(double xnu);

/*outsum sum outward continuum beams */
void outsum(double *outtot, double *outin, double *outout);

/**DrvContPump local continuum pumping rate radiative transfer for all lines 
 \param *t
 \param DopplerWidth
*/ 
double DrvContPump(transition * t, realnum DopplerWidth );

/**cont_gaunt_calc do table look up of gaunt factor 
\param temp
\param z
\param photon
*/
double cont_gaunt_calc( double, double , double );

EXTERN struct t_continuum {
	/** this is information needed to set the energy binning,
	 * full continuum is described by series of ranges where resolution is
	 * constant over that range */
	realnum *filbnd, 

	  *fildel, 

	  *filres;

	long int *ifill0, 
	  /**number of ranges entered for this continuum source*/
	  nrange; 

	/** each of these is the upper bound of an energy band,
	 * the first lowest bound is the low-energy limit of the code */
	double *StoredEnergy,
		/** the resolution, dE/E for each band */
		*StoredResolution;

	/** the number of bands read in */
	long int nStoredBands;

	/** factor to reset continuum resolution set in continuum_mesh.ini,
	 * default is unity, reset with set resolution command */
	double ResolutionScaleFactor;

	/** flag saying that parts of continuum are zero */
	bool lgCon0,
	  lgCoStarInterpolationCaution;

	/** TotalLumin is total intensity in incident continuum erg cm-2 s-1 */
	double TotalLumin, 
	  totlsv;

	/** the incident continuum at Hb and La */
	realnum cn4861, 
	  cn1216, 
	  sv4861, 
	  sv1216;

	realnum 
		fluxv,
		fbeta;

	/** these are number, labels, and bounds of continuum bands
	 * they are specified in continuum_bands.ini in the data dir */
	long int nContBand;
	char **chContBandLabels;
	realnum *ContBandWavelength;
	long int *ipContBandLow , *ipContBandHi;
	/** these are fractions of first and last bin to include in the 
	 * band */
	realnum *BandEdgeCorrLow , *BandEdgeCorrHi;

	/** this is highest energy where k-shell opacities are counted
	 * can be adjusted with "set kshell" command */
	long int KshellLimit;
	realnum EnergyKshell;

	/** the md5sum of the continuum_mesh.ini file, this will be used
	 * to check the energy mesh in grain opacity files, etc. */
	string mesh_md5sum;

	t_continuum()
	{
		nrange = 0;
		mesh_md5sum = MD5datafile( "continuum_mesh.ini" );
	}

}	continuum;

#endif /* CONTINUUM_H_ */
