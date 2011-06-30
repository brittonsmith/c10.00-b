/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef ABUND_H_
#define ABUND_H_


 /**
  AbundancesSet sets initial abundances after parameters are entered by reading input
 */ 
void AbundancesSet(void);

/**
  AbundancesPrt print all abundances, both gas phase and grains
 */ 
void AbundancesPrt( void );

 /**
  AbundancesZero set initial abundances for different mixes
 */ 
void AbundancesZero(void);

/**
  generate abundance set from Fred Hamann's starburst evolution grid 
 \param chCard
*/
class Parser;
void abund_starburst(Parser &p);

 /**
  AbundancesTable interpolate on table of points to do 'element table' command, 
  \param  r0 
  \param  depth 
  \param  iel
 */ 
double AbundancesTable(double r0, 
  double depth, 
  long int iel);

/** abund.h */
EXTERN struct t_abund {

	/** logical flag saying whether to include this element in save output for AGN tables */
	bool lgAGN[LIMELM];

	realnum SolarSave[LIMELM], 
	  OldSolar84[LIMELM], 
	  anova[LIMELM], 
	  apn[LIMELM], 
	  ahii[LIMELM], 
	  camern[LIMELM], 
	  aprim[LIMELM], 
	  aism[LIMELM],
	  aCrab[LIMELM];

	bool lgAbnSolar;

	bool lgElmONapn[LIMELM], 
		lgElmONahii[LIMELM], 
		lgElmONaism[LIMELM],
		lgElmONaCrab[LIMELM];

	/** solar abundances for the current calculation */
	realnum solar[LIMELM];

	/**lgAbunTabl says whether this element is to have its abundance
	 *determined from a table (true) or stored constant (false)
	 *set true with element table command */
	bool lgAbunTabl[LIMELM], 

	  /** lgAbTaDepth says whether depth or radius, true is depth */
	  lgAbTaDepth[LIMELM], 

	  /** general flag saying this option turned on */
	  lgAbTaON;

#	define	LIMTABD	500

	/**AbTabFac abundances for element table*/
	realnum AbTabFac[LIMTABD][LIMELM], 

	/**AbTabRad depth scale 
	 *parameters for dlaw table command*/
	  AbTabRad[LIMTABD][LIMELM];

	long int nAbunTabl;

	/** indices so that abundances can be in any order */
	long int ipSolar[LIMELM], 
	  npSolar;

	/** scale factors to alter abundances of elements, set with element scale */
	realnum ScaleElement[LIMELM];

	/** Depletion is set of stored scale factors for depletion of general ism */
	realnum Depletion[LIMELM], 

	/** depset is unity unless depletion is used */
	  depset[LIMELM];

	/** lgDepln is true if depln used */
	bool lgDepln;

	/** scale factor for metals, set with metals command	 */
	realnum ScaleMetals;

	}	abund;



#endif /* ABUND_H_ */
