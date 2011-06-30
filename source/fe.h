/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef FE_H_
#define FE_H_

/**\file fe.h - vars that are for Fe only */

/** number of levels in model of Fe III */
#define NLFE3  14

/** number of levels in model of Fe VII */
#define	NLFE7	8

/** number of levels in model of Fe 10 Fe X */
#define NLFE10 5

/** number of levels in model of Fe 11 Fe XI */
#define NLFE11 5

/** number of levels in model of Fe 13 Fe XIII */
#define NLFE13 5

EXTERN struct t_fe {

	/** total cooling due to 12 level model FeIV atom, and three uv lines
	 * 4P5/2 - 6S5/2, 4P3/2 - 6S5/2, and 4D+ - 6S
	 * Fe4CoolTot is total cooling due to entire atom, also three selectec lines */
	double Fe4CoolTot, 
	  fe42836, 
	  fe42829, 
	  fe42567, 
	  fe40401, 
	  fe41207, 
	  fe41206, 
	  fe41106, 
	  fe41007, 
	  fe41008, 
	  fe40906;

	/** cooling, wavelengths, and emission in Fe 3 multi level atom */
	double Fe3CoolTot ,
	  **Fe3_wl , **Fe3_emiss;

	/** calculted in Fe 7 multi level atom */
	double Fe7CoolTot ,
	  **Fe7_wl , **Fe7_emiss;

	/** cooling, wavelengths, and emission in Fe 11 multi level atom */
	double Fe11CoolTot ,
	  **Fe11_wl , **Fe11_emiss;

	/** cooling, wavelengths, and emission in Fe 13 multi level atom */
	double Fe13CoolTot ,
	  **Fe13_wl , **Fe13_emiss;

	/** fekhot, fekcld, number of photons in hot and cold iron, per unit vol */
	realnum fekhot, 
	  fekcld;

	/** Fe Ka from iron in grains */
	realnum fegrain;

	/** uv pumping of fe coronal lines */
	long int ipfe10;
	realnum pfe10, 
	  pfe11a, 
	  pfe11b, 
	  pfe14;

	}	fe;


#endif /* FE_H_ */
