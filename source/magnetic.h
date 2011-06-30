/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef MAGNETIC_H_
#define MAGNETIC_H_

/* magnetic.h */

/**ParseMagnet parse magnetic field command  
\param *chCard
*/
class Parser;
void ParseMagnet(Parser &p );

/** Magnetic_init initialize magnetic field parameters */
void Magnetic_init(void);

/**Magnetic_reinit - reinitialized magnetic field at start of new iteration */
void Magnetic_reinit(void);

/**Magnetic_evaluate evaluate some parameters to do with magnetic field */
void Magnetic_evaluate(void);

EXTERN struct t_magnetic {

	/** this says mag fields are turned on */
	bool lgB;

	/** enthalpy density associated with field */
	double EnthalpyDensity;

	/** this is magnetic pressure at current location - this can be negative for ordered field */
	double pressure;

	/** energy density at current location, this is positive */
	double energydensity;

	}	magnetic;


#endif /* MAGNETIC_H_ */
