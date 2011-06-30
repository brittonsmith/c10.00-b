/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef HE_H_
#define HE_H_

/** he.h saves vars that are only for he */
EXTERN struct t_he {

	/** array index wihtin continuum array for Bowen 374A resonance line */
	long int ip374, 
	  ip660;

	/** these remember total He0 dest rate relative to exits from triplets */
	/** largest fraction of He0 dest due to ionization from 23S */
	double frac_he0dest_23S, 
	 /** fraction of that max due to photoioniztion */
	  frac_he0dest_23S_photo;
	/** zone where this max happened */
	long nzone;

	}	he;


#endif /* HE_H_ */
