/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef NITRO_H_
#define NITRO_H_

EXTERN struct t_nitro {

	/** this is fraction of 2D excitations that emit a photon in the 5200 multiplet
	 * this is used for collisional quenching in prt lines routine */
	double quench_5200,

		/** evaluated in cool_nitr, used in prt lines */
		xN5200, 
		xN5198, 
		xN3467, 
		xN3466,
		xN10398,
		xN10397,
		xN10408,
		xN10407,
		c5755, 
		c6584, 

		/** ratio of radiative to total decays from n=3 of NII */
		xN2_A3_tot, 

		//! pump of [NI] 5199 by FUV intercombination lines, only an estimate
		pump5199,

		//! recombination pumping of [NI] 5199, only an estimate
		rec5199;

	}	nitro;


#endif /* NITRO_H_ */
