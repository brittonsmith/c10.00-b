/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef CARB_H_
#define CARB_H_

/** pmp2s.h */
EXTERN struct t_carb {
	double p1909, 
	  p2326, 
	  c8727, 
  	  c9850;

	/** correction for depopulation of excited state of CI */
	realnum r9850;

	}	carb;


#endif /* CARB_H_ */
