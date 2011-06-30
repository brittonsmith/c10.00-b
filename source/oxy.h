/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef OXY_H_
#define OXY_H_

/** oxy.h */
EXTERN struct t_oxy {
	realnum poiii2, 
	  d5007r, 
	  d5007t, 
	  poiii3, 
	  poiii2Max, 
	  poiii3Max, 
	  d4363, 
	  r5007Max, 
	  r4363Max, 
	  poiexc, 
	  poimax, 
	  d6300;

	double p1666, 
	  p1401;

	/** photoionization rate for OII, producing 1665 */
	realnum AugerO3;

	/** lines produced by relaxation following photoionization */
	realnum s3727, 
	  s7325;

	/** array indices for photoionization into excited state of oii */
	long int i2d, 
	  i2p;

	realnum o3cs12, 
	  o3cs13, 
	  o3ex23, 
	  o3br32, 
	  o3enro;

	}	oxy;


#endif /* OXY_H_ */
