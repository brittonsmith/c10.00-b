/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef MEWECOEF_H_
#define MEWECOEF_H_

/* mewecoef.h */
/**
 *     these are expansions for collision strengths
 *     210 is total number of coef from mewe plus gaetz and salpeter
 */
EXTERN struct t_MeweCoef {
	realnum g[210][4];
	}	MeweCoef;


#endif /* MEWECOEF_H_ */
