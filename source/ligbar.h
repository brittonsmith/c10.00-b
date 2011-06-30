/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef LIGBAR_H_
#define LIGBAR_H_

 
 /**ligbar obtain collision strength for any Li-sequence line 
 \param ized 
 \param *t2s2p
 \param *t2s3p
 \param *cs2s2p
 \param *cs2s3p
 */
void ligbar(long, transition * t2s2p , transition * t2s3p ,double*,double*);

#endif /* LIGBAR_H_ */
