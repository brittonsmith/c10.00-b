/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef HELIKE_CS_H_
#define HELIKE_CS_H_


/**HeCollid evaluate collisional rates 
\param nelem
*/
void HeCollid( long int nelem);

/**HeCSInterp interpolate on He1 collision strengths 
\param nelem
\param ipHi
\param ipLo
\param Collider
*/
realnum HeCSInterp( long int nelem,
				 long int ipHi,
				 long int ipLo,
				 long int Collider );

/**AtomCSInterp do the atom	
\param nelem
\param ipHi
\param ipLo
\param *factor
\param **where
\param Collider
*/
realnum AtomCSInterp( long int nelem,
				   long int ipHi,
				   long int ipLo,
				   realnum *factor,
				   const char **where,
				   long int Collider );

/**IonCSInterp do the ions	
\param nelem
\param ipHi
\param ipLo
\param *factor
\param **where
\param Collider
*/
realnum IonCSInterp( long int nelem,
				  long int ipHi,
				  long int ipLo,
				  realnum *factor,
				  const char **where,
				  long int Collider );

/* Three different collision treatments, based on
 * Seaton 1962;
 * Pengelly and Seaton 1964; and
 * Vrinceanu and Flannery 2001.
 */

/**CS_l_mixing - find rate for l-mixing collisions by protons, for neutrals 
  based on Seaton 1962
\param ipISO
\param nelem
\param ipLo
\param ipHi
\param temp
\param Collider
*/
double CS_l_mixing_S62(
	long int ipISO,
	long int nelem,
	long int ipLo,
	long int ipHi,
	double temp,
	long int Collider );

/**CS_l_mixing_PS64 Collision treatment based on Pengelly and Seaton 1964
\param nelem, the chemical element, 1 for He
\param tau,
\param target_charge,
\param n,
\param l,
\param gHi,
\param Collider
*/
double CS_l_mixing_PS64(
	long int nelem,
	double tau,
	double target_charge,
	long int n,
	long int l,
	double gHi,
	long int Collider);

/**CS_l_mixing_VF01 Collision treatment based on Vrinceanu and Flannery 2001
\param ipISO
\param nelem
\param n
\param l
\param lp
\param s
\param temp
\param Collider
*/
double CS_l_mixing_VF01(
	long int ipISO,
	long int nelem,
	long int n,
	long int l,
	long int lp,
	long int s,
	double temp,
	long int Collider );


#endif /* HELIKE_CS_H_ */
