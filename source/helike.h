/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef HELIKE_H_
#define HELIKE_H_


/** \verbatim
  Set this flag to one of the following values
 	0	don't print
 	1	print As
 	2	print only forbidden As
 	3	print Es
 	4	print threshold photoionization cross-sections
 	5	print radiative recombination coefficients.	
        6   print photoionization cross-section grids. \endverbatim <BR>
  Out files are "Helike___.txt" with the blank being, respectively,
 *	"As","fAs","Es","ThPCS","RR", and "PCSgrid."	*/
#define	printflag	(0)

/** the magic number for the table of He-like transition probabilities, YYMMDD */
#define TRANSPROBMAGIC	(60725)
/** the magic number for the table of He1 effective collision strengths, YYMMDD */
#define COLLISMAGIC		(30915)
/** the magic number for the table of He1 case A and B emissivities, YYMMDD */
#define CASEABMAGIC		(51214)

/**HeCollidSetup read in helium collision data files */
void HeCollidSetup( void );

/**helike_energy get helike level energy in cm-1 
\param nelem
\param ipLev 
*/
double helike_energy(long nelem, long ipLev );

/**helike_quantum_defect get quantum defect for helium-like level
\param nelem
\param ipLev
*/
double helike_quantum_defect( long int nelem, long int ipLev );

/**helike_transprob get transition probability for helium-like transition [s-1]
\param nelem
\param ipHi
\param ipLo
*/
realnum helike_transprob( long nelem, long ipHi, long ipLo );

/**AGN_He1_CS routine to save table needed for AGN3 - collision strengths of HeI 
\param *ioPun
*/
void AGN_He1_CS( FILE *ioPun );

#endif /* HELIKE_H_ */
