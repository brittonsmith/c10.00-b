/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef MOLE_CO_PRIV_H_
#define MOLE_CO_PRIV_H_

/* mole_co_priv.h */

#include "hash.h"

extern struct mole_priv_s {
	hashtab *spectab, *reactab, *elemtab;
} mole_priv;

EXTERN struct t_coreactions {
	struct COmole_rate_s **list;
	long int n;
} coreactions;

#define MAXREACTANTS 3
#define MAXPRODUCTS  4

/* Structure containing reaction data */
struct COmole_rate_s {
	int index;
	char *label;
	int nreactants, nrates, nproducts,photon;
	struct molecule *reactants[MAXREACTANTS];
	struct molecule *rate_species[MAXREACTANTS];
	struct molecule *products[MAXPRODUCTS];
	double rk, reduced_mass, a, b, c;
	double (*fun)(struct COmole_rate_s *rate);
};

enum {CHARS_ELEMENT=6};
extern struct chem_element_s {
	int ipCl; /* Index of element in external arrays */
	int
	  ipMl,   /* Index of atomic species in molecule arrays */
		ipMlP,
		ipZ;  /* Index of + ions in molecule arrays */
	char chName[CHARS_ELEMENT]; /* Chemical symbols for elements */
} **chem_element;

extern int32 *ipiv;
extern realnum *tot_ion;
/** CO_step
*/
extern void CO_step(void);
/**CO_solve fills in matrix for heavy elements molecular routines 
\param *lgNegPop set true if we found neg pops
\param *lgZerPop set true if we tried to compute the pops, but some were zero 
*/
extern void CO_solve(
	bool *lgNegPop, 
	bool *lgZerPop );

#endif /* MOLE_CO_PRIV_H_ */
