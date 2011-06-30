/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef ATOMS_H_
#define ATOMS_H_

 /**
  AtomSeqBeryllium compute level populations and emissivity for Be-sequence ions
  \param cs12 
  \param cs13 
  \param cs23 
  \param t 
  \param a30
 */ 
void AtomSeqBeryllium(double cs12, 
  double cs13, 
  double cs23, 
  transition * t, 
  double a30);

 /**
 AtomSeqBoron compute cooling from 5-level boron sequence model atom
  \param  t21 
  \param  t31
  \param  t41 
  \param  t32 
  \param  t42 
  \param  t52 
  \param  cs51
  \param  cs43
  \param  cs53
  \param  cd54
  \param  pump_rate pump rate due to UV permitted lines
  \param chLabel string used to identify calling program in case of error 

 */ 
void AtomSeqBoron(
  transition * t21, 
  transition * t31, 
  transition * t41, 
  transition * t32, 
  transition * t42, 
  transition * t52, 
  double cs51,
  double cs43,
  double cs53,
  double cd54,
  double pump_rate ,
  const char *chLabel
);

/**atom_level2 do level population and cooling for two level atom 
 \param t
*/
void atom_level2(transition * t );

 /**
  atom_level3 compute three level atom, 10, 21, and 20 are line 
  \param t10 
  \param t21 
  \param t20
 */ 
void atom_level3(transition * t10, 
  transition * t21, 
  transition * t20);

 /**
  atom_pop2 do level population for simple two level atom, no radiative transfer
  \param  omega 
  \param  g1 
  \param  g2 
  \param  a21 
  \param  bltz 
  \param  abund
 */ 
double atom_pop2(double omega, 
  double g1, 
  double g2, 
  double a21, 
  double bltz, 
  double abund);

/**
atom_pop3 return value is population for 3-level atom, cm^-3 
\param  g1 statictical weights of level 1
\param  g2 statictical weights of level 2
\param  g3 statictical weights of level 3
\param  o12 collision strengths between three levels
\param  o13 collision strengths between three levels
\param  o23 collision strengths between three levels
\param  a21 transition probabilities between three levels 
\param  a31 transition probabilities between three levels
\param  a32 transition probabilities between three levels
\param  Tex12 excitation energy in Kelvin 
\param  Tex23 excitation energy in Kelvin 
\param  *pop2 returned population of level 2, cm^-3 
\param  abund incoming total abundance of ion 
\param  gam2 possible photodestruction of level 2, normally 0
\param  r12 excitation rates (s-1) due to "other" processes
\param  r13 excitation rates (s-1) due to "other" processes
 */ 
double atom_pop3(
	double g1, double g2, double g3, 
	double o12, double o13, double o23, 
	double a21, double a31, double a32, 
	double Tex12, double Tex23, 
	realnum *pop2, 
	double abund, 
	double gam2,
	double r12, 
	double r13 );

/**
atom_pop5 do populations and cooling for five level atom
\param  g[] 
\param  ex[] 
\param  cs12 
\param  cs13 
\param  cs14 
\param  cs15 
\param  cs23
\param  cs24 
\param  cs25 
\param  cs34 
\param  cs35 
\param  cs45 
\param  a21 
\param  a31 
\param  a41 
\param  a51 
\param  a32 
\param  a42 
\param  a52 
\param  a43 
\param  a53 
\param  a54 
\param  p[] 
\param  abund
\param cooling
\param cooling derivative
\param pump12
\param pump13
\param pump14
\param pump15
*/ 
void atom_pop5(const double g[], 
			   const double ex[], 
			   double cs12, 
			   double cs13, 
			   double cs14, 
			   double cs15, 
			   double cs23, 
			   double cs24, 
			   double cs25, 
			   double cs34, 
			   double cs35, 
			   double cs45, 
			   double a21, 
			   double a31, 
			   double a41, 
			   double a51, 
			   double a32, 
			   double a42, 
			   double a52, 
			   double a43, 
			   double a53, 
			   double a54, 
			   double p[], 
			   realnum abund,
			   double *Cooling, 
			   double *CoolingDeriv,
			   double pump12,
			   double pump13,
			   double pump14,
			   double pump15
			   );

/**
atom_levelN - compute populations of arbitrary n-level atom
\param nlev nlev is the number of levels to compute
\param abund ABUND is total abundance of species, used for nth equation
\param g[] G(ndim) is stat weight of levels
\param ex[] EX(ndim) is excitation potential of levels, either wn or deg K 0 for first one, NOT d(ENER), but energy rel to ground 
\param chExUnits this is 'K' for above ex[] as Kelvin deg, is 'w' for wavenumbers
\param pops[] populations of each level as deduced here
\param depart[] departure coefficient derived here
\param AulEscp net transition rate, A * esc prob, s-1
\param col_str col str rom up to low
\param AulDest AulDest(ilo,ihi) is destruction rate, from up to low, A * dest prob, [s-1], asserts confirm that ihi,lo is zero
\param AulPump AulPump(lo, hi) is pumping rate, A * occ num, (hi,lo) must be zero, [s-1]  
\param CollRate collision rates, evaluated here and returned for cooling by calling function, unless following flag is true.  
	If true then calling function has already filled in these rates.  CollRate[i][j] is rate from i to j 
\param create this is an additional creation rate, normally zero, units cm-3 s-1
\param destroy[] this is an additional destruction rate to continuum, normally zero, units s-1 
\param lgCollRateDone flag saying whether CollRate already done, or we need to do it here
\param cooltl total cooling, set here but nothing done with it
\param coolder derivative of cooling, set here but nothing done with it
\param chLabel string used to identify calling program in case of error
\param lgNegPop lgNegPop flag indicating what we have done positive if negative populations occurred zero if normal calculation done negative if too cold (for some atoms other routine will be called in this case)
\param lgZeroPop true if populations are zero, either due to zero abundance of very low temperature
\param lgDeBug option to print matrices for debugging
\post atoms.PopLevels[n], atoms.DepLTELevels[n] are set lines added to ots array
*/ 
void atom_levelN(
	long int nlev, 
	realnum abund, 
	const double g[], 
	const double ex[], 
	char chExUnits,
	double pops[], 
	double depart[],
	double ***AulEscp, 
	double ***col_str, 
	double ***AulDest, 
	double ***AulPump, 
	double ***CollRate,
	const double create[] ,
	const double destroy[] ,
	bool lgCollRateDone,
	double *cooltl, 
	double *coolder, 
	const char *chLabel, 
	int *nNegPop,
	bool *lgZeroPop ,
	bool lgDeBug);

/**atom_oi drive the solution of OI level populations, Ly-beta pumping 
\param coloi
*/
void atom_oi_calc(double *coloi);

/** number of levels in OI atom */
const int N_OI_LEVELS = 6;
const long LIMLEVELN = 20L;

EXTERN struct t_atoms {

	/** photo excitation of excited state of NI */
	realnum p2nit, d5200r;

	/** collisional rates */
	double c12, c13;

	/** array indices for level energies for OI bowen problem,
	 * used to generate Boltzmann factors */
	long int ipoiex[5];

	/** number of negative OI level populations in current calculation */
	long int nNegOI;

	/** populations from OI fluorescence problem */
	double popoi[N_OI_LEVELS];

	realnum pmpo51, pmpo15;

	/** excited state of Mg+ */
	realnum xMg2Max,
		/** its population */
		popmg2;

	/**
	 * this stores most recently evaluated level populations
	 * resulting from the leven family of routines
	 * PopLevels is population (cm^-3) of the levels
	 * DepLevels is lte departure coef
	 */
	double PopLevels[LIMLEVELN+1], 
	  DepLTELevels[LIMLEVELN+1];

	}	atoms;

#endif /* ATOMS_H_ */
