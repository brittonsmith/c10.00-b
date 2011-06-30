/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef LINES_SERVICE_H_
#define LINES_SERVICE_H_

void linadd(
  double xInten,
  realnum wavelength,
  const char *chLab,
  char chInfo ,	
  const char *chComment );

/*outline_base - adds line photons to reflin and outlin */
void outline_base(double dampXvel, double damp, bool lgTransStackLine, long int ip, double phots, realnum inwd,
						double nonScatteredFraction);

/*outline_base_bin - adds line photons to bins of reflin and outlin */
void outline_base_bin(bool lgTransStackLine, long int ip, double phots, realnum inwd,
						double nonScatteredFraction);

/** put forbidden line into stack, using index derived below 
\param xInten - local emissivity per unit vol
\param wavelength wavelength Angstroms
\param *chLab string label for ion
\param ipnt offset of line in continuum mesh
\param chInfo character type of entry for line - 'c' cooling, 'h' heating, 'i' info only, 'r' recom line
\param lgOutToo should line be included in outward beam?
\param *chComment string explaining line 
*/
void lindst(double xInten,
  realnum wavelength, 
  const char *chLab, 
  long int ipnt, 
  char chInfo, 
  bool lgOutToo,
  const char *chComment);

/** put forbidden line into stack, using index derived below
\param dampXvel - damping constant times Doppler velocity
\param damp - damping constant
\param xInten - local emissivity per unit vol
\param wavelength wavelength Angstroms
\param *chLab string label for ion
\param ipnt offset of line in continuum mesh
\param chInfo character type of entry for line - 'c' cooling, 'h' heating, 'i' info only, 'r' recom line
\param lgOutToo should line be included in outward beam?
\param *chComment string explaining line
*/
void lindst(double dampXvel,
  double damp,
  double xInten,
  realnum wavelength,
  const char *chLab,
  long int ipnt,
  char chInfo,
  bool lgOutToo,
  const char *chComment);

class transition;

/** put forbidden line into stack, using index derived below 
\param t -- pointer to transition data
\param wavelength wavelength Angstroms
\param *chLab string label for ion
\param ipnt offset of line in continuum mesh
\param chInfo character type of entry for line - 'c' cooling, 'h' heating, 'i' info only, 'r' recom line
\param lgOutToo should line be included in outward beam?
\param *chComment string explaining line 
*/
void lindst(
  transition *t,
  const char *chLab, 
  char chInfo, 
  bool lgOutToo,
  const char *chComment);

/** absorption due to continuous opacity 
\param emissivity [erg cm-3 s-1] in inward direction 
\param emissivity [erg cm-3 s-1] in outward direction 
\param array index for continuum frequency 
*/
double emergent_line( 
	 /* emissivity [erg cm-3 s-1] in inward direction */
	 double emissivity_in , 
	 /* emissivity [erg cm-3 s-1] in outward direction */
	 double emissivity_out , 
	 /* array index for continuum frequency */
	 long int ipCont );

/**PntForLine generate pointer for forbidden line 
\param wavelength wavelength of line in Angstroms 
\param *chLabel label for the line
\param *ipnt this is array index on the f, not c scale,
   for the continuum cell holding the line
*/
void PntForLine(double wavelength, 
  const char *chLabel, 
  long int *ipnt);

/**GetGF convert Einstein A into oscillator strength 
\param eina
\param enercm
\param gup
*/
double GetGF(double eina, 
	  double enercm, 
	  double gup);

/**eina convert a gf into an Einstein A 
\param gf
\param enercm
\param gup
*/
double eina(double gf, 
	  double enercm, 
	  double gup);

/**abscf convert gf into absorption coefficient 
\param gf
\param enercm
\param gl
*/
double abscf(double gf, 
	  double enercm, 
	  double gl);

/** setting true will use low-density Lyman branching ratios */
#define LOWDEN_LYMAN 0

/**RefIndex calculates the index of refraction of air using the line energy in wavenumbers,
 * used to convert vacuum wavelengths to air wavelengths. 
 \param EnergyWN
 */
double RefIndex(double EnergyWN);


/**WavlenErrorGet - given the real wavelength in A for a line
 * routine will find the error expected between the real 
 * wavelength and the wavelength printed in the output, with 4 sig figs,
  \param wavelength
  \return function returns difference between exact and 4 sig fig wl, so 
  we have found correct line is fabs(d wl) < return 
  */
realnum WavlenErrorGet( realnum wavelength );

/**EmLineJunk set all elements of emission struc to dangerous values 
\param *t
*/
void EmLineJunk( emission * t );

/**CollisionJunk set all elements of emission struc to dangerous values  
\param *t
*/
void CollisionJunk( collision * t );

/**StateJunk set all elements of quantumState struc to dangerous values  
\param *t
*/
void StateJunk( quantumState * t );

/**EmLineZero set all elements of emission struc to zero 
\param *t
*/
void EmLineZero( emission * t );

/**CollisionZero set all elements of collision struc to zero 
\param *t
*/
void CollisionZero( collision * t );

/**StateZero set all elements of quantumState struc to zero 
\param *t
*/
void StateZero( quantumState * t );

/** convert down coll rate back into electron cs in case other parts of code need this for reference 
\param gHi - stat weight of upper level 
\param rate - deexcitation rate, units s-1 
*/
double ConvRate2CS( realnum gHi , realnum rate );

/** convert collisional deexcitation cross section for into collision strength 
\param CrsSectCM2 - the cross section
\param gLo - statistical weight of lower level of transition
\param E_ProjectileRyd - initial projectile energy in Rydbergs
\param reduced_mass_grams - reduced mass MpMt/(Mp+Mt) of projectile-target system
*/
double ConvCrossSect2CollStr( double CrsSectCM2, double gLo, double E_ProjectileRyd, double reduced_mass_grams );

/**totlin sum total intensity of cooling, recombination, or intensity lines 
\param chInfo chInfor is 1 char,<BR>
	'i' information, <BR>
	'r' recombination or <BR>
	'c' collision
*/
double totlin(
	int chInfo);


/**FndLineHt search through line heat arrays to find the strongest heat source 
\param *level
\param *ipStrong  this is the index of the strongest line in the array on the c scale 
\param *Strong
*/
void FndLineHt(long int *level, 
  /**/
  long int *ipStrong, 
  double *Strong);

/**AddState2Stack add generic quantum state to GenericStates and return pointer to that state. */
quantumState *AddState2Stack( void );

#endif /* LINES_SERVICE_H_ */
