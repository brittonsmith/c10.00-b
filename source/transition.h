/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland 
and
 * others.  For conditions of distribution and use see copyright notice in licen
se.txt */

#ifndef TRANSITION_H_
#define TRANSITION_H_

/*Generalized structure used to hold the transition information for both atoms/ions and molecules*/
class transition
{
 public:
	transition()
	{
		Junk();
	}

	quantumState *Lo, *Hi;
	emission *Emis;
	collision Coll;

	/* Possible linked-list optimization */
	/*transition *nextemis, *nextcoll; */

	/** wavelentgh, usually in Angstroms, used for printout, can be any units */
	realnum WLAng;

	/** transition energy in degrees kelvin*/
	realnum EnergyK;

	/** transition energy in ergs */
	realnum EnergyErg;

	/** transition energy in wavenumbers */
	realnum EnergyWN;

	/** index for line within continuum array,
	 * this is on the f, not c, scale,
	 * negative ipCont means this is not a radiative transition, 
	 * and is used as a sentinel */
	long ipCont;

   /**TransitionJunk set all elements of transition struc to dangerous values  
		\param *t
	*/
	void Junk( void );

	/**TransitionZero set all elements of transition struc to zero 
		\param *t
	*/
	void Zero( void );

   /**outline - adds line photons to reflin and outlin */
	void outline( double nonScatteredFraction, bool lgDoChecks );

   /**outline_resonance - adds line photons to reflin and outlin,
	   setting nonScatteredFraction as default for resonance lines */
	void outline_resonance(  );

};

/** enter lines into the line storage array, called once per zone for each line
\param xInten xInten - local emissivity per unit vol, no fill fac
\param wavelength lam integer wavelength
\param *chLab string label for ion
\param chInfo character type of entry for line - 'c' cooling, 'h' heating, 'i' info only, 'r' recom line
\param *chComment string explaining line 
*/

/**PutLine enter local line intensity into the intensity stack for eventual printout 
\param *t transition structure for line
\param *chComment a description of the line
*/
void PutLine(const transition *t, const char *chComment);

/**PutLine enter local line intensity into the intensity stack for eventual printout 
\param *t transition structure for line
\param *chComment a description of the line
\param *chLabel the line label
*/
void PutLine(const transition *t, const char *chComment, const char *chLabel);

/**TexcLine derive excitation temperature of line from contents of line array 
\param *t
*/
double TexcLine(const transition *t);

/**DumpLine print various information about an emission line vector, used in debugging 
\param *t
*/
void DumpLine(const transition *t);

/** returns fraction of populations the produce emission 
\param *t
*/
double emit_frac(const transition *t);

/** generate null terminated line label from contents of line trans array 
\param *t
*/
void chIonLbl(char*, const transition *t);

/**chLineLbl use information in line transfer arrays to generate a line label<BR>
 this label is null terminated 
 \param *t
 */
char* chLineLbl(const transition *t);

/**PutCS enter a collision strength into an individual line struc 
\param cs
\param *t  the line struc 
*/
void PutCS(double cs, 
  transition * t);

/**GenerateTransitionConfiguration - given transition *t, writes a label
 * t->Lo->chConfig - t->Hi->chConfig (i.e., 2^3S - 2^3P)
 \param t
 \char *chComment
 */
void GenerateTransitionConfiguration( const transition *t, char *chComment );

/**OccupationNumberLine - derive the photon occupation number at line center for any line 
\param *t
*/
double OccupationNumberLine(const transition *t);

/**PutExtra enter and 'extra' intensity source for some line 
\param Extra
*/
void PutExtra(double Extra);

/** convert down coll rate back into electron cs in case other parts of code need this for reference 
\param *t - line struct collision strength is stored in t->cs 
\param rate - deexcitation rate, units s-1 
*/
void LineConvRate2CS( transition * t , realnum rate );

/**lgTauGood returns true is we have good (positive) outward optical depths
 * not true if we have overrun optical depth scale from previous iteration
\param *t
*/
inline bool lgTauGood( transition* t )
{
	// first iteration only use inward optical depths so scale good
	return ( iteration == 1 || 
		// maser - optical depths also ok (but bizarre) 
		t->Emis->TauIn <= 0. || 
		// TauIn < TauTot means outward optical depth is positive, so OK
		t->Emis->TauIn < t->Emis->TauTot );
}

/**MakeCS compute collision strength by g-bar approximations 
\param *t
*/
void MakeCS(transition * t );

/** \todo 2 bring these two together. */
/**AddLine2Stack add generic emission line to GenericLines and return pointer to that state. */
emission *AddLine2Stack( bool lgRadiativeTrans  );
emission *AddLine2Stack( realnum Aul, transition *trans );

#endif // _TRANSITION_H_
