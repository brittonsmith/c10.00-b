/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef PRT_H_
#define PRT_H_

/**PrtZone print out individual zone results */
void PrtZone(void);

/**PrtContinuum print information about continuum if requested with PRINT CONTINUUM command */
void PrtContinuum(void);

/**PrtComment analyze model, generating comments on its features */
void PrtComment(void);

/**PrtFinal create final pages of printout, emission line intensities, etc */
void PrtFinal(void);

/**prt_wl write wavelength to io 
\param *io
\param wavelength
*/
void prt_wl( 
	FILE *io , 
	realnum wavelength );

/**sprt_wl write wavelength to string - must be kept parallel with prt_wl 
\param *chString
\param wl
*/
void sprt_wl( 
	char *chString , 
	realnum wl );

/**PrtHeader print large block of incident continuum numbers at start, 
 just after echoing input commands */
void PrtHeader(void);

/**prt_LineLabels save all labels and wavelengths for emission line array 
\param io file handle to write output
\param lgPrintAll print all if true, if false then do not print parts of 
 transferred lines
*/
void prt_LineLabels(
	FILE * io,
	bool lgPrintAll
	);

/**prtmet print all line optical depths at end of iteration */
void prtmet(void);

/**prme print heavy element line optical depths at end of calculation 
\param chDoIt
\param t
*/
void prme(
  const char *chDoIt, 
  transition * t);

/**PrtMeanIon print mean ionization fractions for all elements,
 * output will go to stream pointed to by argument  
 * chTyp is either 'i' or 't' for mean ionization or temperature 
 \param chType
 \param lgDensity true include density, false do not
 */
void PrtMeanIon( char chType , 
			bool lgDensity,
			FILE *);

/**PrtLineSum parse print line sum command to enter set of lines into sum  
\param chDo the job to do, either " SUM" or "READ"
*/
double PrtLineSum(void);

/**PrtLinePres print line radiation pressures for current conditions 
 * output goes top openned file handle */
void PrtLinePres(FILE *ioPRESSURE);

/**PrtColumns print column densities of all elements 
\param ioMEAN this is stream used for io, is stdout when called by final,
       is save unit when save output generated
*/
void PrtColumns(
	 FILE *ioMEAN );

/** DatabasePrintReference print some database references */
void DatabasePrintReference( void );

/**PrtAllTau master routine controlling printout of optical depths at
 end of calculation */
void PrtAllTau(void);

EXTERN struct t_prt {

	/** lgSortLines, option to sort lines by wavelength- print sort
	 command */
	bool lgSortLines;

	/** if above is set, then one of the following must also be set,
	 * say whether to sort by wavelength or intensity */
	bool lgSortLineWavelength , lgSortLineIntensity;

	/** lower and upper wavelength bounds for printed spectrum,
	 * range option on print sort command */
	realnum wlSort1 , wlSort2;

	/** print hydrogenic level populations, 
	 * set with print hydrogenic command
	bool lgPrintHLevPops; */

	/** print column densities */
	bool lgPrintColumns;

	/** should we print execution time?  normally true, but set false
	 * with no times command so that different runs can compare exactly */
	bool lgPrintTime;

	/** print ages command tells code to print various timescales */
	bool lgPrnAges;

	/** option to print maser lines (true) normally false
	 * print maser turns on */
	bool lgPrtMaser;

	/** lgPrtTau tells whether to print line optical depths */
	bool lgPrtTau;

	/** lgPrintFluxEarth says to print flux of lines at Earth, 
	 * if luminosity can be predicted */
	bool lgPrintFluxEarth;

	/** print line surface brightness command, units either sr or sq arcsec,
	 * default is SR, set to arcsec with arcsec option */
	bool lgSurfaceBrightness , lgSurfaceBrightness_SR;

	/** PrtTauFnt is smallest line optical depth to print */
	realnum PrtTauFnt;

	/** these are various contributors to the line output,
	 * and are changed with the print line or print continuum commands
	 * in prtfinal code uses these to make a final filter over what lines
	 * will be printed */
	bool lgPrnPump, 
	  lgPrnHeat, 
	  lgPrnColl, 
	  lgPrnInwd, 
	  lgPrnDiff;

	/** print predictions from collapsed levels of iso sequences,
	 * print line iso collapsed */
	bool lgPrnIsoCollapsed;

	/** flag set with print continuum command */
	bool lgPrtCont;

	/* flag set with print continuum index command, to identify all lines
	 * that lie within a continuum cell */
	bool lgPrtContIndices;
	/* these are lower and upper limits to the energy range in Rydbergs.
	 * they are the first and second number on the command line, lower and
	 * upper bounds of the code are used if not specified */
	realnum lgPrtContIndices_lo_E , 
		lgPrtContIndices_hi_E;

	/** flags for determining what is included in nFnu */
	bool lgSourceReflected;
	bool lgSourceTransmitted;
	bool lgDiffuseInward;
	bool lgDiffuseOutward;

	/** flag set with print departure coefficients */
	bool lgPrtBN;

	/** if true then print only last iteration */
	bool lgPrtLastIt;

	/** flag set with print short command */
	bool lgPrtShort;

	/** lgOnlyZone set with print only zones */
	bool lgOnlyZone;
	/** lgOnlyHead set with print only header */
	bool lgOnlyHead;

	/** lgPrtStart is option to start printout at certain zone */
	bool lgPrtStart;

	/**nstart is zone number, set with print start command */
	long int nstart;

	/** flag to turn on printout of heat sources */
	bool lgPrintHeating;

	/** flag set with print array command to print ionization recombination arrays */
	bool lgPrtArry[LIMELM];

	/** logical lgFaintOn normally true, says to not print very faint lines
	  set false with print faint off command
	 lines fainter than TooFaint will not be printed.  This is set in 
	 zerologic and reset with print line faint command  */
	realnum TooFaint;
	bool lgFaintOn;

	/** flag set true if print faint command entered,
	 * only used to not override it with print short */
	bool lgFntSet;

	/** these implement the print line cell commmand, 
	 * flag saying to do this */
	bool lgPrnLineCell;
	/** the cell number, on the physics scale, counts from 1, 
	 * will print labels of all lines that lie within that cell */
	long int nPrnLineCell;

	/** option to print main block of lines as a single column
	 * instead of the normal array.  if true then usual array */
	bool lgPrtLineArray;

	/** printing as a column also has an option to print linear quantity 
	 * in exponential format */
	bool lgPrtLineLog;

	/** flag set by print line cumulative command, also print large set of
	 * emission line integrated intensities over time depend model */
	bool lgPrintLineCumulative;

	/** quantities to do with radiation field and printed in header */
	realnum qx, 
	  powion, 
	  xpow, 
	  pbal, 
	  q, 
	  qgam, 
	  pradio, 
	  fx1ryd;
	long int ipeak;
	realnum GammaLumin;

	long int nzdump;

	/** print citations command flag */
	bool lgPrtCitations;

	}	prt;



#endif /* PRT_H_ */
