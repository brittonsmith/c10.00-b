/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef WARNINGS_H_
#define WARNINGS_H_

/* warnings.h */

/**wcnint initialize stack or warnings, cautions, notes */
void wcnint(void);

/**warnin enter warnings at the end of the calculations into large stack 
\param *chLine
*/
void warnin(char *chLine);

/**notein enter a note about calculation into comment array 
\param *chLine
*/
void notein(char *chLine);

/**bangin called by routine comment to enter surprise into comment stack 
\param *chLine
*/
void bangin(char *chLine);

/**caunin called by comment to enter caution into comment stack 
\param *chLine
*/
void caunin(char *chLine);

/** this is the limit to the number of warnings, cautions, and
 * notes that we can save */
#define	LIMWCN	2000

EXTERN struct t_warnings {

	/** this are counters for the number of warnings,
	 * cautions, notes and surprises in the calculation*/
	long int nwarn, 
	  ncaun, 
	  nnote, 
	  nbang;

	/** a comment about the geometry after model stops */
	char chRgcln[2][INPUT_LINE_LENGTH];

	/** these are the strings that contain the warnings, cautions,
	 * and notes about the calculation */
	char chWarnln[LIMWCN][INPUT_LINE_LENGTH], 
	  chCaunln[LIMWCN][INPUT_LINE_LENGTH], 
	  chNoteln[LIMWCN][INPUT_LINE_LENGTH], 
	  chBangln[LIMWCN][INPUT_LINE_LENGTH];

	/** flags set if warnings or cautions present */
	bool lgWarngs, 
	  lgCautns;

	}	warnings;


#endif /* WARNINGS_H_ */
