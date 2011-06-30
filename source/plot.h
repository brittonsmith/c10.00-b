/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef PLOT_H_
#define PLOT_H_

/**plot master routine to generate some sort of plot 
\param *chCall
*/
void plot(const char *chCall);

/** limit to the number of plots we can store at one time */
#define	NDPLOT	10L

EXTERN struct t_plotCom {

	/** which type of plot */
	char chPType[NDPLOT][5];

	/**lgPlotON is flag set when plot turned on */
	bool lgPlotON;

	realnum pltxmn[NDPLOT], 
	  pltxmx[NDPLOT];

	/** number of plots */
	long int nplot;

	/** flag set with plot trace command  */
	bool lgPltTrace[NDPLOT];
	}	plotCom;


#endif /* PLOT_H_ */
