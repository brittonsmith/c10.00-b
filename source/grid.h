/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef GRID_H_
#define GRID_H_

#include "optimize.h"

/**called by cdDrive, this returns 0 if things went ok, 1 for disaster*/
bool grid_do(void);

/**gridXspec
\param xc[]
\param nInterpVars
*/
void gridXspec(realnum *, long);

/**GridGatherInCloudy */
void GridGatherInCloudy( void );

const int NUM_OUTPUT_TYPES = 11;

EXTERN struct t_grid
{
	vector<realnum> Energies;
	multi_arr<realnum,3> Spectra;
	char **paramNames;
	long *paramMethods;
	realnum **paramRange;
	realnum **paramData;
	realnum **interpParameters;

	realnum paramLimits[LIMPAR][2];
	realnum paramIncrements[LIMPAR];
	bool lgLinearIncrements[LIMPAR];

	/** set true when grid command entered */
	bool lgGrid, 
		lgGridDone,
		lgStrictRepeat;

	/** number of grid commands entered */
	long int nGridCommands;

	long nintparm,
		naddparm,
		numEnergies,
		numParamValues[LIMPAR],
		totNumModels;

	bool lgOutputTypeOn[NUM_OUTPUT_TYPES];

	long ipLoEnergy, ipHiEnergy;
	realnum LoEnergy_keV, HiEnergy_keV;

} grid;

#endif /* GRID_H_ */
