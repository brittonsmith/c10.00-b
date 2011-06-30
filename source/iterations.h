/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef ITERATIONS_H_
#define ITERATIONS_H_

/**IterStart, set and save values of many variables at start of iteration */
void IterStart(void);

/**IterRestart, reset values of many variables at start of iteration */
void IterRestart(void);

/** close out this iteration */
void IterEnd(void);

/**iter_end_check called by Cloudy after each zone to determine whether iteration is complete
 * returns true if iteration is complete, false if not */
int iter_end_check(void);

EXTERN struct t_iterations {

	/** these are the variables that control how many iterations are
	 * to be done, and number of the current iteration
	 * itermx is number of iterations to perform, set with iterate command
	 * upper limit is parameter variable ItrDim */
	long int itermx;

	/** amount of space that has been allocated for max iterations */
	long int iter_malloc;

	/**number of zones to print on each iteration*/
	long int *IterPrnt/**[ITR DIM]*/;

	/** this is false on any but the last iteration
	 * set true in startr if iter > itermx */
	bool lgLastIt;

	/** flag saying that another iteration is needed */
	bool lgIterAgain;

	/** next three implement set coverage command to limit iterations and zones */
	bool lgConverge_set;
	long int lim_zone;
	long int lim_iter;

	}	iterations;

#endif /* ITERATIONS_H_ */
