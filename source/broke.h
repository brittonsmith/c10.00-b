/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef BROKE_H_
#define BROKE_H_

/** broke.h */
EXTERN struct t_broke {
	/** logical flag saying that the code is broken, set by calling broken(); 
	 * causes a warning to be printed at the end of the calculation.  prototype
	 * is in cddefines.h */
	bool lgBroke;
	/** flag set with call to fixit, says that code needs attention, but not
	 * broken.  only causes a caution. */
	bool lgFixit;

	/** says that new code is in place that needs to be checked */
	bool lgCheckit;

	}	broke;


#endif /* BROKE_H_ */
