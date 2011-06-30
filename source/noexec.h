/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef NOEXEC_H_
#define NOEXEC_H_

/* noexec.h */

/**
 *flag saying no to compute model, just generate continuum and stop
 *set with cdNoex call
 */
EXTERN struct t_noexec {
	bool lgNoExec;
	}	noexec;



#endif /* NOEXEC_H_ */
