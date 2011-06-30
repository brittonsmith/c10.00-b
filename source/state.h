/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef STATE_H_
#define STATE_H_


/**state_get_put get or save state - job is either "get" or "put" */
void state_get_put( const char chJob[] );

EXTERN struct t_state {

	/** file pointer for the files to get or put the state 
	 * set in parse_state */
	char chPutFilename[INPUT_LINE_LENGTH] , 
		chGetFilename[INPUT_LINE_LENGTH];

	/** flags saying whether to get or put state 
	 * set in parse_state */
	bool lgGet_state , 
		lgPut_state;

	/** option to put all iterations, include ALL on state put */
	bool lgPutAll;

	/** print keyword on state command to turn on printout */
	bool lgState_print;

	}	state;


#endif /* STATE_H_ */
