/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef INPUT_H_
#define INPUT_H_

/* input.h */

/** limit to number of line images that can be read in */
#define	NKRD	4000

/** lgInputComment - parse comment - check if argument is comment string, 
 * either upper or lower case -
 * returns true if line is a comment, false if not 
 * a comment is any line starting with "C ", *, %, //, or # 
 \param *chLine the input line string
 */
bool lgInputComment( 
  const char *chLine );

/** input_readvector: read n numbers from the file chFile and store them in vector[] */
void input_readvector(const char* chFile, /**< file name to read from */
		      double vector[],    /**< vector[n] - the numbers that were read from the input line(s) */
		      long n,             /**< number of elements in vector[] that we need to read */
		      bool* lgEOF);       /**< was EOF reached before enough numbers were read? */

EXTERN struct t_input {

	/** 
	 * we will save the original (not caped) image of the line here 
	 */
	char chCardSav[NKRD][INPUT_LINE_LENGTH], 

	/** 
	 * title entered with the title command 
	 */
	  chTitle[INPUT_LINE_LENGTH];

	/**
	 * delimiter character for file names, / for *nix, \\ for win
	 */
	char chDelimiter[3];

	long int 
	  /** one less than the total number of lines read in with cdRead */
	  nSave,

	  /** this points to the command we are now parsing, within the stack of commands */
	  nRead,		

	  /** number of init commands saved */
	  nSaveIni,	

	  /** +/-1, says whether to increment or decrement nRead, since init commands
		* are at the bottom of the stack and we read backwards */
	  iReadWay,	

	  /** saves current value of nRead, while parsing init commands */
	  nReadSv;	

	/** this is set true if underscore present in input stream, which was
	 * set to space */
	bool lgUnderscoreFound;

	/** this is set true if left or right bracket, [ or ], present in input stream, which was
	 * set to space */
	bool lgBracketFound;

	/** set true with no buffering command, used to print comment at end */
	bool lgSetNoBuffering;

/** get the next input command off the command stack 
 * if more then copy into chCard and set lgEOF false,
 * if all command processed then set lgEOF true 
 \param *chCard the input line string
 \param *lgEOF true if hit end of file
 */
	// friend CodeSmell
	private:
	friend class Parser;
	void readarray(
		char *chCard, 
		bool *lgEOF);

	public:
	void echo( FILE *ipOUT);

	/** called when 'init' command hit, to reset counters for
	 * placing line images within the storage array */
	void init(void);


	}	input;



#endif /* INPUT_H_ */
