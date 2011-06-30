/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef LINES_H_
#define LINES_H_


/**lines main routine to put emission line intensities into line stack */
void lines(void);

/** general information at start of lines */
void lines_general(void);

/** the hydrogenic iso-sequence */
void lines_hydro(void);

/** information about grains */
void lines_grains(void);

/**lines_setup convert level 1 and level 2 line parameters and pointers 
 * into internal form used by code */
void lines_setup(void);

/** enter all continua */
void lines_continuum(void);

/** enter all molecules into emission line stack */
void lines_molecules(void);

/** enter all helium iso seq into emission line stack */
void lines_helium(void);

/**lines_lv1_li_ne place lines of elements lithium through neon into lines storage stack */
void lines_lv1_li_ne(void);

/**lines_lv1_na_ar place lines of elements sodium through argon into lines storage stack */
void lines_lv1_na_ar(void);

/**lines_lv1_k_zn place lines of elements potassium and heavier into lines storage stack */
void lines_lv1_k_zn(void);

/** routine to stuff comments into the stack of comments,
 * return is index to this comment */
long int StuffComment( const char * chComment );

/** lines_table invoked by table lines command, check if we can find all lines in a given list
 * returns 0 if ok, n is n lines not found */
int lines_table(void);

#define NHOLDCOMMENTS 100

/** this struc is different from following since they are only pointer here, will be malloced 
 * to form a large array after number of lines is counted, but this is the final form */
EXTERN struct t_LineSave {
	/** number of emission lines in main stack */

	/** nsum is current index, to last line entered in the array */
	/** linadd increments nsum before doing anything else */
	long int nsum;

	/** nsumAllocated is number of lines allocated in emission line stack
	 * must not change between iterations or grid points */
	long int nsumAllocated;

	/** index to number of comments printed within the block of lines */
	long int nComment;

	/** there are three types of calls to lines() 
	 * ipass = -1, first call, only count number off lines
	 * ipass =  0, second pass, save labels and wavelengths
	 * ipass =  1, integrate intensity*/
	long int ipass;

	/** holds comment strings associated with various blocks of output lines */
	char chHoldComments[NHOLDCOMMENTS][INPUT_LINE_LENGTH];

	/** NormWL is array index for emission line on normalize command */
	long int ipNormWavL;

	/** WavLNorm is wavelength of emission line on normalize command */
	realnum WavLNorm;

	/** the associated uncertainty in the wavelength */
	realnum errorwave;

	/** number of significant figures for lines
	 * this affects all aspects of reading and writing lines */
	long int sig_figs;

	/** ScaleNormLine is the scale factor for its appearance */
	double ScaleNormLine;

	/** chNormLab is optional label */
	char chNormLab[5];

	/** flag saying whether norm has been set */
	bool lgNormSet;

	/** save rec coefficient data for recombination lines of C, N, O */
	realnum RecCoefCNO[4][471];

}	LineSave;

/** this struc is different from above since only pointer here, will be malloced 
 * to form a large array after number of lines is counted.<BR>
 * these are the main line save arrays */
typedef struct t_tag_LineSv {

	/** one char saying whether heat 'h', cooling 'c', information, 'i' */
	char chSumTyp;

	/** the four char string label for the line */
	char chALab[5];

	/** integrated intensity of the line, 
	 * [0] is intrinsic, 
	 * [1] emergent 
	 * [2] is intrinsic, 
	 * [3] emergent 
	 */
	double SumLine[4];

	/** the emissivity, per unit vol, for current conditions, */
	double emslin[2];

	/** the wavelength of the line */
	realnum wavelength;

	/** comment describing the line */
	const char *chComment;

} LinSv;

extern LinSv *LineSv, *LineSvSortWL;


#endif /* LINES_H_ */
