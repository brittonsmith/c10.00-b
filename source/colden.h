/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef COLDEN_H_
#define COLDEN_H_

/* colden.h */
/** number of column densities to remember */
#define	NCOLD	10

/** total hydrogen column density, all forms, xxx + 2*H2 */
#define	ipCOL_HTOT	0
/** H- H minus column denisty*/
#define	ipCOL_HMIN	1
/** column density of H2g */
#define	ipCOL_H2g	2
/** column density of H2s */
#define	ipCOL_H2s	3
/** H2+ */
#define	ipCOL_H2p	4
/** H0 */
#define	ipCOL_H0	5
/** HeH+ */
#define	ipCOL_HeHp	6
/** H+ */
#define	ipCOL_Hp	7
/** H3+ */
#define	ipCOL_H3p	8
/** column density in electrons */
#define	ipCOL_elec	9

EXTERN struct t_colden {

	/**save total column densities in various species for this and 
	 * previous iteration, to check whether
	 * it changed by too much (a bad sign)
	 * column densities, mostly h molecule related */
	realnum colden[NCOLD], 
		/** the previous iteration's coluumn density */
	  colden_old[NCOLD];

	/** integral of n(H2) / v(H2) */
	realnum coldenH2_ov_vel;

	/** integral of ne np over radius */
	double dlnenp;

	/** integral of ne n(He+) over radius */
	double dlnenHep;

	/** integral of ne n(He++) over radius */
	double dlnenHepp;

	/** integral of ne n(C+) over radius */
	double dlnenCp;

	/** integral of n(H0) / Tspin */
	double H0_ov_Tspin;

	/** integral of n(OH) / Tkin */
	double OH_ov_Tspin; 

	/** pops and column density for CII atom */
	realnum C2Pops[5],
		C2Colden[5];

	/** pops and column density for upper level of CIII] atom */
	realnum C3Pops[4],
		C3Colden[4];

	/** pops and column density for SiII atom */
	realnum Si2Pops[5],
		Si2Colden[5];

	/** pops and column density for CI atom */
	realnum C1Pops[3],
		C1Colden[3];

	/** pops and column density for OI atom */
	realnum O1Pops[3],
		O1Colden[3];

	/** He I 23S */
	double He123S;

	/** variables to do with mean mass per particle over model,
	 * called \<Mol\> in final print out, used to get Jeans mass */
	realnum rjnmin, 
	  ajmmin;
	realnum TotMassColl, 
	  tmas, 
	  wmas;

	/** column densities in the lower and upper level of 1s of H0 */
	double H0_21cm_upper;
	double H0_21cm_lower;

	}	colden;

#endif /* COLDEN_H_ */
