/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* HyperfineCreat establish space for hf arrays, reads atomic data from hyperfine.dat */
/* HyperfineCS - returns collision strengths for hyperfine struc transitions */
/*H21cm computes rate for H 21 cm from upper to lower excitation by atomic hydrogen */ 
/*h21_t_ge_20 compute rate for H 21 cm from upper to lower excitation by atomic hydrogen */ 
/*h21_t_lt_20 compute rate for H 21 cm from upper to lower excitation by atomic hydrogen */ 
/*H21cm_electron compute H 21 cm rate from upper to lower excitation by electrons - call by CoolEvaluate */
/*H21cm_H_atom - evaluate H atom spin changing collision rate, called by CoolEvaluate */
/*H21cm_proton - evaluate proton spin changing H atom collision rate, */
#include "cddefines.h"
#include "conv.h"
#include "lines_service.h"
#include "phycon.h"
#include "dense.h"
#include "rfield.h"
#include "taulines.h"
#include "iso.h"
#include "trace.h"
#include "hyperfine.h"
#include "physconst.h"  

/* H21_cm_pops - fine level populations for 21 cm with Lya pumping included 
 * called in CoolEvaluate */
void H21_cm_pops( void )
{
	/*atom_level2( HFLines[0] );*/
	/*return;*/
	/*
	things we know on entry to this routine:
	total population of 2p: StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH2p].Pop
	total population of 1s: StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop
	continuum pumping rate (lo-up) inside 21 cm line: HFLines[0].pump
	upper to lower collision rate inside 21 cm line: HFLines[0].cs*dense.cdsqte
	occupation number inside Lya: OccupationNumberLine( &Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s] )

	level populations (cm-3) must be computed:
	population of upper level of 21cm: HFLines[0].Hi->Pop
	population of lower level of 21cm: HFLines[0].Lo->Pop
	stimulated emission corrected population of lower level: HFLines[0].Emis->PopOpc
	*/

	double x,
		PopTot;
	double a32,a31,a41,a42,a21, occnu_lya ,
		rate12 , rate21 , pump12 , pump21 , coll12 , coll21,
		texc , occnu_lya_23 , occnu_lya_13,occnu_lya_24,occnu_lya_14, texc1, texc2;

	PopTot = StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop;
	/* population can be zero in certain tests where H is turned off,
	 * also if initial solver does not see any obvious source of ionization
	 * also possible to set H0 density to zero with element ionization command,
	 * as is done in func_set_ion test case */
	if( PopTot <0 )
		TotalInsanity();
	else if( PopTot == 0 )
	{
		/*return after zeroing local variables */
		HFLines[0].Hi->Pop = 0.;
		HFLines[0].Lo->Pop = 0.;
		HFLines[0].Emis->PopOpc = 0.;
		HFLines[0].Emis->phots = 0.;
		HFLines[0].Emis->xIntensity = 0.;
		HFLines[0].Emis->ColOvTot = 0.;
		hyperfine.Tspin21cm = 0.;
		return;
	}

	a31 = 2.08e8;   /* Einstein co-efficient for transition 1p1/2 to 0s1/2 */
	a32 = 4.16e8;   /* Einstein co-efficient for transition 1p1/2 to 1s1/2 */
	a41 = 4.16e8;   /* Einstein co-efficient for transition 1p3/2 to 0s1/2 */
	a42 = 2.08e8;   /* Einstein co-efficient for transition 1p3/2 to 1s1/2 */
	/* These A values are determined from eqn. 17.64 of "The theory of Atomic structure
	 * and Spectra" by R. D. Cowan 
	 * A hyperfine level has degeneracy Gf=(2F + 1)
	 * a2p1s = 6.24e8;  Einstein co-efficient for transition 2p to 1s */
	a21 = 2.85e-15; /* Einstein co-efficient for transition 1s1/2 to 0s1/2 */

	/* above is spontaneous rate - the net rate is this times escape and destruction
	 * probabilities */
	a21 *= (HFLines[0].Emis->Pdest + HFLines[0].Emis->Pesc + HFLines[0].Emis->Pelec_esc);
	ASSERT( a21>0. );

	/* hyperfine.lgLya_pump_21cm is option to turn off Lya pump
	 * of 21 cm, with no 21cm lya pump command - note that this
	 * can be negative if Lya mases - can occur during search phase */
	occnu_lya = OccupationNumberLine( &Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s] ) *
		hyperfine.lgLya_pump_21cm;
	if( occnu_lya<0 )
	{
		static bool lgCommentDone = false;
		/* Lya is masing - could be due to very bad solution in search phase */
		if( !conv.lgSearch && !lgCommentDone )
		{
			fprintf(ioQQQ,
			"NOTE Lya masing will invert 21 cm, occupation number set zero\n");
			lgCommentDone = true;
		}
		occnu_lya = 0.;
	}

	/* Lya occupation number for the hyperfine levels 0S1/2 and 1S1/2 are different*/
	texc = TexcLine( &Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s] );
	/* Energy difference between 2p1/2 and 2p3/2 taken from NSRDS */
	if( texc > 0. )
	{
		/* convert to Boltzmann factor, which will applied to occupation 
		 * number of higher energy transition */
		texc1 = sexp(0.068/texc);
		texc2 = sexp(((82259.272-82258.907)*T1CM)/texc);
	}
	else
	{
		texc1 = 0.;
		texc2 = 0.;
	}

	/* the continuum within Lya seen by the two levels is not exactly the same brightness.  They
	 * differ by the exp when Lya is on Wein tail of black body, which must be true if 21 cm is important */

	occnu_lya_23 = occnu_lya;
	occnu_lya_13 = occnu_lya*texc1;
    occnu_lya_14 = occnu_lya_13*texc2;
    occnu_lya_24 = occnu_lya*texc2;

	/* this is the 21 cm upward continuum pumping rate [s-1] for the attenuated incident and
	 * local continuum and including line optical depths */
	pump12 = HFLines[0].Emis->pump;
	pump21 = pump12 * HFLines[0].Lo->g / HFLines[0].Hi->g;

	/* collision rates s-1 within 1s,
	 * were multiplied by collider density when evaluated in CoolEvaluate */
	/* ContBoltz is Boltzmann factor for wavelength of line */
	ASSERT( HFLines[0].Coll.col_str>0. );
	coll12 = HFLines[0].Coll.col_str*dense.cdsqte/HFLines[0].Lo->g*rfield.ContBoltz[HFLines[0].ipCont-1];
	coll21 = HFLines[0].Coll.col_str*dense.cdsqte/HFLines[0].Hi->g;

	/* set up rate (s-1) equations
	 * all process out of 1 that eventually go to 2 */
	rate12 = 
		/* collision rate (s-1) from 1 to 2 */
		coll12 + 
		/* direct external continuum pumping (s-1) in 21 cm line - usually dominated by CMB */
		pump12 +
		/* pump rate (s-1) up to 3, times fraction that decay to 2, hence net 1-2 */
		3.*a31*occnu_lya_13 *a32/(a31+a32)+
		/* pump rate (s-1) up to 4, times fraction that decay to 2, hence net 1-2 */
		/* >>chng 05 apr 04, GS, degeneracy corrected from 6 to 3 */
		3.*a41*occnu_lya_14 *a42/(a41+a42);

	/* set up rate (s-1) equations
	 * all process out of 2 that eventually go to 1 */
	/* spontaneous + induced 2 -> 1 by external continuum inside 21 cm line */
	/* >>chng 04 dec 03, do not include spontaneous decay, for numerical stability */
	rate21 = 
		/* collisional deexcitation */
		coll21 +
		/* net spontaneous decay plus external continuum pumping in 21 cm line */
		pump21 +
		/* rate from 2 to 3 time fraction that go back to 1, hence net 2 - 1 */
		/* >>chng 05 apr 04,GS, degeneracy corrected from 2 to unity */
		occnu_lya_23*a32 * a31/(a31+a32)+
		occnu_lya_24*a42*a41/(a41+a42);

	/* x = HFLines[0].Hi->Pop/HFLines[0].Lo->Pop */
	x = rate12 / SDIV(a21 + rate21);

	/* the Transitions term is the total population of 1s */
	HFLines[0].Hi->Pop = (x/(1.+x))* PopTot;
	HFLines[0].Lo->Pop = (1./(1.+x))* PopTot;
	ASSERT( HFLines[0].Hi->Pop >0. );
	ASSERT( HFLines[0].Lo->Pop >0. );

	/* the population with correction for stimulated emission */
	HFLines[0].Emis->PopOpc = HFLines[0].Lo->Pop*((3*rate21- rate12) + 3*a21)/SDIV(3*(a21+ rate21));

	/* number of escaping line photons, used elsewhere for outward beam */
	HFLines[0].Emis->phots = HFLines[0].Hi->Pop * HFLines[0].Emis->Aul *
		(HFLines[0].Emis->Pesc + HFLines[0].Emis->Pelec_esc);
	ASSERT( HFLines[0].Emis->phots >= 0. );
	/* intensity of line */
	HFLines[0].Emis->xIntensity = HFLines[0].Emis->phots*HFLines[0].EnergyErg;

	/* ratio of collisional to total (collisional + pumped) excitation */
	HFLines[0].Emis->ColOvTot = coll12 / rate12;

	/* finally save the spin temperature */
	if( HFLines[0].Hi->Pop > SMALLFLOAT )
	{
		hyperfine.Tspin21cm = TexcLine( &HFLines[0] );
		/* this line must be non-zero - it does strongly mase in limit_compton_hi_t sim -
		 * in that sim pop ratio goes to unity for a float and TexcLine ret zero */
		if( hyperfine.Tspin21cm == 0. )
			hyperfine.Tspin21cm = phycon.te;
	}
	else
	{
		hyperfine.Tspin21cm = phycon.te;
	}

	return;
}

/*H21cm_electron computes rate for H 21 cm from upper to lower excitation by electrons - call by CoolEvaluate
 * >>refer	H1	cs	Smith, F.J., 1966, Planet. Space Sci 14, 929 */
double H21cm_electron( double temp )
{
	double hold;
	temp = MIN2(1e4 , temp );
	/* following fit is from */
	/* >>refer	H1	21cm	Liszt, H., 2001, A&A, 371, 698 */

	hold = -9.607 + log10( sqrt(temp)) * sexp( pow(log10(temp) , 4.5 ) / 1800. );
	hold = pow(10.,hold );
	return( hold );
}

/* computes rate for H 21 cm from upper to lower excitation by atomic hydrogen 
 * from 
 * >>refer	H1	cs	Allison, A.C., & Dalgarno A., 1969, ApJ 158, 423 */
/* the following is the best current survey of 21 cm excitation */
/* >>refer	H1	21cm	Liszt, H., 2001, A&A, 371, 698 */
#if 0
STATIC double h21_t_ge_20( double temp )
{
	double y;
	double x1,
		teorginal = temp;
	/* data go up to 1,000K must not go above this */
	temp = MIN2( 1000.,temp );
	x1 =1.0/sqrt(temp);
	y =-21.70880995483007-13.76259674006133*x1;
	y = exp(y);

	/* >>chng 02 feb 14, extrapolate above 1e3 K as per Liszt 2001 recommendation 
	 * page 699 of */
	/* >>refer	H1	21cm	Liszt, H., 2001, A&A, 371, 698 */
	if( teorginal > 1e3 )
	{
		y *= pow(teorginal/1e3 , 0.33 );
	}

	return( y );
}

/* this branch for T < 20K, data go down to 1 K */
STATIC double h21_t_lt_20( double temp )
{
	double y;
	double x1;

	/* must not go below 1K */
	temp = MAX2( 1., temp );
	x1 =temp*log(temp);
	y =9.720710314268267E-08+6.325515312006680E-08*x1;
	return(y*y);
}
#endif

/* >> chng 04 dec 15, GS. The fitted rate co-efficients (cm3s-1) in the temperature range 1K to 300K is from
 * >>refer	H1	cs	Zygelman, B. 2005, ApJ preprint doi:10.1086/'427682". 
 * The rate is 4/3 times the Dalgarno (1969) rate for the 
 temperature range 300K to 1000K. Above 1000K, the rate is extrapolated according to Liszt 2001.*/
STATIC double h21_t_ge_10( double temp )
{
	double y;
	double x1,x2,x3,
	teorginal = temp;
	/* data go up to 300K  */
	temp = MIN2( 300., temp );
	x1 =temp;
	y =1.4341127e-9+9.4161077e-15*x1-9.2998995e-9/(log(x1))+6.9539411e-9/sqrt(x1)+1.7742293e-8*(log(x1))/pow2(x1);
	if( teorginal > 300. )
	{
		/* data go up to 1000*/
		x3 = MIN2( 1000., teorginal );
		x2 =1.0/sqrt(x3);
		y =-21.70880995483007-13.76259674006133*x2;
		y = 1.236686*exp(y);

	}
	if( teorginal > 1e3 )
	{
		/*data go above 1000*/
		y *= pow(teorginal/1e3 , 0.33 );
	}
	return( y );
}
/* this branch for T < 10K, data go down to 1 K */
STATIC double h21_t_lt_10( double temp )
{
	double y;
	double x1;

	/* must not go below 1K */
	temp = MAX2(1., temp );
	x1 =temp;
	y =8.5622857e-10+2.331358e-11*x1+9.5640586e-11*pow2(log(x1))-4.6220869e-10*sqrt(x1)-4.1719545e-10/sqrt(x1);
	return(y);
}

/*H21cm_H_atom - evaluate H atom spin changing H atom collision rate, 
 * called by CoolEvaluate 
 * >>refer	H1	cs	Allison, A.C. & Dalgarno, A., 1969, ApJ 158, 423 
 */
double H21cm_H_atom( double temp )
{
	double hold;
	if( temp >= 10. )
	{
		hold = h21_t_ge_10( temp );
	}
	else
	{
		hold = h21_t_lt_10( temp );
	}

	return hold;
}

/*H21cm_proton - evaluate proton spin changing H atom collision rate, 
* called by CoolEvaluate */
double H21cm_proton( double temp )
{
	/*>>refer	21cm	p coll	Furlanetto, S.R. & Furlanetto, M.R. 2007,
	 *>>refcon	MNRAS, doi:10.1111/j.1365-2966.2007.11921.x 
	 * previously had used proton rate, which is 3.2 times H0 rate according to
	 *>>refer	21cm	cs	Liszt, H. A&A, 371, 698 */
	/* fit to table 1 of first paper */
	/*--------------------------------------------------------------*
	TableCurve Function: c:\storage\litera~1\21cm\atomic~1\p21cm.c Jun 20, 2007 3:37:50 PM
	proton coll deex
	X= temperature (K)
	Y= rate coefficient (1e-9 cm3 s-1)
	Eqn# 4419  y=a+bx+cx^2+dx^(0.5)+elnx/x
	r2=0.9999445384690351
	r2adj=0.9999168077035526
	StdErr=5.559328579039901E-12
	Fstat=49581.16793656295
	a= 9.588389834316704E-11
	b= -5.158891920816405E-14
	c= 5.895348443553458E-19
	d= 2.05304960232429E-11
	e= 9.122617940315725E-10
	*--------------------------------------------------------------*/

	double hold;
	/* only fit this range, did not include T = 1K point which 
	 * causes an inflection */
	temp = MAX2( 2. , temp );
	temp = MIN2( 2e4 , temp );

	/* within range of fitted rate coefficients */
	double x1,x2,x3,x4;
	x1 = temp;
	x2 = temp*temp;
	x3 = sqrt(temp);
	x4 = log(temp)/temp;
	hold =9.588389834316704E-11 - 5.158891920816405E-14*x1
		+5.895348443553458E-19*x2 + 2.053049602324290E-11*x3
		+9.122617940315725E-10*x4;

	return hold;
}

/* 
 * HyperfineCreate, HyperfineCS written July 2001
 * William Goddard for Gary Ferland
 * This code calculates line intensities for known
 * hyperfine transitions.
 */

/* two products, the transition structure HFLines, which contains all information for the lines,
 * and nHFLines, the number of these lines.  
 *
 * these are in taulines.h
 *
 * info to create them contained in hyperfine.dat
 *
 * abundances of nuclei are also in hyperfine.dat, stored in 
 */

/* Ion contains twelve varying temperatures, specified above, used for */
/* calculating collision strengths.									   */	
typedef struct 
{
	double strengths[12];
} Ion;

static	Ion *Strengths;

/* HyperfineCreate establish space for hf arrays, reads atomic data from hyperfine.dat */
void HyperfineCreate(void)
{
	FILE *ioDATA;
	char chLine[INPUT_LINE_LENGTH];
	bool lgEOL;
	/*double c, h, k, N, Ne, q12, q21, upsilon, x;*/
	realnum spin, wavelength;
	long int i, j, mass, nelec, ion, nelem;

	DEBUG_ENTRY( "HyperfineCreate()" );

	/* list of ion collision strengths for the temperatures listed in table */
	/* HFLines containing all the data in Hyperfine.dat, and transition is		*/
	/* defined in cddefines.h												*/

	/*transition *HFLines;*/

	/* get the line data for the hyperfine lines */
	if( trace.lgTrace )
		fprintf( ioQQQ," Hyperfine opening hyperfine.dat:");

	ioDATA = open_data( "hyperfine.dat", "r" );

	/* first line is a version number and does not count */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " Hyperfine could not read first line of hyperfine.dat.\n");
		cdEXIT(EXIT_FAILURE);
	}
	/* count how many lines are in the file, ignoring all lines
	 * starting with '#' */
	nHFLines = 0;
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
	{
		/* we want to count the lines that do not start with #
		 * since these contain data */
		if( chLine[0] != '#')
			++nHFLines;
	}

	/* MALLOC the transition HFLines array */
	HFLines = (transition *)MALLOC( (size_t)(nHFLines)*sizeof(transition) );
	/* initialize array to impossible values to make sure eventually done right */
	for( i=0; i< nHFLines; ++i )
	{
		HFLines[i].Junk();
		HFLines[i].Hi = AddState2Stack();
		HFLines[i].Lo = AddState2Stack();
		HFLines[i].Emis = AddLine2Stack( true );
	}

	Strengths = (Ion *)MALLOC( (size_t)(nHFLines)*sizeof(Ion) );
	hyperfine.HFLabundance = (realnum *)MALLOC( (size_t)(nHFLines)*sizeof(realnum) );

	/* now rewind the file so we can read it a second time*/
	if( fseek( ioDATA , 0 , SEEK_SET ) != 0 )
	{
		fprintf( ioQQQ, " Hyperfine could not rewind hyperfine.dat.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* check that magic number is ok, read the line */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " Hyperfine could not read first line of hyperfine.dat.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* check that magic number is ok, scan it in */
	i = 1;
	nelem = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
	nelec = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
	ion = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);

	{
		/* the following is the set of numbers that appear at the start of hyperfine.dat 01 07 12 */ 
		static int iYR=2,  iMN=10, iDY=28;
		if( ( nelem != iYR ) || ( nelec != iMN ) || ( ion != iDY ) )
		{
			fprintf( ioQQQ, 
				" Hyperfine: the version of hyperfine.dat in the data directory is not the current version.\n" );
			fprintf( ioQQQ, 
				" I expected to find the number %i %i %i and got %li %li %li instead.\n" ,
				iYR, iMN , iDY ,
				nelem , nelec , ion );
			fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* 
	 * scan the string taken from Hyperfine.dat, parsing into
	 * needed variables.
	 * nelem is the atomic number.
	 * IonStg is the ionization stage.  Atom = 1, Z+ = 2, Z++ = 3, etc.
	 * Aul is used to find the einstein A in the function GetGF.
	 * most of the variables are floats.
	 */

	/* this will count the number of lines */
	j = 0;

	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
	{
		/* only look at lines without '#' in first col */
		/* make sure line in data file is < 140 char   */
		if( chLine[0] != '#')
		{
			double Aul;
			ASSERT( j < nHFLines );

			double help[3];
			sscanf(chLine, "%li%i%le%i%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le", 
			       &mass, 
			       &HFLines[j].Hi->nelem,
			       &help[0],
			       &HFLines[j].Hi->IonStg,
			       &help[1],
			       &help[2],
			       &Aul,
			       &Strengths[j].strengths[0],
			       &Strengths[j].strengths[1],
			       &Strengths[j].strengths[2],
			       &Strengths[j].strengths[3],
			       &Strengths[j].strengths[4],
			       &Strengths[j].strengths[5],
			       &Strengths[j].strengths[6],
			       &Strengths[j].strengths[7],
			       &Strengths[j].strengths[8],
			       &Strengths[j].strengths[9],
			       &Strengths[j].strengths[10],
			       &Strengths[j].strengths[11]);
			spin = (realnum)help[0];
			hyperfine.HFLabundance[j] = (realnum)help[1];
			wavelength = (realnum)help[2];
			HFLines[j].Emis->Aul = (realnum)Aul;
			HFLines[j].Emis->damp = 1e-20f;
			HFLines[j].Hi->g = (realnum) (2*(spin + .5) + 1);
			HFLines[j].Lo->g = (realnum) (2*(spin - .5) + 1);
			HFLines[j].WLAng = wavelength * 1e8f;
			HFLines[j].EnergyWN = (realnum) (1. / wavelength);
			HFLines[j].Emis->gf = (realnum)(GetGF(HFLines[j].Emis->Aul, HFLines[j].EnergyWN, HFLines[j].Hi->g));
			/*fprintf(ioQQQ,"HFLinesss1\t%li\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
				j,HFLines[j].opacity , HFLines[j].Emis->gf , HFLines[j].Emis->Aul , HFLines[j].EnergyWN,HFLines[j].Lo->g);*/
			HFLines[j].Lo->nelem = HFLines[j].Hi->nelem;
			HFLines[j].Lo->IonStg = HFLines[j].Hi->IonStg;


			ASSERT(HFLines[j].Emis->gf > 0.);
			/* increment the counter */
			j++;
		}

	}
	ASSERT( j==nHFLines );

	/* closing the file */
	fclose(ioDATA);

#	if 0
	/* for debugging and developing only */
	/* calculating the luminosity for each isotope */
	for(i = 0; i < nHFLines; i++)
	{
		N = dense.xIonDense[HFLines[i].Hi->nelem-1][HFLines[i].Hi->IonStg-1];
		Ne = dense.eden;

		h = 6.626076e-27;			/* erg * sec */
		c = 3e10;					/* cm / sec	 */
		k = 1.380658e-16;			/* erg / K   */

		upsilon = HyperfineCS(i);
				/*statistical weights must still be identified */
		q21 = COLL_CONST * upsilon / (phycon.sqrte * HFLines[i].Hi->g);

		q12 = HFLines[i].Hi->g/ HFLines[i].Lo->g * q21 * exp(-1 * h * c * HFLines[i].EnergyWN / (k * phycon.te)); 

		x = Ne * q12 / (HFLines[i].Emis->Aul * (1 + Ne * q21 / HFLines[i].Aul));
		HFLines[i].xIntensity = N * HFLines[i].Emis->Aul * x / (1.0 + x) * h * c / (HFLines[i].WLAng / 1e8);

	}
#	endif

	return;
}


/*HyperfineCS returns interpolated collision strength for element nelem and ion ion */
double HyperfineCS( long i )
{

	/*Search HFLines to find i of ion 		*/
	int /*i = 0,*/ j = 0;
#	define N_TE_TABLE 12
	double 	slope, upsilon, TemperatureTable[N_TE_TABLE] = {.1e6, .15e6, .25e6, .4e6, .6e6, 
		1.0e6, 1.5e6, 2.5e6, 4e6, 6e6, 10e6, 15e6};

	DEBUG_ENTRY( "HyperfineCS()" );

	/*while(i < nHFLines+1 && HFLines[i].Hi->nelem != nelem && HFLines[i].Hi->IonStg != ion)
		++i;*/
	ASSERT( i >= 0. && i <= nHFLines );

	/*j = temperature	*/
	/* calculate actual cooling rate for isotope.								*/
	/* The temperature-dependent collision strength must first be calculated.	*/
	/* phycon.te is compared to the first, last and intermediate values of strengths[]*/
	/* and interpolated by taking the log of the collision strength				*/
	/* and temperature.															*/

	if( phycon.te <= TemperatureTable[0])
	{
		/* temperature below bounds of table */
		upsilon = Strengths[i].strengths[0];
	}
	else if(  phycon.te >= TemperatureTable[N_TE_TABLE-1])
	{
		/* temperature above bounds of table */
		upsilon = Strengths[i].strengths[N_TE_TABLE-1];
	}
	else
	{
		j = 1;
		/* want Table[j-1] < te < Table[j] */
		/*while ( phycon.te >= TemperatureTable[j] && j < N_TE_TABLE )*/
		/*while ( TemperatureTable[j] < phycon.te  && j < N_TE_TABLE )*/
		while ( j < N_TE_TABLE && TemperatureTable[j] < phycon.te )
			j++;

		ASSERT( j >= 0 && j < N_TE_TABLE);
		ASSERT( (TemperatureTable[j-1] <= phycon.te ) && (TemperatureTable[j] >= phycon.te) );

		if( fp_equal( phycon.te, TemperatureTable[j] ) )
		{
			upsilon = Strengths[i].strengths[j];
		}
		/*phycon.te must be less than TemperatureTable[j], greater than TemperatureTable[j-1][0]*/
		else if( phycon.te < TemperatureTable[j])
		{
			slope = (log10(Strengths[i].strengths[j-1]) - log10(Strengths[i].strengths[j])) /
				  (log10(TemperatureTable[j-1]) - log10(TemperatureTable[j]));

			upsilon = (log10((phycon.te)) - (log10(TemperatureTable[j-1])-6.))*slope + log10(Strengths[i].strengths[j-1]);
			upsilon = pow(10., upsilon);
		}
		else 
		{
			upsilon = Strengths[i].strengths[j-1];
		}
	}

	return upsilon;
}
