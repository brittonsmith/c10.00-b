/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*tfidle update some temperature dependent variables */
/*tauff compute optical depth where cloud is thin to free-free and plasma freq */
#include "cddefines.h"
#include "physconst.h"
#include "opacity.h"
#include "iso.h"
#include "dense.h"
#include "phycon.h"
#include "stopcalc.h"
#include "continuum.h"
#include "trace.h"
#include "rfield.h"
#include "doppvel.h"
#include "radius.h"
#include "wind.h"
#include "thermal.h"
#include "conv.h"

/*tauff compute optical depth where cloud is thin to free-free and plasma freq */
STATIC void tauff(void);
/* On first run, fill GauntFF with gaunt factors	*/
STATIC void FillGFF(void);
/* Interpolate on GauntFF to calc gaunt at current temp, phycon.te	*/
STATIC realnum InterpolateGff( long charge , double ERyd );
STATIC int LinterpTable(realnum **t, realnum *v, long int lta, long int ltb, realnum x, realnum *a, long int *pipx);
STATIC int LinterpVector(realnum **t, realnum *v, long lta , long ltb, realnum *yy , long ny, realnum **a);
STATIC void fhunt(realnum *xx, long int n, realnum x, long int *j);

/**tfidle update some temperature dependent variables 
\param lgForceUpdate option to force update of all variables 
*/
STATIC void tfidle(
			bool lgForceUpdate);

static long lgGffNotFilled = true;

const long N_TE_GFF = 41;
static long N_PHOTON_GFF; /* will cover full energy range full range in one-tenth dec steps	*/
static realnum ***GauntFF;
static realnum **GauntFF_T;
/* the array of logs of temperatures at which GauntFF is defined */
static realnum TeGFF[N_TE_GFF];
/* the array of logs of u at which GauntFF is defined	*/
static realnum *PhoGFF;

/**TempChange change kinetic temperature, calls tfidle
*/
void TempChange(
			 double TempNew ,
			 /* option to force update of all variables */
			 bool lgForceUpdate)
{

	DEBUG_ENTRY( "TempChange()" );

	/* set new temperature */
	if( TempNew > phycon.TEMP_LIMIT_HIGH )
	{
		/* temp is too high */
		fprintf(ioQQQ," PROBLEM DISASTER - the kinetic temperature, %.3eK,"
			" is above the upper limit of the code, %.3eK.\n",
			TempNew , phycon.TEMP_LIMIT_HIGH );
		fprintf(ioQQQ," This calculation is aborting.\n Sorry.\n");

		TempNew = phycon.TEMP_LIMIT_HIGH*0.99999;
		lgAbort = true;
	}
	else if( TempNew < phycon.TEMP_LIMIT_LOW )
	{
		/* temp is too low */
		fprintf(ioQQQ," PROBLEM DISASTER - the kinetic temperature, %.3eK,"
			" is below the lower limit of the code, %.3eK.\n",
			TempNew , phycon.TEMP_LIMIT_LOW );
		fprintf(ioQQQ," This calculation is aborting.\n Sorry.\n");
		TempNew = phycon.TEMP_LIMIT_LOW*1.00001;
		lgAbort = true;
	}
	else if( TempNew < StopCalc.TeFloor )
	{
		/* temperature floor option  -
		 * go to constant temperature calculation if temperature
		 * falls below floor */
		thermal.lgTemperatureConstant = true;
		thermal.ConstTemp = (realnum)StopCalc.TeFloor;
		phycon.te = thermal.ConstTemp;
		/*fprintf(ioQQQ,"DEBUG TempChange hit temp floor, setting const temp to %.3e\n",
			phycon.te );*/
	}
	else
	{
		/* temp is within range */
		phycon.te = TempNew;
	}

	/* now update related variables */
	tfidle( lgForceUpdate );
	return;
}
/**TempChange change kinetic temperature, calls tfidle
 * but does not update extensive variables or check for temperature floor,
 * intended for use by routines that are sanity checks rather than real calculation */
void TempChange(
				double TempNew )
{

	DEBUG_ENTRY( "TempChange()" );

	/* set new temperature */
	if( TempNew > phycon.TEMP_LIMIT_HIGH )
	{
		/* temp is too high */
		fprintf(ioQQQ," PROBLEM DISASTER - the kinetic temperature, %.3eK,"
			" is above the upper limit of the code, %.3eK.\n",
			TempNew , phycon.TEMP_LIMIT_HIGH );
		fprintf(ioQQQ," This calculation is aborting.\n Sorry.\n");

		TempNew = phycon.TEMP_LIMIT_HIGH*0.99999;
		lgAbort = true;
	}
	else if( TempNew < phycon.TEMP_LIMIT_LOW )
	{
		/* temp is too low */
		fprintf(ioQQQ," PROBLEM DISASTER - the kinetic temperature, %.3eK,"
			" is below the lower limit of the code, %.3eK.\n",
			TempNew , phycon.TEMP_LIMIT_LOW );
		fprintf(ioQQQ," This calculation is aborting.\n Sorry.\n");
		TempNew = phycon.TEMP_LIMIT_LOW*1.00001;
		lgAbort = true;
	}
	else
	{
		/* temp is within range */
		phycon.te = TempNew;
	}

	/* now update related variables */
	tfidle( false );
	return;
}

void tfidle(
	/* option to force update of all variables */
	bool lgForceUpdate)
{
	static double tgffused=-1., 
	  tgffsued2=-1.;
	static long int nff_defined=-1;
	static long maxion = 0, oldmaxion = 0;
	static double ttused = 0.;
	static bool lgZLogSet = false;
	bool lgGauntF;
	long int ion;
	long int i,
	  nelem,
	  if1,
		ipTe,
		ret;
	double fac,  
	  fanu;

	DEBUG_ENTRY( "tfidle()" );

	/* called with lgForceUpdate true in zero.c, when we must update everything */
	if( lgForceUpdate )
	{
		ttused = -1.;
		tgffused = -1.;
		tgffsued2 = -1.;
	}

	/* check that eden not negative */
	if( dense.eden <= 0. )
	{
		fprintf( ioQQQ, "tfidle called with a zero or negative electron density,%10.2e\n", 
		  dense.eden );
		TotalInsanity();
	}

	/* check that temperature not negative */
	if( phycon.te <= 0. )
	{
		fprintf( ioQQQ, "tfidle called with a negative electron temperature,%10.2e\n", 
		  phycon.te );
		TotalInsanity();
	}

	/* one time only, set up array of logs of charge squared */
	if( !lgZLogSet )
	{
		for( nelem=0; nelem<LIMELM; ++nelem )
		{
			/* this array is used to modify the log temperature array
			 * defined below, for hydrogenic species of charge nelem+1 */
			phycon.sqlogz[nelem] = log10( POW2(nelem+1.) );
		}
		lgZLogSet = true;
	}

	if( ! fp_equal( phycon.te, ttused ) )
	{
		ttused = phycon.te;
		thermal.te_update = phycon.te;
		/* current temperature in various units */
		phycon.te_eV = phycon.te/EVDEGK;
		phycon.te_ryd = phycon.te/TE1RYD;
		phycon.te_wn = phycon.te / T1CM;

		phycon.tesqrd = POW2(phycon.te);
		phycon.sqrte = sqrt(phycon.te);
		thermal.halfte = 0.5/phycon.te;
		thermal.tsq1 = 1./phycon.tesqrd;
		phycon.te32 = phycon.te*phycon.sqrte;
		phycon.teinv = 1./phycon.te;

		phycon.alogte = log10(phycon.te);
		phycon.alnte = log(phycon.te);

		phycon.telogn[0] = phycon.alogte;
		for( i=1; i < 7; i++ )
		{
			phycon.telogn[i] = phycon.telogn[i-1]*phycon.telogn[0];
		}

		phycon.te10 = pow(phycon.te,0.10);
		phycon.te20 = phycon.te10 * phycon.te10;
		phycon.te30 = phycon.te20 * phycon.te10;
		phycon.te40 = phycon.te30 * phycon.te10;
		phycon.te70 = phycon.sqrte * phycon.te20;
		phycon.te90 = phycon.te70 * phycon.te20;

		phycon.te01 = pow(phycon.te,0.01);
		phycon.te02 = phycon.te01 * phycon.te01;
		phycon.te03 = phycon.te02 * phycon.te01;
		phycon.te04 = phycon.te02 * phycon.te02;
		phycon.te05 = phycon.te03 * phycon.te02;
		phycon.te07 = phycon.te05 * phycon.te02;

		phycon.te001 = pow(phycon.te,0.001);
		phycon.te002 = phycon.te001 * phycon.te001;
		phycon.te003 = phycon.te002 * phycon.te001;
		phycon.te004 = phycon.te002 * phycon.te002;
		phycon.te005 = phycon.te003 * phycon.te002;
		phycon.te007 = phycon.te005 * phycon.te002;
		/*>>>chng 06 June 30 -Humeshkar Nemala*/
		phycon.te0001 = pow(phycon.te ,0.0001);
		phycon.te0002 = phycon.te0001 * phycon.te0001;
		phycon.te0003 = phycon.te0002 * phycon.te0001;
		phycon.te0004 = phycon.te0002 * phycon.te0002;
		phycon.te0005 = phycon.te0003 * phycon.te0002;
		phycon.te0007 = phycon.te0005 * phycon.te0002;

	}

	/* >>>chng 99 nov 23, removed this line, so back to old method of h coll */
	/* used for hydrogenic collisions */
	/* 
	 * following electron density has approximate correction for neutrals
	 * corr of hi*1.7e-4 accounts for col ion by HI; 
	 * >>refer	H0	correction for collisional contribution		Drawin, H.W. 1969, Zs Phys 225, 483.
	 * also quoted in Dalgarno & McCray 1972
	 * extensive discussion of this in 
	 *>>refer	H0	collisions	Lambert, D.L. 
	 * used EdenHCorr instead
	 * edhi = eden + hi * 1.7e-4
	 */
	dense.EdenHCorr = dense.eden + 
		/* dense.HCorrFac is unity by default and changed with the set HCOR command */
		dense.xIonDense[ipHYDROGEN][0]*1.7e-4 * dense.HCorrFac;

	/*>>chng 93 jun 04,
	 * term with hi added June 4, 93, to account for warm pdr */
	/* >>chng 05 jan 05, Will Henney noticed that 1.e-4 used here is not same as
	 * 1.7e-4 used for EdenHCorr, which had rewritten the expression.
	 * change so that edensqte uses EdenHCorr rather than reevaluating */
	/*dense.edensqte = ((dense.eden + dense.xIonDense[ipHYDROGEN][0]/1e4)/phycon.sqrte);*/
	dense.edensqte = dense.EdenHCorr/phycon.sqrte;
	dense.cdsqte = dense.edensqte*COLL_CONST;
	dense.SqrtEden = sqrt(dense.eden);

	/* rest have to do with radiation field which may not be defined yet */
	if( !lgRfieldMalloced )
	{
		return;
	}

	/* correction factors for induced recombination, 
	 * also used as Boltzmann factors
	 * check for zero is because ContBoltz is zeroed out in initialization
	 * of code, its possible this is a constant density grid of models
	 * in which the code is called as a subroutine */
	/* >>chng 01 aug 21, must also test on size of continuum nflux because 
	 * conintitemp can increase nflux then call this routine, although 
	 * temp may not have changed */
	if( ! fp_equal(tgffused, phycon.te)
		|| rfield.ContBoltz[0] <= 0. || nff_defined<rfield.nflux )
	{
		tgffused = phycon.te;
		fac = TE1RYD/phycon.te;
		i = 0;
		fanu = fac*rfield.anu[i];
		/* SEXP_LIMIT is 84 in cddefines.h - it is the -ln of smallest number
		 * that sexp can handle, and is used elsewhere in the code.  
		 * atom_level2 uses ContBoltz to see whether the rates will be significant.
		 * If the numbers did not agree then this test would be flawed, resulting in
		 * div by zero */
		while( i < rfield.nupper && fanu < SEXP_LIMIT )
		{
			rfield.ContBoltz[i] = exp(-fanu);
			++i;
			/* this is Boltzmann factor for NEXT cell */
			fanu = fac*rfield.anu[i];
		}
		/* ipMaxBolt is number of cells defined, so defined up through ipMaxBolt-1 */
		rfield.ipMaxBolt = i;

		/* zero out remainder */
		/* >>chng 01 apr 14, upper limit has been ipMaxBolt+1, off by one */
		for( i=rfield.ipMaxBolt; i < rfield.nupper; i++ )
		{
			rfield.ContBoltz[i] = 0.;
		}
	}

	/* find frequency where thin to bremsstrahlung or plasma frequency */
	tauff();

	oldmaxion = maxion;
	maxion = 0;

	/* Find highest maximum stage of ionization	*/
	for( nelem = 0; nelem < LIMELM; nelem++ )
	{
		maxion = MAX2( maxion, dense.IonHigh[nelem] );
	}

	/* reevaluate if temperature or number of cells has changed */
	if( !fp_equal(phycon.te,tgffsued2) || 
	    /* this test - reevaluate if upper bound of defined values is
	     * above nflux, the highest point.  This only triggers in
	     * large grids when continuum source gets harder */
	    nff_defined < rfield.nflux  ||
	    /* this occurs when now have more stages of ionization than in previous defn */
	    maxion > oldmaxion )
	{
		static bool lgFirstRunDone = false;
		long lowion;
		/* >>chng 02 jan 10, only reevaluate needed states of ionization */
		if( fp_equal(phycon.te,tgffsued2) && nff_defined >= rfield.nflux && 
		    maxion > oldmaxion )
		{
			/* this case temperature did not change by much, but upper
			 * stage of ionization increased.  only need to evaluate
			 * stages that have not been already done */
			lowion = oldmaxion;
		}
		else
		{
			/* temperature changed - do them all */
			lowion = 1;
		}

		/* if1 will certainly be set to some positive value in gffsub */
		if1 = 1;

		/* chng 02 may 16, by Ryan...one gaunt factor array for all charges	*/
		/* First index is EFFECTIVE CHARGE!	*/
		/* highest stage of ionization is LIMELM, 
		 * index goes from 1 to LIMELM */

		nff_defined = rfield.nflux;

		ASSERT( if1 >= 0 );

		tgffsued2 = phycon.te;
		lgGauntF = true;

		/* new gaunt factors	*/
		if( lgGffNotFilled )
		{
			FillGFF();
		}

		if( lgFirstRunDone == false )
		{
			maxion = LIMELM;
			lgFirstRunDone = true;
		}

		/* >> chng 03 jan 23, rjrw -- move a couple of loops down into
		 * subroutines, and make those subroutines generic
		 * (i.e. dependences only on arguments, may be useful elsewhere...) */

		/* Make Gaunt table for new temperature */
		ipTe = -1;
		for( ion=1; ion<=LIMELM; ion++ )
		{
			/* Interpolate for all tabulated photon energies at this temperature */
			ret = LinterpTable(GauntFF[ion], TeGFF, N_PHOTON_GFF, N_TE_GFF, (realnum)phycon.alogte, GauntFF_T[ion], &ipTe);
			if(ret == -1) 
			{
				fprintf(ioQQQ," LinterpTable for GffTable called with te out of bounds \n");
				cdEXIT(EXIT_FAILURE);
			}
		}

		/* Interpolate for all ions at required photon energies -- similar
		 * to LinterpTable, but not quite similar enough... */
		/* >>chng 04 jun 30, change nflux to nupper */
		ret = LinterpVector(GauntFF_T+lowion, PhoGFF, maxion-lowion+1, N_PHOTON_GFF,
			rfield.anulog, rfield.nupper,/*rfield.nflux,*/ rfield.gff + lowion); 
		if(ret == -1) 
		{
			fprintf(ioQQQ," LinterpVector for GffTable called with photon energy out of bounds \n");
			cdEXIT(EXIT_FAILURE);
		}
	}
	else
	{
		/* this is flag that would have been set in gffsub, and
		 * printed in debug statement below.  We are not evaluating
		 * so set to -1 */
		if1 = -1;
		lgGauntF = false;
	}

	if( trace.lgTrace && trace.lgTrGant )
	{
		fprintf( ioQQQ, "     tfidle; gaunt factors?" );
		fprintf( ioQQQ, "%2c", TorF(lgGauntF) );

		fprintf( ioQQQ, "%2f g 1 2=%10.2e%9.2ld flag%2f guv(1,n)%10.2e\n", 
		  rfield.gff[1][0], rfield.gff[1][iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][2]-1],
		  if1, rfield.gff[1][iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][2]], 
		  rfield.gff[1][rfield.nflux-1] );
	}
	return;
}

/*tauff compute optical depth where cloud is thin to free-free and plasma freq */
STATIC void tauff(void)
{
	realnum fac;

	/* simply return if space not yet allocated */
	if( !lgOpacMalloced )
		return;

	DEBUG_ENTRY( "tauff()" );

	/* routine sets variable ipEnergyBremsThin, index for lowest cont cell that is optically thin */
	/* find frequency where continuum thin to free-free */
	while( rfield.ipEnergyBremsThin < rfield.nflux && 
		opac.TauAbsGeo[0][rfield.ipEnergyBremsThin] >= 1. )
	{
		++rfield.ipEnergyBremsThin;
	}

	/* TFF will be frequency where cloud becomes optically thin to bremsstrahlung
	 * >>chng 96 may 7, had been 2, change as per Kevin Volk bug report */
	if( rfield.ipEnergyBremsThin > 1 && opac.TauAbsGeo[0][rfield.ipEnergyBremsThin] > 0.001 )
	{
		/* tau can be zero when plasma frequency is within energy grid, */
		fac = (1.f - opac.TauAbsGeo[0][rfield.ipEnergyBremsThin-1])/(opac.TauAbsGeo[0][rfield.ipEnergyBremsThin] - 
		  opac.TauAbsGeo[0][rfield.ipEnergyBremsThin-1]);
		fac = MAX2(fac,0.f);
		rfield.EnergyBremsThin = rfield.anu[rfield.ipEnergyBremsThin-1] + rfield.widflx[rfield.ipEnergyBremsThin-1]*fac;
	}
	else
	{
		rfield.EnergyBremsThin = 0.f;
	}

	/* evaluate the plasma frequency one time per zone to avoid PF moving across
	 * a line during convergence loops */
	if( nzone != rfield.nZonePlsFrqEval )
	{
		rfield.nZonePlsFrqEval = nzone;
		rfield.plsfrq = (realnum)((ELEM_CHARGE_ESU/sqrt(PI*ELECTRON_MASS)/FR1RYD)*sqrt(dense.eden));

		if( rfield.ipPlasma > 0 )
		{
			/* increase index for plasma frequency until within proper cell */
			while( rfield.plsfrq > rfield.anu[rfield.ipPlasma]+rfield.widflx[rfield.ipPlasma]/2. )
				++rfield.ipPlasma;

			/* decrease index for plasma frequency until within proper cell */
			while( rfield.ipPlasma>2 && rfield.plsfrq < rfield.anu[rfield.ipPlasma]-rfield.widflx[rfield.ipPlasma]/2. )
				--rfield.ipPlasma;
		}

		/* also remember the largest plasma frequency we encounter */
		rfield.plsfrqmax = MAX2(rfield.plsfrqmax, rfield.plsfrq);

		/* is plasma frequency within energy grid? */
		if( rfield.plsfrq > rfield.anu[0] )
		{
			rfield.lgPlasNu = true;
		}
	}

	/* did not include plasma freq before
	 * function returns larger of these two frequencies */
	rfield.EnergyBremsThin = MAX2(rfield.plsfrq,rfield.EnergyBremsThin);

	/* now increment ipEnergyBremsThin still further, until above plasma frequency */
	while( rfield.ipEnergyBremsThin < rfield.nflux && 
		rfield.anu[rfield.ipEnergyBremsThin] <= rfield.EnergyBremsThin )
	{
		++rfield.ipEnergyBremsThin;
	}
	return;
}

realnum GetDopplerWidth( realnum massAMU )
{
	ASSERT( massAMU > 0. );
	// force a fairly conservative upper limit
	ASSERT( massAMU < 10000. );
	
	/* usually TurbVel =0, reset with turbulence command
	 * cm/s here, but was entered in km/s with command */
	double turb2 = POW2(DoppVel.TurbVel);

	/* this is option to dissipate the turbulence.  DispScale is entered with
	 * dissipate keyword on turbulence command.  The velocity is reduced here,
	 * by an assumed exponential scale, and also adds heat */
	if( DoppVel.DispScale > 0. )
	{
		/* square of exp depth dependence */
		turb2 *= sexp( 2.*radius.depth / DoppVel.DispScale );
	}

	/* in case of D-Critical flow include initial velocity as
	 * a component of turbulence */
	if( ! ( wind.lgBallistic() || wind.lgStatic() ) )
	{
		turb2 += POW2(wind.windv0);
	}

	realnum width = (realnum)sqrt(2.*BOLTZMANN/ATOMIC_MASS_UNIT*phycon.te/massAMU+turb2);
	ASSERT( width > 0.f );
	return width;
}

realnum GetAveVelocity( realnum massAMU )
{
#if 0
	/* usually TurbVel =0, reset with turbulence command
	 * cm/s here, but was entered in km/s with command */
	double turb2 = POW2(DoppVel.TurbVel);

	/* this is option to dissipate the turbulence.  DispScale is entered with
	 * dissipate keyword on turbulence command.  The velocity is reduced here,
	 * by an assumed exponential scale, and also adds heat */
	if( DoppVel.DispScale > 0. )
	{
		/* square of exp depth dependence */
		turb2 *= sexp( 2.*radius.depth / DoppVel.DispScale );
	}

	/* in case of D-Critical flow include initial velocity as
	 * a component of turbulence */
	if( ! ( wind.lgBallistic() || wind.lgStatic() ) )
	{
		turb2 += POW2(wind.windv0);
	}
#endif

	/* this is average (NOT rms) particle speed for Maxwell distribution, Mihalas 70, 9-70 */
	fixit();  // turbulence was included here for molecules but not ions.  Now neither. Resolve.
	return (realnum)sqrt(8.*BOLTZMANN/PI/ATOMIC_MASS_UNIT*phycon.te/massAMU);
}

STATIC void FillGFF( void )
{

	long i,i1,i2,i3,j,charge,GffMAGIC = 100804;
	double Temp, photon;
	bool lgEOL;

	DEBUG_ENTRY( "FillGFF()" );

#	define chLine_LENGTH 1000
	char chLine[chLine_LENGTH];

	FILE *ioDATA;

	for( i = 0; i < N_TE_GFF; i++ )
	{
		TeGFF[i] = 0.25f*i;
	}
	/* >>chng 06 feb 14, assert thrown at T == 1e10 K, Ryan Porter proposes this fix */
	TeGFF[N_TE_GFF-1] += 0.01f;

	/* number of photon energies */
	N_PHOTON_GFF = (int)( 20. * ( log10(rfield.egamry) - log10(rfield.emm) ) ) + 20;

	PhoGFF = (realnum*)MALLOC( sizeof(realnum)*(unsigned)N_PHOTON_GFF );

	for( i = 0; i< N_PHOTON_GFF; i++ )
	{
		PhoGFF[i] = 0.05f*i + log10(rfield.emm) - 0.5;
	}

	GauntFF = (realnum***)MALLOC( sizeof(realnum**)*(unsigned)(LIMELM+2) );
	for( i = 1; i <= LIMELM; i++ )
	{
		GauntFF[i] = (realnum**)MALLOC( sizeof(realnum*)*(unsigned)N_PHOTON_GFF );
		for( j = 0; j < N_PHOTON_GFF; j++ )
		{
			GauntFF[i][j] = (realnum*)MALLOC( sizeof(realnum)*(unsigned)N_TE_GFF );
		}
	}

	GauntFF_T = (realnum**)MALLOC( sizeof(realnum*)*(unsigned)(LIMELM+2) );
	for( i = 1; i <= LIMELM; i++ )
	{
		GauntFF_T[i] = (realnum*)MALLOC( sizeof(realnum)*(unsigned)N_PHOTON_GFF );
	}

	if( !rfield.lgCompileGauntFF )
	{
		if( trace.lgTrace )
			fprintf( ioQQQ," FillGFF opening gauntff.dat:");

		try
		{
			ioDATA = open_data( "gauntff.dat", "r" );
		}
		catch( cloudy_exit )
		{
			fprintf( ioQQQ, " Defaulting to on-the-fly computation, ");
			fprintf( ioQQQ, "but the code runs much faster if you compile gauntff.dat!\n");
			ioDATA = NULL;
		}

		if( ioDATA == NULL )
		{
			/* Do on the fly computation of Gff's instead.	*/
			for( charge=1; charge<=LIMELM; charge++ )
			{
				for( i=0; i<N_PHOTON_GFF; i++ )
				{
					photon = pow((realnum)10.f,PhoGFF[i]);
					for(j=0; j<N_TE_GFF; j++)
					{
						Temp = pow((realnum)10.f,TeGFF[j]);
						GauntFF[charge][i][j] = (realnum)cont_gaunt_calc( Temp, (double)charge, photon );
					}
				}
			}
		}
		else 
		{
			/* check that magic number is ok */
			if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
			{
				fprintf( ioQQQ, " FillGFF could not read first line of gauntff.dat.\n");
				cdEXIT(EXIT_FAILURE);
			}
			i = 1;
			i1 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
			i2 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
			i3 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);

			if( i1 !=GffMAGIC || i2 !=N_PHOTON_GFF || i3 !=N_TE_GFF )
			{
				fprintf( ioQQQ, 
					" FillGFF: the version of gauntff.dat is not the current version.\n" );
				fprintf( ioQQQ, 
					" FillGFF: I expected to find the numbers  %li %li %li and got %li %li %li instead.\n" ,
					GffMAGIC ,
					N_PHOTON_GFF,
					N_TE_GFF,
					i1 , i2 , i3 );
				fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
				fprintf( ioQQQ, 
					" FillGFF: please recompile the data file with the COMPILE GAUNT command.\n" );
				cdEXIT(EXIT_FAILURE);
			}

			/* now read in the data */
			for( charge = 1; charge <= LIMELM; charge++ )
			{
				for( i = 0; i<N_PHOTON_GFF; i++ )
				{
					/* get next line image */
					if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
					{
						fprintf( ioQQQ, " FillGFF could not read first line of gauntff.dat.\n");
						cdEXIT(EXIT_FAILURE);
					}
					/* each line starts with charge and energy level ( index in rfield ) */
					i3 = 1;
					i1 = (long)FFmtRead(chLine,&i3,INPUT_LINE_LENGTH,&lgEOL);
					i2 = (long)FFmtRead(chLine,&i3,INPUT_LINE_LENGTH,&lgEOL);
					/* check that these numbers are correct */
					if( i1!=charge || i2!=i )
					{
						fprintf( ioQQQ, " FillGFF detected insanity in gauntff.dat.\n");
						fprintf( ioQQQ, 
							" FillGFF: please recompile the data file with the COMPILE GAUNT command.\n" );
						cdEXIT(EXIT_FAILURE);
					}

					/* loop over temperatures to produce array of free free gaunt factors	*/
					for(j = 0; j < N_TE_GFF; j++)
					{
						GauntFF[charge][i][j] = (realnum)FFmtRead(chLine,&i3,chLine_LENGTH,&lgEOL);

						if( lgEOL )
						{
							fprintf( ioQQQ, " FillGFF detected insanity in gauntff.dat.\n");
							fprintf( ioQQQ, 
								" FillGFF: please recompile the data file with the COMPILE GAUNT command.\n" );
							cdEXIT(EXIT_FAILURE);
						}
					}
				}

			}

			/* check that magic number is ok */
			if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
			{
				fprintf( ioQQQ, " FillGFF could not read first line of gauntff.dat.\n");
				cdEXIT(EXIT_FAILURE);
			}
			i = 1;
			i1 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
			i2 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
			i3 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);

			if( i1 !=GffMAGIC || i2 !=N_PHOTON_GFF || i3 !=N_TE_GFF )
			{
				fprintf( ioQQQ, 
					" FillGFF: the version of gauntff.dat is not the current version.\n" );
				fprintf( ioQQQ, 
					" FillGFF: I expected to find the numbers  %li %li %li and got %li %li %li instead.\n" ,
					GffMAGIC ,
					N_PHOTON_GFF,
					N_TE_GFF,
					i1 , i2 , i3 );
				fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
				fprintf( ioQQQ, 
					" FillGFF: please recompile the data file with the COMPILE GAUNT command.\n" );
				cdEXIT(EXIT_FAILURE);
			}

			/* close the data file */
			fclose( ioDATA );
		}
		if( trace.lgTrace )
			fprintf( ioQQQ,"  - opened and read ok.\n");
	}
	else
	{
		/* option to create table of gaunt factors,
		 * executed with the compile gaunt command */
		FILE *ioGFF;

		/*GffMAGIC the magic number for the table of recombination coefficients */
		ioGFF = open_data( "gauntff.dat" , "w", AS_LOCAL_ONLY );
		fprintf(ioGFF,"%li\t%li\t%li\tGFF free-free gaunt factors, created by COMPILE GAUNT command, with %li energy levels and %li temperatures.\n",
			GffMAGIC ,
			N_PHOTON_GFF,
			N_TE_GFF,
			N_PHOTON_GFF,
			N_TE_GFF );

		for( charge = 1; charge <= LIMELM; charge++ )
		{
			for( i=0; i < N_PHOTON_GFF; i++ )
			{
				fprintf(ioGFF, "%li\t%li", charge, i );
				/* loop over temperatures to produce array of gaunt factors	*/
				for(j = 0; j < N_TE_GFF; j++)
				{
					/* Store gaunt factor in N_TE_GFF half dec steps */
					Temp = pow((realnum)10.f,TeGFF[j]);
					photon = pow((realnum)10.f,PhoGFF[i]);
					GauntFF[charge][i][j] = (realnum)cont_gaunt_calc( Temp, (double)charge, photon );
					fprintf(ioGFF, "\t%f", GauntFF[charge][i][j] );
				}
				fprintf(ioGFF, "\n" );
			}
		}

		/* end the file with the same information */
		fprintf(ioGFF,"%li\t%li\t%li\tGFF free-free gaunt factors, created by COMPILE GAUNT command, with %li energy levels and %li temperatures.\n",
			GffMAGIC ,
			N_PHOTON_GFF,
			N_TE_GFF,
			N_PHOTON_GFF,
			N_TE_GFF );

		fclose( ioGFF );

		fprintf( ioQQQ, "FillGFF: compilation complete, gauntff.dat created.\n" );
	}

	lgGffNotFilled = false;

	/* We have already checked the validity of cont_gaunt_calc in sanitycheck.c.
	 * Now we check to see if the InterpolateGff routine also works correctly.	*/
	{
		double gaunt, error;
		double tempsave = phycon.te;
		long logu, loggamma2;

		for( logu=-4; logu<=4; logu++)
		{
			for(loggamma2=-4; loggamma2<=4; loggamma2++)
			{
				double SutherlandGff[9][9]=
				{	{5.5243, 5.5213, 5.4983, 5.3780, 5.0090, 4.4354, 3.8317, 3.2472, 2.7008},
					{4.2581, 4.2577, 4.2403, 4.1307, 3.7816, 3.2436, 2.7008, 2.2126, 1.8041},
					{3.0048, 3.0125, 3.0152, 2.9434, 2.6560, 2.2131, 1.8071, 1.4933, 1.2771},
					{1.8153, 1.8367, 1.8880, 1.9243, 1.7825, 1.5088, 1.2886, 1.1507, 1.0747},
					{0.8531, 0.8815, 0.9698, 1.1699, 1.2939, 1.1988, 1.1033, 1.0501, 1.0237},
					{0.3101, 0.3283, 0.3900, 0.5894, 0.9725, 1.1284, 1.0825, 1.0419, 1.0202},
					{0.1007, 0.1080, 0.1335, 0.2281, 0.5171, 0.9561, 1.1065, 1.0693, 1.0355},
					{0.0320, 0.0344, 0.0432, 0.0772, 0.1997, 0.5146, 0.9548, 1.1042, 1.0680},
					{0.0101, 0.0109, 0.0138, 0.0249, 0.0675, 0.1987, 0.5146, 0.9547, 1.1040}};

				phycon.te = (TE1RYD/pow(10.,(double)loggamma2));
				phycon.alogte = log10(phycon.te);
				double ERyd = pow(10.,(double)(logu-loggamma2));
				if( ERyd > rfield.emm && ERyd < rfield.egamry )
				{
					gaunt = InterpolateGff( 1, ERyd );
					error = fabs( gaunt - SutherlandGff[logu+4][loggamma2+4] ) /gaunt;
					if( error>0.05 )
					{
						fprintf(ioQQQ," PROBLEM DISASTER tfidle found insane gff. log(u) %li, log(gamma2) %li, error %.3e\n",
							logu, loggamma2, error );
						cdEXIT(EXIT_FAILURE);
					}
				}
			}
		}
		phycon.te = tempsave;
		phycon.alogte = log10(phycon.te);
	}

	return;
}

/* Interpolate Gff at some temperature */
STATIC realnum InterpolateGff( long charge , double ERyd )
{
	double GauntAtLowerPho, GauntAtUpperPho;
	static long int ipTe=-1, ipPho=-1;
	double gaunt = 0.;
	double xmin , xmax;
	long i;

	DEBUG_ENTRY( "InterpolateGff()" );

	if( ipTe < 0 )
	{
		/* te totally unknown */
		if( phycon.alogte < TeGFF[0] || phycon.alogte > TeGFF[N_TE_GFF-1] )
		{
			fprintf(ioQQQ," InterpolateGff called with te out of bounds \n");
			cdEXIT(EXIT_FAILURE);
		}
		/* now search for temperature */
		for( i=0; i<N_TE_GFF-1; ++i )
		{
			if( phycon.alogte > TeGFF[i] && phycon.alogte <= TeGFF[i+1] )
			{
				/* found the temperature - use it */
				ipTe = i;
				break;
			}
		}
		ASSERT( (ipTe >=0) && (ipTe < N_TE_GFF-1)  );

	}
	else if( phycon.alogte < TeGFF[ipTe] )
	{
		/* temp is too low, must also lower ipTe */
		ASSERT( phycon.alogte > TeGFF[0] );
		/* decrement ipTe until it is correct */
		while( phycon.alogte < TeGFF[ipTe] && ipTe > 0)
			--ipTe;
	}
	else if( phycon.alogte > TeGFF[ipTe+1] )
	{
		/* temp is too high */
		ASSERT( phycon.alogte <= TeGFF[N_TE_GFF-1] );
		/* increment ipTe until it is correct */
		while( phycon.alogte > TeGFF[ipTe+1] && ipTe < N_TE_GFF-1)
			++ipTe;
	}

	ASSERT( (ipTe >=0) && (ipTe < N_TE_GFF-1)  );

	/* ipTe should now be valid */
	ASSERT( phycon.alogte >= TeGFF[ipTe] && phycon.alogte <= TeGFF[ipTe+1] && ipTe < N_TE_GFF-1 );

	/***************/
	/* This bit is completely analogous to the above, but for the photon vector instead of temp.	*/
	if( ipPho < 0 )
	{
		if( log10(ERyd) < PhoGFF[0] || log10(ERyd) > PhoGFF[N_PHOTON_GFF-1] )
		{
			fprintf(ioQQQ," InterpolateGff called with photon energy out of bounds \n");
			cdEXIT(EXIT_FAILURE);
		}
		for( i=0; i<N_PHOTON_GFF-1; ++i )
		{
			if( log10(ERyd) > PhoGFF[i] && log10(ERyd) <= PhoGFF[i+1] )
			{
				ipPho = i;
				break;
			}
		}
		ASSERT( (ipPho >=0) && (ipPho < N_PHOTON_GFF-1)  );

	}
	else if( log10(ERyd) < PhoGFF[ipPho] )
	{
		ASSERT( log10(ERyd) >= PhoGFF[0] );
		while( log10(ERyd) < PhoGFF[ipPho] && ipPho > 0)
			--ipPho;
	}
	else if( log10(ERyd) > PhoGFF[ipPho+1] )
	{
		ASSERT( log10(ERyd) <= PhoGFF[N_PHOTON_GFF-1] );
		while( log10(ERyd) > PhoGFF[ipPho+1] && ipPho < N_PHOTON_GFF-1)
			++ipPho;
	}
	ASSERT( (ipPho >=0) && (ipPho < N_PHOTON_GFF-1)  );
	ASSERT( log10(ERyd)>=PhoGFF[ipPho] 
		&& log10(ERyd)<=PhoGFF[ipPho+1] && ipPho<N_PHOTON_GFF-1 );

	/* Calculate the answer...must interpolate on two variables.
	 * First interpolate on T, at both the lower and upper photon energies.
	 * Then interpolate between these results for the right photon energy.	*/

	GauntAtLowerPho = (phycon.alogte - TeGFF[ipTe]) / (TeGFF[ipTe+1] - TeGFF[ipTe]) *
		(GauntFF[charge][ipPho][ipTe+1] - GauntFF[charge][ipPho][ipTe]) + GauntFF[charge][ipPho][ipTe];

	GauntAtUpperPho = (phycon.alogte - TeGFF[ipTe]) / (TeGFF[ipTe+1] - TeGFF[ipTe]) *
		(GauntFF[charge][ipPho+1][ipTe+1] - GauntFF[charge][ipPho+1][ipTe]) + GauntFF[charge][ipPho+1][ipTe];

	gaunt = (log10(ERyd) - PhoGFF[ipPho]) / (PhoGFF[ipPho+1] - PhoGFF[ipPho]) * 
		(GauntAtUpperPho - GauntAtLowerPho) + GauntAtLowerPho;

	xmax = MAX4( GauntFF[charge][ipPho][ipTe+1], GauntFF[charge][ipPho+1][ipTe+1],
		GauntFF[charge][ipPho][ipTe], GauntFF[charge][ipPho+1][ipTe] );
	ASSERT( gaunt <= xmax );

	xmin = MIN4( GauntFF[charge][ipPho][ipTe+1], GauntFF[charge][ipPho+1][ipTe+1],
		GauntFF[charge][ipPho][ipTe], GauntFF[charge][ipPho+1][ipTe] );
	ASSERT( gaunt >= xmin );

	ASSERT( gaunt > 0. );

	return (realnum)gaunt;
}

/* Interpolate in table t[lta][ltb], with physical values for the
	 second index given by v[ltb], for values x, and put results in
	 a[lta]; store the index found if that's useful; assumes v[] is
	 sorted */
STATIC int LinterpTable(realnum **t, realnum *v, long int lta, long int ltb, realnum x, realnum *a, long int *pipx)
{
	long int ipx=-1;
	realnum frac;
	long i;

	DEBUG_ENTRY( "LinterpTable()" );

	if(pipx != NULL)
		ipx = *pipx;

	fhunt (v,ltb,x,&ipx); 		/* search for index */
	if(pipx != NULL)
		*pipx = ipx;

	if( ipx == -1 || ipx == ltb )
	{
		return -1;
	}

	ASSERT( (ipx >=0) && (ipx < ltb-1)  );
	ASSERT( x >= v[ipx] && x <= v[ipx+1]);

	frac = (x - v[ipx]) / (v[ipx+1] - v[ipx]);
	for( i=0; i<lta; i++ )
	{
		a[i] = frac*t[i][ipx+1]+(1.f-frac)*t[i][ipx];
	}

	return 0;
}

/* Interpolate in table t[lta][ltb], with physical values for the second index given by v[ltb],
	 for values yy[ny], and put results in a[lta][ly] */
STATIC int LinterpVector(realnum **t, realnum *v, long lta , long ltb, realnum *yy, long ly, realnum **a)
{
	realnum yl, frac;
	long i, j, n;

	DEBUG_ENTRY( "LinterpVector()" );

	if( yy[0] < v[0] || yy[ly-1] > v[ltb-1] )
	{
		return -1;
	}

	n = 0;
	yl = yy[n];
	for(j = 1; j < ltb && n < ly; j++) {
		while (yl < v[j] && n < ly) {
			frac = (yl-v[j-1])/(v[j]-v[j-1]);
			for(i = 0; i < lta; i++)
				a[i][n] = frac*t[i][j]+(1.f-frac)*t[i][j-1];
			n++;
			if(n == ly)
				break;
			ASSERT( yy[n] >= yy[n-1] );
			yl = yy[n];
		}
	}

	return 0;
}
STATIC void fhunt(realnum *xx, long int n, realnum x, long int *j)
{
	/*lint -e731 boolean argument to equal / not equal */
	long int jl, jm, jh, in;
	int up;

	jl = *j;
	up = (xx[n-1] > xx[0]);
	if(jl < 0 || jl >= n) 
	{
		jl = -1;
		jh = n;
	} 
	else 
	{
		in = 1;
		if((x >= xx[jl]) == up) 
		{
			if(jl == n-1) 
			{
				*j = jl;
				return;
			}
			jh = jl + 1;
			while ((x >= xx[jh]) == up)
			{
				jl = jh;
				in += in;
				jh += in;
				if(jh >= n)
				{
					jh = n;
					break;
				}
			}
		}
		else
		{
			if(jl == 0)
			{
				*j = -1;
				return;
			}
			jh = jl--;
			while ((x < xx[jl]) == up)
			{
				jh = jl;
				in += in;
				jl -= in;
				if(jl <= 0)
				{
					jl = 0;
					break;
				}
			}
		}
	}
	while (jh-jl != 1)
	{
		jm = (jh+jl)/2;
		if((x > xx[jm]) == up)
			jl = jm;
		else
			jh = jm;
	}
	*j = jl;
	return;
}
	/*lint +e731 boolean argument to equal / not equal */
