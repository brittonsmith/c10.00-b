/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*GetGF convert Einstein A into oscillator strength */
/*abscf convert gf into absorption coefficient */
/*RefIndex calculates the index of refraction of air using the line energy in wavenumbers,
 * used to convert vacuum wavelengths to air wavelengths. */
/*eina convert a gf into an Einstein A */
/*WavlenErrorGet - find difference between two wavelengths */
/*linadd enter lines into the line storage array, called once per zone */
/*lindst add local line intensity to line luminosity stack */
/*PntForLine generate pointer for forbidden line */
/*EmLineJunk set all elements of transition struc to dangerous values */
/*EmLineZero set all elements of transition struc to zero */
/*totlin sum total intensity of cooling, recombination, or intensity lines */
/*FndLineHt search through line heat arrays to find the strongest heat source */
/*ConvRate2CS convert down coll rate back into electron cs in case other parts of code need this for reference */
#include "cddefines.h"
#include "lines_service.h"
#include "dense.h"
#include "geometry.h"
#include "hydrogenic.h"
#include "ipoint.h"
#include "iso.h"
#include "lines.h"
#include "opacity.h"
#include "physconst.h"
#include "radius.h"
#include "rfield.h"
#include "rt.h"
#include "taulines.h"

STATIC void normalizeProfile( double *profile, long numPoints, realnum *eArray, double E_centroid );

/*eina convert a gf into an Einstein A */
double eina(double gf,
	  double enercm, 
	  double gup)
{
	double eina_v;

	DEBUG_ENTRY( "eina()" );

	/* derive the transition prob, given the following
	 * call to function is gf, energy in cm^-1, g_up
	 * gf is product of g and oscillator strength
	 * eina = ( gf / 1.499e-8 ) / (wl/1e4)**2 / gup  */
	eina_v = (gf/gup)*TRANS_PROB_CONST*POW2(enercm);
	return( eina_v );
}

/*GetGF convert Einstein A into oscillator strength */
double GetGF(double trans_prob, 
	  double enercm, 
	  double gup)
{
	double GetGF_v;

	DEBUG_ENTRY( "GetGF()" );

	ASSERT( enercm > 0. );
	ASSERT( trans_prob > 0. );
	ASSERT( gup > 0.);

	/* derive the transition prob, given the following
	 * call to function is gf, energy in cm^-1, g_up
	 * gf is product of g and oscillator strength
	 * trans_prob = ( GetGF/gup) / 1.499e-8 / ( 1e4/enercm )**2 */
	GetGF_v = trans_prob*gup/TRANS_PROB_CONST/POW2(enercm);
	return( GetGF_v );
}

/*abscf convert gf into absorption coefficient */
double abscf(double gf, 
	  double enercm, 
	  double gl)
{
	double abscf_v;

	DEBUG_ENTRY( "abscf()" );

	ASSERT(gl > 0. && enercm > 0. && gf > 0. );

	/* derive line absorption coefficient, given the following:
	 * gf, enercm, g_low
	 * gf is product of g and oscillator strength */
	abscf_v = 1.4974e-6*(gf/gl)*(1e4/enercm);
	return( abscf_v );
}

/*RefIndex calculates the index of refraction of air using the line energy in wavenumbers,
 * used to convert vacuum wavelengths to air wavelengths. */
double RefIndex(double EnergyWN )
{
	double RefIndex_v, 
	  WaveMic, 
	  xl, 
	  xn;

	DEBUG_ENTRY( "RefIndex()" );

	ASSERT( EnergyWN > 0. );

	/* the wavelength in microns */
	WaveMic = 1.e4/EnergyWN;

	/* only do index of refraction if longward of 2000A */
	if( WaveMic > 0.2 )
	{
		/* longward of 2000A
		 * xl is 1/WaveMic^2 */
		xl = 1.0/WaveMic/WaveMic;
		/* use a formula from 
		 *>>refer	air	index refraction	Allen, C.W. 1973, Astrophysical quantities, 
		 *>>refercon	3rd Edition (AQ), p.124 */
		xn = 255.4/(41. - xl);
		xn += 29498.1/(146. - xl);
		xn += 64.328;
		RefIndex_v = xn/1.e6 + 1.;
	}
	else
	{
		RefIndex_v = 1.;
	}
	ASSERT( RefIndex_v >= 1. );
	return( RefIndex_v );
}

/*WavlenErrorGet - given the real wavelength in A for a line
 * routine will find the error expected between the real 
 * wavelength and the wavelength printed in the output, with 4 sig figs,
 * function returns difference between exact and 4 sig fig wl, so 
 * we have found correct line is fabs(d wl) < return */
realnum WavlenErrorGet( realnum wavelength )
{
	double a;
	realnum errorwave;

	DEBUG_ENTRY( "WavlenErrorGet()" );

	ASSERT( LineSave.sig_figs <= 6 );

	if( wavelength > 0. )
	{
		/* normal case, positive (non zero) wavelength */
		a = log10( wavelength+FLT_EPSILON);
		a = floor(a);
	}
	else
	{
		/* might be called with wl of zero, this is that case */
		/* errorwave = 1e-4f; */
		a = 0.;
	}

	errorwave = 5.f * (realnum)pow( 10., a - (double)LineSave.sig_figs );
	return errorwave;
}

/*linadd enter lines into the line storage array, called once per zone for each line*/
STATIC void lincom(
  double xInten,	/* xInten - local emissivity per unit vol, no fill fac */
  realnum wavelength,	/* realnum wavelength */
  const char *chLab,/* string label for ion */
  // ipnt offset of line in continuum mesh
  long int ipnt, 
  char chInfo,		/* character type of entry for line - given below */
					/* 'c' cooling, 'h' heating, 'i' info only, 'r' recom line */
  // *chComment string explaining line 
  const char *chComment,
  // lgAdd says whether we've come in via linadd (true) or lindst (false)
  bool lgAdd)
{
	DEBUG_ENTRY( "lincom()" );

	/* main routine to actually enter lines into the line storage array
	 * called at top level within routine lines
	 * called series of times in routine PutLine for lines transferred
	 */
	
	/* three values, -1 is just counting, 0 if init, 1 for calculation */
	if( LineSave.ipass > 0 )
	{
		/* >>chng 06 feb 08, add test on xInten positive, no need to evaluate
		 * for majority of zero */
		if (lgAdd || xInten > 0.)
		{
			/* not first pass, sum lines only
			 * total emission from vol */
			/* LineSave.ipass > 0, integration across simulation, sum lines only 
			 * emissivity, emission per unit vol, for this zone */
			LineSv[LineSave.nsum].SumLine[0] += xInten*radius.dVeffAper;
			/* local emissivity in line */
			/* integrated intensity or luminosity, the emissivity times the volume */
			LineSv[LineSave.nsum].emslin[0] = xInten;
		}

		if (lgAdd)
		{
			if (wavelength > 0 )
			{
				/* no need to increment or set [1] version since this is called with no continuum
				 * index, don't know what to do */
				/* only put informational lines, like "Q(H) 4861", in this stack
				 * many continua have a wavelength of zero and are proper intensities,
				 * it would be wrong to predict their transferred intensity */
				LineSv[LineSave.nsum].emslin[1] = LineSv[LineSave.nsum].emslin[0];
				LineSv[LineSave.nsum].SumLine[1] = LineSv[LineSave.nsum].SumLine[0]; 
			}
		}
		else
		{
			if ( xInten > 0. && ipnt <= rfield.nflux )
			{
				/* emergent_line accounts for destruction by absorption outside
				 * the line-forming region */
				const double saveemis = emergent_line( 
					xInten*rt.fracin , xInten*(1.-rt.fracin) , ipnt );
				LineSv[LineSave.nsum].emslin[1] = saveemis;
				LineSv[LineSave.nsum].SumLine[1] += saveemis*radius.dVeffAper; 
			}
		}
	}
		
	else if( LineSave.ipass == 0 )
	{
		/* first call to stuff lines in array, confirm that label is one of
		 * the four correct ones */
		ASSERT( (chInfo == 'c') || (chInfo == 'h') || (chInfo == 'i') || (chInfo == 'r' ) );
		/* then save it into array */
		LineSv[LineSave.nsum].chSumTyp = (char)chInfo;
		LineSv[LineSave.nsum].emslin[0] = 0.;
		LineSv[LineSave.nsum].emslin[1] = 0.;
		LineSv[LineSave.nsum].chComment = chComment;
		/* check that null is correct, string overruns have 
		 * been a problem in the past */
		ASSERT( strlen( chLab )<5 );
		strcpy( LineSv[LineSave.nsum].chALab, chLab );
		
		if (lgAdd)
		{			
			LineSv[LineSave.nsum].wavelength = wavelength;			
		}
		else
		{
			// number of lines ok, set parameters for first pass 
			// negative wavelengh means it is just label, possibly not correct
			LineSv[LineSave.nsum].wavelength = fabs(wavelength);
			LineSv[LineSave.nsum].SumLine[0] = 0.;
			LineSv[LineSave.nsum].SumLine[1] = 0.;
			
			// check that line wavelength and continuum index agree to some extent
			// this check cannot be very precise because some lines have 
			// "wavelengths" that are set by common usage rather than the correct
			// wavelength derived from energy and index of refraction of air
			ASSERT( ipnt > 0 );
#		ifndef NDEBUG		
			double error = MAX2(0.1*rfield.AnuOrg[ipnt-1] , rfield.widflx[ipnt-1] );
			ASSERT( wavelength<=0 ||
					  fabs( rfield.AnuOrg[ipnt-1] - RYDLAM / wavelength) < error );
#		endif
		}
	}
		
	/* increment the line counter */
	++LineSave.nsum;
	
	/* routine can be called with negative LineSave.ipass, in this case
	 * we are just counting number of lines for current setup */
}

/*linadd enter lines into the line storage array, called once per zone for each line*/
void linadd(
  double xInten,	/* xInten - local emissivity per unit vol, no fill fac */
  realnum wavelength,	/* realnum wavelength */
  const char *chLab,/* string label for ion */
  char chInfo,		/* character type of entry for line - given below */
					/* 'c' cooling, 'h' heating, 'i' info only, 'r' recom line */
  const char *chComment )
{
	DEBUG_ENTRY( "linadd()" );
	
	// Values added to get common interface with lindst
	const long int ipnt = LONG_MAX;
	
	lincom( xInten, wavelength, chLab, ipnt, chInfo, chComment, true );
}


/*emergent_line find emission from surface of cloud after correcting for
 * extinction due to continuous opacity for inward & outward directed emission */
double emergent_line( 
	/* emiemission in inward direction */
	double emissivity_in , 
	/* emission in outward direction */
	double emissivity_out , 
	/* array index for continuum frequency on fortran scale */
	long int ipCont )
{

	double emergent_in , emergent_out;
	long int i = ipCont-1;

	DEBUG_ENTRY( "emergent_line()" );

	ASSERT( i >= 0 && i < rfield.nupper-1 );

	/* do nothing if first iteration since we do not know the outward-looking
	 * optical depths.  In version C07.02.00 we assumed an infinite optical
	 * depth in the outward direction, which would be appropriate for a 
	 * HII region on the surface of a molecular cloud.  This converged onto
	 * the correct solution in later iterations, but on the first iteration
	 * this underestimated total emission if the infinite cloud were not
	 * present.  With C07.02.01 we make no assuptions about what is in the
	 * outward direction and simply use the local emission. 
	 * Behavior is unchanged on later iterations */
	if( iteration == 1 )
	{
		/* first iteration - do not know outer optical depths so assume very large 
		 * optical depths */
		emergent_in = emissivity_in*opac.E2TauAbsFace[i];
		emergent_out = emissivity_out;
	}
	else
	{
		if( geometry.lgSphere )
		{
			/* second or later iteration in closed or spherical geometry */
			/* inwardly directed emission must get to central hole then across entire
			 * far side of shell */
			emergent_in = emissivity_in  * opac.E2TauAbsFace[i] *opac.E2TauAbsTotal[i];

			/* E2 is outwardly directed emission to get to outer edge of cloud */
			emergent_out = emissivity_out * opac.E2TauAbsOut[i];
		}
		else
		{
			/* open geometry in second or later iteration, outer optical depths are known 
			 * this is light emitted into the outer direction and backscattered
			 * into the inner */
			double reflected = emissivity_out * opac.albedo[i] * (1.-opac.E2TauAbsOut[i]);
			/* E2 is to get to central hole */
			emergent_in = (emissivity_in + reflected) * opac.E2TauAbsFace[i];
			/* E2 is to get to outer edge */
			emergent_out = (emissivity_out - reflected) * opac.E2TauAbsOut[i];
		}
	}
	/* return the net emission that makes it to the surface */
	return( emergent_in + emergent_out );
}

/* outline_base - calls outline_base_bin after deciding whether to add impulse or resolved line */
void outline_base(double dampXvel, double damp, bool lgTransStackLine, long int ip, double phots, realnum inwd,
						double nonScatteredFraction)
{
	DEBUG_ENTRY( "outline_base()" );

#define DO_PROFILE false 

	if( !DO_PROFILE )
		outline_base_bin(lgTransStackLine, ip, phots, inwd, nonScatteredFraction);
	else
	{
		ASSERT( damp > 0. );
		double LineWidth = dampXvel/damp;
		LineWidth = MIN2( 0.1 * SPEEDLIGHT, LineWidth );
		double sigma = (LineWidth/SPEEDLIGHT);
		long ip3SigmaRed = ipoint( MAX2( rfield.emm, rfield.anu[ip] - 3.*sigma*rfield.anu[ip] ) );
		long ip3SigmaBlue = ipoint( MIN2( rfield.egamry, rfield.anu[ip] + 3.*sigma*rfield.anu[ip] ) );
		ASSERT( ip3SigmaBlue >= ip3SigmaRed );
		long numBins = ip3SigmaBlue - ip3SigmaRed + 1;

		if( numBins < 3 )
			outline_base_bin(lgTransStackLine, ip, phots, inwd, nonScatteredFraction);
		else
		{
			valarray<double> profile(numBins);

			for( long ipBin=ip3SigmaRed; ipBin<=ip3SigmaBlue; ipBin++ )
			{
				double x = (rfield.anu[ip] - rfield.anu[ipBin])/rfield.anu[ip]/sigma;
				profile[ipBin-ip3SigmaRed] = vfun( damp, x );
			}

			normalizeProfile( &profile[0], numBins, rfield.anu+ip3SigmaRed, rfield.anu[ip] );

			for( long ipBin=ip3SigmaRed; ipBin<=ip3SigmaBlue; ipBin++ )
				outline_base_bin(lgTransStackLine, ipBin, phots*profile[ipBin-ip3SigmaRed], inwd, nonScatteredFraction);
		}
	}
}

/*outline_base_bin - adds line photons to bins of reflin and outlin */
void outline_base_bin(bool lgTransStackLine, long int ip, double phots, realnum inwd,
						double nonScatteredFraction)
{
	DEBUG_ENTRY( "outline_base_bin()" );

	if (lgTransStackLine)
	{
		rfield.DiffuseLineEmission[ip] += 
			(realnum)phots;

		/* the reflected part */
		rfield.reflin[0][ip] +=
			(realnum)(inwd*phots*radius.BeamInIn);
		
		/* inward beam that goes out since sphere set */
		rfield.outlin[0][ip] +=
			(realnum)(inwd*phots*radius.BeamInOut*opac.tmn[ip]*nonScatteredFraction);
		
		/* outward part */
		rfield.outlin[0][ip] +=
			(realnum)((1.-inwd)*phots*radius.BeamOutOut*opac.tmn[ip]*nonScatteredFraction);
	}
	else
	{
		rfield.reflin[0][ip] +=
			(realnum)(phots*radius.dVolReflec);

		rfield.outlin[0][ip] +=
			(realnum)(phots*radius.dVolOutwrd*opac.ExpZone[ip]);
	}
}

STATIC void normalizeProfile( double *profile, long numPoints, realnum *eArray, double E_centroid )
{
	double psum = 0.;
	double pEsum = 0.;
	ASSERT( numPoints >= 1 );
	if( numPoints==1 )
	{
		profile[0] = 1.;
		return;
	}

	for( long i=0; i<numPoints; i++ )
	{
		pEsum += eArray[i] * profile[i];
		psum += profile[i];
	}

	for( long i=0; i<numPoints; i++ )
	{
		ASSERT( pEsum > SMALLFLOAT );
		profile[i] *= E_centroid/pEsum;
	}

	for( long i=0; i<numPoints; i++ )
	{
		pEsum += eArray[i] * profile[i];
		psum += profile[i];
	}

	return;
}

/*lindst add line with destruction and outward */
void lindst(
  // xInten - local emissivity per unit vol
  double xInten, 
  // wavelength of line in Angstroms
  realnum wavelength, 
  // *chLab string label for ion
  const char *chLab, 
  // ipnt offset of line in continuum mesh
  long int ipnt, 
  // chInfo character type of entry for line - 'c' cooling, 'h' heating, 'i' info only, 'r' recom line
  char chInfo, 
  // lgOutToo should line be included in outward beam?
  bool lgOutToo,
  // *chComment string explaining line 
  const char *chComment )
{
	DEBUG_ENTRY( "lindst()" );

	// do not add information lines to outward beam
	ASSERT( !lgOutToo || chInfo!='i' );

	lincom(xInten, wavelength, chLab, ipnt, chInfo, chComment, false );

	if( LineSave.ipass > 0 )
	{
		/* >>chng 06 feb 08, add test on xInten positive, no need to evaluate
		 * for majority of zero */
		if (lgOutToo && xInten > 0.)
		{
			/* add line to outward beam 
			 * there are lots of lines that are sums of other lines, or
			 * just for info of some sort.  These have flag lgOutToo false.
			 * Note that the EnergyRyd variable only has a rational
			 * value if PntForLine was called just before this routine - in all
			 * cases where this did not happen the flag is false. */
			const bool lgTransStackLine = false;
			const long int ip = ipnt - 1;
			const double phots = xInten/(rfield.anu[ipnt-1]*EN1RYD);
			const realnum inwd = (realnum)(1.0-(1.+geometry.covrt)/2.);
			const double nonScatteredFraction = 1.;

			outline_base_bin(lgTransStackLine, ip, phots, inwd, nonScatteredFraction);
		}
	}
}

/*lindst add line with destruction and outward */
void lindst(
  double dampXvel,
  double damp,
  // xInten - local emissivity per unit vol
  double xInten,
  // wavelength of line in Angstroms
  realnum wavelength,
  // *chLab string label for ion
  const char *chLab,
  // ipnt offset of line in continuum mesh
  long int ipnt,
  // chInfo character type of entry for line - 'c' cooling, 'h' heating, 'i' info only, 'r' recom line
  char chInfo,
  // lgOutToo should line be included in outward beam?
  bool lgOutToo,
  // *chComment string explaining line
  const char *chComment )
{
	DEBUG_ENTRY( "lindst()" );

	// do not add information lines to outward beam
	ASSERT( !lgOutToo || chInfo!='i' );

	lincom(xInten, wavelength, chLab, ipnt, chInfo, chComment, false );

	if( LineSave.ipass > 0 )
	{
		/* >>chng 06 feb 08, add test on xInten positive, no need to evaluate
		 * for majority of zero */
		if (lgOutToo && xInten > 0.)
		{
			/* add line to outward beam
			 * there are lots of lines that are sums of other lines, or
			 * just for info of some sort.  These have flag lgOutToo false.
			 * Note that the EnergyRyd variable only has a rational
			 * value if PntForLine was called just before this routine - in all
			 * cases where this did not happen the flag is false. */
			const bool lgTransStackLine = false;
			const long int ip = ipnt - 1;
			const double phots = xInten/(rfield.anu[ipnt-1]*EN1RYD);
			const realnum inwd = (realnum)(1.0-(1.+geometry.covrt)/2.);
			const double nonScatteredFraction = 1.;

			outline_base(dampXvel, damp, lgTransStackLine, ip, phots, inwd, nonScatteredFraction);
		}
	}
}

/*lindst add line with destruction and outward */
void lindst(
	transition *t,
  // *chLab string label for ion
  const char *chLab, 
  // chInfo character type of entry for line - 'c' cooling, 'h' heating, 'i' info only, 'r' recom line
  char chInfo, 
  // lgOutToo should line be included in outward beam?
  bool lgOutToo,
  // *chComment string explaining line 
  const char *chComment )
{
	DEBUG_ENTRY( "lindst()" );

	lindst(  t->Emis->dampXvel, t->Emis->damp, t->Emis->xIntensity, t->WLAng, chLab, t->ipCont, chInfo,
			 lgOutToo, chComment );

}

/*PntForLine generate pointer for forbidden line */
void PntForLine(
  /* wavelength of transition in Angstroms */
  double wavelength, 
  /* label for this line */
  const char *chLabel,
  /* this is array index on the f, not c scale,
   * for the continuum cell holding the line */
  long int *ipnt)
{
	/* 
	 * maximum number of forbidden lines - this is a good bet since
	 * new lines do not go into this group, and lines are slowly 
	 * moving to level 1 
	 */
	const int MAXFORLIN = 1000;
	static long int ipForLin[MAXFORLIN]={0};

	/* number of forbidden lines entered into continuum array */
	static long int nForLin;

	DEBUG_ENTRY( "PntForLine()" );

	/* must be 0 or greater */
	ASSERT( wavelength >= 0. );

	if( wavelength == 0. )
	{
		/* zero is special flag to initialize */
		nForLin = 0;
	}
	else
	{

		if( LineSave.ipass > 0 )
		{
			/* not first pass, sum lines only */
			*ipnt = ipForLin[nForLin];
		}
		else if( LineSave.ipass == 0 )
		{
			/* check if number of lines in arrays exceeded */
			if( nForLin >= MAXFORLIN )
			{
				fprintf( ioQQQ, "PROBLEM %5ld lines is too many for PntForLine.\n", 
				  nForLin );
				fprintf( ioQQQ, " Increase the value of maxForLine everywhere in the code.\n" );
				cdEXIT(EXIT_FAILURE);
			}

			/* ipLineEnergy will only put in line label if nothing already there */
			const double EnergyRyd = RYDLAM/wavelength;
			ipForLin[nForLin] = ipLineEnergy(EnergyRyd,chLabel , 0);
			*ipnt = ipForLin[nForLin];
		}
		else
		{
			/* this is case where we are only counting lines */
			*ipnt = 0;
		}
		++nForLin;
	}
	return;
}

/*EmLineJunk set all elements of transition struc to dangerous values */
void EmLineJunk( emission *t )
{

	DEBUG_ENTRY( "EmLineJunk()" );

	/* optical depth in continuum to ill face */
	t->TauCon = -FLT_MAX;

	/* inward and total line optical depths */
	t->TauIn = -FLT_MAX;
	t->TauTot = -FLT_MAX;

	/* type of redistribution function, */
	t->iRedisFun = INT_MIN;

	/* array offset for line within fine opacity */
	t->ipFine = -10000;

	/* inward fraction */
	t->FracInwd = -FLT_MAX;

	/* continuum pumping rate */
	t->pump = -FLT_MAX;

	/* line intensity */
	t->xIntensity = -FLT_MAX;

	/* number of photons emitted per sec in the line */
	t->phots = -FLT_MAX;

	/* gf value */
	t->gf = -FLT_MAX;

	/* escape and destruction probs */
	t->Pesc = -FLT_MAX;
	t->Pdest = -FLT_MAX;
	t->Pelec_esc = -FLT_MAX;

	/* damping constant, and number related to it */
	t->dampXvel = -FLT_MAX;
	t->damp = -FLT_MAX;

	/* ratio of collisional to radiative excitation*/
	t->ColOvTot = -FLT_MAX;

	/* auto-ionization fraction */
	t->AutoIonizFrac = -FLT_MAX;

	/* line opacity */
	t->opacity = -FLT_MAX;

	t->PopOpc = -FLT_MAX;

	/* transition prob, Einstein A upper to lower */
	t->Aul = -FLT_MAX;

	/* ots rate */
	t->ots = -FLT_MAX;
	return;
}

/*CollisionJunk set all elements of transition struc to dangerous values */
void CollisionJunk( collision * t )
{

	DEBUG_ENTRY( "CollisionJunk()" );

	/** collision rate coefficient, [cm^3 s-1], upper to lower */
	t->ColUL = -FLT_MAX;

	/* Coll->cooling and Coll->heating due to collisional excitation */
	t->cool = -FLT_MAX;
	t->heat = -FLT_MAX;

	/* collision strengths for transition */
	t->col_str = -FLT_MAX;

	for( long i=0; i<ipNCOLLIDER; i++ )
		t->col_stri[i] = -FLT_MAX;

	return;
}

/*StateJunk set all elements of transition struc to dangerous values */
void StateJunk( quantumState * t )
{

	DEBUG_ENTRY( "StateJunk()" );

	t->chLabel[0] = '\0';

	/** statistical weight [dimensionless] */
	t->g = -FLT_MAX;

	/** population of state [cm-3] */
	t->Pop = -FLT_MAX;

	 /** ion stage of element, 1 for atom, 2 ion, etc */
	t->IonStg = -10000;

	 /** atomic number of element, 1 for H, 2 for He, etc */
	t->nelem = -10000;

	/* species - initially point to NULL */
	t->sp = NULL;

	return;
}

/*EmLineZero zeros out the emission line structure */
void EmLineZero( emission * t )
{

	DEBUG_ENTRY( "EmLineZero()" );

	/* total optical depth in all overlapping lines to illuminated face,
	 * used for pumping */
	t->TauCon = opac.taumin;

	/* inward and total line optical depths */
	/* >>chng 03 feb 14, from 0 to opac.taumin */
	t->TauIn = opac.taumin;

	/* total optical depths */
	t->TauTot = 1e20f;

	/* inward fraction */
	/* >>chng 03 feb 14, from 0 to 1 */
	t->FracInwd = 1.;

	/* continuum pumping rate */
	t->pump = 0.;

	/* line intensity */
	t->xIntensity = 0.;

	/* number of photons emitted per sec in the line */
	t->phots = 0.;

	/* escape and destruction probs */
	/* >>chng 03 feb 14, change from 0 to 1 */
	t->Pesc = 1.;
	t->Pdest = 0.;
	t->Pelec_esc = 0.;

	/* ratio of collisional to radiative excitation*/
	t->ColOvTot = 1.;

	/* pop that enters net opacity */
	t->PopOpc = 0.;

	/* ots rate */
	t->ots = 0.;
	return;
}

/*CollisionZero zeros out the structure */
void CollisionZero( collision * t )
{

	DEBUG_ENTRY( "CollisionZero()" );

	t->ColUL = 0.;
	/* Coll->cooling and Coll->heating due to collisional excitation */
	t->cool = 0.;
	t->heat = 0.;
	return;
}

/*StateZero zeros out the structure */
void StateZero( quantumState * t )
{

	DEBUG_ENTRY( "StateZero()" );

	/** population of state [cm-3] */
	t->Pop = 0.;
	return;
}

/*ConvRate2CS convert down coll rate back into electron cs in case other parts of code need this for reference */
double ConvRate2CS( realnum gHi , realnum rate )
{

	double cs;

	DEBUG_ENTRY( "ConvRate2CS()" );

	/* return is collision strength, convert from collision rate from 
	 * upper to lower, this assumes pure electron collisions, but that will
	 * also be assumed by anything that uses cs, for self-consistency */
	cs = rate * gHi / dense.cdsqte;

	/* change assert to non-negative - there can be cases (Iin H2) where cs has
	 * underflowed to 0 on some platforms */
	ASSERT( cs >= 0. );
	return cs;
}

/*ConvCrossSect2CollStr convert collisional deexcitation cross section for into collision strength */
double ConvCrossSect2CollStr( double CrsSectCM2, double gLo, double E_ProjectileRyd, double reduced_mass_grams )
{
	double CollisionStrength;

	DEBUG_ENTRY( "ConvCrossSect2CollStr()" );

	ASSERT( CrsSectCM2 >= 0. );
	ASSERT( gLo >= 0. );
	ASSERT( E_ProjectileRyd >= 0. );
	ASSERT( reduced_mass_grams >= 0. );

	CollisionStrength = CrsSectCM2 * gLo * E_ProjectileRyd / (PI*BOHR_RADIUS_CM*BOHR_RADIUS_CM);

	// this part is being tested.
#if 0
	CollisionStrength *= reduced_mass_grams / ELECTRON_MASS;
#endif

	ASSERT( CollisionStrength >= 0. );
	return CollisionStrength;
}

/*totlin sum total intensity of cooling, recombination, or intensity lines */
double totlin(
	/* chInfor is 1 char, 
	'i' information, 
	'r' recombination or 
	'c' collision */
	int chInfo)
{
	long int i;
	double totlin_v;

	DEBUG_ENTRY( "totlin()" );

	/* routine goes through set of entered line
	 * intensities and picks out those which have
	 * types agreeing with chInfo.  Valid types are
	 * 'c', 'r', and 'i'
	 *begin sanity check */
	if( (chInfo != 'i' && chInfo != 'r') && chInfo != 'c' )
	{
		fprintf( ioQQQ, " TOTLIN does not understand chInfo=%c\n", 
		  chInfo );
		cdEXIT(EXIT_FAILURE);
	}
	/*end sanity check */

	/* now find sum of lines of type chInfo */
	totlin_v = 0.;
	for( i=0; i < LineSave.nsum; i++ )
	{
		if( LineSv[i].chSumTyp == chInfo )
		{
			totlin_v += LineSv[i].SumLine[0];
		}
	}
	return( totlin_v );
}


/*FndLineHt search through line heat arrays to find the strongest heat source */
void FndLineHt(long int *level, 
  /* this is the index of the strongest line in the array on the c scale */
  long int *ipStrong, 
  double *Strong)
{
	long int i; 

	DEBUG_ENTRY( "FndLineHt()" );

	*Strong = 0.;
	*level = 0;

	/* do the level 1 lines, 0 is dummy line, <=nLevel1 is correct for c scale */
	for( i=1; i <= nLevel1; i++ )
	{
		/* check if a line was the major heat agent */
		if( TauLines[i].Coll.heat > *Strong )
		{
			*ipStrong = i;
			*level = 1;
			*Strong = TauLines[i].Coll.heat;
		}
	}

	/* now do the level 2 lines */
	for( i=0; i < nWindLine; i++ )
	{
		if( TauLine2[i].Hi->IonStg < TauLine2[i].Hi->nelem+1-NISO )
		{
			/* check if a line was the major heat agent */
			if( TauLine2[i].Coll.heat > *Strong )
			{
				*ipStrong = i;
				*level = 2;
				*Strong = TauLine2[i].Coll.heat;
			}
		}
	}

	/* now do the hyperfine structure lines */
	for( i=0; i < nHFLines; i++ )
	{
		/* check if a line was the major heat agent */
		if( HFLines[i].Coll.heat > *Strong )
		{
			*ipStrong = i;
			*level = 3;
			*Strong = HFLines[i].Coll.heat;
		}
	}

	/* lines from external databases */
	for( i=0; i <linesAdded2; i++)
	{
		/* check if a line was the major heat agent */
		if( dBaseLines[i].tran->Coll.heat > *Strong )
		{
			*ipStrong = i;
			*level = 4;
			*Strong = dBaseLines[i].tran->Coll.heat;
		}
	}

	fixit();  // all other line stacks need to be included here.
	// can we just sweep over line stack?  Is that ready yet?

	return;
}

quantumState *AddState2Stack( void )
{
	DEBUG_ENTRY( "AddState2Stack()" );

	ASSERT( !lgStatesAdded );

	currentState = new quantumState;

	StateJunk( currentState );

	if( statesAdded == 0 )
	{
		GenericStates = currentState;
		GenericStates->next = NULL;
		lastState = GenericStates;
	}
	else
	{
		StateZero( currentState );
		lastState->next = currentState;
		lastState = lastState->next;
	}

	statesAdded++;

	return currentState;
}

emission *AddLine2Stack( bool lgRadiativeTrans )
{
	DEBUG_ENTRY( "AddLine2Stack()" );

	if( !lgRadiativeTrans )
	{
		return &DummyEmis;
	}
	else
	{
		ASSERT( lgLinesAdded == false );

		currentLine = new emission;

		EmLineJunk( currentLine );

		if( linesAdded == 0 )
		{
			GenericLines = currentLine;
			GenericLines->next = NULL;
			lastLine = GenericLines;
		}
		else
		{
			/* 
			\todo 2 Does doing EmLineZero here defeat the purpose of EmLineJunk? 
			* maybe we should pass full set of Emis components, fill everything in 
			* here, and THEN use EmLineZero?  */
			EmLineZero( currentLine );

			lastLine->next = currentLine;
			lastLine = lastLine->next;
		}

		linesAdded++;
		return currentLine;
	}
}
