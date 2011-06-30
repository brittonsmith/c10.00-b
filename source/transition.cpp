/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland 
and
 * others.  For conditions of distribution and use see copyright notice in licen
se.txt */

#include "cddefines.h"
#include "transition.h"
#include "version.h"
#include "dense.h"
#include "elementnames.h"
#include "lines.h"
#include "lines_service.h"
#include "opacity.h"
#include "phycon.h"
#include "radius.h"
#include "rfield.h"
#include "rt.h"
#include "taulines.h"
#include "conv.h"

/*outline - adds line photons to reflin and outlin */
/*PutLine enter local line intensity into the intensity stack for eventual printout */
/*PutExtra enter and 'extra' intensity source for some line */
/*DumpLine print various information about an emission line vector, 
 * used in debugging, print to std out, ioQQQ */
/*TexcLine derive excitation temperature of line from contents of line array */
/*TransitionZero zeros out TransitionZero */
/*LineConvRate2CS convert down coll rate back into electron cs in case other parts of code need this for reference */
/*OccupationNumberLine - derive the photon occupation number at line center for any line */
/*MakeCS compute collision strength by g-bar approximations */
/*gbar1 compute g-bar collision strength using Mewe approximations */
/*gbar0 compute g-bar gaunt factor for neutrals */
/*emit_frac returns fraction of populations the produce emission */
/*chIonLbl use information in line array to generate a null terminated ion label in "Fe 2" */
/*chLineLbl use information in line transfer arrays to generate a line label */
/*PutCS enter a collision strength into an individual line vector */

INSTANTIATE_MULTI_ARR( transition, lgBOUNDSCHECKVAL );

/*outline - adds line photons to reflin and outlin */
void transition::outline_resonance( )
{
	bool lgDoChecks = true;
	outline(Emis->ColOvTot, lgDoChecks);
}

/*outline - adds line photons to reflin and outlin */
void transition::outline( double nonScatteredFraction, bool lgDoChecks  )
{
	long int ip = ipCont-1;

	DEBUG_ENTRY( "transition::outline()" );

	if ( 0 && lgDoChecks)
	{
		// nothing to do if very small photon flux or below plasma frequency
		ASSERT( Emis->phots >= 0. );
		if( Emis->phots < SMALLFLOAT || EnergyErg / EN1RYD <= rfield.plsfrq )
			return;

		ASSERT( Emis->FracInwd >= 0. );
		ASSERT( radius.BeamInIn >= 0. );
		ASSERT( radius.BeamInOut >= 0. );
		double Ptot = Emis->Pesc + Emis->Pelec_esc;
		double PhotEmit = Emis->Aul*Ptot*Hi->Pop;
		// do not assert accuracy if close to bounds of fp precision
		double error = MAX2( SMALLFLOAT*1e3 , 3e-1*PhotEmit );
		// do not assert if do not have valid solution
		bool lgGoodSolution = conv.lgConvEden && conv.lgConvIoniz &&
			conv.lgConvPops && conv.lgConvPres && conv.lgConvTemp;
		// see ticket #135 rt inconsistent results
		// this assert trips on a regular basis.  The fix is to rewrite
		// the RT and level population solvers so that the level population
		// and OTS rates are done simultaneously rather than sequentially
		// For now do not throw the assert on a release version - we know
		// about this problem
		ASSERT( t_version::Inst().lgRelease || !lgGoodSolution || 
			Ptot < 1e-3 ||  fp_equal_tol( Emis->phots, PhotEmit , error ));
	}

	bool lgTransStackLine = true;
	outline_base(Emis->dampXvel, Emis->damp, lgTransStackLine, ip, Emis->phots, Emis->FracInwd,
					 nonScatteredFraction);
}

/*emit_frac returns fraction of populations the produce emission */
double emit_frac(const transition *t)
{
	DEBUG_ENTRY( "emit_frac()" );

	ASSERT( t->ipCont > 0 );

	/* collisional deexcitation and destruction by background opacities
	 * are loss of photons without net emission */
	double deexcit_loss = t->Coll.col_str * dense.cdsqte + t->Emis->Aul*t->Emis->Pdest;
	/* this is what is observed */
	double rad_deexcit = t->Emis->Aul*(t->Emis->Pelec_esc + t->Emis->Pesc);
	return rad_deexcit/(deexcit_loss + rad_deexcit);
}

/*DumpLine print various information about an emission line vector, 
 * used in debugging, print to std out, ioQQQ */
void DumpLine(const transition *t)
{
	char chLbl[110];

	DEBUG_ENTRY( "DumpLine()" );

	ASSERT( t->ipCont > 0 );

	/* routine to print contents of line arrays */
	strcpy( chLbl, "DEBUG ");
	strcat( chLbl, chLineLbl(t) );

	fprintf( ioQQQ, 
		"%10.10s Te%.2e eden%.1e CS%.2e Aul%.1e Tex%.2e cool%.1e het%.1e conopc%.1e albdo%.2e\n", 
	  chLbl, 
	  phycon.te, 
	  dense.eden, 
	  t->Coll.col_str, 
	  t->Emis->Aul, 
	  TexcLine(t), 
	  t->Coll.cool, 
	  t->Coll.heat ,
	  opac.opacity_abs[t->ipCont-1],
	  opac.albedo[t->ipCont-1]);

	fprintf( ioQQQ, 
		"Tin%.1e Tout%.1e Esc%.1e eEsc%.1e DesP%.1e Pump%.1e OTS%.1e PopL,U %.1e %.1e PopOpc%.1e\n", 
	  t->Emis->TauIn, 
	  t->Emis->TauTot, 
	  t->Emis->Pesc, 
	  t->Emis->Pelec_esc, 
	  t->Emis->Pdest, 
	  t->Emis->pump, 
	  t->Emis->ots, 
	  t->Lo->Pop, 
	  t->Hi->Pop ,
	  t->Emis->PopOpc );
	return;
}


/*OccupationNumberLine - derive the photon occupation number at line center for any line */
double OccupationNumberLine(const transition *t)
{
	double OccupationNumberLine_v;

	DEBUG_ENTRY( "OccupationNumberLine()" );

	ASSERT( t->ipCont > 0 );

	/* routine to evaluate line photon occupation number - 
	 * return negative number if line is maser */
	if( fabs(t->Emis->PopOpc) > SMALLFLOAT )
	{
		/* the lower population with correction for stimulated emission */
		OccupationNumberLine_v = ( t->Hi->Pop / t->Hi->g ) /
			( t->Emis->PopOpc / t->Lo->g )  *
			(1. - t->Emis->Pesc);
	}
	else
	{
		OccupationNumberLine_v = 0.;
	}
	/* return value is not guaranteed to be positive - negative if
	 * line mases */
	return( OccupationNumberLine_v );
}

/*TexcLine derive excitation temperature of line from contents of line array */
double TexcLine(const transition *t)
{
	double TexcLine_v;

	DEBUG_ENTRY( "TexcLine()" );

	/* routine to evaluate line excitation temp using contents of line array
	 * */
	if( t->Hi->Pop * t->Lo->Pop > 0. )
	{
		TexcLine_v = ( t->Hi->Pop / t->Hi->g )/( t->Lo->Pop / t->Lo->g );
		TexcLine_v = log(TexcLine_v);
		/* protect against infinite temp limit */
		if( fabs(TexcLine_v) > SMALLFLOAT )
		{
			TexcLine_v = - t->EnergyK / TexcLine_v;
		}
	}
	else
	{
		TexcLine_v = 0.;
	}
	return( TexcLine_v );
}

/*chIonLbl use information in line array to generate a null terminated ion label in "Fe 2" */
void chIonLbl(char *chIonLbl_v, const transition *t)
{
	/*static char chIonLbl_v[5];*/

	DEBUG_ENTRY( "chIonLbl()" );

	/* function to use information within the line array
	 * to generate an ion label, giving element and
	 * ionization stage
	 * */
	if( t->Hi->nelem < 0 )
	{
		/* this line is to be ignored */
		if( t->Hi->chLabel[0]=='\0' )
			strcpy( chIonLbl_v, "Dumy" );
		else 
			strcpy( chIonLbl_v, t->Hi->chLabel );
	}
	else if( t->Hi->nelem-1 >= LIMELM )
	{
		// old 12CO and 13CO lines are disabled, we should not be trying to access those strings
		ASSERT( t->Hi->nelem-1 == LIMELM+2 );
		/* >>chng 02 may 15, from to chElementNameShort to go mole right */
		strcpy( chIonLbl_v , elementnames.chElementNameShort[t->Hi->nelem-1] );

		/* chIonStage is four char null terminated string, starting with "_1__" 
		strcat( chIonLbl_v, "CO");*/
	}

	else
	{
		ASSERT( t->Hi->nelem >= 1 );
		/* ElmntSym.chElementSym is null terminated, 2 ch + null, string giving very
		 * short form of element name */
		strcpy( chIonLbl_v , elementnames.chElementSym[t->Hi->nelem -1] );

		/* chIonStage is four char null terminated string, starting with "_1__" */
		strcat( chIonLbl_v, elementnames.chIonStage[t->Hi->IonStg-1]);
	}
	/* chIonLbl is four char null terminated string */
	return/*( chIonLbl_v )*/;
}

/*chLineLbl use information in line transfer arrays to generate a line label */
/* this label is null terminated */
/* ContCreatePointers has test this with full range of wavelengths */
char* chLineLbl(const transition *t)
{
	static char chLineLbl_v[11];

	DEBUG_ENTRY( "chLineLbl()" );

	if( t->Hi->nelem < 1 && t->Hi->IonStg < 1 )
	{
		ASSERT( t->Hi->sp != NULL );
		ASSERT( t->Hi->sp->lgMolecular );
		return ( t->Hi->chLabel );
	}

	ASSERT( t->Hi->nelem >= 1 );
	ASSERT( t->Hi->IonStg >= 1 && t->Hi->IonStg <= t->Hi->nelem + 1 );

	/* function to use information within the line array
	 * to generate a line label, giving element and
	 * ionization stage
	 * rhs are set in large block data */

	/* NB this function is profoundly slow due to sprintf statement
	 * also - it cannot be evaluated within a write statement itself*/

	if( t->WLAng > (realnum)INT_MAX )
	{
		sprintf( chLineLbl_v, "%2.2s%2.2s%5i%c", 
			elementnames.chElementSym[t->Hi->nelem -1], 
			elementnames.chIonStage[t->Hi->IonStg-1], 
		   (int)(t->WLAng/1e8), 'c' );
	}
	else if( t->WLAng > 9999999. )
	{
		/* wavelength is very large, convert to centimeters */
		sprintf( chLineLbl_v, "%2.2s%2.2s%5.2f%c", 
			elementnames.chElementSym[t->Hi->nelem -1], 
			elementnames.chIonStage[t->Hi->IonStg-1], 
			t->WLAng/1e8, 'c' );
	}
	else if( t->WLAng > 999999. )
	{
		/* wavelength is very large, convert to microns */
		sprintf( chLineLbl_v, "%2.2s%2.2s%5i%c", 
			elementnames.chElementSym[t->Hi->nelem -1], 
			elementnames.chIonStage[t->Hi->IonStg-1], 
			(int)(t->WLAng/1e4), 'm' );
	}
	else if( t->WLAng > 99999. )
	{
		/* wavelength is very large, convert to microns */
		sprintf( chLineLbl_v, "%2.2s%2.2s%5.1f%c", 
			elementnames.chElementSym[t->Hi->nelem -1], 
			elementnames.chIonStage[t->Hi->IonStg-1], 
			t->WLAng/1e4, 'm' );
	}
	else if( t->WLAng > 9999. )
	{
		sprintf( chLineLbl_v, "%2.2s%2.2s%5.2f%c", 
			elementnames.chElementSym[ t->Hi->nelem -1], 
			elementnames.chIonStage[t->Hi->IonStg-1], 
		   t->WLAng/1e4, 'm' );
	}
	else if( t->WLAng >= 100. )
	{
		sprintf( chLineLbl_v, "%2.2s%2.2s%5i%c", 
			elementnames.chElementSym[ t->Hi->nelem -1], 
			elementnames.chIonStage[t->Hi->IonStg-1], 
		   (int)t->WLAng, 'A' );
	}
	/* the following two formats should be changed for the
	 * wavelength to get more precision */
	else if( t->WLAng >= 10. )
	{
		sprintf( chLineLbl_v, "%2.2s%2.2s%5.1f%c", 
			elementnames.chElementSym[ t->Hi->nelem -1], 
			elementnames.chIonStage[t->Hi->IonStg-1], 
		   t->WLAng, 'A' );
	}
	else
	{
		sprintf( chLineLbl_v, "%2.2s%2.2s%5.2f%c", 
			elementnames.chElementSym[ t->Hi->nelem -1], 
			elementnames.chIonStage[t->Hi->IonStg-1], 
		   t->WLAng, 'A' );
	}
	/* make sure that string ends with null character - this should
	 * be redundant */
	ASSERT( chLineLbl_v[10]=='\0' );
	return( chLineLbl_v );
}

/*PutCS enter a collision strength into an individual line vector */
void PutCS(double cs, 
  transition *t)
{
	DEBUG_ENTRY( "PutCS()" );

	/* collision strength must be non-negative */
	ASSERT( cs > 0. );

	t->Coll.col_str = (realnum)cs;

	return;
}

void GenerateTransitionConfiguration( const transition *t, char *chComment )
{
	strcpy( chComment, t->Lo->chConfig );
	strcat( chComment, " - " );
	strcat( chComment, t->Hi->chConfig );
	return;
}

static realnum ExtraInten;

/*PutLine enter local line intensity into the intensity stack for eventual printout */
STATIC void PutLine_base(const transition *t, const char *chComment, const char *chLabelTemp, bool lgLabel)
{
	DEBUG_ENTRY( "PutLine_base()" );

	char chLabel[5];
	double xIntensity,
		other,
		xIntensity_in;
		
	/* routine to use line array data to generate input
	 * for emission line array */
	ASSERT( t->ipCont > 0. );
		
	if (lgLabel == true)
	{
		strncpy( chLabel, chLabelTemp, 4 );
		chLabel[4] = '\0';
	}

	/* if ipass=0 then we must generate label info since first pass
	 * gt.0 then only need line intensity data */
	if( LineSave.ipass == 0 )
	{
		if (lgLabel == false)
		{
			/* these variables not used by linadd if ipass ne 0 */
			/* chIonLbl is function that generates a null terminated 4 char string, of form "C  2" */
			chIonLbl(chLabel, t);
		}
		xIntensity = 0.;
	}
	else
	{
		/* both the counting and integrating modes comes here */
		/* not actually used so set to safe value */
		chLabel[0] = '\0';
		/* total line intensity or luminosity 
		 * these may not be defined in initial calls so define here */
		xIntensity = t->Emis->xIntensity + ExtraInten;
	}

	/* initial counting case, where ipass == -1, just ignored above, call linadd below */
	
	/* ExtraInten is option that allows extra intensity (i.e., recomb)
	 * to be added to this line  with Call PutExtra( exta )
	 * in main lines this extra
	 * contribution must be identified explicitly */
	ExtraInten = 0.;
	/*linadd(xIntensity,wl,chLabel,'i');*/
	/*lindst add line with destruction and outward */
	rt.fracin = t->Emis->FracInwd;
	lindst(xIntensity, 
			 t->WLAng, 
			 chLabel, 
			 t->ipCont, 
			 /* this is information only - has been counted in cooling already */
			 'i', 
			 /* do not add to outward beam, also done separately */
			 false,
			 chComment);
	rt.fracin = 0.5;
		
	/* inward part of line - do not move this away from previous lines
	 * since xIntensity is used here */
	xIntensity_in = xIntensity*t->Emis->FracInwd;
	ASSERT( xIntensity_in>=0. );
	linadd(xIntensity_in,t->WLAng,"Inwd",'i',chComment);
	
	/* cooling part of line */
	other = t->Coll.cool;
	linadd(other,t->WLAng,"Coll",'i',chComment);
	
	/* fluorescent excited part of line */
	double radiative_branching;
	enum { lgNEW = true };
	if (lgNEW)
	{
		// Improved two-level version of radiative branching ratio
		const double AulEscp = t->Emis->Aul*(t->Emis->Pesc+t->Emis->Pelec_esc);
		// Would be better to include all outward transition processes from the
		// line, to cater for the general non-two-level case
		const double sinkrate = AulEscp+t->Emis->Aul*t->Emis->Pdest+t->Coll.ColUL;
		if (sinkrate > 0.0) 
		{
			radiative_branching = AulEscp/sinkrate;
		}
		else
		{
			radiative_branching = 0.;
		}
	}
	else
	{
		// This is the excitation ratio not the de-excitation ratio according
		// to its specification
		radiative_branching = (1.-t->Emis->ColOvTot);
	}

	other = t->Lo->Pop * t->Emis->pump * radiative_branching * t->EnergyErg;
	linadd(other,t->WLAng,"Pump",'i',chComment);
		

	/* heating part of line */
	other = t->Coll.heat;
	linadd(other,t->WLAng,"Heat",'i',chComment);
}

/*PutLine enter local line intensity into the intensity stack for eventual printout */
void PutLine(const transition *t, const char *chComment, const char *chLabelTemp)
{	
	const bool lgLabel = true;
	DEBUG_ENTRY( "PutLine()" );
	PutLine_base(t, chComment, chLabelTemp, lgLabel);
}

/*PutLine enter local line intensity into the intensity stack for eventual printout */
void PutLine(const transition *t,
	     const char *chComment)
{
	const bool lgLabel = false;
	const char *chLabelTemp = NULL;
	DEBUG_ENTRY( "PutLine()" );
	PutLine_base(t, chComment, chLabelTemp, lgLabel);
}

/* ==================================================================== */
/*PutExtra enter and 'extra' intensity source for some line */
void PutExtra(double Extra)
{

	DEBUG_ENTRY( "PutExtra()" );

	ExtraInten = (realnum)Extra;
	return;
}

void transition::Junk( void )
{

	DEBUG_ENTRY( "transition::Junk()" );

		/* wavelength, usually in A, used for printout */
	WLAng = -FLT_MAX;

		/* transition energy in degrees kelvin*/
	EnergyK = -FLT_MAX;
		/* transition energy in ergs */
	EnergyErg = -FLT_MAX;
		/* transition energy in wavenumbers */
	EnergyWN = -FLT_MAX;

	/* array offset for radiative transition within continuum array 
	 * is negative if transition is non-radiative. */
	ipCont = -10000;

	CollisionJunk( &Coll );

	/* set these equal to NULL first. That will cause the code to crash if
	 * the variables are ever used before being deliberately set. */
	Emis = &DummyEmis;
	
	Lo = NULL;
	Hi = NULL;
	return;
}

/*TransitionZero zeros out transition array at start of calculation, sets
 * optical depths to initial values */
void transition::Zero( void )
{

	DEBUG_ENTRY( "transition::Zero()" );

	CollisionZero( &Coll );

	StateZero( Lo );
	StateZero( Hi );
	EmLineZero( Emis );

	return;
}

/*LineConvRate2CS convert down coll rate back into electron cs in case other parts of code need this for reference */
void LineConvRate2CS( transition* t , realnum rate )
{

	DEBUG_ENTRY( "LineConvRate2CS()" );

	/* return is collision strength, convert from collision rate from 
	 * upper to lower, this assumes pure electron collisions, but that will
	 * also be assumed by anything that uses cs, for self-consistency */
	t->Coll.col_str = rate * t->Hi->g / (realnum)dense.cdsqte;

	/* change assert to non-negative - there can be cases (Iin H2) where cs has
	 * underflowed to 0 on some platforms */
	ASSERT( t->Coll.col_str >= 0. );
	return;
}

/*gbar0 compute g-bar gaunt factor for neutrals */
STATIC void gbar0(double ex, 
	  realnum *g)
{
	double a, 
	  b, 
	  c, 
	  d, 
	  y;

	DEBUG_ENTRY( "gbar0()" );

	/* written by Dima Verner
	 *
	 * Calculation of the effective Gaunt-factor by use of 
	 * >>refer	gbar	cs	Van Regemorter, H., 1962, ApJ 136, 906
	 * fits for neutrals
	 *  Input parameters: 
	 * ex - energy ryd - now K
	 * t  - temperature in K
	 *  Output parameter:
	 * g  - effective Gaunt factor
	 * */

	/* y = ex*157813.7/te */
	y = ex/phycon.te;
	if( y < 0.01 )
	{
		*g = (realnum)(0.29*(log(1.0+1.0/y) - 0.4/POW2(y + 1.0))/exp(y));
	}
	else
	{
		if( y > 10.0 )
		{
			*g = (realnum)(0.066/sqrt(y));
		}
		else
		{
			a = 1.5819068e-02;
			b = 1.3018207e00;
			c = 2.6896230e-03;
			d = 2.5486007e00;
			d = log(y/c)/d;
			*g = (realnum)(a + b*exp(-0.5*POW2(d)));
		}
	}
	return;
}

/*gbar1 compute g-bar collision strength using Mewe approximations */
STATIC void gbar1(double ex, 
	  realnum *g)
{
	double y;

	DEBUG_ENTRY( "gbar1()" );

	/*	*written by Dima Verner
	 *
	 * Calculation of the effective Gaunt-factor by use of 
	 * >>refer	gbar	cs	Mewe,R., 1972, A&A 20, 215
	 * fits for permitted transitions in ions MgII, CaII, FeII (delta n = 0)
	 * Input parameters: 
	 * ex - excitation energy in Ryd - now K
	 * te  - temperature in K
	 * Output parameter:
	 * g  - effective Gaunt factor
	 */

	/* y = ex*157813.7/te */
	y = ex/phycon.te;
	*g = (realnum)(0.6 + 0.28*(log(1.0+1.0/y) - 0.4/POW2(y + 1.0)));
	return;
}

/*MakeCS compute collision strength by g-bar approximations */
void MakeCS(transition *t)
{
	long int ion;
	double Abun, 
	  cs;
	realnum
	  gbar;

	DEBUG_ENTRY( "MakeCS()" );

	/* routine to get cs from various approximations */

	/* check if abundance greater than 0 */
	ion = t->Hi->IonStg;

	Abun = dense.xIonDense[ t->Hi->nelem -1 ][ ion-1 ];
	if( Abun <= 0. )
	{
		gbar = 1.;
	}
	else
	{
		/* check if neutral or ion */
		if( ion == 1 )
		{
			/* neutral - compute gbar for eventual cs */
			gbar0(t->EnergyK,&gbar);
		}
		else
		{
			/* ion - compute gbar for eventual cs */
			gbar1(t->EnergyK,&gbar);
		}
	}

	/* above was g-bar, convert to cs */
	cs = gbar*(14.5104/WAVNRYD)*t->Emis->gf/t->EnergyWN;

	/* stuff the cs in place */
	t->Coll.col_str = (realnum)cs;
	return;
}

