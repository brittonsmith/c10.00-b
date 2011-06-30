/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ContCreateMesh calls fill to set up continuum energy mesh if first call, 
 * otherwise reset to original mesh */
/*fill define the continuum energy grid over a specified range */
/*ChckFill perform sanity check confirming that the energy array has been properly filled */
/*rfield_opac_malloc MALLOC space for opacity arrays */
/*read_continuum_mesh read the continuum definition from the file continuum_mesh.ini */
#include "cddefines.h"
#include "rfield.h"
#include "iterations.h"
#include "physconst.h"
#include "dense.h"
#include "trace.h"
#include "opacity.h"
#include "ipoint.h"
#include "geometry.h"
#include "continuum.h"

/* read the continuum definition from the file continuum_mesh.ini */
STATIC void read_continuum_mesh( void );

/*fill define the continuum energy grid over a specified range */
STATIC void fill(double fenlo, 
  double fenhi, 
  double resolv, 
  long int *n0, 
  long int *ipnt,
  /* this says only count, do not fill */
  bool lgCount );

/*rfield_opac_malloc MALLOC space for opacity arrays */
STATIC void rfield_opac_malloc(void);

/*ChckFill perform sanity check confirming that the energy array has been properly filled */
STATIC void ChckFill(void);

void ContCreateMesh(void)
{
	long int 
	  i, 
	  ipnt, 
	  n0;

	/* flag to say whether pointers have ever been evaluated */
	static bool lgPntEval = false;

	DEBUG_ENTRY( "ContCreateMesh()" );

	/* lgPntEval is local static variable defined false when defined. 
	 * it is set true below, so that pointers only created one time in the
	 * history of this coreload. */
	if( lgPntEval )
	{
		if( trace.lgTrace )
		{
			fprintf( ioQQQ, " ContCreateMesh called, not evaluating.\n" );
		}
		/* now save current form of energy array */
		for( i=0; i < rfield.nupper; i++ )
		{
			rfield.anu[i] = rfield.AnuOrg[i];
			rfield.anu2[i] = rfield.anu[i]*rfield.anu[i];
		}
		return;
	}
	else
	{
		if( trace.lgTrace )
		{
			fprintf( ioQQQ, " ContCreateMesh called first time.\n" );
		}
		lgPntEval = true;
	}

	/* set size of arrays that continuum iteration information */

	/* read in the continuum mesh resolution definition */
	/* >>chng 01 sep 29, add external file "continuum_mesh.ini" with fill parameters */
	read_continuum_mesh();

	/* fill in continuum with freq points
	 * arg are range pointer, 2 energy limits, resolution
	 * first argument is lower energy of the continuum range to be defined
	 * second argument is upper range; both are in Rydbergs
	 * third number is the relative energy resolution, dnu/nu,
	 * for that range of the continuum
	 * last two numbers are internal book keeping
	 * N0 is number of energy cells used so far
	 * IPNT is a counter for the number of fills used
	 *
	 * if this is changed, then also change warning in GetTable about using
	 * transmitted continuum - it says version number where continuum changed
	 * */
	n0 = 1;
	ipnt = 0;
	/* this is number of ranges that will be introduced below*/
	continuum.nrange = 0;

	/* ================================================================ */
	/* NB - this block must agree exactly with the one that follows     */
	n0 = 1;
	ipnt = 0;
	/* this is number of ranges that will be introduced below*/
	continuum.nrange = 0;
	continuum.StoredEnergy[continuum.nStoredBands-1] = rfield.egamry;

	fill(rfield.emm, continuum.StoredEnergy[0] , continuum.StoredResolution[0],&n0,&ipnt,true );
	for(i=1; i<continuum.nStoredBands; ++i )
	{
		fill(continuum.StoredEnergy[i-1]      , 
			continuum.StoredEnergy[i],  
			continuum.StoredResolution[i],
			&n0,&ipnt,true);
	}
	/* ================================================================ */

	/* at this point debugger shows that anu and widflx are defined 
	 * up through n0-2 - due to c offset -1 from Fortran! */
	rfield.nupper = n0 - 1;
	/* there must be a cell above nflux for us to pass unity through the vol integrator */
	if( rfield.nupper >= NCELL )
	{
		fprintf(ioQQQ," Currently the arrays that hold interpolated tables can only hold %i points.\n",NCELL);
		fprintf(ioQQQ," This continuum mesh really needs to have %li points.\n",rfield.nupper);
		fprintf(ioQQQ," Please increase the value of NCELL in rfield.h and recompile.\n Sorry.");
		cdEXIT(EXIT_FAILURE);
	}

	/*>>chng 04 oct 10, from nflux = nupper to nupper-1 since vectors must malloc to nupper, but
	 * will address [nflux] for unit continuum test */
	rfield.nflux = rfield.nupper-1;

	/* allocate space for continuum arrays within rfield.h and opacity arrays in opacity.h
	 * sets lgRfieldMalloced true */
	rfield_opac_malloc();

	/* geometry.nend_max is largest number of zones needed on any iteration,
	 * will use it to malloc arrays that save source function as function of zone */
	/* now change all limits, for all iterations, to this value */
	geometry.nend_max = geometry.nend[0];
	for( i=1; i < iterations.iter_malloc; i++ )
	{
		geometry.nend_max = MAX2( geometry.nend_max , geometry.nend[i] );
	}
	/* nend_max+1 because search phase is zone 0, first zone at illumin face is 1 */
	rfield.ConEmitLocal = (realnum**)MALLOC( (size_t)(geometry.nend_max+1)*sizeof(realnum *)  );
	rfield.ConSourceFcnLocal = (realnum**)MALLOC( (size_t)(geometry.nend_max+1)*sizeof(realnum *)  );
	for( i=0;i<(geometry.nend_max+1); ++i )
	{
		rfield.ConEmitLocal[i] = (realnum*)MALLOC( (size_t)rfield.nupper*sizeof(realnum)  );
		rfield.ConSourceFcnLocal[i] = (realnum*)MALLOC( (size_t)rfield.nupper*sizeof(realnum)  );
	}
	for( i=0;i<(geometry.nend_max+1); ++i )
	{
		for( long j=0; j<rfield.nupper; ++j)
		{
			rfield.ConSourceFcnLocal[i][j] = 1.;
		}
	}

	/* ================================================================ */
	n0 = 1;
	ipnt = 0;

	/* this is number of ranges that will be introduced below*/
	continuum.nrange = 0;

	/* the default array values are set in continuum_mesh.ini */
	fill(rfield.emm, continuum.StoredEnergy[0] , continuum.StoredResolution[0] , 
		&n0,&ipnt,false);
	for(i=1; i<continuum.nStoredBands; ++i )
	{
		fill(continuum.StoredEnergy[i-1]      , 
			continuum.StoredEnergy[i],  
			continuum.StoredResolution[i],
			&n0,&ipnt,false);
	}

	/* ================================================================ */

	/* fill in the false highest cell used for unit verification */
	rfield.widflx[rfield.nupper] = rfield.widflx[rfield.nupper-1];
	rfield.anu[rfield.nupper] = rfield.anu[rfield.nupper-1] + 
		rfield.widflx[rfield.nupper];

	/* there must be a cell above nflux for us to pass unity through the vol integrator
	 * as a sanity check.  assert that this is true so we will crash if ever changed 
	ASSERT( rfield.nupper +1 <= rfield.nupper );*/

	/* this is done here when the space is first allocated, 
	 * then done on every subsequent initialization in zero.c */
	rfield_opac_zero( 0 , rfield.nupper );

	/* this is a sanity check for results produced above by fill */
	ChckFill();

	/* now fix widflx array so that it is correct */
	for( i=1; i<rfield.nupper-1; ++i )
	{
		rfield.widflx[i] = ((rfield.anu[i+1] - rfield.anu[i]) + (rfield.anu[i] - 
			rfield.anu[i-1]))/2.f;
	}

	ipnt = 0;
	/* now save current form of array, and define some quantities related to it */
	for( i=0; i < rfield.nupper; i++ )
	{
		double alf , bet;

		rfield.AnuOrg[i] = rfield.anu[i];
		rfield.anusqr[i] = (realnum)sqrt(rfield.AnuOrg[i]);
		/* following are Compton exchange factors from Tarter */
		/* this code also appears in highen, but coef needed before that routine called. */
		alf = 1./(1. + rfield.anu[i]*(1.1792e-4 + 7.084e-10*rfield.anu[i]));
		bet = 1. - alf*rfield.anu[i]*(1.1792e-4 + 2.*7.084e-10*rfield.anu[i])/4.;
		rfield.csigh[i] = (realnum)(alf*rfield.anu[i]*rfield.anu[i]*3.858e-25);
		rfield.csigc[i] = (realnum)(alf*bet*rfield.anu[i]*3.858e-25);
		rfield.anu2[i] = rfield.anu[i]*rfield.anu[i];

		/* >>chng 05 feb 28, add transmission and mapping coef */
		/* map these coarse continua into fine continuum grid */
		if( rfield.anu[i] < rfield.fine_ener_lo || rfield.anu[i]>rfield.fine_ener_hi )
		{
			/* 0 (false) says not defined */
			rfield.ipnt_coarse_2_fine[i] = 0;
		}
		else
		{
			if( ipnt==0 )
			{
				/* this is the first one that maps onto the fine array */
				rfield.ipnt_coarse_2_fine[i] = 0;
				ipnt = 1;
			}
			else
			{
				/* find first fine frequency that is greater than this coarse value */
				while (ipnt < rfield.nfine_malloc && rfield.fine_anu[ipnt]<rfield.anu[i] )
				{
					++ipnt;
				}
				rfield.ipnt_coarse_2_fine[i] = ipnt;
			}
		}
		/*fprintf(ioQQQ," coarse %li nu= %.3e points to fine %li nu=%.3e\n",
			i, rfield.anu[i] , rfield.ipnt_coarse_2_fine[i] , rfield.fine_anu[rfield.ipnt_coarse_2_fine[i]] );*/
		rfield.trans_coef_zone[i] = 1.f;
		rfield.trans_coef_total[i] = 1.f;
	}
	return;
}

/*fill define the continuum energy grid over a specified range, called by ContCreateMesh */
STATIC void fill(
  /* lower bounds to this energy range */
  double fenlo, 
  /* upper bounds to this continuum range */
  double fenhi, 
  /* relative energy resolution */
  double resolv, 
  /* starting index within frequency grid */
  long int *n0, 
  /* which energy band this is */
  long int *ipnt,  
  /* this says only count, do not fill */
  bool lgCount )
{
	long int i, 
	  nbin;
	realnum widtot;
	double aaa , bbb;

	DEBUG_ENTRY( "fill()" );

	ASSERT( fenlo>0. && fenhi>0. && resolv>0. );

	/* this is the number of cells needed to fill the array with numbers at the requested resolution */
	nbin = (long int)(log(10.)*log10(fenhi/fenlo)/resolv + 1);

	if( lgCount )
	{
		/* true means only count number of cells, don't do anything */
		*n0 += nbin;
		return;
	}

	if( *ipnt > 0 && fabs(1.-fenlo/continuum.filbnd[*ipnt]) > 1e-4 )
	{
		fprintf( ioQQQ, " FILL improper bounds.\n" );
		fprintf( ioQQQ, " ipnt=%3ld fenlo=%11.4e filbnd(ipnt)=%11.4e\n", 
		  *ipnt, fenlo, continuum.filbnd[*ipnt] );
		cdEXIT(EXIT_FAILURE);
	}

	ASSERT( *ipnt < continuum.nStoredBands );

	continuum.ifill0[*ipnt] = *n0 - 1;
	continuum.filbnd[*ipnt] = (realnum)fenlo;
	continuum.filbnd[*ipnt+1] = (realnum)fenhi;

	/* this is the number of cells needed to fill the array with numbers 
	nbin = (long int)(log(10.)*log10(fenhi/fenlo)/resolv + 1);*/
	continuum.fildel[*ipnt] = (realnum)(log10(fenhi/fenlo)/nbin);

	if( continuum.fildel[*ipnt] < 0.01 )
	{
		continuum.filres[*ipnt] = (realnum)(log(10.)*continuum.fildel[*ipnt]);
	}
	else
	{
		continuum.filres[*ipnt] = (realnum)((pow(10.,2.*continuum.fildel[*ipnt]) - 1.)/2./
			pow((realnum)10.f,continuum.fildel[*ipnt]));
	}

	if( (*n0 + nbin-2) > rfield.nupper )
	{
		fprintf( ioQQQ, " Fill would need %ld cells to get to an energy of %.3e\n", 
		  *n0 + nbin, fenhi );
		fprintf( ioQQQ, " This is a major logical error in fill.\n");
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	widtot = 0.;
	for( i=0; i < nbin; i++ )
	{
		bbb = continuum.fildel[*ipnt]*((realnum)(i) + 0.5);
		aaa = pow( 10. , bbb );

		rfield.anu[i+continuum.ifill0[*ipnt]] = (realnum)(fenlo*aaa);

		rfield.widflx[i+continuum.ifill0[*ipnt]] = rfield.anu[i+continuum.ifill0[*ipnt]]*
		  continuum.filres[*ipnt];

		widtot += rfield.widflx[i+continuum.ifill0[*ipnt]];
	}

	*n0 += nbin;
	if( trace.lgTrace && (trace.lgConBug || trace.lgPtrace) )
	{
		fprintf( ioQQQ, 
			" FILL range%2ld from%10.3e to%10.3eR in%4ld cell; ener res=%10.3e WIDTOT=%10.3e\n", 
		  *ipnt, 
		  rfield.anu[continuum.ifill0[*ipnt]] - rfield.widflx[continuum.ifill0[*ipnt]]/2., 
		  rfield.anu[continuum.ifill0[*ipnt]+nbin-1] + rfield.widflx[continuum.ifill0[*ipnt]+nbin-1]/2., 
		  nbin, 
		  continuum.filres[*ipnt], 
		  widtot );

		fprintf( ioQQQ, " The requested range was%10.3e%10.3e The requested resolution was%10.3e\n", 
		  fenlo, fenhi, resolv );
	}

	/* nrange is number of ranges  */
	*ipnt += 1;
	continuum.nrange = MAX2(continuum.nrange,*ipnt);
	return;
}

/*ChckFill perform sanity check confirming that the energy array has been properly filled */
STATIC void ChckFill(void)
{
	bool lgFail;
	long int i, 
	  ipnt;
	double energy;

	DEBUG_ENTRY( "ChckFill()" );

	ASSERT( rfield.anu[0] >= rfield.emm*0.99 );
	ASSERT( rfield.anu[rfield.nupper-1] <= rfield.egamry*1.01 );

	lgFail = false;
	for( i=0; i < continuum.nrange; i++ )
	{
		/* test middle of energy bound */
		energy = (continuum.filbnd[i] + continuum.filbnd[i+1])/2.;
		ipnt = ipoint(energy);
		if( energy < rfield.anu[ipnt-1] - rfield.widflx[ipnt-1]*0.5 )
		{
			fprintf( ioQQQ, " ChckFill middle test low fail\n" );
			lgFail = true;
		}

		/* >>chng 02 jul 16, add second test - when "set resol 10" used,
		 * very large values of cell width, combined with fact that cells
		 * are log increasing, causes problem.  */
		else if( (energy > rfield.anu[ipnt-1] + rfield.widflx[ipnt-1]*0.5) &&
			( energy > rfield.anu[ipnt] - rfield.widflx[ipnt]*0.5 ) )
		{
			fprintf( ioQQQ, " ChckFill middle test high fail\n" );
			lgFail = true;
		}

		/* test near low bound */
		energy = continuum.filbnd[i]*0.99 + continuum.filbnd[i+1]*0.01;
		ipnt = ipoint(energy);
		if( energy < rfield.anu[ipnt-1] - rfield.widflx[ipnt-1]*0.5 )
		{
			fprintf( ioQQQ, " ChckFill low test low fail\n" );
			lgFail = true;
		}

		else if( energy > rfield.anu[ipnt-1] + rfield.widflx[ipnt-1]* 0.5 )
		{
			fprintf( ioQQQ, " ChckFill low test high fail\n" );
			lgFail = true;
		}

		/* test near high bound */
		energy = continuum.filbnd[i]*0.01 + continuum.filbnd[i+1]*0.99;
		ipnt = ipoint(energy);

		if( energy < rfield.anu[ipnt-1] - rfield.widflx[ipnt-1]*0.5 )
		{
			fprintf( ioQQQ, " ChckFill high test low fail\n" );
			lgFail = true;
		}
		/* >>chng 02 jul 16, add second test - when "set resol 10" used,
		 * very large values of cell width, combined with fact that cells
		 * are log increasing, causes problem.  */
		else if( (energy > rfield.anu[ipnt-1] + rfield.widflx[ipnt-1]*0.5) &&
			( energy > rfield.anu[ipnt] - rfield.widflx[ipnt]*0.5 ) )
		{
			fprintf( ioQQQ, " ChckFill high test high fail\n" );
			lgFail = true;
		}
	}

	if( lgFail )
	{
		cdEXIT(EXIT_FAILURE);
	}
	return;
}

/* MALLOC arrays within rfield */
STATIC void rfield_opac_malloc(void)
{
	long i;

	DEBUG_ENTRY( "rfield_opac_malloc()" );

	/* allocate one more than we use for the unit integration,
	 * will back up at end of routine */
	++rfield.nupper;

	/* >>chng 03 feb 12, add fine mesh fine grid fine opacity array to keep track of line overlap */
	/** \todo	3	consider making the fine opacity array a double.  with a float, the opacity 
	 * itself often becomes a denormalized number, it then becomes significant when multiplied
	 * by dr - can cause numerical noise.  this is why the coarse opacity array is a double */

	/* frequency range in Rydberg needed for all resonance lines */
	rfield.fine_ener_lo = rfield.emm;
	rfield.fine_ener_hi = 1500.f;

	/* set resolution of fine continuum mesh. 
	 * rfield.fine_opac_velocity_width is width per cell, cm/s 
	 * choose width so that most massive species (usually Fe) is well resolved
	 * 
	 * rfield.fine_opac_nelem is the most massive (hence sharpest line)
	 * we will worry about.  By default this is iron but can be changed
	 * with SET FINE CONTINUUM command 
	 * 
	 * TeLowestFineOpacity of 1e4 K is temperature were line width is 
	 * evaluated.  Tests were done using the stop temperature in its place
	 * Te below 1e4 K made fine opacity grid huge 
	 * do not let temp get higher than 1e4 either - code run with stop temp 10 set
	 * stop temp of 1e10K and assert thrown at line 204 of cont_createpointers.c 
	 * simply use 1e4 K as a characteristic temperature */
	/** \todo	1	set temp of 1e4K will be too coarse a line for PDRs where
	 * H2 line overlap is very important */
	double TeLowestFineOpacity = 1e4;
	rfield.fine_opac_velocity_width = 
		(realnum)sqrt(2.*BOLTZMANN/ATOMIC_MASS_UNIT*TeLowestFineOpacity/
		dense.AtomicWeight[rfield.fine_opac_nelem] ) / 
		/* we want fine_opac_nresolv continuum elements across this line
		 * default is 1, changed with SET FINE CONTINUUM command */
		rfield.fine_opac_nresolv;

	/* we are at first zone so velocity shift is zero */
	rfield.ipFineConVelShift = 0;

	/* dimensionless resolution, dE/E, this is used in ipoint to get offset in find mesh */
	rfield.fine_resol = rfield.fine_opac_velocity_width / SPEEDLIGHT;

	/* the number of cells needed */
	rfield.nfine_malloc = (long)(log10( rfield.fine_ener_hi / rfield.fine_ener_lo ) / log10( 1. + rfield.fine_resol ) );
	if( rfield.nfine_malloc <= 0 )
		TotalInsanity();
	rfield.nfine = rfield.nfine_malloc;

	/* this is the fine opacity array to ghost the main low-resolution array */
	rfield.fine_opac_zone = (realnum *)MALLOC(sizeof(realnum)*(unsigned)rfield.nfine_malloc );
	memset(rfield.fine_opac_zone , 0 , (unsigned long)rfield.nfine_malloc*sizeof(realnum) );

	/* this is the fine total optical array to ghost the main low-resolution array */
	rfield.fine_opt_depth = (realnum *)MALLOC(sizeof(realnum)*(unsigned)rfield.nfine_malloc );
	memset(rfield.fine_opt_depth , 0 , (unsigned long)rfield.nfine_malloc*sizeof(realnum) );

	rfield.fine_anu = (realnum *)MALLOC(sizeof(realnum)*(unsigned)rfield.nfine_malloc );

	/* now fill in energy array */
	ASSERT( rfield.fine_ener_lo > 0. && rfield.fine_resol > 0 );
	for( i=0;i<rfield.nfine_malloc; ++i )
	{
		rfield.fine_anu[i] = rfield.fine_ener_lo * (realnum)pow( (1.+rfield.fine_resol), (i+1.) );
	}
	/* done with fine array */

	/* used to count number of lines per cell */
	rfield.line_count = (long *)MALLOC(sizeof(long)*(unsigned)NCELL );
	for( i=0; i<rfield.nupper; ++i)
	{
		rfield.line_count[i] = 0;
	}
	rfield.anu = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.AnuOrg = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.widflx = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.anulog = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.anusqr = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.anu2 = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.anu3 = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.flux_beam_time = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.flux_isotropic = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.flux_beam_const = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.flux_accum = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.ExtinguishFactor = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.convoc = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.OccNumbBremsCont = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.OccNumbIncidCont = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.OccNumbDiffCont = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.OccNumbContEmitOut = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.ConInterOut = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.SummedCon = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.SummedDif = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.SummedDifSave = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.SummedOcc = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.ConOTS_local_photons = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.DiffuseEscape = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.TotDiff2Pht = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.ConOTS_local_OTS_rate = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.otslin = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.otscon = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.outlin_noplot = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.flux_beam_const_save = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.flux_time_beam_save = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.flux_isotropic_save = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.trans_coef_zone = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.trans_coef_total = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.DiffuseLineEmission = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.ipnt_coarse_2_fine = (long int*)MALLOC((size_t)(rfield.nupper*sizeof(long int)) );

	/* possibly save cumulative flux */
	rfield.flux = (realnum**)MALLOC((size_t)(2*sizeof(realnum*)) );
	rfield.ConEmitReflec = (realnum**)MALLOC((size_t)(2*sizeof(realnum*)) );
	rfield.ConEmitOut = (realnum**)MALLOC((size_t)(2*sizeof(realnum*)) );
	rfield.ConRefIncid = (realnum**)MALLOC((size_t)(2*sizeof(realnum*)) );
	rfield.flux_total_incident = (realnum**)MALLOC((size_t)(2*sizeof(realnum*)) );
	rfield.reflin = (realnum**)MALLOC((size_t)(2*sizeof(realnum*)) );
	rfield.outlin = (realnum**)MALLOC((size_t)(2*sizeof(realnum*)) );

	for( i=0; i<2; ++i )
	{
		rfield.flux[i] = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
		rfield.ConEmitReflec[i] = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
		rfield.ConEmitOut[i] = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
		rfield.ConRefIncid[i] = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
		rfield.flux_total_incident[i] = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
		rfield.reflin[i] = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
		rfield.outlin[i] = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	}
	// the cumulative (time integral) emission
	memset(rfield.flux[1] , 0 , (unsigned long)rfield.nupper*sizeof(realnum) );
	memset(rfield.ConEmitReflec[1] , 0 , (unsigned long)rfield.nupper*sizeof(realnum) );
	memset(rfield.ConEmitOut[1] , 0 , (unsigned long)rfield.nupper*sizeof(realnum) );
	memset(rfield.ConRefIncid[1] , 0 , (unsigned long)rfield.nupper*sizeof(realnum) );
	memset(rfield.flux_total_incident[1] , 0 , (unsigned long)rfield.nupper*sizeof(realnum) );
	memset(rfield.reflin[1] , 0 , (unsigned long)rfield.nupper*sizeof(realnum) );
	memset(rfield.outlin[1] , 0 , (unsigned long)rfield.nupper*sizeof(realnum) );

	/* chng 02 may 16, by Ryan...added array for gaunt factors for ALL charges, malloc here.	*/
	/* First index is EFFECTIVE CHARGE MINUS ONE!	*/
	rfield.gff = (realnum**)MALLOC((size_t)((LIMELM+1)*sizeof(realnum*)) );
	for( i = 1; i <= LIMELM; i++ )
	{
		rfield.gff[i] = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	}

	rfield.csigh = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	rfield.csigc = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );

	rfield.comdn = (double*)MALLOC((size_t)(rfield.nupper*sizeof(double)) );
	rfield.comup = (double*)MALLOC((size_t)(rfield.nupper*sizeof(double)) );
	rfield.ContBoltz = (double*)MALLOC((size_t)(rfield.nupper*sizeof(double)) );

	/*realnum rfield.otssav[NC_ELL][2];*/
	rfield.otssav = (realnum**)MALLOC((size_t)(rfield.nupper*sizeof(realnum*)));
	for( i=0; i<rfield.nupper; ++i)
	{
		rfield.otssav[i] = (realnum*)MALLOC(2*sizeof(realnum));
	}


	/* char rfield.chLineLabel[NLINES][5];*/
	rfield.chLineLabel = (char**)MALLOC((size_t)(rfield.nupper*sizeof(char*)));
	rfield.chContLabel = (char**)MALLOC((size_t)(rfield.nupper*sizeof(char*)));

	/* now allocate all the labels for each of the above lines */
	for( i=0; i<rfield.nupper; ++i)
	{
		rfield.chLineLabel[i] = (char*)MALLOC(5*sizeof(char));
		rfield.chContLabel[i] = (char*)MALLOC(5*sizeof(char));
	}

	opac.TauAbsFace = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	memset( opac.TauAbsFace , 0 , rfield.nupper*sizeof(realnum) );

	opac.TauScatFace = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	opac.E2TauAbsFace = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	opac.E2TauAbsTotal = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	opac.TauAbsTotal = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	opac.E2TauAbsOut = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	opac.ExpmTau = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );
	opac.tmn = (realnum*)MALLOC((size_t)(rfield.nupper*sizeof(realnum)) );

	opac.opacity_abs = (double*)MALLOC((size_t)(rfield.nupper*sizeof(double)) );
	opac.opacity_abs_savzon1 = (double*)MALLOC((size_t)(rfield.nupper*sizeof(double)) );
	opac.OldOpacSave = (double*)MALLOC((size_t)(rfield.nupper*sizeof(double)) );
	opac.opacity_sct = (double*)MALLOC((size_t)(rfield.nupper*sizeof(double)) );
	opac.albedo = (double*)MALLOC((size_t)(rfield.nupper*sizeof(double)) );
	opac.opacity_sct_savzon1 = (double*)MALLOC((size_t)(rfield.nupper*sizeof(double)) );
	opac.OpacStatic = (double*)MALLOC((size_t)(rfield.nupper*sizeof(double)) );
	opac.FreeFreeOpacity = (double*)MALLOC((size_t)(rfield.nupper*sizeof(double)) );
	opac.ExpZone = (double*)MALLOC((size_t)(rfield.nupper*sizeof(double)) );

	opac.TauAbsGeo = (realnum**)MALLOC((size_t)(2*sizeof(realnum *)) );
	opac.TauScatGeo = (realnum**)MALLOC((size_t)(2*sizeof(realnum *)) );
	opac.TauTotalGeo = (realnum**)MALLOC((size_t)(2*sizeof(realnum *)) );

	for( i=0; i<2; ++i)
	{
		opac.TauAbsGeo[i] = (realnum*)MALLOC(rfield.nupper*sizeof(realnum));
		opac.TauScatGeo[i] = (realnum*)MALLOC(rfield.nupper*sizeof(realnum));
		opac.TauTotalGeo[i] = (realnum*)MALLOC(rfield.nupper*sizeof(realnum));
	}

	/* fix allocate trick for one more than we use for the unit integration */
	--rfield.nupper;

	/* say that space exists */
	lgRfieldMalloced = true;
	return;
}


/* read the continuum definition from the file continuum_mesh.ini */
STATIC void read_continuum_mesh( void )
{
	FILE *ioDATA;
	char chLine[INPUT_LINE_LENGTH];
	long i;
	bool lgEOL;
	long i1 , i2 , i3;

	DEBUG_ENTRY( "read_continuum_mesh()" );

	if( trace.lgTrace )
		fprintf( ioQQQ," read_continuum_mesh opening continuum_mesh.ini:");

	ioDATA = open_data( "continuum_mesh.ini", "r" );

	/* first line is a version number and does not count */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " read_continuum_mesh could not read first line of continuum_mesh.ini.\n");
		cdEXIT(EXIT_FAILURE);
	}
	/* count how many lines are in the file, ignoring all lines
	 * starting with '#' */
	continuum.nStoredBands = 0;
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
	{
		/* we want to count the lines that do not start with #
		 * since these contain data */
		if( chLine[0] != '#')
			++continuum.nStoredBands;
	}

	/* we now have number of lines containing pairs of bounds,
	 * allocate space for the arrays we will need */
	continuum.filbnd = 
		((realnum *)MALLOC( (size_t)(continuum.nStoredBands+1)*sizeof(realnum )));
	continuum.fildel = 
		((realnum *)MALLOC( (size_t)(continuum.nStoredBands+1)*sizeof(realnum )));
	continuum.filres = 
		((realnum *)MALLOC( (size_t)(continuum.nStoredBands+1)*sizeof(realnum )));
	continuum.ifill0 = 
		((long *)MALLOC( (size_t)(continuum.nStoredBands+1)*sizeof(long )));
	continuum.StoredEnergy = 
		((double *)MALLOC( (size_t)(continuum.nStoredBands+1)*sizeof(double )));
	continuum.StoredResolution = 
		((double *)MALLOC( (size_t)(continuum.nStoredBands+1)*sizeof(double )));

	/* now rewind the file so we can read it a second time*/
	if( fseek( ioDATA , 0 , SEEK_SET ) != 0 )
	{
		fprintf( ioQQQ, " read_continuum_mesh could not rewind continuum_mesh.ini.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* check that magic number is ok */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " read_continuum_mesh could not read first line of continuum_mesh.ini.\n");
		cdEXIT(EXIT_FAILURE);
	}

	i = 1;
	/* continuum mesh magic number */
	i1 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
	i2 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
	i3 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);

	bool lgResPower;

	/* the following is the set of numbers that appear at the start of continuum_mesh.ini */
	if( i1 == 1 && i2 == 9 && i3 == 29 )
		// old version of the file (c08 and older), this has pairs: upper limit freq range, resolution
		// this format is still supported to accomodate users with existing continuum_mesh.ini files.
		lgResPower = false;
	else if( i1 == 10 && i2 == 8 && i3 == 8 )
		// new version of the file (c10 and newer), this has pairs: upper limit freq range, resolving power
		// resolving power = 1./resolution
		lgResPower = true;
	else
	{
		fprintf( ioQQQ, 
			" read_continuum_mesh: the version of continuum_mesh.ini is not supported.\n" );
		fprintf( ioQQQ, 
			" I found version number %li %li %li.\n" ,
			 i1 , i2 , i3 );
		fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
		cdEXIT(EXIT_FAILURE);
	}

	/* this starts at 1 not 0 since zero is reserved for the
	 * dummy line */
	continuum.nStoredBands = 0;
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
	{
		/* only look at lines without '#' in first col */
		if( chLine[0] != '#')
		{
			i = 1;
			continuum.StoredEnergy[continuum.nStoredBands] = FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
			continuum.StoredResolution[continuum.nStoredBands] = FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);

			// continuum energy could be 0 to indicate low or high energy bounds of code
			// but none can be negative
			if( continuum.StoredEnergy[continuum.nStoredBands]<0. ||
					continuum.StoredResolution[continuum.nStoredBands]<=0. )
			{
				fprintf(ioQQQ, "DISASTER PROBLEM continuum_mesh.ini has a non-positive number.\n");
				cdEXIT(EXIT_FAILURE);
			}

			// convert resolving power (entered quantity) into resolution
			if( lgResPower )
				continuum.StoredResolution[continuum.nStoredBands] = 1./
					continuum.StoredResolution[continuum.nStoredBands];

			/* this is option to rescale resolution with set resolution command */
			continuum.StoredResolution[continuum.nStoredBands] *= continuum.ResolutionScaleFactor;

			++continuum.nStoredBands;
		}
	}

	fclose( ioDATA );

	/* now verify continuum grid is ok - first are all values but the last positive? */
	for( i=1; i<continuum.nStoredBands-1; ++i )
	{
		if( continuum.StoredEnergy[i-1]>=continuum.StoredEnergy[i] )
		{
			fprintf( ioQQQ, 
				" read_continuum_mesh: The continuum definition array energies must be in increasing order.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}
	if( continuum.StoredEnergy[continuum.nStoredBands-1]!=0 )
	{
		fprintf( ioQQQ, 
			" read_continuum_mesh: The last continuum array energies must be zero.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return;
}

/*rfield_opac_zero zero out rfield arrays between certain limits */
void rfield_opac_zero( 
					  /* index for first element in arrays to be set to zero */
					  long lo , 
					  /* array index for highest element to be set */
					  long ihi )
{
	long int i;

	/* >>chng 01 aug 19, space not allocated yet,
	* following code must also be present in contcreatemesh where
	* space allocated for the first time */
	if( lgRfieldMalloced )
	{
		unsigned long n=(unsigned long)(ihi-lo+1);
		memset(&rfield.OccNumbDiffCont[lo]      , 0 , n*sizeof(realnum) );
		memset(&rfield.OccNumbContEmitOut[lo]   , 0 , n*sizeof(realnum) );
		memset(&rfield.ContBoltz[lo]            , 0 , n*sizeof(double) );
		/*>>chng 06 aug 15, this is now 2D array, saving diffuse continuum
		* over all zones for use in exact RT */
		/*memset(&rfield.ConEmitLocal[lo]         , 0 , n*sizeof(realnum) );*/
		memset(&rfield.ConEmitReflec[0][lo]        , 0 , n*sizeof(realnum) );
		memset(&rfield.ConEmitOut[0][lo]           , 0 , n*sizeof(realnum) );
		memset(&rfield.reflin[0][lo]               , 0 , n*sizeof(realnum) );
		memset(&rfield.ConRefIncid[0][lo]          , 0 , n*sizeof(realnum) );
		memset(&rfield.SummedCon[lo]            , 0 , n*sizeof(realnum) );
		memset(&rfield.OccNumbBremsCont[lo]     , 0 , n*sizeof(realnum) );
		memset(&rfield.convoc[lo]               , 0 , n*sizeof(realnum) );
		memset(&rfield.flux[0][lo]                 , 0 , n*sizeof(realnum) );
		memset(&rfield.flux_total_incident[0][lo]  , 0 , n*sizeof(realnum) );
		memset(&rfield.flux_beam_const_save[lo] , 0 , n*sizeof(realnum) );
		memset(&rfield.flux_time_beam_save[lo]  , 0 , n*sizeof(realnum) );
		memset(&rfield.flux_isotropic_save[lo]  , 0 , n*sizeof(realnum) );
		memset(&rfield.SummedOcc[lo]            , 0 , n*sizeof(realnum) );
		memset(&rfield.SummedDif[lo]            , 0 , n*sizeof(realnum) );
		memset(&rfield.flux_accum[lo]           , 0 , n*sizeof(realnum) );
		memset(&rfield.otslin[lo]               , 0 , n*sizeof(realnum) );
		memset(&rfield.otscon[lo]               , 0 , n*sizeof(realnum) );
		memset(&rfield.ConInterOut[lo]          , 0 , n*sizeof(realnum) );
		memset(&rfield.outlin[0][lo]               , 0 , n*sizeof(realnum) );
		memset(&rfield.outlin_noplot[lo]        , 0 , n*sizeof(realnum) );
		memset(&rfield.ConOTS_local_OTS_rate[lo], 0 , n*sizeof(realnum) );
		memset(&rfield.ConOTS_local_photons[lo] , 0 , n*sizeof(realnum) );
		memset(&opac.OldOpacSave[lo]            , 0 , n*sizeof(double) );
		memset(&opac.opacity_abs[lo]            , 0 , n*sizeof(double) );
		memset(&opac.opacity_sct[lo]            , 0 , n*sizeof(double) );
		memset(&opac.albedo[lo]                 , 0 , n*sizeof(double) );
		memset(&opac.FreeFreeOpacity[lo]        , 0 , n*sizeof(double) );

		/* these are not defined on first iteration */
		memset( &opac.E2TauAbsTotal[lo]        , 0 , n*sizeof(realnum) );
		memset( &opac.E2TauAbsOut[lo]          , 0 , n*sizeof(realnum) );
		memset( &opac.TauAbsTotal[lo]          , 0 , n*sizeof(realnum) );

		for( i=lo; i <= ihi; i++ )
		{
			opac.TauTotalGeo[0][i] = opac.taumin;
			opac.TauAbsGeo[0][i] = opac.taumin;
			opac.TauScatGeo[0][i] = opac.taumin;
			opac.tmn[i] = 1.;
			opac.ExpZone[i] = 1.;
			opac.E2TauAbsFace[i] = 1.;
			opac.ExpmTau[i] = 1.;
			opac.OpacStatic[i] = 1.;
		}
		/* also zero out fine opacity fine grid fine mesh array */
		memset(rfield.fine_opac_zone , 0 , (unsigned long)rfield.nfine_malloc*sizeof(realnum) );
		/* also zero out fine opacity array */
		memset(rfield.fine_opt_depth , 0 , (unsigned long)rfield.nfine_malloc*sizeof(realnum) );
	}
	return;
}
