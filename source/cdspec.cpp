/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*cdSPEC returns the spectrum needed for Keith Arnaud's XSPEC */
/*Spec_cont called by cdSPEC to generate actual spectrum */
#include "cddefines.h"
#include "cddrive.h"
#include "physconst.h"
#include "geometry.h"
#include "radius.h"
#include "rfield.h"
#include "opacity.h"
#include "grid.h"

/* 
 * this routine returns the spectrum needed for Keith Arnaud's XSPEC
 * X-Ray analysis code.  It should be called after cdDrive has successfully
 * computed a model.  the calling routine must ensure that the vectors
 * have enough space to store the resulting spectrum, 
 * given the bounds and energy resolution 
 */

/* array index within energy array in Cloudy grids -  this has file scope */
static long int iplo , iphi;

/* energies of lower and upper bound of XSPEC cell we are considering */
static double Elo , Ehi;

void cdSPEC( 
	/* option - the type of spectrum to be returned
	 * 1	the incident continuum 4\pi nuJ_nu, , erg cm-2 s-1
	 *
	 * 2	the attenuated incident continuum, same units
	 * 3	the reflected continuum, same units
	 *
	 * 4	diffuse continuous emission outward direction
	 * 5	diffuse continuous emission, reflected
	 *
	 * 6	collisional+recombination lines, outward
	 * 7	collisional+recombination lines, reflected
	 * 
	 * 8	pumped lines, outward <= not implemented
	 * 9	pumped lines, reflected <= not implemented
	 *
	 *		all lines and continuum emitted by the cloud assume full coverage of 
	 *		continuum source */
	int nOption ,

	/* the energy of the lower edge of each cell 
	 * (in Ryd to be consistent with all the rest of Cloudy) */
	double EnergyLow[] , 

	/* the number of cells + 1*/
	long int nEnergy ,

	/* the returned spectrum, same size is two energy spectra (see option), returns nEnergy -1 pts */
	double ReturnedSpectrum[] )

{
	/* this pointer will bet set to one of the cloudy continuum arrays */
	realnum *cont , 
		refac;
	long int ncell , j;

	/* flag saying whether we will need to free cont at the end */
	bool lgFREE;

	DEBUG_ENTRY( "cdSPEC()" );

	ASSERT( nEnergy <= rfield.nflux );

	if( nOption == 1 )
	{
		/* this is the incident continuum, col 2 of save continuum command */
		cont = rfield.flux_total_incident[0];
		lgFREE = false;
	}
	else if( nOption == 2 )
	{
		/* the attenuated transmitted continuum, no diffuse emission,
		 * col 3 of save continuum command */
		cont = rfield.flux[0];
		lgFREE = false;
	}
	else if( nOption == 3 )
	{
		/* reflected incident continuum, col 6 of save continuum command */
		lgFREE = false;
		cont = rfield.ConRefIncid[0];
	}
	else if( nOption == 4 )
	{
		/* diffuse continuous emission outward direction  */
		cont = (realnum*)MALLOC( sizeof(realnum)*(size_t)rfield.nupper );

		/* need to free the vector once done */
		lgFREE = true;
		refac = (realnum)radius.r1r0sq*geometry.covgeo;
		for( j=0; j<rfield.nflux; ++j)
		{
			cont[j] = rfield.ConEmitOut[0][j]*refac;
		}
	}
	else if( nOption == 5 )
	{
		/* reflected diffuse continuous emission */
		cont = (realnum*)MALLOC( sizeof(realnum)*(size_t)rfield.nupper );
		/* need to free the vector once done */
		lgFREE = true;
		refac = (realnum)radius.r1r0sq*geometry.covgeo;
		for( j=0; j<rfield.nflux; ++j)
		{
			cont[j] = rfield.ConEmitReflec[0][j]*refac;
		}
	}
	else if( nOption == 6 )
	{
		/* all outward lines,   */
		cont = (realnum*)MALLOC( sizeof(realnum)*(size_t)rfield.nupper );
		/* need to free the vector once done */
		lgFREE = true;
		/* correct back to inner radius */
		refac = (realnum)radius.r1r0sq*geometry.covgeo;
		for( j=0; j<rfield.nflux; ++j)
		{
			/* units of lines here are to cancel out with tricks applied to continuum cells
			 * when finally extracted below */
			cont[j] = rfield.outlin[0][j] *rfield.widflx[j]/rfield.anu[j]*refac;
		}
	}
	else if( nOption == 7 )
	{
		/* all reflected lines */
		if( geometry.lgSphere )
		{
			refac = 0.;
		}
		else
		{
			refac = 1.;
		}

		cont = (realnum*)MALLOC( sizeof(realnum)*(size_t)rfield.nupper );
		/* need to free the vector once done */
		lgFREE = true;
		for( j=0; j<rfield.nflux; ++j)
		{
			/* units of lines here are to cancel out with tricks applied to continuum cells
			 * when finally extracted below */
			cont[j] = rfield.reflin[0][j] *rfield.widflx[j]/rfield.anu[j]*refac;
		}
	}
	else
	{
		fprintf(ioQQQ," cdSPEC called with impossible nOption (%i)\n", nOption);
		cdEXIT(EXIT_FAILURE);
	}

	/* this is array index used in Spec_cont */
	iplo = 0;
	iphi = 0;
	/* now generate the continua */
	for( ncell = 0; ncell < nEnergy-1; ++ncell )
	{
		/* Spec_cont will find the spectrum between these two energies */
		Elo = EnergyLow[ncell];
		Ehi = EnergyLow[ncell+1];
		ReturnedSpectrum[ncell] = cont[ncell] * EN1RYD * rfield.anu2[ncell] / rfield.widflx[ncell];
	}

	/* need to free the vector once done if this flag was set */
	if( lgFREE )
	{
		free(cont);
	}
	return;
}


/* returns in units photons/cm^2/s/bin */
void cdSPEC2( 
	/* option - the type of spectrum to be returned
	 * 1	the incident continuum 4\pi nuJ_nu, , erg cm-2 s-1
	 *
	 * 2	the attenuated incident continuum, same units
	 * 3	the reflected continuum, same units
	 *
	 * 4	diffuse emission, lines + continuum, outward
	 * 5	diffuse emission, lines + continuum, reflected
	 *
	 * 6	diffuse continuous emission outward direction
	 * 7	diffuse continuous emission, reflected
	 *
	 * 8	total transmitted, lines and continuum
	 * 9	total reflected, lines and continuum
	 *
	 *10    exp(-tau) to the illuminated face
	 *
	 *		all lines and continuum emitted by the cloud assume full coverage of 
	 *		continuum source */
	int nOption ,

	/* the number of cells */
	long int nEnergy,

	/* the index of the lowest and highest energy bounds to use. */
	long ipLoEnergy,
	long ipHiEnergy,

	/* the returned spectrum, same size is two energy spectra (see option), returns nEnergy -1 pts */
	realnum ReturnedSpectrum[] )

{
	realnum refac;

	DEBUG_ENTRY( "cdSPEC2()" );

	ASSERT( ipLoEnergy >= 0 );
	ASSERT( ipLoEnergy < ipHiEnergy );
	ASSERT( ipHiEnergy < rfield.nupper );
	ASSERT( nEnergy == (ipHiEnergy-ipLoEnergy+1) );
	ASSERT( nEnergy >= 2 );

	ASSERT( nOption <= NUM_OUTPUT_TYPES );

	for( long i = 0; i < nEnergy; i++ )
	{
		long j = ipLoEnergy + i;

		if( j >= rfield.nflux )
		{
			ReturnedSpectrum[i] = SMALLFLOAT;
			continue;
		}
		
		if( nOption == 0 )
		{
			/* the attenuated incident continuum */
			realnum flxatt = rfield.flux[0][j]*
				(realnum)radius.r1r0sq * rfield.trans_coef_total[j];

			/* the outward emitted continuum */
			realnum conem = (rfield.ConEmitOut[0][j] + rfield.outlin[0][j])*
				(realnum)radius.r1r0sq*geometry.covgeo;

			/* the reflected continuum */
			realnum flxref = rfield.ConRefIncid[0][j] + rfield.ConEmitReflec[0][j] +
				rfield.reflin[0][j];

			ReturnedSpectrum[i] = flxatt + conem + flxref;
		}
		else if( nOption == 1 )
		{
			/* this is the incident continuum, col 2 of save continuum command */
			ReturnedSpectrum[i] = rfield.flux_total_incident[0][j];
		}
		else if( nOption == 2 )
		{
			/* the attenuated transmitted continuum, no diffuse emission,
			 * col 3 of save continuum command */
			ReturnedSpectrum[i] = rfield.flux[0][j]*
				(realnum)radius.r1r0sq * rfield.trans_coef_total[j];
		}
		else if( nOption == 3 )
		{
			/* reflected incident continuum, col 6 of save continuum command */
			ReturnedSpectrum[i] = rfield.ConRefIncid[0][j];
		}
		else if( nOption == 4 )
		{
			/* all outward diffuse emission */
			/* correct back to inner radius */
			refac = (realnum)radius.r1r0sq*geometry.covgeo;
			ReturnedSpectrum[i] = (rfield.outlin[0][j]+rfield.ConEmitOut[0][j])*refac;
		}
		else if( nOption == 5 )
		{
			/* all reflected diffuse emission */
			if( geometry.lgSphere )
			{
				refac = 0.;
			}
			else
			{
				refac = 1.;
			}

			ReturnedSpectrum[i] = (rfield.reflin[0][j]+rfield.ConEmitReflec[0][j])*refac;
		}
		else if( nOption == 6 )
		{
			/* all outward line emission */
			/* correct back to inner radius */
			refac = (realnum)radius.r1r0sq*geometry.covgeo;
			ReturnedSpectrum[i] = rfield.outlin[0][j]*refac;
		}
		else if( nOption == 7 )
		{
			/* all reflected line emission */
			if( geometry.lgSphere )
			{
				refac = 0.;
			}
			else
			{
				refac = 1.;
			}

			ReturnedSpectrum[i] = rfield.reflin[0][j]*refac;
		}
		else if( nOption == 8 )
		{
			/* total transmitted continuum */
			/* correct back to inner radius */
			refac = (realnum)radius.r1r0sq*geometry.covgeo;
			ReturnedSpectrum[i] = (rfield.ConEmitOut[0][j]+ rfield.outlin[0][j])*refac
				+ rfield.flux[0][j]*(realnum)radius.r1r0sq*rfield.trans_coef_total[j];
		}
		else if( nOption == 9 )
		{
			/* total reflected continuum */
			ReturnedSpectrum[i] = rfield.ConRefIncid[0][j] + rfield.ConEmitReflec[0][j] +
				rfield.reflin[0][j];
		}
		else if( nOption == 10 )
		{
			/* this is exp(-tau) */
			/* This is the TOTAL attenuation in both the continuum and lines.  
			 * Jon Miller discovered that the line attenuation was missing in 07.02 */
			ReturnedSpectrum[i] = opac.ExpmTau[j]*rfield.trans_coef_total[j];
		}
		else
		{
			fprintf(ioQQQ," cdSPEC called with impossible nOption (%i)\n", nOption);
			cdEXIT(EXIT_FAILURE);
		}

		ASSERT( ReturnedSpectrum[i] >=0.f );
	}

	return;
}
