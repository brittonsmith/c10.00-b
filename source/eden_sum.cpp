/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "hmi.h"
#include "trace.h"
#include "grainvar.h"
#include "rfield.h"
#include "mole.h"
#include "dense.h"
#include "taulines.h"

// eden_sum: sum all contributions to the electron density, sets variable dense.EdenTrue
// called by ConvBase - ConvEdenIoniz actually updates the electron density dense.eden
// here we allow dense.EdenTrue to become negative, ConvEdenIoniz will deal with this
// and will also assure that dense.eden is always positive. this routine always returns 0
int eden_sum(void)
{
	DEBUG_ENTRY( "eden_sum()" );

	/* this variable is set with the set eden command, 
	 * is supposed to override physical electron density */
	if( dense.EdenSet > 0.f )
	{
		dense.EdenTrue = dense.EdenSet;
		dense.eden_from_metals = 1.;

		if( trace.lgTrace || trace.lgESOURCE )
			fprintf( ioQQQ, "     eden_sum zn: %.2f eden set to: %.4e\n", fnzone, dense.EdenSet );
	}
	else
	{
		/* EdenExtra is normally zero, set with EDEN command, to add extra e- */
		dense.EdenTrue = dense.EdenExtra;

		/* sum over all ions */
		double eden_ions[LIMELM];
		double sum_all_ions = 0.;
		double sum_metals = 0.;
		for( long nelem=ipHYDROGEN; nelem < LIMELM; nelem++ )
		{
			eden_ions[nelem] = 0.;
			for( long ion=1; ion <= nelem+1; ion++ )
				eden_ions[nelem] += ion*dense.xIonDense[nelem][ion];

			sum_all_ions += eden_ions[nelem];
			if( nelem >= ipLITHIUM )
				sum_metals += eden_ions[nelem];
		}
		dense.EdenTrue += sum_all_ions;

		/* add electrons from all H molecules, had just been H- */
		double hmole_eden = 0.;
		for( long i=0; i < N_H_MOLEC; i++ )
		{
			/* hmi.nElectron is zero for H+ since counted as an ion, -1 for H-, etc */
			hmole_eden += hmi.Hmolec[i]*hmi.nElectron[i];
		}
		dense.EdenTrue += hmole_eden;

		/* electrons contributed by heavy molecules */
		co.comole_eden = 0.;
		for( long i=0; i < mole.num_comole_calc; i++ )
		{
			if( COmole[i]->n_nuclei != 1 )
				co.comole_eden += COmole[i]->hevmol*COmole[i]->nElec;
		}
		dense.EdenTrue += co.comole_eden;

		/* gv.lgGrainElectrons - should grain electron source/sink be included in overall electron sum?
		 * default is true, set false with no grain electrons command */
		dense.EdenTrue += gv.TotalEden*gv.lgGrainElectrons;

		/* fraction of electrons from ions heavier than helium */
		dense.eden_from_metals = safe_div( sum_metals, dense.EdenTrue, 1. );

		if( trace.lgTrace || trace.lgESOURCE )
		{
			fprintf( ioQQQ, 
				 "     eden_sum zn: %.2f current: %.4e new true: %.4e ions: %.4e comole: %.4e"
				 " hmole: %.4e grain: %.4e extra: %.4e LaOTS: %.4e\n",
				 fnzone ,
				 dense.eden , 
				 dense.EdenTrue , 
				 sum_all_ions ,
				 co.comole_eden ,
				 hmole_eden ,
				 gv.TotalEden*gv.lgGrainElectrons,
				 dense.EdenExtra ,
				 rfield.otslin[Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].ipCont-1] );

			if( trace.lgNeBug )
			{
				for( long nelem=ipHYDROGEN; nelem < LIMELM; nelem++ )
				{
					if( nelem == 0 )
						fprintf( ioQQQ, "      eden_sum H -Ne:" );
					else if( nelem == 10 )
						fprintf( ioQQQ, "      eden_sum Na-Ca:" );
					else if( nelem == 20 )
						fprintf( ioQQQ, "      eden_sum Sc-Zn:" );
					fprintf( ioQQQ, " %.4e", eden_ions[nelem] );
					if( nelem%10 == 9 )
						fprintf( ioQQQ, "\n" );
				}
			}
		}
	}

	/* case where electron density is set with set eden command, make sure we use it */
	ASSERT( dense.EdenSet <= 0.f || fp_equal((realnum)dense.EdenTrue, dense.EdenSet) );

	return 0;
}
