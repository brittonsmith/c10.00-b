/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CO_step fills in matrix for heavy elements molecular routines */
#include "cddefines.h"
#include "taulines.h"
#include "dense.h"
#include "ionbal.h"
#include "thermal.h"
#include "phycon.h"
#include "hmi.h"
#include "dynamics.h"
#include "conv.h"
#include "trace.h"
#include "coolheavy.h"
#include "timesc.h"
#include "thirdparty.h"
#include "mole.h"
#include "mole_co_priv.h"
/* Nick Abel between July and October of 2003 assisted Dr. Ferland in improving the heavy element 
 * molecular network in Cloudy. Before this routine would predict negative abundances if 
 * the fraction of carbon in the form of molecules came close to 100%. A reorganizing of 
 * the reaction network detected several bugs.  Treatment of "coupled reactions",
 * in which both densities in the reaction rate were being predicted by Cloudy, were also 
 * added.  Due to these improvements, Cloudy can now perform calculations
 * where 100% of the carbon is in the form of CO without predicting negative abundances
 *
 * Additional changes were made in November of 2003 so that our reaction 
 * network would include all reactions from the TH85 paper.  This involved 
 * adding silicon to the chemical network.  Also the reaction rates were
 * labeled to make identification with the reaction easier and the matrix 
 * elements of atomic C, O, and Si are now done in a loop, which makes 
 * the addition of future chemical species (like N or S) easy.
 * */

void CO_solve(
	/* set true if we found neg pops */
	bool *lgNegPop, 
	/* set true if we tried to compute the pops, but some were zero */
	bool *lgZerPop )
{

	int32 merror;
	long int i, j, k, n,
		nelem , ion , ion2;
	double
		co_denominator;
	realnum cartot_mol, oxytot_mol;
	struct chem_element_s *element;
	double sum;

	DEBUG_ENTRY( "CO_solve()" );

	CO_step();                  /* Calculate the matrix elements */

	/* Ugly hack to deal with species which have become uncoupled */
	for(i=0;i<mole.num_comole_calc;i++)
	{
		sum = 0.;
		for(j=0;j<mole.num_comole_calc;j++)
		{
			sum = sum+fabs(mole.c[i][j]);
		}
		if(sum < SMALLFLOAT && fabs(mole.b[i]) < SMALLFLOAT)
		{
			mole.c[i][i] = 1.;
		}
	}

	cartot_mol = dense.xMolecules[ipCARBON] + findspecies("C")->hevmol + findspecies("C+")->hevmol;
	oxytot_mol = dense.xMolecules[ipOXYGEN] + findspecies("O")->hevmol + findspecies("O+")->hevmol;

	ASSERT( cartot_mol >= 0. && oxytot_mol >= 0.);

	*lgNegPop = false; 
	*lgZerPop = false;

	/* zero out molecular charge transfer rates */
	for(nelem=ipLITHIUM; nelem < LIMELM; ++nelem)
	{
		for(ion=0; ion < nelem+2; ++ion)
		{
			/*zero out the arrays */				
			mole.sink[nelem][ion]     = 0.;
			mole.source[nelem][ion]   = 0.;

			for(ion2=0; ion2< nelem+2; ++ion2)
			{
				mole.xMoleChTrRate[nelem][ion][ion2] = 0.;
			}
		}
	}


	/* >>06 jun 29, add advective terms here */
	/* Add rate terms for dynamics to equilibrium, makes c[][] non-singular 
	 * do before cross talk with heavy elements is done to not double count charge
	 * transfer */
	if( iteration > dynamics.n_initial_relax+1 && dynamics.lgAdvection && ! dynamics.lgEquilibrium && dynamics.Rate != 0.0 ) 
	{
		for(i=0;i<mole.num_comole_calc;i++) 
		{
			if(COmole[i]->n_nuclei > 1) {
				mole.c[i][i] -= dynamics.Rate;
				/* this is the net gain of the species due to advection -
				 * in this form we have the net gain as the difference between
				 * current species flowing out of this region, downstream, and
				 * new advected material coming into this region from upstream */
				mole.b[i] -= (COmole[i]->hevmol*dynamics.Rate-dynamics.CO_molec[i]);
			}
		}
	}

	/* >>chng Oct. 21, 2004 -- Terms that contribute to the ionization balance of C, O, S, Si, and N
	will now be inserted directly into the ion solver.  Before, we just corrected the recombination rate
	by scaling saying that IONIZATION_RATE = RECOMBINATION_RATE * [n(X+)/n(X0)].  This code follows the logic
	of hmole_step, written by Robin Williams. */


	/* sink terms should include only terms that either form or destroy an atomic or singly ionized
		element in the network, but not both.  This means that terms that destroy X at the expense of 
		forming X+ can't be included.  The rate of destruction of X or X+ is equal to the formation of 
		the inverse.  We therefore add this rate to sink, which gets rid of this dependence */

	/* The following code is all generalized.  After all the molecules, the total number being mole.num_heavy_molec,
	   there are 2*mole.num_elements remaining.  The first set, equal to mole.num_elements, are the ionized elemental
	   species.  The second set, also equal to mole.num_elements, are the atomic species.  They are in order such that

	   array number for X = array number for (X+) + mole.num_elements

	   In the future, if one wants to add another element to the network, such as Na the macros must be changed.
	   Also, the array number for Na and Na+ must obey the equation above.  If this is done, then the sources,
	   sinks, and molecular recombination terms are automatically added!!! */

	/* COmole[i]->nelem_hevmol[j] is just the dominant element for a given molecule.  For an element this
	   string is just the element itself!! */

	/* >> chng 2006 Sep 02 rjrw: Change to using element properties rather than properties of elements as molecules, for ease of
	   generalization -- ordering of atoms and ions no longer hardwired */

	for(i=0;i<mole.num_elements; ++i)
	{ 
		element = chem_element[i];
		mole.sink[element->ipCl][0] -= (mole.c[element->ipMl][element->ipMl] + mole.c[element->ipMl][element->ipMlP])*dense.lgElmtOn[element->ipCl];
		mole.sink[element->ipCl][1] -= (mole.c[element->ipMlP][element->ipMlP] + mole.c[element->ipMlP][element->ipMl] )*dense.lgElmtOn[element->ipCl];
	}

	/* source terms must be multiplied by the density of the corresponding matrix element */

	for(j=0;j < mole.num_elements; ++j)
	{
		element = chem_element[j];
		for(i=0; i < mole.num_comole_calc; ++i)
		{
			mole.source[element->ipCl][1] += COmole[i]->hevmol*mole.c[i][element->ipMlP]*dense.lgElmtOn[element->ipCl];
			mole.source[element->ipCl][0] += COmole[i]->hevmol*mole.c[i][element->ipMl]*dense.lgElmtOn[element->ipCl];
		}
	}

	/* subtract out diagonal terms, they are already in sink.  Also subtract recombination terms,
	   as these are done by mole.xMoleChTrRate */

	for(i=0;i < mole.num_elements; ++i)
	{ 
		element = chem_element[i];
		mole.source[element->ipCl][1]   -= ( mole.c[element->ipMlP][element->ipMlP]*COmole[element->ipMlP]->hevmol + 
			mole.c[element->ipMl][element->ipMlP]*COmole[element->ipMl]->hevmol)*dense.lgElmtOn[element->ipCl];

		mole.source[element->ipCl][0]   -= ( mole.c[element->ipMl][element->ipMl]*COmole[element->ipMl]->hevmol + 
			mole.c[element->ipMlP][element->ipMl]*COmole[element->ipMlP]->hevmol)*dense.lgElmtOn[element->ipCl];

		/* The source terms, as they are right now, include negative contributions.  This is because the 
		linearization "trick" creates source terms.   Take for example the reaction:
		C+        H2O      =>        HCO+   H          
		This reaction destroys C+, but in the act of linearizing, there will be the following term:
		mole.c[ipH2O][ipCP]
		This is only a numerical "trick" and not a real matrix term.  All terms like this
		have to be removed to get the "true" source terms.  Fortunately, the sum of all these terms
		is just mole.b[i].  To remove these terms, we have to add mole.b[i] if the overall contribution from these terms
		was negative, subtract if their contribution was positive.  This is done by subtracting mole.b[i] 
		from the current value of source. */

		mole.source[element->ipCl][1] = mole.source[element->ipCl][1] - mole.b[element->ipMlP];
		mole.source[element->ipCl][0] = mole.source[element->ipCl][0] - mole.b[element->ipMl];

	}

	/* negative source terms are actually destruction mechanisms, so should go into the sink vector.
	   source and sinks have different units, so if source is negative divide by the density of
	   the corresponding element to get same units.  For example, if:

	   mole.source[ipCARBON][0]

	   is negative, we divide by dense.xIonDense[ipCARBON][0] to get rate in same units as mole.sink */

	for(i=2; i < LIMELM; ++i)
	{
		for(j=0; j < 2; ++j)
		{
			/*if source term is negative, make it a sink and then set to zero*/ 
			if(mole.source[i][j] < 0)
			{
				mole.sink[i][j] -= (mole.source[i][j]/SDIV(dense.xIonDense[i][j]));
				mole.source[i][j] = 0;
			}
		}
	}
	/* These are rates of recombination for atomic and singly ionized species that are due to molecular processes
	   This will be added to ion_solver to get correct ionization balance */

	for(i=0;i < mole.num_elements; ++i)
	{ 
		element = chem_element[i];
		mole.xMoleChTrRate[element->ipCl][1][0] = (realnum)(mole.c[element->ipMlP][element->ipMl] - 
				ionbal.RateRecomTot[element->ipCl][0])*dense.lgElmtOn[element->ipCl];

		mole.xMoleChTrRate[element->ipCl][0][1] =  (realnum)(mole.c[element->ipMl][element->ipMlP] - 
				ionbal.RateIonizTot(element->ipCl,0))*dense.lgElmtOn[element->ipCl];
	}

	/* If rate for going from 1-0 or 0-1 is negative then it should be added to the inverse process */
	for(i=2; i < LIMELM; ++i)
	{
		for(j=0; j < 2; ++j)
		{
			for(k=0; k< 2; ++k)
			{
			/*only possible charge transfers are 0-1 or 1-0 */		
				if(j != k)
				{
					if( mole.xMoleChTrRate[i][j][k] < 0. )
					{
						mole.xMoleChTrRate[i][k][j] -= mole.xMoleChTrRate[i][j][k];
						mole.xMoleChTrRate[i][j][k] = 0.;
					}
				}
			}
		}
	}

	for(n=0;n<mole.num_elements;n++) 
	{
		tot_ion[n] = dense.gas_phase[chem_element[n]->ipCl];
		for( i=2; i < chem_element[n]->ipCl+2; i++ ) 
		{
			tot_ion[n] -= dense.xIonDense[chem_element[n]->ipCl][i];
		}
	}


	/* at this point the cartot_mol, sum of molecules and atom/first ion,
	 * should be equal to the gas_phase minus double and higher ions */
	/*fprintf(ioQQQ," dbuggggas\t %.2f\t%f\t%f\n",fnzone,
		cartot_mol/cartot_ion,
		oxytot_mol/oxytot_ion);*/

	/* these will be used in population equation in case of homogeneous equation */


	/* rjrw 2006 Aug 08: messing with the matrix like this will likely break advection --
		 was done correctly in hmole */

	if( iteration > dynamics.n_initial_relax+1 && dynamics.lgAdvection) 
		fixit();

	for(n=0;n<mole.num_elements;n++) 
	{
		mole.b[chem_element[n]->ipMl] = tot_ion[n];  
	}

	/* <<chng 03 Nov 11--Nick Abel,  Set up the non-zero matrix elements that go into the conservation equations
	   for atomic C, O, and Si, this is now set up by looping over the atomic species in
	   co.h and setting the number of atoms of C, O, and Si equal to the variable findspecies("CARBON")->nElem, 
	   findspecies("OXYGEN")->nElem, and findspecies("SILICON")->nElem respectively.  For every element (excluding Hydrogen) not in
	   the network an if statement will be needed */

	/* loop over last mole.num_elements elements in co vector - these are atoms of the mole.num_elements */
	for( j=0; j < mole.num_comole_calc; j++ )	
	{
		for(i=0;i<mole.num_elements;i++) 
		{
			element = chem_element[i];
			mole.c[j][element->ipMl] = COmole[j]->nElem[element->ipCl];
		}
	}

	/*--------------------------------------------------------------------
		* */
	/*printf( " Here are all matrix elements\n" );

	for(i=0; i < mole.num_comole_calc; i++)

	{


			for(j=0; j < mole.num_comole_calc; j++)
			{
				printf( "%4.4s", COmole[i]->label );
				printf( "%4.4s\t", COmole[j]->label );
				printf( "%.8e,\n", mole.c[j][i] );
			}
	}
	printf( "b's are:\n" );
	for(i=0; i < mole.num_comole_calc; i++)
	{
		printf( "%.8e,\n", mole.b[i] );
	}*/


	if( trace.lgTrace && trace.lgTr_CO_Mole )
	{
		fprintf( ioQQQ, " COMOLE matrix\n           " );

#		if 0
		for( i=0; i < MIN2(mole.num_comole_calc,8); i++ )
		{
			fprintf( ioQQQ, "%4.4s", COmole[i]->label );
		}
		fprintf( ioQQQ, "     \n" );

		for( j=0; j < mole.num_comole_calc; j++ )
		{
			fprintf( ioQQQ, " %4.4s", COmole[j]->label );
			fprintf( ioQQQ, " " );
			for( i=0; i < MIN2(mole.num_comole_calc,8); i++ )
			{
				fprintf( ioQQQ, "%9.1e", mole.c[i][j] );
			}
			fprintf( ioQQQ, "\n" );
		}
#		endif

		fprintf( ioQQQ, " COMOLE matrix\n           " );

		for( i=0; i < mole.num_comole_calc; i++ )
		{
			fprintf( ioQQQ, "%4.4s", COmole[i]->label );
		}
		fprintf( ioQQQ, "     \n" );

		for( j=0; j < mole.num_comole_calc; j++ )
		{
			fprintf( ioQQQ, " %4.4s", COmole[j]->label );
			fprintf( ioQQQ, " " );
			for( i=0; i < mole.num_comole_calc; i++ )
			{
				/*fprintf( ioQQQ, "%9.1e", mole.c[i][j] );*/
				fprintf( ioQQQ, ", %.15e", mole.c[i][j] );
			}
			fprintf( ioQQQ, "\n" );
		}

		fprintf( ioQQQ, " COMOLE b\n           " );
		for( j=0; j < mole.num_comole_calc; j++ )
		{
			fprintf( ioQQQ, " %4.4s", COmole[j]->label );
			fprintf( ioQQQ, ", %.15e\n",mole.b[j] );
		}

#if 0
		fprintf( ioQQQ, 
			" COMOLE H2 den:%10.3e, H2,3+=%10.2e%10.2e Carb sum=%10.3e Oxy sum=%10.3e\n", 
			hmi.H2_total, 
			hmi.Hmolec[ipMH2p], 
			hmi.Hmolec[ipMH3p], 
			mole.c[mole.num_heavy_molec][findspecies("C")->index], 
						 mole.c[mole.num_heavy_molec][findspecies("O")->index] );

#endif
	}

	/* remember longest molecule destruction timescales */
	for( j=0; j < mole.num_comole_calc; j++ )
	{
		if( -mole.c[j][j] > SMALLFLOAT )
		{
			/* this is rate CO is destroyed, equal to formation rate in equilibrium */
			timesc.AgeCOMoleDest[j] = -1./mole.c[j][j];
			/* moved to radinc */
			/*timesc.BigCOMoleForm = (realnum)MAX2(timesc.BigCOMoleForm,timesc.AgeCOMoleForm);*/
		}
		else
		{
			timesc.AgeCOMoleDest[j] = 0.;
		}
	}

	/* copy to amat, saving c for later use */
	for( j=0; j < mole.num_comole_calc; j++ )
	{
		for( i=0; i < mole.num_comole_calc; i++ )
		{
			mole.amat[i][j] = mole.c[i][j];
		}
	}

	merror = 0;

	getrf_wrapper(mole.num_comole_calc,mole.num_comole_calc,mole.amat[0],mole.num_comole_calc, ipiv,&merror);

	if( merror != 0 )
	{
		fprintf( ioQQQ, " CO_solve getrf_wrapper error\n" );
		cdEXIT(EXIT_FAILURE);
	}

	getrs_wrapper('N',mole.num_comole_calc,1,mole.amat[0],mole.num_comole_calc,ipiv,mole.b,mole.num_comole_calc,&merror);

	if( merror != 0 )
	{
		fprintf( ioQQQ, " CO_solve: dgetrs finds singular or ill-conditioned matrix\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* check for negative populations, which happens when 100% co */
	*lgNegPop = false;
	*lgZerPop = false;
	co_denominator = 0;
#	if 0
	if(fnzone > 1.02)
	{
			for( i=0; i < mole.num_comole_calc; i++)
			{			

				for( j=61; j < 62; j++  )
				{

					printf("ZONE\t%.5e\tSPECIES_1\t%s\tSPECIES_2\t%s\tDEST_RATE\t%.3e\tbvec\t%.3e\n",
						fnzone, COmole[i]->label,COmole[j]->label,mole.c[i][j],mole.b[i]);

				}
			}
	}
#	endif
	for( i=0; i < mole.num_comole_calc; i++ )
	{			
		/* fprintf(stderr,"%ld [%d] %g\n",i,mole.num_comole_calc,mole.b[i]); */
		if( mole.b[i] < 0. )
		{

			/* >>chng 03 sep 11, the following*/
			/* comparison between the solution vector from 32 bit machines vs
			 * maple on a pc shows that very small negative numbers are produced by
			 * the linear algebra package used by the code, but is positive with
			 * the math package.  allow very small negative numbers to occur,
			 * and silently reset to a postive number. 
			 *
			 * >>chng 04 feb 25
			 * Here we want to check if the negative abundance is a significant
			 * fraction of any of the species in the network.  The code below
			 * checks to see which elements the molecule in question is made of,
			 * then finds which one has the least abundance.  The one with
			 * the least abundance will then be used as the divisor in checking 
			 * if the negative abundance is too large to be ignored.*/

			/*  >>chng 04 apr 21, when an element in the chemical network is 
			 *  turned off, the molecular abundance of a species containing that
			 *  element was not zero but rather a small number of order 1e-20. When 
			 *  these abundances went negative, the code thought the abundances were
			 *  important because it checks the ratio of the species to its gas phase
			 *  abundance.  When the gas phase abundance is zero, the denominator is set 
			 *  to SMALLFLOAT, which is ~1e-35.  This led to a ratio of molecule to gas
			 *	phase of ~1e15, clearly unphysical.  Here the value of co_denominator 
			 *  will be set to a high value if the gas phase abundance is zero.*/

			if( dense.lgElmtOn[COmole[i]->nelem_hevmol] )
			{
				co_denominator = dense.gas_phase[COmole[i]->nelem_hevmol];
			}
			else
			{
				/* >>chng 04 apr 20, set to zero if element is off */
				co_denominator = 1e10;
			}


			/* we must have a positive value by now, unlesss element is turned off */
			ASSERT( co_denominator > SMALLFLOAT || !dense.lgElmtOn[COmole[i]->nelem_hevmol] );

			/* >>chng 04 feb 28, threshold from -1e-10 to -1e-6, the 
			 * level of roundoff in a realnum */
			/*if( mole.b[i] / MAX2 (co_denominator, SMALLFLOAT) > -1e-10) */
			/*if( mole.b[i] / SDIV(co_denominator) > -1e-6 ) */
			/* >>chng 04 oct 31, change limit from -1e-6 to -1e-5, produced only
			 * a few comments in func_map.in */
			if( mole.b[i] / SDIV(co_denominator) > -1e-5 ) 
			{
				/* this case, it is only slightly negative relative to
				 * the parent species - multiply by -1 to make positive
				 * and press on */
				mole.b[i] *= -1.;
			}
			else
			{
				/*>>chng 04 mar 06 press on */
				*lgNegPop = true;
				/* 05 jul 16, there had been a bug such that CO was always evaluated as long as
				 * the temperature was below 20,000K.  Once fixed, one sim turned on CO mid way
				 * into the H+ zone and had neg pops during first trys - change check on whether it
				 * is to be commented on only if this is not first zone with CO soln */
				/*if( nzone>0 )*/
				if( nzone>co.co_nzone )
				{
					static long int nzFail=-2;
					long int nzFail2 = (long)fnzone*100;
					/* during a map we may be in search phase but not in first zone */
					if( !conv.lgSearch )
						fprintf(ioQQQ," PROBLEM ");
					fprintf(ioQQQ,
						" CO_solve neg pop for species %li %s, value is %.2e rel value is %.2e zone %.2f Te %.4e Search?%c\n",
						i , 
						COmole[i]->label , 
						/* abs value */
						mole.b[i], 
						/* relative value */
						mole.b[i] / SDIV(co_denominator) ,
						fnzone,
						phycon.te,TorF( conv.lgSearch ) );
					/* if well beyond search phase and still having problems, announce them 
					 * only announce one failure per sweep through solvers by ConvFail */
					if( nzFail2 !=nzFail )
					{
						nzFail = nzFail2;
						ConvFail( "pops" , "CO");
					}
				}
				mole.b[i] *= -1;
				/*>>chng 04 jan 24 set to zero instead of *-1,
				 * h2_cr.in co density went to infinite, with some negative
				 * var getting to - inf 
				 *>>chng 04 nov 03, remove this - not needed h2_cr works without
				mole.b[i] = SMALLFLOAT;*/
			} 
		}
		else if( mole.b[i] == 0. )
		{
			/* this is not used for anything in calling routine and 
			 * could be cleaned up */
			*lgZerPop = true;
			/* >>chng 04 feb 28, zero pop not really a problem in very deep neutral
			 * gas - this happens */
#			if 0
			if( /*nzone>0 ||*/ CODEBUG>1 )
			{
				fprintf(ioQQQ," PROBLEM ");
				fprintf(ioQQQ,
					" CO_solve zero pop for species %li %s, value is %.2e zone %li\n",
					i , COmole[i]->label , mole.b[i], nzone);
			}
#			endif
			mole.b[i] = SMALLFLOAT;
		}
		COmole[i]->hevmol = (realnum)mole.b[i];
		/* check for NaN */
		ASSERT( !isnan( COmole[i]->hevmol ) );
	}
	/* >>chng 04 mar 06 pass negative pop as a problem, but not a fatal one,
	 * also do not call this during 0th zone when search for conditions
	 * is underway */
	/* 05 jul 16, there had been a bug such that CO was always evaluated as long as
	 * the temperature was below 20,000K.  Once fixed, on sim turned on CO mid way
	 * into the H+ zone and had neg pops during first trys - change check on whether it
	 * is to be commented on only if this is not first zone with CO soln */
	/*if( *lgNegPop && (nzone>0 &&!conv.lgSearch)  )*/
	if( *lgNegPop && (nzone>co.co_nzone &&!conv.lgSearch)  )
	{
		conv.lgConvPops = false;
		fprintf(ioQQQ," CO network negative population occurred, Te=%.4e, calling ConvFail. ",
			phycon.te);
		fprintf( ioQQQ, " CO/C=%9.1e", findspecies("CO")->hevmol/SDIV(dense.gas_phase[ipCARBON]) );
		fprintf( ioQQQ, "\n" );
		/*ConvFail("pops", "CO");*/
		*lgNegPop = false;
	}
	/* if negative pops were present need to renorm b to get
	 * proper sum rule */
	if( 0 && *lgNegPop )
	{

		for(j=0; j<2; ++j )
		{
			double sumcar = 0.;
			double sumoxy = 0.;
			double renorm;
			for( i=0; i < mole.num_comole_calc; i++ )
			{
				/* this case, different molecules */
				sumcar += COmole[i]->hevmol*COmole[i]->nElem[ipCARBON];
				sumoxy += COmole[i]->hevmol*COmole[i]->nElem[ipOXYGEN];
			}
			renorm = (tot_ion[0] + tot_ion[1]) / SDIV( sumcar+sumoxy);
			if(j)
				fprintf(ioQQQ,"\t%f\n", renorm);
			else
				fprintf(ioQQQ,"renormco\t%.2f\t%f", fnzone, renorm);
			for( i=0; i < mole.num_comole_calc; i++ )
			{
				COmole[i]->hevmol *= (realnum)renorm;
				/* check for NaN */
				ASSERT( !isnan( COmole[i]->hevmol ) );
			}
		}
	}

	if( merror != 0 )
	{
		fprintf( ioQQQ, " COMOLE matrix inversion error, MERROR=%5ld zone=%5ld\n", 
				 (long)merror, nzone );
		ShowMe();
		fprintf( ioQQQ, " Product matrix\n           " );

		for( i=0; i < MIN2(mole.num_comole_calc,8); i++ )
		{
			fprintf( ioQQQ, "%4.4s", COmole[i]->label );
		}
		fprintf( ioQQQ, "     \n" );

		for( j=0; j < mole.num_comole_calc; j++ )
		{
			fprintf( ioQQQ, " %4.4s", COmole[j]->label );
			fprintf( ioQQQ, " " );

			for( i=0; i < MIN2(mole.num_comole_calc,8); i++ )
			{
				fprintf( ioQQQ, "%9.1e", mole.amat[i][j]*
						 COmole[i]->hevmol );
			}
			fprintf( ioQQQ, "\n" );
		}

		if( mole.num_comole_calc >= 9 )
		{
			fprintf( ioQQQ, " COMOLE matrix\n           " );
			for( i=8; i < mole.num_comole_calc; i++ )
			{
				fprintf( ioQQQ, "%4.4s", COmole[i]->label );
			}
			fprintf( ioQQQ, "     \n" );

			for( j=0; j < mole.num_comole_calc; j++ )
			{
				fprintf( ioQQQ, " %4.4s", COmole[j]->label );
				fprintf( ioQQQ, " " );
				for( i=8; i < mole.num_comole_calc; i++ )
				{
					fprintf( ioQQQ, "%9.1e", 
							 mole.amat[i][j]* COmole[i]->hevmol );
				}
				fprintf( ioQQQ, "\n" );
			}
		}

		fprintf( ioQQQ, " Mole dens:" );
		for( j=0; j < mole.num_comole_calc; j++ )
		{
			fprintf( ioQQQ, " %4.4s:%.2e", COmole[j]->label
					 , COmole[j]->hevmol );
		}
		fprintf( ioQQQ, " \n" );

		cdEXIT(EXIT_FAILURE);
	}

	if( trace.lgTrace && trace.lgTr_CO_Mole )
	{
		fprintf( ioQQQ, " Product matrix\n           " );

		for( i=0; i < MIN2(mole.num_comole_calc,8); i++ )
		{
			fprintf( ioQQQ, "%4.4s", COmole[i]->label );
		}
		fprintf( ioQQQ, "     \n" );

		for( j=0; j < mole.num_comole_calc; j++ )
		{
			fprintf( ioQQQ, " %4.4s", COmole[j]->label );
			fprintf( ioQQQ, " " );
			for( i=0; i < MIN2(mole.num_comole_calc,8); i++ )
			{
				fprintf( ioQQQ, "%9.1e", mole.amat[i][j]*
						 COmole[i]->hevmol );
			}
			fprintf( ioQQQ, "\n" );

		}

		if( mole.num_comole_calc >= 9 )
		{
			fprintf( ioQQQ, " COMOLE matrix\n           " );
			for( i=8; i < mole.num_comole_calc; i++ )
			{
				fprintf( ioQQQ, "%4.4s", COmole[i]->label );
			}
			fprintf( ioQQQ, "     \n" );

			for( j=0; j < mole.num_comole_calc; j++ )
			{
				fprintf( ioQQQ, " %4.4s", COmole[j]->label );
				fprintf( ioQQQ, " " );

				for( i=8; i < mole.num_comole_calc; i++ )
				{
					fprintf( ioQQQ, "%9.1e", mole.amat[i][j]* COmole[i]->hevmol );
				}
				fprintf( ioQQQ, "\n" );
			}
		}

		fprintf( ioQQQ, " Mole dens:" );
		for( j=0; j < mole.num_comole_calc; j++ )
		{
			fprintf( ioQQQ, " %4.4s:%.2e", COmole[j]->label
					 , COmole[j]->hevmol );
		}
		fprintf( ioQQQ, " \n" );
	}

	/* heating due to CO photodissociation */
	co.CODissHeat = (realnum)(CO_findrate("PHOTON,CO=>C,O")*1e-12);

	thermal.heating[0][9] = co.CODissHeat;

	/* now set total density of each element locked in gas phase */
	for( i=0;i<LIMELM; ++i )
	{
		dense.xMolecules[i] = 0.;
	}

	/* total number of H per unit vol in molecules,
	 * of course not including H0/H+ */
	for(i=0;i<N_H_MOLEC;i++) 
	{
		dense.xMolecules[ipHYDROGEN] += hmi.Hmolec[i]*hmi.nProton[i];
	}
	dense.xMolecules[ipHYDROGEN] -= (hmi.Hmolec[ipMH] + hmi.Hmolec[ipMHp]);

	/* >>chng 02 sep 05, dense.xMolecules[[ipHYDROGEN] is the part of H 
	 * that is not computed in hmole
	 * add in gas phase abundances locked in molecules
	 * H is special since treated in ion_solver with all hmole molecules
	 * xMolecules are the densities of these species that are done in co,
	 * this sum is only up to mole.num_heavy_molec and so does not include the atoms/ions */
	dense.H_sum_in_CO = 0;
	for(i=0; i<mole.num_comole_calc; ++i)
	{
		if(COmole[i]->n_nuclei <= 1)
			continue;
		/* special sum of H in CO network */
		dense.H_sum_in_CO += COmole[i]->nElem[ipHYDROGEN]*COmole[i]->hevmol;

		for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
		{
			if( mole.lgElem_in_chemistry[nelem] )
			{
				if( COmole[i]->hevmol > BIGFLOAT )
				{
					fprintf(ioQQQ, "PROBLEM DISASTER mole_co_solve has found "
						"a CO-network molecule with an insane abundance.\n");
					fprintf(ioQQQ, "Species number %li with abundance %.2e\n",
						i , COmole[i]->hevmol);
					TotalInsanity();
				}
				else if( conv.lgSearch && 
					COmole[i]->nElem[nelem]*COmole[i]->hevmol > dense.gas_phase[nelem])
				{
					if( !conv.lgSearch )
						fprintf(ioQQQ,
							"PROBLEM COmole limit trip call %li Seatch? %c nMole %li "
							"n(mole) %.2e n(gas) %.2e)\n", 
							conv.nTotalIoniz,
							TorF(conv.lgSearch),i , 
							COmole[i]->nElem[nelem]*COmole[i]->hevmol,
							dense.gas_phase[nelem]);
					COmole[i]->hevmol = dense.gas_phase[nelem]/COmole[i]->nElem[nelem]/2.f;
				}
				dense.xMolecules[nelem] +=COmole[i]->nElem[nelem]*COmole[i]->hevmol;
			}
		}
	}

	/* This takes the average of an element's atomic and singly ionized density from ion_solver 
	   and comole and sets the solution to the updated COmole[i]->hevmol value.  The macro mole.num_heavy_molec
	   is the number of heavy element molecules in the network.  In addition there are 2*mole.num_elements
	   in the network, half atomic elements and half ionized elements.  Finally, the total number of species
	   in the network, including elements, equals mole.num_comole_calc.  The following will first take the average of 
	   the ionized species (between mole.num_heavy_molec and mole.num_heavy_molec+mole.num_elements) and the other statement will
	   take the average of the atomic species (between mole.num_heavy_molec+mole.num_elements and mole.num_comole_calc).  
	   This is generalized, so that if one wants to insert, for example, a sodium chemistry, once the macro's
	   are updated this if statement immediately works */

	/* >> chng 2006 Sep 02 rjrw: use vectors rather than offsets above to index, for ease of generalization */

	for(i=0;i<mole.num_elements; ++i)
	{ 
		element = chem_element[i];
		/*printf("MOL\t%3e\tION\t%e\tSPECIES\t%s\n", COmole[i]->hevmol,dense.xIonDense[COmole[i]->nelem_hevmol][1],COmole[i].label);*/
		COmole[element->ipMlP]->hevmol = ((dense.xIonDense[element->ipCl][1]+ COmole[element->ipMlP]->hevmol)/2)*dense.lgElmtOn[element->ipCl];
		/*printf("MOL\t%3e\tION\t%e\tSPECIES\t%s\n", COmole[i]->hevmol,dense.xIonDense[COmole[i]->nelem_hevmol][0],COmole[i].label);*/
		COmole[element->ipMl]->hevmol = ((dense.xIonDense[element->ipCl][0]+ COmole[element->ipMl]->hevmol)/2)*dense.lgElmtOn[element->ipCl];

		/* check for NaN */
		ASSERT( !isnan( COmole[element->ipMlP]->hevmol ) );
		ASSERT( !isnan( COmole[element->ipMl]->hevmol ) );
	}
	
	findspecies("^13CO")->hevmol = findspecies("CO")->hevmol / co.C12_C13_isotope_ratio;

	/* check whether ion and chem solvers agree yet */
	if( dense.lgElmtOn[ipCARBON] &&
		fabs(dense.xIonDense[ipCARBON][0]- findspecies("C")->hevmol)/SDIV(dense.gas_phase[ipCARBON]) >0.001 )
	{
		conv.lgConvIoniz = false;
		sprintf( conv.chConvIoniz, "CO C0 con");
		conv.BadConvIoniz[0] = dense.xIonDense[ipCARBON][0];
		conv.BadConvIoniz[1] = findspecies("C")->hevmol;
	}
	else if( dense.lgElmtOn[ipCARBON] &&
		fabs(dense.xIonDense[ipCARBON][1]- findspecies("C+")->hevmol)/SDIV(dense.gas_phase[ipCARBON]) >0.001 )
	{
		conv.lgConvIoniz = false;
		sprintf( conv.chConvIoniz, "CO C1 con");
		conv.BadConvIoniz[0] = dense.xIonDense[ipCARBON][1];
		conv.BadConvIoniz[1] = findspecies("C+")->hevmol;
	}
	else if( dense.lgElmtOn[ipOXYGEN] &&
		fabs(dense.xIonDense[ipOXYGEN][0]- findspecies("O")->hevmol)/SDIV(dense.gas_phase[ipOXYGEN]) >0.001 )
	{
		conv.lgConvIoniz = false;
		sprintf( conv.chConvIoniz, "CO O0 con");
		conv.BadConvIoniz[0] = dense.xIonDense[ipOXYGEN][0];
		conv.BadConvIoniz[1] = findspecies("O")->hevmol;
	}
	else if( dense.lgElmtOn[ipOXYGEN] &&
		fabs(dense.xIonDense[ipOXYGEN][1]- findspecies("O+")->hevmol)/SDIV(dense.gas_phase[ipOXYGEN]) >0.001 )
	{
		conv.lgConvIoniz = false;
		sprintf( conv.chConvIoniz, "CO O1 con");
		conv.BadConvIoniz[0] = dense.xIonDense[ipOXYGEN][1];
		conv.BadConvIoniz[1] = findspecies("O+")->hevmol;
	}

	/* now update ionization distribution of the elements we just did */
	for(i=0;i<mole.num_elements; ++i)
	{ 
		element = chem_element[i];
		dense.xIonDense[element->ipCl][0] = COmole[element->ipMl]->hevmol*dense.lgElmtOn[element->ipCl]; 
		dense.xIonDense[element->ipCl][1] = COmole[element->ipMlP]->hevmol*dense.lgElmtOn[element->ipCl]; 
		dense.IonLow[element->ipCl] = 0;
		dense.IonHigh[element->ipCl] = MAX2( dense.IonHigh[element->ipCl] , 1 );
	}

	/* if populations not conserved then not converged */
#	define EPS_MOLE	0.1
	if( dense.lgElmtOn[ipHYDROGEN] &&
		dense.xMolecules[ipHYDROGEN] > dense.gas_phase[ipHYDROGEN]*(1.+EPS_MOLE) )
	{
		conv.lgConvIoniz = false;
		sprintf( conv.chConvIoniz, "COcon%2i",ipHYDROGEN );
		conv.BadConvIoniz[0] = dense.xMolecules[ipHYDROGEN];
		conv.BadConvIoniz[1] = dense.gas_phase[ipHYDROGEN];
	}
	else if( dense.lgElmtOn[ipCARBON] &&
		dense.xMolecules[ipCARBON] > dense.gas_phase[ipCARBON]*(1.+EPS_MOLE) )
	{
		conv.lgConvIoniz = false;
		sprintf( conv.chConvIoniz, "COcon%2i",ipCARBON );
		conv.BadConvIoniz[0] = dense.xMolecules[ipCARBON];
		conv.BadConvIoniz[1] = dense.gas_phase[ipCARBON];
	}
	else if( dense.lgElmtOn[ipOXYGEN] &&
		dense.xMolecules[ipOXYGEN] > dense.gas_phase[ipOXYGEN]*(1.+EPS_MOLE) )
	{
		conv.lgConvIoniz = false;
		sprintf( conv.chConvIoniz, "COcon%2i",ipOXYGEN );
		conv.BadConvIoniz[0] = dense.xMolecules[ipOXYGEN];
		conv.BadConvIoniz[1] = dense.gas_phase[ipOXYGEN];
	}
	else if( dense.lgElmtOn[ipSILICON] &&
		dense.xMolecules[ipSILICON] > dense.gas_phase[ipSILICON]*(1.+EPS_MOLE) )
	{
		conv.lgConvIoniz = false;
		sprintf( conv.chConvIoniz, "COcon%2i",ipSILICON );
		conv.BadConvIoniz[0] = dense.xMolecules[ipSILICON];
		conv.BadConvIoniz[1] = dense.gas_phase[ipSILICON];
	}
	else if( dense.lgElmtOn[ipSULPHUR] &&
		dense.xMolecules[ipSULPHUR] > dense.gas_phase[ipSULPHUR]*(1.+EPS_MOLE) )
	{
		conv.lgConvIoniz = false;
		sprintf( conv.chConvIoniz, "COcon%2i",ipSULPHUR );
		conv.BadConvIoniz[0] = dense.xMolecules[ipSULPHUR];
		conv.BadConvIoniz[1] = dense.gas_phase[ipSULPHUR];
	}
	else if( dense.lgElmtOn[ipNITROGEN] &&
		dense.xMolecules[ipNITROGEN] > dense.gas_phase[ipNITROGEN]*(1.+EPS_MOLE) )
	{
		conv.lgConvIoniz = false;
		sprintf( conv.chConvIoniz, "COcon%2i",ipNITROGEN );
		conv.BadConvIoniz[0] = dense.xMolecules[ipNITROGEN];
		conv.BadConvIoniz[1] = dense.gas_phase[ipNITROGEN];
	}

	else if( dense.lgElmtOn[ipCHLORINE] &&
		dense.xMolecules[ipCHLORINE] > dense.gas_phase[ipCHLORINE]*(1.+EPS_MOLE) )
	{
		conv.lgConvIoniz = false;
		sprintf( conv.chConvIoniz, "COcon%2i",ipCHLORINE );
		conv.BadConvIoniz[0] = dense.xMolecules[ipCHLORINE];
		conv.BadConvIoniz[1] = dense.gas_phase[ipCHLORINE];
	}

	/* the total N photo dissociation rate, cm-3 s-1 */
	co.nitro_dissoc_rate = (realnum)(CO_dissoc_rate("N"));
	return;

/*lint +e550 */
/*lint +e778 constant expression evaluates to 0 in operation '-' */
}

#	undef EPS_MOLE
