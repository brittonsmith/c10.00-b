/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CO_step fills in matrix for heavy elements molecular routines */
#include "cddefines.h"
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
/* Robin Williams in August 2006 onwards reorganized the coding to cut down repetitions.  
 * This isolated several further bugs, and allows a sigificant number of lines of
 * code to be eliminated.  The balance of S2/S2+ amd ClO/ClO+ seems highly sensitive
 * (with small log scale results varying significantly if the order of arithmetic
 * operations is changed) -- I suspect this may imply a bug somewhere.
 * */
/*lint -e778 constant expression evaluatess to 0 in operation '-' */
/*=================================================================*/


void CO_step(void)
{
	long int i, j, n, nt, ratei, ratej;
	struct COmole_rate_s *rate;
	double rate_tot, rate_deriv[MAXREACTANTS], rated, rk, rate_bval;

	DEBUG_ENTRY("CO_step()");
	/* zero out array used for formation rates */
	for( i=0; i < mole.num_comole_calc; i++ )
	{
		mole.b[i] = 0.;
	}
	for( j=0; j < mole.num_comole_tot; j++ )
	{
		for( i=0; i < mole.num_comole_calc; i++ )
		{
			mole.c[j][i] = 0.;
		}
	}


	/* Call all the routines that set up the matrix 
	 * CO_solve will call this routine, therefore all other matrix elements are
	 * included here so that, when CO_solve is called, everything is accounted for */	

	/* All now cross-validated with new treatment, switching causes only v. minor
	 * differences in results */
	/* Revised molecular network implementation */
	/* Fills matrix elements for passive species -- these can be used to
		 derive sources & sinks resulting from this part of the network */
	nt = coreactions.n;
	for(n=0;n<nt;n++) 
	{
		rate = coreactions.list[n];
		rk = rate->rk;
		for(i=0;i<rate->nrates;i++)
		{
			rate_deriv[i] = rk;
			for(j=0;j<rate->nrates;j++)
			{
				if(i!=j)
				{
					rate_deriv[i] *= rate->rate_species[j]->hevmol;
				}
			}
		}

		if(rate->nreactants != 1) 
		{
			rate_tot = rate_deriv[0]*rate->rate_species[0]->hevmol;
			rate_bval = (rate->nreactants-1)*rate_tot;
			for(i=0;i<rate->nreactants;i++)
			{
				ratei = rate->reactants[i]->index;
				mole.b[ratei] -= rate_bval;
			}

			for(i=0;i<rate->nproducts;i++)
			{
				ratei = rate->products[i]->index;
				mole.b[ratei] += rate_bval;
			}
		}

		for(j=0;j<rate->nrates;j++)
		{
			ratej = rate->rate_species[j]->index;
			rated = rate_deriv[j];
			for(i=0;i<rate->nreactants;i++)
			{
				mole.c[ratej][rate->reactants[i]->index] -= rated;
			}
			for(i=0;i<rate->nproducts;i++)
			{
				mole.c[ratej][rate->products[i]->index] += rated;
			}
		}
	}

}
