/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*hmole_step do a step in chemical network */
/*hmirat compute radiative association rate for H- */
/* >> chng 02 nov 7 rjrw, Mole Moreliano:
 *   changes to linearized iterative form */
/* from Robin Williams:
The process for these kind of problems seems to be pretty uniform:
switch printsol on to check which terms in the chemical matrix change,
and next to switch on the prints in the matrix assembly which apply to
the species involved to find what reactions are involved.  It's a bit
of a pain grepping down to find the 47th reaction, so I guess some
kind of naming scheme for the reactions may come in handy (I'd thought
about generating the in and out vectors from a text string, e.g. "H +
H => H2";-), but you'd have to verify uniqueness).
*/
/*lint -e778 const express evaluates to 0 */
/*lint -e725 expect positive indentation */
#include "cddefines.h"
#include "physconst.h"
#include "iso.h"
#include "atmdat.h"
#include "grainvar.h"
#include "ionbal.h"
#include "dense.h"
#include "secondaries.h"
#include "mole.h"
#include "mole_co_priv.h"
#include "opacity.h"
#include "rfield.h"
#include "thermal.h"
#include "timesc.h"
#include "trace.h"
#include "phycon.h"
#include "doppvel.h"
#include "thirdparty.h"
#include "gammas.h"
#include "h2.h"
#include "dynamics.h"
#include "conv.h"
#include "radius.h"
#include "hextra.h" 
#include "hmi.h"
#include "taulines.h"

/* HP cc cannot compile following except in -O1 mode */
#if defined(__HP_aCC)
#	pragma OPT_LEVEL 1
#endif

#define ABSLIM  1e-12

/* Calculate number of elements in an integer vector */
#define INTSZ(a)     (sizeof(a)/sizeof(int))

struct Hmole_rate_s {
	int index;
	int nreactants, nrates, nproducts;
	int reactants[MAXREACTANTS];
	int rate_species[MAXREACTANTS];
	int products[MAXPRODUCTS];
	double rk; 
	struct Hmole_rate_s *next;
};
typedef struct Hmole_rate_s reaction;

/* Generate new element for reaction list */
reaction *newreaction(
	/* reaction index, a number to reference current reaction */
	int rindex, 
	/* vector of ints that are incoming species, the reactants */
	int *in, 
	/* number of reactants */
	int nin, 
	/* same for products */
	int *out, 
	/* number of products */
	int nout, 
	/* rate determining species if they are not the same as the reactants,
	 * these are non-null in cases where a part of the chemical network
	 * acts as a catalyst in the reaction 
	 * as H2g + H2g -> H2g + H2* or where */
	int *rate, 
	/* number of these */
	int nrate)
{
	static reaction *list = NULL, *r;
	static int poolsize=1, index = 0;
	int i;

	/* this routine is called only to initialize structure with information
	 * on the reactants and products - this does not change, so no need
	 * to do this but once.  In later cases only the rate coefficients are
	 * updated in hmole_step 
	 * in hmole_stop on later call linked list is incremented upon each
	 * reaction, and rate coefficient for current temperature is stored 
	 * in r->rk 
	 * r is pointer to structure of reaction information including
	 * r->next, which is the next reaction in the linked list */
	/* fprintf(ioQQQ,"New reaction %d %d %d\n",rindex,nin,nout); */

	/* default assumption for chemical kinetics,
	 * rate is determined by all incoming species that are in *in 
	 * when not null then the rate determining species are not the
	 * same as *in */
	if(rate == NULL) 
	{
		rate=in;
		nrate=nin;
	}

	/* space for the linked list "list" */
	/*lint -e701 shift left of signed quantity */
	if(list == NULL || index == poolsize) 
	{
		poolsize <<=1;
		list = ((reaction *)MALLOC( (size_t)poolsize*sizeof(reaction) ));
		index = 0;
	}
	/*lint +e701 */

	/* fprintf(ioQQQ,"Getting element %d+1 of %d\n",index,poolsize); */
	r = list+index;
	index++;
	r->next = NULL;
	r->index = rindex;
	ASSERT(nin <= MAXREACTANTS && nout <= MAXPRODUCTS && nrate <= MAXREACTANTS);

	r->nreactants =	nin;
	r->nrates =	nrate;
	r->nproducts  =	nout;

	/* incoming reactants */
	for(i=0; i<r->nreactants; i++)
		r->reactants[i] = in[i];

	/* rate determining species */
	for(i=0; i<r->nrates; i++)
		r->rate_species[i] = rate[i];

	/* outgoing products */
	for(i=0; i<r->nproducts; i++)
		r->products[i] = out[i];

	return r;
}

// icc 10.0 miscompiles this routine with higher optimization
#if defined (__ICC) && defined(__amd64)
#pragma optimization_level 1
#endif

/*hmole_step do a step in chemical network */
void hmole_step(int *nFixup, double *error)
{
	enum {PRINTSOL = false};

	int32 ipiv[N_H_MOLEC];

	long int i, 
	  ipConserve, 
	  j, 
	  limit ,
	  mol;

	int printsol = PRINTSOL;

	bool lgNegPop;
	int iworst;
	/* >>chng 05 jul 31, from float to double, since very nearly 1 for H - route */
	double frac_H2star_grains,
		frac_H2star_hminus; 
	double frac_H2_grains,
		frac_H2_hminus; 
	double sum_first_ions;

	double 
		bhneut, 
		Hneut,
		c3bod,
		cionhm, 
		corr, 
		H2star_deexcit,
		deexc_htwo,
		deexc_hneut,
		desh2p, 
		etmp,
		eh3p_3h,
		Boltz_fac_H2_H2star,
		fhneut, 
		gamheh, 
		h1fnd, 
		h1rat, 
		h2pcin,  
		h2phhp, 
		h2pion, 
		H2star_excit ,
		radath, 
		hmihph2p,
		h2phmhhh,
		h2crphh,
		h2scrphh,
		h2crphphm,
		h2scrphphm,
		h2crphpeh,
		h2scrphpeh,
		h2crh2pe,
		h2crhphe,
		h2scrhphe,
		h2scrh2pe,
		h2pehh,
		h3ph2ph,
		hphmhhpe,
		h2hmhh2e,
		hehmeheh,
		hephmhhe,
		fracneg,
		fracnegtmp,
		fracnegfac,
		sum_H0_Hp,
		conserve,
		rate,
		rk,
		rated,
		rate_deriv[MAXREACTANTS],
		sinkrate[MAXREACTANTS],
		T_ortho_para_crit ,
		TStancil;


	static double 
	  gamtwo, 
	  h2phet, 
	  proton_sum_old,
	  proton_sum_new;

	static double /*amat[N_H_MOLEC][N_H_MOLEC], */
	  b2pcin, 
	  *amat=NULL,
	  *bvec=NULL/*[N_H_MOLEC]*/, 
	  *Hmolec_old=NULL/*[N_H_MOLEC]*/, 
	  **c=NULL/*[N_H_MOLEC+1][N_H_MOLEC+1]*/, 
	  plte;
	static double oatomic = -1., oion = -1.;
	/*static long int iter_eval = -2;*/

	/* if this is still true then must create space for arrays */
	static bool lgMustMalloc = true;

	static reaction *rlist = NULL;
	reaction *r;
	long int rindex, ratei, ratej;
	
	realnum AveVelH = GetAveVelocity( dense.AtomicWeight[ipHYDROGEN] );
	realnum AveVelH2 = GetAveVelocity( 2.f*dense.AtomicWeight[ipHYDROGEN] );

	DEBUG_ENTRY( "hmole_step()" );

	if( lgMustMalloc )
	{
		/* on very first call must create space */
		lgMustMalloc = false;

		bvec = ((double*)MALLOC( (size_t)N_H_MOLEC*sizeof(double) ));
		Hmolec_old = ((double*)MALLOC( (size_t)N_H_MOLEC*sizeof(double) ));
		amat = ((double*)MALLOC( (size_t)(N_H_MOLEC*N_H_MOLEC)*sizeof(double) ));
		c = ((double**)MALLOC( (size_t)N_H_MOLEC*sizeof(double *) ));
		for( i=0; i<N_H_MOLEC; ++i )
		{
			/* this is the Jacobian array, MALLOC sets it to NaN, first reagents
			 * will be filled in, then Jacobian array set in loop below.  Search
			 * for Jacobian array */
			c[i] = ((double*)MALLOC( (size_t)N_H_MOLEC*sizeof(double) ));
		}
	}

	/* Assume no error for cases with abundances set */
	*error = 0;

	/* there are two "no molecules" options, the no co, which turns off everything,
	 * and the no n2, which only turns off the h2.  in order to not kill the co
	 * part we still need to compute the hydrogen network here, and then set h2 to
	 * small values */
	if( hmi.lgNoH2Mole )
	{
		dense.xMolecules[ipHYDROGEN] = 0.;

		/* these are the molecular species */
		for(mol=0; mol<N_H_MOLEC; ++mol) 
		{
			hmi.Hmolec[mol] = 0.;
		}
		hmi.Hmolec[ipMH] = dense.xIonDense[ipHYDROGEN][0];
		hmi.Hmolec[ipMHp] = dense.xIonDense[ipHYDROGEN][1];
		hmi.H2_total = 0.;
		/* this is where the transition struc expects to find the H2 abundance */
		//dense.xIonDense[LIMELM+2][0] = 0.;
		hmi.hmihet = 0.;
		hmi.h2plus_heat = 0.;
		hmi.H2Opacity = 0.;
		hmi.hmicol = 0.;
		hmi.hmidep = 1.;
		hmi.rh2dis = 0.;
		hmi.HalphaHmin = 0.;

		hmi.HeatH2Dish_used = 0.; 
		hmi.HeatH2Dish_BigH2 = 0.;
		hmi.HeatH2Dish_TH85 = 0.;
		hmi.HeatH2Dish_BD96 = 0.;
		hmi.HeatH2Dish_BHT90 = 0.;
		hmi.HeatH2Dish_ELWERT = 0.;

		/** HeatH2Dexc_used is heating due to collisional deexcitation of vib-excited 
		* H2 actually used */
		hmi.HeatH2Dexc_used = 0.;
		hmi.HeatH2Dexc_BigH2 = 0.;
		hmi.HeatH2Dexc_TH85 = 0.;
		hmi.HeatH2Dexc_BD96 = 0.;
		hmi.HeatH2Dexc_BHT90 = 0.;
		hmi.HeatH2Dexc_ELWERT = 0.;

		/** these are derivative wrt temp for collisional processes within X */
		hmi.deriv_HeatH2Dexc_used = 0.;
		hmi.deriv_HeatH2Dexc_BigH2 = 0.;
		hmi.deriv_HeatH2Dexc_TH85 = 0.;
		hmi.deriv_HeatH2Dexc_BD96 = 0.;
		hmi.deriv_HeatH2Dexc_BHT90 = 0.;
		hmi.deriv_HeatH2Dexc_ELWERT = 0.;
		return;
	}

	/* option to force H2 abundance, for testing h2 molecules,
	 * hmi.H2_frac_abund_set is fraction in molecules that is set by 
	 * set h2 fraction command */
	if( hmi.H2_frac_abund_set > 0. )
	{
		for(mol=0;mol<N_H_MOLEC;mol++) 
		{
			hmi.Hmolec[mol] = 0.;
		}
		/* >>chng 03 jul 19, from 0 to SMALLFLOAT, to pass asserts in ConvBase,
		 * problem is that ion range has not been reset for hydrogen */
		dense.xIonDense[ipHYDROGEN][0] = dense.xIonDense[ipHYDROGEN][1] = 
			2.f*SMALLFLOAT*dense.gas_phase[ipHYDROGEN];
		/* put it all in the ground state */
		hmi.Hmolec[ipMH2g] = (realnum)(dense.gas_phase[ipHYDROGEN] * hmi.H2_frac_abund_set);
		hmi.Hmolec[ipMH2s] = 0.;

		hmi.H2_total = hmi.Hmolec[ipMH2g] + hmi.Hmolec[ipMH2s];
		/* first guess at ortho and para densities */
		h2.ortho_density = 0.75*hmi.H2_total;
		h2.para_density = 0.25*hmi.H2_total;

		hmi.hmihet = 0.;
		hmi.h2plus_heat = 0.;
		hmi.H2Opacity = 0.;
		hmi.hmicol = 0.;
		hmi.HeatH2Dish_TH85 = 0.;
		hmi.HeatH2Dexc_TH85 = 0.;
		hmi.deriv_HeatH2Dexc_TH85 = 0.;
		hmi.hmidep = 1.;
		hmi.HalphaHmin = 0.;

		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			gv.bin[nd]->rate_h2_form_grains_used = 0.;
		}
		return;
	}

	/* update these two to current values of atomic and ionized hydrogen density */
	hmi.Hmolec[ipMH] = dense.xIonDense[ipHYDROGEN][0];
	hmi.Hmolec[ipMHp] = dense.xIonDense[ipHYDROGEN][1];
	/* now copy all of H moles into old array */
	for(mol=0;mol<N_H_MOLEC;mol++) 
	{
		Hmolec_old[mol] = hmi.Hmolec[mol];
	}
	for(i=0; i<MAXREACTANTS; ++i )
	{
		rate_deriv[i] = 0.;
		sinkrate[i] = 0.;
	}

	/* collisional ionization of H-, rate from Janev, Langer et al. */
	if( phycon.te < 3074. )
	{
		cionhm = 1.46e-32*(powi(phycon.te,6))*phycon.sqrte*hmi.exphmi;
	}
	else if( phycon.te >= 3074. && phycon.te < 30000. )
	{
		cionhm = 5.9e-19*phycon.tesqrd*phycon.sqrte*phycon.te05;
	}
	else
	{
		cionhm = 1.54e-7;
	}

	/* H2 formation on grains;
	 * rate from 
	 * >>refer	H2	grain formation	Hollenbach, D., & McKee, C.F., 1979, ApJS, 41, 555 eq 3.4 3.8 */
	if( gv.lgDustOn() )
	{

#		ifndef IGNORE_QUANTUM_HEATING
		/* hmole is called before grains, so assure that all the grain stuff is properly initialized */
		GrainDrive();
#		endif

		/* these are rates (s-1) H2 will be deactivated by collisions with grains 
		 * will be incremented below 
		 * H2 ortho - para conversion on grain surface */
		hmi.rate_grain_h2_op_conserve = 0.;
		/* rate (s-1) v=0, J=1 level goes to 0 */
		hmi.rate_grain_h2_J1_to_J0 = 0.;

		/* loop over all grain species */
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
#			ifndef IGNORE_QUANTUM_HEATING
			long k, qnbin;
			double *qtemp, *qprob;
			bool lgUseQHeat = gv.lgGrainPhysicsOn && gv.bin[nd]->lgQHeat;
#			endif
			/* >>chng 02 feb 15, removed check tedust > 1.01, change in GrainsInit
			 * guarantees that all relevant parameters are initialized, PvH */

			/* sticking probability, 2H + grain equation 3.7 of
			 * >>refer	grain	phys	Hollenbach, D.J., & McKee, C.F., 1979, ApJS, 41, 555,
			 * fraction of H impacts on grain surface that stick */
			/* this sticking probability is used for both HM79 and CT02 */
			double sticking_probability_H = 1./(1. + 0.04*sqrt(gv.bin[nd]->tedust+phycon.te) + 
				0.002*phycon.te + 8e-6*phycon.te*phycon.te);

#			ifndef IGNORE_QUANTUM_HEATING
			/* >>chng 04 feb 21, included quantum heating in calculation of formation rate, PvH */
			if( lgUseQHeat )
			{
				qtemp = (double*)MALLOC((size_t)(NQGRID*sizeof(double)));
				qprob = (double*)MALLOC((size_t)(NQGRID*sizeof(double)));

				qheat(qtemp,qprob,&qnbin,nd);

				if( gv.bin[nd]->lgUseQHeat )
				{
					ASSERT( qnbin > 0 );
				}
				else
				{
					qnbin = 1;
					qprob[0] = 1.;
					qtemp[0] = gv.bin[nd]->tedust;
				}

				gv.bin[nd]->rate_h2_form_grains_HM79 = 0.;

				for( k=0; k < qnbin; k++ )
				{
					/* fraction of impacts that produce H2 before evaporation from grain surface.
					 * this is equation 3.4 of
					 * >>refer	grain	phys	Hollenbach, D.J., & McKee, C.F., 1979, ApJS, 41, 555
					 * 1e4 is ratio of total absorption sites to appropriate sites 
					 * 920 is D_H and chosen to get f_a = 0.5 at 100 K.
					 * factor of 0.6252 needed to obtain std ism rate to be 3e-17 at 100 K,
					 * the value deduced by
					 * >>refer	H2	grain physics	Jura, M., 1974, ApJ, 197, 581 */
					double conversion_efficiency_HM79 = 1/(1. + 1e4*sexp(920./qtemp[k]));
					sticking_probability_H = 1./(1. + 0.04*sqrt(qtemp[k]+phycon.te) + 
								     0.002*phycon.te + 8e-6*phycon.te*phycon.te);

					gv.bin[nd]->rate_h2_form_grains_HM79 += qprob[k] * sticking_probability_H *
						conversion_efficiency_HM79;
				}

				/* NB IntArea is total, not projected, area, must div by 4 */
				/* gv.bin[nd]->rate_h2_form_grains_HM79 has units s^-1 since gv.bin[nd]->cnv_H_pCM3 has units cm-3 */
				/* cnv_H_pCM3 converts <unit>/H (default depletion) -> <unit>/cm^3 (actual depletion), units are cm-3 */
				gv.bin[nd]->rate_h2_form_grains_HM79 *= 0.5 * AveVelH *
					gv.bin[nd]->IntArea/4. * gv.bin[nd]->cnv_H_pCM3;

				ASSERT( gv.bin[nd]->rate_h2_form_grains_HM79 > 0. );
			}
			else
#			endif
			{
				/* fraction of impacts that produce H2 before evaporation from grain surface.
				* this is equation 3.4 of
				* >>refer	grain	phys	Hollenbach, D.J., & McKee, C.F., 1979, ApJS, 41, 555
				* 1e4 is ratio of total absorption sites to appropriate sites 
				* 920 is D_H and chosen to get f_a = 0.5 at 100 K.
				* factor of 0.6252 needed to obtain std ism rate to be 3e-17 at 100 K,
				* the value deduced by
				* >>refer	H2	grain physics	Jura, M., 1974, ApJ, 197, 581 */
				double conversion_efficiency_HM79 = 1/(1. + 1e4*sexp(920./gv.bin[nd]->tedust));

				/* NB IntArea is total area per H for default abundances, not projected area, must div by 4 
				 * units s^-1 since gv.bin[nd]->cnv_H_pCM3 has units H cm-3 
				 * final units are cm s-1*/
				gv.bin[nd]->rate_h2_form_grains_HM79 = 0.5 * AveVelH * gv.bin[nd]->IntArea/4. * 
					/* cnv_H_pCM3 converts <unit>/H (default depletion) -> <unit>/cm^3 (actual depletion), units are cm-3 */
					gv.bin[nd]->cnv_H_pCM3 * sticking_probability_H * conversion_efficiency_HM79;
				ASSERT( gv.bin[nd]->rate_h2_form_grains_HM79 > 0. );
			}

#			ifndef IGNORE_QUANTUM_HEATING
			if( lgUseQHeat )
			{
				/* H2 formation on grains from 
				 * >>refer	H2	form	Cazaux, S., & Tielens, A.G.G.M., 2002, ApJ, 575, L29 */
				/* number of monolayers per second - only affects efficiency at very low or high temperatures */
				double f = 1e-10;
				/* equation 17 
				double sqrt_term = POW2( 1. + sqrt( (10000.-200.)/(600.-200.) ) );*/
				double sqrt_term = 35.399494936611667;

				gv.bin[nd]->rate_h2_form_grains_CT02 = 0.;

				for( k=0; k < qnbin; k++ )
				{
					double beta_alpha = 0.25 * sqrt_term *sexp(200./qtemp[k] );
					/* equation 16 */
					double xi =  1./ (1. + 1.3e13*sexp(1.5*1e4/qtemp[k])*sqrt_term/(2.*f) );
					/* expression for beta comes from just after equation 5 */
					double beta = 3e12 * sexp( 320. / qtemp[k] );
					/* recombination efficiency given by their equation 15, they call
					 * this epsilon_H2 */
					double recombination_efficiency_CT02 = xi / (1. + 0.005*f/2./SDIV(beta) + beta_alpha );
					sticking_probability_H = 1./(1. + 0.04*sqrt(qtemp[k]+phycon.te) + 
								     0.002*phycon.te + 8e-6*phycon.te*phycon.te);

					/* printf( " k %ld Td %.6e re*sp %.6e\n", k, qtemp[k], recombination_efficiency_CT02* */
					/* sticking_probability_H ); */

					gv.bin[nd]->rate_h2_form_grains_CT02 += qprob[k] * sticking_probability_H *
						recombination_efficiency_CT02;
				}

				/* gv.bin[nd]->IntArea integrated grain surface area Int(4pi*a^2), normalized per H, in cm^2/H,
				 * so x/4 is projected area of circle */
				/* gv.bin[nd]->cnv_H_pCM3 is H density [cm-3] times grain depletion factor */
				/* gv.bin[nd]->rate_h2_form_grains_CT02 units s-1 */
				gv.bin[nd]->rate_h2_form_grains_CT02 *= 0.5 * AveVelH *
					gv.bin[nd]->IntArea/4. * gv.bin[nd]->cnv_H_pCM3;

				ASSERT( gv.bin[nd]->rate_h2_form_grains_CT02 > 0. );

				free(qtemp);
				free(qprob);
			}
			else
#			endif
			{
				/* H2 formation on grains from 
				 * >>refer	H2	form	Cazaux, S., & Tielens, A.G.G.M., 2002, ApJ, 575, L29 */
				/* number of monolayers per second - only affects efficiency at very low or high temperatures */
				double f = 1e-10;
				/* equation 17 
				double sqrt_term = POW2( 1. + sqrt( (10000.-200.)/(600.-200.) ) );*/
				double sqrt_term = 35.399494936611667;
				double beta_alpha = 0.25 * sqrt_term *sexp(200./gv.bin[nd]->tedust );
				/* equation 16 */
				double xi =  1./ (1. + 1.3e13*sexp(1.5*1e4/ gv.bin[nd]->tedust )*sqrt_term/(2.*f) );
				/* expression for beta comes from just after equation 5 */
				double beta = 3e12 * sexp( 320. / gv.bin[nd]->tedust );
				/* recombination efficiency given by their equation 15, they call
				 * this epsilon_H2 */
				double recombination_efficiency_CT02 = xi / (1. + 0.005*f/2./SDIV(beta) + beta_alpha );

				/* gv.bin[nd]->IntArea integrated grain surface area Int(4pi*a^2), normalized per H, in cm^2/H,
				 * so x/4 is projected area of circle */
				/* gv.bin[nd]->cnv_H_pCM3 is H density [cm-3] times grain depletion factor */
				/* units s-1 */
				gv.bin[nd]->rate_h2_form_grains_CT02 = 0.5 * AveVelH * gv.bin[nd]->IntArea/4. * 
					gv.bin[nd]->cnv_H_pCM3 * sticking_probability_H * recombination_efficiency_CT02;
				ASSERT( gv.bin[nd]->rate_h2_form_grains_CT02 > 0. );
			}

#			ifndef	IGNORE_QUANTUM_HEATING
			/* reset sticking probability for code below */
			sticking_probability_H = 1./(1. + 0.04*sqrt(gv.bin[nd]->tedust+phycon.te) + 
				0.002*phycon.te + 8e-6*phycon.te*phycon.te);
#			endif

			/* rate (s-1) all H2 v,J levels go to 0 or 1, preserving nuclear spin */
			/* ortho to para on grain surfaces, taken from 
			 *>refer	H2	sticking	Le Bourlot, J., 2000, A&A, 360, 656-662 
			 * >chng 05 apr 30, GS, hmi.H2_total/dense.gas_phase[ipHYDROGEN] is removed
			 * This is used in h2.c.
			 * NB IntArea is total are per H, not projected area, must div by 4 
			 * gv.bin[nd]->cnv_H_pCM3 has units H cm-3 to product with above
			 * is cm2/H H/cm3 or cm-1 or an opacity
			 * multiply by velocity of H2, cm s-1, so product 
			 * hmi.rate_grain_h2_op_conserve has units s^-1  */
			hmi.rate_grain_h2_op_conserve += AveVelH2 * gv.bin[nd]->IntArea/4. *
				gv.bin[nd]->cnv_H_pCM3 * sticking_probability_H;

			/* ortho to para on grain surfaces, taken from 
			 *>refer	H2	sticking	Le Bourlot, J., 2000, A&A, 360, 656-662 
			 * For all grain temperatures, this process corresponds to high J going to
			 * either 0 or 1 preserving nuclear spin.  All ortho go to 1 and para go to 0.
			 * When the dust temperature is below Tcrit all 1 go to 0 and so all J go to 0.

			 * this temperature depends on grain composition, discussion left column of page 657,
			 * this is for a bare grain */
			 /** \todo	2	- put in actual composition dependent Tad - this is only valid 
			 * for bare surfaces - not ice - for ice Tad is 555K 
			 * hmi.Tad is binding energy expressed as a temperature 
			 * note that hmi.Tad is set to 800. in zero 
			 * tau_nu the first equation in section 2.5
			 * equation one paragraph before equation 2 
			 * at low grain temperatures all end in para, J=0 */

			/* AveVEl[LIMELM+2] is average speed of H2 molecules 
			 * for now assume that sticking probability for H2 on the grain is equal to
			 * that for H 
			 * efficiency factor efficiency_opr is vary fast function of t dust - 
			 * large at low Td and small at Td > T_ortho_para_crit
			 * start evaluating just above the critical temperature 
			 * T_ortho_para_crit this is roughly 24.345 K,GS */
			T_ortho_para_crit = 2. * hmi.Tad / log( POW2(60. *1.1e11)*hmi.Tad); 
			if( gv.bin[nd]->tedust < T_ortho_para_crit )
			{
				double efficiency_opr = sexp(60.*1.1e11*sqrt(hmi.Tad)*sexp(hmi.Tad/gv.bin[nd]->tedust));
				/* rate (s-1) all v,J levels go to 0, regardless of nuclear spin
				 * see above discussion for how units work out */
				hmi.rate_grain_h2_J1_to_J0 += AveVelH2 * gv.bin[nd]->IntArea/4. * 
					gv.bin[nd]->cnv_H_pCM3 * sticking_probability_H * efficiency_opr;
			}
		}
		/*fprintf(ioQQQ," H2 grain form rate HM79 %.2e  %.2e CT02 %.2e  %.2e O-P grn %.2e %.2e\n", 
			gv.bin[nd]->rate_h2_form_grains_HM79 , 
			gv.bin[nd]->rate_h2_form_grains_HM79 ,
			gv.bin[nd]->rate_h2_form_grains_CT02 , 
			gv.bin[nd]->rate_h2_form_grains_CT02 , 
			hmi.rate_grain_h2_J1_to_J0,
			hmi.rate_h2_allX_2_J1_grains
			);*/
		/* options to turn off grain collision with atom h2 collisions grains off command */
		hmi.rate_grain_h2_op_conserve *= mole.lgH2_grain_deexcitation;
		hmi.rate_grain_h2_J1_to_J0 *= mole.lgH2_grain_deexcitation;

	}
	else
	{
		/* grains are not enabled, set these to zero */
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			gv.bin[nd]->rate_h2_form_grains_CT02 = 0.;
			gv.bin[nd]->rate_h2_form_grains_HM79 = 0.;
		}
		/* rate all H2 goes to either 0 or 1 depending on ortho/para */
		hmi.rate_grain_h2_op_conserve = 0.;
		/* at low temp, rate all H2 goes to J=0 */
		hmi.rate_grain_h2_J1_to_J0 = 0.;
	}

	/* the H2 catalysis rate on grains that is actually used in calculations
	 * hmi.ScaleJura is scale factor set with set Jura scale command 
	 * units are s-1 
	 * default is 'C' Cazaux & Tielens */
	gv.rate_h2_form_grains_used_total = 0.;
	for( size_t nd=0; nd < gv.bin.size(); nd++ )
	{
		if( hmi.chJura == 'C' )
		{
			/* use the new rate by 
			 * >>refer	H2	form	Cazaux, S., & Tielens, A.G.G.M., 2002, ApJ, 575, L29 
			 * units are s-1*/
			gv.bin[nd]->rate_h2_form_grains_used = 
				gv.bin[nd]->rate_h2_form_grains_CT02*hmi.ScaleJura;
			gv.rate_h2_form_grains_used_total += gv.bin[nd]->rate_h2_form_grains_used;
		}
		else if( hmi.chJura == 'T' )
		{
			/* rate from Hollenbach & McKee 1979  */
			gv.bin[nd]->rate_h2_form_grains_used = 
				gv.bin[nd]->rate_h2_form_grains_HM79*hmi.ScaleJura;
			gv.rate_h2_form_grains_used_total += gv.bin[nd]->rate_h2_form_grains_used;
		}
		else if( hmi.chJura == 'S' )
		{
			/* H2 formation rate from 
			 * >>refer	H2	form	Sternberg, A. & Neufeld, D.A. 1999, ApJ, 516, 371 */
			gv.bin[nd]->rate_h2_form_grains_used = 
				3e-18 * phycon.sqrte / gv.bin.size() * dense.gas_phase[ipHYDROGEN]*hmi.ScaleJura;
			/* this is simple rate from Sternberg & Neufeld 99 */
			gv.rate_h2_form_grains_used_total += gv.bin[nd]->rate_h2_form_grains_used;
		}
		/*>>chng 07 jan 10, this had been C for constant, and so could never have been triggered.
		 * caught by robin Williams, fixed by nick Abel, error was in sense that any set jura rate
		 * would use Cazaux & Tielens */
		else if( hmi.chJura == 'F' )
		{
			/* command "set H2 rate" - enters log of Jura rate - C for constant,
			 * no dependence on grain properties */
			gv.bin[nd]->rate_h2_form_grains_used = hmi.rate_h2_form_grains_set*dense.gas_phase[ipHYDROGEN] / gv.bin.size();
			gv.rate_h2_form_grains_used_total += gv.bin[nd]->rate_h2_form_grains_used;
		}
	}
	ASSERT( gv.rate_h2_form_grains_used_total >= 0. );

#	ifndef IGNORE_QUANTUM_HEATING
	printf( " fnzone %.2f H2 rate %.4e\n", fnzone, gv.rate_h2_form_grains_used_total );
#	endif

	/* >>chng 03 sep 09, get ratio of excited to ground state H2 */
	if( h2.lgH2ON  && hmi.lgBigH2_evaluated && hmi.lgH2_Chemistry_BigH2 )
	{
		frac_H2star_grains = hmi.H2star_forms_grains / 
			SDIV(hmi.H2star_forms_grains+hmi.H2_forms_grains);

		frac_H2star_hminus = hmi.H2star_forms_hminus / 
			SDIV(hmi.H2star_forms_hminus+hmi.H2_forms_hminus);

		frac_H2_grains = hmi.H2_forms_grains / 
			SDIV(hmi.H2star_forms_grains+hmi.H2_forms_grains);

		frac_H2_hminus = hmi.H2_forms_hminus / 
			SDIV(hmi.H2star_forms_hminus+hmi.H2_forms_hminus);
		
		/* option print statement for above */
		/*printf( "DEBUG H2s frac grain %.3e f(H2g) %.3e ",frac_H2star_grains ,
			hmi.H2_forms_grains/SDIV(hmi.H2star_forms_grains+hmi.H2_forms_grains) );
		printf( " H2s frac h- %.3e f(H2g) %.3e\n",frac_H2star_hminus ,
			hmi.H2_forms_hminus/SDIV(hmi.H2star_forms_hminus+hmi.H2_forms_hminus));*/
	}
	else
	{
		/* the large H2 molecule was not evaluated, so we can't use exact
		 * branching ratios.  These are the distribution fractions for around 500K */
		/*These depend on temperature and distribution function and the definition of ENERGY_H2_STAR.
		  So reset the values properly*/
		/* >>chng 05 jul 13, TE, with the new definition of H2s these are both 1 */
		/* >>chng 05 jul 31, activate above print, rest for current 0.5 ev defn */
		frac_H2star_grains = 0.9416;
		frac_H2star_hminus = 0.999995062;
		frac_H2_grains = 0.0584;
		frac_H2_hminus = 4.938e-6;
	}

	/* print rate coefficient */
	/*fprintf(ioQQQ," total grain h2 form rate %.3e\n",gv.rate_h2_form_grains_used_total);*/

	/* collisional dissociation, rate from 
	 * >>refer	H2	collisional dissociation	Dove, J.E., and Mandy, M. E., 1986, ApJ, 311, L93.
	 * corr is correction for approach to high density limit
	 * H2 + H => 3H - rate very uncertain */
	corr = MIN2(6.,14.44-phycon.alogte*3.08);

	if(corr > 0.)
		corr = pow(10.,corr*Hmolec_old[ipMH]/(Hmolec_old[ipMH]+1.6e4));
	else
		corr = 1.;
	/* must kill H- when UMIST is in place since they do not consider it */
	hmi.rh2dis = (realnum)(1.55e-8/phycon.sqrte*sexp(65107./phycon.te)* corr)*co.lgUMISTrates;

	/* old hminus rate Hollenbach & McKee 1979
	 *>>chng 98 jan 02, from 2.12e4 to 2.123e4 */
	/*hmi.bh2h2p = 1.8e-12f*phycon.sqrte*phycon.te10/phycon.te01*2.f/16.f;
	hmi.rh2h2p = 1.8e-12*phycon.sqrte*phycon.te10/phycon.te01*sexp(2.123e4/
	  phycon.te);*/

	/* forward and back reactions for H2+ + H+ <=> H2+ + H */
	/*>>chng 02 oct 25, update rate from above (Hollenbach & McKee 1979) to
	 * >>refer	H2	form	Karpas, Z., Anicich, V., & Huntress, W.T. 1979, J Chem Phys, 70, 2877 
	 * following is from this paper.\:
	We note that the application of detailed balance is only strictly valid for 
	state-to-state reactions, i.e., when the , J level of the reactant and product 
	molecules are known. While typical laboratory conditions are such that  = 0 and J 
	is likely to be small for the reactant, the product , J is usually unknown. 
	For example, Krsti (2002) finds that for the reverse of reaction (1), the product 
	H2 is primarily formed into  = 4, not  = 0. Therefore, estimation of reaction (1) 
	by the application of detailed balance to the measured rate coefficient for the 
	reverse reaction gives the rate coefficient for H2( = 4), which can be as much as 
	an order of magnitude larger than for H2( = 0). Rate coefficients that are estimated 
	by detailed balance are therefore suspect. 
	*/

	/* H2+ + H => H2 + H+ */
	hmi.bh2h2p = 6.4e-10f;
#	if 0
	/* H2  +  H+  =>  H2+  +  H
	 * >>chng 04 Feb 24, get back reaction from above */
	 if(hmi.rel_pop_LTE_H2g != 0.)
		hmi.rh2h2p = hmi.bh2h2p*hmi.rel_pop_LTE_H2p/hmi.rel_pop_LTE_H2g;
	 else
		hmi.rh2h2p = 0;
#	endif
	/*>> refer	H2	chem	Savin, D.W., Krstic, P.S., Haiman, Z., & Stancil, P.C., 2004,
	 *>>refercon	ApJL, 606L, 167, astro-ph/0404288 */
	if(phycon.te<=3.e4 )
	{
		/* this is lower bound to their temperature */
		double teused = MAX2(100., phycon.te );
		double telog = log(teused);
		hmi.rh2h2p = sexp(2.123715e4/teused)*
			(((((((3.5311932e-13*telog-1.8171411e-11)*telog+3.9731542e-10)*telog-4.781372e-9)*telog+3.4172805e-8)*telog-1.4491368e-7)*telog+3.3735382e-7)*telog-3.3232183e-7);
		/* option to kill process when Leiden hacks are in place */
		hmi.rh2h2p = hmi.rh2h2p*co.lgUMISTrates;
	}
	else
		hmi.rh2h2p= 0;

	/* >>chng 05 aug 05 NPA comment The UMIST rate uses a different photorate for H2+.  We do it better.  
	Therefore, for this case we use the UMIST rate only when the UMIST hack is set.
	We cannot actually turn off all H- reactions because to do so would cause the matrix solver 
	to crash.  Therefore a couple of reactions still exist, but do not affect the Leiden models */

	t_phoHeat photoHeat;

	/* H2+  +  HNU  =>  H+  + H */
	gamtwo = GammaK(opac.ih2pnt[0],opac.ih2pnt[1],opac.ih2pof,1.,&photoHeat);

	/* this only occurs when set units rates is entered */
	if(!co.lgUMISTrates)
		gamtwo = 5.7e-10*hmi.UV_Cont_rel2_Habing_TH85_face*(realnum)sexp((1.9*rfield.extin_mag_V_point))/1.66f;

	/*GammaPrt(opac.ih2pnt[0],opac.ih2pnt[1],opac.ih2pof,ioQQQ,gamtwo,0.01);*/

	h2phet = photoHeat.HeatNet;

	/* >> chng 02 nov 15 rjrw: ionization fractions to multiply c[ipHo][*] terms
	 * as b[ipMHo] contains _both_ H0 and H+ */
	/* sum_H0_Hp = ((double)dense.xIonDense[ipHYDROGEN][0])+((double)dense.xIonDense[ipHYDROGEN][1]); */

	sum_H0_Hp = Hmolec_old[ipMH]+Hmolec_old[ipMHp];

	rindex = 0;
	r = rlist;
	/* Special case, put null reaction at head of list */
	if(r == NULL) 
	{
		int in[]={-1},out[]={-1};
		r = rlist = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	rindex++;

	/*-------------------------------------------------------------------- */

	/* H- H minus hminus balance equations
	 * (IHMI,IPMHO) == processes making H- from Ho =+sign
	 * radiative attachment: HI + NE => H-
	 * H + e -> H- + hnu */
	/* Use Newton-Raphson step to improve solution, so bvec[] contains reaction rates
	 * and c[][] components of the Jacobian of the rates */

	/* This block adds a reaction H => H- to the stack if it wasn't
	 * there already.
	 *
	 * >>>> ONLY CHANGE the elements of the in[] and out[] vectors and
	 * the rate constant, keep the rest fixed for all reactions 
	 * */
	if(r->next == NULL) {
		int in[]={ipMH},out[]={ipMHm};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;

	/* >>chng 05 aug 05, NPA comment:  The Leiden comparison does not consider H-.  Therefore the UMIST hack 
	   is used to turn off H- so that it is never important */

	r->rk = (hmi.hminus_rad_attach + hmi.HMinus_induc_rec_rate)*co.lgUMISTrates;

	/* >>chng 02 oct 29, add these two chemical processes */
	/* H- + H+ => H2+ + e
	 * equation (H6) from 
	 * >>refer	H2	chemistry	Galli,D., & Palla, F. 1998, A&A, 335, 403-420 
	 * hmihph2p = 6.9e-9f*(Tg)^(-0.35) for Tg<=8000
	 * hmihph2p = 6.9e-9f*(Tg)^(-0.9) for Tg>=8000  */
	/* >>chng 02 nov 07 rjrw, include H+ ion density in rate constant */
	if(phycon.te <= 7891.)
	{
	    /*hmihph2p = 6.9e-9*pow(phycon.te , -0.35);*/
	    hmihph2p = 6.9e-9 / (phycon.te30 * phycon.te05);
	}
	else 
	{
		/* >>chng 02 nov 18, had typo for leading coefficient here */
		/*hmihph2p = 9.6e-7*pow(phycon.te , -0.9);*/
		hmihph2p = 9.6e-7 / phycon.te90;
	}

	if(r->next == NULL) {
		int in[]={ipMHm,ipMHp},out[]={ipMH2p};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;

	/* >>chng 05 aug 05, NPA comment:The Leiden comparison does not consider H-.  
	Therefore the UMIST hack is used to turn off H- so that it is never important */

	hmihph2p = hmihph2p*co.lgUMISTrates;
	r->rk = hmihph2p;

	/* >>chng 03 feb 6 */
	/* H2+ + H- => H2 + H
	 * equation (32) from
	 * >>refer	H2+	k	Stancil, P.C, & Lepp, S, & Dalgarno, A. 1998,ApJ, 509, 1-10
	 * h2phmh2h = 1.4e-7f*pow(phycon.te/300.0, -0.5) */
	/* >>chng 03 sep 01, rm the pow function */
	/*h2phmh2h = 1.4e-7f*pow(phycon.te/300 , -0.5);*/
	/* the fits in this paper cannot be used below 10K or above 1e4K.  Limit
	 * the range of evaluation */
	TStancil = MIN2(phycon.te, 1e4 );
	TStancil = MAX2( TStancil , 10. );

	/** \todo	1	equivalent reaction for H2* is not included in chemistry, Big h2 does not include this reaction, what to do? GS */	
	/* >>chng 05 aug 05, NPA comment:The Leiden comparison does not consider H-.  
	Therefore the UMIST hack is used to turn off H- so that it is never important */

	hmi.h2phmh2h = 1.4e-7*co.lgUMISTrates*17.305/phycon.sqrte;
	if(r->next == NULL) {
		int in[]={ipMH2p,ipMHm},out[]={ipMH2g,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = hmi.h2phmh2h;

	/* >>chng 03 feb 7 */
	/* H2+ + H- => H + H + H
	 * equation (33) from
	 * >>refer	H2+	k	Stancil, P.C, & Lepp, S, & Dalgarno, A. 1998,ApJ, 509, 1-10
	 * h2phmhhh = 1.4e-7f*pow(phycon.te/300.0, -0.5) */
	/* >>chng 03 sep 01, rm the pow function */
	/*h2phmhhh = 1.4e-7f*pow(phycon.te/300 , -0.5);*/
	h2phmhhh = 1.4e-7f*17.3205/phycon.sqrte;
	if(r->next == NULL) {
		int in[]={ipMH2p,ipMHm},out[]={ipMH,ipMH,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	/* UMIST Leiden does not include H- so must kill H- reactions */
	r->rk = h2phmhhh*co.lgUMISTrates;

	/* >>chng 03 sep 30 */
	/* H+ + H- => H + H+ + e
	 * >>refer	H-	k	Paolo Lenzuni, David F. Chernoff, Edwin E. Salpeter, 1991, ApJS, 76, 759L  (table 5)
	 * hphmhhpe = 4.75e-30*pow(phycon.te,3.1); */
	hphmhhpe = 4.75e-30*pow(phycon.te,3.1);
	if(r->next == NULL) {
		int in[]={ipMHp,ipMHm},out[]={ipMH,ipMHp};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	/* UMIST Leiden does not include H- so must kill H- reactions */
	r->rk =hphmhhpe*co.lgUMISTrates;

	/** \todo	1	equivalent reaction for H2* is not included in chemistry, Big h2 does not include this reaction, what to do? GS */
	/* >>chng 03 sep 30 */
	/* H2 + H- => H + H2 + e
	 * >>refer	H-	k	Paolo Lenzuni, David F. Chernoff, Edwin E. Salpeter, 1991, ApJS, 76, 759L (table 5)
	 * h2hmhh2e = 6.74e-17*pow(phycon.te,2)*sexp(19870/phycon.te); */
	/* UMIST Leiden does not include H- so must kill H- reactions */
	h2hmhh2e = 6.74e-17*co.lgUMISTrates*phycon.tesqrd*sexp(19870/phycon.te);
	if(r->next == NULL) {
		int in[]={ipMH2g,ipMHm},out[]={ipMH,ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =h2hmhh2e;

	/* (IHMI,IHMI) = processes destroying H- =-sign
	 * photodissociation, H- + H NU => H + NE */

	if(r->next == NULL) {
		int in[]={ipMHm},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = hmi.HMinus_photo_rate;

	/* mutual neutralization with heavies, rate from Dalgarno and McCray
	 * all charged ions contribute equally,
	 * H- + A+ => H + A */
	/* >>chng 04 feb 19, do actual sum of first ions rather than following kludge */
	/* find the sum of all single ion densities for species heavier than helium */
	sum_first_ions = 0.;
	for( i=ipLITHIUM; i < LIMELM; i++ )
	{
		sum_first_ions += dense.xIonDense[i][1];
	}

	{
		/* this debug print statement compares H2 formation through grn vs H- */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC && nzone>140 )
		{
			fprintf(ioQQQ,"DEBUG sumfirstions\t%.2f\t%.4e\t%.4e\t%.4e",
				fnzone,phycon.te,
				sum_first_ions,
				sum_first_ions/dense.eden);
			for( i=ipLITHIUM; i < LIMELM; i++ )
			{
				if( dense.xIonDense[i][1]/sum_first_ions >0.1 )
					fprintf(ioQQQ,"\t%li\t%.3e",
					i,dense.xIonDense[i][1]/sum_first_ions);
			}
			fprintf(ioQQQ,"\n");
		}
	}

	hmi.hmin_ct_firstions = 4e-6f/(realnum)phycon.sqrte*atmdat.lgCTOn;

	if(r->next == NULL) {
		int in[]={ipMHm},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = hmi.hmin_ct_firstions*sum_first_ions*co.lgUMISTrates;

	/* electron collisional ionization of H- */
	cionhm *= dense.eden;

	if(r->next == NULL) {
		int in[]={ipMHm},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = cionhm*co.lgUMISTrates;

	/* inverse process; three body rec */
	c3bod = cionhm*(hmi.rel_pop_LTE_Hmin*dense.eden);

	if(r->next == NULL) {
		int in[]={ipMH},out[]={ipMHm};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = c3bod*co.lgUMISTrates;

	/* form molecular hydrogen from H minus,
	 * associative detachment:  H- + H => H2 + E */
	/* make H2 from H- 
	 * associative detachment; H- + H => H2: 
	 * >>referold	H2	rates	Browne & Dalgarno J PHys B 2, 885 */
	/* rate coefficient from 
	 * >>refer	H2	form	Launay, J.R., Le Dourneuf, M., & Zeippen, C.J., 
	 * >>refercon	1991, A&A, 252, 842-852*/
	/* >>chng 02 oct 17, temp dependent fit to rate, updated reference,
	 * about 40% larger than before */
	{
		double y , x;
		x = MAX2(10., phycon.te );
		x = MIN2(1e4, x );
		y=545969508.1323510+x*71239.23653059864;
		hmi.assoc_detach = 1./y;
	}

	/* >>chng 02 nov 7 rjrw, example case of 2-body process */
	/* this one is into ground H2 */
	if(r->next == NULL) {
		int in[]={ipMH, ipMHm},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = hmi.assoc_detach*co.lgUMISTrates*frac_H2_hminus;

	/* >>chng 03 sep 10, multiply above by correction for excited state,
	 * add below reaction for population of excited state */
	/* this one is into excited H2 */
	if(r->next == NULL) {
		int in[]={ipMH, ipMHm},out[]={ipMH2s};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = hmi.assoc_detach*frac_H2star_hminus*co.lgUMISTrates;

	{
		/* this debug print statement compares H2 formation through grn vs H- */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC && nzone>140 )
		{
			fprintf(ioQQQ," debuggggrn grn\t%.2f\t%.3e\t%.3e\tfrac\t%.3e\tH-\t%.3e\t%.3e\tfrac\t%.3e\t%.3e\t%.3e\t%.3e\n",
				fnzone ,
				gv.rate_h2_form_grains_used_total , 
				hmi.H2_forms_grains+hmi.H2star_forms_grains ,
				frac_H2star_grains,
				hmi.assoc_detach*dense.xIonDense[ipHYDROGEN][0]*hmi.Hmolec[ipMHm],
				hmi.H2star_forms_hminus+hmi.H2_forms_hminus,
				frac_H2star_hminus,
				hmi.assoc_detach,dense.xIonDense[ipHYDROGEN][0],hmi.Hmolec[ipMHm]
				);
		}
	}

	/* convert H2 into H- 
	 * the back reaction, H2(grnd) + e => H- + Ho */
	if( hmi.rel_pop_LTE_H2g > 0. )
	{ 
		hmi.assoc_detach_backwards_grnd = hmi.assoc_detach*hmi.rel_pop_LTE_Hmin/hmi.rel_pop_LTE_H2g*
			dense.eden*co.lgUMISTrates;
	}
	else
	{
		hmi.assoc_detach_backwards_grnd = 0.;
	}

	/* convert H2 into H- 
	 * the back reaction, H2(exct) + e => H- + Ho */
	if( hmi.rel_pop_LTE_H2s > 0. )
	{ 

		hmi.assoc_detach_backwards_exct = hmi.assoc_detach*hmi.rel_pop_LTE_Hmin/hmi.rel_pop_LTE_H2s*
			dense.eden*co.lgUMISTrates; 
	}
	else
	{
		hmi.assoc_detach_backwards_exct = 0.;
	}

	{
		/* often the H- route is the most efficient formation mechanism for H2,
		 * will be through rate called Hmolec_old[ipMH]*hmi.assoc_detach
		 * this debug print statement is to trace h2 oscillations */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC && nzone>140/*&& iteration > 1*/)
		{
			/* rapid increase in H2 density caused by rapid increase in hmi.rel_pop_LTE_H2g */
			fprintf(ioQQQ,"hmi.assoc_detach_backwards_grnd\t%.2f\t%.5e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n", 
				/* total forward rate */
				fnzone,
				phycon.te, 
				dense.eden,
				/* rate H- + H => H2 + E */
				hmi.assoc_detach,
				hmi.assoc_detach_backwards_grnd, 
				hmi.assoc_detach_backwards_exct, 
				hmi.hminus_rad_attach,
				hmi.HMinus_induc_rec_rate,
				/* H0 */
				hmi.Hmolec[ipMH],
				/* H+ */
				hmi.Hmolec[ipMHp],
				/* H- */
				hmi.Hmolec[ipMHm],
				hmi.H2_total,
				hmi.rel_pop_LTE_Hmin,
				hmi.rel_pop_LTE_H2g,
				hmi.rel_pop_LTE_H2s
				);
		}
	}

	/* >>chng 03 sep 11, resolve H2 and H2*, use fraction determined above */
	if(r->next == NULL) {
		int in[]={ipMH2g},out[]={ipMH,ipMHm};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;

	/* >>chng 05 aug 05, NPA comment:The Leiden comparison does not consider H-.  
	Therefore the UMIST hack is used to turn off H- so that it is never important */

	/* >>chng 05 oct 03, TE, rearrange to get the correct H2 destruction file */ 
	hmi.assoc_detach_backwards_grnd *= (frac_H2_hminus * co.lgUMISTrates);
	r->rk = hmi.assoc_detach_backwards_grnd;

	/* >>chng 03 sep 11, resolve H2 and H2*, add new destruction process for H2* */
	if(r->next == NULL) {
		int in[]={ipMH2s},out[]={ipMH,ipMHm};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	/* >>chng 04 jan 28, had wrong Boltzmann factor for this reaction,
	 * fixed by Gargi Shaw */

	/* >>chng 05 aug 05, NPA comment:The Leiden comparison does not consider H-.  
	Therefore the UMIST hack is used to turn off H- so that it is never important */

	/* >>chng 05 oct 03, TE, rearrange to get the correct H2 destruction file */ 
	hmi.assoc_detach_backwards_exct *= (frac_H2star_hminus * co.lgUMISTrates);
	r->rk = hmi.assoc_detach_backwards_exct;

	/*#	define Hneut 7e-8*/
	/* >>chng 05 sept 12 - NPA.  change rate for mutual neutralization of H-
	 * and H+ to the rate from 
	 * >>refer H-	mutual neut	Lepp, S., Stancil, P.C. & Dalgarno, A. 2002, J. Phys. B, 35, R57 */
	if( phycon.te < 14125. )
	{
		/* the fit in Lepp et al. explodes at high temperature,
		 * Te = 14,125 is the temp where the rates reaches its lowest value */
		Hneut = 1.4e-7*pow(phycon.te/300,-0.487)*exp(phycon.te/29300);
	}
	else
	{
		Hneut = 3.4738192887404660e-008;
	}
	/* mutual neut, mostly into n=3; rates from Janev et al
	 * H- + H+ => H + H(n=3) */
	/** \todo	2	process is net source term for H(n=3) states, must be added in */
	fhneut = Hmolec_old[ipMHp]*Hneut; /* dense.xIonDense[ipHYDROGEN][1]*7e-8; */

	if(r->next == NULL) {
		int in[]={ipMHm,ipMHp},out[]={ipMH,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;

	/* >>chng 05 aug 05, NPA comment:The Leiden comparison does not consider H-.  
	Therefore the UMIST hack is used to turn off H- so that it is never important */

	r->rk = Hneut*co.lgUMISTrates;

	/* back reaction from excited state H */
	if( phycon.te > 1000. )
	{
		/* HBN(3,1) is defined; when <HydTempLimit then set to 1 */
		bhneut = (Hneut*hmi.rel_pop_LTE_Hmin*dense.eden)*iso.DepartCoef[ipH_LIKE][ipHYDROGEN][3];
	}
	else
	{
		bhneut = 0.;
	}

	/* mutual neut, mostly into n=3; rates from Janev et al
	 * H + H(n=3) => H- + H+ */
	/** \todo	2	process is net ionization term for H(n=3) states */
	/* this is the back reaction, forming H- from Ho */

	if(r->next == NULL) {
		int in[]={ipMH,ipMH},out[]={ipMHm,ipMHp}, ratesp[]={ipMH,ipMHp};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),ratesp,INTSZ(ratesp));
	}
	r = r->next;
	rindex++;
	r->rk = bhneut*co.lgUMISTrates; 
	bhneut *= Hmolec_old[ipMHp];

	/*--------------------------------------------------------------------
	 *
	 * molecular hydrogen H2 Htwo balance equation
	 * (IPMH2,IPMHO)==create H2 from Ho =+ */

	/* H2 formation on grains */
	/* >>chng 01 jan 05, remove from matrix part and add hden to hmi.rate_h2_form_grains_used, */
	/* the large molecule keeps explicit track of the fraction that goes into 
	 * excited vs ground H2.  Use that ratio if H2 turned on, else use an
	 * estimate of it */
	/* The reaction rate is only proportional to one of the ipMH, due to
	 * surface saturation (?) */

#	define	CATALYST	true
	if( CATALYST )
	{
		/* this is the method used by the code for most of its history.  The grain
		 * is only a catalytic agent, and so the rate goes as the square of the
		 * incoming H0 density */
		/* This goes to excited H2 */
		if(r->next == NULL) {
			int in[]={ipMH,ipMH},out[]={ipMH2s},ratesp[]={ipMH};
			r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),ratesp,INTSZ(ratesp));
		}
		r = r->next;
		rindex++;
		r->rk = gv.rate_h2_form_grains_used_total*frac_H2star_grains;

		/* >>chng 03 sep 10, multiply above by correction for excited state,
		* add below reaction for population of excited state */

		/* This goes to ground H2 */
		if(r->next == NULL) {
			int in[]={ipMH,ipMH},out[]={ipMH2g},ratesp[]={ipMH};
			r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),ratesp,INTSZ(ratesp));
		}
		r = r->next;
		rindex++;
		r->rk = gv.rate_h2_form_grains_used_total*frac_H2_grains;
	}

	else
	{
		/* >>chng 03 nov 25, go to this formalism */
		/* the grain is not a true catalyst, but rather a target loaded with H atoms
		 * ready to react.  So the rate is the number of these grains, times their
		 * cross section, times the number of incident H atoms.  The number of grains
		 * is replaced with the total hydrogen density, which is not in the network
		 * but is a constant */
		/* This goes to excited H2 */
		if(r->next == NULL) {
			int in[]={ipMH,ipMH},out[]={ipMH2s},ratesp[]={ipMH};
			r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),ratesp,INTSZ(ratesp));
		}
		r = r->next;
		rindex++;
		if( dense.xIonDense[ipHYDROGEN][0]>0 )
			r->rk = gv.rate_h2_form_grains_used_total*frac_H2star_grains*
				dense.gas_phase[ipHYDROGEN] / dense.xIonDense[ipHYDROGEN][0];
		else
			r->rk = 0.;

		/* >>chng 03 sep 10, multiply above by correction for excited state,
		* add below reaction for population of excited state */

		/* This goes to ground H2 */
		if(r->next == NULL) {
			int in[]={ipMH,ipMH},out[]={ipMH2g},ratesp[]={ipMH};
			r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),ratesp,INTSZ(ratesp));
		}
		r = r->next;
		rindex++;
		if( dense.xIonDense[ipHYDROGEN][0]>0 )
			r->rk = gv.rate_h2_form_grains_used_total*frac_H2_grains*
				dense.gas_phase[ipHYDROGEN]/dense.xIonDense[ipHYDROGEN][0];
		else
			r->rk = 0.;
	}

	/* excited atom radiative association,
	 * H(n=2) + H(n=1) => H2 + hnu
	 * written as H(n=1)*pop ratio + H(n=1) -> H2 + hnu but ratio of pops is in
	 * terms of pop to ion, so initial term is
	 * n(H+) * pop2ion, 
	 * >>refer	H2	rates	Latter, W.B., & Black, J.H., 1991, Ap.J. 372, 161 */
	/* hmi.radasc = ((StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH2p].Pop + StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH2s].Pop)*dense.xIonDense[ipHYDROGEN][1])*3e-14; */

	/** \todo	1	equivalent reaction for H2* is not included in chemistry, Big h2 does not include this reaction, what to do? GS */
	if( dense.xIonDense[ipHYDROGEN][1] > 0. )
		hmi.radasc = (StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH2p].Pop + 
			StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH2s].Pop)/dense.xIonDense[ipHYDROGEN][1]*3e-14;
	else
		hmi.radasc = 0.;
	/* >>chng 02 nov 7 rjrw: correct for n^2 behaviour w.r.t. H 
	   >>chng 02 nov 7 rjrw, correct stoichiometry */

	/* Possible that changing to a rate proportional to ipMHp would be more consistent */
	if(r->next == NULL) {
		int in[]={ipMH,ipMH},out[]={ipMH2g},ratesp[]={ipMH,ipMHp};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),ratesp,INTSZ(ratesp));
	}
	r = r->next;
	rindex++;
	r->rk = hmi.radasc*co.lgUMISTrates;  
	hmi.radasc *= Hmolec_old[ipMHp]; /* why this is reset here? GS*/

	/* photo-destroy H2 */
	/* >>chng 00 nov 25 factor of 0.1, assume pump is total, and 10% destroy H2 is 21*/
	if(r->next == NULL) {
		int in[]={ipMH2g},out[]={ipMH,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	/* >>chng 03 mar 07, had factor of 0.1 for branching ratio from H2** to H+H, 
	 * but branching is now already included */
	/*r->rk = hmi.H2_Solomon_dissoc_rate_used*0.1;*/
	r->rk = hmi.H2_Solomon_dissoc_rate_used_H2g;

	/* >>chng 03 sep 11, add this process */
	/* photo-destroy H2* by Solomon process at same rate as H2ground dissociation,
		see above eqn A12 in TH85 */
	/* >>chng 00 nov 25 factor of 0.1, assume pump is total, and 10% destroy H2 */
	if(r->next == NULL) {
		int in[]={ipMH2s},out[]={ipMH,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	/* >>chng 03 mar 07, had factor of 0.1 for branching ratio from H2** to H+H, 
	 * but branching is now already included */
	/*r->rk = hmi.H2_Solomon_dissoc_rate_used*0.1; is #22*/
	r->rk = hmi.H2_Solomon_dissoc_rate_used_H2s;

	/* H2 + H+ => H3+ HNU
	 * equation H21 from 
	 * >>refer	H2	chemistry	Galli,D., & Palla, F. 1998, A&A, 335, 403-420 */
	/* >>chng 02 nov 07 rjrw, include H+ ion density in rate constant */

	/* >>chng 05 aug 05, NPA comment.  This reaction is not in UMIST, therefore I turned it 
	   off when comparing to the other codes */
	/** \todo	1	equivalent reaction for H2* is not included in chemistry, Big h2 does not include this reaction, what to do? GS */

	hmi.h2hph3p = 1.0e-16f*co.lgUMISTrates;

	if(r->next == NULL) {
		int in[]={ipMH2g,ipMHp},out[]={ipMH3p};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = hmi.h2hph3p;

	/* collisional dissociation, rate from 
	 * >>refer	H2	collisional dissociation	Dove, J.E., and Mandy, M. E., 1986, ApJ, 311, L93.
	 * H_2 + H => 2H + H 	 
	 * >>chng 02 nov 7 rjrw, correct stoichiometry */

	/* Rate is catalyzed by an additional H */
	if(r->next == NULL) {
		int in[]={ipMH2g},out[]={ipMH,ipMH},ratesp[]={ipMH,ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),ratesp,INTSZ(ratesp));
	}
	r = r->next;
	rindex++;
	hmi.rh2dis *= co.lgUMISTrates;
	r->rk = hmi.rh2dis;

	/* >>chng 04 apr 21 */
	/* 2H + H2 => H2 + H2
	 * equation (5) from 
	 * >>refer	H2	chemistry Palla, F., Salpeter, E.E., & Stahler, S.W., 1983, ApJ,271, 632-641
	 * bh2h22hh2= 5.5e-29/(8*phycon.te) */
	/** \todo	1	equivalent reaction for H2* is not included in chemistry, Big h2 does not include this reaction, what to do? GS */
	hmi.bh2h22hh2 = 5.5e-29*co.lgUMISTrates/(8.*phycon.te);

	if(r->next == NULL) {
		int in[]={ipMH,ipMH,ipMH2g},out[]={ipMH2g,ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = hmi.bh2h22hh2;

	/* >>chng 04 apr 21 */
	/* H2 + H2 => 2H + H2
	 * equation (5) from 
	 * >>refer	H2	chemistry Palla, F., Salpeter, E.E., & Stahler, S.W., 1983, ApJ,271, 632-641
	 * h2h22hh2 = bh2h22hh2/hmi.rel_pop_LTE_H2g */
	/** \todo	1	equivalent reaction for H2* is not included in chemistry, Big h2 does not include this reaction, what to do? GS */
	if( hmi.rel_pop_LTE_H2g > 0. )
	{
		hmi.h2h22hh2 = hmi.bh2h22hh2/hmi.rel_pop_LTE_H2g;
	}
	else
	{
		hmi.h2h22hh2 =0.;
	}

	if(r->next == NULL) {
		int in[]={ipMH2g,ipMH2g},out[]={ipMH,ipMH,ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	hmi.h2h22hh2 *= co.lgUMISTrates;
	r->rk = hmi.h2h22hh2;

	/* back rate, three body recombination, 2H + S => H_2 + S 	 
	 * >>chng 02 nov 7 rjrw: correct for n^2 behaviour w.r.t. H 
	 * >>chng 02 nov 7 rjrw, correct stoichiometry 
	 * >>chng 02 nov 7 rjrw, correct for n^3 behaviour w.r.t. H !! */
	/* hmi.bh2dis = hmi.rh2dis*hmi.rel_pop_LTE_H2g*dense.xIonDense[ipHYDROGEN][0]*dense.xIonDense[ipHYDROGEN][0]; */
	/** \todo	1	equivalent reaction for H2* is not included in chemistry, Big h2 does not include this reaction, what to do? GS */
	hmi.bh2dis = hmi.rh2dis*hmi.rel_pop_LTE_H2g*co.lgUMISTrates*hmi.H2_formation_scale;	

	if(r->next == NULL) {
		int in[]={ipMH,ipMH},out[]={ipMH2g},ratesp[]={ipMH,ipMH,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),ratesp,INTSZ(ratesp));
	}
	r = r->next;
	rindex++;
	r->rk = hmi.bh2dis;

	hmi.bh2dis = hmi.rh2dis*hmi.rel_pop_LTE_H2g*Hmolec_old[ipMH]*Hmolec_old[ipMH]*co.lgUMISTrates;

	/* H2 + HNU=>  H2+ + E
	 * photoionization by hard photons, crossection=3*HI - in molecular environments this is
	 * only the highest energy photons that have penetrated the H+ and H0 regions*/
	/* following copies from opacity_addtotal line 353 */
	/** \todo	0	update photoelectric opacity for H2 to include real cross sections and energies.
	 * this is not a higher priority because when H2 is formed there can be very little ionizing
	 * radiation. this process must be trivial compared with the Solomon process 
	 * following reference gives cross section for all energies
	 * >>refer	H2	photo cs	Yan, M., Sadeghpour, H.R., & Dalgarno, A., 1998, ApJ, 496, 1044 
	 * Wilms, J., Allen, A., & McCray, R. 2000, ApJ, 542, 914 */
	/* >>chng 02 jan 16, approximate inclusion of H_2 photoelectric opacity */
	/* include H_2 in total photoelectric opacity */
	/* set lower and upper limits to this range */
	/*hmi.H2_photoionize_rate = iso.gamnc[ipH_LIKE][ipHYDROGEN][ipH1s];
	fprintf(ioQQQ,"DEBUG H2 photo\t%.3e", hmi.H2_photoionize_rate );*/
	/** \todo	0	must include heating, Compton ionization */
	/* >>chng 05 nov 24, evaluate real photo rate,
	 * had used H0 rate - photo heating had not been included */
	{
		static long int nzone_eval = -1, iteration_evaluated=-1;
		/* must reevaluate During search phase */
		if( ( nzone_eval!=nzone || iteration_evaluated!=iteration ) || !nzone )
		{
			/* generally not important, do one time per zone */
			hmi.H2_photoionize_rate = 
				GammaK(opac.ipH2_photo_thresh , 
				rfield.nupper,
				opac.ipH2_photo_opac_offset,1.,
				&photoHeat)*
				ionbal.lgPhotoIoniz_On +
				/* Compton recoil ionization - we include this in the H2 photoionization
				 * rate but not the heating rate - factor of two since assume 2H
				 * is same as two H0 at such high energies */
				 2.*ionbal.CompRecoilIonRate[ipHYDROGEN][0];

			/* photo heating - this has units s-1 - needs H2 density
			 * to become vol heat rate */
			hmi.H2_photo_heat_soft = photoHeat.HeatLowEnr * ionbal.lgPhotoIoniz_On;
			hmi.H2_photo_heat_hard = photoHeat.HeatHiEnr * ionbal.lgPhotoIoniz_On;
			nzone_eval = nzone;
			iteration_evaluated = iteration;
		}
	}

	/*fprintf(ioQQQ,"\t %.3e\n", hmi.H2_photoionize_rate );*/

	/* cosmic rays predominantly H2 + cr -> H2+ + e, as per table 10 of TH85 */
	/* >>chng 00 nov 28, factor of 0.93 from
	 >>refer	cosmic ray	ionization rate	Tielens, A.G.G.M., & Hollenbach, D., 1985, ApJ, 291, 722
	 * also cosmic rays producing secondary ionization csupra */
	/* >>chng 00 nov 28, factor of 0.93 from
	 >>refer	cosmic ray	ionization rate	Maloney, P.R., Hollenbach, D., & Tielens, A. G. G. M., 1998, ApJ, 466, 561
	 */
	/* >>chng 04jan 26, assume ion(H2) = 2x ion(H), H ion rate of
	 * 2.5e-17 s-1, as per
	 * >>refer	cosmic ray	ionization	Williams, J.P., Bergin, E.A., Caselli, P., 
	 * >>refercon	Myers, P.C., & Plume, R. 1998, ApJ, 503, 689 
	 * so H2 secondary ionization rate is 5e-17 s-1 */
	if(r->next == NULL) {
		int in[]={ipMH2g},out[]={ipMH2p};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	/* >> chng 05 jul 07, TE, rename to get correct h2 destruction file */
	/* ratio of H2 to H cr rates from table 10 of tielens & hollenbach 1985 */
	/* >> chng 05 aug 05, NPA comment.  Our definition of the cosmic ray reaction, 
	   we include the factor hmi.H2_photoionize_rate.  The Leiden comparison wanted a constant cosmic 
	   ray rate.  Therefore if the UMIST rate is set we use a constant 4.4e-17 ionization
	   rate.  Otherwise we just use what Cloudy naturally does */
	if(co.lgUMISTrates)
	{
		h2crh2pe = hmi.H2_photoionize_rate + secondaries.csupra[ipHYDROGEN][0]*2.02;
	}

	else
	{
		h2crh2pe = 4.4e-17;
	}

	r->rk = h2crh2pe;

	/* >>chng 04 apr 22, add H2 + cr -> H+ H + e, TH85 table 10 */
	if(r->next == NULL) {
		int in[]={ipMH2g},out[]={ipMH,ipMHp};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;

	/* >> chng 05 jul 07, TE, rename to get correct h2 destruction file */
	/* ratio of H2 to H cr rates from table 10 of tielens & hollenbach 1985 */
	/* >> chng 05 aug 05, NPA comment.  The Leiden comparison wanted a constant cosmic 
	   ray rate.  Therefore if the UMIST rate is set we use a constant 1e-19 ionization
	   rate.  Otherwise we just use what Cloudy naturally does */

	if(co.lgUMISTrates)
	{
		h2crhphe = secondaries.csupra[ipHYDROGEN][0]*0.0478;
	}
	else
	{
		h2crhphe =  1e-19;
	}
	r->rk = h2crhphe;

	/* >> chng 05 sep 26, TE, include the same reaction for H2s */
	/* H2s + CR -> H+ H + e, TH85 table 10 */
	if(r->next == NULL) {
		int in[]={ipMH2s},out[]={ipMH,ipMHp};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	if(co.lgUMISTrates)
	{
		h2scrhphe = secondaries.csupra[ipHYDROGEN][0]*0.0478;
	}
	else
	{
		h2scrhphe =  1e-19;
	}
	r->rk = h2scrhphe;


	/* >> chng 05 jul 07, TE, rename to get correct h2 destruction file */
	/* >>chng 05 jul 01,GS */
	/* H2s + CRP => H2+ + e; 
	 * Cosmic ray ionization of H2s added*/
	if(r->next == NULL) {
		int in[]={ipMH2s},out[]={ipMH2p};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;

	if( co.lgUMISTrates )
	{
		/* not using UMIST - do the general case */
		h2scrh2pe = hmi.H2_photoionize_rate + secondaries.csupra[ipHYDROGEN][0]*2.02;
	}
	else
	{
		/* use UMIST - this is from Sternberg email defining Leiden meeting */
		h2scrh2pe = 4.4e-17;
	}
	r->rk = h2scrh2pe;

	/* >>chng 03 apr 11 */
	/* H2 + CRP => H + H; CRP=Cosmic Ray proton 
	 * equation (3643) from
	 * >>refer	H2	k	Millar, T.J. et.al, 1997,A&AS, 121, 139
	 * h2crphh = 1.3e-18f */ 
	if(r->next == NULL) {
		int in[]={ipMH2g},out[]={ipMH,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;

	/* >>chng 05 jun 16, GS, use the rate from big H2 network */
	if( h2.lgH2ON  && hmi.lgBigH2_evaluated && hmi.lgH2_Chemistry_BigH2 )
	{
		h2crphh = hmi.H2_tripletdissoc_H2g;
	}
	else
	{
		h2crphh = secondaries.x12tot*3.;
	}

	/* co.lgUMISTrates is set false with the set Leiden hack command, which also
	 * sets their standard cosmic ray rates */
	/* >> chng 05 aug 05, NPA comment.  The Leiden comparison wanted a constant cosmic 
	   ray rate.  Therefore if the UMIST rate is set we use a constant 5e-18 ionization
	   rate.  Otherwise we just use what Cloudy naturally does */
	if(!co.lgUMISTrates)
		 h2crphh = 5e-18;

	r->rk = h2crphh;

	/* >>chng 05 jun 16, GS, use the rate from big H2 network, small network does not have this rate */
	if( h2.lgH2ON  && hmi.lgBigH2_evaluated && hmi.lgH2_Chemistry_BigH2 )
	{
		h2scrphh = hmi.H2_tripletdissoc_H2s;
	}
	else
	{
		h2scrphh = secondaries.x12tot*3.;
	}

	if(r->next == NULL) {
		int in[]={ipMH2s},out[]={ipMH,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = h2scrphh;

	/* >>chng 03 apr 11 */
	/* H2 + CRP => H+ + H_; CRP=Cosmic Ray Proton
	 * equation (3644) from
	 * >>refer	H2	k	Millar, T.J., et.al, 1997,A&AS, 121, 139
	 * h2crphphm = 3.9e-21 */
	/* >> chng 05 aug 05, NPA comment.  Turn off H- for the Leiden comparison */
	 h2crphphm = 3.9e-21 * hextra.cryden_ov_background * co.lgUMISTrates;

	if(r->next == NULL) {
		int in[]={ipMH2g},out[]={ipMHp,ipMHm};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;

	r->rk = h2crphphm;

	/* >> chng 05 sep 26, TE, include the same reaction for H2s */
	/* H2s + CRP => H+ + H_; CRP=Cosmic Ray Proton */
	 h2scrphphm = 3.9e-21 * hextra.cryden_ov_background * co.lgUMISTrates;

	if(r->next == NULL) {
		int in[]={ipMH2s},out[]={ipMHp,ipMHm};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;

	r->rk = h2scrphphm;

	/* >>chng 03 apr 11 */
	/* H2 + CRP => H+ + H + e; CRP=Cosmic Ray Proton
	 * equation (3641) from
	 * >>refer	H2	k	Millar, T.J., et.al, 1997,A&AS, 121, 139
	 * h2crphpeh = 2.2e-19f */ 
	/* >> chng 05 aug 05, NPA comment.  Amiel Sternberg said to not consider this process
	      in the benchmark calculations.  Therefore, the UMIST rate is used to turn this reaction off */
	 h2crphpeh = 2.2e-19 * hextra.cryden_ov_background * co.lgUMISTrates;

	if(r->next == NULL) {
		int in[]={ipMH2g},out[]={ipMHp,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = h2crphpeh;

	/* >> chng 05 sep 26, TE, include the same reaction for H2s */
	/* H2s + CRP => H+ + H + e; CRP=Cosmic Ray Proton */
	 h2scrphpeh = 2.2e-19 * hextra.cryden_ov_background * co.lgUMISTrates;

	if(r->next == NULL) {
		int in[]={ipMH2s},out[]={ipMHp,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = h2scrphpeh;

	/*>>chng 05 jul 01, GS, h2s ionization by cosmic ray added*/
	/* >>chng 05 jun 29, TE, used in new save H2 destruction file*/
	hmi.CR_reac_H2g = h2crh2pe + h2crhphe + h2crphh + h2crphphm + h2crphpeh;
	hmi.CR_reac_H2s = h2scrh2pe + h2scrhphe + h2scrphh + h2scrphphm + h2scrphpeh;

	/* >>chng 03 apr 11 */
	/* H3+ + H-=> H2 + H2; 
	 * equation (5,table 9) from
	 * >>refer	H3+	k	Maloney et.al, 1996,ApJ, 466, 561
	 * h3phm2h2 = 1.3e-7f*pow(phycon.te/1000., -0.5) */ 
	/** \todo	1	equivalent reaction for H2* is not included in chemistry, Big h2 does not include this reaction, what to do? GS */
	/* >> chng 05 aug 05, NPA comment.  Turn off H- for the Leiden comparison */
	hmi.h3phm2h2 = 1.3e-7 / (phycon.sqrte/31.62278) * co.lgUMISTrates;/*pow(phycon.te/1000., -0.5);*/
	if(r->next == NULL) {
		int in[]={ipMH3p,ipMHm},out[]={ipMH2g,ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = hmi.h3phm2h2;

	/* >>chng 03 sep 30 */
	/* H3+ + HNU=> H2+ + H; 
	 * equation (table 5) from
	 * >>refer	H2	dissoc	Tielens, A.G.G.M., & Hollenbach, D., 1985, ApJ, 291, 722  
	 * h3ph2ph = 7.9e-9*hmi.UV_Cont_rel2_Habing_TH85_depth;*/
	/* >>chng 04 jun 13 --  update this rate to match that in the UMIST database */

	/* >> chng 05 aug 05, NPA comment.  This is one of the few instances where 
	I actually updated a rate permanently.  Originally this rate was taken from TH85.  
	I changed this to the UMIST rate.  The UMIST database gives the reference for their 
	rate as van Dishoeck and Black, 1987.  Therefore, since it comes from such a regarded 
	source and after the TH85 paper, I think this is appropriate.  The UMIST hack is used 
	to extinguish the continuum by exp(-a*AV) instead of through our radiative transfer 
	solution, which is what the other codes in the comparison did. */

	if(co.lgUMISTrates)
	{
		h3ph2ph = 5.0e-13*hmi.UV_Cont_rel2_Habing_TH85_depth/1.66f;
	}
	else
	{
		h3ph2ph = 5.0e-13*hmi.UV_Cont_rel2_Habing_TH85_face*(realnum)sexp((2.3*rfield.extin_mag_V_point))/1.66f;
	}

	if(r->next == NULL) {
		int in[]={ipMH3p},out[]={ipMH2p,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =h3ph2ph;

	/* >>chng 03 sep 30 */
	/* H3+ + HNU=> H2 + H+; 
	 * equation (table 5) from
	 * >>refer	H2	dissoc	Tielens, A.G.G.M., & Hollenbach, D., 1985, ApJ, 291, 722  
	 * h3ph2hp = 2.0e-8*hmi.UV_Cont_rel2_Habing_TH85_depth;*/
	/* >>chng 04 jun 13 --  update this rate to match that in the UMIST database */	
	/* >> chng 05 aug 05, NPA comment.  This is one of the few instances where 
	I actually updated a rate permanently.  Originally this rate was taken from TH85.  
	I changed this to the UMIST rate.  The UMIST database gives the reference for their 
	rate as van Dishoeck and Black, 1987.  Therefore, since it comes from such a regarded 
	source and after the TH85 paper, I think this is appropriate.  The UMIST hack is used 
	to extinguish the continuum by exp(-a*AV) instead of through our radiative transfer 
	solution, which is what the other codes in the comparison did. */
	/** \todo	1	equivalent reaction for H2* is not included in chemistry, Big h2 does not include this reaction, what to do? reverse of this reaction i not in detailed balance,why? GS */
	if(co.lgUMISTrates)
	{
		hmi.h3ph2hp = 5.0e-13*hmi.UV_Cont_rel2_Habing_TH85_depth/1.66f;
	}
	else
	{
		hmi.h3ph2hp = 5.0e-13*hmi.UV_Cont_rel2_Habing_TH85_face*sexp((1.8*rfield.extin_mag_V_point))/1.66;
	}

	if(r->next == NULL) {
		int in[]={ipMH3p},out[]={ipMH2g,ipMHp};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =hmi.h3ph2hp;

	/* >> chng 02 nov 15 rjrw: multiply c[ipMHo][*] terms by ionization fraction
	 * as b[ipMHo] contains _both_ H0 and H+ */
	/* H2*  +  H+  =>  H2+  +  H */
	/* >> chng 05 jul 14, TE, 
	 * to maintain detailed balance with bh2h2p, only consider H2s*/
	/* >> chng 05 sept 28, GS, 
	 * H2g  +  H+  =>  H2+  +  H 
	 *  bh2h2p and rh2h2p are not in detailed balance,astro-ph/0404288*/
	if(r->next == NULL) {
		int in[]={ipMH2g,ipMHp},out[]={ipMH,ipMH2p};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = hmi.rh2h2p;

	/* (3,IPMH2P) == destroy H2+ = -sign
	 * H + H2+ => H+ + H2 */
	/* >>chng 02 nov 7 rjrw, remove the destruction rate 
	 *	 c[ipMHo][ipMHo] += -hmi.bh2h2p*hmi.Hmolec[ipMH2p]; 
	 * twice -- reaction changes state of H within single [H0,H+] `species' */

	if(r->next == NULL) 
	{
	/* >> chng 05 jul 13, TE, 
	 * this process populates v=4,no J information assume into J=0 -> H2s not H2g */
		int in[]={ipMH,ipMH2p},out[]={ipMHp,ipMH2s};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = hmi.bh2h2p;

	/** \todo	2	this rate drives numerical instability in such models as secondary1 and 2.in */
	/* this rate couples H2+ and H3+, and tends to destabilize the matrix in both highly
	 * ionized and fully molecular conditions.  Setting this to zero had no effect - the th85
	 * predictions were identical.  
	 *
	 */
	/* H + H3+ => H2 + H2+ */
	/** \todo	1	equivalent reaction for H2* is not included in chemistry, Big h2 does not include this reaction, what to do? GS */
	/* >> chng 05 aug 05, NPA comment.  This rate is not in UMIST.  Therefore 
	      it is turned off for the Leiden comparison */

	hmi.h3ph2p = HMRATE(2.08e-9,0.,1.88e4)*co.lgUMISTrates;

	if(r->next == NULL) 
	{
		int in[]={ipMH,ipMH3p},out[]={ipMH2g,ipMH2p};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = hmi.h3ph2p;

	/* >>chng 03 feb 7 */
		/* H3+ + H- => H2 + H + H
	 * equation (50) from
	 * >>refer	H-	k	Stancil, P.C, & Lepp, S, & Dalgarno, A. 1998,ApJ, 509, 1-10
	 * h3phmh2hh = 2.3e-7f*pow(phycon.te/300.0, -0.5) */
	/* >> chng 05 aug 05, NPA comment.  This rate is not in UMIST.  Therefore it 
	   is turned off for the Leiden comparison */
	/** \todo	1	equivalent reaction for H2* is not included in chemistry, Big h2 does not include this reaction, what to do? GS */
	hmi.h3phmh2hh = 2.3e-7f*pow(phycon.te/300 , -0.5)*co.lgUMISTrates;
	if(r->next == NULL) 
	{
		int in[]={ipMH3p,ipMHm},out[]={ipMH2g,ipMH,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = hmi.h3phmh2hh;

	/* H2 + H3+ => H2 + H2+ + H */
	/* >> chng 05 aug 05, NPA comment.  This rate is not in UMIST.  Therefore it 
	is turned off for the Leiden comparison */
	/** \todo	2	equivalent reaction for H2* is not included in chemistry, Big h2 does not include this reaction, what to do? GS */
	hmi.h3petc = HMRATE(3.41e-11,0.5,7.16e4)*co.lgUMISTrates;

	if(r->next == NULL) 
	{
		int in[]={ipMH2g,ipMH3p},out[]={ipMH2g,ipMH2p,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = hmi.h3petc;

	/* H2 + H3+ => H2 + H+ + H2 */
	/* >> chng 05 aug 05, NPA comment.  This rate is not in UMIST.  Therefore it 
	is turned off for the Leiden comparison */
	/** \todo	2	equivalent reaction for H2* is not included in chemistry, Big h2 does not include this reaction, what to do? GS */
	hmi.h32h2 = HMRATE(3.41e-11,0.5,5.04e4)*co.lgUMISTrates;

	if(r->next == NULL) {
		int in[]={ipMH2g,ipMH3p},out[]={ipMHp,ipMH2g,ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = hmi.h32h2;

	/* e + H3+ => H2 + H */
	/* e + H3+ => 3H was supposed to be included in this rate, 
	 * and stoichiometric factor 2* on sink rate seemed wrong */
	/* >>chng 03 feb 10, increase rate by factor of 13.6 to agree with
	 * >>refer	H3+	DR	McCall, B.J., et al. 2003, Nature, in press (astro-ph 0302106)*/
	/* >>chng 03 feb 13, extra 0.2 since 20% of these go to H2 + H, Stancil private comm */
	/* >>chng 04 apr 22 , update the next two rates to match that of:
	   >>refer	H3+	k	Stancil, P. C., Lepp, S., and Dalgarno, A 509, 1-10;
	   *>>refercon	Table 1, reactions #48 and #49 */
	/** \todo	1	equivalent reaction for H2* is not included in chemistry,
	 * Big h2 does not include this reaction, what to do? GS */
	/* >>chng 06 jan 23, Stancil's rate is rescaled by 2.25 to match McCall's rate, GS*/
	/* >>chng 07 jan 05, as per GS discussions, USE_MCCALL was commented out,
	 * turn back on.  Should not have been commented out */
#	define USE_MCCALL
#	ifdef USE_MCCALL
#	define	FACTOR	2.25
#	else
#	define	FACTOR	1.0
#	endif
	hmi.eh3_h2h = HMRATE(4.00e-8/FACTOR,-0.5,0.)*dense.eden;

	/* >> chng 05 aug 05, NPA comment.  Rate we use and UMIST uses is different.  If UMIST 
	   hack is on then we use UMIST, otherwise we use our rate   */
	if(!co.lgUMISTrates)
		hmi.eh3_h2h = HMRATE(2.5e-8,-0.3,0.)*dense.eden;


	if(r->next == NULL) {
		int in[]={ipMH3p},out[]={ipMH,ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = hmi.eh3_h2h;

	/* e + H3+ => 3H */
	/* >>chng 06 jan 23, Stancil's rate is rescaled by 2.25 to match McCall's rate, GS*/
	eh3p_3h = HMRATE(1.6e-7/FACTOR,-0.5,0.)*dense.eden;
#	undef FACTOR
#	ifdef USE_MCCALL
#	undef USE_MCCALL
#	endif

		/* >> chng 05 aug 05, NPA comment.  Rate we use and UMIST uses is different.  If UMIST 
	   hack is on then we use UMIST, otherwise we use our rate   */
	if(!co.lgUMISTrates)
		eh3p_3h = HMRATE(7.5e-8,-0.3,0.)*dense.eden;

	/* >>chng 03 feb 10, increase rate by factor of 13.6 to agree with
	 * >>refer	H3+	DR	McCall, B.J., et al. 2003, Nature, in press (astro-ph 0302106)*/
	/* >>chng 03 feb 13, extra 0.8 since 80% of these go to 3H, Stancil private comm */

	if(r->next == NULL) {
		int in[]={ipMH3p},out[]={ipMH,ipMH,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = eh3p_3h;

	if( (trace.lgTrace && trace.lgTr_H2_Mole) )
	{
		if( hmi.H2_rate_destroy > SMALLFLOAT )
		{
			fprintf( ioQQQ, 
			  " H2 destroy rate=%.2e DIS;%.3f bat;%.3f h2dis;%.3f hmi.H2_photoionize_rate;%.3f h2h2p;%.3f E-h;%.3f hmi.h2hph3p;%.3f sec;%.3f\n", 
			  hmi.H2_rate_destroy, 
			  hmi.H2_Solomon_dissoc_rate_used_H2g / hmi.H2_rate_destroy, 
			  hmi.assoc_detach_backwards_grnd / hmi.H2_rate_destroy, 
			  hmi.rh2dis*dense.xIonDense[ipHYDROGEN][0] / hmi.H2_rate_destroy, 
			  hmi.H2_photoionize_rate / hmi.H2_rate_destroy, 
			  hmi.rh2h2p*dense.xIonDense[ipHYDROGEN][1] / hmi.H2_rate_destroy, 
			  hmi.eh2hh /hmi.H2_rate_destroy, 
			  hmi.h2hph3p / hmi.H2_rate_destroy ,
			  secondaries.csupra[ipHYDROGEN][0]*2.02 / hmi.H2_rate_destroy
			  );
		}
		else
		{
			fprintf( ioQQQ, " Destroy H2: rate=0\n" );
		}
	}

	/*------------------------------------------------------------------- */

	/* h2plus H2+ balance equations */

	/** \todo	2	must add process H2+ + H- => H2 + H, Dalgarno&Lepp 87 */
	/* >>refer	H2+	chemistry	Dalgarno, A., & Lepp, S., 1987, in Astrochemistry, eds. 
	 * >>refercon	M.S. Vardya & S.P. Tarafar, Reidel, Dordrecht, p 109 */
	/* rate = 5e-7 * sqrt(100. / phycon.te); */

	/** \todo	2	put in H2+ + gamma => H + H+ */
	/* >>refer	H2+	chemistry	Stancil, P.C., 1994, ApJ, 430, 360 */
	/* cross section is log10( cs_25) = -1.6547717e6 + 1.8660333e5 ln(nu) - 7.8986431e3*ln(nu)^2
	 * 148.73693 * ln(nu)^3 - 1.0513032*ln(nu)^4 */

	/* make H2+ from Ho
	 * H+  + H  =>  H2+ + HNU
	 * approximation was from Kurucz thesis, not meant for hot gas 
	 * >>chng 02 nov 7 rjrw, stoichiometric factor */
	radath = MAX2(0.,2.325*MIN2(5000.,phycon.te)-1375.)*1e-20;

	/* >> chng 05 aug 05, NPA comment.  Rate we use and UMIST uses is different.  If UMIST 
	   hack is on then we use UMIST, otherwise we use our rate   */

	if( !co.lgUMISTrates)
		radath = HMRATE(5.3e-19,1.85,0);

	if( r->next == NULL) {
		int in[]={ipMH,ipMHp},out[]={ipMH2p};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	/* Printf("O: %ld %g %g\n",rindex,radath*dense.xIonDense[ipHYDROGEN][1]*dense.xIonDense[ipHYDROGEN][0],
		 radath*dense.xIonDense[ipHYDROGEN][1]); */
	rindex++;
	r->rk = radath;

	/* H2+  +  H+  => H + H+ + H+; Janev et al. 3.2.6 */
	/* >>chng 02 nov 7 rjrw, stoichiometric factor */

	/* >> chng 05 aug 05, NPA comment.  This reaction is not in UMIST, so turn if off 
	   for the comparison */

	h2pion = 2.4e-27*POW3(phycon.te)*co.lgUMISTrates;

	if(r->next == NULL) {
		int in[]={ipMH2p},out[]={ipMH,ipMHp},ratesp[]={ipMHp,ipMH2p};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),ratesp,INTSZ(ratesp));
	}
	r = r->next;
	rindex++;
	r->rk = h2pion;

	/* H2+  +  E  => H + H+ + e-; Janev et al. */
	/* >>chng 02 nov 7 rjrw, stoichiometric factor */

	/* >> chng 05 aug 05, NPA comment.  This reaction is not in UMIST, so turn if off 
	 for the comparison */
	h2pcin = 2e-7*sexp(30720./phycon.te)*dense.eden*co.lgUMISTrates;

	if(r->next == NULL) {
		int in[]={ipMH2p},out[]={ipMH,ipMHp};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = h2pcin;

	/* >>chng 04 jul 06 -- NPA, include H ionization/recombination due
	to O-H charge transfer in the molecular solver.  */
	/* O + H+ => H + O+ */
	if( iteration==1 && conv.lgSearch )
	{
		/* during search for first iteration do not do ct */
		oatomic = 0.;
		oion = 0.; 
	}
	else
	{
		/* same zone and iteration, take mean of old and new abund */
#		define OLD_FRAC	0.0
		oatomic = oatomic*OLD_FRAC + dense.xIonDense[ipOXYGEN][0]*(1.-OLD_FRAC); 
		oion = oion*OLD_FRAC + dense.xIonDense[ipOXYGEN][1]*(1.-OLD_FRAC); 
		/* ionbal.lgHO_ct_chem is normally 1 set to 0 with command
		 * set HO charge transfer ionization, in which case we do H O
		 * charge transfer in ionization solver */
		oatomic *= ionbal.lgHO_ct_chem;
		oion *= ionbal.lgHO_ct_chem;
	}
	/*	oatomic = dense.xIonDense[ipOXYGEN][0];
		oion = dense.xIonDense[ipOXYGEN][1]; */

	if(r->next == NULL) {
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}

	r = r->next;
	rindex++;
	/*r->rk = atmdat.HCharExcIonOf[ipOXYGEN][0]*dense.xIonDense[ipOXYGEN][0]; */
	r->rk = atmdat.HCharExcIonOf[ipOXYGEN][0]*oatomic; 

	/* O+ + H => H+ + O */
	if(r->next == NULL) {
		int in[]={ipMH},out[]={ipMHp};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	/*r->rk = atmdat.HCharExcRecTo[ipOXYGEN][0]*dense.xIonDense[ipOXYGEN][1]; */
	r->rk = atmdat.HCharExcRecTo[ipOXYGEN][0]*oion; 

	/* >>chng 03 aug 22 */
	/* H2+ + e=> H + H; 
	 * equation (6,table 4) from
	 * >>refer	H2	l	Maloney et.al, 1996,ApJ, 466, 561
	 * h2pehh = 2.8e-8f*pow(phycon.te/1000., -0.37) */
	/* >>chng 03 sep 01, rm the pow function */
	/*h2pehh = 2.8e-8f*pow(phycon.te/1000., -0.37);*/
	h2pehh = 2.8e-8*12.882/(phycon.te30*phycon.te07)*dense.eden;

	/* >> chng 05 aug 05, NPA comment.  Rate we use and UMIST uses is different.  
	   If UMIST hack is on then we use UMIST, otherwise we use our rate */

	if(!co.lgUMISTrates)
		h2pehh = HMRATE(1.6e-8,-0.43,0);

	if(r->next == NULL) {
		int in[]={ipMH2p},out[]={ipMH,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =h2pehh;

	/* back reaction, H + H+ + e => h2+ + e */
	/* >>chng 05 aug 02, rm eden umist terms since now in above 
	 * note on units - this will go into rk as a two-body reaction,
	 * so we need to multiply by an extra eden before feeding into the matrix
	 * this is why the eden is left in the h2pcin from above, the 
	 * real back rate coefficient for a three body process would be this
	 * divided by eden */
	b2pcin = h2pcin*hmi.rel_pop_LTE_H2p;
	/* this is the hot reaction at high densities */

	if(r->next == NULL) {
		int in[]={ipMH,ipMHp},out[]={ipMH2p};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = b2pcin;

	/* H2+  +  HNU  =>  H+  + H */

	if(r->next == NULL) {
		int in[]={ipMH2p},out[]={ipMH,ipMHp};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = gamtwo;

	/*>>KEYWORD	H2+ photoionization
	 * photoionization by hard photons, crossection =H0 in high-energy limit 
	 * H2+ + hnu -> H + H+ 
	 * one electron system */
	if(r->next == NULL) {
		int in[]={ipMH2p},out[]={ipMH,ipMHp};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;

	/* >> chng 05 aug 05, NPA comment.  This reaction is not in UMIST, for the case 
	 * of hard photons.  Turn if off for the comparison. 
	 * >>chng 05 nov 27, factor of two had been in front of H photo rate 
	 * by analogy with high-energy limit for H2 - but this is a one-electron
	 * system so not appropriate 
	 * note that iso.gamnc include bound Compton ionization */
	/*r->rk = 2.*iso.gamnc[ipH_LIKE][ipHYDROGEN][ipH1s]*co.lgUMISTrates;*/
	r->rk = iso.gamnc[ipH_LIKE][ipHYDROGEN][ipH1s]*co.lgUMISTrates;

	/* H2 + H2+ => H + H3+ */
	/** \todo	1	equivalent reaction for H2* is not included in chemistry, 
	 * Big h2 does not include this reaction, what to do? GS */
	hmi.h2ph3p = 1.40e-9*(1. - sexp(9940./phycon.te));

	if(!co.lgUMISTrates)
		hmi.h2ph3p = 2.08e-9;

	if(r->next == NULL) {
		int in[]={ipMH2g,ipMH2p},out[]={ipMH,ipMH3p};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = hmi.h2ph3p;

	/* destroy H2+ via H2+ + H2 => H + H+ + H2 */
	/** \todo	2	equivalent reaction for H2* is not included in chemistry, Big h2 does not include this reaction, what to do? GS */			
	/* >> chng 05 aug 05, NPA comment.  Rate we use and UMIST uses is different.  
	   If UMIST hack is on then we use UMIST, otherwise we use our rate */

	h2phhp = 2.41e-12*phycon.sqrte*sexp(30720./phycon.te)*co.lgUMISTrates;

	if(r->next == NULL) {
		int in[]={ipMH2g,ipMH2p},out[]={ipMH,ipMHp,ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = h2phhp;

	/*------------------------------------------------------------------ */

	/* H3+ balance equations*/

	/* photoionization by hard photons, crossection =2*HI (wild guess)
	 * -- rjrw: where do they go??? 
	 * -- H3+ + hv => H2+ + H+ + e, best guess (P. Stancil, priv comm) */

	if(r->next == NULL) {
		int in[]={ipMH3p},out[]={ipMH2p,ipMHp};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;

	/* >> chng 05 aug 05, NPA comment.  This reaction is not in UMIST, for the case 
	   of hard photons.  Turn if off for the comparison. */
	r->rk = 2.*iso.gamnc[ipH_LIKE][ipHYDROGEN][ipH1s]*co.lgUMISTrates;

	/*------------------------------------------------------------------ */

	/* vib excited H2, called H2* balance equations, these closely follow
	 * >>refer	h2	fits	Tielens, A.G.G.M., & Hollenbach, D., 1985a, ApJ 291, 722 */
	/* population of vib-excited H2, from discussion on pp 736-737 of TH85 */

	/* deexcitation rate from upper level, H2* => H2 */

	/* deexc_hneut is H2* + H0 -> H2g + H 
	 * deexc_htwo  is H2* + H2 -> H2g + H2 */
	Boltz_fac_H2_H2star = 1.*sexp( 30172./phycon.te);
	if( h2.lgH2ON  && hmi.lgBigH2_evaluated && hmi.lgH2_Chemistry_BigH2 )
	{
		deexc_htwo = hmi.Average_collH2_deexcit;
		deexc_hneut = hmi.Average_collH_deexcit;
	}
	else
	{
		deexc_htwo = (1.4e-12*phycon.sqrte * sexp( 18100./(phycon.te + 1200.) ))/6.;
		deexc_hneut =  (1e-12*phycon.sqrte * sexp(1000./phycon.te ))/6.;
	}
	/* total H2* -> H2g rate, s-1 */
	H2star_deexcit = hmi.H2_total*deexc_htwo + hmi.Hmolec[ipMH] * deexc_hneut;

	/* H2g -> H2s */
	if( h2.lgH2ON  && hmi.lgBigH2_evaluated && hmi.lgH2_Chemistry_BigH2 )
	{
		H2star_excit = hmi.Average_collH2_excit *hmi.H2_total + 
			hmi.Average_collH_excit*hmi.Hmolec[ipMH];
	}
	else
	{
		H2star_excit = Boltz_fac_H2_H2star * H2star_deexcit;
	}

	/* depopulate H2_star, 2e-7 is spontaneous deexcitation rate,
	 * which also appears in lines where intensity of vib lines is entered into line stack */
	/* H2* + H2g -> H2g + H2g */

	if(r->next == NULL) {
		int in[]={ipMH2s},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;

	/* >>chng 05 jul 11, TE, rename to use in save file*/
	/* >>chng 05 jul 9, GS, use average A calculated from Big H2 */
	if( h2.lgH2ON  && hmi.lgBigH2_evaluated && hmi.lgH2_Chemistry_BigH2 )
	{
		hmi.h2s_sp_decay = hmi.Average_A;
	}
	else
	{
		hmi.h2s_sp_decay = 2e-7;
	}
	r->rk = hmi.h2s_sp_decay;


	if(r->next == NULL) {
		int in[]={ipMH2s},out[]={ipMH2g},ratesp[]={ipMH2s,ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),ratesp,INTSZ(ratesp));
	}
	r = r->next;
	rindex++;
	r->rk = deexc_htwo;

	if(r->next == NULL) {
		int in[]={ipMH2s},out[]={ipMH2g},ratesp[]={ipMH2s,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),ratesp,INTSZ(ratesp));
	}
	r = r->next;
	rindex++;
	r->rk = deexc_hneut;

	/* collisional excitation of vib from ground, 
	 * H2g + H2g -> H2* + H2g, so H2g acts as catalyst
	 * stat weight of ground 1, excit 6, as per TH discussion
	 * this must normally be zero */
	/* H2 producing H2_star */
	/* >>chng 03 sep 11, had been 6, changed to 1 */
	/*Boltz_fac_H2_H2star = 1.*sexp( hmi.H2_BigH2_H2s_av * T1CM / phycon.te);*/
	/* total excitation rate to H2*, s-1, NB - this is also used in the cooling - heating
	 * rate below */
	/*H2star_excit = Boltz_fac_H2_H2star * H2star_deexcit;*/

	if(r->next == NULL) {
		int in[]={ipMH2g},out[]={ipMH2s},ratesp[]={ipMH2g,ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),ratesp,INTSZ(ratesp));
	}
	r = r->next;
	rindex++;
	/* >>chng 05 jul 10, GS, use average collisional rate calculated from Big H2 */
	if( h2.lgH2ON  && hmi.lgBigH2_evaluated && hmi.lgH2_Chemistry_BigH2 )
	{
		r->rk = hmi.Average_collH2_excit;
	}
	else
	{
		r->rk = deexc_htwo*Boltz_fac_H2_H2star;
	}

	if(r->next == NULL) {
		int in[]={ipMH2g},out[]={ipMH2s},ratesp[]={ipMH,ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),ratesp,INTSZ(ratesp));
	}
	r = r->next;
	rindex++;
	/* >>chng 05 jul 10, GS, use average collisional rate calculated from Big H2 */
	if( h2.lgH2ON  && hmi.lgBigH2_evaluated && hmi.lgH2_Chemistry_BigH2 )
	{
		r->rk = hmi.Average_collH_excit;
	}
	else
	{
		r->rk = deexc_hneut*Boltz_fac_H2_H2star;
	}

	/* >>chng 03 aug 28 */
	/* H2* + H => H + H + H
	 * equation  from table 9
	 * >>refer	H2*	k	Tielens, A.G.G.M., & Hollenbach, D., 1985, ApJ, 291, 722*/

	/*	hmi.h2sh = HMRATE(9.8e-12,0.5,2.7e4);*/
	/* >>chng 05 jul 19, TE, update to UMIST rate */
	/* >>chng 05 mar 18, TE, used in new save H2 destruction file*/
	hmi.h2sh = HMRATE(4.67e-7,-1.,5.5e4);

	if(r->next == NULL) {
		int in[]={ipMH2s,ipMH},out[]={ipMH,ipMH,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = hmi.h2sh;


	/*  The following is a general prescription on how to generate chemical rates for H2* from H2.  It
	    is an e-mail from Phillip Stancil received 05 aug 04.  
		The topic of the e-mail involved the reaction H2 + H2 => H2 + 2H and 
		getting reaction rates for H2* + H2 -> H2 + 2H and H2* + H2* -> H2/H2* + 2H from the ground state
		reaction rate, which has a strong temperature dependence.  E-mail inserted by TE, NPA on 05 aug 05 */

	/***********************************************************************************************
		Consider that both H2's are in arbitrary v,j levels and after the collision
		they are in arbitrary levels with the constraint that at least on of  
		the H2's been dissociated or both are dissociated, but both cannot be excited to 
		bound levels (better to use actual data for those cases).

		H2(vj) + H2(v'j') -> H2(v''j'') + H2(v'''j''')


		Then an approximate rate coefficient would be

		k_vj,v'j'->v''j'',v'''j''' = 10^-11 * exp(-beta/kT)

		where beta = (E_vj - E_v''j'') + (E_v'j' - E_v'''j'''),

		but if beta<0, set beta=0. So that the rate never becomes greater than ~10^-11 cm^3/s.
		Here the energies are dissociation energies (i.e., 4.478 eV for vj=00  and 0 for a
		dissociation state).

		So, looking at limits, if vj=v'j'=v''j''=00 and v'''j'''=dissociation  
		state (15,0, say), then we get

		10^-11*exp(-4.478 eV/kT)

		[ H2(0,0) + H2(0,0) -> H2(0,0) + 2H  ]

		if vj=v'j'=00 and v''j''=v'''j'''=dissociation, i.e. 4H is the product, then we get

		10^-11*exp(-2*4.478/kT)

		[ H2(0,0) + H2(0,0) -> 2H + 2H]


		if vj=v'j'=00, v''j'' is an excited bound state, and v'''j''' is  
		dissociation, the number in parentheses is between 4.478 and 2*4.478 
		(i.e., dissociation plus excitation)

		[ H2(0,0) + H2(0,0) -> H2*(v,j) + 2H ]

		However, say we have a case where both H2 in the H2* (1.88 eV) level.

		H2* + H2* -> H2(0,0) + 2H

		then beta = (1.88-4.478)+(1.88-0)=-0.718 eV, so take beta=0 to give k=10^-11.

		This argument considers only asymptotic energies.
		So, it neglects selection rules (which would be needed if you considered state-to-state
		reactions, but not the two-level approximation) and any effect rovibrational overlaps.

		Phillip Stancil

	**************************************************************************************/


	/* >>chng 03 aug 28 */
	/* H2* + H2 => H2 + H + H
	 * equation from table 9
	 * >>refer	H2*	k	Tielens, A.G.G.M., & Hollenbach, D., 1985, ApJ, 291, 722-746*/

	/* >>chng 05 jul 19, TE, update to UMIST rate */
	/* >>chng 05 mar 18, TE, used in new save H2 destruction file*/
	/*	hmi.h2sh2g = HMRATE(9.8e-12,0.5,2.7e4);*/
	/* hmi.h2sh2g = HMRATE(1e-8,0.,8.41e4); */
	/* >>chng 05 aug 05, TE, update to the scheme of Phillip Stancil */
	hmi.h2sh2g = HMRATE(1e-11,0.,2.18e4);
	if(r->next == NULL) {
		int in[]={ipMH2s,ipMH2g},out[]={ipMH2g,ipMH,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	/* >>chng 05 jul 20, GS, TE */
	if( h2.lgH2ON  && hmi.lgBigH2_evaluated && hmi.lgH2_Chemistry_BigH2 )
	{
		hmi.h2sh2g = hmi.Average_collH2s_dissoc; 
	}

	r->rk = hmi.h2sh2g;

	/* >>chng 03 aug 28 */
	/* H2* + H2* => H2 + H + H
	 * equation from table 9
	 * >>refer	H2*	k	Tielens, A.G.G.M., & Hollenbach, D., 1985, ApJ, 291, 722-746*/
	/* >>chng 05 mar 18, TE, used in new save H2 destruction file*/
	/*	hmi.h2sh2sh2g2h = HMRATE(9.8e-12,0.5,0.);*/
	/* >>chng 05 jul 19, TE, update to UMIST rate */
	/* hmi.h2sh2sh2g2h = HMRATE(1e-8,0.,0.); */
	/* >>chng 05 aug 05, TE, update to the scheme of Phillip Stancil */
	hmi.h2sh2sh2g2h = HMRATE(1e-11,0.,0.);
	if(r->next == NULL) {
		int in[]={ipMH2s,ipMH2s},out[]={ipMH2g,ipMH,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	/* >>chng 05 jul 20, GS, TE */
	if( h2.lgH2ON  && hmi.lgBigH2_evaluated && hmi.lgH2_Chemistry_BigH2 )
	{
		hmi.h2sh2sh2g2h = hmi.Average_collH2s_dissoc;
	}

	r->rk = hmi.h2sh2sh2g2h;


	/* >>chng 05 jul 21, TE, GS */
	/* H2* + H2* => H2* + H + H */
	/* hmi.h2sh2sh2s2h = HMRATE(1e-8,0.,0.);*/
	/* >>chng 05 aug 05, TE, update to the scheme of Phillip Stancil */
	hmi.h2sh2sh2s2h = HMRATE(1e-11,0.,2.18e4);
	if(r->next == NULL) {
		int in[]={ipMH2s,ipMH2s},out[]={ipMH2s,ipMH,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	if( h2.lgH2ON  && hmi.lgBigH2_evaluated && hmi.lgH2_Chemistry_BigH2 )
	{
		hmi.h2sh2sh2s2h = hmi.Average_collH2s_dissoc;
	}

	r->rk = hmi.h2sh2sh2s2h;


	/* >>03 may 22, hmi.H2_Solomon_dissoc_rate_used is now rate electronic excited states
	 * decay into X continuum, not the rate that elec excited states are excited.
	 * assuming that 10% of excitations lead to ionization, this rate, for relax
	 * into ro-vib excited states, is 10x the Solomon rate */
	/* assume that 0.9 of H2 dissociations lead to H2_star,
	 * H2 + 0.9*hmi.H2_Solomon_dissoc_rate_used => h2_star */
	if(r->next == NULL) {
		int in[]={ipMH2g},out[]={ipMH2s};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	/* >>chng 03 may 22, see comment just above */
	/*r->rk = 0.9*hmi.H2_Solomon_dissoc_rate_used;
	r->rk = 10.*hmi.H2_Solomon_dissoc_rate_used;*/
	/* >>chng 03 sep 11, use real evaluated rate */
	r->rk = hmi.H2_H2g_to_H2s_rate_used;

	/* rate of photodissoc of vib-excit H2, A12 of TH85 */
	if(r->next == NULL) {
		int in[]={ipMH2s},out[]={ipMH,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = hmi.H2_photodissoc_used_H2s;


	/* >>chng 05 mar 24, TE, include continuum photodissociation from H2g*/
	if(r->next == NULL) {
		int in[]={ipMH2g},out[]={ipMH,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = hmi.H2_photodissoc_used_H2g;

	/* >>chng 05 oct 04, TE, rename and include in save H2 destruction
	 * >>chng 05 sept 30, GS, include collisional dissociation by electron, H2g + e = 2H + e*/
	if(r->next == NULL) {
		int in[]={ipMH2g},out[]={ipMH,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	hmi.h2ge2h = 1e-14*sexp(4.478*EVDEGK/phycon.te)*dense.eden;
	r->rk = hmi.h2ge2h;

	/* >>chng 05 oct 04, TE, rename and include in save H2 destruction
	 * >>chng 05 sept 30, GS, include collisional dissociation by electron, H2s + e = 2H + e*/	
	if(r->next == NULL) {
		int in[]={ipMH2s},out[]={ipMH,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	hmi.h2se2h = 1e-14*sexp(1.978*EVDEGK/phycon.te)*dense.eden;
	r->rk = hmi.h2se2h;

	/*---------------------------------------------------------------- */

	/* He H+ formation rates taken from Flower+Roueff, Black */

	/* He+ + H => HeH+
	 * radiative association from 
	 * >>refer	heh+	rate	Zygelman, B., and Dalgarno, A. 1990, ApJ 365, 239 */

	if(r->next == NULL) {
		int in[]={ipMH},out[]={ipMHeHp};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = 1e-15*dense.xIonDense[ipHELIUM][1];

	/* >>chng 03 sep 30 */
	/* back reaction of He+ + H=> He + H+	*/
	/* He + H+=> He+ + H;   
	 * bhephhphe = hephhphe*dense.xIonDense[ipHELIUM][1]*dense.xIonDense[ipHYDROGEN][0]/dense.xIonDense[ipHYDROGEN][1];*/
	/*bhephhphe = hephhphe*dense.xIonDense[ipHELIUM][1]*dense.xIonDense[ipHYDROGEN][0]/SDIV(dense.xIonDense[ipHYDROGEN][1]);
	if(r->next == NULL) {
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = bhephhphe;*/

	/* >>chng 03 sep 30 */
	/* He + H-=> He + H + e; 
	 * equation (table 5) from
	 * >>refer	H-	k	Paolo Lenzuni, David F. Chernoff, Edwin E. Salpeter, 1991, ApJS, 76, 759L  
	 * hehmeheh = 4.1e-17*pow(phycon.te,2)*sexp(19870/phycon.te)*dense.xIonDense[ipHELIUM][0];*/
	hehmeheh = 4.1e-17*phycon.tesqrd*sexp(19870/phycon.te)*dense.xIonDense[ipHELIUM][0];
	if(r->next == NULL) {
		int in[]={ipMHm},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =hehmeheh;

	/* >>chng 03 sep 30 */
	/* He+ + H2=> He + H+ + H; 
	 * equation (table 6) from
	 * >>refer	H2	dissoc	Tielens, A.G.G.M., & Hollenbach, D., 1985, ApJ, 291, 722  
	 * heph2hpheh = 1.5e-13*dense.xIonDense[ipHELIUM][1];*/
	/*hmi.rheph2hpheh = 1.5e-13f;*/
	/* >>chng 04 jun 10 */
	/* use the rate for this reaction as defined in UMIST */

	/* >> chng 05 aug 05, NPA comment.  This reaction can be an important He+ destruction term 
	   deep in molecular clouds.  This reaction is only slightly different from TH85, and has a temperature dependence.
	   Therefore this rate was switched to UMIST*/
	/** \todo	2	equivalent reaction for H2* is not included in chemistry, Big h2 does not include this reaction, what to do? GS */
	hmi.rheph2hpheh = (realnum)HMRATE(3.7e-14, 0., 35);
	if(r->next == NULL) {
		int in[]={ipMH2g},out[]={ipMHp,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = hmi.rheph2hpheh*dense.xIonDense[ipHELIUM][1];

	/* >>chng 04 jun 10 -- This reaction was not originally included.  It is in
	 * the UMIST database and seems to be an important He+ destruction mechanism */
	/* He+ + H2=> He + H2+  
	 * use the rate for this reaction as defined in UMIST */
	/** \todo	2	equivalent reaction for H2* is not included in chemistry, Big h2 does not include this reaction, what to do? GS */
	hmi.heph2heh2p = (realnum)HMRATE(7.2e-15, 0., 0);
	if(r->next == NULL) {
		int in[]={ipMH2g},out[]={ipMH2p};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = hmi.heph2heh2p*dense.xIonDense[ipHELIUM][1];

	/* He + H+ => HeH+ */	
	if(r->next == NULL) {
		int in[]={ipMHp},out[]={ipMHeHp};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = 1e-20*dense.xIonDense[ipHELIUM][0];

	/* H2+ + HE => HEH+ + H0 */
	if(r->next == NULL) {
		int in[]={ipMH2p},out[]={ipMH,ipMHeHp};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;

	/* >> chng 05 aug 05, NPA comment.  Turn off HeH+ for the Leiden comparison.  Can't
	   turn if off completely or else the matrix is unstable.  Just make sure that 
	   the abundance never affects the other molecules. */

	r->rk = 3e-10*exp(-6717./phycon.te)*dense.xIonDense[ipHELIUM][0]*co.lgUMISTrates;

	/* photodissociation through 1.6->2.3 continuum */

	/* why is this in a look instead of GammaK?
	 * to fix must set opacities into stack */
	gamheh = 0.;
	limit = MIN2(hmi.iheh2-1 , rfield.nflux );
	for( i=hmi.iheh1-1; i < limit; i++ )
	{
		gamheh += rfield.flux[0][i] + rfield.ConInterOut[i]+ rfield.outlin[0][i] + rfield.outlin_noplot[i];
	}
	gamheh *= 4e-18;

	/* hard radiation */
	gamheh += 3.*iso.gamnc[ipH_LIKE][ipHYDROGEN][ipH1s];

	/* recombination, HeH+  +  e => He + H */
	gamheh += dense.eden*1e-9;

	if(r->next == NULL) {
		int in[]={ipMHeHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = gamheh;

	/* HeH+  +  H => H2+  + He */
	if(r->next == NULL) {
		int in[]={ipMH,ipMHeHp},out[]={ipMH2p};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;

	/* >> chng 05 aug 05, NPA comment.  Turn off HeH+ for the Leiden comparison.  Can't
	   turn if off completely or else the matrix is unstable.  Just make sure that 
	   the abundance never affects the other molecules. */

	r->rk = 1e-10*co.lgUMISTrates;

	/* >>chng 03 sep 30 */
	/* HeH+  +  H2 => H3+  + He 
	 * equation (He13) from
	 * >>refer	HeH+	k	Galli, D., & Palla, F. 1998, A&A,335,403-420 
	 * hehph2h3phe = 1.3e-9;*/
	/*UMIST: 1.5e-9*/

	/* >> chng 05 aug 05, NPA comment.  Turn off HeH+ for the Leiden comparison.  Can't
	   turn if off completely or else the matrix is unstable.  Just make sure that 
	   the abundance never affects the other molecules. */
	/** \todo	2	equivalent reaction for H2* is not included in chemistry, Big h2 does not include this reaction, what to do? GS */

	hmi.hehph2h3phe = 1.3e-9*co.lgUMISTrates;
	if(r->next == NULL) {
		int in[]={ipMH2g,ipMHeHp},out[]={ipMH3p};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = hmi.hehph2h3phe;

	/* >>chng 03 sep 30 */
	/* He+  +  H- => H  + He 
	 * equation (20) from
	 * >>refer	H-	k	Stancil, P.C., Lepp, S., Dalgarno, A.1998, ApJ,509,1-10
	 * hephmhhe = 2.32e-7*pow(phycon.te/300,-.52)*sexp(phycon.te/-22400.)*dense.xIonDense[ipHELIUM][1];*/
	/* >>chng 03 oct 22, only add for Te < 1e5 - caused exception as T -> inf */
	hephmhhe = 2.32e-7*pow(phycon.te/300,-.52)*exp(TStancil/22400.)*dense.xIonDense[ipHELIUM][1];
	if(r->next == NULL) {
		int in[]={ipMHm},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =hephmhhe;

#	define CO_ON	true
	/*>>chng 06 jul 02, do not include cross terms with co molecules when advection on -
	 * this causes the H sum to be incorrect by as much as a few percent - no problem is
	 * co network if not on when dynamics is done */
	if( CO_ON && !dynamics.lgAdvection )
	{
	/* >>chng 04 may 26, NPA add this code */
	/* Many reactions in co.c also will create or destroy molecules that are predicted 
	   in hmole.  These reactions need to be inserted into this solver to obtain a more
	   accurate solution. These are the rates for reactions in co.c that create or destroy
	   species that are predicted in hmole_step.c.*/

	co.H_CH_C_H_H				= CO_findrk("H,CH=>C,H,H")*findspecies("CH")->hevmol;
	co.H_OH_O_H_H				= CO_findrk("H,OH=>O,H,H")*findspecies("OH")->hevmol;
	co.H_H2O_OH_H_H				= CO_findrk("H,H2O=>OH,H,H")*findspecies("H2O")->hevmol;
	co.H_COP_CO_HP				= CO_findrk("H,CO+=>CO,H+")*findspecies("CO+")->hevmol;
	co.H_CH_C_H2				= CO_findrk("H,CH=>C,H2")*findspecies("CH")->hevmol;
	co.H_CHP_CP_H2				= CO_findrk("H,CH+=>C+,H2")*findspecies("CH+")->hevmol;
	co.H_CH2_CH_H2				= CO_findrk("H,CH2=>CH,H2")*findspecies("CH2")->hevmol;
	co.H_CH3P_CH2P_H2			= CO_findrk("H,CH3+=>CH2+,H2")*findspecies("CH3+")->hevmol;
	co.H_OH_O_H2				= CO_findrk("H,OH=>O,H2")*findspecies("OH")->hevmol;
	co.H_H2O_OH_H2				= CO_findrk("H,H2O=>OH,H2")*findspecies("H2O")->hevmol;
	co.Hminus_HCOP_CO_H2	    = CO_findrk("H-,HCO+=>CO,H2")*findspecies("HCO+")->hevmol;
	co.Hminus_H3OP_H2O_H2	    = CO_findrk("H-,H3O+=>H2O,H2")*findspecies("H3O+")->hevmol;
	co.Hminus_H3OP_OH_H2_H	    = CO_findrk("H-,H3O+=>OH,H2,H")*findspecies("H3O+")->hevmol;
	co.HP_CH_CHP_H			    = CO_findrk("H+,CH=>CH+,H")*findspecies("CH")->hevmol;
	co.HP_CH2_CH2P_H		    = CO_findrk("H+,CH2=>CH2+,H")*findspecies("CH2")->hevmol;
	co.HP_H2O_H2OP_H		    = CO_findrk("H+,H2O=>H2O+,H")*findspecies("H2O")->hevmol;
	co.HP_O2_O2P_H			    = CO_findrk("H+,O2=>O2+,H")*findspecies("O2")->hevmol;
	co.HP_OH_OHP_H			    = CO_findrk("H+,OH=>OH+,H")*findspecies("OH")->hevmol;
	co.HP_SiO_SiOP_H		    = CO_findrk("H+,SiO=>SiO+,H")*findspecies("SiO")->hevmol;
	co.HP_CH2_CHP_H2		    = CO_findrk("H+,CH2=>CH+,H2")*findspecies("CH2")->hevmol;
	co.HP_SiH_SiP_H2		    = CO_findrk("H+,SiH=>Si+,H2")*findspecies("SiH")->hevmol;
	co.H2_CHP_CH2P_H		    = CO_findrk("H2,CH+=>CH2+,H")*findspecies("CH+")->hevmol;
	co.H2_CH2P_CH3P_H		    = CO_findrk("H2,CH2+=>CH3+,H")*findspecies("CH2+")->hevmol;
	co.H2_OHP_H2OP_H		    = CO_findrk("H2,OH+=>H2O+,H")*findspecies("OH+")->hevmol;
	co.H2_H2OP_H3OP_H		    = CO_findrk("H2,H2O+=>H3O+,H")*findspecies("H2O+")->hevmol;
	co.H2_COP_HCOP_H		    = CO_findrk("H2,CO+=>HCO+,H")*findspecies("CO+")->hevmol;
	co.H2_OP_OHP_H			    = CO_findrk("H2,O+=>OH+,H")*findspecies("O+")->hevmol;
	co.H2_SiOP_SiOHP_H		    = CO_findrk("H2,SiO+=>SiOH+,H")*findspecies("SiO+")->hevmol;
	co.H2_C_CH_H				= CO_findrk("H2,C=>CH,H")*findspecies("C")->hevmol;
	co.H2_CP_CHP_H				= CO_findrk("H2,C+=>CH+,H")*findspecies("C+")->hevmol;
	co.H2_CH_CH2_H				= CO_findrk("H2,CH=>CH2,H")*findspecies("CH")->hevmol;
	co.H2_OH_H2O_H				= CO_findrk("H2,OH=>H2O,H")*findspecies("OH")->hevmol;
	co.H2_O_OH_H				= CO_findrk("H2,O=>OH,H")*findspecies("O")->hevmol;
	co.H2_CH_C_H2_H				= CO_findrk("H2,CH=>C,H2,H")*findspecies("CH")->hevmol;
	co.H2_OH_O_H2_H				= CO_findrk("H2,OH=>O,H2,H")*findspecies("OH")->hevmol;
	co.H2_H2O_OH_H2_H			= CO_findrk("H2,H2O=>OH,H2,H")*findspecies("H2O")->hevmol;
	co.H2_O2_O_O_H2				= CO_findrk("H2,O2=>O,O,H2")*findspecies("O2")->hevmol;
	co.H2s_CH_C_H2_H			= CO_findrk("H2*,CH=>C,H2,H")*findspecies("CH")->hevmol*hmi.lgLeiden_Keep_ipMH2s;
	co.H2s_OH_O_H2_H			= CO_findrk("H2*,OH=>O,H2,H")*findspecies("OH")->hevmol*hmi.lgLeiden_Keep_ipMH2s;
	co.H2s_H2O_OH_H2_H			= CO_findrk("H2*,H2O=>OH,H2,H")*findspecies("H2O")->hevmol*hmi.lgLeiden_Keep_ipMH2s;
	co.H2s_O2_O_O_H2			= CO_findrk("H2*,O2=>O,O,H2")*findspecies("O2")->hevmol*hmi.lgLeiden_Keep_ipMH2s;
	co.H2P_C_CHP_H				= CO_findrk("H2+,C=>CH+,H")*findspecies("C")->hevmol;
	co.H2P_CH_CH2P_H			= CO_findrk("H2+,CH=>CH2+,H")*findspecies("CH")->hevmol;
	co.H2P_CH2_CH3P_H			= CO_findrk("H2+,CH2=>CH3+,H")*findspecies("CH2")->hevmol;
	co.H2P_OH_H2OP_H			= CO_findrk("H2+,OH=>H2O+,H")*findspecies("OH")->hevmol;
	co.H2P_H2O_H3OP_H			= CO_findrk("H2+,H2O=>H3O+,H")*findspecies("H2O")->hevmol;
	co.H2P_CO_HCOP_H			= CO_findrk("H2+,CO=>HCO+,H")*findspecies("CO")->hevmol;
	co.H2P_O_OHP_H				= CO_findrk("H2+,O=>OH+,H")*findspecies("O")->hevmol;
	co.H2P_CH_CHP_H2			= CO_findrk("H2+,CH=>CH+,H2")*findspecies("CH")->hevmol;
	co.H2P_CH2_CH2P_H2			= CO_findrk("H2+,CH2=>CH2+,H2")*findspecies("CH2")->hevmol;
	co.H2P_CO_COP_H2			= CO_findrk("H2+,CO=>CO+,H2")*findspecies("CO")->hevmol;
	co.H2P_H2O_H2OP_H2			= CO_findrk("H2+,H2O=>H2O+,H2")*findspecies("H2O")->hevmol;
	co.H2P_O2_O2P_H2			= CO_findrk("H2+,O2=>O2+,H2")*findspecies("O2")->hevmol;
	co.H2P_OH_OHP_H2			= CO_findrk("H2+,OH=>OH+,H2")*findspecies("OH")->hevmol;
	co.H3P_C_CHP_H2				= CO_findrk("H3+,C=>CH+,H2")*findspecies("C")->hevmol;
	co.H3P_CH_CH2P_H2			= CO_findrk("H3+,CH=>CH2+,H2")*findspecies("CH")->hevmol;
	co.H3P_CH2_CH3P_H2			= CO_findrk("H3+,CH2=>CH3+,H2")*findspecies("CH2")->hevmol;
	co.H3P_OH_H2OP_H2			= CO_findrk("H3+,OH=>H2O+,H2")*findspecies("OH")->hevmol;
	co.H3P_H2O_H3OP_H2			= CO_findrk("H3+,H2O=>H3O+,H2")*findspecies("H2O")->hevmol;
	co.H3P_CO_HCOP_H2			= CO_findrk("H3+,CO=>HCO+,H2")*findspecies("CO")->hevmol;
	co.H3P_O_OHP_H2				= CO_findrk("H3+,O=>OH+,H2")*findspecies("O")->hevmol;
	co.H3P_SiH_SiH2P_H2			= CO_findrk("H3+,SiH=>SiH2+,H2")*findspecies("SiH")->hevmol;
	co.H3P_SiO_SiOHP_H2			= CO_findrk("H3+,SiO=>SiOH+,H2")*findspecies("SiO")->hevmol;
	co.H2s_CH_CH2_H				= CO_findrk("H2*,CH=>CH2,H")*findspecies("CH")->hevmol*hmi.lgLeiden_Keep_ipMH2s;
	co.H2s_O_OH_H				= CO_findrk("H2*,O=>OH,H")*findspecies("O")->hevmol*hmi.lgLeiden_Keep_ipMH2s;
	co.H2s_OH_H2O_H				= CO_findrk("H2*,OH=>H2O,H")*findspecies("OH")->hevmol*hmi.lgLeiden_Keep_ipMH2s;
	co.H2s_C_CH_H				= CO_findrk("H2*,C=>CH,H")*findspecies("C")->hevmol*hmi.lgLeiden_Keep_ipMH2s;
	co.H2s_CP_CHP_H				= CO_findrk("H2*,C+=>CH+,H")*findspecies("C+")->hevmol*hmi.lgLeiden_Keep_ipMH2s;
	co.H_CH3_CH2_H2				= CO_findrk("H,CH3=>CH2,H2")*findspecies("CH3")->hevmol;
	co.H_CH4P_CH3P_H2			= CO_findrk("H,CH4+=>CH3+,H2")*findspecies("CH4+")->hevmol;
	co.H_CH5P_CH4P_H2			= CO_findrk("H,CH5+=>CH4+,H2")*findspecies("CH5+")->hevmol;
	co.H2_CH2_CH3_H				= CO_findrk("H2,CH2=>CH3,H")*findspecies("CH2")->hevmol;
	co.H2_CH3_CH4_H				= CO_findrk("H2,CH3=>CH4,H")*findspecies("CH3")->hevmol;
	co.H2_CH4P_CH5P_H			= CO_findrk("H2,CH4+=>CH5+,H")*findspecies("CH4+")->hevmol;
	co.H2s_CH2_CH3_H			= CO_findrk("H2*,CH2=>CH3,H")*findspecies("CH2")->hevmol*hmi.lgLeiden_Keep_ipMH2s;
	co.H2s_CH3_CH4_H			= CO_findrk("H2*,CH3=>CH4,H")*findspecies("CH3")->hevmol*hmi.lgLeiden_Keep_ipMH2s;
	co.H2P_CH4_CH3P_H2			= CO_findrk("H2+,CH4=>CH3+,H2,H")*findspecies("CH4")->hevmol;
	co.H2P_CH4_CH4P_H2			= CO_findrk("H2+,CH4=>CH4+,H2")*findspecies("CH4")->hevmol;
	co.H2P_CH4_CH5P_H			= CO_findrk("H2+,CH4=>CH5+,H")*findspecies("CH4")->hevmol;
	/* >>chng 05 jul 15, TE, change from findspecies("CH")->hevmol to findspecies("CH3")->hevmol */
	co.H3P_CH3_CH4P_H2			= CO_findrk("H3+,CH3=>CH4+,H2")*findspecies("CH3")->hevmol;
	co.H3P_CH4_CH5P_H2			= CO_findrk("H3+,CH4=>CH5+,H2")*findspecies("CH4")->hevmol;
	co.HP_CH3_CH3P_H			= CO_findrk("H+,CH3=>CH3+,H")*findspecies("CH3")->hevmol;
	co.HP_CH4_CH3P_H2			= CO_findrk("H+,CH4=>CH3+,H2")*findspecies("CH4")->hevmol;
	co.HP_CH4_CH4P_H			= CO_findrk("H+,CH4=>CH4+,H")*findspecies("CH4")->hevmol;


	/*  >>chng 04 jul 14 -- NPA. This is the contribution that the heavy element 
	 *  network makes to the hydrogen chemistry due to reactions with N and S. */

	co.H2_N_NH_H           = CO_findrk("H2,N=>NH,H")*findspecies("N")->hevmol;
	co.H2_NH_NH2_H         = CO_findrk("H2,NH=>NH2,H")*findspecies("NH")->hevmol;
	co.H2_NH2_NH3_H        = CO_findrk("H2,NH2=>NH3,H")*findspecies("NH2")->hevmol;
	co.H2_CN_HCN_H         = CO_findrk("H2,CN=>HCN,H")*findspecies("CN")->hevmol;
	co.HP_HNO_NOP_H2       = CO_findrk("H+,HNO=>NO+,H2")*findspecies("HNO")->hevmol;
	co.HP_HS_SP_H2         = CO_findrk("H+,HS=>S+,H2")*findspecies("HS")->hevmol;
	co.H_HSP_SP_H2         = CO_findrk("H,HS+=>S+,H2")*findspecies("HS+")->hevmol;
	co.H2P_N_NHP_H         = CO_findrk("H2+,N=>NH+,H")*findspecies("N")->hevmol;
	co.H2_NP_NHP_H         = CO_findrk("H2,N+=>NH+,H")*findspecies("N+")->hevmol;
	co.H2_NHP_N_H3P        = CO_findrk("H2,NH+=>N,H3+")*findspecies("NH+")->hevmol;
	co.H2P_NH_NH2P_H       = CO_findrk("H2+,NH=>NH2+,H")*findspecies("NH")->hevmol;
	co.H2_NHP_NH2P_H       = CO_findrk("H2,NH+=>NH2+,H")*findspecies("NH+")->hevmol;
	co.H2_NH2P_NH3P_H      = CO_findrk("H2,NH2+=>NH3+,H")*findspecies("NH2+")->hevmol;
	co.H2_NH3P_NH4P_H      = CO_findrk("H2,NH3+=>NH4+,H")*findspecies("NH3+")->hevmol;
	co.H2P_CN_HCNP_H       = CO_findrk("H2+,CN=>HCN+,H")*findspecies("CN")->hevmol;
	co.H2_CNP_HCNP_H       = CO_findrk("H2,CN+=>HCN+,H")*findspecies("CN+")->hevmol;
	co.H2P_NO_HNOP_H       = CO_findrk("H2+,NO=>HNO+,H")*findspecies("NO")->hevmol;
	co.H2_SP_HSP_H         = CO_findrk("H2,S+=>HS+,H")*findspecies("S+")->hevmol;
	co.H2_CSP_HCSP_H       = CO_findrk("H2,CS+=>HCS+,H")*findspecies("CS+")->hevmol;
	co.H3P_NH_NH2P_H2      = CO_findrk("H3+,NH=>NH2+,H2")*findspecies("NH")->hevmol;
	co.H3P_NH2_NH3P_H2     = CO_findrk("H3+,NH2=>NH3+,H2")*findspecies("NH2")->hevmol;
	co.H3P_NH3_NH4P_H2     = CO_findrk("H3+,NH3=>NH4+,H2")*findspecies("NH3")->hevmol;
	co.H3P_CN_HCNP_H2      = CO_findrk("H3+,CN=>HCN+,H2")*findspecies("CN")->hevmol;
	co.H3P_NO_HNOP_H2      = CO_findrk("H3+,NO=>HNO+,H2")*findspecies("NO")->hevmol;
	co.H3P_S_HSP_H2        = CO_findrk("H3+,S=>HS+,H2")*findspecies("S")->hevmol;
	co.H3P_CS_HCSP_H2      = CO_findrk("H3+,CS=>HCS+,H2")*findspecies("CS")->hevmol;
	co.H3P_NO2_NOP_OH_H2   = CO_findrk("H3+,NO2=>NO+,OH,H2")*findspecies("NO2")->hevmol;
	co.HP_NH_NHP_H         = CO_findrk("H+,NH=>NH+,H")*findspecies("NH")->hevmol;
	co.HP_NH2_NH2P_H       = CO_findrk("H+,NH2=>NH2+,H")*findspecies("NH2")->hevmol;
	co.HP_NH3_NH3P_H       = CO_findrk("H+,NH3=>NH3+,H")*findspecies("NH3")->hevmol;
	co.H_CNP_CN_HP         = CO_findrk("H,CN+=>CN,H+")*findspecies("CN+")->hevmol;
	co.HP_HCN_HCNP_H       = CO_findrk("H+,HCN=>HCN+,H")*findspecies("HCN")->hevmol;
	co.H_HCNP_HCN_HP       = CO_findrk("H,HCN+=>HCN,H+")*findspecies("HCN+")->hevmol;
	co.H_N2P_N2_HP         = CO_findrk("H,N2+=>N2,H+")*findspecies("N2+")->hevmol;
	co.HP_NO_NOP_H         = CO_findrk("H+,NO=>NO+,H")*findspecies("NO")->hevmol;
	co.HP_HS_HSP_H         = CO_findrk("H+,HS=>HS+,H")*findspecies("HS")->hevmol;
	co.HP_SiN_SiNP_H       = CO_findrk("H+,SiN=>SiN+,H")*findspecies("SiN")->hevmol;
	co.HP_CS_CSP_H         = CO_findrk("H+,CS=>CS+,H")*findspecies("CS")->hevmol;
	co.HP_NS_NSP_H         = CO_findrk("H+,NS=>NS+,H")*findspecies("NS")->hevmol;
	co.HP_SO_SOP_H         = CO_findrk("H+,SO=>SO+,H")*findspecies("SO")->hevmol;
	co.HP_OCS_OCSP_H       = CO_findrk("H+,OCS=>OCS+,H")*findspecies("OCS")->hevmol;
	co.HP_S2_S2P_H         = CO_findrk("H+,S2=>S2+,H")*findspecies("S2")->hevmol;
	co.H2P_NH_NHP_H2       = CO_findrk("H2+,NH=>NH+,H2")*findspecies("NH")->hevmol;
	co.H2P_NH2_NH2P_H2     = CO_findrk("H2+,NH2=>NH2+,H2")*findspecies("NH2")->hevmol;
	co.H2P_NH3_NH3P_H2     = CO_findrk("H2+,NH3=>NH3+,H2")*findspecies("NH3")->hevmol;
	co.H2P_CN_CNP_H2       = CO_findrk("H2+,CN=>CN+,H2")*findspecies("CN")->hevmol;
	co.H2P_HCN_HCNP_H2     = CO_findrk("H2+,HCN=>HCN+,H2")*findspecies("HCN")->hevmol;
	co.H2P_NO_NOP_H2       = CO_findrk("H2+,NO=>NO+,H2")*findspecies("NO")->hevmol;
	/*>>chng 05 jul 22, TE, included reaction which contribute to the molecular hydrogen network */
	/* that were in mole_co_step, but not in mole_h_step*/
	co.HP_C2_C2P_H		   = CO_findrk("H+,C2=>C2+,H")*findspecies("C2")->hevmol;
	co.H2P_C2_C2P_H2	   = CO_findrk("H2+,C2=>C2+,H2")*findspecies("C2")->hevmol;	
	co.Hminus_NH4P_NH3_H2  = CO_findrk("H-,NH4+=>NH3,H2")*findspecies("NH4+")->hevmol;
	co.Hminus_NP_N_H       = CO_findrk("H-,N+=>N,H")*findspecies("N+")->hevmol;

	/* >> chng 05 jul 15, TE, make the chlorine reactions a member of mole.h */
	/* >>chng 05 mar 23 -- NPA. Add Chlorine chemical reactions to hydrogen ionization balance */
	/* >>chng 05 jul 21 -- NPA. TE noticed that two reactions that contribute to the molecular hydrogen network,
	 * and in mole_co_step.c, were not included in h_step.  These reactions are now inserted.*/
	/* >>chng 05 aug 1 -- NPA. Gargi Shaw noticed that the reaction H2 + S => HS + H was not in the CO network
	 * This has been changed.  This reaction also affects H2 network so must go here also .*/

	co.H2_ClP_HClP_H = CO_findrk("H2,Cl+=>HCl+,H")*findspecies("Cl+")->hevmol;
	co.H2_HClP_H2ClP_H = CO_findrk("H2,HCl+=>H2Cl+,H")*findspecies("HCl+")->hevmol;
	co.H3P_Cl_HClP_H2 = CO_findrk("H3+,Cl=>HCl+,H2")*findspecies("Cl")->hevmol;
	co.H3P_HCl_H2ClP_H2 = CO_findrk("H3+,HCl=>H2Cl+,H2")*findspecies("HCl")->hevmol;
	co.HP_HCl_HClP_H = CO_findrk("H+,HCl=>HCl+,H")*findspecies("HCl")->hevmol;
	/* >>chng 06 feb 06 rjrw -- was multiplied by HS not S */	
	co.H2_S_HS_H = CO_findrk("H2,S=>HS,H")*findspecies("S")->hevmol;
	/* >>chng 05 sept 14 - NPA, add HNC and HCNH+ to hmole_step */
	co.HP_HNC_HCN_HP = CO_findrk("H+,HNC=>HCN,H+")*findspecies("HNC")->hevmol;
	co.H_HNC_HCN_H = CO_findrk("H,HNC=>HCN,H")*findspecies("HNC")->hevmol;
	/* >>chng 06 feb 06 rjrw -- was duplicated */
	/* co.H_HNC_HCN_H = CO_findrk("H,HNC=>HCN,H")*findspecies("HNC")->hevmol; */
	co.H2_HCNP_HCNHP_H = CO_findrk("H2,HCN+=>HCNH+,H")*findspecies("HCN+")->hevmol;
	/* >>chng 06 Sep 03 rjrw HNC -> HCN in COmole index for consistency */
	co.H3P_HCN_HCNHP_H2 = CO_findrk("H3+,HCN=>HCNH+,H2")*findspecies("HCN")->hevmol;
	/* >>chng 05 nov 17, TE added following reaction, the same rate as H2g is assumed*/
	co.H2s_OP_OHP_H	= CO_findrk("H2*,O+=>OH+,H")*findspecies("O+")->hevmol;



	/* >>chng 05 dec 01 - NPA.  Add reactions of C3, C3+, C2H, C2H+, C3H, C3H+, C2H2,
	 * C2H2+, and C2H3+ with hydrogen molecules to the overall balance */
	co.HP_C2H2_C2H2P_H = CO_findrk("H+,C2H2=>C2H2+,H")*findspecies("C2H2")->hevmol;
	co.HP_C2H2_C2HP_H2 = CO_findrk("H+,C2H2=>C2H+,H2")*findspecies("C2H2")->hevmol;
	co.HP_C3H_C3HP_H = CO_findrk("H+,C3H=>C3H+,H")*findspecies("C3H")->hevmol;
	co.HP_C3H_C3P_H2 = CO_findrk("H+,C3H=>C3+,H2")*findspecies("C3H")->hevmol;
	co.H2P_C2H_C2H2P_H = CO_findrk("H2+,C2H=>C2H2+,H")*findspecies("C2H")->hevmol;
	co.H2P_C2H2_C2H2P_H2 = CO_findrk("H2+,C2H2=>C2H2+,H2")*findspecies("C2H2")->hevmol;
	co.H3P_C2H_C2H2P_H2 = CO_findrk("H3+,C2H=>C2H2+,H2")*findspecies("C2H")->hevmol;
	co.H3P_C3_C3HP_H2 = CO_findrk("H3+,C3=>C3H+,H2")*findspecies("C3")->hevmol;
	co.H2_C2HP_C2H2P_H = CO_findrk("H2,C2H+=>C2H2+,H")*findspecies("C2H+")->hevmol;
	co.H2_C3P_C3HP_H = CO_findrk("H2,C3+=>C3H+,H")*findspecies("C3+")->hevmol;
	co.H_C2H3P_C2H2P_H2 = CO_findrk("H,C2H3+=>C2H2+,H2")*findspecies("C2H3+")->hevmol;
	co.H3P_C2H2_C2H3P_H2 = CO_findrk("H3+,C2H2=>C2H3+,H2")*findspecies("C2H2")->hevmol;
	co.H2P_C2H2_C2H3P_H = CO_findrk("H2+,C2H2=>C2H3+,H")*findspecies("C2H2")->hevmol;
	co.HP_C3_C3P_H = CO_findrk("H+,C3=>C3+,H")*findspecies("C3")->hevmol;
	co.HP_C2H_C2HP_H = CO_findrk("H+,C2H=>C2H+,H")*findspecies("C2H")->hevmol;
	co.H2P_C2_C2HP_H = CO_findrk("H2+,C2=>C2H+,H")*findspecies("C2")->hevmol;
	co.H2P_C2H_C2HP_H2 = CO_findrk("H2+,C2H=>C2H+,H2")*findspecies("C2H")->hevmol;
	co.H3P_C2_C2HP_H2 = CO_findrk("H3+,C2=>C2H+,H2")*findspecies("C2")->hevmol;
	co.H2_C2P_C2HP_H = CO_findrk("H2,C2+=>C2H+,H")*findspecies("C2+")->hevmol;
	co.HP_C2H_C2P_H2 = CO_findrk("H+,C2H=>C2+,H2")*findspecies("C2H")->hevmol;
	/* >>chng 13 Apr 2006, Add N2H+ to chemistry, which 
	 * should improve modeling of nitrogen chemistry, in 
	 * particular NH and N2 */
	co.N2_H3P_N2HP_H2 = CO_findrk("N2,H3+=>N2H+,H2")*findspecies("N2")->hevmol;
	if(r->next == NULL) 
	{
		int in[]={ipMH},out[]={ipMH,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = co.H_CH_C_H_H;

	if(r->next == NULL)
	{
		int in[]={ipMH},out[]={ipMH,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H_OH_O_H_H;

	if(r->next == NULL)
	{
		int in[]={ipMH},out[]={ipMH,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H_H2O_OH_H_H;

	if(r->next == NULL)
	{
		int in[]={ipMH},out[]={ipMHp};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H_COP_CO_HP;

	if(r->next == NULL)
	{
		int in[]={ipMH},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H_CH_C_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H_CHP_CP_H2;

	if(r->next == NULL) 
	{
		int in[]={ipMH},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H_CH2_CH_H2;

	if(r->next == NULL) 
	{
		int in[]={ipMH},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H_CH3P_CH2P_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H_OH_O_H2;

	if(r->next == NULL) 
	{
		int in[]={ipMH},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H_H2O_OH_H2;


	if(r->next == NULL) 
	{
		int in[]={ipMHm},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.Hminus_HCOP_CO_H2;

	if(r->next == NULL) 
	{
		int in[]={ipMHm},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.Hminus_H3OP_H2O_H2;

	if(r->next == NULL)
	{
		int in[]={ipMHm},out[]={ipMH2g, ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.Hminus_H3OP_OH_H2_H;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.HP_CH_CHP_H;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.HP_CH2_CH2P_H;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.HP_H2O_H2OP_H;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.HP_O2_O2P_H;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.HP_OH_OHP_H;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.HP_SiO_SiOP_H;

	if(r->next == NULL) 
	{
		int in[]={ipMHp},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.HP_CH2_CHP_H2;

	if(r->next == NULL) 
	{
		int in[]={ipMHp},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.HP_SiH_SiP_H2;

	if(r->next == NULL) 
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2_CHP_CH2P_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2_CH2P_CH3P_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2_OHP_H2OP_H;

	if(r->next == NULL) 
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2_H2OP_H3OP_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2_COP_HCOP_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2_OP_OHP_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2_SiOP_SiOHP_H;

	if(r->next == NULL) 
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2_C_CH_H;

	if(r->next == NULL) 
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2_CP_CHP_H;

	if(r->next == NULL) 
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2_CH_CH2_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2_OH_H2O_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2_O_OH_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH2g,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2_CH_C_H2_H;

	if(r->next == NULL) 
	{
		int in[]={ipMH2g},out[]={ipMH2g, ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2_OH_O_H2_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH2g, ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2_H2O_OH_H2_H;

	if(r->next == NULL) 
	{
		int in[]={ipMH2g},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2_O2_O_O_H2;

	/* >>chng 05 mar 18, by TE, this was in twice, had little effect since high
	 * temperature barrier */
	/*if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2_C_CH_H;*/

	if(r->next == NULL) 
	{
		int in[]={ipMH2s},out[]={ipMH2g,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2s_CH_C_H2_H;

	if(r->next == NULL) 
	{
		int in[]={ipMH2s},out[]={ipMH2g,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2s_OH_O_H2_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2s},out[]={ipMH2g,ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2s_H2O_OH_H2_H;

	if(r->next == NULL) 
	{
		int in[]={ipMH2s},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2s_O2_O_O_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH2p},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2P_C_CHP_H;

	if(r->next == NULL) 
	{
		int in[]={ipMH2p},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2P_CH_CH2P_H;

	if(r->next == NULL) {
		int in[]={ipMH2p},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2P_CH2_CH3P_H;

	if(r->next == NULL) {
		int in[]={ipMH2p},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2P_OH_H2OP_H;

	if(r->next == NULL) {
		int in[]={ipMH2p},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2P_H2O_H3OP_H;

	if(r->next == NULL) {
		int in[]={ipMH2p},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2P_CO_HCOP_H;

	if(r->next == NULL) {
		int in[]={ipMH2p},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2P_O_OHP_H;

	if(r->next == NULL) {
		int in[]={ipMH2p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2P_CH_CHP_H2;

	if(r->next == NULL) {
		int in[]={ipMH2p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2P_CH2_CH2P_H2;

	if(r->next == NULL) {
		int in[]={ipMH2p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2P_CO_COP_H2;

	if(r->next == NULL) {
		int in[]={ipMH2p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2P_H2O_H2OP_H2;

	if(r->next == NULL) {
		int in[]={ipMH2p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2P_O2_O2P_H2;

	if(r->next == NULL) {
		int in[]={ipMH2p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2P_OH_OHP_H2;

	if(r->next == NULL) {
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H3P_C_CHP_H2;

	if(r->next == NULL) {
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H3P_CH_CH2P_H2;

	if(r->next == NULL) {
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H3P_CH2_CH3P_H2;

	if(r->next == NULL) {
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H3P_OH_H2OP_H2;

	if(r->next == NULL) {
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H3P_H2O_H3OP_H2;

	if(r->next == NULL) {
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H3P_CO_HCOP_H2;

	if(r->next == NULL) 
	{
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H3P_O_OHP_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H3P_SiH_SiH2P_H2;

	if(r->next == NULL) {
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H3P_SiO_SiOHP_H2;

	if(r->next == NULL) 
	{
		int in[]={ipMH2s},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2s_CH_CH2_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2s},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2s_O_OH_H;

	if(r->next == NULL) 
	{
		int in[]={ipMH2s},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2s_OH_H2O_H;


	if(r->next == NULL) 
	{
		int in[]={ipMH2s},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2s_C_CH_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2s},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2s_CP_CHP_H;

	if(r->next == NULL) 
	{
		int in[]={ipMH},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H_CH3_CH2_H2;

	if(r->next == NULL) 
	{
		int in[]={ipMH},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H_CH4P_CH3P_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H_CH5P_CH4P_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2_CH2_CH3_H;


	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2_CH3_CH4_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2_CH4P_CH5P_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2s},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2s_CH2_CH3_H;

	if(r->next == NULL) 
	{
		int in[]={ipMH2s},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2s_CH3_CH4_H;

	if(r->next == NULL) 
	{
		int in[]={ipMH2p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2P_CH4_CH3P_H2;

	if(r->next == NULL) 
	{
		int in[]={ipMH2p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2P_CH4_CH4P_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH2p},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2P_CH4_CH5P_H;

	if(r->next == NULL)
	{
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H3P_CH3_CH4P_H2;

	if(r->next == NULL) 
	{
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H3P_CH4_CH5P_H2;

	if(r->next == NULL) {
		int in[]={ipMHp},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.HP_CH4_CH3P_H2;

	if(r->next == NULL) {
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.HP_CH4_CH4P_H;

	if(r->next == NULL) {
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2_ClP_HClP_H;

	if(r->next == NULL) {
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2_HClP_H2ClP_H;

	if(r->next == NULL) {
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H3P_Cl_HClP_H2;

	if(r->next == NULL) {
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H3P_HCl_H2ClP_H2;

	if(r->next == NULL) {
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.HP_HCl_HClP_H;
	/* >>chng 05 aug 02, NA add following reaction */

	if(r->next == NULL) {
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2_S_HS_H;


	/* The following terms are for the reactions with the rates defined above.  The 
	 * actual reaction is self-contained in the rate label. */

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2_N_NH_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2_NH_NH2_H;


	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2_NH2_NH3_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2_CN_HCN_H;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.HP_HNO_NOP_H2;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.HP_HS_SP_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H_HSP_SP_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH2p},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2P_N_NHP_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2_NP_NHP_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH3p};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2_NHP_N_H3P;

	if(r->next == NULL)
	{
		int in[]={ipMH2p},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2P_NH_NH2P_H;


	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2_NHP_NH2P_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2_NH2P_NH3P_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2_NH3P_NH4P_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2p},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2P_CN_HCNP_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2_CNP_HCNP_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2p},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2P_NO_HNOP_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2_SP_HSP_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2_CSP_HCSP_H;


	if(r->next == NULL)
	{
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H3P_NH_NH2P_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H3P_NH2_NH3P_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H3P_NH3_NH4P_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H3P_CN_HCNP_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H3P_NO_HNOP_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H3P_S_HSP_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H3P_CS_HCSP_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H3P_NO2_NOP_OH_H2;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.HP_NH_NHP_H;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.HP_NH2_NH2P_H;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.HP_NH3_NH3P_H;


	if(r->next == NULL)
	{
		int in[]={ipMH},out[]={ipMHp};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H_CNP_CN_HP;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.HP_HCN_HCNP_H;

	if(r->next == NULL)
	{
		int in[]={ipMH},out[]={ipMHp};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H_HCNP_HCN_HP;

	if(r->next == NULL)
	{
		int in[]={ipMH},out[]={ipMHp};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H_N2P_N2_HP;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.HP_NO_NOP_H;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.HP_HS_HSP_H;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.HP_SiN_SiNP_H;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.HP_CS_CSP_H;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.HP_NS_NSP_H;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.HP_SO_SOP_H;


	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.HP_OCS_OCSP_H;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.HP_S2_S2P_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2P_NH_NHP_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH2p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2P_NH2_NH2P_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH2p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2P_NH3_NH3P_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH2p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2P_CN_CNP_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH2p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2P_HCN_HCNP_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH2p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2P_NO_NOP_H2;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.HP_C2_C2P_H;


	if(r->next == NULL)
	{
		int in[]={ipMH2p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2P_C2_C2P_H2;

	if(r->next == NULL)
	{
		int in[]={ipMHm},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.Hminus_NH4P_NH3_H2;

	if(r->next == NULL)
	{
		int in[]={ipMHm},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.Hminus_NP_N_H;

	/* >>chng 05 sept 14 - NPA.  Reactions involving HNC or HCNH+ that affect hydrogen chemistry */

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMHp};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.HP_HNC_HCN_HP;

	if(r->next == NULL)
	{
		int in[]={ipMH},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H_HNC_HCN_H;

	/* >>chng 06 feb 06 rjrw -- was duplicated */
	/* if(r->next == NULL)
	{
		int in[]={ipMH},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H_HNC_HCN_H; */

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H2_HCNP_HCNHP_H;

	if(r->next == NULL)
	{
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =	co.H3P_HCN_HCNHP_H2;

	/* >>chng 05 nov 17, TE  added following reaction*/
	if(r->next == NULL)
	{
		int in[]={ipMH2s},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.H2s_OP_OHP_H;

	/* >>chng 05 dec 01 - NPA.  Add reactions for some more complex molecules
	 * into the hydrogen chemistry */
	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = co.HP_C2H2_C2H2P_H;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = co.HP_C2H2_C2HP_H2;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = co.HP_C3H_C3HP_H;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk =co.HP_C3H_C3P_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH2p},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = co.H2P_C2H_C2H2P_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = co.H2P_C2H2_C2H2P_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = co.H3P_C2H_C2H2P_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = co.H3P_C3_C3HP_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = co.H2_C2HP_C2H2P_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = co.H2_C3P_C3HP_H;

	if(r->next == NULL)
	{
		int in[]={ipMH},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = co.H_C2H3P_C2H2P_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = co.H3P_C2H2_C2H3P_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH2p},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = co.H2P_C2H2_C2H3P_H;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = co.HP_C3_C3P_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2p},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = co.H2P_C2_C2HP_H;

	if(r->next == NULL)
	{
		int in[]={ipMH2p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = co.H2P_C2H_C2HP_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = co.H3P_C2_C2HP_H2;

	if(r->next == NULL)
	{
		int in[]={ipMH2g},out[]={ipMH};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = co.H2_C2P_C2HP_H;

	if(r->next == NULL)
	{
		int in[]={ipMHp},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = co.HP_C2H_C2P_H2;


	/* >>chng 13 Apr 2006, Add N2H+ to chemistry, which 
	 * should improve modeling of nitrogen chemistry, in 
	 * particular NH and N2 */

	if(r->next == NULL)
	{
		int in[]={ipMH3p},out[]={ipMH2g};
		r->next = newreaction(rindex,in,INTSZ(in),out,INTSZ(out),NULL,0);
	}
	r = r->next;
	rindex++;
	r->rk = co.N2_H3P_N2HP_H2;
	}
	else
	{
		fixit();/* rm macro set above to hose co chem */
	}

	/* >>chng 06 jun 30, these are not used anywhere else and should be removed from co header
	 * this couples H and CO chem but was unstable.  ignore their contribution to H
	 * chemistry for now - OK because CO moles have little effect on H solution? */
#	if 0
	/*  The following reactions in co.c include a formation or destruction
	    process, but not both, for one of the molecules in mole_h_step.c.
		Because of this, the following rates have to be put into the bvec[i]
		part of the matrix. */

	/* Zero out the CO contribution to the sources and sinks in hmole */


	for(i=0; i<N_H_MOLEC; ++i)
	{
		co.hydro_source[i] = 0;
		co.hydro_sink[i] = 0;
	}

	/* C + H3OP =  HCOP + H2 */   

	co.C_H3OP_HCOP_H2_1 =  HMRATE(1.0e-11,0,0)*findspecies("H3O+")->hevmol*findspecies("C")->hevmol;
	co.hydro_source[ipMH2g] += (realnum)co.C_H3OP_HCOP_H2_1;


	/* C + OH =  CO + H  */  

	co.C_OH_CO_H_1 = HMRATE(1.1e-10,0.5,0)*findspecies("OH")->hevmol*findspecies("C")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.C_OH_CO_H_1;

	/* CP + OH =>  CO + HP */   

	co.CP_OH_CO_HP_1 = HMRATE(7.7e-10,0,0)*findspecies("OH")->hevmol*findspecies("C+")->hevmol;
	co.hydro_source[ipMHp] += (realnum)co.CP_OH_CO_HP_1;

	/* CP + H2O =>  HCOP + H  */  

	co.CP_H2O_HCOP_H_1 = HMRATE(9.0e-10,0,0)*findspecies("H2O")->hevmol*findspecies("C+")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.CP_H2O_HCOP_H_1;

	/* CP + OH =>  COP + H  */  

	co.CP_OH_COP_H_1 = HMRATE(7.7e-10,0,0)*findspecies("OH")->hevmol*findspecies("C+")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.CP_OH_COP_H_1;

	/* O + CH =>  CO + H  */  

	co.O_CH_CO_H_1 = HMRATE(6.6e-11,0,0)*findspecies("CH")->hevmol*findspecies("O")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.O_CH_CO_H_1;

	/* O + CHP =>  COP + H  */  

	co.O_CHP_COP_H_1 = HMRATE(3.5e-10,0,0)*findspecies("CH+")->hevmol*findspecies("O")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.O_CHP_COP_H_1;

	/* O + CH2 =>  CO + H + H  */  

	co.O_CH2_CO_H_H_1 = HMRATE(1.33e-10,0,0)*findspecies("CH2")->hevmol*findspecies("O")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.O_CH2_CO_H_H_1;

	/* O + CH2 =>  CO + H2 */  

	co.O_CH2_CO_H2_1 = HMRATE(8.0e-11,0,0)*findspecies("CH2")->hevmol*findspecies("O")->hevmol;
	co.hydro_source[ipMH2g] += (realnum)co.O_CH2_CO_H2_1;

	/* O + CH2P =>  HCOP + H  */  

	co.O_CH2P_HCOP_H_1 = HMRATE(7.5e-10,0,0)*findspecies("CH2+")->hevmol*findspecies("O")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.O_CH2P_HCOP_H_1;

	/* O + CH3P =>  HCOP + H2  */   

	co.O_CH3P_HCOP_H2_1	= HMRATE(4.0e-10,0,0)*findspecies("CH3+")->hevmol*findspecies("O")->hevmol;
	co.hydro_source[ipMH2g] += (realnum)co.O_CH3P_HCOP_H2_1;

	/* O + H2OP =>  O2P + H2  */   

	co.O_H2OP_O2P_H2_1 = HMRATE(4.0e-11,0,0)*findspecies("H2O+")->hevmol*findspecies("O")->hevmol;
	co.hydro_source[ipMH2g] += (realnum)co.O_H2OP_O2P_H2_1;

	/* O + OH =>  O2 + H  */  

	co.O_OH_O2_H_1 = HMRATE(4.34e-11,-0.5,30)*findspecies("OH")->hevmol*findspecies("O")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.O_OH_O2_H_1;

	/* O + OHP =>  O2P + H  */  

	co.O_OHP_O2P_H_1 = HMRATE(7.1e-10,0,0)*findspecies("OH+")->hevmol*findspecies("O")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.O_OHP_O2P_H_1;

	/* O + SiH =>  SiO + H  */  

	co.O_SiH_SiO_H_1 = HMRATE(4.0e-11,0.5,0)*findspecies("SiH")->hevmol*findspecies("O")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.O_SiH_SiO_H_1;

	/* O + SiH2P =>  SiOHP + H  */  

	co.O_SiH2P_SiOHP_H_1 = HMRATE(6.3e-10,0,0)*findspecies("SiH2+")->hevmol*findspecies("O")->hevmol;	
	co.hydro_source[ipMH] += (realnum)co.O_SiH2P_SiOHP_H_1;

	/* OP + CH =>  COP + H  */  

	co.OP_CH_COP_H_1 = HMRATE(3.5e-10,0,0)*findspecies("CH")->hevmol*findspecies("O+")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.OP_CH_COP_H_1;

	/* OP + OH =>  O2P + H  */  

	co.OP_OH_O2P_H_1 = HMRATE(3.6e-10,0,0)*findspecies("OH")->hevmol*findspecies("O+")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.OP_OH_O2P_H_1;

	/* Si + OH =>  SiO + H  */  

	co.Si_OH_SiO_H_1 = HMRATE(2.0e-10,0.5,0)*findspecies("OH")->hevmol*findspecies("Si")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.Si_OH_SiO_H_1;

	/* SiP + H2O =>  SiOHP + H  */  

	co.SiP_H2O_SiOHP_H_1 = HMRATE(2.3e-10,0,0)*findspecies("H2O")->hevmol*findspecies("Si+")->hevmol;			
	co.hydro_source[ipMH] += (realnum)co.SiP_H2O_SiOHP_H_1;

	/* SiP + OH =>  SiOP + H  */  

	co.SiP_OH_SiOP_H_1 = HMRATE(6.3e-10,0,0)*findspecies("OH")->hevmol*findspecies("Si+")->hevmol;		
	co.hydro_source[ipMH] += (realnum)co.SiP_OH_SiOP_H_1;

	/* CHP + H2O =>  HCOP +  H2  */   

	co.CHP_H2O_HCOP_H2_1 = HMRATE(2.9e-9,0,0)*findspecies("H2O")->hevmol*findspecies("CH+")->hevmol;			
	co.hydro_source[ipMH2g] += (realnum)co.CHP_H2O_HCOP_H2_1;

	/* CHP + OH =>  COP + H2  */   

	co.CHP_OH_COP_H2_1 = HMRATE(7.5e-10,0,0)*findspecies("OH")->hevmol*findspecies("CH+")->hevmol;		
	co.hydro_source[ipMH2g] += (realnum)co.CHP_OH_COP_H2_1;

	/* H + C =>  CH + nu  */

	co.H_C_CH_nu = HMRATE(0.00000000000000001,0,0)*dense.xIonDense[ipHYDROGEN][0]*findspecies("C")->hevmol;
	co.hydro_sink[ipMH] += (realnum)co.H_C_CH_nu;

	/* H + CP =>  CHP + nu  */

	co.H_CP_CHP_nu = HMRATE(1.7e-17,0,0)*findspecies("C+")->hevmol*dense.xIonDense[ipHYDROGEN][0];	
	co.hydro_sink[ipMH] += (realnum)co.H_CP_CHP_nu;

	/* H + OH =>  H2O + nu  */

	co.H_OH_H2O_nu = HMRATE(5.26E-18,-5.22,90)*findspecies("OH")->hevmol*dense.xIonDense[ipHYDROGEN][0];	
	co.hydro_sink[ipMH] += (realnum)co.H_OH_H2O_nu;

	/* Hminus + CH =>  CH2 + e  */

	co.Hminus_CH_CH2_e = HMRATE(1.0e-10,0,0)*findspecies("CH")->hevmol*hmi.Hmolec[ipMHm];		
	co.hydro_sink[ipMHm] += (realnum)co.Hminus_CH_CH2_e;

	/* Hminus + C =>  CH + e  */

	co.Hminus_C_CH_e = HMRATE(1.0e-9,0,0)*findspecies("C")->hevmol*hmi.Hmolec[ipMHm];	
	co.hydro_sink[ipMHm] += (realnum)co.Hminus_C_CH_e;

	/* Hminus + OH =>  H2O + e  */

	co.Hminus_OH_H2O_e = HMRATE(1.0e-10,0,0)*findspecies("OH")->hevmol*hmi.Hmolec[ipMHm];			
	co.hydro_sink[ipMHm] += (realnum)co.Hminus_OH_H2O_e;

	/* Hminus + O =>  OH + e  */

	co.Hminus_O_OH_e = HMRATE(1.0e-9,0,0)*findspecies("O")->hevmol*hmi.Hmolec[ipMHm];				
	co.hydro_sink[ipMHm] += (realnum)co.Hminus_O_OH_e;

	/* H2 + C =>  CH2 + nu  */

	co.H2_C_CH2_nu = HMRATE(1.0e-17,0,0)*findspecies("C")->hevmol*hmi.Hmolec[ipMH2g];		
	co.hydro_sink[ipMH2g] += (realnum)co.H2_C_CH2_nu;

	/* H2 + CP =>  CH2P + nu  */

	co.H2_CP_CH2P_nu = HMRATE(4.0e-16,-0.2,0)*findspecies("C+")->hevmol*hmi.Hmolec[ipMH2g];			
	co.hydro_sink[ipMH2g] += (realnum)co.H2_CP_CH2P_nu;

	/* H2 + SiP =>  SiH2P + nu  */

	co.H2_SiP_SiH2P_nu = HMRATE(3.0e-18,0,0)*findspecies("Si+")->hevmol*hmi.Hmolec[ipMH2g];	
	co.hydro_sink[ipMH2g] += (realnum)co.H2_SiP_SiH2P_nu;

	/* H2 + O2 =>  OH + OH */

	co.H2_O2_OH_OH = HMRATE(3.16e-10,0,21890)*findspecies("O2")->hevmol*hmi.Hmolec[ipMH2g];
	co.hydro_sink[ipMH2g] += (realnum)co.H2_O2_OH_OH;

	/* HeP + CH =>  CP + He + H  */

	co.HeP_CH_CP_He_H = HMRATE(1.1e-9,0,0)*findspecies("CH")->hevmol*dense.xIonDense[ipHELIUM][1];	
	co.hydro_source[ipMH] += (realnum)co.HeP_CH_CP_He_H;

	/* HeP + CH2 =>  CHP + He + H  */

	co.HeP_CH2_CHP_He_H = HMRATE(7.5e-10,0,0)*findspecies("CH2")->hevmol*dense.xIonDense[ipHELIUM][1];	
	co.hydro_source[ipMH] += (realnum)co.HeP_CH2_CHP_He_H;

	/* HeP + OH =>  OP + He + H  */

	co.HeP_OH_OP_He_H = HMRATE(1.1e-9,0,0)*findspecies("OH")->hevmol*dense.xIonDense[ipHELIUM][1];	
	co.hydro_source[ipMH] += (realnum)co.HeP_OH_OP_He_H;

	/* HeP + H2O  OHP + He + H  */

	co.HeP_H2O_OHP_He_H = HMRATE(2.86e-10,0,0)*findspecies("H2O")->hevmol*dense.xIonDense[ipHELIUM][1];	
	co.hydro_source[ipMH] += (realnum)co.HeP_H2O_OHP_He_H;

	/* HeP + SiH =>  SiP + He + H  */

	co.HeP_SiH_SiP_He_H = HMRATE(1.8e-9,0,0)*findspecies("SiH")->hevmol*dense.xIonDense[ipHELIUM][1];	
	co.hydro_source[ipMH] += (realnum)co.HeP_SiH_SiP_He_H;

	/* HeP + H2O =>  OH + He + HP */

	co.HeP_H2O_OH_He_HP = HMRATE(2.04e-10,0,0)*findspecies("H2O")->hevmol*dense.xIonDense[ipHELIUM][1];	
	co.hydro_source[ipMHp] += (realnum)co.HeP_H2O_OH_He_HP;

	/* HeP + CH2 =>  CP + He + H2  */

	co.HeP_CH2_CP_He_H2 = HMRATE(7.5e-10,0,0)*findspecies("CH2")->hevmol*dense.xIonDense[ipHELIUM][1];	
	co.hydro_source[ipMH2g] += (realnum)co.HeP_CH2_CP_He_H2;

	/* crnu + CH =>  C + H  */

	co.crnu_CH_C_H = findspecies("CH")->hevmol*secondaries.csupra[ipHYDROGEN ][0] * 2. * 756;
	co.hydro_source[ipMH] += (realnum)co.crnu_CH_C_H;

	/* crnu + CHP =>  CP + H  */

	co.crnu_CHP_CP_H = findspecies("CH+")->hevmol*secondaries.csupra[ipHYDROGEN][0] * 2. * 183;	
	co.hydro_source[ipMH] += (realnum)co.crnu_CHP_CP_H;

	/* crnu + H2O =>  OH + H  */

	co.crnu_H2O_OH_H = findspecies("H2O")->hevmol*secondaries.csupra[ipHYDROGEN][0] * 2. * 979;	
	co.hydro_source[ipMH] += (realnum)co.crnu_H2O_OH_H;

	/* crnu + OH =>  O + H  */

	co.crnu_OH_O_H = findspecies("OH")->hevmol*secondaries.csupra[ipHYDROGEN][0] * 2. *522;
	co.hydro_source[ipMH] += (realnum)co.crnu_OH_O_H;

	/* crnu + SiH =>  Si + H  */

	co.crnu_SiH_Si_H = findspecies("SiH")->hevmol*secondaries.csupra[ipHYDROGEN][0] * 2. *500;	
	co.hydro_source[ipMH] += (realnum)co.crnu_SiH_Si_H;

	/* nu + CH =>  C + H  */

	co.nu_CH_C_H = HMRATE(8.6e-10,0,0)*(hmi.UV_Cont_rel2_Habing_TH85_depth/1.66)*findspecies("CH")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.nu_CH_C_H;

	/* nu + CHP =>  CP + H  */

	co.nu_CHP_CP_H = HMRATE(2.5e-10,0,0)*(hmi.UV_Cont_rel2_Habing_TH85_depth/1.66)*findspecies("CH+")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.nu_CHP_CP_H;

	/* nu + CH2 =>  CH + H  */

	co.nu_CH2_CH_H = HMRATE(7.2e-10,0,0)*(hmi.UV_Cont_rel2_Habing_TH85_depth/1.66)*findspecies("CH2")->hevmol;	
	co.hydro_source[ipMH] += (realnum)co.nu_CH2_CH_H;

	/* nu + CH2P =>  CHP + H  */

	co.nu_CH2P_CHP_H = HMRATE(1.7e-9,0,0)*(hmi.UV_Cont_rel2_Habing_TH85_depth/1.66)*findspecies("CH2+")->hevmol;		
	co.hydro_source[ipMH] += (realnum)co.nu_CH2P_CHP_H;

	/* nu + CH3P =>  CH2P + H  */

	co.nu_CH3P_CH2P_H = HMRATE(1.0e-9,0,0)*(hmi.UV_Cont_rel2_Habing_TH85_depth/1.66)*findspecies("CH3+")->hevmol;		
	co.hydro_source[ipMH] += (realnum)co.nu_CH3P_CH2P_H;

	/* nu + CH3P =>  CHP + H2  */

	co.nu_CH3P_CHP_H2 = HMRATE(1.0e-9,0,0)*(hmi.UV_Cont_rel2_Habing_TH85_depth/1.66)*findspecies("CH3+")->hevmol;			
	co.hydro_source[ipMH2g] += (realnum)co.nu_CH3P_CHP_H2;

	/* nu + H2O =>  OH + H  */

	co.nu_H2O_OH_H = HMRATE(5.9e-10,0,0)*(hmi.UV_Cont_rel2_Habing_TH85_depth/1.66)*findspecies("H2O")->hevmol;		
	co.hydro_source[ipMH] += (realnum)co.nu_H2O_OH_H;

	/* nu + OH =>  O + H  */

	co.nu_OH_O_H = HMRATE(3.5e-10,0,0)*(hmi.UV_Cont_rel2_Habing_TH85_depth/1.66)*findspecies("OH")->hevmol;			
	co.hydro_source[ipMH] += (realnum)co.nu_OH_O_H;

	/* nu + OHP =>  O + HP */

	co.nu_OHP_O_HP = HMRATE(1.0e-12,0,0)*(hmi.UV_Cont_rel2_Habing_TH85_depth/1.66)*findspecies("OH+")->hevmol;		
	co.hydro_source[ipMHp] += (realnum)co.nu_OHP_O_HP;

	/* nu + SiH =>  Si + H  */

	co.nu_SiH_Si_H = HMRATE(2.8e-9,0,0)*(hmi.UV_Cont_rel2_Habing_TH85_depth/1.66)*findspecies("SiH")->hevmol;	
	co.hydro_source[ipMH] += (realnum)co.nu_SiH_Si_H;

	/* e + CHP =>  C + H  */

	co.e_CHP_C_H = HMRATE(1.5e-7,-0.42,0)*dense.eden*findspecies("CH+")->hevmol;	
	co.hydro_source[ipMH] += (realnum)co.e_CHP_C_H;

	/* e + CH2P =>  CH + H  */

	co.e_CH2P_CH_H = HMRATE(1.6e-7,-0.6,0)*dense.eden*findspecies("CH2+")->hevmol;	
	co.hydro_source[ipMH] += (realnum)co.e_CH2P_CH_H;

	/* e + CH2P =>  C + H + H  */

	co.e_CH2P_C_H_H = HMRATE(4.03e-7,-0.6,0)*dense.eden*findspecies("CH2+")->hevmol;	
	co.hydro_source[ipMH] += (realnum)co.e_CH2P_C_H_H;

	/* e + CH2P =>  C + H2  */

	co.e_CH2P_C_H2 = HMRATE(7.68e-8,-0.6,0)*dense.eden*findspecies("CH2+")->hevmol;	
	co.hydro_source[ipMH2g] += (realnum)co.e_CH2P_C_H2;

	/* e + CH3P =>  C + H2 + H  */

	co.e_CH3P_C_H2_H = HMRATE(1.05e-7,-0.5,0)*dense.eden*findspecies("CH3+")->hevmol;		
	co.hydro_source[ipMH] += (realnum)co.e_CH3P_C_H2_H;
	co.hydro_source[ipMH2g] += (realnum)co.e_CH3P_C_H2_H;

	/* e + CH3P =>  CH2 + H  */

	co.e_CH3P_CH2_H	= HMRATE(1.4e-7,-0.5,0)*dense.eden*findspecies("CH3+")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.e_CH3P_CH2_H;

	/* e + CH3P =>  CH + H + H  */

	co.e_CH3P_CH_H_H = HMRATE(5.6e-8,-0.5,0)*dense.eden*findspecies("CH3+")->hevmol;	
	co.hydro_source[ipMH] += (realnum)co.e_CH3P_CH_H_H;

	/* e + CH3P =>  CH + H2  */

	co.e_CH3P_CH_H2	= HMRATE(4.9e-8,-0.5,0)*dense.eden*findspecies("CH3+")->hevmol;
	co.hydro_source[ipMH2g] += (realnum)co.e_CH3P_CH_H2;

	/* e + H2OP =>  OH + H  */

	co.e_H2OP_OH_H = HMRATE(7.92e-8,-0.5,0)*dense.eden*findspecies("H2O+")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.e_H2OP_OH_H;

	/* e + H2OP =>  O + H + H  */

	co.e_H2OP_O_H_H = HMRATE(2.45e-7,-0.5,0)*dense.eden*findspecies("H2O+")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.e_H2OP_O_H_H;

	/* e + H2OP =>  O + H2  */

	co.e_H2OP_O_H2 = HMRATE(3.6e-8,-0.5,0)*dense.eden*findspecies("H2O+")->hevmol;
	co.hydro_source[ipMH2g] += (realnum)co.e_H2OP_O_H2;

	/* e + H3OP =>  H2O + H  */

	co.e_H3OP_H2O_H = HMRATE(1.08e-7,-0.5,0)*dense.eden*findspecies("H3O+")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.e_H3OP_H2O_H;

	/* e + H3OP =>  OH + H + H  */

	co.e_H3OP_OH_H_H = HMRATE(2.58e-7,-0.5,0)*dense.eden*findspecies("H3O+")->hevmol;	
	co.hydro_source[ipMH] += (realnum)co.e_H3OP_OH_H_H;

	/* e + H3OP =>  OH + H2  */

	co.e_H3OP_OH_H2	= HMRATE(6.45e-8,-0.5,0)*dense.eden*findspecies("H3O+")->hevmol;
	co.hydro_source[ipMH2g] += (realnum)co.e_H3OP_OH_H2;

	/* e + H3OP =>  O + H2 + H  */

	co.e_H3OP_O_H2_H = HMRATE(5.59e-9,-0.5,0)*dense.eden*findspecies("H3O+")->hevmol;	
	co.hydro_source[ipMH] += (realnum)co.e_H3OP_O_H2_H;

	/* e + HCOP =>  CO + H  */

	co.e_HCOP_CO_H = HMRATE(1.1e-7,-1,0)*dense.eden*findspecies("HCO+")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.e_HCOP_CO_H;

	/* e + OHP =>  O + H  */

	co.e_OHP_O_H = HMRATE(3.75e-8,-0.5,0)*dense.eden*findspecies("OH+")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.e_OHP_O_H;

	/* e + SiH2P =>  SiH + H  */

	co.e_SiH2P_SiH_H = HMRATE(1.5e-7,-0.5,0)*dense.eden*findspecies("SiH2+")->hevmol;	
	co.hydro_source[ipMH] += (realnum)co.e_SiH2P_SiH_H;

	/* e + SiH2P =>  Si + H + H  */

	co.e_SiH2P_Si_H_H = HMRATE(2.0e-7,-0.5,0)*dense.eden*findspecies("SiH2+")->hevmol;	
	co.hydro_source[ipMH] += (realnum)co.e_SiH2P_Si_H_H;

	/* e + SiH2P =>  Si + H2  */

	co.e_SiH2P_Si_H2 = HMRATE(1.5e-7,-0.5,0)*dense.eden*findspecies("SiH2+")->hevmol;	
	co.hydro_source[ipMH2g] += (realnum)co.e_SiH2P_Si_H2;

	/* e + SiOHP =>  SiO + H  */

	co.e_SiOHP_SiO_H = HMRATE(1.5e-7,-0.5,0)*dense.eden*findspecies("SiOH+")->hevmol;	
	co.hydro_source[ipMH] += (realnum)co.e_SiOHP_SiO_H;

	/* H2 + CH =>  CH3 + nu  */

	co.H2_CH_CH3_nu	= HMRATE(5.09E-18,-0.71,11.6)*findspecies("H2")->hevmol*findspecies("CH")->hevmol;
	co.hydro_sink[ipMH2g] += (realnum)co.H2_CH_CH3_nu;

	/* H2 + CH3P =>  CH5P + nu  */

	co.H2_CH3P_CH5P_nu = HMRATE(1.3e-15,-1,0)*findspecies("CH3+")->hevmol*hmi.Hmolec[ipMH2g];	
	co.hydro_sink[ipMH2g] += (realnum)co.H2_CH3P_CH5P_nu;

	/* H2s + CH =>  CH3 + nu  */

	co.H2s_CH_CH3_nu = HMRATE(5.09E-18,-0.71,0)*findspecies("CH")->hevmol*hmi.Hmolec[ipMH2s]*hmi.lgLeiden_Keep_ipMH2s;	
	co.hydro_sink[ipMH2s] += (realnum)co.H2s_CH_CH3_nu;

	/* Hminus + CH2 =>  CH3 + e  */

	co.Hminus_CH2_CH3_e	= HMRATE(1.0e-9,0,0)*findspecies("CH2")->hevmol*hmi.Hmolec[ipMHm];	
	co.hydro_sink[ipMHm] += (realnum)co.Hminus_CH2_CH3_e;

	/* Hminus + CH3 =>  CH4 + e  */

	co.Hminus_CH3_CH4_e	= HMRATE(1.0e-9,0,0)*findspecies("CH3")->hevmol*hmi.Hmolec[ipMHm];	
	co.hydro_sink[ipMHm] += (realnum)co.Hminus_CH3_CH4_e;

	/* nu + CH3 =>  CH2 + H  */

	co.nu_CH3_CH2_H	= HMRATE(2.5e-10,0,0)*(hmi.UV_Cont_rel2_Habing_TH85_depth/1.66)*findspecies("CH3")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.nu_CH3_CH2_H;

	/* nu + CH3 =>  CH + H2  */

	co.nu_CH3_CH_H2	= HMRATE(2.5e-10,0,0)*(hmi.UV_Cont_rel2_Habing_TH85_depth/1.66)*findspecies("CH3")->hevmol;
	co.hydro_source[ipMH2g] += (realnum)co.nu_CH3_CH_H2;

	/* nu + CH4 =>  CH3 + H  */

	co.nu_CH4_CH3_H	= HMRATE(2.2e-10,0,0)*(hmi.UV_Cont_rel2_Habing_TH85_depth/1.66)*findspecies("CH4")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.nu_CH4_CH3_H;

	/* nu + CH4 =>  CH2 + H2  */

	co.nu_CH4_CH2_H2 = HMRATE(9.8e-10,0,0)*(hmi.UV_Cont_rel2_Habing_TH85_depth/1.66)*findspecies("CH4")->hevmol;	
	co.hydro_source[ipMH2g] += (realnum)co.nu_CH4_CH2_H2;

	/* nu + CH4 =>  CH + H2  */

	co.nu_CH4_CH_H2	= HMRATE(2.2e-10,0,0)*(hmi.UV_Cont_rel2_Habing_TH85_depth/1.66)*findspecies("CH4")->hevmol;
	co.hydro_source[ipMH2g] += (realnum)co.nu_CH4_CH_H2;

	/* crnu + CH3 =>  CH2 + H  */

	co.crnu_CH3_CH2_H = findspecies("CH3")->hevmol*secondaries.csupra[ipHYDROGEN][0] * 2.*500;	
	co.hydro_source[ipMH] += (realnum)co.crnu_CH3_CH2_H;

	/* crnu + CH3 =>  CH + H2  */

	co.crnu_CH3_CH_H2 = findspecies("CH3")->hevmol*secondaries.csupra[ipHYDROGEN][0] * 2.*500;	
	co.hydro_source[ipMH2g] += (realnum)co.crnu_CH3_CH_H2;

	/* crnu + CH4 =>  CH2 + H2  */

	co.crnu_CH4_CH2_H2 = findspecies("CH4")->hevmol*secondaries.csupra[ipHYDROGEN][0] * 2.*2272;	
	co.hydro_source[ipMH2g] += (realnum)co.crnu_CH4_CH2_H2;

	/* e + CH5P =>  CH3 + H2  */

	co.e_CH5P_CH3_H2 = HMRATE(5.5e-7,-0.3,0)*dense.eden*findspecies("CH5+")->hevmol;	
	co.hydro_source[ipMH2g] += (realnum)co.e_CH5P_CH3_H2;

	/* e + CH5P =>  CH4 + H  */

	co.e_CH5P_CH4_H	= HMRATE(5.5e-7,-0.3,0)*dense.eden*findspecies("CH5+")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.e_CH5P_CH4_H;

	/* e + CH4P =>  CH3 + H  */

	co.e_CH4P_CH3_H	= HMRATE(1.75e-7,-0.5,0)*dense.eden*findspecies("CH4+")->hevmol;
	co.hydro_source[ipMH] += (realnum)co.e_CH4P_CH3_H;

	/* e + CH4P =>  CH2 + H + H  */

	co.e_CH4P_CH2_H_H = HMRATE(1.75e-7,-0.5,0)*dense.eden*findspecies("CH4+")->hevmol;	
	co.hydro_source[ipMH] += (realnum)co.e_CH4P_CH2_H_H;
#	endif
	/******* END OF REACTIONS ********/
#	if 0
	if( iteration>1 )
	{
		r = rlist;
		i = 0;
		while (r->next != NULL)
		{
			++i;
			r = r->next;
			fprintf(ioQQQ,"DEBUG r\t%li\t%.3e\n", i, r->rk);
		}
	}
#	endif

	/* Generate chemical error vector and Jacobian array from reaction list */
	for( i=0; i < N_H_MOLEC; i++ )
	{
		bvec[i] = 0.;
		for( j=0; j < N_H_MOLEC; j++ )
		{
			c[j][i] = 0.;
		}
	}
	/* Subtotal rates for H_0 and H^+ within ipMHo */
	for(i=0;i<2;i++) 
	{
		mole.source[ipHYDROGEN][i] = mole.sink[ipHYDROGEN][i] = 0.;
	}

	/* reinitialize linked list and move through it */
	/* set up Jacobian matrix */
	r = rlist;

	/*  Comments about this line of code, made by Nick Abel.  This section of the code
		does the same thing as in co.c, where it decouples products of two densities
		that the code is trying to predict. Reactions of the form n(X)*n(Y) cannot be 
		treated properly in a linear solver otherwise.

		This section of the code starts off with a while
		statement that loops over all reactions in the network (hence the r->next variable).
		For each reaction, the first thing that is done is to state that for the reaction
		of interest, the rate coefficient k is equal to r->rk, which is the stored value from
		when the reaction was originally set up above.  After this a loop over all reactants is 
		performed.  This is actually a double loop.  If the number of reactants is greater than 
		one, then this code will generate (via the i!=j if statement) two "decoupled" products of 
		the rate coefficient times the previous solution for the density.  This is explained in co.c 
		and in that code it is just the reaction rates that end in _1 or _2. 

		After this loop the	variable "rate" is multiplied by the the density that was not multiplied
		in the loop.  This sets up "rate" to be what co.c calls "bvec" for that reaction. The bvec's 
		are needed because in the process of "decoupling" a reaction, their is a leftover term of
		the form k*n(X)old*n(Y)old which has to go into the solution vector for the matrix equation Ax = b.

		The next two loops set up the bvec so that it goes with the proper species.  Also, in the
		case of H and H+, these rates are saved so that they can be fed into the main ionization 
		solver.  This allows Cloudy to account for formation and destruction processes
		for Hydrogen that are due to reaction with molecules.  

		The last for statement inside the while loop takes the rates that were calculated previously 
		and stores them in the appropriate part of the matrix A. */

	/* rates complete - all reactions have been stored, now fill in the c[][] matrix */
	while (r->next != NULL)
	{
		r = r->next;
		/* >>chng 04 feb 05, this was option to cut chemistry short in testing
		if(r->index == rindex)
			break;*/
		rk = r->rk;

		/* if this blows, rk is NaN */
		/*ASSERT( rk == rk );*/
		ASSERT( !isnan( rk ) );

		/* There's an O(n) algorithm for this -- but it doesn't improve
		 * things unless nreactants is >= 4...!*/
		/* loop over all rate determining species */
		for(i=0;i<r->nrates;i++)
		{
			rate_deriv[i] = rk;
			for(j=0;j<r->nrates;j++)
			{
				/* Hmolec_old was previous abundance,
				 * rate_deriv[i] is derivative of rates coefficient by species r->rate_species[i] */
				if(i!=j)
				{
					rate_deriv[i] *= Hmolec_old[r->rate_species[j]];
					/* if this blows, rate_deriv[i] is NaN */
					/*ASSERT( rate_deriv[i] == rate_deriv[i] );*/
					ASSERT( !isnan( rate_deriv[i] ) );
				}
			}
		}

		/* this is total rate, rate_deriv times old population */
		rate = rate_deriv[0]*Hmolec_old[r->rate_species[0]];

		/* is this blows, rate is NaN */
		ASSERT( !isnan( rate ) );

		/* Get sink terms (i.e. rate/abundance) to include in ionization ladders */
		for(i=0;i < r->nreactants;i++)
		{
			int ok = 0;
			for(j=0;j < r->nrates && !ok;j++)
			{
				if(r->rate_species[j] == r->reactants[i]) 
				{
					sinkrate[i] = rate_deriv[j];
					ok = true;
				}
			}
			if(!ok) 
			{
				/* Odd, the rate didn't depend on one of the species it used
				 * at all!  An alternative way of getting sink rate is
				 *
				 * sinkrate[i] = rate/Hmolec_old[r->reactants[i]]; 
				 *
				 * but this uses the possibly zero Hmolec_old, and is prone to underflow of rate.
				 * */
				fprintf(ioQQQ,"A chemical rate in hmole was independent of the species it used\n");
				fprintf(ioQQQ,"This probably shouldn't happen (so you shouldn't see this message).\n");
				cdEXIT(EXIT_FAILURE);
			}
		}

		/* if(nzone == 416)
			 fprintf(ioQQQ,"Adding reaction %d rate %g\n",r->index,rate); */

		/* loop over all reactions, find total consumption rate,
		 * also keep track of rates that use up H0 or H+ */
		for(i=0;i<r->nreactants;i++)
		{
			ratei = r->reactants[i];
			bvec[ratei] -= rate;
			/*if((nzone == 421 || nzone == 422) && ratei == ipMHm)
					fprintf(ioQQQ,"snk %s %d %g\n",hmi.chLab[ratei],r->index,rate);*/
			/* mole.sink[ipHYDROGEN] is how chemical reaction network reacts with ionization
			 * network, so this keeps track of total rates */
			if(ratei == ipMH || ratei == ipMHp)
				mole.sink[ipHYDROGEN][ratei] += sinkrate[i];
		}

		/* loop over all reactions, find total production rate,
		 * also keep track of rates that produce H0 or H+ */
		for(i=0;i<r->nproducts;i++)
		{
			ratei = r->products[i];
			bvec[ratei] += rate;
			/*if((nzone == 421 || nzone == 422) && ratei == ipMHm)
					fprintf(ioQQQ,"src %s %d %g\n",hmi.chLab[ratei],r->index,rate); */
			if(ratei == ipMH || ratei == ipMHp)
			{
				mole.source[ipHYDROGEN][ratei] += rate;

				/* confirm mole.source[ipHYDROGEN][ratei] is valid float */
				ASSERT( !isnan( mole.source[ipHYDROGEN][ratei] ) );
			}
		}

		/* The first thing that must be said about the for statements below is that 
		 * it is inside a while statement (starting on line 4151).  This while loops over
		 * all reactions. Each individual reaction is stored with the pointer r->next,
		 * which goes from 1 to the number of reactions in the hmole_step.  Also, for
		 * each reaction, the code keeps track of the number of products and reactants
		 * in reaction r->next, what each product or reactant is 
		 * (H2, H2*, H3+, et cetera), and the reaction rate coefficient.
		 *
		 *	For example, if the first reaction is H2 + H+ => H2+ + H, then:		 
		 *
		 *	r->next = 1 (this is the first reaction)
		 *	r->nreactants[i] = 2 (the number of reactants equals two)
		 *	r->nproducts[i] = 2  (the number of products equals two)
		 *	r->nrates =  (unless the reactant is repeated, such as H2 + H2 => H2* + H2, 
		 *	            then r->nrates is always the number of reactants (see line 107)
		 *
		 *	rate_deriv[j] = reaction rate determined from lines (4161-4173), 
		 *				this is the product of the rate coefficient k and the density of
		 *				one of the reactants.
		 *
		 *	r->rate_species[j], r->reactants[i], and r->products[i] is the species
		 *				corresponding to each product or reactant.
		 *
		 *
		 *
		 *	So H2 + H+ => H2+ + H does the following in the code below.
		 *
		 *
		 *	1)  Since nrates = 2, the for loop goes over all reactants, 1 and 2. 
		 *		The first reactant considered is 1
		 *	2)  ratej is set equal to H2, and rated is equal to k*[density of H+]
		 *	3)  The second for statement  loops over all reactants, and puts
		 *		fills in some matrix elements:
		 *
		 *		c[ipMH2g][ipMH2g] -= k*[density of H+]
		 *		c[ipMH2g][ipMHp] -=k*[density of H+]
		 *
		 *	4)  Now the third for statement fills in some more reactions,
		 *		involving the products:
		 *
		 *		c[ipMH2g][ipMH2p] += k*[density of H+]
		 *		c[ipMH2g][ipMH] += k*[density of H+]
		 *
		 *	5)  At this state, we go back up to the first for statement,
		 *		and the reactant is changed from 1 to 2 (H+)
		 *	6)  Also ratej is set to H+ and rated is now k*[density of H2]
		 *	7)  Some more matrix elements are filled in:
		 *
		 *		c[ipMHp][ipMH2g] -= k*[density of H2]
		 *		c[ipMHp][ipMHp] -= k*[density of H2]
		 *
		 *		c[ipMHp][ipMH2p] += k*[density of H2]
		 *		c[ipMHp][ipMH] += k*[density of H2]
		 *
		 *		This is the more elegant way of linearizing the series of non-linear
		 *		equations in the molecular network, incorporated by Robin Williams.
		 *		co.c does the same thing more explicitly, but also takes up way too 
		 *		much space. */


		/* fill Jacobian rate matrix */
		for(j=0;j<r->nrates;j++)
		{
			ratej = r->rate_species[j];
			rated = rate_deriv[j];
			for(i=0;i<r->nreactants;i++)
			{
				c[ratej][r->reactants[i]] -= rated;
			}
			for(i=0;i<r->nproducts;i++)
			{
				c[ratej][r->products[i]] += rated;
			}
		}
	}
	/* the c[][] matrix Jacobian array has now been filled with all reagents */

	/* save rate H2 is destroyed units s-1 */
	/* >>chng 05 mar 18, TE, add terms - 
		total destruction rate is: dest_tot = n_H2g/n_H2tot * dest_H2g + n_H2s/n_H2tot * dest_H2s */
	/* as reactions that change H2s to H2g and vice versa are not counted destruction processes, the terms c[ipMH2g][ipMH2s] *
	   and c[ipMH2s][ipMH2g], which have a different sign than [ipMH2g][ipMH2g] and [ipMH2s][ipMH2s], have to be added	*/
	if( hmi.H2_total>SMALLFLOAT )
		hmi.H2_rate_destroy = (hmi.Hmolec[ipMH2g] * (-c[ipMH2g][ipMH2g]-c[ipMH2g][ipMH2s]) +
			hmi.Hmolec[ipMH2s] * (-c[ipMH2s][ipMH2s]-c[ipMH2s][ipMH2g])) / hmi.H2_total;
	else
		hmi.H2_rate_destroy = hmi.H2_photodissoc_used_H2s;

	{
		/* following should be set true to print populations */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			if( DEBUG_LOC && (nzone > 570) ) 
			{
				printsol = 1;
				fprintf(ioQQQ,"Temperature %g\n",phycon.te);
				fprintf(ioQQQ," Net mol ion rate [%g %g] %g\n",mole.source[ipHYDROGEN][1],mole.sink[ipHYDROGEN][1],
								mole.source[ipHYDROGEN][1]-mole.sink[ipHYDROGEN][1]*Hmolec_old[ipMHp]);
			}
		}
	}

	/* save total H2P destruction rate for possible later printout:
	 * NB this must come last */
	desh2p = -c[ipMH2p][ipMH2p]; 

	/* Check that matrix and vector generated in above loops make sense */
	/*if(!defined(NDEBUG))	*/
	/* in std C NDEBUG is a macro set at compile time */
	/** \todo	0	had to comment following test out - NA change to hmole caused
	 * massive prints */
#	if 0
/*#	if !defined(NDEBUG)*/
#	ifndef NDEBUG
	{
		double total, mtotal;
		for(i=0;i<N_H_MOLEC;i++) 
		{
			total = 0.;
			for( j=0;j<N_H_MOLEC;j++) 
			{
				total += c[i][j]*hmi.nProton[j];
			}
			if( fabs(total) > 1e-5*fabs(c[i][i]*hmi.nProton[i])) 
			{
				fprintf(ioQQQ,"PROBLEM Subtotal1 %.2e\n",fabs(total)/fabs(c[i][i]*hmi.nProton[i]));
				fprintf(ioQQQ,"Species %li Total %g Diag %g\n",i,total,c[i][i]*hmi.nProton[i]);
			}
			else if( fabs(total) > 1e-6*fabs(c[i][i]*hmi.nProton[i]) && phycon.te< 1e6 ) 
			{
				fprintf(ioQQQ,"NOTE Subtotal1 %.2e Te=%.4e\n",
					fabs(total)/fabs(c[i][i]*hmi.nProton[i]),phycon.te);
				fprintf(ioQQQ,"Species %li Total %g Diag %g\n",i,total,c[i][i]*hmi.nProton[i]);
			}
		}
		total = mtotal = 0.;
		for(j=0;j<N_H_MOLEC;j++) 
		{ 
			total += bvec[j]*hmi.nProton[j]; 
			mtotal += fabs(bvec[j]*hmi.nProton[j]); 
		}
		if(fabs(total) > 1e-30 && fabs(total) > 1e-10*rtot) 
		{ 
			fprintf(ioQQQ,"PROBLEM Subtotal2 %.2e\n",fabs(total)/mtotal);
			fprintf(ioQQQ,"RHS Total %g cf %g\n",total,mtotal);
		} 
		else if(fabs(total) > 1e-7*mtotal)  
		{
			fprintf(ioQQQ,"WARNING zone %li Hmole RHS conservation error %.2e of %.2e\n",nzone,total,mtotal);
			fprintf(ioQQQ,"(may be due to high rate equilibrium reactions)\n");
		}
	}
#		endif
#	endif


#define MOLMIN  1
#define N_H_MAT (N_H_MOLEC-MOLMIN)
	/* Will collapse ipMH and ipMHp into single species, as don't include
	 * all ionizations and recombinations here */
	/* last test - do not include advection if we have overrun the radius scale 
	 * of previous iteration */
	/* >>chng 06 mar 17, comment out test on old full depth - keep old solution if overrun scale */
	if( iteration >= dynamics.n_initial_relax+1  && dynamics.lgAdvection && !dynamics.lgEquilibrium
		&& dynamics.Rate != 0. ) 
	{
		/* Don't use conservation form in matrix solution */
		ipConserve = -1; 
		/* Add rate terms for dynamics to equilibrium, makes c[][] non-singular */
		for(i=0;i<N_H_MOLEC;i++) 
		{
			c[i][i] -= dynamics.Rate;
			bvec[i] -= (Hmolec_old[i]*dynamics.Rate-dynamics.H2_molec[i]);
		}

		/* Dynamics implies conservation of advected material */
		proton_sum_old = 0.;
		for(i=0; i<N_H_MOLEC;i++)
		{
			proton_sum_old += hmi.nProton[i]*dynamics.H2_molec[i]/dynamics.Rate;
		}

		/* bring H+ and H0 together since their ratio is set in H atom solver,
		 * we determine sum of two here */
		for(i=0;i<N_H_MOLEC;i++) 
		{
			c[ipMHp][i] = (Hmolec_old[ipMH]*c[ipMH][i]+Hmolec_old[ipMHp]*c[ipMHp][i])/SDIV(sum_H0_Hp);
			c[ipMH][i] = 0.;
		}
		for(i=1;i<N_H_MOLEC;i++) 
		{
			c[i][ipMHp] += c[i][ipMH];
			c[i][ipMH] = 0.;
		}
		bvec[ipMHp] += bvec[ipMH];
		bvec[ipMH] = 0.;
		Hmolec_old[ipMHp] += Hmolec_old[ipMH];
		Hmolec_old[ipMH] = 0.;
	}
	else
	{
		/* usual branch, no advection */
		/* bring H+ and H0 together since their ratio is set in H atom solver,
		 * we determine sum of two here */
		for(i=0;i<N_H_MOLEC;i++) 
		{
			/* >>chng 04 feb 04, sum_H0_Hp goes to zero when no ionization,
			 * add test on SMALLFLOAT */
			if( sum_H0_Hp > SMALLFLOAT )
				c[ipMHp][i] = (Hmolec_old[ipMH]*c[ipMH][i]+Hmolec_old[ipMHp]*c[ipMHp][i])/sum_H0_Hp;
			c[ipMH][i] = 0.;
		}
		Hmolec_old[ipMHp] += Hmolec_old[ipMH];
		bvec[ipMH] = Hmolec_old[ipMH] = 0.;
		ipConserve = ipMHp;
		/* For Newton-Raphson method, want the change in populations to be zero,
		 * so the conserved component must also be zero */
		bvec[ipConserve] = 0.;

		/* proton_sum_old is the sum of all protons in H-bearing molecules */
		proton_sum_old = 0.;
		for(i=MOLMIN;i<N_H_MOLEC;i++) 
		{
			c[i][ipConserve] = hmi.nProton[i];
			proton_sum_old += hmi.nProton[i]*Hmolec_old[i];
		}
	}

	{
		/* following should be set true to print populations */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			/* these are the raw results */
			fprintf( ioQQQ, " HMOLE h2 %.2e h2* %.2e\n" , Hmolec_old[ipMH2g] ,Hmolec_old[ipMH2s] );
		}
	}

	/*------------------------------------------------------------------ */
	if(printsol || (trace.lgTrace && trace.lgTr_H2_Mole ))
	{

		/*

		[0][0]  [0][1]  [0][2]  [0][3]  [0][4]  [0][5]
		[1][0]  [1][1]  [1][2]  [1][3]  [1][4]  [1][5]
		[2][0]  [2][1]  [2][2]  [2][3]  [2][4]  [2][5]
		[3][0]  [3][1]  [3][2]  [3][3]  [3][4]  [3][5]
		[4][0]  [4][1]  [4][2]  [4][3]  [4][4]  [4][5]
		[5][0]  [5][1]  [5][2]  [5][3]  [5][4]  [5][5]

		[ipMHo][ipMHo]  [ipMHo][ipMHm]  [ipMHo][ipMH2g]  [ipMHo][ipMH2p]  [ipMHo][ipMH3p]  [ipMHo][ipMH2s]
		[ipMHm][ipMHo] [ipMHm][ipMHm] [ipMHm][ipMH2g] [ipMHm][ipMH2p] [ipMHm][ipMH3p] [ipMHm][ipMH2s]
		[ipMH2g][ipMHo]  [ipMH2g][ipMHm]  [ipMH2g][ipMH2g]  [ipMH2g][ipMH2p]  [ipMH2g][ipMH3p]  [ipMH2g][ipMH2s]
		[ipMH2p][ipMHo] [ipMH2p][ipMHm] [ipMH2p][ipMH2g] [ipMH2p][ipMH2p] [ipMH2p][ipMH3p] [ipMH2p][ipMH2s]
		[ipMH3p][ipMHo] [ipMH3p][ipMHm] [ipMH3p][ipMH2g] [ipMH3p][ipMH2p] [ipMH3p][ipMH3p] [ipMH3p][ipMH2s]
		[ipMH2s][ipMHo] [ipMH2s][ipMHm] [ipMH2s][ipMH2g] [ipMH2s][ipMH2p] [ipMH2s][ipMH3p]  [ipMH2s][ipMH2s]

		*/

		fprintf(ioQQQ,"       MOLE old abundances\t%.2f",fnzone);
		for( i=0; i<N_H_MOLEC; i++ )
			fprintf(ioQQQ,"\t%.2e", Hmolec_old[i] );
		fprintf(ioQQQ,"\n" );

		/* print the full matrix */
		fprintf( ioQQQ, "                ");
		for( i=MOLMIN; i < N_H_MOLEC; i++ )
		{
			fprintf( ioQQQ, "      %s", hmi.chLab[i] );
		}
		fprintf( ioQQQ, "   bvec \n" );

		for( i=MOLMIN; i < N_H_MOLEC; i++ )
		{
			fprintf( ioQQQ, "       MOLE%2ld %s", i-MOLMIN ,hmi.chLab[i] );
			for( j=MOLMIN; j < N_H_MOLEC; j++ )
			{
				fprintf( ioQQQ, "%10.2e", c[j][i] );
			}
			fprintf( ioQQQ, "%10.2e", bvec[i] );
			fprintf( ioQQQ, "\n" );
		}
	}

	/*fprintf(ioQQQ,"DEBUG %.2e %.2e %.2e %.2e\n", 
		c[ipMH][ipMH2g] , c[ipMH][ipMH2s],
		c[ipMH2g][ipMH] , c[ipMH2s][ipMH]);
	fprintf(ioQQQ,"DEBUG %.2e %.2e %.2e %.2e\n\n", 
		c[ipMHp][ipMH2g] , c[ipMHp][ipMH2s],
		c[ipMH2g][ipMHp] , c[ipMH2s][ipMHp]);*/
	/* establish local timescale for H2 molecule destruction */
	if( -c[ipMH2g][ipMH2g] > SMALLFLOAT )
	{
		/* units are now seconds */
		timesc.time_H2_Dest_here = -1./c[ipMH2g][ipMH2g];
	}
	else
	{
		timesc.time_H2_Dest_here = 0.;
	}

	/* local timescale for H2 formation 
	 * both grains and radiative attachment */
	if( dense.xIonDense[ipHYDROGEN][0]>0 )
		timesc.time_H2_Form_here = gv.rate_h2_form_grains_used_total * 
			/* this corrects for fact that we the timescale for H2 to form from an atomic gas.
			 * The rate becomes very small when gas is fully molecular, and ratio of total hydrogen
			 * to atomic hydrogen corrections for this. */
			dense.gas_phase[ipHYDROGEN]/dense.xIonDense[ipHYDROGEN][0];
	else
		timesc.time_H2_Form_here = 0.;
	timesc.time_H2_Form_here += hmi.hminus_rad_attach;
	/* timescale is inverse of this rate */
	if( timesc.time_H2_Form_here > SMALLFLOAT )
	{
		/* units are now seconds */
		timesc.time_H2_Form_here = 1./timesc.time_H2_Form_here;
	}
	else
	{
		timesc.time_H2_Form_here = 0.;
	}

#	ifdef MAT
#		undef MAT
#	endif
#	define MAT(a,I_,J_)	(*((a)+(I_)*(N_H_MAT)+(J_)))

	/* copy contents over to 1D array */
	for( j=0; j < N_H_MAT; j++ )
	{
		for( i=0; i < N_H_MAT; i++ )
		{
			MAT(amat,i,j) = c[i+MOLMIN][j+MOLMIN];
		}
	}

	if(printsol)
	{
		double total=0;
		fprintf(ioQQQ,"Zone %.2f input:\n",fnzone);
		for( i=MOLMIN; i < N_H_MOLEC; i++ )
		{
				fprintf(ioQQQ,"%15.7g\t",Hmolec_old[i]);
				total += hmi.nProton[i]*Hmolec_old[i];
		}
		fprintf(ioQQQ,"sum = %15.7g\n",total);
	}

	int32 merror1 = 0;
	int32 merror2 = 0;

	/* now invert the matrix */
	getrf_wrapper(N_H_MAT,N_H_MAT,(double*)amat,N_H_MAT,ipiv,&merror1);
	getrs_wrapper('N',N_H_MAT,1,(double*)amat,N_H_MAT,ipiv,bvec+MOLMIN,N_H_MAT,&merror2);

	if( merror1 != 0 || merror2 != 0 )
	{
		fprintf( ioQQQ, "PROBLEM hmole_step: dgetrs finds singular or ill-conditioned matrix\n" );
		cdEXIT(EXIT_FAILURE);
	}

	if(printsol)
	{
		double total=0;
		fprintf(ioQQQ,"solution:\n");
		for( i=MOLMIN; i < N_H_MOLEC; i++ )
		{
			fprintf(ioQQQ,"%15.7g\t",bvec[i]);
			total += hmi.nProton[i]*bvec[i];
		}
		fprintf(ioQQQ,"sum = %15.7g\n",total);
	}

	*error = 0;
	/* loop starts from MOLMIN=1 rather than zero since
	 * H0 and H+ rates have been collapsed into one, since that solution
	 * comes from H atom solver.  
	 *
	 * bvec is (old - new) solutions coming into this routine
	 * loops converts bvec to new density */
	for( i=MOLMIN; i < N_H_MOLEC; i++ )
	{
		/* Smooth the error mode tailoff */
		etmp = bvec[i]/(ABSLIM+Hmolec_old[i]);

		if(printsol)
			fprintf(ioQQQ,"%15.7g\t",etmp);

		/* square of change in abundance of this species, in this iteration */
		*error += etmp*etmp;
		/* change bvec from being the difference into being the new value
		 * bvec is now new density */
		bvec[i] = Hmolec_old[i]-bvec[i];
	}
	/* bvec is now the density */
	*error = sqrt(*error)/N_H_MAT;

	if(printsol)
	{
		double total=0;
		fprintf(ioQQQ,"err = %15.7g\n",*error);
		/* fprintf(ioQQQ,"derived:\n"); */
		for( i=MOLMIN; i < N_H_MOLEC; i++ )
		{
				fprintf(ioQQQ,"%15.7g\t",bvec[i]);
				total += hmi.nProton[i]*bvec[i];
		}
		fprintf(ioQQQ,"sum = %15.7g\n",total);
	}

	proton_sum_new = 0.;
	/* check for negative populations and do proton sum */
	lgNegPop = false;
	fracneg = 0.;
	fracnegfac = 0.;
	iworst = -1;
	for( i=MOLMIN; i < N_H_MOLEC; i++ )
	{
		if( bvec[i] < 0. ) 
		{
			lgNegPop = true;
		}
		/* largest relative change in solution for neg soln */
		if( Hmolec_old[i]> 0. )
			fracnegtmp = -bvec[i]/SDIV(Hmolec_old[i]);
		else
			fracnegtmp = 0.;
		/* this can only occur for negative solutions since fracneg starts
		 * as zero */
		if(fracnegtmp > fracneg) 
		{
			fracneg = fracnegtmp;
			fracnegfac = (0.5*Hmolec_old[i]-bvec[i])/(Hmolec_old[i]-bvec[i]);
			iworst = i;
		}
		/* sum total number of protons used - hmi.nProton is number of protons in species bvec[i] */
		proton_sum_new += hmi.nProton[i] * bvec[i];
	}

	/* this is difference between number of protons in hmi.Hmolec upon entry into this routine
	 * and number of protons we found here */
	conserve = (proton_sum_old - proton_sum_new)/SDIV(proton_sum_old);
	/* >>chng 06 jun 29, from conserve < 1e-8 to twice FLT_EPSILON - the CO network now includes
	 * part of the H - the old upstream fraction of H in CO molecules is likely different from
	 * the current fractions.  the CO chem is only solved to a certain error, should not
	 * demand higher accuracy than this 
	 * the factor 10.*FLT_EPSILON also appears in ion_solver in total H conservation */
	/*if( fabs(conserve) > 1e-8 )*/
	if( fabs(conserve) > 10.*FLT_EPSILON )
		fprintf(ioQQQ,"PROBLEM hmoleee zn %li proton_sum_old %.8e, proton_sum_new %.8e n(H) %.8e (old-new)/old %.3e nH-old %.3e nH-new %.3e\n", 
		nzone , 
		proton_sum_old , 
		proton_sum_new , 
		dense.gas_phase[ipHYDROGEN] , 
		conserve ,
		dense.gas_phase[ipHYDROGEN]-proton_sum_old,
		dense.gas_phase[ipHYDROGEN]-proton_sum_new);

#	if 0
	/* NDEBUG is set by the compiler to indicate that a debugging mode
	 * has not been specified.  */
#	ifndef NDEBUG
	/*if(NDEBUG)*/
	{
		fprintf( ioQQQ, "Zone %li\n",nzone);
		for( i=MOLMIN; i < N_H_MOLEC; i++ )
		{
			fprintf(ioQQQ," %s %.2e", hmi.chLab[i], Hmolec_old[i]);
		}
		fprintf( ioQQQ, " =>\n" );
		for( i=MOLMIN; i < N_H_MOLEC; i++ )
		{
			fprintf(ioQQQ," %s %.2e", hmi.chLab[i], bvec[i]);
		}
		fprintf( ioQQQ, "\n" );
	}
#	endif
#	endif

	if(lgNegPop)
	{
#		ifndef NDEBUG
		/* very common to obtain negative solution on very first try - 
		 * don't print in this case */
		if(*nFixup )
		{
			fprintf( ioQQQ, " PROBLEM  hmole_step finds negative H molecule, in zone %.2f.\n",fnzone );
			fprintf( ioQQQ, " Worst is species %d -ve by fraction %g.\n",iworst,fracneg );
			fprintf( ioQQQ, " The populations follow:\n");
			for( i=MOLMIN; i < N_H_MOLEC; i++ )
			{
				fprintf(ioQQQ," %s %.2e", hmi.chLab[i], bvec[i]);
			}
			fprintf( ioQQQ, "\n" );
		}
#		endif

		/* Fix negative abundance -- assume the new solution is better in some ways */
		{
			double total=0., ntotal=0., ratio;
			enum {FIXUPTYPE = 1};

			++*nFixup;

			if(FIXUPTYPE == 1) {
				for( i=MOLMIN; i < N_H_MOLEC; i++ )
				{
					total += hmi.nProton[i]*bvec[i];
					if(bvec[i] < 0) 
					{
						ntotal += hmi.nProton[i]*bvec[i];
						bvec[i] = 0.;
					}
				}
				ratio = total/(total-ntotal);
				for( i=MOLMIN; i < N_H_MOLEC; i++ )
				{
					bvec[i] *= ratio;
				}
			}
			else if(FIXUPTYPE == 2) 
			{
				for( i=MOLMIN; i < N_H_MOLEC; i++ )
				{
					bvec[i] = fracnegfac*Hmolec_old[i]+(1-fracnegfac)*bvec[i];
				}
			}

#			ifndef NDEBUG
			/*if(NDEBUG)*/
			/* very common to obtain negative solution on very first try - 
			 * don't print in this case */
			if( *nFixup>1 )
			{
				fprintf(ioQQQ," FIXUP taken %i time%s.\n\n", *nFixup, (*nFixup == 1)?"":"s");
				fprintf( ioQQQ, " Initially:\n");
				for( i=MOLMIN; i < N_H_MOLEC; i++ )
				{
					fprintf(ioQQQ," %s %.2e", hmi.chLab[i], Hmolec_old[i]);
				}
				fprintf( ioQQQ, "\n" );
				fprintf( ioQQQ, " Changed to:\n");
				for( i=MOLMIN; i < N_H_MOLEC; i++ )
				{
					fprintf(ioQQQ," %s %.2e", hmi.chLab[i], bvec[i]);
				}
				fprintf( ioQQQ, "\n" );
			}
#			endif
		}
	}

	/* put derived abundances back into appropriate molecular species,
	 * bvec[ipMHp] is still sum of H0 and H+ from chemistry,
	 * now split up using ratio found in H atom solver */
	h1fnd = bvec[ipMHp];
	/* >>chng 04 feb 04, add SMALLFLOAT to protect against fully molecular limit */
	h1rat = h1fnd/SDIV(sum_H0_Hp);
	/* put back into proper places in solution vector so following loops work
	 * as expected */
	bvec[ipMH] = dense.xIonDense[ipHYDROGEN][0] * h1rat;
	bvec[ipMHp] = dense.xIonDense[ipHYDROGEN][1] * h1rat;
	/* ASSERT(fabs(bvec[ipMH]+bvec[ipMHp]-h1fnd) < 1e-12 * h1fnd); */

	if(fabs(bvec[ipMH]+bvec[ipMHp]-h1fnd) >= 1e-12 * h1fnd) 
	{
		static bool lgPrint=true;
		fprintf(ioQQQ,"PROBLEM h1fnd residual error, bvec[ipMH]=%g [ipMHp}=%g" 
			" h1fnd=%g h1rat=%g bvec[ipMH]+bvec[ipMHp]-h1fnd=%g\n",
			bvec[ipMH],bvec[ipMHp],h1fnd,h1rat,bvec[ipMH]+bvec[ipMHp]-h1fnd);
		/* nearly all cases of this problem are due to ZERO H ionization rate - this can't happen if
		 * cosmic rays are present */
		if( lgPrint )
		{
			fprintf(ioQQQ," This problem is usually caused by little or no sources of ionization.\n");
			fprintf(ioQQQ," Is the simulation physically motivated?\n");
			if( hextra.cryden==0. && lgPrint )
			{
				fprintf(ioQQQ,"PROBLEM h1fnd - no cosmic rays are present - is this physical?\n");
				fprintf(ioQQQ,"PROBLEM h1fnd - Consider including the COSMIC RAY BACKGROUND command.\n");
			}
			lgPrint = false;
		}
	}

	/* copy new solutions over the hmi.Hmolec */
	for(mol=0;mol<N_H_MOLEC;mol++)
	{
		hmi.Hmolec[mol] = (realnum) bvec[mol];
		ASSERT( hmi.Hmolec[mol] < MAX_DENSITY );
	}

	dense.xIonDense[ipHYDROGEN][0] = (realnum) bvec[ipMH];
	dense.xIonDense[ipHYDROGEN][1] = (realnum) bvec[ipMHp];

	/* total H2 - all forms */
	hmi.H2_total = hmi.Hmolec[ipMH2s] + hmi.Hmolec[ipMH2g];
	/* first guess at ortho and para densities */
	h2.ortho_density = 0.75*hmi.H2_total;
	h2.para_density = 0.25*hmi.H2_total;

	/* NB the first index must be kept parallel with nelem and ionstag in
	 * H2Lines transition struc,
	 * since that struc expects to find the abundances here */
	/* >>chng 04 feb 19, had been ipMH2g, chng to total */
	//dense.xIonDense[LIMELM+2][0] = hmi.H2_total;

	/* identify dominant H2 formation process */
	{
		/* following should be set true to identify H2 formation and destruction processes */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC && (nzone>50) )
		{
			double createsum ,create_from_Hn2 , create_3body_Ho, create_h2p, 
				create_h3p, create_h3pe, create_grains, create_hminus;
			double destroysum, destroy_hm ,destroy_soloman ,destroy_2h ,destroy_hp,
				destroy_h,destroy_hp2,destroy_h3p;

			/* H(n=2) + H(n=1) -> H2 */
			create_from_Hn2 = hmi.radasc*dense.xIonDense[ipHYDROGEN][0];
			/* 3H => H2 + H */
			create_3body_Ho = hmi.bh2dis*dense.xIonDense[ipHYDROGEN][0];
			/* H2+ + H => H2 + H+ */
			create_h2p = hmi.bh2h2p*dense.xIonDense[ipHYDROGEN][0]*Hmolec_old[ipMH2p];
			/* H + H3+ => H2 + H2+ */
			create_h3p = hmi.h3ph2p*dense.xIonDense[ipHYDROGEN][0]*hmi.Hmolec[ipMH3p];
			/* e + H3+ => H2 + H */
			create_h3pe = hmi.eh3_h2h*dense.eden * hmi.Hmolec[ipMH3p];
			/* from grains */
			create_grains = gv.rate_h2_form_grains_used_total;
			/* from H- */
			create_hminus = Hmolec_old[ipMH]*hmi.assoc_detach*hmi.Hmolec[ipMHm];

			createsum = create_from_Hn2 + create_3body_Ho + create_h2p +
				create_h3p + create_h3pe + create_grains + create_hminus;

			fprintf(ioQQQ,"H2 create zone\t%.2f \tsum\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
				fnzone,
				createsum,
				create_hminus   / createsum,
				create_from_Hn2 / createsum,
				create_3body_Ho / createsum,
				create_h2p      / createsum,
				create_h3p      / createsum,
				create_h3pe     / createsum,
				create_grains   / createsum	);

			/* all h2 molecular hydrogen destruction processes */
			/* >>chng 04 jan 28, had wrong Boltzmann factor,
			 * caught by gargi Shaw */
			destroy_hm = hmi.assoc_detach_backwards_grnd+hmi.assoc_detach_backwards_exct;
			/*destroy_hm2 = eh2hhm;*/
			destroy_soloman = hmi.H2_Solomon_dissoc_rate_used_H2g;
			destroy_2h = hmi.eh2hh;
			destroy_hp = hmi.h2hph3p*dense.xIonDense[ipHYDROGEN][1];
			destroy_h = hmi.rh2dis*dense.xIonDense[ipHYDROGEN][0];
			destroy_hp2 = hmi.rh2h2p*dense.xIonDense[ipHYDROGEN][1];
			destroy_h3p = hmi.h3petc * hmi.Hmolec[ipMH3p];
			destroysum = destroy_hm + /*destroy_hm2 +*/ destroy_soloman + destroy_2h + 
				destroy_hp+ destroy_h+ destroy_hp2+ destroy_h3p;

			fprintf(ioQQQ,"H2 destroy\t%.3f \t%.2e\tsum\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
				fnzone,
				destroysum,
				destroy_hm / destroysum ,
				destroy_soloman / destroysum ,
				destroy_2h / destroysum ,
				destroy_hp / destroysum ,
				destroy_h / destroysum ,
				destroy_hp2 / destroysum ,
				destroy_h3p / destroysum );

		}
	}

	{
		/* following should be set true to identify H- formation and destruction processes */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC && (nzone>140)/**/ )
		{
			double create_from_Ho,create_3body_Ho,create_batach,destroy_photo,
				destroy_coll_heavies,destroy_coll_electrons,destroy_Hattach,destroy_fhneut,
				destsum , createsum;

			create_from_Ho = (hmi.hminus_rad_attach + hmi.HMinus_induc_rec_rate);
			create_3body_Ho = c3bod;
			/* total formation is sum of g and s attachment */
			create_batach = hmi.assoc_detach_backwards_grnd + hmi.assoc_detach_backwards_exct;
			destroy_photo = hmi.HMinus_photo_rate;
			destroy_coll_heavies = hmi.hmin_ct_firstions*sum_first_ions;
			destroy_coll_electrons = cionhm;
			destroy_Hattach = Hmolec_old[ipMH]*hmi.assoc_detach;
			destroy_fhneut = fhneut;

			destsum = destroy_photo + destroy_coll_heavies + destroy_coll_electrons + 
				destroy_Hattach + destroy_fhneut;
			fprintf(ioQQQ,"H- destroy zone\t%.2f\tTe\t%.4e\tsum\t%.2e\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", 
				fnzone,
				phycon.te,
				destsum,
				destroy_photo/destsum , 
				destroy_coll_heavies/destsum,
				destroy_coll_electrons/destsum, 
				destroy_Hattach/destsum,
				destroy_fhneut/destsum );

			createsum = create_from_Ho+create_3body_Ho+create_batach;
			fprintf(ioQQQ,"H- create\t%.2f\tTe\t%.4e\tsum\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
				fnzone,
				phycon.te,
				createsum,
				dense.eden,
				create_from_Ho/createsum,
				create_3body_Ho/createsum,
				create_batach/createsum);
		}
	}

	/* rate H-alpha is created by H- ct */
	hmi.HalphaHmin = (realnum)(fhneut*hmi.Hmolec[ipMHm]);

	/* heating due to H2 dissociation */
	if( hmi.lgNoH2Mole )
	{
		hmi.HeatH2Dish_TH85 = 0.;
		hmi.HeatH2Dexc_TH85 = 0.;
		hmi.deriv_HeatH2Dexc_TH85 = 0.;
	}
	else
	{
		/* H2 photodissociation heating, eqn A9 of Tielens & Hollenbach 1985a */
		/*hmi.HeatH2Dish_TH85 = (realnum)(1.36e-23*hmi.Hmolec[ipMH2g]*h2esc*hmi.UV_Cont_rel2_Habing_TH85_depth);*/
		/* >>chng 04 feb 07, more general to express heating in terms of the assumed
		 * photo rates - the 0.25 was obtained by inverting A8 & A9 of TH85 to find that
		 * there are 0.25 eV per dissociative pumping, ie, 10% of total 
		 * this includes both H2g and H2s - TH85 say just ground but they include
		 * process for both H2 and H2s - as we did above - both must be in
		 * heating term */
		/* >>chng 05 mar 11, TE, old had used H2_Solomon_dissoc_rate_used, which was only
		 * H2g.  in regions where Solomon process is fast, H2s has a large population
		 * and the heating rate was underestimated. */
		/* >>chng 05 jun 23, 
		 * >>chng 05 dec 05, TE, modified to approximate the heating better for the
		 * new approximation */
		/* >>chng 00 nov 25, explicitly break out this heat agent */
		/* 2.6 eV of heat per deexcitation, consider difference
		* between deexcitation (heat) and excitation (cooling) */
		/* >>chng 04 jan 27, code moved here and simplified */
		/* >>chng 05 jul 10, GS*/ 
		/*  average collisional rate for H2* to H2g calculated from big H2, GS*/

		/* TH85 dissociation heating - this is ALWAYS defined for reference
		 * may be output for comparison with other rates*/
		hmi.HeatH2Dish_TH85 = 0.25 * EN1EV *
			(hmi.H2_Solomon_dissoc_rate_used_H2g * hmi.Hmolec[ipMH2g] +
			 hmi.H2_Solomon_dissoc_rate_used_H2s * hmi.Hmolec[ipMH2s]);

		/* TH85 deexcitation heating */
		hmi.HeatH2Dexc_TH85 = (hmi.Hmolec[ipMH2s]*H2star_deexcit - hmi.Hmolec[ipMH2g]*H2star_excit) * 4.17e-12;
		/* this is derivative wrt temperature, only if counted as a cooling source */
		hmi.deriv_HeatH2Dexc_TH85 = (realnum)(MAX2(0.,-hmi.HeatH2Dexc_TH85)* ( 30172. * thermal.tsq1 - thermal.halfte ) );

		if( hmi.chH2_small_model_type == 'H' )
		{
			/* Burton et al. 1990 */
			hmi.HeatH2Dish_BHT90 = 0.25 * EN1EV *
				(hmi.H2_Solomon_dissoc_rate_used_H2g * hmi.Hmolec[ipMH2g] +
				hmi.H2_Solomon_dissoc_rate_used_H2s * hmi.Hmolec[ipMH2s]);

			/* Burton et al. 1990 heating */
			hmi.HeatH2Dexc_BHT90 = (hmi.Hmolec[ipMH2s]*H2star_deexcit - hmi.Hmolec[ipMH2g]*H2star_excit) * 4.17e-12;
			/* this is derivative wrt temperature, only if counted as a cooling source */
			hmi.deriv_HeatH2Dexc_BHT90 = (realnum)(MAX2(0.,-hmi.HeatH2Dexc_BHT90)* ( 30172. * thermal.tsq1 - thermal.halfte ) );
		}
		else if( hmi.chH2_small_model_type == 'B')
		{
			/* Bertoldi & Draine */
			hmi.HeatH2Dish_BD96 = 0.25 * EN1EV *
				(hmi.H2_Solomon_dissoc_rate_used_H2g * hmi.Hmolec[ipMH2g] +
				hmi.H2_Solomon_dissoc_rate_used_H2s * hmi.Hmolec[ipMH2s]);
			/* Bertoldi & Draine heating, same as TH85 */
			hmi.HeatH2Dexc_BD96 = (hmi.Hmolec[ipMH2s]*H2star_deexcit - hmi.Hmolec[ipMH2g]*H2star_excit) * 4.17e-12;
			/* this is derivative wrt temperature, only if counted as a cooling source */
			hmi.deriv_HeatH2Dexc_BD96 = (realnum)(MAX2(0.,-hmi.HeatH2Dexc_BD96)* ( 30172. * thermal.tsq1 - thermal.halfte ) );
		}
		else if(hmi.chH2_small_model_type == 'E')
		{
			/* heating due to dissociation of H2
			 * >>chng 05 oct 19, TE, define new approximation for the heating due to the destruction of H2
			 *	use this approximation for the specified cloud parameters, otherwise
			 * use limiting cases for 1 <= G0, G0 >= 1e7, n >= 1e7, n <= 1 */

			double log_density, 
					f1, f2,f3, f4, f5;
			static double log_G0_face = -1;
			static double k_f4;


			/* test for G0 
			 * this is a constant so only do it in zone 0 */
			if( !nzone )
			{
				if(hmi.UV_Cont_rel2_Draine_DB96_face <= 1.) 
				{ 
					log_G0_face = 0.;
				}
				else if(hmi.UV_Cont_rel2_Draine_DB96_face >= 1e7) 
				{ 
					log_G0_face = 7.;
				}
				else 
				{ 
					log_G0_face = log10(hmi.UV_Cont_rel2_Draine_DB96_face); 
				}
				/*>>chng 06 oct 24 TE change Go face for spherical geometry*/
				log_G0_face /= radius.r1r0sq;
			}
			/* test for density */
			if(dense.gas_phase[ipHYDROGEN] <= 1.) 
			{ 
				log_density = 0.; 
			}
			else if(dense.gas_phase[ipHYDROGEN] >= 1e7) 
			{ 
				log_density = 7.; 
			}
			else 
			{ 
				log_density = log10(dense.gas_phase[ipHYDROGEN]); 
			}

			f1 = 0.15 * log_density + 0.75;
			f2 = -0.5 * log_density + 10.;

			hmi.HeatH2Dish_ELWERT = 0.25 * EN1EV *  f1 * 
				(hmi.H2_Solomon_dissoc_rate_used_H2g * hmi.Hmolec[ipMH2g] +
				 hmi.H2_Solomon_dissoc_rate_used_H2s * hmi.Hmolec[ipMH2s] ) + 
				f2 * secondaries.x12tot * EN1EV * hmi.H2_total;

			/*fprintf( ioQQQ, "f1: %.2e, f2: %.2e,heat Solomon: %.2e",f1,f2,hmi.HeatH2Dish_TH85);*/


			/* heating due to collisional deexcitation within X of H2
			* >>chng 05 oct 19, TE, define new approximation for the heating due to the destruction of H2
			*	use this approximation for the specified cloud parameters, otherwise
			* use limiting cases for 1 <= G0, G0 >= 1e7, n >= 1e7, n <= 1 */

			/* set value of k_f4 by testing on value of G0 */
			if(hmi.UV_Cont_rel2_Draine_DB96_face <= 1.) 
			{ 
				log_G0_face = 0.;
			}
			else if(hmi.UV_Cont_rel2_Draine_DB96_face >= 1e7) 
			{ 
				log_G0_face = 7.;
			}
			else 
			{ 
				log_G0_face = log10(hmi.UV_Cont_rel2_Draine_DB96_face); 
			}
			/* 06 oct 24, TE introduce effects of spherical geometry */
			log_G0_face /= radius.r1r0sq;

			/* terms only dependent on G0_face */
			k_f4 = -0.25 * log_G0_face + 1.25;

			/* test for density */
			if(dense.gas_phase[ipHYDROGEN] <= 1.) 
			{ 
				log_density = 0.; 
				f4 = 0.; 
			}
			else if(dense.gas_phase[ipHYDROGEN] >= 1e7) 
			{ 
				log_density = 7.; 
				f4 = pow2(k_f4) * pow( 10. , 2.2211 * log_density - 29.8506);
			}
			else 
			{ 
				log_density = log10(dense.gas_phase[ipHYDROGEN]); 
				f4 = pow2(k_f4) * pow( 10., 2.2211 * log_density - 29.8506);
			}

			f3 = MAX2(0.1, -4.75 * log_density + 24.25);
			f5 = MAX2(1.,0.95 * log_density - 1.45) * 0.2 * log_G0_face;

			hmi.HeatH2Dexc_ELWERT = (hmi.Hmolec[ipMH2s]*H2star_deexcit - hmi.Hmolec[ipMH2g]*H2star_excit) * 4.17e-12 * f3 + 
				f4 * (hmi.Hmolec[ipMH]/dense.gas_phase[ipHYDROGEN]) +
				f5 * secondaries.x12tot * EN1EV * hmi.H2_total;

			if(log_G0_face == 0.&& dense.gas_phase[ipHYDROGEN] > 1.) 
				hmi.HeatH2Dexc_ELWERT *= 0.001 / dense.gas_phase[ipHYDROGEN];

			/* >>chng 06 oct 24, TE introduce effects of spherical geometry */
			/*if(radius.depth/radius.rinner >= 1.0) */
			hmi.HeatH2Dexc_ELWERT /= POW2(radius.r1r0sq);

			/* this is derivative wrt temperature, only if counted as a cooling source */
			hmi.deriv_HeatH2Dexc_ELWERT = (realnum)(MAX2(0.,-hmi.HeatH2Dexc_ELWERT)* ( 30172. * thermal.tsq1 - thermal.halfte ) );

			/*fprintf( ioQQQ, "\tf3: %.2e, f4: %.2e, f5: %.2e, heat coll dissoc: %.2e\n",f3,f4,f5,hmi.HeatH2Dexc_TH85);*/
		}
		/* end Elwert branch for photo rates */
		else
			TotalInsanity();

		if( h2.lgH2ON  && hmi.lgBigH2_evaluated && hmi.lgH2_Chemistry_BigH2 )
		{
				deexc_htwo = hmi.Average_collH2_deexcit;
				deexc_hneut = hmi.Average_collH_deexcit;
		}
		else
		{
			deexc_htwo = (1.4e-12*phycon.sqrte * sexp( 18100./(phycon.te + 1200.) ))/6.;
			deexc_hneut =  (1e-12*phycon.sqrte * sexp(1000./phycon.te ))/6.;
		}

		H2star_deexcit = hmi.H2_total*deexc_htwo + hmi.Hmolec[ipMH] * deexc_hneut;

		if( h2.lgH2ON  && hmi.lgBigH2_evaluated && hmi.lgH2_Chemistry_BigH2 )
		{
			H2star_excit = hmi.Average_collH2_excit *hmi.H2_total + hmi.Average_collH_excit*hmi.Hmolec[ipMH];
		}
		else
		{
			H2star_excit = Boltz_fac_H2_H2star * H2star_deexcit;
		}

		/* Leiden hacks try to turn off H2*, which is all unphysical.  do not include heating
		 * due to H2 deexcitation since H2* is bogus */
		/* >>chng 05 aug 12, do not turn off vibrational heating when Leiden hack is in place
		 * other codes included heating but did not include H2s on chemistry 
		hmi.HeatH2Dexc_TH85 *= hmi.lgLeiden_Keep_ipMH2s;*/
		/*fprintf(ioQQQ,
			"DEBUG hmole H2 deex heating:\t%.2f\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
			fnzone,
			hmi.HeatH2Dexc_TH85,
			thermal.htot,
			hmi.Hmolec[ipMH2s],
			H2star_deexcit, 
			hmi.Hmolec[ipMH2g],
			H2star_excit,
			Hmolec_old[ipMH2g],
			Hmolec_old[ipMH],
			phycon.te);*/
	}

	{
		/* following should be set true to print populations */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			/* these are the raw results */
			fprintf( ioQQQ, " HMOLE raw; hi\t%.2e" , dense.xIonDense[ipHYDROGEN][0]);
			for( i=0; i < N_H_MOLEC; i++ )
			{
				fprintf( ioQQQ, "\t%s\t%.2e", hmi.chLab[i], bvec[i] );
			}
			fprintf( ioQQQ, " \n" );
		}
	}

	if( trace.lgTrace && trace.lgTr_H2_Mole )
	{
		/* these are the raw results */
		fprintf( ioQQQ, "\n raw; " );
		for( i=0; i < N_H_MOLEC; i++ )
		{
			fprintf( ioQQQ, " %s:%.2e", hmi.chLab[i], bvec[i] );
		}
		fprintf( ioQQQ, " \n" );
	}

	/* >>chng 05 jul 11, TE, each term must have unit cm-3 s-1,
	 * more terms added */
	hmi.H2_rate_create = gv.rate_h2_form_grains_used_total  * dense.xIonDense[ipHYDROGEN][0] + 
		hmi.assoc_detach * hmi.Hmolec[ipMH] * hmi.Hmolec[ipMHm] +
		hmi.bh2dis * dense.xIonDense[ipHYDROGEN][0] + 
		hmi.bh2h2p * dense.xIonDense[ipHYDROGEN][0] * hmi.Hmolec[ipMH2p] + 
		hmi.radasc * dense.xIonDense[ipHYDROGEN][0] + 
		hmi.h3ph2p * dense.xIonDense[ipHYDROGEN][0] * hmi.Hmolec[ipMH3p] + 
		hmi.h2phmh2h * hmi.Hmolec[ipMH2p] * hmi.Hmolec[ipMHm] +
		hmi.bh2h22hh2 * 2 * dense.xIonDense[ipHYDROGEN][0] * hmi.Hmolec[ipMH2g] +
		hmi.h3phmh2hh * hmi.Hmolec[ipMH3p] * hmi.Hmolec[ipMHm] +
		hmi.h3phm2h2 * hmi.Hmolec[ipMH3p] * hmi.Hmolec[ipMHm] +
		hmi.h32h2 * hmi.Hmolec[ipMH2p] * hmi.Hmolec[ipMH3p] +
		hmi.eh3_h2h * hmi.Hmolec[ipMH3p] +
		hmi.h3ph2hp * hmi.Hmolec[ipMH3p]+
		co.H_CH_C_H2 * dense.xIonDense[ipHYDROGEN][0] +
		co.H_CHP_CP_H2 * dense.xIonDense[ipHYDROGEN][0] +
		co.H_CH2_CH_H2 * dense.xIonDense[ipHYDROGEN][0] +
		co.H_CH3P_CH2P_H2 * dense.xIonDense[ipHYDROGEN][0] +
		co.H_OH_O_H2 * dense.xIonDense[ipHYDROGEN][0] +
		co.Hminus_HCOP_CO_H2 * hmi.Hmolec[ipMHm] +
		co.Hminus_H3OP_H2O_H2 * hmi.Hmolec[ipMHm] +
		co.Hminus_H3OP_OH_H2_H * hmi.Hmolec[ipMHm] +
		co.HP_CH2_CHP_H2 * hmi.Hmolec[ipMHp] +
		co.HP_SiH_SiP_H2 * hmi.Hmolec[ipMHp] +
		co.H2P_CH_CHP_H2 * hmi.Hmolec[ipMH2p] + 
		co.H2P_CH2_CH2P_H2 * hmi.Hmolec[ipMH2p] + 
		co.H2P_CO_COP_H2 * hmi.Hmolec[ipMH2p] + 
		co.H2P_H2O_H2OP_H2 * hmi.Hmolec[ipMH2p] + 
		co.H2P_O2_O2P_H2 * hmi.Hmolec[ipMH2p] + 
		co.H2P_OH_OHP_H2 * hmi.Hmolec[ipMH2p] + 
		co.H3P_C_CHP_H2 * hmi.Hmolec[ipMH3p] + 
		co.H3P_CH_CH2P_H2 * hmi.Hmolec[ipMH3p] + 
		co.H3P_CH2_CH3P_H2 * hmi.Hmolec[ipMH3p] + 
		co.H3P_OH_H2OP_H2 * hmi.Hmolec[ipMH3p] + 
		co.H3P_H2O_H3OP_H2 * hmi.Hmolec[ipMH3p] + 
		co.H3P_CO_HCOP_H2 * hmi.Hmolec[ipMH3p] + 
		co.H3P_O_OHP_H2 * hmi.Hmolec[ipMH3p] + 
		co.H3P_SiH_SiH2P_H2 * hmi.Hmolec[ipMH3p] + 
		co.H3P_SiO_SiOHP_H2	* hmi.Hmolec[ipMH3p] +
		co.H_CH3_CH2_H2 * dense.xIonDense[ipHYDROGEN][0] +
		co.H_CH4P_CH3P_H2 * dense.xIonDense[ipHYDROGEN][0] +
		co.H_CH5P_CH4P_H2 * dense.xIonDense[ipHYDROGEN][0] +
		co.H2P_CH4_CH3P_H2 * hmi.Hmolec[ipMH2p] + 
		co.H2P_CH4_CH4P_H2 * hmi.Hmolec[ipMH2p] + 
		co.H3P_CH3_CH4P_H2 * hmi.Hmolec[ipMH3p] + 
		co.H3P_CH4_CH5P_H2 * hmi.Hmolec[ipMH3p] + 
		co.HP_CH4_CH3P_H2 * hmi.Hmolec[ipMHp] +	
		co.HP_HNO_NOP_H2 * hmi.Hmolec[ipMHp] +
		co.HP_HS_SP_H2 * hmi.Hmolec[ipMHp] +
		co.H_HSP_SP_H2 * dense.xIonDense[ipHYDROGEN][0] +
		co.H3P_NH_NH2P_H2 * hmi.Hmolec[ipMH3p] + 
		co.H3P_NH2_NH3P_H2 * hmi.Hmolec[ipMH3p] + 
		co.H3P_NH3_NH4P_H2 * hmi.Hmolec[ipMH3p] + 
		co.H3P_CN_HCNP_H2 * hmi.Hmolec[ipMH3p] + 
		co.H3P_NO_HNOP_H2 * hmi.Hmolec[ipMH3p] + 
		co.H3P_S_HSP_H2 * hmi.Hmolec[ipMH3p] + 
		co.H3P_CS_HCSP_H2 * hmi.Hmolec[ipMH3p] + 
		co.H3P_NO2_NOP_OH_H2 * hmi.Hmolec[ipMH3p] + 
		co.H2P_NH_NHP_H2 * hmi.Hmolec[ipMH2p] + 
		co.H2P_NH2_NH2P_H2 * hmi.Hmolec[ipMH2p] + 
		co.H2P_NH3_NH3P_H2 * hmi.Hmolec[ipMH2p] + 
		co.H2P_CN_CNP_H2 * hmi.Hmolec[ipMH2p] + 
		co.H2P_HCN_HCNP_H2 * hmi.Hmolec[ipMH2p] + 
		co.H2P_NO_NOP_H2 * hmi.Hmolec[ipMH2p] +
		co.H3P_Cl_HClP_H2 * hmi.Hmolec[ipMH3p]+
		co.H3P_HCl_H2ClP_H2 * hmi.Hmolec[ipMH3p]+
		co.H2P_C2_C2P_H2 * hmi.Hmolec[ipMH2p]+	
		co.Hminus_NH4P_NH3_H2 * hmi.Hmolec[ipMHm]+
		co.H3P_HCN_HCNHP_H2 * hmi.Hmolec[ipMH3p];

	/* option to print rate H2 forms */
	/* trace.lgTr_H2_Mole is trace molecules option,
	 * save htwo */
	if( (trace.lgTrace && trace.lgTr_H2_Mole) )
	{

		if( hmi.H2_rate_create > SMALLFLOAT )
		{
			fprintf( ioQQQ, 
				" Create H2, rate=%10.2e grain;%5.3f hmin;%5.3f bhedis;%5.3f h2+;%5.3f hmi.radasc:%5.3f hmi.h3ph2p:%5.3f hmi.h3petc:%5.3f\n", 
			  hmi.H2_rate_create, 
			  gv.rate_h2_form_grains_used_total/hmi.H2_rate_create, 
			  hmi.assoc_detach * hmi.Hmolec[ipMH] * hmi.Hmolec[ipMHm] /hmi.H2_rate_create, 
			  hmi.bh2dis * dense.xIonDense[ipHYDROGEN][0]/hmi.H2_rate_create, 
			  hmi.bh2h2p*dense.xIonDense[ipHYDROGEN][0]*hmi.Hmolec[ipMH2p]/hmi.H2_rate_create, 
			  hmi.radasc*dense.xIonDense[ipHYDROGEN][0]/hmi.H2_rate_create, 
			  hmi.h3ph2p*hmi.Hmolec[ipMH3p]/hmi.H2_rate_create, 
			  hmi.h3petc*hmi.Hmolec[ipMH3p]/hmi.H2_rate_create );
		}
		else
		{
			fprintf( ioQQQ, " Create H2, rate=0\n" );
		}
	}

	/* this is H2+ */
	if( trace.lgTrace && trace.lgTr_H2_Mole )
	{
		rate = hmi.rh2h2p*hmi.Hmolec[ipMH2g]*dense.xIonDense[ipHYDROGEN][1] + b2pcin*dense.xIonDense[ipHYDROGEN][1]*dense.eden*dense.xIonDense[ipHYDROGEN][0] + 
		  hmi.h3ph2p*hmi.Hmolec[ipMH3p] + hmi.h3petc*hmi.Hmolec[ipMH3p];
		if( rate > 1e-25 )
		{
			fprintf( ioQQQ, " Create H2+, rate=%10.2e hmi.rh2h2p;%5.3f b2pcin;%5.3f hmi.h3ph2p;%5.3f hmi.h3petc+;%5.3f\n", 
			  rate, hmi.rh2h2p*dense.xIonDense[ipHYDROGEN][1]*hmi.Hmolec[ipMH2g]/rate, 
				b2pcin*dense.xIonDense[ipHYDROGEN][1]*dense.xIonDense[ipHYDROGEN][1]*
			  dense.eden/rate, hmi.h3ph2p*hmi.Hmolec[ipMH3p]/rate, hmi.h3petc*hmi.Hmolec[ipMH3p]/
			  rate );
		}
		else
		{
			fprintf( ioQQQ, " Create H2+, rate=0\n" );
		}
	}

	double denom = (double)dense.xIonDense[ipHYDROGEN][0]*dense.eden*hmi.rel_pop_LTE_Hmin;
	if( hmi.Hmolec[ipMHm] > 0. && denom > 0. )
		hmi.hmidep = (double)hmi.Hmolec[ipMHm]/ denom;
	else
		hmi.hmidep = 1.;

	/* this will be net volume heating rate, photo heat - induc cool */
	hmi.hmihet = hmi.HMinus_photo_heat*hmi.Hmolec[ipMHm] - hmi.HMinus_induc_rec_cooling;
	hmi.h2plus_heat = h2phet*hmi.Hmolec[ipMH2p];

	/* departure coefficient for H2 defined rel to n(1s)**2
	 * (see Littes and Mihalas Solar Phys 93, 23) */
	plte = (double)dense.xIonDense[ipHYDROGEN][0] * hmi.rel_pop_LTE_H2g * (double)dense.xIonDense[ipHYDROGEN][0];
	if( plte > 0. )
	{
		hmi.h2dep = hmi.Hmolec[ipMH2g]/plte;
	}
	else
	{
		hmi.h2dep = 1.;
	}

	/* departure coefficient of H2+ defined rel to n(1s) n(p)
	 * sec den was HI before 85.34 */
	plte = (double)dense.xIonDense[ipHYDROGEN][0]*hmi.rel_pop_LTE_H2p*(double)dense.xIonDense[ipHYDROGEN][1];
	if( plte > 0. )
	{
		hmi.h2pdep = hmi.Hmolec[ipMH2p]/plte;
	}
	else
	{
		hmi.h2pdep = 1.;
	}

	/* departure coefficient of H3+ defined rel to N(H2+) N(p) */
	if( hmi.rel_pop_LTE_H3p > 0. )
	{
		hmi.h3pdep = hmi.Hmolec[ipMH3p]/hmi.rel_pop_LTE_H3p;
	}
	else
	{
		hmi.h3pdep = 1.;
	}


	if( trace.lgTrace && trace.lgTr_H2_Mole )
	{
		fprintf( ioQQQ, " HMOLE, Dep Coef, H-:%10.2e H2:%10.2e H2+:%10.2e\n", 
		  hmi.hmidep, hmi.h2dep, hmi.h2pdep );
		fprintf( ioQQQ, "     H- creat: Rad atch%10.3e Induc%10.3e bHneut%10.2e 3bod%10.2e b=>H2%10.2e N(H-);%10.2e b(H-);%10.2e\n", 
		  hmi.hminus_rad_attach, hmi.HMinus_induc_rec_rate, bhneut, c3bod, hmi.assoc_detach_backwards_grnd, hmi.Hmolec[ipMHm], hmi.hmidep );

		fprintf( ioQQQ, "     H- destr: Photo;%10.3e mut neut%10.2e e- coll ion%10.2e =>H2%10.2e x-ray%10.2e p+H-%10.2e\n", 
		  hmi.HMinus_photo_rate, hmi.hmin_ct_firstions*sum_first_ions, cionhm, 
		  Hmolec_old[ipMH]*hmi.assoc_detach, iso.gamnc[ipH_LIKE][ipHYDROGEN][ipH1s], 
		  fhneut );
		fprintf( ioQQQ, "     H- heating:%10.3e Ind cooling%10.2e Spon cooling%10.2e\n", 
		  hmi.hmihet, hmi.HMinus_induc_rec_cooling, hmi.hmicol );
	}

	/* identify creation and destruction processes for H2+ */
	if( trace.lgTrace && trace.lgTr_H2_Mole )
	{
		rate = desh2p;
		if( rate != 0. )
		{
			fprintf( ioQQQ, 
				" Destroy H2+: rate=%10.2e e-;%5.3f phot;%5.3f hard gam;%5.3f H2col;%5.3f h2phhp;%5.3f pion;%5.3f bh2h2p:%5.3f\n", 
			  rate, h2pcin*dense.eden/rate, gamtwo/rate, 2.*iso.gamnc[ipH_LIKE][ipHYDROGEN][ipH1s]/
			  rate, hmi.h2ph3p/rate, h2phhp/rate, h2pion/rate, hmi.bh2h2p*
			  dense.xIonDense[ipHYDROGEN][0]/rate );

			rate *= hmi.Hmolec[ipMH2p];
			if( rate > 0. )
			{
				fprintf( ioQQQ, 
					" Create  H2+: rate=%.2e HII HI;%.3f Col H2;%.3f HII H2;%.3f HI HI;%.3f\n", 
				  rate, 
				  radath*dense.xIonDense[ipHYDROGEN][1]*dense.xIonDense[ipHYDROGEN][0]/rate, 
				  (hmi.H2_photoionize_rate + secondaries.csupra[ipHYDROGEN][0]*2.02)*hmi.Hmolec[ipMH2g]/rate, 
				  hmi.rh2h2p*dense.xIonDense[ipHYDROGEN][1]*hmi.Hmolec[ipMH2g]/rate, 
				  b2pcin*dense.xIonDense[ipHYDROGEN][0]*dense.xIonDense[ipHYDROGEN][1]*dense.eden/rate );
			}
			else
			{
				fprintf( ioQQQ, " Create  H2+: rate= is zero\n" );
			}
		}
	}

	{
		/* following should be set true to print populations */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			fprintf(ioQQQ,"hmole bugg\t%.3f\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
				fnzone,
				iso.gamnc[ipH_LIKE][ipHYDROGEN][0],
				hmi.HMinus_photo_rate,
				hmi.Hmolec[ipMH2g] , 
				hmi.Hmolec[ipMHm] ,
				dense.xIonDense[ipHYDROGEN][1]);
		}
	}
	return;
}
#if defined(__HP_aCC)
#pragma OPTIMIZE OFF
#pragma OPTIMIZE ON
#endif
/*lint +e778 const express eval to 0 */
/*lint +e725 expect positive indentation */
