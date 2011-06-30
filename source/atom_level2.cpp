/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*atom_level2 do level population and cooling for two level atom,
 * side effects:
 * set elements of transition struc
 * cooling via 	CoolAdd( chLab, (long)t->WLAng , t->cool);
 * cooling derivative */
#include "cddefines.h"
#include "phycon.h"
#include "transition.h"
#include "dense.h"
#include "rfield.h"
#include "thermal.h"
#include "cooling.h"
#include "atoms.h"

void atom_level2(transition * t)
{
	char chLab[5];
	long int ion, 
	  ip, 
	  nelem;
	double AbunxIon, 
	  a21, 
	  boltz, 
	  col12, 
	  col21, 
	  coolng, 
	  g1, 
	  g2, 
	  omega, 
	  pfs1,
	  pfs2,
	  r, 
	  rate12, 
	  ri21;

	DEBUG_ENTRY( "atom_level2()" );

	/* result is density (cm-3) of excited state times A21
	 * result normalized to N1+N2=ABUND
	 * routine increments dCooldT, call CoolAdd
	 * CDSQTE is EDEN / SQRTE * 8.629E-6
	 */

	ion = t->Hi->IonStg;
	nelem = t->Hi->nelem;

	/* dense.xIonDense[nelem][i] is density of ith ionization stage (cm^-3) */
	AbunxIon = dense.xIonDense[nelem-1][ion-1];

	/* continuum pointer */
	ip = t->ipCont;

	/* approximate Boltzmann factor to see if results zero */
	boltz = rfield.ContBoltz[ip-1];

	/* collision strength for this transition, omega is zero for hydrogenic
	 * species which are done in special hydro routines */
	omega = t->Coll.col_str;

	/* ROUGH check whether upper level has significant population,*/
	r = (boltz*dense.cdsqte + t->Emis->pump)/(dense.cdsqte + t->Emis->Aul);

	/* following first test needed for 32 bit cpu on search phase
	 * >>chng 96 jul 02, was 1e-30 crashed on dec, change to 1e-25
	 * if( AbunxIon.lt.1e-25 .or. boltz.gt.30. ) then
	 * >>chng 96 jul 11, to below, since can be strong pumping when
	 * Boltzmann factor essentially zero */
	/* omega in following is zero for hydrogenic species, since done
	 * in hydro routines, so this should cause us to quit on this test */
	/* >>chng 99 nov 29, from 1e-25 to 1e-30 to keep same result for
	 * very low density models, where AbunxIon is very small but still significant*/
	/*if( omega*AbunxIon < 1e-25 || r < 1e-25 )*/
	if( omega*AbunxIon < 1e-30 || r < 1e-25 )
	{
		/* put in pop since possible just too cool */
		t->Lo->Pop = AbunxIon;
		t->Emis->PopOpc = AbunxIon;
		t->Hi->Pop = 0.;
		t->Emis->xIntensity = 0.;
		t->Coll.cool = 0.;
		t->Emis->ots = 0.;
		t->Emis->phots = 0.;
		t->Emis->ColOvTot = 0.;
		t->Coll.heat = 0.;
		/* level populations */
		atoms.PopLevels[0] = AbunxIon;
		atoms.PopLevels[1] = 0.;
		atoms.DepLTELevels[0] = 1.;
		atoms.DepLTELevels[1] = 0.;
		return;
	}

	/* net rate down A21*(escape + destruction) */
	a21 = t->Emis->Aul*(t->Emis->Pesc+ t->Emis->Pdest + t->Emis->Pelec_esc);

	/* put null terminated line label into chLab using line array*/
	chIonLbl(chLab,t);

	/* statistical weights of upper and lower levels */
	g1 = t->Lo->g;
	g2 = t->Hi->g;

	/* now get real Boltzmann factor */
	boltz = t->EnergyK/phycon.te;

	ASSERT( boltz > 0. );
	boltz = sexp(boltz);

	ASSERT( g1 > 0. && g2 > 0. );

	/* this lacks the upper statistical weight */
	col21 = dense.cdsqte*omega;
	/* upward coll rate s-1 */
	col12 = col21/g1*boltz;
	/* downward coll rate s-1 */
	col21 /= g2;

	/* rate 1 to 2 is both collisions and pumping */
	/* the total excitation rate from lower to upper, collisional and pumping */
	rate12 = col12 + t->Emis->pump;

	/* induced emissions down */
	ri21 = t->Emis->pump*g1/g2;

	/* this is the ratio of lower to upper level populations */
	r = (a21 + col21 + ri21)/rate12;

	/* upper level pop */
	pfs2 = AbunxIon/(r + 1.);
	atoms.PopLevels[1] = pfs2;
	t->Hi->Pop = pfs2;

	/* pop of ground */
	pfs1 = pfs2*r;
	atoms.PopLevels[0] = pfs1;

	/* compute ratio Aul/(Aul+Cul) */
	/* level population with correction for stimulated emission */
	t->Lo->Pop = atoms.PopLevels[0];

	t->Emis->PopOpc = (atoms.PopLevels[0] - atoms.PopLevels[1]*g1/g2 );

	/* departure coef of excited state rel to ground */
	atoms.DepLTELevels[0] = 1.;
	if( boltz > 1e-20 && atoms.PopLevels[1] > 1e-20 )
	{
		/* this line could have zero boltz factor but radiatively excited
		 * dec alpha does not obey () in fast mode?? */
		atoms.DepLTELevels[1] = (atoms.PopLevels[1]/atoms.PopLevels[0])/
		  (boltz*g2/g1);
	}
	else
	{
		atoms.DepLTELevels[1] = 0.;
	}

	/* number of escaping line photons, used elsewhere for outward beam */
	t->Emis->phots = t->Emis->Aul*(t->Emis->Pesc + t->Emis->Pelec_esc)*pfs2;

	/* intensity of line */
	t->Emis->xIntensity = t->Emis->phots*t->EnergyErg;
	//double plower = AbunxIon - pfs2;

	/* ratio of collisional to total (collisional + pumped) excitation */
	t->Emis->ColOvTot = col12/rate12;

	/* two cases - collisionally excited (usual case) or 
	 * radiatively excited - in which case line can be a heat source
	 * following are correct heat exchange, they will mix to get correct deriv 
	 * the sum of heat-cool will add up to EnrUL - EnrLU - this is a trick to
	 * keep stable solution by effectively dividing up heating and cooling,
	 * so that negative cooling does not occur */

	//double Enr12 = plower*col12*t->EnergyErg;
	//double Enr21 = pfs2*col21*t->EnergyErg;

	/* energy exchange due to this level
	 * net cooling due to excit minus part of de-excit -
	 * note that ColOvTot cancels out in the sum heat - cool */
	//coolng = Enr12 - Enr21*t->Emis->ColOvTot;

	/* this form of coolng is guaranteed to be non-negative */
	coolng = t->EnergyErg*AbunxIon*col12*(a21 + ri21)/(a21 + col21 + ri21 + rate12);
	ASSERT( coolng >= 0. );

	t->Coll.cool = coolng;

	/* net heating is remainder */
	t->Coll.heat = t->EnergyErg*AbunxIon*col21*(t->Emis->pump)/(a21 + col21 + ri21 + rate12);

	/* expression pre jul 3 95, changed for case where line heating dominates
	 * coolng = (plower*col12 - pfs2*col21)*t->t(ipLnEnrErg)
	 * t->t(ipLnCool) = cooling */

	/* add to cooling - heating part was taken out above,
	 * and is not added in here - it will be added to thermal.heating[0][22]
	 * in CoolSum */
	CoolAdd( chLab, t->WLAng , t->Coll.cool);

	/* derivative of cooling function */
	thermal.dCooldT += coolng * (t->EnergyK * thermal.tsq1 - thermal.halfte );
	return;
}
