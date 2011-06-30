/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*atom_level3 compute three level atom, 10, 21, and 20 are line */
#include "cddefines.h"
#include "phycon.h"
#include "physconst.h"
#include "dense.h"
#include "transition.h"
#include "thermal.h"
#include "cooling.h"
/* NB - use the nLev3Fail failures mode indicators when debugging this routine */
#include "atoms.h"
#include "rfield.h"

void atom_level3(transition * t10, 
  transition * t21, 
  transition * t20)
{
	char chLab[5], 
	  chLab10[11];
	double AbunxIon, 
	  a, 
	  a10, 
	  a20, 
	  a21, 
	  b, 
	  beta, 
	  bolt01, 
	  bolt02, 
	  bolt12, 
	  c, 
	  ener10, 
	  ener20, 
	  ener21, 
	  g0, 
	  g010, 
	  g020, 
	  g1, 
	  g110, 
	  g121, 
	  g2, 
	  g220, 
	  g221, 
	  o10, 
	  o20, 
	  o21, 
	  p0, 
	  p1, 
	  p2, 
	  pump01, 
	  pump02, 
	  pump12;

	int nLev3Fail;

	double TotCool, 
	  TotHeat, 
	  TotInten, 
	  alpha, 
	  alpha1, 
	  alpha2, 
	  c01, 
	  c02, 
	  c10, 
	  c12, 
	  c20, 
	  c21, 
	  cnet01, 
	  cnet02, 
	  cnet12, 
	  cool01, 
	  cool02, 
	  cool12, 
	  heat10, 
	  heat20, 
	  heat21, 
	  hnet01, 
	  hnet02, 
	  hnet12, 
	  pump10, 
	  pump20, 
	  pump21, 
	  r01, 
	  r02, 
	  r10, 
	  r12, 
	  r20, 
	  r21, 
	  temp01, 
	  temp02, 
	  temp12;

	DEBUG_ENTRY( "atom_level3()" );

	/* compute three level atom, 10, 21, and 20 are line
	 * arrays for 10, 21, and 20 transitions.
	 * one can be a dummy */
	/* >>chng 96 dec 06, to double precision due to round off problems below */

	/* generalized three level atom for any ion
	 * sum of three levels normalized to total abundance
	 *
	 * stat weights of all three lines
	 * sanity check will confirm ok */
	g010 = t10->Lo->g;
	g110 = t10->Hi->g;

	g121 = t21->Lo->g;
	g221 = t21->Hi->g;

	g020 = t20->Lo->g;
	g220 = t20->Hi->g;

	/* these are statistical weights */
	if( g010 > 0. )
	{
		g0 = g010;
	}

	else if( g020 > 0. )
	{
		g0 = g020;
	}

	else
	{
		g0 = -1.;
		strcpy( chLab10, chLineLbl(t10) );
		fprintf( ioQQQ, " PROBLEM atom_level3: insane stat weights g0 :%10.2e%10.2e %10.10s\n", 
		  g010, g020, chLab10 );
		TotalInsanity();
	}

	if( g110 > 0. )
	{
		g1 = g110;
	}

	else if( g121 > 0. )
	{
		g1 = g121;
	}

	else
	{
		g1 = -1.;
		strcpy( chLab10, chLineLbl(t10) );
		fprintf( ioQQQ, " PROBLEM atom_level3: insane stat weights g1 :%10.2e%10.2e %10.10s\n", 
		  g010, g020, chLab );
		TotalInsanity();
	}

	if( g220 > 0. )
	{
		g2 = g220;
	}
	else if( g221 > 0. )
	{
		g2 = g221;
	}

	else
	{
		g2 = -1.;
		strcpy( chLab10, chLineLbl(t20) );
		fprintf( ioQQQ, " PROBLEM atom_level3: insane stat weights g2 :%10.2e%10.2e %10.10s\n", 
		  g010, g020, chLab10 );
		TotalInsanity();
	}

	/* abundances from the elements grid
	 * one of these must be a true line */
	if( g010 > 0. )
	{
		/* put null terminated line label into chLab using line array*/
		chIonLbl(chLab, t10);
		AbunxIon = dense.xIonDense[ t10->Hi->nelem -1][t10->Hi->IonStg-1];
	}

	else if( g121 > 0. )
	{
		/* put null terminated line label into chLab using line array*/
		chIonLbl(chLab, t21);
		AbunxIon = dense.xIonDense[t21->Hi->nelem -1][t21->Hi->IonStg-1];
	}

	else
		/* this cannot possibly happen */
	{
		chLab[0] = ' ';
		AbunxIon = 0.;
		fprintf( ioQQQ, " PROBLEM atom_level3: insanity at g010 g121 branch \n" );
		TotalInsanity();
	}

	a = t10->EnergyK*phycon.teinv;
	b = t21->EnergyK*phycon.teinv;
	c = t20->EnergyK*phycon.teinv;

	if( c == 0. )
	{
		c = a + b;
	}

	/* if still neg at end, then success!, so possible to
	 * to check why zero returned */
	nLev3Fail = -1;

	/* if two of the lines are below the plasma frequency, hand the third in atom_level2 */
	if( ( t10->EnergyErg/EN1RYD < rfield.plsfrq && t20->EnergyErg/EN1RYD < rfield.plsfrq ) )
	{
		atom_level2( t21 );
		return;
	}
	else if( ( t10->EnergyErg/EN1RYD < rfield.plsfrq && t21->EnergyErg/EN1RYD < rfield.plsfrq ) )
	{
		atom_level2( t20 );
		return;
	}
	else if( ( t20->EnergyErg/EN1RYD < rfield.plsfrq && t21->EnergyErg/EN1RYD < rfield.plsfrq ) )
	{
		atom_level2( t10 );
		return;
	}

	/** \todo	2	test on c checks whether collisions are possible at this temperature,
	 * should add photo excitation */
	if( AbunxIon <= 1e-30 || c > 60. )
	{
		nLev3Fail = 0;

		/* all populations are zero */
		atoms.PopLevels[0] = AbunxIon;
		atoms.PopLevels[1] = 0.;
		atoms.PopLevels[2] = 0.;

		/** \todo	2	these pops ARE NOT defined below */
		atoms.DepLTELevels[0] = 1.;
		atoms.DepLTELevels[1] = 0.;
		atoms.DepLTELevels[2] = 0.;

		/* level populations */
		t21->Emis->PopOpc = 0.;
		t10->Emis->PopOpc = AbunxIon;
		t20->Emis->PopOpc = AbunxIon;
		t21->Lo->Pop = 0.;
		t10->Lo->Pop = AbunxIon;
		t20->Lo->Pop = AbunxIon;
		t21->Hi->Pop = 0.;
		t10->Hi->Pop = 0.;
		t20->Hi->Pop = 0.;

		/* line heating */
		t20->Coll.heat = 0.;
		t21->Coll.heat = 0.;
		t10->Coll.heat = 0.;

		/* intensity of line */
		t21->Emis->xIntensity = 0.;
		t10->Emis->xIntensity = 0.;
		t20->Emis->xIntensity = 0.;

		/* line cooling */
		t20->Coll.cool = 0.;
		t21->Coll.cool = 0.;
		t10->Coll.cool = 0.;

		/* local ots rates */
#		if 0
		/* >>chng 03 oct 04, move to RT_OTS */
		t20->ots = 0.;
		t21->ots = 0.;
		t10->ots = 0.;
#		endif

		/* number of photons in line zero */
		t21->Emis->phots = 0.;
		t10->Emis->phots = 0.;
		t20->Emis->phots = 0.;

		/* ratio coll over total excitation */
		t21->Emis->ColOvTot = 0.;
		t10->Emis->ColOvTot = 0.;
		t20->Emis->ColOvTot = 0.;

		/* add zero to cooling */
		CoolAdd(chLab, t21->WLAng ,0.);
		CoolAdd(chLab, t10->WLAng ,0.);
		CoolAdd(chLab, t20->WLAng ,0.);
		return;
	}

	/* collision strengths */
 	o10 = t10->Coll.col_str;
	o21 = t21->Coll.col_str;
	o20 = t20->Coll.col_str;

	/* begin sanity checks, check statistic weights, 
	 * first check is protection against dummy lines */
	ASSERT( g010*g020 == 0. || fp_equal( g010, g020 ) );

	ASSERT( g110*g121 == 0. || fp_equal( g110, g121 ) );

	ASSERT( g221*g220 == 0. || fp_equal( g221, g220 ) );

	/* both abundances must be same, equal abundance
	 * dense.xIonDense(nelem,i) is density of ith ionization stage (cm^-3) */
	ASSERT( (t10->Hi->IonStg*t21->Hi->IonStg == 0) || (t10->Hi->IonStg == t21->Hi->IonStg));

	ASSERT( (t20->Hi->IonStg*t21->Hi->IonStg == 0) || (t20->Hi->IonStg == t21->Hi->IonStg ) );

	ASSERT( (t10->Hi->nelem * t21->Hi->nelem == 0) || (t10->Hi->nelem == t21->Hi->nelem) );

	ASSERT( (t20->Hi->nelem * t21->Hi->nelem == 0) || (t20->Hi->nelem == t21->Hi->nelem) );

	ASSERT( o10 > 0. && o21 > 0. && o20 > 0. );

	/*end sanity checks */

	/* net loss of line escape and destruction */
	a21 = t21->Emis->Aul * (t21->Emis->Pesc+ t21->Emis->Pelec_esc + t21->Emis->Pdest);
	a10 = t10->Emis->Aul * (t10->Emis->Pesc+ t10->Emis->Pelec_esc + t10->Emis->Pdest);
	a20 = t20->Emis->Aul * (t20->Emis->Pesc+ t20->Emis->Pelec_esc + t20->Emis->Pdest);

	/* find energies of all transitions - one line could be a dummy
	 * also find Boltzmann factors */
	if( t10->Emis->Aul == 0. )
	{
		ener20 = t20->EnergyErg;
		ener21 = t21->EnergyErg;
		ener10 = ener20 - ener21;
		bolt12 = exp(-t21->EnergyK/phycon.te);
		bolt02 = exp(-t20->EnergyK/phycon.te);
		bolt01 = bolt02/bolt12;
		temp12 = t21->EnergyK;
		temp02 = t20->EnergyK;
		temp01 = temp02 - temp12;
	}

	else if( t21->Emis->Aul == 0. )
	{
		ener10 = t10->EnergyErg;
		ener20 = t20->EnergyErg;
		ener21 = ener20 - ener10;
		bolt01 = exp(-t10->EnergyK/phycon.te);
		bolt02 = exp(-t20->EnergyK/phycon.te);
		bolt12 = bolt02/bolt01;
		temp02 = t20->EnergyK;
		temp01 = t10->EnergyK;
		temp12 = temp02 - temp01;
	}

	else if( t20->Emis->Aul == 0. )
	{
		ener10 = t10->EnergyErg;
		ener21 = t21->EnergyErg;
		ener20 = ener21 + ener10;
		bolt01 = exp(-t10->EnergyK/phycon.te);
		bolt12 = exp(-t21->EnergyK/phycon.te);
		bolt02 = bolt01*bolt12;
		temp01 = t10->EnergyK;
		temp12 = t21->EnergyK;
		temp02 = temp01 + temp12;
	}

	else
	{
		/* all lines are ok */
		ener10 = t10->EnergyErg;
		ener20 = t20->EnergyErg;
		ener21 = t21->EnergyErg;
		bolt01 = exp(-t10->EnergyK/phycon.te);
		bolt12 = exp(-t21->EnergyK/phycon.te);
		bolt02 = bolt01*bolt12;
		temp02 = t20->EnergyK;
		temp01 = t10->EnergyK;
		temp12 = t21->EnergyK;
	}

	/* check all energies positive */
	ASSERT( ener10 > 0. && ener20 > 0. && ener21 > 0. );

	/* check if energy order is ok */
	ASSERT( ener10 < ener20 && ener21 < ener20 );

	/* check if energy scale is ok */
	ASSERT( fabs((ener10+ener21)/ener20-1.) < 1e-4 );

	pump01 = t10->Emis->pump;
	pump10 = pump01*g0/g1;
	pump12 = t21->Emis->pump;
	pump21 = pump12*g1/g2;
	pump02 = t20->Emis->pump;
	pump20 = pump02*g0/g2;

	/* cdsqte is 8.629E-6 / sqrte * eden */
	c01 = o10*bolt01*dense.cdsqte/g0;
	r01 = c01 + pump01;
	c10 = o10*dense.cdsqte/g1;
	r10 = c10 + a10 + pump10;
	c20 = o20*dense.cdsqte/g2;
	r20 = c20 + a20 + pump20;
	c02 = o20*bolt02*dense.cdsqte/g0;
	r02 = c02 + pump02;
	c12 = o21*bolt12*dense.cdsqte/g1;
	r12 = c12 + pump12;
	c21 = o21*dense.cdsqte/g2;
	r21 = c21 + a21 + pump21;

	alpha1 = (double)(AbunxIon)*(r01+r02)/(r10+r01+r02);
	alpha2 = (double)(AbunxIon)*(r01)/(r10+r12+r01);
	alpha = alpha1 - alpha2;

	/*  1( DBLE(r01+r02)/DBLE(r10+r01+r02) - DBLE(r01)/DBLE(r10+r12+r01) )
	 * beta is factor with n2 */
	beta = (r21 - r01)/(r10 + r12 + r01) + (r20 + r01 + r02)/(r10 + 
	  r01 + r02);

	if( alpha/MAX2(alpha1,alpha2) < 1e-11 )
	{
		/* this catches both negative and round off */
		p2 = 0.;
		alpha = 0.;
		nLev3Fail = 1;
	}

	else
	{
		p2 = alpha/beta;
	}
	atoms.PopLevels[2] = p2;

	if( alpha < 0. || beta < 0. )
	{
		fprintf( ioQQQ, " atom_level3: insane n2 pop alf, bet, p2=%10.2e%10.2e%10.2e %10.10s t=%10.2e\n", 
		  alpha, beta, p2, chLab, phycon.te );
		fprintf( ioQQQ, " gs are%5.1f%5.1f%5.1f\n", g0, g1, 
		  g2 );
		fprintf( ioQQQ, " Bolts are%10.2e%10.2e%10.2e\n", 
		  bolt01, bolt12, bolt02 );
		fprintf( ioQQQ, " As are%10.2e%10.2e%10.2e\n", a10, 
		  a21, a20 );
		fprintf( ioQQQ, " Energies are%10.2e%10.2e%10.2e\n", 
		  ener10, ener21, ener20 );
		fprintf( ioQQQ, " 2 terms, dif of alpha are%15.6e%15.6e\n", 
		  (r01 + r02)/(r10 + r01 + r02), r01/(r10 + r12 + r01) );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	alpha = (double)(AbunxIon)*(r01+r02) - (double)(p2)*(r20+r01+r02);
	/* guard against roundoff - this should really have been zero
	 * >>chng 96 nov 14, protection against round-off to zero
	 * >>chng 96 dec 03, made r01, etc, double, and changed limit to 1e-9 */
	if( fabs(alpha)/(MAX2(AbunxIon*(r01+r02),p2*(r20+r01+r02))) < 1e-9 )
	{
		alpha = 0.;
		nLev3Fail = 2;
	}

	beta = r10 + r01 + r02;
	p1 = alpha/beta;
	atoms.PopLevels[1] = p1;

	if( p1 < 0. )
	{
		if( p1 > -1e-37 )
		{
			/* slightly negative solution, probably just round-off, zero it */
			p1 = 0.;
			atoms.PopLevels[1] = p1;
			nLev3Fail = 3;
		}

		else
		{
			/* very negative solution, better punt */
			fprintf( ioQQQ, " atom_level3: insane n1 pop alf, bet, p1=%10.2e%10.2e%10.2e %10.10s%5f\n", 
			  alpha, beta, p1, chLab,  t10->WLAng );
			fprintf( ioQQQ, " local electron density and temperature were%10.2e%10.2e\n", 
			  dense.eden, phycon.te );
			ShowMe();
			cdEXIT(EXIT_FAILURE);
		}
	}

	p0 = AbunxIon - p2 - p1;

	/* population of lowest level */
	atoms.PopLevels[0] = p0;
	if( p0 <= 0. )
	{
		fprintf( ioQQQ, " atom_level3: insane n0 pop p1, 2, abun=%10.2e%10.2e%10.2e \n", 
		  p1, p2, AbunxIon );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	/* level populations for line opacities */
	t21->Lo->Pop = p1;
	t10->Lo->Pop = p0;
	t20->Lo->Pop = p0;
	t21->Emis->PopOpc = (p1 - p2*g1/g2);
	t10->Emis->PopOpc = (p0 - p1*g0/g1);
	t20->Emis->PopOpc = (p0 - p2*g0/g2);
	t21->Hi->Pop = p2;
	t10->Hi->Pop = p1;
	t20->Hi->Pop = p2;

	/* line emission - net emission in line */
	t21->Emis->phots = t21->Emis->Aul * (t21->Emis->Pesc + t21->Emis->Pelec_esc)*p2;
	t21->Emis->xIntensity = t21->Emis->phots * t21->EnergyErg;

	t20->Emis->phots = t20->Emis->Aul * (t20->Emis->Pesc + t20->Emis->Pelec_esc)*p2;
	t20->Emis->xIntensity = t20->Emis->phots * t20->EnergyErg;

	t10->Emis->phots = t10->Emis->Aul * (t10->Emis->Pesc + t10->Emis->Pelec_esc)*p1;
	t10->Emis->xIntensity = t10->Emis->phots * t10->EnergyErg;

#	if 0
	/* >>chng 03 oct 04, move to RT_OTS */
	t21->ots = p2 * t21->Emis->Aul * t21->Emis->Pdest;
	t20->ots = p2 * t20->Emis->Aul * t20->Emis->Pdest;
	t10->ots = p2 * t10->Emis->Aul * t10->Emis->Pdest;
	/*  now add thess lines to ots field, routine works on f not c scale */
	RT_OTS_AddLine( t21->ots , t21->ipCont);
	RT_OTS_AddLine( t20->ots , t20->ipCont);
	RT_OTS_AddLine( t10->ots , t10->ipCont);
#	endif

	/* total intensity used to divide line up - one may be 0 */
	/* >>chng 99 nov 30, rewrite algebra so double prec throughout,
	 * for very low density models */
	/*TotInten = t21->Emis->xIntensity + t20->xIntensity + t10->xIntensity;*/
	TotInten = t21->Emis->phots * (double)t21->EnergyErg 
		+ t20->Emis->phots * (double)t20->EnergyErg + 
		t10->Emis->phots * (double)t10->EnergyErg;

	/* fraction that was collisionally excited */
	if( r12 > 0. )
	{
		t21->Emis->ColOvTot = c12/r12;
	}
	else
	{
		t21->Emis->ColOvTot = 0.;
	}

	if( r01 > 0. )
	{
		t10->Emis->ColOvTot = c01/r01;
	}
	else
	{
		t10->Emis->ColOvTot = 0.;
	}

	if( r02 > 0. )
	{
		t20->Emis->ColOvTot = c02/r02;
	}
	else
	{
		t20->Emis->ColOvTot = 0.;
	}

	/* heating or cooling due to each line */
	heat20 = p2*c20*ener20;
	cool02 = p0*c02*ener20;
	heat21 = p2*c21*ener21;
	cool12 = p1*c12*ener21;
	heat10 = p1*c10*ener10;
	cool01 = p0*c01*ener10;

	/* two cases - collisionally excited (usual case) or 
	 * radiatively excited - in which case line can be a heat source
	 * following are correct heat exchange, they will mix to get correct deriv 
	 * the sum of heat-cool will add up to EnrUL - EnrLU - this is a trick to
	 * keep stable solution by effectively dividing up heating and cooling,
	 * so that negative cooling does not occur */

	/* now get net heating or cooling */
	cnet02 = cool02 - heat20*t20->Emis->ColOvTot;
	hnet02 = heat20 *(1. - t20->Emis->ColOvTot);
	cnet12 = cool12 - heat21*t21->Emis->ColOvTot;
	hnet12 = heat21 *(1. - t21->Emis->ColOvTot);
	cnet01 = cool01 - heat10*t10->Emis->ColOvTot;
	hnet01 = heat10 *(1. - t10->Emis->ColOvTot);

	/*TotalCooling = p0*(c01*ener10 + c02*ener20) + p1*c12*ener21 -
		(p2*(c21*ener21 + c20*ener20)  + p1*c10*ener10);*/
	/* >>chng 96 nov 22, very dense cool models, roundoff error
	 *could cause [OI] 63 mic to be dominant heating, cooling, or
	 *just zero
	 * >>chng 96 dec 06, above change caused o iii 1666 cooling to
	 *   be zeroed when important in a model in kirk's grid.  was at 1e-6,
	 *   set to 1e-7
	 * >>chng 96 dec 17, from 1e-7 to 1e-8, caused temp fail */
	/* >>chng 99 nov 29, had been 1e-30 to prevent div by zero (?),
	 * change to dble max since 1e-30 was very large compared to
	 * cooling for 1e-10 cm-3 den models */
	/*if( fabs(cnet01/MAX2(1e-30,cool01)) < 1e-8 )*/

	/* >>chng 02 jan 28, min from 1e-8 to 1e-10, conserve.in had massive
	 * temp failures when no molecules turned on, due to this tripping */
	/*if( fabs(cnet01/MAX2(DBL_MIN,cool01)) < 1e-8 )*/
	if( fabs(cnet01/MAX2(DBL_MIN,cool01)) < 1e-10 )
	{
		nLev3Fail = 4;
		cnet02 = 0.;
		hnet02 = 0.;
		cnet12 = 0.;
		hnet12 = 0.;
		cnet01 = 0.;
		hnet01 = 0.;
	}

	TotCool = cnet02 + cnet12 + cnet01;
	TotHeat = hnet02 + hnet12 + hnet01;


	if( TotInten > 0. )
	{
		cool02 = TotCool * t20->Emis->phots * (double)t20->EnergyErg/TotInten;
		cool12 = TotCool * t21->Emis->phots * (double)t21->EnergyErg/TotInten;
		cool01 = TotCool * t10->Emis->phots * (double)t10->EnergyErg/TotInten;
		heat20 = TotHeat * t20->Emis->phots * (double)t20->EnergyErg/TotInten;
		heat21 = TotHeat * t21->Emis->phots * (double)t21->EnergyErg/TotInten;
		heat10 = TotHeat * t10->Emis->phots * (double)t10->EnergyErg/TotInten;
		t20->Coll.cool = cool02;
		t21->Coll.cool = cool12;
		t10->Coll.cool = cool01;
		t20->Coll.heat = heat20;
		t21->Coll.heat = heat21;
		t10->Coll.heat = heat10;
	}
	else
	{
		nLev3Fail = 5;
		cool02 = 0.;
		cool12 = 0.;
		cool01 = 0.;
		heat20 = 0.;
		heat21 = 0.;
		heat10 = 0.;
		t20->Coll.cool = 0.;
		t21->Coll.cool = 0.;
		t10->Coll.cool = 0.;
		t20->Coll.heat = 0.;
		t21->Coll.heat = 0.;
		t10->Coll.heat = 0.;
	}

	/* add cooling due to each line,
	 * heating broken out above, will be added to thermal.heating[0][22] in CoolEvaluate*/
	/* >>chng 99 nov 30, rewrite algebra to keep precision on very low density models*/
	CoolAdd(chLab, t21->WLAng ,cool12);
	CoolAdd(chLab, t10->WLAng ,cool01);
	CoolAdd(chLab, t20->WLAng ,cool02);

	/* derivative of cooling function
	 * dC/dT = Cooling * ( T(excitation)/T_e^2 - 1/(2T) )
	 * in following I assume that a 1-2 exciation will have the 0-2
	 * exponential in the dcdt term = NOT A TYPO */

	thermal.dCooldT += t10->Coll.cool*(temp01*thermal.tsq1 - thermal.halfte) + 
	  (t20->Coll.cool + t21->Coll.cool)*(temp02*thermal.tsq1 - thermal.halfte);
	/* two t20->t's above are not a typo!
	 * */
	{
		enum{DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			fprintf(ioQQQ,"atom_level3 nLev3Fail %i\n", 
				nLev3Fail );
		}
	}
	return;
}
