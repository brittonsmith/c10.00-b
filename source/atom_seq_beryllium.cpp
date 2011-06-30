/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*AtomSeqBeryllium compute level populations and emissivity for Be-sequence ions */
#include "cddefines.h"
#include "phycon.h"
#include "transition.h"
#include "thirdparty.h"
#include "dense.h"
#include "cooling.h"
#include "thermal.h"
#include "atoms.h"

void AtomSeqBeryllium(double cs12, 
  double cs13, 
  double cs23, 
  transition * t, 
  double a30)
{

	char chLab[5];

	bool lgNegPop;

	int32 ipiv[4], nerror;
	long int i, j;

	double AbunxIon, 
	  Enr01, 
	  Enr10, 
	  a20, 
	  boltz, 
	  c01, 
	  c02, 
	  c03, 
	  c10, 
	  c12, 
	  c13, 
	  c20, 
	  c21, 
	  c23, 
	  c30, 
	  c31, 
	  c32, 
	  coolng, 
	  cs1u, 
	  excit, 
	  r01, 
	  r02, 
	  r03, 
	  r10, 
	  r12, 
	  r13, 
	  r20, 
	  r21, 
	  r23, 
	  r30, 
	  r31, 
	  r32, 
	  ri02, 
	  ri20, 
	  tot20;

	double amat[4][4], 
	  bvec[4], 
	  zz[5][5];

	DEBUG_ENTRY( "AtomSeqBeryllium()" );

	/* function returns total emission in both components of line
	 * destruction by background opacity computed and added to otslin field
	 * here,
	 * total cooling computed and added via CoolAdd
	 *
	 * finds population of four level Be-like atom
	 *
	 * level 1 = ground, J=0
	 * levels 2,3,4 are J=0,1,2, OF 3P.
	 * levels 3-1 is the fast one, j=1 to ground j=0
	 *
	 * cs1u is coll str to all J sub-levels of 3P; C.S. goes as 2J+1
	 * when routine exits the collision strength in the fast line is 1/3 of the entry value,
	 * so that save lines data does give the right critical density */

	/* AtomSeqBeryllium(cs12,cs13,cs23,t->t,a30,chLab)
	 * dense.xIonDense[Hi.nelem,i] is density of ith ionization stage (cm^-3) */
	AbunxIon = dense.xIonDense[t->Hi->nelem -1][t->Hi->IonStg-1];

	/* put null terminated line label into chLab using line array*/
	chIonLbl(chLab, t);
	boltz = t->EnergyK/phycon.te;

	/* set the cs before if below, since we must reset to line cs in all cases */
	cs1u = t->Coll.col_str;
	/* >>chng 01 sep 09, reset cs coming in, to propoer cs for the ^3P_1 level only
	 * so that critical density is printed correctly with save lines data command */
	t->Coll.col_str /= 3.f;

	/* low density cutoff to keep matrix happy */
	if( AbunxIon <= 1e-20 || boltz > 30. )
	{
		/* put in pop since possible just too cool */
		t->Lo->Pop = AbunxIon;
		t->Emis->PopOpc = AbunxIon;
		t->Hi->Pop = 0.;
		t->Emis->xIntensity = 0.;
		t->Coll.cool = 0.;
		t->Emis->phots = 0.;
		/*>>chng 03 oct 04, move to RT_OTS */
		/*t->ots = 0.;*/
		t->Emis->ColOvTot = 0.;
		t->Coll.heat = 0.;
		/* level populations */
		atoms.PopLevels[0] = AbunxIon;
		atoms.PopLevels[1] = 0.;
		atoms.PopLevels[2] = 0.;
		atoms.PopLevels[3] = 0.;
		atoms.DepLTELevels[0] = 1.;
		atoms.DepLTELevels[1] = 0.;
		atoms.DepLTELevels[2] = 0.;
		atoms.DepLTELevels[3] = 0.;
		CoolAdd(chLab, t->WLAng ,0.);
		return;
	}
	excit = exp(-boltz);

	/* these must be the statistical weights */
	ASSERT( t->Lo->g == 1. );
	ASSERT( t->Hi->g == 3. );

	/* collision strength must be positive */
	ASSERT( cs1u > 0. );

	/* incuded rates for fastest transition */
	ri02 = t->Emis->pump;

	/* back reaction has ratio of stat weights */
	ri20 = t->Emis->pump*1./3.;

	/* net rate out of level 3, with destruction */
	a20 = t->Emis->Aul*(t->Emis->Pesc + t->Emis->Pelec_esc + t->Emis->Pdest);
	tot20 = a20 + ri20;

	/* rates between j=0, lowest 3P level,
	 * 1/9 is ratio of level stat weight to term stat weight */
	c10 = cs1u*dense.cdsqte/9.;
	c01 = c10*excit;
	r01 = c01;
	r10 = c10;

	/* stat weights canceled out here */
	c20 = c10;
	c02 = c01*3.;
	r02 = c02 + ri02;
	r20 = c20 + tot20;

	c30 = c10;
	c03 = c01*5.;
	r30 = c30 + a30;
	r03 = c03;

	c21 = cs12*dense.cdsqte/3.;
	c12 = c21*3.;
	r12 = c12;
	r21 = c21;

	c31 = cs13*dense.cdsqte/5.;
	c13 = c31*5.;
	r31 = c31;
	r13 = c13;

	c32 = cs23*dense.cdsqte/5.;
	c23 = c32*1.667;
	r32 = c32;
	r23 = c23;

	/* set up matrix */
	for( i=0; i <= 3; i++ )
	{
		/* first equation will be sum to abund */
		zz[i][0] = 1.;
		zz[4][i] = 0.;
	}

	/* zz(0,4) = AbunxIon */
	zz[4][0] = 1.;

	/* ground level 0 is implicit
	 * level 1 balance equation */
	zz[0][1] = -r01;
	zz[1][1] = r10 + r12 + r13;
	zz[2][1] = -r21;
	zz[3][1] = -r31;

	/* level 2 balance equation */
	zz[0][2] = -r02;
	zz[1][2] = -r12;
	zz[2][2] = r20 + r21 + r23;
	zz[3][2] = -r32;

	/* level 3 balance equation */
	zz[0][3] = -r03;
	zz[1][3] = -r13;
	zz[2][3] = -r23;
	zz[3][3] = r30 + r31 + r32;

	/* this one may be more robust */
	for( j=0; j <= 3; j++ )
	{
		for( i=0; i <= 3; i++ )
		{
			amat[i][j] = zz[i][j];
		}
		bvec[j] = zz[3+1][j];
	}

	nerror = 0;

  	getrf_wrapper(4, 4, (double*)amat, 4, ipiv, &nerror);
	getrs_wrapper('N', 4, 1, (double*)amat, 4, ipiv, bvec, 4, &nerror);

	if( nerror != 0 )
	{
		fprintf( ioQQQ, " AtomSeqBeryllium: dgetrs finds singular or ill-conditioned matrix\n" );
		cdEXIT(EXIT_FAILURE);
	}
	/* now put results back into z so rest of code treates only
		* one case - as if matin1 had been used */
	for( i=0; i <= 3; i++ )
	{
		zz[3+1][i] = bvec[i];
	}

	lgNegPop = false;
	for( i=0; i <= 3; i++ )
	{
		atoms.PopLevels[i] = zz[4][i]*AbunxIon;
		if( atoms.PopLevels[i] < 0. )
			lgNegPop = true;
	}

	/* check for negative level populations, this would be a major error */
	if( lgNegPop )
	{
		fprintf( ioQQQ, " AtomSeqBeryllium finds non-positive pop,=" );
		for( i=0; i <= 3; i++ )
		{
			fprintf( ioQQQ, "%g ", atoms.PopLevels[i] );
		}
		fprintf( ioQQQ, "%s \n", chLab );
		fprintf( ioQQQ, " te=%g  total abund=%g  boltz=%g \n", 
		  phycon.te, AbunxIon, boltz );
		cdEXIT(EXIT_FAILURE);
	}

	/* convert level populations over to departure coeficients relative to ground */
	atoms.DepLTELevels[0] = 1.;
	atoms.DepLTELevels[1] = (atoms.PopLevels[1]/atoms.PopLevels[0])/
	  excit;
	atoms.DepLTELevels[2] = (atoms.PopLevels[2]/atoms.PopLevels[0])/
	  (excit*3.);
	atoms.DepLTELevels[3] = (atoms.PopLevels[3]/atoms.PopLevels[0])/
	  (excit*5.);

	/* compute ratio Aul/(Aul+Cul) */
	t->Emis->ColOvTot = c02/r02;

	t->Lo->Pop = AbunxIon;

	/* >>chng 96 sep 12, ipLnPopu had not been set before */
	t->Hi->Pop = atoms.PopLevels[2];

	t->Emis->PopOpc = AbunxIon - atoms.PopLevels[2]*1./3.;

	/* this will be escaping part of line
	 * number of escaping line photons, used elsewhere for outward beam */
	t->Emis->phots = atoms.PopLevels[2] * t->Emis->Aul * ( t->Emis->Pesc + t->Emis->Pelec_esc );

	t->Emis->xIntensity = t->Emis->phots * t->EnergyErg;

	/* two cases - collisionally excited (usual case) or 
	 * radiatively excited - in which case line can be a heat source
	 * following are correct heat exchange, they will mix to get correct deriv 
	 * the sum of heat-cool will add up to EnrUL - EnrLU - this is a trick to
	 * keep stable solution by effectively dividing up heating and cooling,
	 * so that negative cooling does not occur */
	Enr01 = atoms.PopLevels[0]*(c01 + c02 + c03)*t->EnergyErg;
	Enr10 = (atoms.PopLevels[1]*c10 + atoms.PopLevels[2]*c20 + 
	  atoms.PopLevels[3]*c30)*t->EnergyErg;

	/* net cooling due to excit minus part of de-excit */
	t->Coll.cool = Enr01 - Enr10*t->Emis->ColOvTot;

	/* net heating is remainder */
	t->Coll.heat = Enr10*(1. - t->Emis->ColOvTot);

	/* put line cooling into cooling stack */
	coolng = t->Coll.cool;
	CoolAdd(chLab, t->WLAng ,coolng);

	/* derivative of cooling function */
	thermal.dCooldT += coolng*(t->EnergyK*thermal.tsq1 - thermal.halfte);

	/* >>chhg 03 oct 04, move to RT_OTS */
	/* destroyed part of line
	dest = atoms.PopLevels[2]*t->Emis->Aul*t->Pdest;
	t->ots = dest; */

	/* now add to ots field 
	 * >>chng 03 oct 03, moved to RT_OTS
	RT_OTS_AddLine(dest , t->ipCont );*/
	return;
}
