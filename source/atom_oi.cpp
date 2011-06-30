/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*atom_oi drive the solution of OI level populations, Ly-beta pumping */
/*oi_level_pops get OI level population with Ly-beta pumping */
#include "cddefines.h"
#include "taulines.h"
#include "doppvel.h"
#include "iso.h"
#include "trace.h"
#include "dense.h"
#include "rt.h"
#include "rfield.h"
#include "phycon.h"
#include "lines_service.h"
#include "thirdparty.h"
#include "atoms.h"

/*oi_level_pops get OI level population with Ly-beta pumping */
STATIC void oi_level_pops(double abundoi, 
			  double *coloi);

/*atom_oi drive the solution of OI level populations, Ly-beta pumping */
void atom_oi_calc(double *coloi)
{
	double esab,
	  eslb, 
	  esoi, 
	  flb, 
	  foi, 
	  opaclb, 
	  opacoi, 
	  xlb, 
	  xoi;
	double aoi = TauLines[ipTO1025].Emis->Aul;
	double alb = Transitions[ipH_LIKE][ipHYDROGEN][ipH3p][ipH1s].Emis->Aul;

	DEBUG_ENTRY( "atom_oi_calc()" );

	fixit(); // ticket #78 refers
	// The code below should be calculating the O I 1025 pumping by H Ly beta, as well as
	// the inverse process (this can become important in hydrogen-deficient environments).
	// It now uses the Elitzur & Netzer (1985, ApJ, 291, 464) theory, which is no longer
	// valid since the line overlap code prevents us from getting at the escape probability
	// of individual lines.

	/* A's from Pradhan; OI pump line; Ly beta, 8446 */

	/* called by LINES before calc really start, protect here
	 * also used for cases where OI not present */
	if( dense.xIonDense[ipOXYGEN][0] <= 0. )
	{
		for( int i=0; i < 6; i++ )
		{
			atoms.popoi[i] = 0.;
		}
		*coloi = 0.;
		atoms.pmpo15 = 0.;
		atoms.pmpo51 = 0.;
		return;
	}

	// line overlap code makes this the escape probability of the combined lines
	esab = TauLines[ipTO1025].Emis->Pelec_esc + TauLines[ipTO1025].Emis->Pesc;

	// these two are no longer correct, the line overlap code makes it impossible
	// to get at the escape probabilities of the individual lines...
	esoi = TauLines[ipTO1025].Emis->Pelec_esc + TauLines[ipTO1025].Emis->Pesc;
	eslb = Transitions[ipH_LIKE][ipHYDROGEN][ipH3p][ipH1s].Emis->Pelec_esc + 
		Transitions[ipH_LIKE][ipHYDROGEN][ipH3p][ipH1s].Emis->Pesc;

	/* all trace output turned on with "trace ly beta command' */
	if( trace.lgTr8446 && trace.lgTrace )
	{
		fprintf( ioQQQ, 
			"       P8446 finds Lbeta, OI widths=%10.2e%10.2e and esc prob=%10.2e%10.2e esAB=%10.2e\n", 
		  GetDopplerWidth(dense.AtomicWeight[ipHYDROGEN]), GetDopplerWidth(dense.AtomicWeight[ipOXYGEN]), eslb, esoi, esab );
	}

	/* find relative opacities for two lines */
	opacoi = 2.92e-9*dense.xIonDense[ipOXYGEN][0]*0.5556/GetDopplerWidth(dense.AtomicWeight[ipOXYGEN]);
	opaclb = 1.22e-8*StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop/GetDopplerWidth(dense.AtomicWeight[ipHYDROGEN]);

	/* these are x sub a (OI) and x sub b (ly beta) defined in Elitz+Netz */
	xoi = opacoi/(opacoi + opaclb);
	xlb = opaclb/(opacoi + opaclb);

	/* find relative line-widths, assume same rest freq */
	foi = MIN2(GetDopplerWidth(dense.AtomicWeight[ipHYDROGEN]),GetDopplerWidth(dense.AtomicWeight[ipOXYGEN]))/GetDopplerWidth(dense.AtomicWeight[ipOXYGEN]);
	flb = MIN2(GetDopplerWidth(dense.AtomicWeight[ipHYDROGEN]),GetDopplerWidth(dense.AtomicWeight[ipOXYGEN]))/GetDopplerWidth(dense.AtomicWeight[ipHYDROGEN])*
		MAX2(0.,1.- TauLines[ipTO1025].Emis->Pelec_esc - TauLines[ipTO1025].Emis->Pesc);

	if( trace.lgTr8446 && trace.lgTrace )
	{
		fprintf( ioQQQ, 
			"       P8446 opac Lb, OI=%10.2e%10.2e X Lb, OI=%10.2e%10.2e FLb, OI=%10.2e%10.2e\n", 
		  opaclb, opacoi, xlb, xoi, flb, foi );
	}

	/* pumping of OI by L-beta - this goes into OI matrix as 1-5 rate
	 * lgInducProcess set false with no induced command, usually true */
	/* Notation: pmpo15 net rate oxygen level 5 populated from 1 */
	if( rfield.lgInducProcess )
	{
		atoms.pmpo15 = (realnum)((flb*alb*StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH3p].Pop*
					 xoi*(1. - esab)/dense.xIonDense[ipOXYGEN][0]));
		/* net decay rate from upper level */
		atoms.pmpo51 = (realnum)(aoi*(1. - (1. - foi)*(1. - esoi) - xoi*(1. - esab)*foi));
	}
	else
	{
		atoms.pmpo15 = 0.;
		atoms.pmpo51 = 0.;
	}

	/* find level populations for OI */
	oi_level_pops(dense.xIonDense[ipOXYGEN][0],coloi);

	/* continuum pumping due to J=1, 0 sub states.
	 * neglect J=2 since L-Beta very optically thick */

	/* lower level populations */
	TauLines[ipT1304].Lo->Pop = atoms.popoi[0];
	TauLines[ipTO1025].Lo->Pop = atoms.popoi[0];
	TauLines[ipT1039].Lo->Pop = atoms.popoi[0];
	TauLines[ipT8446].Lo->Pop = atoms.popoi[1];
	TauLines[ipT4368].Lo->Pop = atoms.popoi[1];
	TauLines[ipTOI13].Lo->Pop = atoms.popoi[2];
	TauLines[ipTOI11].Lo->Pop = atoms.popoi[2];
	TauLines[ipTOI29].Lo->Pop = atoms.popoi[3];
	TauLines[ipTOI46].Lo->Pop = atoms.popoi[4];

	/* upper level populations */
	TauLines[ipT1304].Hi->Pop = atoms.popoi[1];
	TauLines[ipTO1025].Hi->Pop = atoms.popoi[4];
	TauLines[ipT1039].Hi->Pop = atoms.popoi[3];
	TauLines[ipT8446].Hi->Pop = atoms.popoi[2];
	TauLines[ipT4368].Hi->Pop = atoms.popoi[5];
	TauLines[ipTOI13].Hi->Pop = atoms.popoi[3];
	TauLines[ipTOI11].Hi->Pop = atoms.popoi[4];
	TauLines[ipTOI29].Hi->Pop = atoms.popoi[5];
	TauLines[ipTOI46].Hi->Pop = atoms.popoi[5];

	/* opacity factor
	 * t1304(ipLnEmis.PopOpc) = (popoi(1)-popoi(2)*3.0) */
	/** \todo	2	following needed to get badbugs/bug8.in to work */
	TauLines[ipT1304].Emis->PopOpc = atoms.popoi[0];
	TauLines[ipTO1025].Emis->PopOpc = (atoms.popoi[0] - atoms.popoi[4]*0.6);
	TauLines[ipT1039].Emis->PopOpc = (atoms.popoi[0] - atoms.popoi[3]*3.0);

	/* t8446(ipLnEmis.PopOpc) = (popoi(2)-popoi(3)*.33) */
	/** \todo	2	following needed to get badbugs/bug5.in to work */
	TauLines[ipT8446].Emis->PopOpc = 
		(MAX2(0.,atoms.popoi[1]-atoms.popoi[2]*.33));
	TauLines[ipT4368].Emis->PopOpc = (atoms.popoi[1] - atoms.popoi[5]*.33);
	TauLines[ipTOI13].Emis->PopOpc = (atoms.popoi[2] - atoms.popoi[3]*3.0);
	TauLines[ipTOI11].Emis->PopOpc = (atoms.popoi[2] - atoms.popoi[4]*0.6);
	TauLines[ipTOI29].Emis->PopOpc = (atoms.popoi[3] - atoms.popoi[5]*.33);
	TauLines[ipTOI46].Emis->PopOpc = (atoms.popoi[4] - atoms.popoi[5]*1.67);
	return;
}

/*oi_level_pops get OI level population with Ly-beta pumping */
STATIC void oi_level_pops(double abundoi, 
			  double *coloi)
{
	bool lgNegPop;

	long int i, j;

	int32 ipiv[6], ner;

	double a21, 
	  a32, 
	  a41, 
	  a43, 
	  a51, 
	  a53, 
	  a62, 
	  a64, 
	  a65, 
	  c12, 
	  c13, 
	  c14, 
	  c15, 
	  c16, 
	  c21, 
	  c23, 
	  c24, 
	  c25, 
	  c26, 
	  c31, 
	  c32, 
	  c34, 
	  c35, 
	  c36, 
	  c41, 
	  c42, 
	  c43, 
	  c45, 
	  c46, 
	  c51, 
	  c52, 
	  c53, 
	  c54, 
	  c56, 
	  c61, 
	  c62, 
	  c63, 
	  c64, 
	  c65, 
	  cs, 
	  deptoi[6], 
	  e12, 
	  e23, 
	  e34, 
	  e45, 
	  e56, 
	  simple;

	double amat[6][6], 
	  bvec[6], 
	  zz[7][7];

	static double g[6]={9.,3.,9.,3.,15.,9};

	DEBUG_ENTRY( "oilevl()" );

	/* following used for linpac matrix inversion */

	/* compute emission from six level OI atom*/

	/* Boltzmann factors for ContBoltz since collisions not dominant in UV tran
	 * ipoiex is array lof dEnergy for each level, set in DoPoint */
	e12 = rfield.ContBoltz[atoms.ipoiex[0]-1];
	e23 = rfield.ContBoltz[atoms.ipoiex[1]-1];
	e34 = rfield.ContBoltz[atoms.ipoiex[2]-1];
	e45 = rfield.ContBoltz[atoms.ipoiex[3]-1];
	e56 = rfield.ContBoltz[atoms.ipoiex[4]-1];

	/* total rad rates here have dest by background continuum */
	a21 = TauLines[ipT1304].Emis->Aul*(TauLines[ipT1304].Emis->Pdest+ TauLines[ipT1304].Emis->Pesc + TauLines[ipT1304].Emis->Pelec_esc);
	a41 = TauLines[ipT1039].Emis->Aul*(TauLines[ipT1039].Emis->Pdest+ TauLines[ipT1039].Emis->Pesc + TauLines[ipT1039].Emis->Pelec_esc);
	a51 = TauLines[ipTO1025].Emis->Aul*(TauLines[ipTO1025].Emis->Pdest+ TauLines[ipTO1025].Emis->Pesc + TauLines[ipTO1025].Emis->Pelec_esc);
	a51 = atoms.pmpo51;
	a32 = TauLines[ipT8446].Emis->Aul*(TauLines[ipT8446].Emis->Pdest+ TauLines[ipT8446].Emis->Pesc + TauLines[ipT8446].Emis->Pelec_esc);
	a62 = TauLines[ipT4368].Emis->Aul*(TauLines[ipT4368].Emis->Pdest+ TauLines[ipT4368].Emis->Pesc + TauLines[ipT4368].Emis->Pelec_esc);
	a43 = TauLines[ipTOI13].Emis->Aul*(TauLines[ipTOI13].Emis->Pdest+ TauLines[ipTOI13].Emis->Pesc + TauLines[ipTOI13].Emis->Pelec_esc);
	a53 = TauLines[ipTOI11].Emis->Aul*(TauLines[ipTOI11].Emis->Pdest+ TauLines[ipTOI11].Emis->Pesc + TauLines[ipTOI11].Emis->Pelec_esc);
	a64 = TauLines[ipTOI29].Emis->Aul*(TauLines[ipTOI29].Emis->Pdest+ TauLines[ipTOI29].Emis->Pesc + TauLines[ipTOI29].Emis->Pelec_esc);
	a65 = TauLines[ipTOI46].Emis->Aul*(TauLines[ipTOI46].Emis->Pdest+ TauLines[ipTOI46].Emis->Pesc + TauLines[ipTOI46].Emis->Pelec_esc);

	/* even at density of 10^17 excited states not in lte due
	 * to fast transitions down - just need to kill a21 to get to unity at 10^17*/

	/* the 2-1 transition is 1302, cs Wang and McConkey '92 Jphys B 25, 5461 */
	cs = 2.151e-5*phycon.te/phycon.te03;
	PutCS(cs,&TauLines[ipT1304]);

	/* the 5-1 transition is 1027, cs Wang and McConkey '92 Jphys B 25, 5461 */
	cs = 9.25e-7*phycon.te*phycon.te10/phycon.te01/phycon.te01;
	PutCS(cs,&TauLines[ipTO1025]);
	c21 = dense.cdsqte*TauLines[ipT1304].Coll.col_str/g[1];
	c51 = dense.cdsqte*TauLines[ipTO1025].Coll.col_str/g[4];

	/* all following are g-bar approx, g-bar = 0.2 */
	c31 = dense.cdsqte*1.0/g[2];
	PutCS(0.27,&TauLines[ipT1039]);
	c41 = dense.cdsqte*TauLines[ipT1039].Coll.col_str/g[3];
	c61 = dense.cdsqte*1./g[5];

	c12 = c21*g[1]/g[0]*e12;
	c13 = c31*g[2]/g[0]*e12*e23;
	c14 = c41*g[3]/g[0]*e12*e23*e34;
	c15 = c51*g[4]/g[0]*e12*e23*e34*e45;
	c16 = c61*g[5]/g[0]*e12*e23*e34*e45*e56;

	c32 = dense.cdsqte*85./g[2];
	c42 = dense.cdsqte*85./g[3];
	c52 = dense.cdsqte*85./g[4];
	c62 = dense.cdsqte*85./g[5];

	c23 = c32*g[2]/g[1]*e23;
	c24 = c42*g[3]/g[1]*e23*e34;
	c25 = c52*g[4]/g[1]*e23*e34*e45;
	c26 = c62*g[5]/g[1]*e23*e34*e45*e56;

	c43 = dense.cdsqte*70./g[3];
	c53 = dense.cdsqte*312./g[4];
	c63 = dense.cdsqte*1./g[5];

	c34 = c43*g[3]/g[2]*e34;
	c35 = c53*g[4]/g[2]*e34*e45;
	c36 = c63*g[5]/g[2]*e34*e45*e56;

	c54 = dense.cdsqte*50./g[4];
	c64 = dense.cdsqte*415./g[5];

	c45 = c54*g[4]/g[3]*e45;
	c46 = c64*g[5]/g[3]*e45*e56;

	c65 = dense.cdsqte*400./g[5];
	c56 = c65*g[5]/g[4]*e56;

	/** \todo	2	this must have all stimulated emission, pump by cont, etc*/

	/* this is check for whether matrix inversion likely to fail */
	simple = (c16 + atoms.pmpo15)/(c61 + c62 + c64 + a65 + a64 + a62);
	if( simple < 1e-19 )
	{
		atoms.popoi[0] = abundoi;
		for( i=1; i < 6; i++ )
		{
			atoms.popoi[i] = 0.;
		}
		*coloi = 0.;
		return;
	}

	/*--------------------------------------------------------- */

	for( i=0; i < 6; i++ )
	{
		zz[i][0] = 1.0;
		zz[6][i] = 0.;
	}

	/* first equation is sum of populations */
	zz[6][0] = abundoi;

	/* level two, 3s 3So */
	zz[0][1] = -c12;
	zz[1][1] = c21 + c23 + c24 + c25 + c26 + a21;
	zz[2][1] = -c32 - a32;
	zz[3][1] = -c42;
	zz[4][1] = -c52;
	zz[5][1] = -c62 - a62;

	/* level three */
	zz[0][2] = -c13;
	zz[1][2] = -c23;
	zz[2][2] = c31 + c32 + c34 + c35 + c36 + a32;
	zz[3][2] = -c43 - a43;
	zz[4][2] = -c53 - a53;
	zz[5][2] = -c63;

	/* level four */
	zz[0][3] = -c14;
	zz[1][3] = -c24;
	zz[2][3] = -c34;
	zz[3][3] = c41 + c42 + c43 + c45 + c46 + a41 + a43;
	zz[4][3] = -c54;
	zz[5][3] = -c64 - a64;

	/* level five */
	zz[0][4] = -c15 - atoms.pmpo15;
	zz[1][4] = -c25;
	zz[2][4] = -c35;
	zz[3][4] = -c45;
	zz[4][4] = c51 + c52 + c53 + c54 + c56 + a51 + a53;
	zz[5][4] = -c65 - a65;

	/* level six */
	zz[0][5] = -c16;
	zz[1][5] = -c26;
	zz[2][5] = -c36;
	zz[3][5] = -c46;
	zz[4][5] = -c56;
	zz[5][5] = c61 + c62 + c63 + c64 + c65 + a65 + a64 + a62;

	/* this one may be more robust */
	for( j=0; j < 6; j++ )
	{
		for( i=0; i < 6; i++ )
		{
			amat[i][j] = zz[i][j];
		}
		bvec[j] = zz[6][j];
	}

	ner = 0;

  	getrf_wrapper(6, 6, (double*)amat, 6, ipiv, &ner);
	getrs_wrapper('N', 6, 1, (double*)amat, 6, ipiv, bvec, 6, &ner);

	/*DGETRF(6,6,(double*)amat,6,ipiv,&ner);*/
	/*DGETRS('N',6,1,(double*)amat,6,ipiv,bvec,6,&ner);*/
	if( ner != 0 )
	{
		fprintf( ioQQQ, " oi_level_pops: dgetrs finds singular or ill-conditioned matrix\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* now put results back into z so rest of code treates only
		* one case - as if matin1 had been used */
	for( i=0; i < 6; i++ )
	{
		zz[6][i] = bvec[i];
	}

	lgNegPop = false;
	for( i=0; i < 6; i++ )
	{
		atoms.popoi[i] = zz[6][i];
		if( atoms.popoi[i] < 0. )
			lgNegPop = true;
	}

	/* following used to confirm that all dep coef are unity at
	 * density of 1e17, t=10,000, and all A's set to zero */
	if( trace.lgTrace && trace.lgTr8446 )
	{
		deptoi[0] = 1.;
		deptoi[1] = atoms.popoi[1]/atoms.popoi[0]/(g[1]/g[0]*
		  e12);
		deptoi[2] = atoms.popoi[2]/atoms.popoi[0]/(g[2]/g[0]*
		  e12*e23);
		deptoi[3] = atoms.popoi[3]/atoms.popoi[0]/(g[3]/g[0]*
		  e12*e23*e34);
		deptoi[4] = atoms.popoi[4]/atoms.popoi[0]/(g[4]/g[0]*
		  e12*e23*e34*e45);
		deptoi[5] = atoms.popoi[5]/atoms.popoi[0]/(g[5]/g[0]*
		  e12*e23*e34*e45*e56);

		fprintf( ioQQQ, " oilevl finds levl pop" );
		for(i=0; i < 6; i++)
			fprintf( ioQQQ, "%11.3e", atoms.popoi[i] );
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, " oilevl finds dep coef" );
		for(i=0; i < 6; i++)
			fprintf( ioQQQ, "%11.3e", deptoi[i] );
		fprintf( ioQQQ, "\n" );
	}

	/* this happens due to numerical instability in matrix inversion routine */
	if( lgNegPop )
	{
		atoms.nNegOI += 1;

		fprintf( ioQQQ, " OILEVL finds negative population" );
		for(i=0;i < 6; i++)
			fprintf( ioQQQ, "%10.2e", atoms.popoi[i] );
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, " simple 5 =%10.2e\n", simple );

		atoms.popoi[5] = 0.;
		atoms.popoi[4] = abundoi*(c15 + atoms.pmpo15)/(a51 + a53 + c51 + c53);
		atoms.popoi[3] = 0.;
		atoms.popoi[2] = (atoms.popoi[4]*(a53 + c53)) / (a32 + c32);
		atoms.popoi[1] = (atoms.popoi[4]*(a53 + c53) + abundoi*c12) / (a21 + c21);
		atoms.popoi[0] = abundoi;
		/*  write(QQ,'('' OILEVL resets this to simple pop'',1P,6E10.2)')
		 *  1   popoi */
	}

	/* this is total cooling due to model atom, can be neg (heating) */
	*coloi = 
	   (atoms.popoi[0]*c12 - atoms.popoi[1]*c21)*1.53e-11 + 
	   (atoms.popoi[0]*c14 - atoms.popoi[3]*c41)*1.92e-11 + 
	   (atoms.popoi[0]*c15 - atoms.popoi[4]*c51)*1.94e-11 + 
	   (atoms.popoi[1]*c23 - atoms.popoi[2]*c32)*2.36e-12 + 
	   (atoms.popoi[1]*c26 - atoms.popoi[5]*c62)*4.55e-12 + 
	   (atoms.popoi[2]*c35 - atoms.popoi[4]*c53)*1.76e-12 + 
	   (atoms.popoi[2]*c34 - atoms.popoi[3]*c43)*1.52e-12 + 
	   (atoms.popoi[3]*c46 - atoms.popoi[5]*c64)*6.86e-13 + 
	   (atoms.popoi[4]*c56 - atoms.popoi[5]*c65)*4.32e-13;

	return;
}
