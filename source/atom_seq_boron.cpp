/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*AtomSeqBoron compute cooling from 5-level boron sequence model atom */
#include "cddefines.h"
#include "cooling.h"
#include "thermal.h"
#include "dense.h"
#include "atoms.h"
#include "transition.h"

void AtomSeqBoron(
	/* indices for all lines are on the C scale since they will be stuffed into
	 * C arrays.  so, t10 refers to the 2-1 transition */
	transition * t10, 
	transition * t20, 
	transition * t30, 
	transition * t21, 
	transition * t31, 
	transition * t41, 
	double cs40,
	double cs32,
	double cs42,
	double cs43,
	/* pump rate s-1 due to UV permitted lines */
	double pump_rate ,
	/* string used to identify calling program in case of error */
	const char *chLabel
	)
{

	/* this routine has three possible returns:
	 * abundance is zero
	 * too cool for full 5-level atom, but still do ground term 
	 * full solution */

	/* boron sequence is now a five level atom */
#	define	N_SEQ_BORON	5
	static double 
		**AulEscp ,
		**col_str ,
		**AulDest, 
		/* AulPump[low][high] is rate (s^-1) from lower to upper level */
		**AulPump,
		**CollRate,
		*pops,
		*create,
		*destroy,
		*depart,
		/* statistical weight */
		*stat ,
		/* excitation energies in kelvin */
		*excit;

	double b_cooling,
	  dCoolDT;
	double EnrLU, EnrUL;
	realnum abundan;

	static bool lgFirst=true,
		lgZeroPop;
	int i , j;
	int 
		/* flag to signal negative level populations */
		lgNegPop;
		/* flag to turn on debug print in atom_levelN */
	bool lgDeBug;

	DEBUG_ENTRY( "AtomSeqBoron()" );

	if( lgFirst )
	{
		/* will never do this again */
		lgFirst = false;
		/* allocate the 1D arrays*/
		excit = (double *)MALLOC( sizeof(double)*(N_SEQ_BORON) );
		stat = (double *)MALLOC( sizeof(double)*(N_SEQ_BORON) );
		pops = (double *)MALLOC( sizeof(double)*(N_SEQ_BORON) );
		create = (double *)MALLOC( sizeof(double)*(N_SEQ_BORON) );
		destroy = (double *)MALLOC( sizeof(double)*(N_SEQ_BORON) );
		depart = (double *)MALLOC( sizeof(double)*(N_SEQ_BORON) );
		/* create space for the 2D arrays */
		AulPump = ((double **)MALLOC((N_SEQ_BORON)*sizeof(double *)));
		CollRate = ((double **)MALLOC((N_SEQ_BORON)*sizeof(double *)));
		AulDest = ((double **)MALLOC((N_SEQ_BORON)*sizeof(double *)));
		AulEscp = ((double **)MALLOC((N_SEQ_BORON)*sizeof(double *)));
		col_str = ((double **)MALLOC((N_SEQ_BORON)*sizeof(double *)));
		for( i=0; i<(N_SEQ_BORON); ++i )
		{
			AulPump[i] = ((double *)MALLOC((N_SEQ_BORON)*sizeof(double )));
			CollRate[i] = ((double *)MALLOC((N_SEQ_BORON)*sizeof(double )));
			AulDest[i] = ((double *)MALLOC((N_SEQ_BORON)*sizeof(double )));
			AulEscp[i] = ((double *)MALLOC((N_SEQ_BORON)*sizeof(double )));
			col_str[i] = ((double *)MALLOC((N_SEQ_BORON)*sizeof(double )));
		}
	}

	/* total abundance of this species */
	abundan = dense.xIonDense[ t10->Hi->nelem -1][t10->Hi->IonStg-1];

	/** \todo	2	use TransitionZero here */
	if( abundan <= 0. )
	{
		/* this branch, no abundance of ion */
		t10->Lo->Pop = 0.;
		t20->Lo->Pop = 0.;
		t30->Lo->Pop = 0.;
		t21->Lo->Pop = 0.;
		t31->Lo->Pop = 0.;
		t41->Lo->Pop = 0.;

		t10->Emis->PopOpc = 0.;
		t20->Emis->PopOpc = 0.;
		t30->Emis->PopOpc = 0.;
		t21->Emis->PopOpc = 0.;
		t31->Emis->PopOpc = 0.;
		t41->Emis->PopOpc = 0.;

		t10->Hi->Pop = 0.;
		t20->Hi->Pop = 0.;
		t30->Hi->Pop = 0.;
		t21->Hi->Pop = 0.;
		t31->Hi->Pop = 0.;
		t41->Hi->Pop = 0.;

		t10->Emis->xIntensity = 0.;
		t20->Emis->xIntensity = 0.;
		t30->Emis->xIntensity = 0.;
		t21->Emis->xIntensity = 0.;
		t31->Emis->xIntensity = 0.;
		t41->Emis->xIntensity = 0.;

		t10->Coll.cool = 0.;
		t20->Coll.cool = 0.;
		t30->Coll.cool = 0.;
		t21->Coll.cool = 0.;
		t31->Coll.cool = 0.;
		t41->Coll.cool = 0.;

		t10->Emis->phots = 0.;
		t20->Emis->phots = 0.;
		t30->Emis->phots = 0.;
		t21->Emis->phots = 0.;
		t31->Emis->phots = 0.;
		t41->Emis->phots = 0.;

		t10->Emis->ColOvTot = 0.;
		t20->Emis->ColOvTot = 0.;
		t30->Emis->ColOvTot = 0.;
		t21->Emis->ColOvTot = 0.;
		t31->Emis->ColOvTot = 0.;
		t41->Emis->ColOvTot = 0.;

		t10->Coll.heat = 0.;
		t20->Coll.heat = 0.;
		t30->Coll.heat = 0.;
		t21->Coll.heat = 0.;
		t31->Coll.heat = 0.;
		t41->Coll.heat = 0.;

		CoolAdd( chLabel, t10->WLAng , 0.);
		CoolAdd( chLabel, t20->WLAng , 0.);
		CoolAdd( chLabel, t30->WLAng , 0.);
		CoolAdd( chLabel, t21->WLAng , 0.);
		CoolAdd( chLabel, t31->WLAng , 0.);
		CoolAdd( chLabel, t41->WLAng , 0.);

		/* level populations */
		/* LIMLEVELN is the dimension of the atoms vectors */
		ASSERT( N_SEQ_BORON <= LIMLEVELN);
		for( i=0; i < N_SEQ_BORON; i++ )
		{
			atoms.PopLevels[i] = 0.;
			atoms.DepLTELevels[i] = 1.;
		}
		return;
	}

	ASSERT( t10->Coll.col_str > 0.);
	ASSERT( t20->Coll.col_str > 0.);
	ASSERT( t30->Coll.col_str > 0.);
	ASSERT( t21->Coll.col_str > 0.);
	ASSERT( t31->Coll.col_str > 0.);
	ASSERT( t41->Coll.col_str > 0.);
	ASSERT( cs40>0.);
	ASSERT( cs32>0.);
	ASSERT( cs42>0.);
	ASSERT( cs43>0.);

	/* all elements are used, and must be set to zero if zero */
	for( i=0; i < N_SEQ_BORON; i++ )
	{
		create[i] = 0.;
		destroy[i] = 0.;
		for( j=0; j < N_SEQ_BORON; j++ )
		{
			/*data[j][i] = -1e33;*/
			AulEscp[j][i] = 0.;
			AulDest[j][i] = 0.;
			AulPump[j][i] = 0.;
			col_str[j][i] = 0.;
		}
	}

	/* statistical weights */
	stat[0] = t10->Lo->g;
	stat[1] = t10->Hi->g;
	stat[2] = t20->Hi->g;
	stat[3] = t30->Hi->g;
	stat[4] = t41->Hi->g;
	ASSERT( stat[0]>0. && stat[1]>0. &&stat[2]>0. &&stat[3]>0. &&stat[4]>0.);
	ASSERT( fabs(t10->Lo->g/2.-1.) < FLT_EPSILON);
	ASSERT( fabs(t10->Hi->g/4.-1.) < FLT_EPSILON);
	ASSERT( fabs(t20->Lo->g/2.-1.) < FLT_EPSILON);
	ASSERT( fabs(t20->Hi->g/2.-1.) < FLT_EPSILON);
	ASSERT( fabs(t30->Lo->g/2.-1.) < FLT_EPSILON);
	ASSERT( fabs(t30->Hi->g/4.-1.) < FLT_EPSILON);
	ASSERT( fabs(t21->Lo->g/4.-1.) < FLT_EPSILON);
	ASSERT( fabs(t21->Hi->g/2.-1.) < FLT_EPSILON);
	ASSERT( fabs(t31->Lo->g/4.-1.) < FLT_EPSILON);
	ASSERT( fabs(t31->Hi->g/4.-1.) < FLT_EPSILON);
	ASSERT( fabs(t41->Lo->g/4.-1.) < FLT_EPSILON);
	ASSERT( fabs(t41->Hi->g/6.-1.) < FLT_EPSILON);

	/* excitation energy of each level relative to ground, in Kelvin */
	excit[0] = 0.;
	excit[1] = t10->EnergyK;
	excit[2] = t20->EnergyK;
	excit[3] = t30->EnergyK;
	excit[4] = t41->EnergyK + t10->EnergyK;
	ASSERT( excit[1]>0. &&excit[2]>0. &&excit[3]>0. &&excit[4]>0.);

	/* fill in Einstein As, collision strengths, pumping rates */
	AulEscp[1][0] = t10->Emis->Aul*(t10->Emis->Pesc + t10->Emis->Pelec_esc);
	AulDest[1][0] = t10->Emis->Aul*t10->Emis->Pdest;
	col_str[1][0] = t10->Coll.col_str;
	AulPump[0][1] = t10->Emis->pump;

	/* add FUV pump transitions to this pump rate */
	AulPump[0][1] += pump_rate;

	AulEscp[2][0] = t20->Emis->Aul*(t20->Emis->Pesc + t20->Emis->Pelec_esc);
	AulDest[2][0] = t20->Emis->Aul*t20->Emis->Pdest;
	col_str[2][0] = t20->Coll.col_str;
	AulPump[0][2] = t20->Emis->pump;

	AulEscp[3][0] = t30->Emis->Aul*(t30->Emis->Pesc + t30->Emis->Pelec_esc);
	AulDest[3][0] = t30->Emis->Aul*t30->Emis->Pdest;
	col_str[3][0] = t30->Coll.col_str;
	AulPump[0][3] = t30->Emis->pump;

	AulEscp[4][0] = 1e-8;/* made up trans prob */
	AulDest[4][0] = 0.;
	col_str[4][0] = cs40;
	AulPump[0][4] = 0.;

	AulEscp[2][1] = t21->Emis->Aul*(t21->Emis->Pesc + t21->Emis->Pelec_esc);
	AulDest[2][1] = t21->Emis->Aul*t21->Emis->Pdest;
	col_str[2][1] = t21->Coll.col_str;
	AulPump[1][2] = t21->Emis->pump;

	AulEscp[3][1] = t31->Emis->Aul*(t31->Emis->Pesc + t31->Emis->Pelec_esc);
	AulDest[3][1] = t31->Emis->Aul*t31->Emis->Pdest;
	col_str[3][1] = t31->Coll.col_str;
	AulPump[1][3] = t31->Emis->pump;

	AulEscp[4][1] = t41->Emis->Aul*(t41->Emis->Pesc + t41->Emis->Pelec_esc);
	AulDest[4][1] = t41->Emis->Aul*t41->Emis->Pdest;
	col_str[4][1] = t41->Coll.col_str;
	AulPump[1][4] = t41->Emis->pump;

	AulEscp[3][2] = 1e-8;/* made up trans prob */
	AulDest[3][2] = 0.;
	col_str[3][2] = cs32;
	AulPump[2][3] = 0.;

	AulEscp[4][2] = 1e-8;/* made up trans prob */
	AulDest[4][2] = 0.;
	col_str[4][2] = cs42;
	AulPump[2][4] = 0.;

	AulEscp[4][3] = 1e-8;/* made up trans prob */
	AulDest[4][3] = 0.;
	col_str[4][3] = cs43;
	AulPump[3][4] = 0.;

	lgDeBug = false;

	/* lgNegPop positive if negative pops occurred, negative if too cold */
	atom_levelN(N_SEQ_BORON,
		abundan,
		stat,
		excit,
		'K',
		pops,
		depart,
		&AulEscp,
		&col_str,
		&AulDest,
		&AulPump,
		&CollRate,
		create,
		destroy,
		false,/* say atom_levelN should evaluate coll rates from cs */
		&b_cooling,
		&dCoolDT,
		chLabel,
		&lgNegPop,
		&lgZeroPop,
		lgDeBug );/* option to print stuff - set to true for debug printout */

	/* atom_levelN did not evaluate PopLevels, so save pops here */
	/* LIMLEVELN is the dimension of the atoms vectors */
	ASSERT( N_SEQ_BORON <= LIMLEVELN);
	for( i=0; i< N_SEQ_BORON; ++i )
	{
		atoms.PopLevels[i] = pops[i];
		atoms.DepLTELevels[i] = depart[i];
	}
	/* this branch, we have a full valid solution */
	t10->Lo->Pop = pops[0];
	t20->Lo->Pop = pops[0];
	t30->Lo->Pop = pops[0];
	t21->Lo->Pop = pops[1];
	t31->Lo->Pop = pops[1];
	t41->Lo->Pop = pops[1];

	t10->Emis->PopOpc = (pops[0] - pops[1]*t10->Lo->g/t10->Hi->g);
	t20->Emis->PopOpc = (pops[0] - pops[2]*t20->Lo->g/t20->Hi->g);
	t30->Emis->PopOpc = (pops[0] - pops[3]*t30->Lo->g/t30->Hi->g);
	t21->Emis->PopOpc = (pops[1] - pops[2]*t21->Lo->g/t21->Hi->g);
	t31->Emis->PopOpc = (pops[1] - pops[3]*t31->Lo->g/t31->Hi->g);
	t41->Emis->PopOpc = (pops[1] - pops[4]*t41->Lo->g/t41->Hi->g);

	t10->Hi->Pop = pops[1];
	t20->Hi->Pop = pops[2];
	t30->Hi->Pop = pops[3];
	t21->Hi->Pop = pops[2];
	t31->Hi->Pop = pops[3];
	t41->Hi->Pop = pops[4];

	t10->Emis->phots = t10->Emis->Aul*(t10->Emis->Pesc + t10->Emis->Pelec_esc)*pops[1];
	t20->Emis->phots = t20->Emis->Aul*(t20->Emis->Pesc + t20->Emis->Pelec_esc)*pops[2];
	t30->Emis->phots = t30->Emis->Aul*(t30->Emis->Pesc + t30->Emis->Pelec_esc)*pops[3];
	t21->Emis->phots = t21->Emis->Aul*(t21->Emis->Pesc + t21->Emis->Pelec_esc)*pops[2];
	t31->Emis->phots = t31->Emis->Aul*(t31->Emis->Pesc + t31->Emis->Pelec_esc)*pops[3];
	t41->Emis->phots = t41->Emis->Aul*(t41->Emis->Pesc + t41->Emis->Pelec_esc)*pops[4];

	t10->Emis->xIntensity = t10->Emis->phots*t10->EnergyErg;
	t20->Emis->xIntensity = t20->Emis->phots*t20->EnergyErg;
	t30->Emis->xIntensity = t30->Emis->phots*t30->EnergyErg;
	t21->Emis->xIntensity = t21->Emis->phots*t21->EnergyErg;
	t31->Emis->xIntensity = t31->Emis->phots*t31->EnergyErg;
	t41->Emis->xIntensity = t41->Emis->phots*t41->EnergyErg;

	/* ratio of collisional to total excitation */
	t10->Emis->ColOvTot = CollRate[0][1]/SDIV(CollRate[0][1]+t10->Emis->pump);
	t20->Emis->ColOvTot = CollRate[0][2]/SDIV(CollRate[0][2]+t20->Emis->pump);
	t30->Emis->ColOvTot = CollRate[0][3]/SDIV(CollRate[0][3]+t30->Emis->pump);
	t21->Emis->ColOvTot = CollRate[1][2]/SDIV(CollRate[1][2]+t21->Emis->pump);
	t31->Emis->ColOvTot = CollRate[1][3]/SDIV(CollRate[1][3]+t31->Emis->pump);
	t41->Emis->ColOvTot = CollRate[1][4]/SDIV(CollRate[1][4]+t41->Emis->pump);

	/* derivative of cooling function */
	thermal.dCooldT += dCoolDT;

	/* two cases - collisionally excited (usual case) or
	 * radiatively excited - in which case line can be a heat source
	 * following are correct heat exchange, they will mix to get correct derivative
	 * the sum of heat-cool will add up to EnrUL - EnrLU - this is a trick to
	 * keep stable solution by effectively dividing up heating and cooling,
	 * so that negative cooling does not occur */

	EnrLU = t10->Lo->Pop*CollRate[0][1]*t10->EnergyErg;
	EnrUL = t10->Hi->Pop*CollRate[1][0]*t10->EnergyErg;
	/* energy exchange due to this level
	 * net cooling due to excitation minus part of de-excit */
	t10->Coll.cool = EnrLU - EnrUL*t10->Emis->ColOvTot;
	/* net heating is remainder */
	t10->Coll.heat = EnrUL*(1. - t10->Emis->ColOvTot);
	/* add to cooling */
	CoolAdd( chLabel, t10->WLAng , t10->Coll.cool);
	/* derivative of cooling function */
	thermal.dCooldT += t10->Coll.cool * (t10->EnergyK * thermal.tsq1 - thermal.halfte );

	EnrLU = t20->Lo->Pop*CollRate[0][2]*t20->EnergyErg;
	EnrUL = t20->Hi->Pop*CollRate[2][0]*t20->EnergyErg;
	t20->Coll.cool = EnrLU - EnrUL*t20->Emis->ColOvTot;
	t20->Coll.heat = EnrUL*(1. - t20->Emis->ColOvTot);
	/* add to cooling */
	CoolAdd( chLabel, t20->WLAng , t20->Coll.cool);
	thermal.dCooldT += t20->Coll.cool * (t20->EnergyK * thermal.tsq1 - thermal.halfte );

	EnrLU = t30->Lo->Pop*CollRate[0][3]*t30->EnergyErg;
	EnrUL = t30->Hi->Pop*CollRate[3][0]*t30->EnergyErg;
	t30->Coll.cool = EnrLU - EnrUL*t30->Emis->ColOvTot;
	t30->Coll.heat = EnrUL*(1. - t30->Emis->ColOvTot);
	/* add to cooling */
	CoolAdd( chLabel, t30->WLAng , t30->Coll.cool);
	thermal.dCooldT += t30->Coll.cool * (t30->EnergyK * thermal.tsq1 - thermal.halfte );

	EnrLU = t21->Lo->Pop*CollRate[1][2]*t21->EnergyErg;
	EnrUL = t21->Hi->Pop*CollRate[2][1]*t21->EnergyErg;
	t21->Coll.cool = EnrLU - EnrUL*t21->Emis->ColOvTot;
	t21->Coll.heat = EnrUL*(1. - t21->Emis->ColOvTot);
	/* add to cooling */
	CoolAdd( chLabel, t21->WLAng , t21->Coll.cool);
	/* use of 20 is intentional in following - that is Boltzmann factor */
	thermal.dCooldT += t21->Coll.cool * (t20->EnergyK * thermal.tsq1 - thermal.halfte );

	EnrLU = t31->Lo->Pop*CollRate[1][3]*t31->EnergyErg;
	EnrUL = t31->Hi->Pop*CollRate[3][1]*t31->EnergyErg;
	t31->Coll.cool = EnrLU - EnrUL*t31->Emis->ColOvTot;
	t31->Coll.heat = EnrUL*(1. - t31->Emis->ColOvTot);
	/* add to cooling */
	CoolAdd( chLabel, t31->WLAng , t31->Coll.cool);
	/* use of 30 is intentional in following - that is Boltzmann factor */
	thermal.dCooldT += t31->Coll.cool * (t30->EnergyK * thermal.tsq1 - thermal.halfte );

	EnrLU = t41->Lo->Pop*CollRate[1][4]*t41->EnergyErg;
	EnrUL = t41->Hi->Pop*CollRate[4][1]*t41->EnergyErg;
	t41->Coll.cool = EnrLU - EnrUL*t41->Emis->ColOvTot;
	t41->Coll.heat = EnrUL*(1. - t41->Emis->ColOvTot);
	/* add to cooling */
	CoolAdd( chLabel, t41->WLAng , t41->Coll.cool);
	/* use of 41 is intentional in following - that is Boltzmann factor (no 40 here) */
	thermal.dCooldT += t41->Coll.cool * (t41->EnergyK * thermal.tsq1 - thermal.halfte );

	return;
}
