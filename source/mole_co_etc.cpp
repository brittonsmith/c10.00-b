/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CO_Init called from cdInit to initialize co routines */
/*CO_update_chem_rates update rate coefficients, only temp part - in mole_co_etc.c */
#include "cddefines.h"
#include "physconst.h"
#include "mole.h"
#include "mole_co_priv.h"
#include "hmi.h"
#include "rfield.h"
#include "dense.h"
#include "ionbal.h"
#include "grainvar.h"
#include "timesc.h"
/*lint -e778 constant expression evaluates to 0 in operation '-' */

/*CO_update_chem_rates update rate coefficients, only temp part - in mole_co_etc.c 
 * called in conv_base before any chemistry or ionization is done */

enum spectype {MOLECULE, OTHER};
enum molstate {ACTIVE, PASSIVE};

STATIC void newelement(const char label[], int ipion, 
											 int priority);
STATIC struct molecule *newspecies(const char label[7], enum spectype type, 
											 enum molstate state, realnum *location, double frac0);
STATIC struct chem_element_s *findelement(const char buf[]);
STATIC int isactive(data_u *dat);
STATIC int ispassive(data_u *dat);
STATIC int isCOnet(data_u *dat);

struct chem_element_s **chem_element;
/* List of element structures indexed by atom index 
	 -- could use (element_list[nelem] != NULL) for mole.lgElem_in_chemistry[nelem] */
struct chem_element_s *element_list[LIMELM];
int32 *ipiv;
realnum *tot_ion;

struct mole_priv_s mole_priv;
struct molecule null_mole;

/*=================================================================*/
/*CO_Init called from cdInit to initialize CO routines */
void CO_Init(void)
{
	long int i,
		nelem;
	struct molecule *sp;
	static bool lgCO_Init_called=false;

	DEBUG_ENTRY( "CO_Init()" );

	/* these tell the molecular solver what zone and iteration it has
	 * been evaluated on */
	co.co_nzone = -2;
	co.iteration_co = -2;

	/* prevent memory leaks */
	/* \todo	this is a temporary fix for PR14. We should improve the overall design
	 * of this code to prevent valid pointers being overwritten in a second call to CO_Init */
	if( lgCO_Init_called )
	{
		return;
	}

	/* say that we have been called */
	lgCO_Init_called = true;

	mole_priv.spectab = newhash(NULL);
	mole_priv.reactab = newhash(NULL);
	mole_priv.elemtab = newhash(NULL);

	for(nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		element_list[nelem] = NULL;
		mole.lgElem_in_chemistry[nelem] = false;
	}

	/* set up concordance of elemental species to external Cloudy indices */
	/* final number is 'element priority':-
		 rjrw 2006 Aug 11: Nick Abel explains this order is roughly gas phase
		 abundance -- the encoding of atomic species is from the last
		 elements of "enum molecule_codes" in mole.h. Should this make any
		 difference?  In practice it does. */
	newelement("C" ,ipCARBON,3);
	newelement("^13C",ipCARBON,0);
	newelement("O" ,ipOXYGEN,2);
	newelement("Si",ipSILICON,5);
	newelement("N" ,ipNITROGEN,4);
	newelement("S" ,ipSULPHUR,6);
	newelement("Cl",ipCHLORINE,7);
	newelement("H" ,ipHYDROGEN,1);
	newelement("He",ipHELIUM,0);
	newelement("Mg",ipMAGNESIUM,0);
	newelement("Fe",ipIRON,0);

	/* set up properties of molecular species -- chemical formulae,
		 array indices, elementary components (parsed from formula), 
		 status within CO network, location of stored value external 
		 to CO network (must be floating point). */

	/* final arguments are relative to neutral carbon and are the solution to the first zone of the 
	 * TH85 pdr */

	/* Sizes of different parts of network are calculated in newspecies */
	mole.num_comole_calc = mole.num_comole_tot = mole.num_elements = 0;
	null_mole.hevmol = null_mole.hevcol = 0.;  /* Non-molecule to allow valid pointer return inactive species are accessed */
	null_mole.index = -1;

	sp = newspecies("C     ",MOLECULE,ACTIVE,NULL,1.00e+00);
	sp = newspecies("O     ",MOLECULE,ACTIVE,NULL,1.13e+05);
	sp = newspecies("Si    ",MOLECULE,ACTIVE,NULL,3.56e-04);
	sp = newspecies("N     ",MOLECULE,ACTIVE,NULL,2.25e+00);
	sp = newspecies("S     ",MOLECULE,ACTIVE,NULL,6.90e-03);
	sp = newspecies("Cl    ",MOLECULE,ACTIVE,NULL,6.90e-03);
	sp = newspecies("C+    ",MOLECULE,ACTIVE,NULL,6.79e+04);
	sp = newspecies("O+    ",MOLECULE,ACTIVE,NULL,1.58e+01);
	sp = newspecies("Si+   ",MOLECULE,ACTIVE,NULL,1.79e+02);
	sp = newspecies("N+    ",MOLECULE,ACTIVE,NULL,2.62e-10);
	sp = newspecies("S+    ",MOLECULE,ACTIVE,NULL,1.79e+03);
	sp = newspecies("Cl+   ",MOLECULE,ACTIVE,NULL,3.71e-09);
	sp = newspecies("CH    ",MOLECULE,ACTIVE,NULL,1.10e-06);
	sp = newspecies("CH+   ",MOLECULE,ACTIVE,NULL,6.99e-05);
	sp = newspecies("OH    ",MOLECULE,ACTIVE,NULL,2.15e-03);
	sp = newspecies("OH+   ",MOLECULE,ACTIVE,NULL,2.45e-04);
	sp = newspecies("O2    ",MOLECULE,ACTIVE,NULL,1.91e-07);
	sp = newspecies("CO    ",MOLECULE,ACTIVE,NULL,4.05e-05);
	sp = newspecies("^13CO ",MOLECULE,ACTIVE,NULL,4.05e-05);
	sp = newspecies("CO+   ",MOLECULE,ACTIVE,NULL,2.16e-06);
	sp = newspecies("H2O   ",MOLECULE,ACTIVE,NULL,1.60e-11);
	sp = newspecies("H2O+  ",MOLECULE,ACTIVE,NULL,1.27e-10);
	sp = newspecies("O2+   ",MOLECULE,ACTIVE,NULL,1.17e-06);
	sp = newspecies("H3O+  ",MOLECULE,ACTIVE,NULL,3.44e-17);
	sp = newspecies("CH2+  ",MOLECULE,ACTIVE,NULL,3.71e-09);
	sp = newspecies("CH2   ",MOLECULE,ACTIVE,NULL,8.59e-13);
	sp = newspecies("HCO+  ",MOLECULE,ACTIVE,NULL,4.01e-10);
	sp = newspecies("CH3+  ",MOLECULE,ACTIVE,NULL,7.02e-16);
	sp = newspecies("CH3   ",MOLECULE,ACTIVE,NULL,4.59e-19);
	sp = newspecies("CH4   ",MOLECULE,ACTIVE,NULL,9.12e-28);
	sp = newspecies("CH4+  ",MOLECULE,ACTIVE,NULL,1.03e-28);
	sp = newspecies("CH5+  ",MOLECULE,ACTIVE,NULL,8.16e-28);
	sp = newspecies("SiH2+ ",MOLECULE,ACTIVE,NULL,2.07e-13);
	sp = newspecies("SiH   ",MOLECULE,ACTIVE,NULL,1.80e-14);
	sp = newspecies("SiOH+ ",MOLECULE,ACTIVE,NULL,5.39e-15);
	sp = newspecies("SiO   ",MOLECULE,ACTIVE,NULL,3.89e-14);
	sp = newspecies("SiO+  ",MOLECULE,ACTIVE,NULL,2.27e-08);
	sp = newspecies("N2    ",MOLECULE,ACTIVE,NULL,2.46e-17);
	sp = newspecies("N2+   ",MOLECULE,ACTIVE,NULL,2.22e-17);
	sp = newspecies("NO    ",MOLECULE,ACTIVE,NULL,5.99e-12);
	sp = newspecies("NO+   ",MOLECULE,ACTIVE,NULL,1.31e-11);
	sp = newspecies("S2    ",MOLECULE,ACTIVE,NULL,1.21e-11);
	sp = newspecies("S2+   ",MOLECULE,ACTIVE,NULL,1.00e-13);
	sp = newspecies("OCN   ",MOLECULE,ACTIVE,NULL,3.86e-11);
	sp = newspecies("OCN+  ",MOLECULE,ACTIVE,NULL,1.65e-22);
	sp = newspecies("NH    ",MOLECULE,ACTIVE,NULL,1.30e-09);
	sp = newspecies("NH+   ",MOLECULE,ACTIVE,NULL,2.20e-10);
	sp = newspecies("NH2   ",MOLECULE,ACTIVE,NULL,6.78e-20);
	sp = newspecies("NH2+  ",MOLECULE,ACTIVE,NULL,1.71e-16);
	sp = newspecies("NH3   ",MOLECULE,ACTIVE,NULL,1.25e-29);
	sp = newspecies("NH3+  ",MOLECULE,ACTIVE,NULL,3.03e-23);
	sp = newspecies("NH4+  ",MOLECULE,ACTIVE,NULL,1.15e-33);
	sp = newspecies("CN    ",MOLECULE,ACTIVE,NULL,3.38e-11);
	sp = newspecies("CN+   ",MOLECULE,ACTIVE,NULL,5.07e-11);
	sp = newspecies("HCN   ",MOLECULE,ACTIVE,NULL,1.07e-17);
	sp = newspecies("HCN+  ",MOLECULE,ACTIVE,NULL,9.89e-17);
	sp = newspecies("HNO   ",MOLECULE,ACTIVE,NULL,9.09e-21);
	sp = newspecies("HNO+  ",MOLECULE,ACTIVE,NULL,1.91e-18);
	sp = newspecies("HS    ",MOLECULE,ACTIVE,NULL,1.71e-14);
	sp = newspecies("HS+   ",MOLECULE,ACTIVE,NULL,4.83e-13);
	sp = newspecies("CS    ",MOLECULE,ACTIVE,NULL,6.69e-18);
	sp = newspecies("CS+   ",MOLECULE,ACTIVE,NULL,4.12e-11);
	sp = newspecies("NO2   ",MOLECULE,ACTIVE,NULL,2.67e-26);
	sp = newspecies("NO2+  ",MOLECULE,ACTIVE,NULL,2.41e-23);
	sp = newspecies("NS    ",MOLECULE,ACTIVE,NULL,3.12e-11);
	sp = newspecies("NS+   ",MOLECULE,ACTIVE,NULL,6.81e-13);
	sp = newspecies("SO    ",MOLECULE,ACTIVE,NULL,4.00e-15);
	sp = newspecies("SO+   ",MOLECULE,ACTIVE,NULL,1.20e-07);
	sp = newspecies("SiN   ",MOLECULE,ACTIVE,NULL,7.03e-12);
	sp = newspecies("SiN+  ",MOLECULE,ACTIVE,NULL,7.60e-14);
	sp = newspecies("N2O   ",MOLECULE,ACTIVE,NULL,1.05e-11);
	sp = newspecies("HCS+  ",MOLECULE,ACTIVE,NULL,8.91e-17);
	sp = newspecies("OCS   ",MOLECULE,ACTIVE,NULL,8.83e-24);
	sp = newspecies("OCS+  ",MOLECULE,ACTIVE,NULL,1.72e-21);
	sp = newspecies("HNC   ",MOLECULE,ACTIVE,NULL,4.00e-15);
	sp = newspecies("HCNH+ ",MOLECULE,ACTIVE,NULL,4.00e-15);
	sp = newspecies("C2    ",MOLECULE,ACTIVE,NULL,3.71e-09);
	sp = newspecies("C2+   ",MOLECULE,ACTIVE,NULL,8.59e-13);
	sp = newspecies("CCl   ",MOLECULE,ACTIVE,NULL,4.00e-15);
	sp = newspecies("ClO   ",MOLECULE,ACTIVE,NULL,4.00e-15);
	sp = newspecies("HCl+  ",MOLECULE,ACTIVE,NULL,4.00e-15);
	sp = newspecies("HCl   ",MOLECULE,ACTIVE,NULL,4.00e-15);
	sp = newspecies("H2Cl+ ",MOLECULE,ACTIVE,NULL,4.00e-15);
	sp = newspecies("CCl+  ",MOLECULE,ACTIVE,NULL,4.00e-15);
	sp = newspecies("H2CCl+",MOLECULE,ACTIVE,NULL,0.00e+00);
	sp = newspecies("ClO+  ",MOLECULE,ACTIVE,NULL,4.00e-15);
	/* fprintf(stderr,"ETC: %d %d\n",gv.lgDustOn(), mole.lgGrain_mole_deplete); */
	if(gv.lgDustOn() && mole.lgGrain_mole_deplete )
	{
		sp = newspecies("COgrn ",MOLECULE,ACTIVE,NULL,1.00e-15);
		sp = newspecies("H2Ogrn",MOLECULE,ACTIVE,NULL,0.00e+00);
		sp = newspecies("OHgrn ",MOLECULE,ACTIVE,NULL,1.00e-15);
	}
	sp = newspecies("C2H   ",MOLECULE,ACTIVE,NULL,4.00e-15);
	sp = newspecies("C2H+  ",MOLECULE,ACTIVE,NULL,4.00e-15);
	sp = newspecies("C3    ",MOLECULE,ACTIVE,NULL,1.00e-15);
	sp = newspecies("C3+   ",MOLECULE,ACTIVE,NULL,4.00e-15);
	sp = newspecies("C3H+  ",MOLECULE,ACTIVE,NULL,4.00e-15);
	sp = newspecies("C3H   ",MOLECULE,ACTIVE,NULL,4.00e-15);
	sp = newspecies("C2H2  ",MOLECULE,ACTIVE,NULL,4.00e-15);
	sp = newspecies("C2H2+ ",MOLECULE,ACTIVE,NULL,4.00e-15);
	sp = newspecies("C2H3+ ",MOLECULE,ACTIVE,NULL,4.00e-15);
	sp = newspecies("N2H+  ",MOLECULE,ACTIVE,NULL,4.00e-15);
	/* Add passive species to complete network */
	sp = newspecies("e-    ",OTHER,PASSIVE,&(dense.eden_f),0.00e+00);
	sp->nElec = -1;	sp->mole_mass = (realnum)ELECTRON_MASS; /* Augment properties for this non-molecular species */

	sp = newspecies("Fe    ",MOLECULE,PASSIVE,&(dense.xIonDense[ipIRON][0]),0.00e+00);
	sp = newspecies("Fe+   ",MOLECULE,PASSIVE,&(dense.xIonDense[ipIRON][1]),0.00e+00);
	sp = newspecies("H     ",MOLECULE,PASSIVE,&(dense.xIonDense[ipHYDROGEN][0]),0.00e+00);
	sp = newspecies("H-    ",MOLECULE,PASSIVE,&(hmi.Hmolec[ipMHm]),0.00e+00);
	sp = newspecies("H+    ",MOLECULE,PASSIVE,&(dense.xIonDense[ipHYDROGEN][1]),0.00e+00);
	sp = newspecies("H2    ",MOLECULE,PASSIVE,&(hmi.Hmolec[ipMH2g]),0.00e+00);
	sp = newspecies("H2*   ",MOLECULE,PASSIVE,&(hmi.Hmolec[ipMH2s]),0.00e+00);
	sp = newspecies("H2+   ",MOLECULE,PASSIVE,&(hmi.Hmolec[ipMH2p]),0.00e+00);
	sp = newspecies("H3+   ",MOLECULE,PASSIVE,&(hmi.Hmolec[ipMH3p]),0.00e+00);
	sp = newspecies("He    ",MOLECULE,PASSIVE,&(dense.xIonDense[ipHELIUM][0]),0.00e+00);
	sp = newspecies("He+   ",MOLECULE,PASSIVE,&(dense.xIonDense[ipHELIUM][1]),0.00e+00);
	sp = newspecies("Mg    ",MOLECULE,PASSIVE,&(dense.xIonDense[ipMAGNESIUM][0]),0.00e+00);
	sp = newspecies("Mg+   ",MOLECULE,PASSIVE,&(dense.xIonDense[ipMAGNESIUM][1]),0.00e+00);

	/* Create linear list of species and populate it... */
	COmole = (struct molecule **)MALLOC((size_t)mole.num_comole_tot*
																		 sizeof(struct molecule *));

	/* ...first active species */
	i = makeplist(mole_priv.spectab,(void **)COmole,
								mole.num_comole_calc,isactive); 
	ASSERT (i == mole.num_comole_calc); 

	/* ...then passive species at end of list */
	i = makeplist(mole_priv.spectab,(void **)COmole+mole.num_comole_calc,
								mole.num_comole_tot-mole.num_comole_calc,ispassive);
	ASSERT (i == mole.num_comole_tot-mole.num_comole_calc); 

	/* Set molecule indices to order of list just created */
	for(i=0;i<mole.num_comole_tot;i++) 
	{
		COmole[i]->index = i;
	}

	/* Register the atomic ladders which are active in this network */
	for(i=0;i<mole.num_comole_calc;i++) 
	{
		sp = COmole[i];
		if(sp->active && sp->n_nuclei == 1)
		{
			if(sp->nElec == 0) 
			{
				element_list[sp->nelem_hevmol]->ipMl = i;
			}
			else if(sp->nElec == 1)
			{
				element_list[sp->nelem_hevmol]->ipMlP = i;
			}
		}
	}

	/* Create and fill cache for looping on active atomic species */
	chem_element = (struct chem_element_s **)
		MALLOC((size_t)(mole.num_elements)*sizeof(struct chem_element_s *));
	i = makeplist(mole_priv.elemtab,(void **)chem_element,mole.num_elements,isCOnet);
	ASSERT(i == mole.num_elements);


	/* Create other workspace arrays which have sizes given by the species numbers */
	mole.amat = (double **)MALLOC((size_t)mole.num_comole_calc*sizeof(double *));
	mole.amat[0] = (double *)MALLOC((size_t)mole.num_comole_calc*
															 mole.num_comole_calc*sizeof(double));
	for(i=1;i<mole.num_comole_calc;i++)
	{
		mole.amat[i] = mole.amat[i-1]+mole.num_comole_calc;
	}
	mole.b = (double *)MALLOC((size_t)mole.num_comole_calc*sizeof(double));
	mole.c = (double **)MALLOC((size_t)mole.num_comole_tot*sizeof(double *));
	mole.c[0] = (double *)MALLOC((size_t)mole.num_comole_tot*
																	mole.num_comole_calc*sizeof(double));
	for(i=1;i<mole.num_comole_tot;i++)
	{
		mole.c[i] = mole.c[i-1]+mole.num_comole_calc;
	}
	ipiv = (int32 *)MALLOC((size_t)mole.num_comole_calc*sizeof(int32));

	tot_ion = (realnum *)MALLOC((size_t)mole.num_elements*sizeof(realnum));

	return;

}
/*lint +e778 constant expression evaluates to 0 in operation '-' */

/* Fill element linking structure */
STATIC void newelement(const char label[], int ipion, int priority)
{
	struct chem_element_s *element;
	data_u *p;
	int exists;

	DEBUG_ENTRY("newelement()");

	element = (struct chem_element_s *) MALLOC(sizeof(struct chem_element_s));
	element->ipCl = ipion;
	element->ipZ = priority;
	ASSERT ((size_t)strlen(label) < sizeof(element->chName));
	strncpy(element->chName,label,sizeof(element->chName));
	element->ipMl = element->ipMlP = -1; 	/* Chemical network species indices not yet defined */
	p = addentry(element->chName,0,mole_priv.elemtab,&exists);
	p->p = (void *) element;
	element_list[ipion] = element;
}
/* Parse species string to find constituent atoms, charge etc. */
STATIC struct molecule *newspecies(const char label[7], enum spectype type, 
											 enum molstate state, realnum *location, double frac0)
{
	int exists;
	data_u *p;
	char mylab[7], thisel[CHARS_ELEMENT], *s;
	long int i, n, nnuc, nelem, nel;
	struct molecule *mol;
	struct chem_element_s *el, *maxel;

	DEBUG_ENTRY("newspecies()");

	mol = (struct molecule *) MALLOC(sizeof(struct molecule));

	mole.num_comole_tot++;
	if(state == ACTIVE)
		mole.num_comole_calc++;

	mol->pdr_mole_co = (realnum) frac0;
	for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		mol->nElem[nelem] = 0;
	}
	mol->co_save = 0.;

	strncpy(mol->label,label,sizeof(mol->label));
	s = strchr_s(mol->label,' ');
	if(s)
		*s = '\0';

	/* Insert species into hash table */
	p = addentry(mol->label,0,mole_priv.spectab,&exists);
	p->p = (void *) mol;

	if(state == ACTIVE) 
	{
		mol->active = 1;
	} 
	else 
	{
		mol->active = 0;
	}
	mol->location = location;
	mol->n_nuclei = 0;

	if(type == MOLECULE)
	{
		/* Copy to private workspace */
		strncpy(mylab,label,7);

		/* Trailing modifiers */
		mol->nElec = mol->Excit = 0;

		/* Excitation... */
		s = strpbrk(mylab,"*");
		if(s)
		{
			mol->Excit = 1;
			*s = '\0';
		} 

		/* ...Charge */
		s = strpbrk(mylab,"+-");
		if(s)
		{
			if(isdigit(*(s+1))) 
				n = atoi(s+1);
			else
				n = 1;
			if(*s == '+')
				mol->nElec = n;
			else
				mol->nElec = -n;
			*s = '\0';
		}
		/* ...Grain */
		s = strstr_s(mylab,"grn");
		if(s) 
		{
			mol->lgGas_Phase = false;
			*s = '\0';
		} 
		else 
		{
			mol->lgGas_Phase = true;
		}

		/* Now analyse chemical formula */
		nnuc = 0;  /* Keep account of number of atoms contained */
		i = 0;
		maxel = NULL;
		while (mylab[i] != '\0' && mylab[i] != ' ' && mylab[i] != '*') 
		{
			nel = 0;
		
			/* first check for isotope */	
			if( mylab[i] == '^' )
			{
				long A = 0;
				thisel[nel++] = mylab[i++];
				do {    
					A = 10*A+(long int)(mylab[i]-'0');
					thisel[nel++] = mylab[i++];
				} while (isdigit(mylab[i]));
			}
			
			/* Select next element in species, matches regexp [A-Z][a-z]? */
			thisel[nel++] = mylab[i++];
			if(islower(mylab[i])) 
			{
				thisel[nel++] = mylab[i++];
			}
			thisel[nel] = '\0';

			el = findelement(thisel);
			if(el == NULL) 
			{
				fprintf(stderr,"Did not recognize element at %s in \"%s \"[%ld]\n",mylab+i,label,i);
				cdEXIT( EXIT_FAILURE );
			}

			/* Determine 'heaviest' atom in molecular species */
			/* >> chng 06 Sep 25 rjrw: only C, O, Si, N, S, Cl and H count for this */
			if(el->ipZ != 0 && (maxel == NULL || el->ipZ > maxel->ipZ))
			{
				maxel = el;
			}

			if(isdigit(mylab[i])) /* If there is >1 of this atom */
			{
				n = 0;
				do {
					n = 10*n+(long int)(mylab[i]-'0');
					i++;
				} while (isdigit(mylab[i]));
			}
			else
			{
				n = 1;
			}
			mol->nElem[el->ipCl] += n;
			nnuc += n;
		}

		mol->n_nuclei = nnuc;

		ASSERT(mol->active == 0 || maxel != NULL );

		if(maxel != NULL)
			mol->nelem_hevmol = maxel->ipCl;
		else
			mol->nelem_hevmol = -1;

		mol->mole_mass = 0.;
		for(i=0;i<LIMELM;++i)
		{
			mol->mole_mass += 
				mol->nElem[i]*dense.AtomicWeight[i]*(realnum)ATOMIC_MASS_UNIT;
		}

		/* If we've found a neutral atomic species active in this network */
		if(mol->n_nuclei == 1 && mol->active && mol->nElec == 0) 
		{
			mole.num_elements++;
			mole.lgElem_in_chemistry[mol->nelem_hevmol] = true;
		}
	}
	return mol;
}
STATIC int isactive(data_u *dat)
{

	DEBUG_ENTRY("isactive()");

	struct molecule *p = (struct molecule *) (dat->p);

	return p->active;
}
STATIC int ispassive(data_u *dat)
{
	return !((struct molecule *) dat->p)->active;
}
STATIC int isCOnet(data_u *dat)
{
	return (((struct chem_element_s *) dat->p)->ipMl != -1);
}

struct molecule *findspecies(const char buf[])
{
	DEBUG_ENTRY("findspecies()");

	// strip string of the first space and anything after it
	string s;
	for (const char *pb = buf; *pb && *pb != ' '; ++pb)
	{
		s += *pb;
	}

	data_u *p = lookup(s.c_str(),0,mole_priv.spectab);
	
	if(p) 
		return (struct molecule *) p->p;
	else 
		return &null_mole;
}

STATIC struct chem_element_s *findelement(const char buf[])
{
	data_u *p;

	DEBUG_ENTRY("findelement()");

	p = lookup(buf,0,mole_priv.elemtab);

	if(p)  
		return (struct chem_element_s *) p->p;
	else 
		return NULL;
}

struct COmole_rate_s *CO_findrate_s(const char buf[])
{
	data_u *p;

	DEBUG_ENTRY("CO_findrate_s()");

	p = lookup(buf,0,mole_priv.reactab);

	if(p)
		return (struct COmole_rate_s *) p->p;
	else
		return NULL;
}
double CO_findrk(const char buf[])
{
	struct COmole_rate_s *rate;

	DEBUG_ENTRY("CO_findrk()");

	rate = CO_findrate_s(buf);

	if(!rate)
		return 0.;

	/* check for NaN */
	ASSERT( !isnan( rate->rk ) );

	return rate->rk;
}
double CO_findrate(const char buf[])
{
	struct COmole_rate_s *rate;
	double ret;
	int i;

	DEBUG_ENTRY("CO_findrate()");

	rate = CO_findrate_s(buf);
	if(!rate)
	{
		return 0.;
	}

	ret = rate->rk;
	for(i=0;i<rate->nrates;i++)
		ret *= rate->rate_species[i]->hevmol;
	return ret;
}

void CO_update_species_cache(void)
{
	int i;

	DEBUG_ENTRY("CO_update_species_cache()");

	dense.eden_f = (realnum)dense.eden;  /* Need floating point version for compatibility with all other values */

	for(i=0;i<mole.num_comole_tot;i++) 
	{
		if(COmole[i]->location) 
		{
			COmole[i]->hevmol = *(COmole[i]->location);

			/* check for NaN */
			ASSERT( !isnan( COmole[i]->hevmol ) );
			/* check that below limit */
			ASSERT( COmole[i]->hevmol < MAX_DENSITY );
		}
	}
}

/* Calculate rate at which molecular network abstracts species */

/* Need to check rate_species vs reactant behaviour */
double CO_sink_rate(const char chSpecies[7])
{
	int ipthis, i, n, nt;
	double ratev, ratevi;
	struct COmole_rate_s *rate;
	struct molecule *sp;

	DEBUG_ENTRY("CO_sink_rate()");

	sp = findspecies(chSpecies);
	ratev = 0;
	nt = coreactions.n;

	for(n=0;n<nt;n++) 
	{
		rate = coreactions.list[n];
		ipthis = -1;
		for(i=0;i<rate->nrates && ipthis == -1;i++)
		{
			if(rate->rate_species[i] == sp) 
			{
				ipthis = i;
			}
		}
		if(ipthis != -1) {
			ratevi = rate->rk;
			for(i=0;i<rate->nrates;i++)
			{
				if(i!=ipthis)
				{
					ratevi *= rate->rate_species[i]->hevmol;
				}
			}
			ratev += ratevi;
		}
	}	

	return ratev;
}
#ifdef PRINTREACT
STATIC void printreact(struct COmole_rate_s *rate)
{
	int i;

	DEBUG_ENTRY("printreact()");

	for(i=0;i<rate->nreactants;i++) {
		fprintf(stderr,"%s,",rate->reactants[i]->label);
	}
	fprintf(stderr,"=>");
	for(i=0;i<rate->nproducts;i++) {
		fprintf(stderr,"%s,",rate->products[i]->label);
	}
	fprintf(stderr,"\n");

}
#endif
void CO_update_rks(void)
{
	int n, nt;
	struct COmole_rate_s *rate;

	DEBUG_ENTRY("CO_update_rks()");

	nt = coreactions.n;

	for(n=0;n<nt;n++) 
	{
		rate = coreactions.list[n];
		if(rate->fun != NULL) {
			rate->rk = rate->a*rate->fun(rate);
		}
	}	
}

double CO_dissoc_rate(const char chSpecies[7])
{
	int ipthis, i, n, nt;
	double ratev, ratevi;
	struct COmole_rate_s *rate;
	struct molecule *sp;

	DEBUG_ENTRY("CO_dissoc_rate()");

	sp = findspecies(chSpecies);
	ratev = 0;
	nt = coreactions.n;

	for(n=0;n<nt;n++) 
	{
		rate = coreactions.list[n];
		if(rate->photon == -1) 
		{
			ipthis = 0;
			for(i=0;i<rate->nproducts;i++)
			{
				if(rate->products[i] == sp) 
				{
					ipthis++;
				}
			}
			if(ipthis) 
			{
				ratevi = rate->rk;
				for(i=0;i<rate->nrates;i++)
				{
					ratevi *= rate->rate_species[i]->hevmol;
				}
				ratev += ipthis*ratevi;
			}
		}
	}	

	return ratev;
}
double CO_source_rate(const char chSpecies[7])
{
	int ipthis, i, n, nt;
	double ratev, ratevi;
	struct COmole_rate_s *rate;
	struct molecule *sp;

	DEBUG_ENTRY("CO_source_rate()");

	sp = findspecies(chSpecies);
	ratev = 0;
	nt = coreactions.n;

	for(n=0;n<nt;n++) 
	{
		rate = coreactions.list[n];
		ipthis = 0;
		for(i=0;i<rate->nproducts;i++)
		{
			if(rate->products[i] == sp) 
			{
				ipthis++;
			}
		}
		if(ipthis) {
			ratevi = rate->rk;
			for(i=0;i<rate->nrates;i++)
			{
				ratevi *= rate->rate_species[i]->hevmol;
			}
			ratev += ipthis*ratevi;
		}
	}	

	return ratev;
}

/* Punch all rates involving specified species */
void CO_punch_mol(FILE *punit, const char chSpecies[], char header[], double depth)
{
	int n, nt, i, ipthis;
	double ratevi;
	struct COmole_rate_s *rate;
	struct molecule *sp;
	char *s;

	DEBUG_ENTRY("CO_punch_mol()");

	s = header;

	sp = findspecies(chSpecies);
	if(punit == NULL)
	{
		if (sp == &null_mole)
		{
			fprintf(ioQQQ,"Could not find species '%s' to save chemistry rates.\nSorry.\n",
					  chSpecies);
			cdEXIT(EXIT_FAILURE);
		}
		sprintf (s,"#Depth");
		s += strlen(s);
	}
	else
	{
		fprintf (punit,"%.5e",depth);
	}

	nt = coreactions.n;
	for(n=0;n<nt;n++) {
		rate = coreactions.list[n];
		ipthis = 0;
		for(i=0;i<rate->nrates;i++)
		{
			if(rate->rate_species[i] == sp) 
			{
				ipthis++;
			}
		}
		for(i=0;i<rate->nproducts;i++)
		{
			if(rate->products[i] == sp) 
			{
				ipthis++;
			}
		}
		if(ipthis) 
		{
			if(punit == NULL)
			{
				sprintf(s,"\t%s",rate->label);
				s += strlen(s);
			}
			else
			{
				ratevi = rate->rk;
				for(i=0;i<rate->nrates;i++)
				{
					ratevi *= rate->rate_species[i]->hevmol;
				}				
				fprintf(punit,"\t%.3e",ratevi);
			}
		}
	}
	if(punit == NULL)
	{
		sprintf(s,"\n");
	}
	else
	{
		fprintf(punit,"\n");
	}
}	

void CO_zero(void)
{
	long int i;

	static bool lgFirstCall = true;
	static long int num_comole_calc_MALLOC=-1;

	DEBUG_ENTRY("CO_zero()");

	if( lgFirstCall )
	{
		/* do one-time malloc of timesc.AgeCOMoleDest */
		timesc.AgeCOMoleDest = (double*)MALLOC( mole.num_comole_calc*sizeof(double) );

		lgFirstCall = false;
		num_comole_calc_MALLOC = mole.num_comole_calc;
	}
	else if( mole.num_comole_calc>num_comole_calc_MALLOC )
	{
		/* number of species has increased since last time - this can't happen
		 * tsuite / programs / comp4 has 95 first time, 98 second time */
		fprintf(ioQQQ,"DISASTER - the number of species in the CO network has increased.  This is not allowed.\n");
		fprintf(ioQQQ,"This could happen if an element was initially turned off or grains not included, then the element or grains was included.  There are not allowed.\n");
		fprintf(ioQQQ,"Sorry.\n");
		cdEXIT(EXIT_FAILURE);
	}

	for( i=0; i < mole.num_comole_calc; i++ )
	{
		timesc.AgeCOMoleDest[i] = 0.;
	}

	/* largest fraction of atoms in molecules */
	for( i=0; i<mole.num_comole_calc; ++i )
		COmole[i]->xMoleFracMax = 0.;

	/* zero out molecular species */
	for( i=0; i < mole.num_comole_tot; i++ )
	{
		COmole[i]->hevmol = 0.;
		COmole[i]->hevcol = 0.;
	}
}
