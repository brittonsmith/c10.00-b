#include "cddefines.h"
#include "mole.h"
#include "mole_co_priv.h"
#include "hmi.h"
#include "conv.h"
#include "grainvar.h"

#ifdef _MSC_VER
#	pragma warning( disable : 4100 )/* unreferenced formal parameter  */
#endif

/*
 * HOWTO:- add a reaction to the new CO network, as at 2006 December 07.
 *
 * add a line of the form
 *   newreact("O,C2=>CO,C",hmrate,SCALE,B,C); / * Data source * /
 * to this file for the new reaction.  The first argument is the chemical
 * reaction, the second is the function which is used to evaluate the 
 * rate coefficients.
 *
 * SCALE is an overall constant by which the reaction rate is scaled,
 * the remaining arguments constants used by this function.
 *
 *
 * If all the species have previously been defined, and the parser in
 * mole_co_etc.c can understand the reaction string, then that's it.
 * The sources and sinks to other networks will need to be defined, 
 * though.
 *
 */

/* Nick Abel, 06 Nov 27 states: "In all the crnu reactions, the factor of two is the albedo factor 1/(1-albedo)."
	 Can this be obtained from the grain physics? */

/* Structures containing reaction data */
struct t_coreactions coreactions;

STATIC void newreact(const char label[], 
										 double (*fun)(struct COmole_rate_s *rate), double a, double b, double c);
STATIC double hmrate(struct COmole_rate_s *rate);
STATIC double constrate(struct COmole_rate_s *rate);
STATIC double th85rate(struct COmole_rate_s *rate);
STATIC double crnurate(struct COmole_rate_s *rate);
STATIC double co_lnu_c_o_lnu(struct COmole_rate_s *rate);
STATIC double ele_ion_ladder(struct COmole_rate_s *rate);
STATIC double vib_evap(struct COmole_rate_s *rate);
STATIC double th85rate_co(struct COmole_rate_s *rate);
STATIC double grn_abs(struct COmole_rate_s *rate);
STATIC double oh_c2h2_co_ch3(struct COmole_rate_s *rate);
STATIC double h_hnc_hcn_h(struct COmole_rate_s *rate);

static bool lgReactInitialized = false;

void CO_create_react( void )
{
	/* Should adaptively select reactions rather than list them explicitly here */

	/* prevent memory leaks */
	/* \todo	this is a temporary fix for PR14. We should improve the overall design
	 * of this code to prevent valid pointers being overwritten in a second call to CO_create_react */
	if( lgReactInitialized )
		return;

	lgReactInitialized = true;

	long int i;
	/* Initialize number of reactions added to list in coreactions of reactions
	   treated by new methodology, incremented by newreact() */

	DEBUG_ENTRY("CO_create_react()");

	coreactions.n = 0; 

	newreact("N+,(e-)=>N",ele_ion_ladder,1.,0.,0.);/*   */
	newreact("S+,(e-)=>S",ele_ion_ladder,1.,0.,0.);/*   */
	newreact("Cl+,(e-)=>Cl",ele_ion_ladder,1.,0.,0.);/*   */
	newreact("C+,(e-)=>C",ele_ion_ladder,1.,0.,0.);/*   */
	newreact("O+,(e-)=>O",ele_ion_ladder,1.,0.,0.);/*   */
	newreact("Si+,(e-)=>Si",ele_ion_ladder,1.,0.,0.);/*   */
	newreact("C,PHOTON=>C+,e-",ele_ion_ladder,1.,0.,0.);/*   */
	newreact("O,PHOTON=>O+,e-",ele_ion_ladder,1.,0.,0.);/*   */
	newreact("Si,PHOTON=>Si+,e-",ele_ion_ladder,1.,0.,0.);/*   */
	newreact("Cl,PHOTON=>Cl+,e-",ele_ion_ladder,1.,0.,0.);/*   */
	newreact("N,PHOTON=>N+,e-",ele_ion_ladder,1.,0.,0.);/*   */
	newreact("S,PHOTON=>S+,e-",ele_ion_ladder,1.,0.,0.);/*   */
	/* >>chng 06 Feb 28 -- NPA.  Charge transfer between S+, Si+, and C+ and Mg/Fe is 
	 * sometimes important in molecular abundance determinations.  We therefore include them here.  Someday, we may need to 
	 * include Mg and Fe in the molecular network, but for now include the reaction in the existing heavy element
	 * network */
	newreact("S+,Mg=>S,Mg+",hmrate,2.8e-10,0.,0.); /* TH85	 */
	newreact("S+,Fe=>S,Fe+",hmrate,1.8e-10,0.,0.); /* TH85	 */
	newreact("Si+,Mg=>Si,Mg+",hmrate,2.9e-10,0.,0.); /* TH85	 */
	newreact("Si+,Fe=>Si,Fe+",hmrate,1.9e-10,0.,0.); /* TH85	 */
	newreact("C+,Mg=>C,Mg+",hmrate,1.1e-9,0.,0.); /* TH85	 */
	newreact("C+,Fe=>C,Fe+",hmrate,2.6e-9,0.,0.); /* TH85	 */
	if(co.lgUMISTrates) 
		newreact("CO,lnu=>C,O,lnu",co_lnu_c_o_lnu,1.,0.,0.);/*  Inner shell photoionization???  */
	newreact("PHOTON,CO=>C,O",th85rate_co,2.0e-10,3.2,0.);/*  UMIST	 */
	if(!co.lgUMISTrates) 
		newreact("C,CRPHOT=>C+,e-",crnurate,2*510,0.,0.); /*  */
	if( gv.lgDustOn() && mole.lgGrain_mole_deplete )
	{
		/* fprintf(stderr,"REA: %d %d\n",gv.lgDustOn(), mole.lgGrain_mole_deplete); */
		newreact("COgrn=>CO,grn",vib_evap,1.,1210.,0.);/*   */
		newreact("COgrn,CRPHOT=>CO,grn",crnurate,196./2.17,0.,0.); /*  */
		newreact("CO,grn=>COgrn",grn_abs,1.,0.,0.);/*   */	/* >>chng 06 mar 05, include missing factor 8/PI in average speed of molecules, PvH */
		newreact("H2Ogrn=>H2O,grn",vib_evap,1.,1860.,0.);/*   */
		newreact("H2Ogrn,CRPHOT=>H2O,grn",crnurate,0.028/2.17,0.,0.); /*  */
		newreact("H2O,grn=>H2Ogrn",grn_abs,1.,0.,0.);/*   */
		newreact("OHgrn=>OH,grn",vib_evap,1.,1260.,0.);/*   */
		newreact("OHgrn,CRPHOT=>OH,grn",crnurate,126./2.17,0.,0.); /*  */
		newreact("OH,grn=>OHgrn",grn_abs,1.,0.,0.);/*   */
	}
	newreact("H,CO+=>CO,H+",hmrate,7.5e-10,0,0); /* UMIST	 */
	newreact("H-,HCO+=>CO,H2",hmrate,0.00000023,-0.5,0); /* UMIST	 */
	newreact("H2+,CO=>HCO+,H",hmrate,2.16e-9,0,0); /* UMIST	 */
	newreact("H2+,CO=>CO+,H2",hmrate,6.4e-10,0,0); /* UMIST	 */
	newreact("H3+,CO=>HCO+,H2",hmrate,1.7e-9,0,0); /* UMIST	 */
	newreact("He+,CO=>O+,C,He",hmrate,0.00000000000000014,-0.5,0); /* UMIST	 */
	newreact("He+,CO=>O,C+,He",hmrate,1.6e-9,0,0); /* UMIST	 */
	newreact("CRPHOT,CO=>C,O",crnurate,2. * 10,0.,0.); /* UMIST	 */
	newreact("CRP,CO=>CO+,e-",hmrate,0.000000000000000039,0.,0.); /* UMIST	 */
	newreact("C,CO=>C2,O",hmrate,2.94e-11,0.5,58025); /* UMIST	 */
	newreact("O,C2=>CO,C",hmrate,5.e-11,0.5,0); /* UMIST	 */
	newreact("C2,O2=>CO,CO",hmrate,1.83e-12,0,0); /* UMIST	 */
	newreact("C2,O2+=>CO+,CO",hmrate,4.1e-10,0,0); /* UMIST	 */
	newreact("C2,CO+=>CO,C2+",hmrate,8.4e-10,0,0); /* UMIST	 */
	newreact("C2+,O2=>CO+,CO",hmrate,8.e-10,0,0); /* UMIST	 */
	newreact("C,CO+=>CO,C+",hmrate,1.1e-10,0,0); /* UMIST	 */
	newreact("C,HCO+=>CO,CH+",hmrate,1.1e-9,0,0); /* UMIST	 */
	newreact("C,O=>CO,PHOTON",hmrate,2.1E-19,0,0); /* UMIST	 */
	newreact("C,O2=>CO,O",hmrate,3.3e-11,0,0); /* UMIST	 */
	newreact("C,OH=>CO,H",hmrate,1.1e-10,0.5,0); /* UMIST	 */
	newreact("C,SiO+=>Si+,CO",hmrate,1.0e-9,0,0); /* UMIST	 */
	/* >>chng 05 aug 17, NA, need to use UMIST rate for PDR comparison, even though the rate is 
	 * an order of magnitude too small at 50 K */
	if(co.lgProtElim)
	{
		if(!co.lgUMISTrates)
			newreact("C+,OH=>CO,H+",hmrate,7.7e-10,0.,0.); /* UMIST; Dubernet et al. 1992, ApJ, 239, 855	 */
		else
			newreact("C+,OH=>CO,H+",hmrate,2.7e-9,-0.3508,0); /* UMIST; Dubernet et al. 1992, ApJ, 239, 855	 */
	}
	if(!co.lgUMISTrates)
		newreact("C+,OH=>CO+,H",hmrate,7.7e-10,0.,0.);/*  UMIST; Dubernet et al. 1992, ApJ, 239, 855	 */
	else
		newreact("C+,OH=>CO+,H",hmrate,2.7e-9,-0.3508,0); /* UMIST; Dubernet et al. 1992, ApJ, 239, 855	 */
	newreact("C+,SiO=>Si+,CO",hmrate,5.4e-10,0,0); /* UMIST	 */
	newreact("C+,O2=>CO,O+",hmrate,6.2e-10,0,0); /* UMIST	 */
	newreact("O,CH2=>CO,H,H",hmrate,1.33e-10,0,0); /* UMIST	 */
	newreact("O,CH2=>CO,H2",hmrate,8.e-11,0,0); /* UMIST	 */
	newreact("O,CO+=>CO,O+",hmrate,1.4e-10,0,0); /* UMIST	 */
	newreact("O+,CO=>CO+,O",hmrate,4.9e-12,0.5,4580); /* UMIST	 */
	newreact("CH,CO+=>CO,CH+",hmrate,3.2e-10,0,0); /* UMIST	 */
	newreact("CH,HCO+=>CO,CH2+",hmrate,6.3e-10,0,0); /* UMIST	 */
	newreact("CH,O2=>CO,OH",hmrate,2.6e-11,0,0); /* UMIST	 */
	newreact("CH2,CO+=>CO,CH2+",hmrate,4.3e-10,0,0); /* UMIST	 */
	newreact("CH2,HCO+=>CO,CH3+",hmrate,8.6e-10,0,0); /* UMIST	 */
	newreact("CH2,O2=>CO,H2O",hmrate,2.48e-10,-3.3,1443); /* UMIST	 */
	newreact("CO+,O2=>O2+,CO",hmrate,1.2e-10,0,0); /* UMIST	 */
	newreact("H2O,CO+=>CO,H2O+",hmrate,1.72e-9,0,0); /* UMIST	 */
	newreact("H2O,HCO+=>CO,H3O+",hmrate,2.5e-9,0,0); /* UMIST	 */
	newreact("H2O+,CO=>HCO+,OH",hmrate,5.0e-10,0,0); /* UMIST	 */
	newreact("HCO+,SiH=>SiH2+,CO",hmrate,8.7e-10,0,0); /* UMIST	 */
	newreact("HCO+,SiO=>SiOH+,CO",hmrate,7.9e-10,0,0); /* UMIST	 */
	newreact("OH,CO+=>CO,OH+",hmrate,3.1e-10,0,0); /* UMIST	 */
	newreact("OH,HCO+=>CO,H2O+",hmrate,6.2e-10,0,0); /* UMIST	 */
	newreact("OH+,CO=>HCO+,O",hmrate,1.05e-9,0,0); /* UMIST	 */
	newreact("CO+,CH4=>CO,CH4+",hmrate,7.93e-10,0,0); /* UMIST     */
	newreact("CO,CH4+=>HCO+,CH3",hmrate,1.4e-9,0,0); /* UMIST     */
	newreact("CO,CH5+=>HCO+,CH4",hmrate,1.0e-9,0,0); /* UMIST     */
	newreact("C,NO=>CO,N",hmrate,4.65e-11,0,0); /* UMIST   */
	newreact("C,OCN=>CO,CN",hmrate,4.0e-11,0.5,0); /* UMIST   */
	newreact("C,SO=>S,CO",hmrate,7.2e-11,0,0); /* UMIST   */
	newreact("O,CN=>CO,N",hmrate,4.36e-11,0.46,364); /* UMIST   */
	newreact("O,HCN=>CO,NH",hmrate,0.00000000000073,1.14,3742); /* UMIST   */
	newreact("O,OCN=>NO,CO",hmrate,9.43e-11,-0.09,100); /* UMIST   */
	newreact("O,CS=>S,CO",hmrate,2.48e-10,-0.65,783); /* UMIST   */
	newreact("O,OCS=>SO,CO",hmrate,1.6e-11,0,2150); /* UMIST   */
	newreact("OH,HCN=>CO,NH2",hmrate,1.07e-13,0,5892); /* UMIST   */
	newreact("CN,NO=>N2,CO",hmrate,1.79e-10,0,4040); /* UMIST   */
	newreact("CN,O2=>NO,CO",hmrate,0.00000000000053,0,0); /* UMIST   */
	newreact("CO,HS=>OCS,H",hmrate,0.0000000000000595,1.12,8330); /* UMIST   */
	newreact("H+,OCS=>HS+,CO",hmrate,2.1e-9,0,0); /* UMIST   */
	newreact("He+,OCS=>S+,CO,He",hmrate,7.6e-10,0,0); /* UMIST   */
	newreact("C+,SO=>S+,CO",hmrate,2.6e-10,0,0); /* UMIST   */
	newreact("C+,OCS=>CS+,CO",hmrate,1.6e-9,0,0); /* UMIST   */
	newreact("CH+,OCS=>HCS+,CO",hmrate,1.05e-9,0,0); /* UMIST   */
	newreact("N+,CO=>NO+,C",hmrate,1.45e-10,0,0); /* UMIST   */
	newreact("N+,OCS=>S+,CO,N",hmrate,3.08e-10,0,0); /* UMIST   */
	newreact("NH+,CO=>HCO+,N",hmrate,4.41e-10,0,0); /* UMIST   */
	newreact("NH+,CO=>OCN+,H",hmrate,5.39e-10,0,0); /* UMIST   */
	newreact("NH,HCO+=>CO,NH2+",hmrate,6.4e-10,0,0); /* UMIST   */
	newreact("NH2,HCO+=>CO,NH3+",hmrate,8.9e-10,0,0); /* UMIST   */
	newreact("NH3,HCO+=>CO,NH4+",hmrate,2.2e-9,0,0); /* UMIST   */
	newreact("CN+,O2=>NO+,CO",hmrate,8.6e-11,0,0); /* UMIST   */
	newreact("HCN+,CO=>HCO+,CN",hmrate,1.4e-10,0,0); /* UMIST   */
	newreact("CO,HNO+=>NO,HCO+",hmrate,1.0e-10,0,0); /* UMIST   */
	newreact("N2+,OCS=>S+,N2,CO",hmrate,1.04e-9,0,0); /* UMIST   */
	newreact("HCO+,S=>HS+,CO",hmrate,3.3e-10,0,0); /* UMIST   */
	newreact("HCO+,CS=>HCS+,CO",hmrate,1.2e-9,0,0); /* UMIST   */
	newreact("NH,CO+=>CO,NH+",hmrate,3.2e-10,0,0); /* UMIST   */
	newreact("NH2,CO+=>CO,NH2+",hmrate,4.5e-10,0,0); /* UMIST   */
	newreact("NH3,CO+=>CO,NH3+",hmrate,2.02e-9,0,0); /* UMIST   */
	newreact("CN+,CO=>CO+,CN",hmrate,6.3e-10,0,0); /* UMIST   */
	newreact("HCN,CO+=>CO,HCN+",hmrate,3.4e-9,0,0); /* UMIST   */
	newreact("CO,N2+=>N2,CO+",hmrate,7.4e-11,0,0); /* UMIST   */
	newreact("CO+,NO=>NO+,CO",hmrate,3.3e-10,0,0); /* UMIST   */
	newreact("OCN+,e-=>CO,N",hmrate,0.0000003,-0.5,0); /* UMIST   */
	newreact("OCS+,e-=>S,CO",hmrate,0.00000015,-0.5,0); /* UMIST   */
	newreact("CO,S=>OCS,PHOTON",hmrate,1.6e-17,-1.5,0); /* UMIST   */
	newreact("OCS,PHOTON=>S,CO",th85rate,3.7e-9,0.,0.); /* UMIST   */
	newreact("OCS,CRP=>S,CO",crnurate,5360 * 2,0.,0.); /* UMIST   */
	newreact("N+,CO=>CO+,N",hmrate,8.25e-10,0,0); /* UMIST   */
	newreact("CO+,S=>S+,CO",hmrate,1.1e-9,0,0); /* UMIST   */
	newreact("O,CCl=>Cl,CO",hmrate,0.0000000000996,0,0); /* UMIST   */
	newreact("CO,H2Cl+=>HCl,HCO+",hmrate,0.00000000078,0,0); /* UMIST   */
	newreact("HNC,HCO+=>HCNH+,CO",hmrate,0.0000000031,0,0); /*  */
	newreact("HCN,HCO+=>HCNH+,CO",hmrate,0.0000000031,0,0); /*  */
	newreact("O,C2H=>CO,CH",hmrate,0.000000000017,0,0); /*  */
	newreact("C2H,CO+=>CO,C2H+",hmrate,0.00000000039,0,0); /*  */
	newreact("C2,HCO+=>CO,C2H+",hmrate,0.00000000083,0,0); /*  */
	newreact("O2,C3=>CO,C2,O",hmrate,0.000000000001,0,0); /*  */
	newreact("O,C3=>CO,C2",hmrate,0.00000000005,0.5,0); /*  */
	newreact("C2H,HCO+=>CO,C2H2+",hmrate,0.00000000078,0,0); /*  */
	newreact("O,C3H=>CO,C2H",hmrate,0.0000000001,0,250); /*  */
	newreact("O,C2H2=>CO,CH2",hmrate,0.00000000000839,1.03,1197); /*  */
	newreact("OH,C2H2=>CO,CH3",oh_c2h2_co_ch3,6.51E-18,4.,-1006);/*   */
	newreact("HCO+,C3=>C3H+,CO",hmrate,0.000000002,0,0); /*  */
	newreact("H2O,C3H+=>CO,C2H3+",hmrate,0.00000000018,0,0); /*  */
	newreact("C2H2,HCO+=>CO,C2H3+",hmrate,0.0000000014,0,0); /*  */
	newreact("e-,CH+=>C,H",hmrate,0.00000015,-0.42,0); /* UMIST	 */
	newreact("e-,CH2+=>C,H,H",hmrate,0.000000403,-0.6,0); /* UMIST	 */
	newreact("e-,CH3+=>CH,H,H",hmrate,0.000000056,-0.5,0); /* UMIST	 */
	newreact("e-,CH3+=>CH3,PHOTON",hmrate,1.1e-10,-0.5,0); /* UMIST	 */

	newreact("e-,H2O+=>OH,H",hmrate,0.0000000792,-0.5,0); /* UMIST	 */
	newreact("e-,H2O+=>O,H,H",hmrate,0.000000245,-0.5,0); /* UMIST	 */
	newreact("e-,H2O+=>O,H2",hmrate,0.000000036,-0.5,0); /* UMIST	 */
	newreact("e-,H3O+=>H2O,H",hmrate,0.000000108,-0.5,0); /* UMIST	 */
	newreact("e-,H3O+=>OH,H,H",hmrate,0.000000258,-0.5,0); /* UMIST	 */
	newreact("e-,H3O+=>OH,H2",hmrate,0.0000000645,-0.5,0); /* UMIST	 */
	newreact("e-,H3O+=>O,H2,H",hmrate,5.59e-9,-0.5,0); /* UMIST	 */
	newreact("e-,O2+=>O,O",hmrate,0.000000195,-0.7,0); /* UMIST	 */
	newreact("e-,OH+=>O,H",hmrate,0.0000000375,-0.5,0); /* UMIST	 */
	newreact("e-,SiH2+=>SiH,H",hmrate,0.00000015,-0.5,0); /* UMIST	 */
	newreact("e-,SiH2+=>Si,H,H",hmrate,0.0000002,-0.5,0); /* UMIST	 */
	newreact("e-,SiH2+=>Si,H2",hmrate,0.00000015,-0.5,0); /* UMIST	 */
	newreact("e-,SiO+=>Si,O",hmrate,0.0000002,-0.5,0); /* UMIST	 */
	newreact("e-,SiOH+=>SiO,H",hmrate,0.00000015,-0.5,0); /* UMIST	 */
	newreact("e-,SiOH+=>Si,OH",hmrate,0.00000015,-0.5,0); /* UMIST	 */
	newreact("e-,CH5+=>CH3,H2",hmrate,0.00000055,-0.3,0); /* UMIST     */
	newreact("e-,CH5+=>CH4,H",hmrate,0.00000055,-0.3,0); /* UMIST     */
	newreact("e-,CH4+=>CH3,H",hmrate,0.000000175,-0.5,0); /* UMIST     */
	newreact("e-,CH4+=>CH2,H,H",hmrate,0.000000175,-0.5,0); /* UMIST     */
	newreact("C2+,e-=>C,C",hmrate,0.0000003,-0.5,0); /* UMIST	 */
	newreact("NH+,e-=>N,H",hmrate,0.000000043,-0.5,0); /* UMIST   */
	newreact("NH2+,e-=>N,H,H",hmrate,0.000000198,-0.5,0); /* UMIST   */
	newreact("NH2+,e-=>NH,H",hmrate,0.000000102,-0.5,0); /* UMIST   */
	newreact("NH3+,e-=>NH,H,H",hmrate,0.000000155,-0.5,0); /* UMIST   */
	newreact("NH3+,e-=>NH2,H",hmrate,0.000000155,-0.5,0); /* UMIST   */
	newreact("NH4+,e-=>NH2,H,H",hmrate,0.000000286,-0.5,0); /* UMIST   */
	newreact("NH4+,e-=>NH2,H2",hmrate,0.000000137,-0.5,0); /* UMIST   */
	newreact("NH4+,e-=>NH3,H",hmrate,0.000000938,-0.5,0); /* UMIST   */
	newreact("CN+,e-=>N,C",hmrate,0.00000018,-0.5,0); /* UMIST   */
	newreact("HCN+,e-=>CN,H",hmrate,0.0000002,-0.5,0); /* UMIST   */
	newreact("N2+,e-=>N,N",hmrate,0.000000036,-0.42,0); /* UMIST   */
	newreact("NO+,e-=>O,N",hmrate,0.00000043,-0.37,0); /* UMIST   */
	newreact("HNO+,e-=>NO,H",hmrate,0.0000003,-0.5,0); /* UMIST   */
	newreact("HS+,e-=>S,H",hmrate,0.0000002,-0.5,0); /* UMIST   */
	newreact("SiN+,e-=>Si,N",hmrate,0.0000002,-0.5,0); /* UMIST   */
	newreact("CS+,e-=>S,C",hmrate,0.0000002,-0.5,0); /* UMIST   */
	newreact("HCS+,e-=>CS,H",hmrate,0.00000005,-0.75,0); /* UMIST   */
	newreact("NO2+,e-=>NO,O",hmrate,0.0000003,-0.5,0); /* UMIST   */
	newreact("NS+,e-=>S,N",hmrate,0.0000002,-0.5,0); /* UMIST   */
	newreact("SO+,e-=>S,O",hmrate,0.0000002,-0.5,0); /* UMIST   */
	newreact("OCS+,e-=>CS,O",hmrate,0.00000015,-0.5,0); /* UMIST   */
	newreact("S2+,e-=>S,S",hmrate,0.0000002,-0.5,0); /* UMIST   */
	newreact("HCl+,e-=>Cl,H",hmrate,0.0000003,-0.5,0); /* UMIST   */
	newreact("H2Cl+,e-=>Cl,H,H",hmrate,0.00000027,-0.5,0); /* UMIST   */
	newreact("H2Cl+,e-=>HCl,H",hmrate,0.00000003,-0.5,0); /* UMIST   */
	newreact("CCl+,e-=>Cl,C",hmrate,0.0000003,-0.5,0); /* UMIST   */
	newreact("H2CCl+,e-=>CCl,H,H",hmrate,0.0000003,-0.5,0); /* UMIST   */
	newreact("ClO+,e-=>Cl,O",hmrate,0.0000002,-0.5,0); /* UMIST   */
	newreact("HCNH+,e-=>HCN,H",hmrate,0.00000009,-0.5,0); /*  */
	newreact("HCNH+,e-=>CN,H,H",hmrate,0.00000018,-0.5,0); /*  */
	newreact("HCNH+,e-=>HNC,H",hmrate,0.00000009,-0.5,0); /*  */
	newreact("C2H+,e-=>CH,C",hmrate,0.000000135,-0.5,0); /*  */
	newreact("C2H+,e-=>C2,H",hmrate,0.000000116,-0.76,0); /* Herbst */
	newreact("C3+,e-=>C2,C",hmrate,0.0000003,-0.5,0); /*  */
	newreact("C2H2+,e-=>C2,H,H",hmrate,0.00000009,-0.5,0); /*  */
	newreact("C2H2+,e-=>CH,CH",hmrate,0.00000009,-0.5,0); /*  */
	newreact("C2H2+,e-=>C2H,H",hmrate,0.00000029,-0.5,0); /* Herbst */
	newreact("C3H+,e-=>C2H,C",hmrate,0.00000015,-0.5,0); /*  */
	newreact("C3H+,e-=>C3,H",hmrate,0.00000015,-0.5,0); /*  */
	newreact("C2H3+,e-=>C2H2,H",hmrate,0.00000045,-0.5,0); /*  */
	/* >>chng 06 Apr 13, add N2H+ to chemistry.  which 
	 * should improve modeling of nitrogen chemistry, in 
	 * particular NH and N2 */
	newreact("N2H+,e-=>N2,H",hmrate,0.36e-7,-0.51,0); /*   */
	newreact("N2H+,e-=>NH,N",hmrate,0.64e-7,-0.51,0); /*   */
	newreact("NH,PHOTON=>NH+,e-",th85rate,1.0e-11,0.,0.); /* UMIST   */
	newreact("NH,PHOTON=>N,H",th85rate,5.0e-10,0.,0.); /* UMIST   */
	newreact("NH2,PHOTON=>NH,H",th85rate,3.9e-10,0.,0.); /* UMIST   */
	newreact("NH2,PHOTON=>NH2+,e-",th85rate,1.73e-10,0.,0.); /* UMIST   */
	newreact("NH3,PHOTON=>NH3+,e-",th85rate,1.8e-10,0.,0.); /* UMIST   */
	newreact("NH3,PHOTON=>NH,H2",th85rate,3.3e-10,0.,0.); /* UMIST   */
	newreact("NH3,PHOTON=>NH2,H",th85rate,1.1e-9,0.,0.); /* UMIST   */
	newreact("CN,PHOTON=>N,C",th85rate,1.1e-9,0.,0.);/*  UMIST   */	/* >>chng 06 jun 01, comment out "??" - in more than one time*/
	newreact("HCN,PHOTON=>CN,H",th85rate,1.3e-9,0.,0.); /* UMIST   */
	newreact("N2,PHOTON=>N,N",th85rate,2.3e-10,0.,0.); /* UMIST   */
	newreact("NO,PHOTON=>O,N",th85rate,4.3e-10,0.,0.); /* UMIST   */
	newreact("NO,PHOTON=>NO+,e-",th85rate,2.6e-10,0.,0.); /* UMIST   */
	newreact("HNO,PHOTON=>NO,H",th85rate,1.7e-10,0.,0.); /* UMIST   */
	newreact("HS,PHOTON=>S,H",th85rate,9.7e-10,0.,0.); /* UMIST   */
	newreact("HS+,PHOTON=>S+,H",th85rate,3.e-10,0.,0.); /* UMIST   */
	newreact("OCN,PHOTON=>CN,O",th85rate,1.0e-11,0.,0.); /* UMIST   */
	newreact("CS,PHOTON=>S,C",th85rate,9.7e-10,0.,0.); /* UMIST   */
	newreact("CS+,PHOTON=>S+,C",th85rate,2.0e-10,0.,0.); /* UMIST   */
	newreact("CS,PHOTON=>CS+,e-",th85rate,2.0e-10,0.,0.); /* UMIST   */
	newreact("NO2,PHOTON=>NO,O",th85rate,1.2e-9,0.,0.); /* UMIST   */
	newreact("NS,PHOTON=>S,N",th85rate,1.0e-11,0.,0.); /* UMIST   */
	newreact("SO,PHOTON=>S,O",th85rate,3.7e-9,0.,0.); /* UMIST   */
	newreact("SO,PHOTON=>SO+,e-",th85rate,6.e-10,0.,0.); /* UMIST   */
	newreact("OCS,PHOTON=>OCS+,e-",th85rate,4.2e-10,0.,0.); /* UMIST   */
	newreact("NH,CRP=>NH+,e-",crnurate,500 * 2,0.,0.); /* UMIST   */
	newreact("NH,CRP=>N,H",crnurate,500 * 2,0.,0.); /* UMIST   */
	newreact("NH2,CRP=>NH,H",crnurate,81 * 2,0.,0.); /* UMIST   */
	newreact("NH2,CRP=>NH2+,e-",crnurate,610 * 2,0.,0.); /* UMIST   */
	newreact("NH3,CRP=>NH3+,e-",crnurate,543 * 2,0.,0.); /* UMIST   */
	newreact("NH3,CRP=>NH,H2",crnurate,541 * 2,0.,0.); /* UMIST   */
	newreact("NH3,CRP=>NH2,H",crnurate,1325 * 2,0.,0.); /* UMIST   */
	newreact("CN,CRP=>N,C",crnurate,10580 * 2,0.,0.); /* UMIST   */
	newreact("HCN,CRP=>CN,H",crnurate,3114 * 2,0.,0.); /* UMIST   */
	newreact("N2,CRP=>N,N",crnurate,50 * 2,0.,0.); /* UMIST   */
	newreact("NO,CRP=>O,N",crnurate,427 * 2,0.,0.); /* UMIST   */
	newreact("NO,CRP=>NO+,e-",crnurate,430 * 2,0.,0.); /* UMIST   */
	newreact("HNO,CRP=>NO,H",crnurate,1000 * 2,0.,0.); /* UMIST   */
	newreact("HS,CRP=>S,H",crnurate,500 * 2,0.,0.); /* UMIST   */
	newreact("OCN,CRP=>CN,O",crnurate,1500 * 2,0.,0.); /* UMIST   */
	newreact("CS,CRP=>S,C",crnurate,500 * 2,0.,0.); /* UMIST   */
	newreact("CS,CRP=>CS+,e-",crnurate,500 * 2,0.,0.); /* UMIST   */
	newreact("NO2,CRP=>NO,O",crnurate,1500 * 2,0.,0.); /* UMIST   */
	newreact("NS,CRP=>S,N",crnurate,500 * 2,0.,0.); /* UMIST   */
	newreact("SO,CRP=>S,O",crnurate,500 * 2,0.,0.); /* UMIST   */
	newreact("SO,CRP=>SO+,e-",crnurate,500 * 2,0.,0.); /* UMIST   */
	newreact("OCS,CRP=>OCS+,e-",crnurate,1444 * 2,0.,0.); /* UMIST   */
	newreact("PHOTON,CH=>CH+,e-",th85rate,7.6e-10,2.8,0.); /*  UMIST	 */  	/* >>chng 06 jun 01, comment out "??4.2e-10*TH85_field;" - in more than one time*/
	newreact("PHOTON,CH=>C,H",th85rate,8.6e-10,1.2,0.); /*  UMIST	 */  	/* >>chng 06 jun 01, comment out "??1.3e-9*TH85_field;" - in more than one time*/
	newreact("PHOTON,CH+=>C+,H",th85rate,2.5e-10,2.5,0.); /*  UMIST	 */
	newreact("PHOTON,CH2=>CH2+,e-",th85rate,1.0e-9,2.3,0.); /*  UMIST	 */
	newreact("PHOTON,CH2=>CH,H",th85rate,7.2e-10,1.7,0.); /*  UMIST	 */
	newreact("PHOTON,CH2+=>CH+,H",th85rate,1.7e-9,1.7,0.); /*  UMIST	 */
	newreact("PHOTON,CH3+=>CH2+,H",th85rate,1.0e-9,1.7,0.); /*  UMIST	 */
	newreact("PHOTON,CH3+=>CH+,H2",th85rate,1.0e-9,1.7,0.); /*  UMIST	 */
	newreact("PHOTON,H2O=>H2O+,e-",th85rate,3.3e-11,3.9,0.); /*  UMIST	 */
	newreact("PHOTON,H2O=>OH,H",th85rate,5.9e-10,1.7,0.); /*  UMIST	 */
	newreact("PHOTON,O2=>O2+,e-",th85rate,5.6e-11,3.7,0.); /*  UMIST	 */
	newreact("PHOTON,O2=>O,O",th85rate,6.9e-10,1.8,0.); /*  UMIST	 */
	newreact("PHOTON,OH=>OH+,e-",th85rate,1.6e-12,3.1,0.); /*  UMIST	 */
	newreact("PHOTON,OH=>O,H",th85rate,3.5e-10,1.7,0.); /*  UMIST	 */
	newreact("PHOTON,OH+=>O,H+",th85rate,1.e-12,1.8,0.); /*  UMIST	 */
	newreact("PHOTON,SiH=>Si,H",th85rate,2.8e-9,1.1,0.); /*  UMIST	 */
	newreact("PHOTON,SiO=>Si,O",th85rate,1.0e-10,2.3,0.); /*  UMIST	 */
	newreact("HNC,PHOTON=>CN,H",th85rate,1.5e-9,0.,0.); /*  */
	newreact("HNC,CRPHOT=>CN,H",crnurate,2*2986,0.,0.); /*  */
	newreact("HCl,PHOTON=>Cl,H",th85rate,0.00000000011,0.,0.); /* UMIST   */
	newreact("CCl,PHOTON=>Cl,C",th85rate,0.0000000001,0.,0.); /* UMIST   */
	newreact("ClO,PHOTON=>Cl,O",th85rate,0.0000000001,0.,0.);/*  UMIST   */
	newreact("HCl,CRPHOT=>Cl,H",crnurate,2*610,0.,0.); /* UMIST   */
	newreact("CCl,CRPHOT=>Cl,C",crnurate,2*500,0.,0.); /* UMIST   */
	newreact("ClO,CRPHOT=>Cl,O",crnurate,2*500,0.,0.); /* UMIST   */
	newreact("cr,CH=>C,H",crnurate,2. * 756,0.,0.); /* UMIST	 */
	newreact("cr,CH+=>C+,H",crnurate,2. * 183,0.,0.); /* UMIST	 */
	newreact("cr,H2O=>OH,H",crnurate,2. * 979,0.,0.); /* UMIST	 */
	newreact("cr,O2=>O2+,e-",crnurate,2. *88,0.,0.); /* UMIST	 */
	newreact("cr,O2=>O,O",crnurate,2. *730,0.,0.); /* UMIST	 */
	newreact("cr,OH=>O,H",crnurate,2. *522,0.,0.); /* UMIST	 */
	newreact("cr,SiH=>Si,H",crnurate,2. *500,0.,0.); /* UMIST	 */
	newreact("cr,SiO=>Si,O",crnurate,2. *500,0.,0.); /* UMIST	 */
	newreact("PHOTON,CH3=>CH3+,e-",th85rate,1.0e-10,2.1,0.); /*  UMIST     */
	newreact("PHOTON,CH3=>CH2,H",th85rate,2.5e-10,1.9,0.); /*  UMIST     */
	newreact("PHOTON,CH3=>CH,H2",th85rate,2.5e-10,1.9,0.); /*  UMIST     */
	newreact("CRPHOT,CH3=>CH3+,e-",crnurate,2.*500,0.,0.); /* UMIST     */
	newreact("CRPHOT,CH3=>CH2,H",crnurate,2.*500,0.,0.); /* UMIST     */
	newreact("CRPHOT,CH3=>CH,H2",crnurate,2.*500,0.,0.); /* UMIST     */
	newreact("PHOTON,CH4=>CH3,H",th85rate,2.2e-10,2.2,0.); /*  UMIST     */
	newreact("PHOTON,CH4=>CH2,H2",th85rate,9.8e-10,2.2,0.); /*  UMIST     */
	newreact("PHOTON,CH4=>CH,H2,H",th85rate,2.2e-10,2.2,0.); /*  UMIST     */
	newreact("CRPHOT,CH4=>CH2,H2",crnurate,2.*2272,0.,0.); /* UMIST     */
	newreact("C2,PHOTON=>C,C",th85rate,4.7e-11,0.,0.); /* UMIST	 */ 	/*>>chng 05 dec 17 The rate for this reaction comes from Eric Herbst's website */
	newreact("C2,PHOTON=>C2+,e-",th85rate,1e-10,0.,0.); /* UMIST	 */ 	/*>>chng 05 dec 17 The rate for this reaction comes from Eric Herbst's website */
	newreact("C2,cr=>C,C",crnurate,2. *500,0.,0.); /* UMIST	 */
	newreact("C2+,PHOTON=>C+,C",th85rate,1.0e-11,0.,0.); /* UMIST	 */
	newreact("C2H,CRPHOT=>C2,H",crnurate,2*5000,0.,0.); /*  */
	newreact("C2H,PHOTON=>C2,H",th85rate,0.00000000051,0.,0.); /*  */
	newreact("C2H,PHOTON=>C2H+,e-",th85rate,0.00000000001,0.,0.); /*  */
	newreact("C2H,CRPHOT=>C2H+,e-",crnurate,2*5000,0.,0.); /*  */
	newreact("C3,PHOTON=>C2,C",th85rate,0.0000000038,0.,0.); /*  */
	newreact("C3,CRPHOT=>C2,C",crnurate,2*1119,0.,0.); /*  */
	newreact("C2H2,CRPHOT=>C2H2+,e-",crnurate,2*1309,0.,0.); /*  */
	newreact("C2H2,CRPHOT=>C2H,H",crnurate,2*5155,0.,0.); /*  */
	newreact("C3H,CRPHOT=>C3,H",crnurate,2*5000,0.,0.); /*  */
	newreact("C2H2,PHOTON=>C2H2+,e-",th85rate,0.00000000046,0.,0.); /*  */
	newreact("C2H2,PHOTON=>C2H,H",th85rate,0.0000000073,0.,0.); /*  */
	newreact("C3H,PHOTON=>C3,H",th85rate,0.000000001,0.,0.); /*  */
	newreact("H,CH=>C,H,H",hmrate,6.0e-9,0,40200); /* UMIST	 */
	newreact("H,OH=>O,H,H",hmrate,6.0e-9,0,50900); /* UMIST	 */
	newreact("H,H2O=>OH,H,H",hmrate,5.8e-9,0,52900); /* UMIST	 */
	if(!co.lgFederman)
	{
		newreact("O,CH=>CO,H",hmrate,6.6e-11,0.,0.);/*  UMIST */
		newreact("e-,CH2+=>CH,H",hmrate,0.00000016,-0.6,0);/*  UMIST */
		newreact("e-,CH2+=>C,H2",hmrate,0.0000000768,-0.6,0); /* UMIST	 */
		newreact("e-,CH3+=>C,H2,H",hmrate,0.000000105,-0.5,0); /* UMIST; Federman	 */
		newreact("e-,CH3+=>CH2,H",hmrate,0.00000014,-0.5,0); /* UMIST;Federman	 */
		newreact("e-,CH3+=>CH,H2",hmrate,0.000000049,-0.5,0); /* UMIST; Federman	 */
		newreact("e-,CO+=>C,O",hmrate,0.0000002,-0.48,0); /* UMIST; Federman	 */
		newreact("H,CH=>C,H2",hmrate,2.7e-11,0.38,2200); /* UMIST; Federman	 */ 	/* >>chng 06 Feb 27 -- NPA  TH85 has a temperature barrier for this reaction, and the Meudon PDR code datafile does as well.  Use the temperature barrier of 2200 K */
		newreact("H2,CO+=>HCO+,H",hmrate,1.8e-9,0.,0.); /* UMIST;Federman;	 */
		newreact("H2,CH2+=>CH3+,H",hmrate,1.6e-9,0.,0.); /*  UMIST;Federman;*/
		newreact("H2,C+=>CH2+,PHOTON",hmrate,0.0000000000000004,-0.2,0); /* UMIST;Federman	 */
		newreact("N,C2=>CN,C",hmrate,5.e-11,0.5,0); /* UMIST; Federman	 */
		newreact("C+,CH=>C2+,H",hmrate,27e-9,-0.50,0); /* UMIST;Federman */
		newreact("C+,CH=>CH+,C",hmrate,3.8e-10,0.,0.); /* UMIST;Federman	 */
		newreact("O,CH3+=>HCO+,H2",hmrate,4.0e-10,0.,0.); /* UMIST;Federman	 */
		newreact("CH,N=>CN,H",hmrate,1.66e-10,-0.09,0); /* UMIST; Federman  */
		newreact("C+,NH=>CN+,H",hmrate,7.8e-10,0.,0.); /* UMIST; Federman   */
	}
	else
	{
		newreact("O,CH=>CO,H",hmrate,9.5e-11,0.5,0); /* UMIST; Federman	 */
		newreact("e-,CH2+=>CH,H",hmrate,1.25e-7,-0.5,0); /* UMIST; Federman	 */
		newreact("e-,CH2+=>C,H2",hmrate,1.25e-7,-0.5,0); /* UMIST; Federman	 */
		newreact("e-,CH3+=>C,H2,H",hmrate,3.5e-8,-0.5,0); /* UMIST; Federman	 */
		newreact("e-,CH3+=>CH2,H",hmrate,7.0e-8,-0.5,0); /* UMIST;Federman	 */
		newreact("e-,CH3+=>CH,H2",hmrate,2.45e-7,-0.5,0); /* UMIST; Federman	 */
		newreact("e-,CO+=>C,O",hmrate,1.0e-7,-0.46,0); /* UMIST; Federman	 */
		newreact("H,CH=>C,H2",hmrate,5.0e-11,0.5,2200); /* UMIST; Federman	 */ 	/* >>chng 06 Feb 27 -- NPA  TH85 has a temperature barrier for this reaction, and the Meudon PDR code datafile does as well.  Use the temperature barrier of 2200 K */
		newreact("H2,CO+=>HCO+,H",hmrate,1.4e-9,0.,0.); /* UMIST;Federman;	 */
		newreact("H2,CH2+=>CH3+,H",hmrate,1.2e-9,0.,0.); /*  UMIST;Federman;*/
		newreact("H2,C+=>CH2+,PHOTON",hmrate,1.0e-16,0.,0.); /* UMIST;Federman	 */
		newreact("N,C2=>CN,C",hmrate,1.7e-11,0.5,0); /* UMIST; Federman	 */
		newreact("C+,CH=>C2+,H",hmrate,5.4e-10,0.,0.); /* UMIST;Federman */
		newreact("C+,CH=>CH+,C",hmrate,1.5e-9,0.,0.); /* UMIST;Federman	 */
		newreact("O,CH3+=>HCO+,H2",hmrate,4.4e-10,0.,0.); /* UMIST;Federman	 */
		newreact("CH,N=>CN,H",hmrate,2e-11,0.5,0); /* UMIST; Federman  */
		newreact("C+,NH=>CN+,H",hmrate,5.0e-10,0.,0.); /* UMIST; Federman   */
	}

	newreact("e-,HCO+=>CO,H",hmrate,2.4e-7,-0.69,0); /* UMIST; Federman	 */
	newreact("H,CH+=>C+,H2",hmrate,7.5e-10,0,0); /* UMIST	 */
	newreact("H,CH2=>CH,H2",hmrate,6.64e-11,0,370); /* UMIST	 */ 	/* >>chng 06 Feb 27 -- NPA  TH85 has a temperature barrier for this reaction, and the Meudon PDR code datafile does as well.  Use the temperature barrier of 370 K*/	
	newreact("H,CH3+=>CH2+,H2",hmrate,7.0e-10,0,10560); /* UMIST	 */
	newreact("H,OH=>O,H2",hmrate,0.0000000000000699,2.8,1950); /* UMIST	 */
	newreact("H,H2O=>OH,H2",hmrate,1.59e-11,1.2,9610); /* UMIST	 */
	newreact("H,O2=>OH,O",hmrate,2.61e-10,0,8156); /* UMIST	 */
	newreact("H,O2=>O,O,H",hmrate,6.0e-9,0,52300); /* UMIST	 */
	newreact("H,C=>CH,PHOTON",hmrate,1.e-17,0,0); /* UMIST	 */
	newreact("H,C+=>CH+,PHOTON",hmrate,1.7e-17,0,0); /* UMIST	 */
	newreact("H,OH=>H2O,PHOTON",hmrate,5.26E-18,-5.22,90); /* UMIST	 */
	newreact("H,O=>OH,PHOTON",hmrate,9.9E-19,-0.38,0); /* UMIST	 */
	newreact("H-,CH=>CH2,e-",hmrate,1.0e-10,0,0); /* UMIST	 */
	newreact("H-,C=>CH,e-",hmrate,1.0e-9,0,0); /* UMIST	 */
	newreact("H-,OH=>H2O,e-",hmrate,1.0e-10,0,0); /* UMIST	 */
	newreact("H-,O=>OH,e-",hmrate,1.0e-9,0,0); /* UMIST	 */
	newreact("H-,H3O+=>H2O,H2",hmrate,0.00000023,-0.5,0); /* UMIST	 */
	newreact("H-,H3O+=>OH,H2,H",hmrate,0.00000023,-0.5,0); /* UMIST	 */
	newreact("H+,CH=>CH+,H",hmrate,1.9e-9,0,0); /* UMIST	 */
	newreact("H+,CH2=>CH2+,H",hmrate,1.4e-9,0,0); /* UMIST	 */
	newreact("H+,H2O=>H2O+,H",hmrate,6.9e-9,0,0); /* UMIST	 */
	newreact("H+,O2=>O2+,H",hmrate,2.0e-9,0,0); /* UMIST	 */
	newreact("H+,OH=>OH+,H",hmrate,2.1e-9,0,0); /* UMIST	 */
	newreact("H+,SiO=>SiO+,H",hmrate,3.3e-9,0,0); /* UMIST	 */
	newreact("H+,CH2=>CH+,H2",hmrate,1.4e-9,0,0); /* UMIST	 */
	newreact("H+,SiH=>Si+,H2",hmrate,1.7e-9,0,0); /* UMIST	 */
	newreact("H2,C=>CH,H",hmrate,6.64e-10,0,11700); /* UMIST	 */
	newreact("H2,C+=>CH+,H",hmrate,1.0e-10,0,4640); /* UMIST	 */
	newreact("H2,CH=>CH2,H",hmrate,5.46e-10,0,1943); /* UMIST	 */
	newreact("H2,CH+=>CH2+,H",hmrate,1.2e-9,0,0); /* UMIST	 */
	newreact("H2,OH=>H2O,H",hmrate,2.05e-12,1.52,1736); /* UMIST	 */
	newreact("H2,OH+=>H2O+,H",hmrate,1.01e-9,0,0); /* UMIST	 */
	newreact("H2,H2O+=>H3O+,H",hmrate,6.4e-10,0,0); /* UMIST	 */
	newreact("H2,O=>OH,H",hmrate,3.14e-13,2.7,3150); /* UMIST	 */
	newreact("H2,O+=>OH+,H",hmrate,1.7e-9,0,0); /* UMIST	 */
	newreact("H2,SiO+=>SiOH+,H",hmrate,3.2e-10,0,0); /* UMIST	 */
	newreact("H2,CH=>C,H2,H",hmrate,6.0e-9,0,40200); /* UMIST	 */
	newreact("H2,OH=>O,H2,H",hmrate,6.0e-9,0,50900); /* UMIST	 */
	newreact("H2,H2O=>OH,H2,H",hmrate,5.8e-9,0,52900); /* UMIST	 */
	newreact("H2,O2=>O,O,H2",hmrate,6.0e-9,0,52300); /* UMIST	 */
	newreact("H2,O2=>OH,OH",hmrate,3.16e-10,0,21890); /* UMIST	 */
	newreact("H2,C=>CH2,PHOTON",hmrate,1.e-17,0,0); /* UMIST	 */
	newreact("H2,Si+=>SiH2+,PHOTON",hmrate,0.000000000000000003,0,0); /* UMIST	 */
	newreact("H2+,C=>CH+,H",hmrate,2.4e-9,0,0); /* UMIST	 */
	newreact("H2+,CH=>CH2+,H",hmrate,7.1e-10,0,0); /* UMIST	 */
	newreact("H2+,CH2=>CH3+,H",hmrate,1.0e-9,0,0); /* UMIST	 */
	newreact("H2+,OH=>H2O+,H",hmrate,7.6e-10,0,0); /* UMIST	 */
	newreact("H2+,H2O=>H3O+,H",hmrate,3.4e-9,0,0); /* UMIST	 */
	newreact("H2+,O=>OH+,H",hmrate,1.5e-9,0,0); /* UMIST	 */
	newreact("H2+,CH=>CH+,H2",hmrate,7.1e-10,0,0); /* UMIST	 */
	newreact("H2+,CH2=>CH2+,H2",hmrate,1.0e-9,0,0); /* UMIST	 */
	newreact("H2,S=>HS,H",hmrate,1.76e-13,2.88,6126); /* UMIST	 *//* >>chng 05 aug 02, NA added this */
	newreact("H2+,H2O=>H2O+,H2",hmrate,3.9e-9,0,0); /* UMIST	 */
	newreact("H2+,O2=>O2+,H2",hmrate,8.e-10,0,0); /* UMIST	 */
	newreact("H2+,OH=>OH+,H2",hmrate,7.6e-10,0,0); /* UMIST	 */
	newreact("H3+,C=>CH+,H2",hmrate,2.0e-9,0,0); /* UMIST	 */
	newreact("H3+,CH=>CH2+,H2",hmrate,1.2e-9,0,0); /* UMIST	 */
	newreact("H3+,CH2=>CH3+,H2",hmrate,1.7e-9,0,0); /* UMIST	 */
	newreact("H3+,OH=>H2O+,H2",hmrate,1.3e-9,0,0); /* UMIST	 */
	newreact("H3+,H2O=>H3O+,H2",hmrate,5.9e-9,0,0); /* UMIST	 */
	newreact("H3+,O=>OH+,H2",hmrate,8.e-10,0,0); /* UMIST	 */
	newreact("H3+,SiH=>SiH2+,H2",hmrate,2.0e-9,0,0); /* UMIST	 */
	newreact("H3+,SiO=>SiOH+,H2",hmrate,2.0e-9,0,0); /* UMIST	 */
	newreact("He+,CH=>CH+,He",hmrate,5.0e-10,0,0); /* UMIST	 */
	newreact("He+,H2O=>H2O+,He",hmrate,6.05e-11,0,0); /* UMIST	 */
	newreact("He+,O2=>O2+,He",hmrate,3.3e-11,0,0); /* UMIST	 */
	newreact("He+,CH=>C+,H,He",hmrate,1.1e-9,0,0); /* UMIST	 */
	newreact("He+,CH2=>CH+,H,He",hmrate,7.5e-10,0,0); /* UMIST	 */
	newreact("He+,OH=>O+,H,He",hmrate,1.1e-9,0,0); /* UMIST	 */
	newreact("He+,H2O=>OH+,H,He",hmrate,2.86e-10,0,0); /* UMIST	 */
	newreact("He+,SiH=>Si+,H,He",hmrate,1.8e-9,0,0); /* UMIST	 */
	newreact("He+,H2O=>OH,H+,He",hmrate,2.04e-10,0,0); /* UMIST	 */
	newreact("He+,CH2=>C+,H2,He",hmrate,7.5e-10,0,0); /* UMIST	 */
	newreact("He+,O2=>O+,O,He",hmrate,1.0e-9,0,0); /* UMIST	 */
	newreact("He+,SiO=>Si+,O,He",hmrate,8.6e-10,0,0); /* UMIST	 */
	newreact("He+,SiO=>Si,O+,He",hmrate,8.6e-10,0,0); /* UMIST	 */
	newreact("He+,Si=>Si+,He",hmrate,0.0000000033,0,0); /*  */
	newreact("H,CH3=>CH2,H2",hmrate,1.0e-10,0,7600); /* UMIST     */
	newreact("H,CH4+=>CH3+,H2",hmrate,1.0e-11,0,0); /* UMIST     */
	newreact("H,CH5+=>CH4+,H2",hmrate,2.0e-11,0,0); /* UMIST     */
	newreact("H2,CH2=>CH3,H",hmrate,5.18e-11,0.17,6400); /* UMIST     */
	newreact("H2,CH=>CH3,PHOTON",hmrate,5.09E-18,-0.71,11.6); /* UMIST     */
	newreact("H2,CH3=>CH4,H",hmrate,0.0000000000000686,2.74,4740); /* UMIST     */
	newreact("H2,CH4+=>CH5+,H",hmrate,3.3e-11,0,0); /* UMIST     */
	newreact("H2,CH3+=>CH5+,PHOTON",hmrate,0.000000000000013,-1,0); /* UMIST     */
	newreact("H2+,CH4=>CH3+,H2,H",hmrate,2.3e-9,0,0); /* UMIST     */
	newreact("H2+,CH4=>CH4+,H2",hmrate,1.4e-9,0,0); /* UMIST     */
	newreact("H2+,CH4=>CH5+,H",hmrate,1.14e-10,0,0); /* UMIST     */
	newreact("H3+,CH3=>CH4+,H2",hmrate,2.1e-9,0,0); /* UMIST     */
	newreact("H3+,CH4=>CH5+,H2",hmrate,2.4e-9,0,0); /* UMIST     */
	newreact("He+,CH3=>CH+,He,H2",hmrate,1.8e-9,0,0); /* UMIST     */
	newreact("He+,CH4=>CH+,He,H2,H",hmrate,2.4e-10,0,0); /* UMIST     */
	newreact("He+,CH4=>CH2+,He,H2",hmrate,9.5e-10,0,0); /* UMIST     */
	newreact("He+,CH4=>CH3,He,H+",hmrate,4.8e-10,0,0); /* UMIST     */
	newreact("He+,CH4=>CH3+,He,H",hmrate,8.5e-11,0,0); /* UMIST     */
	newreact("He+,CH4=>CH4+,He",hmrate,5.1e-11,0,0); /* UMIST     */
	newreact("H-,CH2=>CH3,e-",hmrate,1.0e-9,0,0); /* UMIST     */
	newreact("H-,CH3=>CH4,e-",hmrate,1.0e-9,0,0); /* UMIST     */
	newreact("H+,CH3=>CH3+,H",hmrate,3.4e-9,0,0); /* UMIST     */
	newreact("H+,CH4=>CH3+,H2",hmrate,2.3e-9,0,0); /* UMIST     */
	newreact("H+,CH4=>CH4+,H",hmrate,1.5e-9,0,0); /* UMIST     */
	newreact("He+,C2=>C+,C,He",hmrate,1.6e-9,0,0); /* UMIST	 */
	if(hmi.lgLeiden_Keep_ipMH2s) 
	{
		newreact("H2*,CH=>C,H2,H",hmrate,6.0e-9,0.,0.); /*  */
		newreact("H2*,OH=>O,H2,H",hmrate,6.0e-9,0.,0.); /*  */
		newreact("H2*,H2O=>OH,H2,H",hmrate,5.8e-9,0.,0.); /*  */
		newreact("H2*,O2=>O,O,H2",hmrate,6.0e-9,0.,0.); /*  */
		newreact("H2*,CH2=>CH3,H",hmrate,5.18e-11,0.17,0); /*  */
		newreact("H2*,CH=>CH3,PHOTON",hmrate,5.09E-18,-0.71,0); /*  */
		newreact("H2*,CH3=>CH4,H",hmrate,0.0000000000000686,2.74,0); /*  */
		newreact("H2*,CH=>CH2,H",hmrate,5.46e-10,0, 0. ); /* Tielens & Hollenbach 1985, Ap. J. 291, 722 Table 9 */
		newreact("H2*,O=>OH,H",hmrate,3.14e-13,2.7, 0. ); /* Tielens & Hollenbach 1985, Ap. J. 291, 722 Table 9 */
		newreact("H2*,OH=>H2O,H",hmrate,2.05e-12,1.52, 0. ); /* Tielens & Hollenbach 1985, Ap. J. 291, 722 Table 9 */
		newreact("H2*,O2=>OH,OH",hmrate,3.16e-10,0, 0. ); /* Tielens & Hollenbach 1985, Ap. J. 291, 722 Table 9 */
		newreact("H2*,C=>CH,H",hmrate,6.64e-10,0, 0. ); /* Tielens & Hollenbach 1985, Ap. J. 291, 722 Table 9 */
		newreact("H2*,C+=>CH+,H",hmrate,1.0e-10,0, 0. ); /* Tielens & Hollenbach 1985, Ap. J. 291, 722 Table 9 */
		newreact("H2*,O+=>OH+,H",hmrate,1.7e-9,0.,0.); /*  */
	}
	newreact("H2,N=>NH,H",hmrate,1.69e-9,0,18095); /* UMIST   */
	newreact("H2,NH=>NH2,H",hmrate,5.96e-11,0,7782); /* UMIST   */
	newreact("H2,NH2=>NH3,H",hmrate,2.05e-15,3.89,1400); /* UMIST   */
	newreact("H2,CN=>HCN,H",hmrate,0.000000000000404,2.87,820); /* UMIST   */
	newreact("H+,HNO=>NO+,H2",hmrate,4.e-9,0,0); /* UMIST   */
	newreact("H+,HS=>S+,H2",hmrate,1.6e-9,0,0); /* UMIST   */
	newreact("H,HS+=>S+,H2",hmrate,1.1e-10,0,0); /* UMIST   */
	newreact("H2+,N=>NH+,H",hmrate,1.9e-9,0,0); /* UMIST   */
	newreact("H2,N+=>NH+,H",hmrate,1.0e-9,0,85); /* UMIST   */
	newreact("H2,NH+=>N,H3+",hmrate,2.25e-10,0,0); /* UMIST   */
	newreact("H2+,NH=>NH2+,H",hmrate,7.6e-10,0,0); /* UMIST   */
	newreact("H2,NH+=>NH2+,H",hmrate,1.28e-9,0,0); /* UMIST   */
	newreact("H2,NH2+=>NH3+,H",hmrate,2.7e-10,0,0); /* UMIST   */
	newreact("H2,NH3+=>NH4+,H",hmrate,0.0000000000002,0,0); /* UMIST   */
	newreact("H2+,CN=>HCN+,H",hmrate,1.2e-9,0,0); /* UMIST   */
	newreact("H2,CN+=>HCN+,H",hmrate,1.0e-9,0,0); /* UMIST   */
	newreact("H2+,NO=>HNO+,H",hmrate,1.1e-9,0,0); /* UMIST   */
	newreact("H2,S+=>HS+,H",hmrate,1.1e-10,0,9860); /* UMIST   */
	newreact("H2,CS+=>HCS+,H",hmrate,4.5e-10,0,0); /* UMIST   */
	newreact("H2,NO2+=>NO+,H2O",hmrate,1.5e-10,0,0); /* UMIST   */
	newreact("H3+,NH=>NH2+,H2",hmrate,1.3e-9,0,0); /* UMIST   */
	newreact("H3+,NH2=>NH3+,H2",hmrate,1.8e-9,0,0); /* UMIST   */
	newreact("H3+,NH3=>NH4+,H2",hmrate,2.7e-9,0,0); /* UMIST   */
	newreact("H3+,CN=>HCN+,H2",hmrate,2.0e-9,0,0); /* UMIST   */
	newreact("H3+,NO=>HNO+,H2",hmrate,1.1e-9,0,0); /* UMIST   */
	newreact("H3+,S=>HS+,H2",hmrate,2.6e-9,0,0); /* UMIST   */
	newreact("H3+,CS=>HCS+,H2",hmrate,2.9e-9,0,0); /* UMIST   */
	newreact("H3+,NO2=>NO+,OH,H2",hmrate,7.0e-10,0,0); /* UMIST   */
	newreact("He+,NH=>N+,He,H",hmrate,1.1e-9,0,0); /* UMIST   */
	newreact("He+,NH2=>N+,He,H2",hmrate,8.e-10,0,0); /* UMIST   */
	newreact("He+,NH2=>NH+,He,H",hmrate,8.e-10,0,0); /* UMIST   */
	newreact("He+,NH3=>NH+,He,H2",hmrate,1.76e-10,0,0); /* UMIST   */
	newreact("He+,NH3=>NH2+,He,H",hmrate,1.76e-9,0,0); /* UMIST   */
	newreact("He+,CN=>N,C+,He",hmrate,8.8e-10,0,0); /* UMIST   */
	newreact("He+,CN=>N+,C,He",hmrate,8.8e-10,0,0); /* UMIST   */
	newreact("He+,HCN=>N,CH+,He",hmrate,6.51e-10,0,0); /* UMIST   */
	newreact("He+,HCN=>N+,CH,He",hmrate,2.17e-10,0,0); /* UMIST   */
	newreact("He+,HCN=>N,C+,He,H",hmrate,7.75e-10,0,0); /* UMIST   */
	newreact("He+,HCN=>CN+,He,H",hmrate,1.46e-9,0,0); /* UMIST   */
	newreact("He+,N2=>N+,N,He",hmrate,9.6e-10,0,0); /* UMIST   */
	newreact("He+,NO=>O+,N,He",hmrate,2.0e-10,0,0); /* UMIST   */
	newreact("He+,NO=>O,N+,He",hmrate,1.4e-9,0,0); /* UMIST   */
	newreact("He+,HNO=>NO+,He,H",hmrate,1.0e-9,0,0); /* UMIST   */
	newreact("He+,HNO=>NO,He,H+",hmrate,1.0e-9,0,0); /* UMIST   */
	newreact("He+,HS=>S+,He,H",hmrate,1.7e-9,0,0); /* UMIST   */
	newreact("He+,OCN=>CN,O+,He",hmrate,3.0e-9,0,0); /* UMIST   */
	newreact("He+,OCN=>CN+,O,He",hmrate,3.0e-9,0,0); /* UMIST   */
	newreact("He+,SiN=>Si+,N,He",hmrate,2.0e-9,0,0); /* UMIST   */
	newreact("He+,N2O=>N2,O+,He",hmrate,2.76e-10,0,0); /* UMIST   */
	newreact("He+,N2O=>N2+,O,He",hmrate,1.24e-9,0,0); /* UMIST   */
	newreact("He+,N2O=>NO,N+,He",hmrate,3.e-10,0,0); /* UMIST   */
	newreact("He+,N2O=>NO+,N,He",hmrate,4.83e-10,0,0); /* UMIST   */
	newreact("He+,CS=>S+,C,He",hmrate,1.3e-9,0,0); /* UMIST   */
	newreact("He+,CS=>S,C+,He",hmrate,1.3e-9,0,0); /* UMIST   */
	newreact("He+,NS=>S,N+,He",hmrate,1.2e-9,0,0); /* UMIST   */
	newreact("He+,NS=>S+,N,He",hmrate,1.2e-9,0,0); /* UMIST   */
	newreact("He+,SO=>S,O+,He",hmrate,8.3e-10,0,0); /* UMIST   */
	newreact("He+,SO=>S+,O,He",hmrate,8.3e-10,0,0); /* UMIST   */
	newreact("He+,OCS=>S,CO+,He",hmrate,7.6e-10,0,0); /* UMIST   */
	newreact("He+,OCS=>CS+,O,He",hmrate,7.6e-10,0,0); /* UMIST   */
	newreact("He+,OCS=>CS,O+,He",hmrate,7.6e-10,0,0); /* UMIST   */
	newreact("He+,S2=>S+,S,He",hmrate,2.0e-9,0,0); /* UMIST   */
	newreact("H+,NH=>NH+,H",hmrate,2.1e-9,0,0); /* UMIST   */
	newreact("H+,NH2=>NH2+,H",hmrate,2.9e-9,0,0); /* UMIST   */
	newreact("H+,NH3=>NH3+,H",hmrate,1.1e-9,0,0); /* UMIST   */
	newreact("H,CN+=>CN,H+",hmrate,1.9e-10,0,0); /* UMIST   */
	newreact("H+,HCN=>HCN+,H",hmrate,0.0000000105,-0.13,0); /* UMIST   */
	newreact("H,HCN+=>HCN,H+",hmrate,3.7e-11,0,0); /* UMIST   */
	newreact("H,N2+=>N2,H+",hmrate,1.2e-10,0,0); /* UMIST   */
	newreact("H+,NO=>NO+,H",hmrate,2.9e-9,0,0); /* UMIST   */
	newreact("H+,HS=>HS+,H",hmrate,1.6e-9,0,0); /* UMIST   */
	newreact("H+,SiN=>SiN+,H",hmrate,3.0e-9,0,0); /* UMIST   */
	newreact("H+,CS=>CS+,H",hmrate,4.9e-9,0,0); /* UMIST   */
	newreact("H+,NS=>NS+,H",hmrate,4.7e-9,0,0); /* UMIST   */
	newreact("H+,SO=>SO+,H",hmrate,3.2e-9,0,0); /* UMIST   */
	newreact("H+,OCS=>OCS+,H",hmrate,2.1e-9,0,0); /* UMIST   */
	newreact("H+,S2=>S2+,H",hmrate,3.0e-9,0,0); /* UMIST   */
	newreact("H2+,NH=>NH+,H2",hmrate,7.6e-10,0,0); /* UMIST   */
	newreact("H2+,NH2=>NH2+,H2",hmrate,2.1e-9,0,0); /* UMIST   */
	newreact("H2+,NH3=>NH3+,H2",hmrate,5.7e-9,0,0); /* UMIST   */
	newreact("H2+,CN=>CN+,H2",hmrate,1.2e-9,0,0); /* UMIST   */
	newreact("H2+,HCN=>HCN+,H2",hmrate,2.7e-9,0,0); /* UMIST   */
	newreact("H2+,NO=>NO+,H2",hmrate,1.1e-9,0,0); /* UMIST   */
	newreact("He+,NH3=>NH3+,He",hmrate,2.64e-10,0,0); /* UMIST   */
	newreact("He+,N2=>N2+,He",hmrate,6.4e-10,0,0); /* UMIST   */
	newreact("H-,N=>NH,e-",hmrate,1.0e-9,0,0); /* UMIST   */
	newreact("H-,NH=>NH2,e-",hmrate,1.0e-10,0,0); /* UMIST   */
	newreact("H-,NH2=>NH3,e-",hmrate,1.0e-9,0,0); /* UMIST   */
	newreact("H-,CN=>HCN,e-",hmrate,1.0e-10,0,0); /* UMIST   */
	newreact("H-,NH4+=>NH3,H2",hmrate,0.00000023,-0.5,0); /* UMIST   */
	newreact("H-,N+=>N,H",hmrate,0.00000023,-0.5,0); /* UMIST   */
	newreact("H2,Cl+=>HCl+,H",hmrate,0.000000001,0,0); /* UMIST   */
	newreact("H2,HCl+=>H2Cl+,H",hmrate,0.0000000013,0,0); /* UMIST   */
	newreact("H3+,Cl=>HCl+,H2",hmrate,0.000000001,0,0); /* UMIST   */
	newreact("H3+,HCl=>H2Cl+,H2",hmrate,0.0000000038,0,0); /* UMIST   */
	newreact("He+,HCl=>Cl+,He,H",hmrate,0.0000000033,0,0); /* UMIST   */
	newreact("He+,CCl=>Cl,C+,He",hmrate,0.0000000033,0,0); /* UMIST   */
	newreact("He+,ClO=>Cl+,O,He",hmrate,0.000000001,0,0); /* UMIST   */
	newreact("H+,HCl=>HCl+,H",hmrate,0.0000000033,1,0); /* UMIST   */
	newreact("H3+,HNC=>HCNH+,H2",hmrate,0.0000000081,0,0); /*  */
	newreact("He+,HNC=>NH+,C,He",hmrate,0.0000000005,0,0); /*  */
	newreact("He+,HNC=>N,C+,He,H",hmrate,0.0000000005,0,0); /*  */
	newreact("He+,HNC=>CN+,He,H",hmrate,0.0000000005,0,0); /*  */
	newreact("H+,HNC=>HCN,H+",hmrate,0.000000001,0,0); /*  */
	newreact("H,HNC=>HCN,H",h_hnc_hcn_h,0.000000000000136,4.48,0);/*   */
	newreact("H2,HCN+=>HCNH+,H",hmrate,0.0000000009,0,0); /*  */
	newreact("H3+,HCN=>HCNH+,H2",hmrate,0.0000000081,0,0); /*  */
	newreact("H,C2=>CH,C",hmrate,4.67e-10,0.5,30450); /* UMIST	 */
	newreact("H+,C2=>C2+,H",hmrate,3.1e-9,0,0); /* UMIST	 */
	newreact("H2+,C2=>C2+,H2",hmrate,1.1e-9,0,0); /* UMIST	 */
	newreact("He+,C2=>C2+,He",hmrate,5.0e-10,0,0); /* UMIST	 */
	newreact("He+,C2H=>C2+,He,H",hmrate,0.000000001,0,0); /*  */
	newreact("He+,C2H=>CH+,C,He",hmrate,0.0000000015,0,0); /*  */
	newreact("He+,C2H=>CH,C+,He",hmrate,0.0000000015,0,0); /*  */
	newreact("H-,C2=>C2H,e-",hmrate,0.00000000235,0,0); /*  */
	newreact("H+,C2H=>C2+,H2",hmrate,0.0000000014,0,0); /*  */
	newreact("H+,C2H=>C2H+,H",hmrate,0.0000000014,0,0); /*  */
	newreact("H2+,C2=>C2H+,H",hmrate,0.0000000015,0,0); /*  */
	newreact("H2+,C2H=>C2H+,H2",hmrate,0.0000000014,0,0); /*  */
	newreact("H3+,C2=>C2H+,H2",hmrate,0.000000002,0,0); /*  */
	newreact("H2,C2+=>C2H+,H",hmrate,0.000000000254,0,0); /*  */
	newreact("H+,C3=>C3+,H",hmrate,0.00000000485,0,0); /*  */
	newreact("He+,C3=>C2,C+,He",hmrate,0.00000000094,0,0); /*  */
	newreact("He+,C2H2=>C2H+,He,H",hmrate,0.000000001,0,0); /*  */
	newreact("He+,C3H=>C3+,He,H",hmrate,0.000000002,0,0); /*  */
	newreact("He+,C2H2=>C2+,He,H2",hmrate,0.000000002,0,0); /*  */
	newreact("He+,C2H2=>CH+,CH,He",hmrate,0.0000000075,0,0); /*  */
	newreact("H-,C2H=>C2H2,e-",hmrate,0.000000000294,0,0); /*  */
	newreact("H+,C2H2=>C2H2+,H",hmrate,0.00000000086,0,0); /*  */
	newreact("H+,C2H2=>C2H+,H2",hmrate,0.00000000077,0,0); /*  */
	newreact("H+,C3H=>C3H+,H",hmrate,0.00000000057,0,0); /*  */
	newreact("H+,C3H=>C3+,H2",hmrate,0.0000000024,0,0); /*  */
	newreact("H2+,C2H=>C2H2+,H",hmrate,0.000000003,0,0); /*  */
	newreact("H2+,C2H2=>C2H2+,H2",hmrate,0.0000000014,0,0); /*  */
	newreact("H3+,C2H=>C2H2+,H2",hmrate,0.000000002,0,0); /*  */
	newreact("H3+,C3=>C3H+,H2",hmrate,0.0000000015,0,0); /*  */
	newreact("H2,C2H+=>C2H2+,H",hmrate,0.000000001,0,0); /*  */
	newreact("H2,C3+=>C3H+,H",hmrate,0.0000000022,0,0); /*  */
	newreact("He+,C2H2=>C2H2+,He",hmrate,0.000000001,0,0); /*  */
	newreact("H,C2H3+=>C2H2+,H2",hmrate,0.0000000000000003,-1,0); /*  */
	newreact("H2+,C2H2=>C2H3+,H",hmrate,0.0000000014,0,0); /*  */
	newreact("H3+,C2H2=>C2H3+,H2",hmrate,0.000000001,0,0); /*  */
	newreact("N2,H3+=>N2H+,H2",hmrate,1.8e-9,0.,0.); /* Bohme et al. 1973     */
	newreact("NH2+,NH2=>NH3+,NH",hmrate,1.0e-9,0,0); /* UMIST   */
	newreact("NH2,OH+=>NH3+,O",hmrate,5.0e-10,0,0); /* UMIST   */
	newreact("NH2+,OH=>H2O+,NH",hmrate,7.1e-10,0,0); /* UMIST   */
	newreact("NH2+,NH3=>NH4+,NH",hmrate,1.61e-9,0,0); /* UMIST   */
	newreact("NH2,NH3+=>NH4+,NH",hmrate,1.0e-11,0,0); /* UMIST   */
	newreact("NH2,CH5+=>NH3+,CH4",hmrate,9.9e-10,0,0); /* UMIST   */
	newreact("NH2+,H2O=>NH3+,OH",hmrate,1.0e-10,0,0); /* UMIST   */
	newreact("NH2+,H2O=>H3O+,NH",hmrate,2.76e-9,0,0); /* UMIST   */
	newreact("NH2,H2O+=>NH3+,OH",hmrate,4.9e-10,0,0); /* UMIST   */
	newreact("NH2+,H2O=>NH4+,O",hmrate,1.45e-10,0,0); /* UMIST   */
	newreact("NH2,H3O+=>H2O,NH3+",hmrate,9.7e-10,0,0); /* UMIST   */
	newreact("NH2,HCN+=>CN,NH3+",hmrate,9.0e-10,0,0); /* UMIST   */
	newreact("NH2,CO+=>HCO+,NH",hmrate,4.5e-10,0,0); /* UMIST   */
	newreact("NH2,HNO+=>NO,NH3+",hmrate,8.8e-10,0,0); /* UMIST   */
	newreact("NH2+,O2=>HNO+,OH",hmrate,2.1e-11,0,0); /* UMIST   */
	newreact("NH2+,S=>HS+,NH",hmrate,4.4e-10,0,0); /* UMIST   */
	newreact("CH4+,NH3=>NH4+,CH3",hmrate,1.15e-9,0,0); /* UMIST   */
	newreact("CH4,NH3+=>NH4+,CH3",hmrate,4.8e-10,0,0); /* UMIST   */
	newreact("CH4,N2+=>N2,CH2+,H2",hmrate,7.e-11,0,0); /* UMIST   */
	newreact("CH4,N2+=>N2,CH3+,H",hmrate,9.3e-10,0,0); /* UMIST   */
	newreact("CH4,HNO+=>NO,CH5+",hmrate,1.0e-10,0,0); /* UMIST   */
	newreact("CH4,S+=>HCS+,H2,H",hmrate,2.0e-11,0,0); /* UMIST   */
	newreact("CH4,CS+=>HCS+,CH3",hmrate,5.0e-10,0,0); /* UMIST   */
	newreact("OH+,NH3=>NH4+,O",hmrate,1.2e-9,0,0); /* UMIST   */
	newreact("OH,NH3+=>NH4+,O",hmrate,7.0e-10,0,0); /* UMIST   */
	newreact("OH+,CN=>HCN+,O",hmrate,1.0e-9,0,0); /* UMIST   */
	newreact("OH,HCN+=>CN,H2O+",hmrate,6.3e-10,0,0); /* UMIST   */
	newreact("OH+,NO=>HNO+,O",hmrate,6.11e-10,0,0); /* UMIST   */
	newreact("OH,HNO+=>NO,H2O+",hmrate,6.2e-10,0,0); /* UMIST   */
	newreact("OH+,S=>HS+,O",hmrate,4.3e-10,0,0); /* UMIST   */
	newreact("OH+,S=>SO+,H",hmrate,4.3e-10,0,0); /* UMIST   */
	newreact("OH,S+=>SO+,H",hmrate,6.1e-10,0,0); /* UMIST   */
	newreact("NH3+,NH3=>NH4+,NH2",hmrate,2.2e-9,0,0); /* UMIST   */
	newreact("NH3,CH5+=>NH4+,CH4",hmrate,2.5e-9,0,0); /* UMIST   */
	newreact("NH3+,H2O=>NH4+,OH",hmrate,1.1e-10,0,0); /* UMIST   */
	newreact("NH3,H2O+=>NH4+,OH",hmrate,9.45e-10,0,0); /* UMIST   */
	newreact("NH3,H3O+=>NH4+,H2O",hmrate,2.2e-9,0,0); /* UMIST   */
	newreact("NH3,CO+=>HCO+,NH2",hmrate,4.12e-11,0,0); /* UMIST   */
	newreact("NH3,HNO+=>NO,NH4+",hmrate,1.1e-9,0,0); /* UMIST   */
	newreact("NH3,HS+=>S,NH4+",hmrate,9.75e-10,0,0); /* UMIST   */
	newreact("NH3,HCS+=>CS,NH4+",hmrate,2.0e-9,0,0); /* UMIST   */
	newreact("CH5+,S=>HS+,CH4",hmrate,1.3e-9,0,0); /* UMIST   */
	newreact("H2O,CN+=>HCN+,OH",hmrate,1.6e-9,0,0); /* UMIST   */
	newreact("H2O,CN+=>HCO+,NH",hmrate,1.6e-10,0,0); /* UMIST   */
	newreact("H2O,HCN+=>CN,H3O+",hmrate,1.8e-9,0,0); /* UMIST   */
	newreact("H2O,HNO+=>NO,H3O+",hmrate,2.3e-9,0,0); /* UMIST   */
	newreact("H2O+,S=>HS+,OH",hmrate,4.3e-10,0,0); /* UMIST   */
	newreact("H2O,HS+=>S,H3O+",hmrate,7.8e-10,0,0); /* UMIST   */
	newreact("H3O+,CS=>HCS+,H2O",hmrate,1.0e-9,0,0); /* UMIST   */
	newreact("CN+,NO=>OCN+,N",hmrate,1.9e-10,0,0); /* UMIST   */
	newreact("CN,HNO+=>NO,HCN+",hmrate,8.7e-10,0,0); /* UMIST   */
	newreact("CN+,O2=>OCN+,O",hmrate,8.6e-11,0,0); /* UMIST   */
	newreact("HCN+,S=>HS+,CN",hmrate,5.7e-10,0,0); /* UMIST   */
	newreact("HNO+,S=>HS+,NO",hmrate,1.1e-9,0,0); /* UMIST   */
	newreact("O2,S+=>SO+,O",hmrate,1.5e-11,0,0); /* UMIST   */
	newreact("O2+,S=>SO+,O",hmrate,5.4e-10,0,0); /* UMIST   */
	newreact("O2,CS+=>OCS+,O",hmrate,1.3e-10,0,0); /* UMIST   */
	newreact("S,SiO+=>SO,Si+",hmrate,1.0e-9,0,0); /* UMIST   */
	newreact("C+,NH3=>NH3+,C",hmrate,5.06e-10,0,0); /* UMIST   */
	newreact("C,CN+=>CN,C+",hmrate,1.1e-10,0,0); /* UMIST   */
	newreact("C,N2+=>N2,C+",hmrate,1.1e-10,0,0); /* UMIST   */
	newreact("C+,NO=>NO+,C",hmrate,5.2e-10,0,0); /* UMIST   */
	newreact("C+,SiN=>SiN+,C",hmrate,1.0e-9,0,0); /* UMIST   */
	newreact("C,CS+=>CS,C+",hmrate,1.6e-9,0,0); /* UMIST   */
	newreact("C+,NS=>NS+,C",hmrate,7.6e-10,0,0); /* UMIST   */
	newreact("C+,SO=>SO+,C",hmrate,2.6e-10,0,0); /* UMIST   */
	newreact("C+,OCS=>OCS+,C",hmrate,4.0e-10,0,0); /* UMIST   */
	newreact("CH,NH2+=>NH2,CH+",hmrate,3.5e-10,0,0); /* UMIST   */
	newreact("CH+,NH3=>NH3+,CH",hmrate,4.59e-10,0,0); /* UMIST   */
	newreact("CH,CN+=>CN,CH+",hmrate,6.4e-10,0,0); /* UMIST   */
	newreact("CH,N2+=>N2,CH+",hmrate,6.3e-10,0,0); /* UMIST   */
	newreact("CH+,NO=>NO+,CH",hmrate,7.6e-10,0,0); /* UMIST   */
	newreact("N+,NH=>NH+,N",hmrate,3.7e-10,0,0); /* UMIST   */
	newreact("N+,NH2=>NH2+,N",hmrate,1.0e-9,0,0); /* UMIST   */
	newreact("N+,NH3=>NH3+,N",hmrate,1.97e-9,0,0); /* UMIST   */
	newreact("N+,CN=>CN+,N",hmrate,1.1e-9,0,0); /* UMIST   */
	newreact("N+,HCN=>HCN+,N",hmrate,1.2e-9,0,0); /* UMIST   */
	newreact("N,N2+=>N2,N+",hmrate,1.0e-11,0,0); /* UMIST   */
	newreact("N+,NO=>NO+,N",hmrate,4.51e-10,0,0); /* UMIST   */
	newreact("N+,OCS=>OCS+,N",hmrate,1.02e-9,0,0); /* UMIST   */
	newreact("CH2,NH2+=>NH2,CH2+",hmrate,4.9e-10,0,0); /* UMIST   */
	newreact("CH2,CN+=>CN,CH2+",hmrate,8.8e-10,0,0); /* UMIST   */
	newreact("CH2,N2+=>N2,CH2+",hmrate,8.7e-10,0,0); /* UMIST   */
	newreact("CH2+,NO=>NO+,CH2",hmrate,4.2e-10,0,0); /* UMIST   */
	newreact("NH,O+=>O,NH+",hmrate,3.6e-10,0,0); /* UMIST   */
	newreact("NH,OH+=>OH,NH+",hmrate,3.6e-10,0,0); /* UMIST   */
	newreact("NH+,NH3=>NH3+,NH",hmrate,1.8e-9,0,0); /* UMIST   */
	newreact("NH+,H2O=>H2O+,NH",hmrate,1.05e-9,0,0); /* UMIST   */
	newreact("NH,CN+=>CN,NH+",hmrate,6.5e-10,0,0); /* UMIST   */
	newreact("NH,N2+=>N2,NH+",hmrate,6.5e-10,0,0); /* UMIST   */
	newreact("NH+,NO=>NO+,NH",hmrate,7.12e-10,0,0); /* UMIST   */
	newreact("NH+,O2=>O2+,NH",hmrate,4.51e-10,0,0); /* UMIST   */
	newreact("NH+,S=>S+,NH",hmrate,6.9e-10,0,0); /* UMIST   */
	newreact("CH3+,NO=>NO+,CH3",hmrate,1.0e-9,0,0); /* UMIST   */
	newreact("O+,NH2=>NH2+,O",hmrate,1.0e-9,0,0); /* UMIST   */
	newreact("O+,NH3=>NH3+,O",hmrate,1.2e-9,0,0); /* UMIST   */
	newreact("O,CN+=>CN,O+",hmrate,6.5e-11,0,0); /* UMIST   */
	newreact("O,HCN+=>HCN,O+",hmrate,6.5e-11,0,0); /* UMIST   */
	newreact("O,N2+=>N2,O+",hmrate,1.0e-11,0,0); /* UMIST   */
	newreact("O+,NO=>NO+,O",hmrate,1.7e-12,0,0); /* UMIST   */
	newreact("O+,OCS=>OCS+,O",hmrate,6.5e-10,0,0); /* UMIST   */
	newreact("NH2,OH+=>OH,NH2+",hmrate,5.0e-10,0,0); /* UMIST   */
	newreact("NH2+,NH3=>NH3+,NH2",hmrate,6.9e-10,0,0); /* UMIST   */
	newreact("NH2,H2O+=>H2O,NH2+",hmrate,4.9e-10,0,0); /* UMIST   */
	newreact("NH2,CN+=>CN,NH2+",hmrate,9.1e-10,0,0); /* UMIST   */
	newreact("NH2,N2+=>N2,NH2+",hmrate,8.9e-10,0,0); /* UMIST   */
	newreact("NH2+,NO=>NO+,NH2",hmrate,7.0e-10,0,0); /* UMIST   */
	newreact("NH2,O2+=>O2,NH2+",hmrate,8.7e-10,0,0); /* UMIST   */
	newreact("NH2+,S=>S+,NH2",hmrate,4.4e-10,0,0); /* UMIST   */
	newreact("CH4+,NH3=>NH3+,CH4",hmrate,1.65e-9,0,0); /* UMIST   */
	newreact("CH4+,OCS=>OCS+,CH4",hmrate,4.2e-10,0,0); /* UMIST   */
	newreact("OH+,NH3=>NH3+,OH",hmrate,1.2e-9,0,0); /* UMIST   */
	newreact("OH,CN+=>CN,OH+",hmrate,6.4e-10,0,0); /* UMIST   */
	newreact("OH,N2+=>N2,OH+",hmrate,6.3e-10,0,0); /* UMIST   */
	newreact("OH+,NO=>NO+,OH",hmrate,3.59e-10,0,0); /* UMIST   */
	newreact("NH3,H2O+=>H2O,NH3+",hmrate,2.21e-9,0,0); /* UMIST   */
	newreact("NH3,HCN+=>HCN,NH3+",hmrate,1.68e-9,0,0); /* UMIST   */
	newreact("NH3,N2+=>N2,NH3+",hmrate,1.9e-9,0,0); /* UMIST   */
	newreact("NH3+,Si=>Si+,NH3",hmrate,1.9e-9,0,0); /* UMIST   */
	newreact("NH3+,NO=>NO+,NH3",hmrate,7.2e-10,0,0); /* UMIST   */
	newreact("NH3,O2+=>O2,NH3+",hmrate,2.0e-9,0,0); /* UMIST   */
	newreact("NH3,S+=>S,NH3+",hmrate,1.44e-9,0,0); /* UMIST   */
	newreact("NH3,HS+=>HS,NH3+",hmrate,5.25e-10,0,0); /* UMIST   */
	newreact("NH3,SO+=>SO,NH3+",hmrate,1.3e-9,0,0); /* UMIST   */
	newreact("H2O,HCN+=>HCN,H2O+",hmrate,1.8e-9,0,0); /* UMIST   */
	newreact("H2O,N2+=>N2,H2O+",hmrate,2.3e-9,0,0); /* UMIST   */
	newreact("H2O+,NO=>NO+,H2O",hmrate,2.7e-10,0,0); /* UMIST   */
	newreact("CN+,HCN=>HCN+,CN",hmrate,1.79e-9,0,0); /* UMIST   */
	newreact("CN,N2+=>N2,CN+",hmrate,1.0e-10,0,0); /* UMIST   */
	newreact("CN+,NO=>NO+,CN",hmrate,5.7e-10,0,0); /* UMIST   */
	newreact("CN+,O2=>O2+,CN",hmrate,2.58e-10,0,0); /* UMIST   */
	newreact("CN+,S=>S+,CN",hmrate,1.1e-9,0,0); /* UMIST   */
	newreact("HCN,N2+=>N2,HCN+",hmrate,3.9e-10,0,0); /* UMIST   */
	newreact("HCN+,NO=>NO+,HCN",hmrate,8.1e-10,0,0); /* UMIST   */
	newreact("HCN+,O2=>O2+,HCN",hmrate,3.2e-10,0,0); /* UMIST   */
	newreact("HCN+,S=>S+,HCN",hmrate,5.7e-10,0,0); /* UMIST   */
	newreact("N2+,NO=>NO+,N2",hmrate,4.4e-10,0,0); /* UMIST   */
	newreact("N2+,O2=>O2+,N2",hmrate,5.e-11,0,0); /* UMIST   */
	newreact("N2+,S=>S+,N2",hmrate,1.1e-9,0,0); /* UMIST   */
	newreact("Si,NO+=>NO,Si+",hmrate,1.6e-9,0,0); /* UMIST   */
	newreact("Si,HS+=>HS,Si+",hmrate,1.4e-9,0,0); /* UMIST   */
	newreact("Si,CS+=>CS,Si+",hmrate,1.5e-10,0,0); /* UMIST   */
	newreact("NO,HNO+=>HNO,NO+",hmrate,7.0e-10,0,0); /* UMIST   */
	newreact("NO,O2+=>O2,NO+",hmrate,4.5e-10,0,0); /* UMIST   */
	newreact("NO,S+=>S,NO+",hmrate,3.7e-10,0,0); /* UMIST   */
	newreact("NO,HS+=>HS,NO+",hmrate,4.5e-10,0,0); /* UMIST   */
	newreact("NO,SiO+=>SiO,NO+",hmrate,7.2e-10,0,0); /* UMIST   */
	newreact("NO,S2+=>S2,NO+",hmrate,5.1e-10,0,0); /* UMIST   */
	newreact("O2+,NO2=>NO2+,O2",hmrate,6.6e-10,0,0); /* UMIST   */
	newreact("S,HS+=>HS,S+",hmrate,9.7e-10,0,0); /* UMIST   */
	newreact("C,N=>CN,PHOTON",hmrate,1.e-17,0,0); /* UMIST   */
	newreact("C,S=>CS,PHOTON",hmrate,4.36E-19,0.22,0); /* UMIST   */
	newreact("C+,S=>CS+,PHOTON",hmrate,3.07E-19,0.15,0); /* UMIST   */
	newreact("N+,N=>N2+,PHOTON",hmrate,3.71E-18,0.24,26.1); /* UMIST   */
	newreact("CH,N+=>N,CH+",hmrate,3.6e-10,0,0); /* UMIST   */
	newreact("CH+,S=>S+,CH",hmrate,4.7e-10,0,0); /* UMIST   */
	newreact("N+,CH2=>CH2+,N",hmrate,1.0e-9,0,0); /* UMIST   */
	newreact("N+,CH4=>CH4+,N",hmrate,2.8e-11,0,0); /* UMIST   */
	newreact("N+,OH=>OH+,N",hmrate,3.7e-10,0,0); /* UMIST   */
	newreact("N+,H2O=>H2O+,N",hmrate,2.8e-9,0,0); /* UMIST   */
	newreact("N+,O2=>O2+,N",hmrate,3.11e-10,0,0); /* UMIST   */
	newreact("OH+,S=>S+,OH",hmrate,4.3e-10,0,0); /* UMIST   */
	newreact("H2O+,S=>S+,H2O",hmrate,4.3e-10,0,0); /* UMIST   */
	newreact("Si,S+=>S,Si+",hmrate,1.6e-9,0,0); /* UMIST   */
	newreact("O2+,S=>S+,O2",hmrate,5.4e-10,0,0); /* UMIST   */
	newreact("O,CCl=>ClO,C",hmrate,0.000000000138,0,16050); /* UMIST   */
	newreact("O,ClO=>Cl,O2",hmrate,0.000000000038,0,0); /* UMIST   */
	newreact("C+,HCl=>CCl+,H",hmrate,0.0000000011,0,0); /* UMIST   */
	newreact("CH3+,HCl=>H2CCl+,H2",hmrate,0.0000000013,0,0); /* UMIST   */
	newreact("H2O,H2Cl+=>HCl,H3O+",hmrate,0.000000002,0,0); /* UMIST   */
	newreact("C+,CCl=>CCl+,C",hmrate,0.000000001,0,0); /* UMIST   */
	newreact("C+,ClO=>ClO+,C",hmrate,0.000000001,0,0); /* UMIST   */
	newreact("O2,Cl+=>Cl,O2+",hmrate,0.00000000046,0,0); /* UMIST   */
	newreact("C,NH2=>HNC,H",hmrate,0.0000000000326,-0.1,-9); /*  */
	newreact("N,CH2=>HNC,H",hmrate,0.0000000000789,0.17,0); /*  */
	newreact("CH+,HNC=>HCNH+,C",hmrate,0.0000000018,0,0); /*  */
	newreact("CH,HCNH+=>HNC,CH2+",hmrate,0.000000000315,0,0); /*  */
	newreact("CH2,HCNH+=>HNC,CH3+",hmrate,0.000000000435,0,0); /*  */
	newreact("NH+,HNC=>HCNH+,N",hmrate,0.0000000018,0,0); /*  */
	newreact("NH2+,HNC=>HCNH+,NH",hmrate,0.0000000012,0,0); /*  */
	newreact("NH2,HCNH+=>HNC,NH3+",hmrate,0.000000000445,0,0); /*  */
	newreact("OH+,HNC=>HCNH+,O",hmrate,0.0000000012,0,0); /*  */
	newreact("NH3,HCNH+=>HNC,NH4+",hmrate,0.0000000011,0,0); /*  */
	newreact("CH5+,HNC=>HCNH+,CH4",hmrate,0.0000000012,0,0); /*  */
	newreact("H2O+,HNC=>HCNH+,OH",hmrate,0.0000000011,0,0); /*  */
	newreact("H3O+,HNC=>HCNH+,H2O",hmrate,0.000000004,0,0); /*  */
	newreact("HCN+,HNC=>HCNH+,CN",hmrate,0.000000001,0,0); /*  */
	newreact("HNC,HNO+=>NO,HCNH+",hmrate,0.00000000099,0,0); /*  */
	newreact("HNC,HS+=>S,HCNH+",hmrate,0.00000000086,0,0); /*  */
	newreact("CH+,HCN=>HCNH+,C",hmrate,0.0000000018,0,0); /*  */
	newreact("CH,HCNH+=>HCN,CH2+",hmrate,0.000000000315,0,0); /*  */
	newreact("N+,CH4=>HCNH+,H,H",hmrate,0.00000000038,0,0); /*  */
	newreact("CH2,HCNH+=>HCN,CH3+",hmrate,0.000000000435,0,0); /*  */
	newreact("NH,CH3+=>HCNH+,H2",hmrate,0.00000000074,0,0); /*  */
	newreact("NH+,HCN=>HCNH+,N",hmrate,0.0000000018,0,0); /*  */
	newreact("NH2+,HCN=>HCNH+,NH",hmrate,0.0000000012,0,0); /*  */
	newreact("NH2,HCNH+=>HCN,NH3+",hmrate,0.000000000445,0,0); /*  */
	newreact("CH4,HCN+=>HCNH+,CH3",hmrate,0.00000000104,0,0); /*  */
	newreact("OH+,HCN=>HCNH+,O",hmrate,0.0000000012,0,0); /*  */
	newreact("NH3,HCN+=>HCNH+,NH2",hmrate,0.00000000084,0,0); /*  */
	newreact("NH3,HCNH+=>HCN,NH4+",hmrate,0.0000000011,0,0); /*  */
	newreact("CH5+,HCN=>HCNH+,CH4",hmrate,0.0000000012,0,0); /*  */
	newreact("H2O+,HCN=>HCNH+,OH",hmrate,0.0000000011,0,0); /*  */
	newreact("H3O+,HCN=>HCNH+,H2O",hmrate,0.000000004,0,0); /*  */
	newreact("HCN+,HCN=>HCNH+,CN",hmrate,0.0000000016,0,0); /*  */
	newreact("HCN,HNO+=>NO,HCNH+",hmrate,0.00000000099,0,0); /*  */
	newreact("HCN,HS+=>S,HCNH+",hmrate,0.00000000086,0,0); /*  */
	newreact("C,CH2=>C2H,H",hmrate,0.00000000005,0.5,0); /*  */
	newreact("O+,C2H=>CO+,CH",hmrate,0.00000000046,0,0); /*  */
	newreact("C2H,CO+=>HCO+,C2",hmrate,0.00000000039,0,0); /*  */
	newreact("C,CH2+=>C2H+,H",hmrate,0.0000000012,0,0); /*  */
	newreact("C,CH3+=>C2H+,H2",hmrate,0.0000000012,0,0); /*  */
	newreact("C2,HCN+=>CN,C2H+",hmrate,0.00000000084,0,0); /*  */
	newreact("C2,HNO+=>NO,C2H+",hmrate,0.00000000082,0,0); /*  */
	newreact("C2H,CN+=>CN,C2H+",hmrate,0.0000000008,0,0); /*  */
	newreact("C2H,N2+=>N2,C2H+",hmrate,0.00000000079,0,0); /*  */
	newreact("C2H+,HCN=>HCNH+,C2",hmrate,0.0000000014,0,0); /*  */
	newreact("C2H+,HNC=>HCNH+,C2",hmrate,0.0000000014,0,0); /*  */
	newreact("C2H+,NO=>NO+,C2H",hmrate,0.00000000012,0,0); /*  */
	newreact("C2H+,S=>S+,C2H",hmrate,0.0000000012,0,0); /*  */
	newreact("CH,C2H+=>C2,CH2+",hmrate,0.00000000032,0,0); /*  */
	newreact("CH2,C2H+=>C2,CH3+",hmrate,0.00000000044,0,0); /*  */
	newreact("CH4,C2+=>C2H+,CH3",hmrate,0.000000000238,0,0); /*  */
	newreact("CH5+,C2=>C2H+,CH4",hmrate,0.00000000095,0,0); /*  */
	newreact("CH+,CH2=>C2H+,H2",hmrate,0.000000001,0,0); /*  */
	newreact("C+,CH2=>C2H+,H",hmrate,0.000000000434,-0.5,0); /*  */
	newreact("C+,CH3=>C2H+,H2",hmrate,0.000000001,0,0); /*  */
	newreact("H2O,C2+=>C2H+,OH",hmrate,0.00000000044,0,0); /*  */
	newreact("H2O+,C2=>C2H+,OH",hmrate,0.00000000047,0,0); /*  */
	newreact("H2O+,C2H=>C2H+,H2O",hmrate,0.00000000044,0,0); /*  */
	newreact("H3O+,C2=>C2H+,H2O",hmrate,0.00000000092,0,0); /*  */
	newreact("N,C2H+=>CN,CH+",hmrate,0.00000000009,0,0); /*  */
	newreact("NH,C2+=>C2H+,N",hmrate,0.00000000033,0,0); /*  */
	newreact("NH2,C2H+=>C2,NH3+",hmrate,0.00000000046,0,0); /*  */
	newreact("NH2+,C2=>C2H+,NH",hmrate,0.00000000097,0,0); /*  */
	newreact("NH3,C2H+=>C2,NH4+",hmrate,0.00000000055,0,0); /*  */
	newreact("NH+,C2=>C2H+,N",hmrate,0.00000000049,0,0); /*  */
	newreact("N+,C2H=>C2H+,N",hmrate,0.00000000095,0,0); /*  */
	newreact("O,C2H+=>HCO+,C",hmrate,0.00000000033,0,0); /*  */
	newreact("OH+,C2=>C2H+,O",hmrate,0.00000000048,0,0); /*  */
	newreact("OH+,C2H=>C2H+,OH",hmrate,0.00000000045,0,0); /*  */
	newreact("O+,C2H=>C2H+,O",hmrate,0.00000000046,0,0); /*  */
	newreact("C+,C2H=>C3+,H",hmrate,0.000000001,0,0); /*  */
	newreact("CH+,C2H=>C3+,H2",hmrate,0.00000000098,0,0); /*  */
	newreact("C2+,C2=>C3+,C",hmrate,0.00000000087,0,0); /*  */
	newreact("CH+,C2=>C3+,H",hmrate,0.000000001,0,0); /*  */
	newreact("CH,C2+=>C3+,H",hmrate,0.00000000032,0,0); /*  */
	newreact("C,C2H+=>C3+,H",hmrate,0.0000000011,0,0); /*  */
	newreact("C,C2H2=>C3H,H",hmrate,0.0000000002,0,0); /*  */
	newreact("C,C2H2+=>C3H+,H",hmrate,0.0000000011,0,0); /*  */
	newreact("C2H,HCN+=>C2H2+,CN",hmrate,0.00000000079,0,0); /*  */
	newreact("C2H2+,HCN=>HCNH+,C2H",hmrate,0.00000000023,0,0); /*  */
	newreact("C2H2+,NO=>NO+,C2H2",hmrate,0.00000000012,0,0); /*  */
	newreact("C2H+,HCN=>C2H2+,CN",hmrate,0.0000000014,0,0); /*  */
	newreact("CH,C2H+=>C3H+,H",hmrate,0.00000000032,0,0); /*  */
	newreact("CH,CH3+=>C2H2+,H2",hmrate,0.00000000071,0,0); /*  */
	newreact("CH2,C2+=>C3H+,H",hmrate,0.00000000045,0,0); /*  */
	newreact("CH2,CH2=>C2H2,H,H",hmrate,0.00000000018,0,400); /*  */
	newreact("CH2,CH2=>C2H2,H2",hmrate,0.00000000263,0,6013); /*  */
	newreact("CH3+,C2=>C3H+,H2",hmrate,0.00000000099,0,0); /*  */
	newreact("CH4,C2H+=>C2H2+,CH3",hmrate,0.000000000374,0,0); /*  */
	newreact("CH4,C2+=>C2H2+,CH2",hmrate,0.000000000182,0,0); /*  */
	newreact("CH4,C2+=>C3H+,H2,H",hmrate,0.000000000196,0,0); /*  */
	newreact("CH4+,C2H2=>C2H2+,CH4",hmrate,0.00000000113,0,0); /*  */
	newreact("CH5+,C2H=>C2H2+,CH4",hmrate,0.0000000009,0,0); /*  */
	newreact("CH+,CH4=>C2H2+,H2,H",hmrate,0.000000000143,0,0); /*  */
	newreact("C+,C2H2=>C3H+,H",hmrate,0.0000000022,0,0); /*  */
	newreact("C+,CH3=>C2H2+,H",hmrate,0.0000000013,0,0); /*  */
	newreact("C+,CH4=>C2H2+,H2",hmrate,0.0000000004,0,0); /*  */
	newreact("H2O,C2H2+=>C2H,H3O+",hmrate,0.00000000022,0,0); /*  */
	newreact("H2O,C3H+=>HCO+,C2H2",hmrate,0.00000000027,0,0); /*  */
	newreact("H2O+,C2H=>C2H2+,OH",hmrate,0.00000000044,0,0); /*  */
	newreact("H2O+,C2H2=>C2H2+,H2O",hmrate,0.0000000019,0,0); /*  */
	newreact("H3O+,C3=>C3H+,H2O",hmrate,0.000000002,0,0); /*  */
	newreact("N,C2H2+=>HCN,CH+",hmrate,0.000000000025,0,0); /*  */
	newreact("NH2,C2H2+=>C2H,NH3+",hmrate,0.00000000045,0,0); /*  */
	newreact("NH2+,C2H=>C2H2+,NH",hmrate,0.00000000091,0,0); /*  */
	newreact("NH3,C2H2+=>C2H,NH4+",hmrate,0.0000000011,0,0); /*  */
	newreact("NH3,C2H2+=>C2H2,NH3+",hmrate,0.0000000021,0,0); /*  */
	newreact("NH3,C3H+=>C3,NH4+",hmrate,0.0000000008,0,0); /*  */
	newreact("NH3,C3H+=>C3H,NH3+",hmrate,0.00000000032,0,0); /*  */
	newreact("NH3+,C2=>C2H2+,NH",hmrate,0.00000000001,0,0); /*  */
	newreact("NH+,C2H=>C2H2+,N",hmrate,0.0000000014,0,0); /*  */
	newreact("NO,C3H+=>C3H,NO+",hmrate,0.00000000013,0,0); /*  */
	newreact("O,C2H2=>C2H,OH",hmrate,0.0000000053,0,8520); /*  */
	newreact("O,C2H2+=>HCO+,CH",hmrate,0.000000000085,0,0); /*  */
	newreact("OH,C2H2=>C2H,H2O",hmrate,0.000000000000105,2.68,6060); /*  */
	newreact("OH+,C2H=>C2H2+,O",hmrate,0.00000000045,0,0); /*  */
	newreact("O+,C2H2=>C2H2+,O",hmrate,0.000000000039,0,0); /*  */
	newreact("C,C2H3+=>C3H+,H2",hmrate,0.000000001,0,0); /*  */
	newreact("C2H,C2H3+=>C2H2+,C2H2",hmrate,0.00000000033,0,0); /*  */
	newreact("CH2,CH3+=>C2H3+,H2",hmrate,0.00000000099,0,0); /*  */
	newreact("CH4,C3H+=>C2H3+,C2H2",hmrate,0.000000000612,0,0); /*  */
	newreact("CH4,HCN+=>C2H3+,NH2",hmrate,0.00000000026,0,0); /*  */
	newreact("CH4+,C2H2=>C2H3+,CH3",hmrate,0.00000000125,0,0); /*  */
	newreact("CH5+,C2H2=>C2H3+,CH4",hmrate,0.0000000016,0,0); /*  */
	newreact("CH+,CH4=>C2H3+,H2",hmrate,0.00000000109,0,0); /*  */
	newreact("C+,CH4=>C2H3+,H",hmrate,0.0000000011,0,0); /*  */
	newreact("H2O,C2H3+=>C2H2,H3O+",hmrate,0.00000000111,0,0); /*  */
	newreact("HCN,C2H3+=>HCNH+,C2H2",hmrate,0.0000000029,0,0); /*  */
	newreact("HNC,C2H3+=>HCNH+,C2H2",hmrate,0.0000000029,0,0); /*  */
	newreact("NH3,C2H3+=>C2H2,NH4+",hmrate,0.0000000025,0,0); /*  */
	newreact("O,C2H3+=>HCO+,CH2",hmrate,0.0000000001,0,0); /*  */
	newreact("C,CH=>C2,H",hmrate,6.59e-11,0,0); /* UMIST	 */
	newreact("C,CN=>C2,N",hmrate,4.98e-10,0,18116); /* UMIST	 */
	newreact("C,CS=>S,C2",hmrate,1.44e-11,0.5,20435); /* UMIST	 */
	newreact("C2,S=>CS,C",hmrate,1.73e-11,0.5,0); /* UMIST	 */
	newreact("NH+,C2=>HCN+,C",hmrate,4.9e-10,0,0); /* UMIST	 */
	newreact("O+,C2=>CO+,C",hmrate,4.8e-10,0,0); /* UMIST	 */
	newreact("C+,S=>C,S+",hmrate,1.0e-9,0,0); /* UMIST	 */
	newreact("C2,S+=>CS+,C",hmrate,8.1e-10,0,0); /* UMIST	 */
	newreact("C,C2+=>C2,C+",hmrate,1.1e-10,0,0); /* UMIST	 */
	newreact("CH,C2+=>C2,CH+",hmrate,3.2e-10,0,0); /* UMIST	 */
	newreact("N+,C2=>C2+,N",hmrate,1.0e-9,0,0); /* UMIST	 */
	newreact("CH2,C2+=>C2,CH2+",hmrate,4.5e-10,0,0); /* UMIST	 */
	newreact("O+,C2=>C2+,O",hmrate,4.8e-10,0,0); /* UMIST	 */
	newreact("NH2,C2+=>C2,NH2+",hmrate,4.6e-10,0,0); /* UMIST	 */
	newreact("OH+,C2=>C2+,OH",hmrate,4.8e-10,0,0); /* UMIST	 */
	newreact("OH,C2+=>C2,OH+",hmrate,6.5e-10,0,0); /* UMIST	 */
	newreact("H2O+,C2=>C2+,H2O",hmrate,4.7e-10,0,0); /* UMIST	 */
	newreact("C2,CN+=>CN,C2+",hmrate,8.5e-10,0,0); /* UMIST	 */
	newreact("C2,N2+=>N2,C2+",hmrate,8.4e-10,0,0); /* UMIST	 */
	newreact("C2+,NO=>NO+,C2",hmrate,3.4e-10,0,0); /* UMIST	 */
	newreact("C2,O2+=>O2,C2+",hmrate,4.1e-10,0,0); /* UMIST	 */
	newreact("C2+,S=>S+,C2",hmrate,5.8e-10,0,0); /* UMIST	 */
	newreact("C,C=>C2,PHOTON",hmrate,4.36E-18,0.35,161.3); /* UMIST	 */
	newreact("C,CH+=>C2+,H",hmrate,1.2e-9,0,0); /* UMIST	 */
	newreact("CH+,CH=>C2+,H2",hmrate,7.4e-10,0,0); /* UMIST	 */
	newreact("N,C2+=>CN,C+",hmrate,4.0e-11,0,0); /* UMIST	 */
	newreact("O,C2+=>CO+,C",hmrate,3.1e-10,0,0); /* UMIST	 */
	newreact("C2+,S=>CS+,C",hmrate,5.8e-10,0,0); /* UMIST	 */
	newreact("C+,C=>C2+,PHOTON",hmrate,4.01E-18,0.17,101.5); /* UMIST	 */
	newreact("C,CH2=>CH,CH",hmrate,2.69e-12,0,23550); /* UMIST	 */
	newreact("C,H2O+=>OH,CH+",hmrate,1.1e-9,0,0); /* UMIST	 */
	newreact("C,H3O+=>HCO+,H2",hmrate,1.0e-11,0,0); /* UMIST	 */
	newreact("C,O2+=>O2,C+",hmrate,5.2e-11,0,0); /* UMIST	 */
	newreact("C,O2+=>CO+,O",hmrate,5.2e-11,0,0); /* UMIST	 */
	newreact("C,OH=>O,CH",hmrate,2.25e-11,0.5,14800); /* UMIST	 */
	newreact("C,OH+=>O,CH+",hmrate,1.2e-9,0,0); /* UMIST	 */
	newreact("C+,CH2=>CH2+,C",hmrate,5.2e-10,0,0); /* UMIST	 */
	newreact("C+,H2O=>HCO+,H",hmrate,9.0e-10,0,0); /* UMIST	 */
	newreact("C+,O=>CO+,PHOTON",hmrate,2.5E-18,0,0); /* UMIST	 */
	newreact("C+,O2=>CO+,O",hmrate,3.8e-10,0,0); /* UMIST	 */
	newreact("O,CH=>OH,C",hmrate,2.52e-11,0,2381); /* UMIST	 */
	newreact("O,CH=>HCO+,e-",hmrate,2.0e-11,0.44,0); /* UMIST	 */
	newreact("O,CH+=>CO+,H",hmrate,3.5e-10,0,0); /* UMIST	 */
	newreact("O,CH2=>OH,CH",hmrate,4.98e-10,0,6000); /* UMIST	 */
	newreact("O,CH2+=>HCO+,H",hmrate,7.5e-10,0,0); /* UMIST	 */
	newreact("O,H2O=>OH,OH",hmrate,1.85e-11,0.95,8571); /* UMIST	 */
	newreact("O,H2O+=>O2+,H2",hmrate,4.0e-11,0,0); /* UMIST	 */
	newreact("O,O=>O2,PHOTON",hmrate,4.9E-20,1.58,0); /* UMIST	 */
	newreact("O,OH=>O2,H",hmrate,4.34e-11,-0.5,30); /* UMIST	 */
	newreact("O,OH+=>O2+,H",hmrate,7.1e-10,0,0); /* UMIST	 */
	newreact("O,Si=>SiO,PHOTON",hmrate,5.52E-18,0.31,0); /* UMIST	 */
	newreact("O,Si+=>SiO+,PHOTON",hmrate,1.e-17,0,0); /* UMIST	 */
	newreact("O,SiH=>SiO,H",hmrate,4.0e-11,0.5,0); /* UMIST	 */
	newreact("O,SiH2+=>SiOH+,H",hmrate,6.3e-10,0,0); /* UMIST	 */
	newreact("O,SiO+=>O2,Si+",hmrate,2.0e-10,0,0); /* UMIST	 */
	newreact("O+,CH=>O,CH+",hmrate,3.5e-10,0,0); /* UMIST	 */
	newreact("O+,CH=>CO+,H",hmrate,3.5e-10,0,0); /* UMIST	 */
	newreact("O+,CH2=>O,CH2+",hmrate,9.7e-10,0,0); /* UMIST	 */
	newreact("O+,H2O=>H2O+,O",hmrate,3.2e-9,0,0); /* UMIST	 */
	newreact("O+,O2=>O2+,O",hmrate,1.9e-11,0,0); /* UMIST	 */
	newreact("O+,OH=>O2+,H",hmrate,3.6e-10,0,0); /* UMIST	 */
	newreact("O+,OH=>OH+,O",hmrate,3.6e-10,0,0); /* UMIST	 */
	newreact("Si,CH+=>Si+,CH",hmrate,2.0e-10,0,0); /* UMIST	 */
	newreact("Si,H2O+=>Si+,H2O",hmrate,3.0e-9,0,0); /* UMIST	 */
	newreact("Si,OH=>SiO,H",hmrate,2.0e-10,0.5,0); /* UMIST	 */
	newreact("Si,O2+=>O2,Si+",hmrate,1.6e-9,0,0); /* UMIST	 */
	newreact("Si+,H2O=>SiOH+,H",hmrate,2.3e-10,0,0); /* UMIST	 */
	newreact("Si+,OH=>SiO+,H",hmrate,6.3e-10,0,0); /* UMIST	 */
	newreact("Si+,O2=>SiO+,O",hmrate,1.e-13,0,0); /* UMIST	 */
	newreact("CH,CO+=>HCO+,C",hmrate,3.2e-10,0,0); /* UMIST	 */
	newreact("CH,H2O+=>H2O,CH+",hmrate,3.4e-10,0,0); /* UMIST	 */
	newreact("CH,H2O+=>OH,CH2+",hmrate,3.4e-10,0,0); /* UMIST	 */
	newreact("CH,H3O+=>H2O,CH2+",hmrate,6.8e-10,0,0); /* UMIST	 */
	newreact("CH,O2+=>O2,CH+",hmrate,3.1e-10,0,0); /* UMIST	 */
	newreact("CH,O2+=>HCO+,O",hmrate,3.1e-10,0,0); /* UMIST	 */
	newreact("CH,OH+=>OH,CH+",hmrate,3.5e-10,0,0); /* UMIST	 */
	newreact("CH,OH+=>O,CH2+",hmrate,3.5e-10,0,0); /* UMIST	 */
	newreact("CH,SiO+=>HCO+,Si",hmrate,5.9e-10,0,0); /* UMIST	 */
	newreact("CH+,H2O=>H3O+,C",hmrate,5.8e-10,0,0); /* UMIST	 */
	newreact("CH+,H2O=>HCO+,H2",hmrate,2.9e-9,0,0); /* UMIST	 */
	newreact("CH+,O2=>HCO+,O",hmrate,9.7e-10,0,0); /* UMIST	 */
	newreact("CH+,OH=>CO+,H2",hmrate,7.5e-10,0,0); /* UMIST	 */
	newreact("CH+,O2=>CO+,OH",hmrate,1.0e-11,0,0); /* UMIST	 */
	newreact("CH2,CO+=>HCO+,CH",hmrate,4.3e-10,0,0); /* UMIST	 */
	newreact("CH2,H2O+=>H2O,CH2+",hmrate,4.7e-10,0,0); /* UMIST	 */
	newreact("CH2,H2O+=>OH,CH3+",hmrate,4.7e-10,0,0); /* UMIST	 */
	newreact("CH2,H3O+=>H2O,CH3+",hmrate,1e-10*9.4,0.,0.); /* UMIST	 */
	newreact("CH2,O2+=>O2,CH2+",hmrate,1e-10*4.3,0.,0.); /* UMIST	 */
	/* Orphaned comment which was sitting above this reaction... */
	/* >>chng 06 Apr 13, add N2H+ to chemistry.  Leading coefficient
	 * is branching ratio for this reaction, which is taken from:
	 * >>refer	mole	Geppert, W. D. et al. 2005, ApJ, 609, 459 */
	newreact("CH2,OH=>H2O,CH",hmrate,1.44e-11,0.5,3000); /* UMIST	 */
	newreact("CH2,OH+=>OH,CH2+",hmrate,1e-10*4.8,0.,0.); /* UMIST	 */
	newreact("CH2,OH+=>O,CH3+",hmrate,1e-10*4.8,0.,0.); /* UMIST	 */
	newreact("CH2+,O2=>HCO+,OH",hmrate,1e-10*9.1,0.,0.); /* UMIST	 */
	newreact("H2O,CO+=>HCO+,OH",hmrate,1e-10*8.84,0.,0.); /* UMIST	 */
	newreact("H2O+,H2O=>H3O+,OH",hmrate,1e-9*2.1,0.,0.); /* UMIST	 */
	newreact("H2O+,O2=>O2+,H2O",hmrate,1e-10*4.6,0.,0.); /* UMIST	 */
	newreact("H3O+,SiH=>SiH2+,H2O",hmrate,1e-10*9.7,0.,0.); /* UMIST	 */
	newreact("H3O+,SiO=>SiOH+,H2O",hmrate,2.0e-9,0.,0.); /* UMIST	 */
	newreact("OH,CO+=>HCO+,O",hmrate,3.1e-10,0,0); /* UMIST	 */
	newreact("OH,H2O+=>H3O+,O",hmrate,6.9e-10,0,0); /* UMIST	 */
	newreact("OH,OH=>H2O,O",hmrate,1.65e-12,1.14,50); /* UMIST	 */
	newreact("OH+,H2O=>H3O+,O",hmrate,1.3e-9,0,0); /* UMIST	 */
	newreact("OH+,H2O=>H2O+,OH",hmrate,1.59e-9,0,0); /* UMIST	 */
	newreact("OH+,O2=>O2+,OH",hmrate,5.9e-10,0,0); /* UMIST	 */
	newreact("OH+,OH=>H2O+,O",hmrate,7.0e-10,0,0); /* UMIST	 */
	newreact("OH+,SiH=>SiH2+,O",hmrate,1.0e-9,0,0); /* UMIST	 */
	newreact("OH+,SiO=>SiOH+,O",hmrate,9.4e-10,0,0); /* UMIST	 */
	newreact("C,CH5+=>CH4,CH+",hmrate,1.2e-9,0,0); /* UMIST     */
	newreact("O,CH4=>OH,CH3",hmrate,2.29e-12,2.2,3820); /* UMIST     */
	newreact("O,CH4+=>OH,CH3+",hmrate,1.0e-9,0,0); /* UMIST     */
	newreact("O,CH5+=>H3O+,CH2",hmrate,2.2e-10,0,0); /* UMIST     */
	newreact("O+,CH4=>OH,CH3+",hmrate,1.1e-10,0,0); /* UMIST     */
	newreact("O+,CH4=>CH4+,O",hmrate,8.9e-10,0,0); /* UMIST     */
	newreact("CH4,CH=>CH3,CH2",hmrate,2.28e-11,0.7,3000); /* UMIST     */
	newreact("CH5+,CH=>CH4,CH2+",hmrate,6.9e-10,0,0); /* UMIST     */
	newreact("CH2,CH2=>CH3,CH",hmrate,4.0e-10,0,5000); /* UMIST     */
	newreact("CH4,CH2=>CH3,CH3",hmrate,7.13e-12,0,5050); /* UMIST     */
	newreact("OH,CH2=>O,CH3",hmrate,1.44e-11,0.5,3000); /* UMIST     */
	newreact("CH5+,CH2=>CH4,CH3+",hmrate,9.6e-10,0,0); /* UMIST     */
	newreact("OH,CH3=>CH4,O",hmrate,3.27e-14,2.2,2240); /* UMIST     */
	newreact("OH,CH3=>H2O,CH2",hmrate,1.2e-10,0,1400); /* UMIST     */
	newreact("H2O,CH3=>OH,CH4",hmrate,2.3e-15,3.47,6681); /* UMIST     */
	newreact("CH3,CH3=>CH4,CH2",hmrate,7.13e-12,0,5052); /* UMIST     */
	newreact("OH,CH4=>H2O,CH3",hmrate,3.77e-13,2.42,1162); /* UMIST     */
	newreact("OH+,CH4=>CH5+,O",hmrate,1.95e-10,0,0); /* UMIST     */
	newreact("OH+,CH4=>H3O+,CH2",hmrate,1.31e-9,0,0); /* UMIST     */
	newreact("H2O+,CH4=>H3O+,CH3",hmrate,1.4e-9,0,0); /* UMIST     */
	newreact("CO+,CH4=>HCO+,CH3",hmrate,4.55e-10,0,0); /* UMIST     */
	newreact("CH4,CH4+=>CH5+,CH3",hmrate,1.5e-9,0,0); /* UMIST     */
	newreact("H2O,CH4+=>H3O+,CH3",hmrate,2.6e-9,0,0); /* UMIST     */
	newreact("O2,CH4+=>O2+,CH4",hmrate,4.0e-10,0,0); /* UMIST     */
	newreact("H2O,CH5+=>H3O+,CH4",hmrate,3.7e-9,0,0); /* UMIST     */
	newreact("CH5+,OH=>H2O+,CH4",hmrate,7.0e-10,0,0); /* UMIST     */
	newreact("C,NH=>N,CH",hmrate,1.73e-11,0.5,4000); /* UMIST   */
	newreact("C,NH=>CN,H",hmrate,1.1e-10,0.5,0); /* UMIST   */
	newreact("C,N2=>CN,N",hmrate,8.69e-11,0,22600); /* UMIST   */
	newreact("C,NO=>CN,O",hmrate,4.8e-11,0,0); /* UMIST   */
	newreact("C,HS=>S,CH",hmrate,1.2e-11,0.58,5880); /* UMIST   */
	newreact("C,HS=>CS,H",hmrate,2.0e-11,0,0); /* UMIST   */
	newreact("C,NS=>S,CN",hmrate,2.0e-11,0.5,0); /* UMIST   */
	newreact("C,NS=>CS,N",hmrate,1.73e-11,0.5,4000); /* UMIST   */
	newreact("C,S2=>CS,S",hmrate,1.73e-11,0.5,0); /* UMIST   */
	newreact("CH,N=>NH,C",hmrate,3.03e-11,0.65,1207); /* UMIST   */
	newreact("CH,N2=>HCN,N",hmrate,0.00000000000056,0.88,10128); /* UMIST   */
	newreact("CH,NO=>HCN,O",hmrate,5.59e-9,0,10814); /* UMIST   */
	newreact("CH,NO=>CN,OH",hmrate,2.32E-26,0,0); /* UMIST   */
	newreact("CH,NO=>OCN,H",hmrate,1.13E-25,0,0); /* UMIST   */
	newreact("CH,HNO=>NO,CH2",hmrate,1.73e-11,0.5,0); /* UMIST   */
	newreact("CH,S=>HS,C",hmrate,1.73e-11,0.5,4000); /* UMIST   */
	newreact("CH,S=>CS,H",hmrate,1.1e-12,0,0); /* UMIST   */
	newreact("N,NH=>N2,H",hmrate,4.98e-11,0,0); /* UMIST   */
	newreact("N,CH3=>HCN,H,H",hmrate,3.32e-13,0,0); /* UMIST   */
	newreact("N,CH3=>HCN,H2",hmrate,1.3e-11,0.5,0); /* UMIST   */
	newreact("N,OH=>O,NH",hmrate,1.88e-11,0.1,10700); /* UMIST   */
	newreact("N,OH=>NO,H",hmrate,5.32e-11,-0.25,0); /* UMIST   */
	newreact("N,CN=>N2,C",hmrate,3.e-10,0,0); /* UMIST   */
	newreact("N,SiH=>SiN,H",hmrate,5.e-11,0.5,0); /* UMIST   */
	newreact("N,NO=>N2,O",hmrate,3.75e-11,0,26); /* UMIST   */
	newreact("N,HNO=>NO,NH",hmrate,2.94e-12,0.5,1000); /* UMIST   */
	newreact("N,HNO=>N2O,H",hmrate,1.43e-12,0.5,1500); /* UMIST   */
	newreact("N,O2=>NO,O",hmrate,2.26e-12,0.86,3134); /* UMIST   */
	newreact("N,HS=>NS,H",hmrate,1.73e-11,0.5,0); /* UMIST   */
	newreact("N,HS=>S,NH",hmrate,1.73e-11,0.5,9060); /* UMIST   */
	newreact("N,CS=>S,CN",hmrate,3.8e-11,0.5,1160); /* UMIST   */
	newreact("N,NS=>S,N2",hmrate,1.73e-11,0.5,0); /* UMIST   */
	newreact("N,SO=>NS,O",hmrate,4.68e-11,0.5,8254); /* UMIST   */
	newreact("N,SO=>S,NO",hmrate,1.73e-11,0.5,750); /* UMIST   */
	newreact("N,S2=>NS,S",hmrate,1.73e-11,0.5,4000); /* UMIST   */
	newreact("CH2,CN=>HCN,CH",hmrate,5.3e-12,0,2500); /* UMIST   */
	newreact("CH2,NO=>HCN,OH",hmrate,0.000000000000832,0,1443); /* UMIST   */
	newreact("CH2,HNO=>NO,CH3",hmrate,1.73e-11,0.5,0); /* UMIST   */
	newreact("NH,NH=>N2,H,H",hmrate,1.16e-9,0,0); /* UMIST   */
	newreact("NH,NH=>N2,H2",hmrate,1.7e-11,0,0); /* UMIST   */
	newreact("NH,O=>OH,N",hmrate,1.16e-11,0,0); /* UMIST   */
	newreact("NH,O=>NO,H",hmrate,1.16e-10,0,0); /* UMIST   */
	newreact("NH,OH=>HNO,H",hmrate,3.32e-11,0,0); /* UMIST   */
	newreact("NH,OH=>NH2,O",hmrate,2.93e-12,0.1,5800); /* UMIST   */
	newreact("NH,OH=>H2O,N",hmrate,3.11e-12,1.2,0); /* UMIST   */
	newreact("NH,H2O=>OH,NH2",hmrate,1.83e-12,1.6,14090); /* UMIST   */
	newreact("NH,CN=>HCN,N",hmrate,2.94e-12,0.5,1000); /* UMIST   */
	newreact("NH,NO=>N2O,H",hmrate,1.16e-10,-1.03,420); /* UMIST   */
	newreact("NH,NO=>N2,O,H",hmrate,5.e-11,0,0); /* UMIST   */
	newreact("NH,NO=>N2,OH",hmrate,1.46e-11,-0.58,37); /* UMIST   */
	newreact("NH,S=>HS,N",hmrate,1.73e-11,0.5,4000); /* UMIST   */
	newreact("NH,S=>NS,H",hmrate,1.73e-11,0.5,0); /* UMIST   */
	newreact("NH,NO2=>HNO,NO",hmrate,5.72e-12,0.5,2500); /* UMIST   */
	newreact("NH,NO2=>N2O,OH",hmrate,1.44e-13,0, 1140); /* UMIST   */
	newreact("CH3,NH3=>CH4,NH2",hmrate,0.0000000000000955,0,4890); /* UMIST   */
	newreact("CH3,CN=>HCN,CH2",hmrate,9.21e-12,0.7,1500); /* UMIST   */
	newreact("CH3,HNO=>NO,CH4",hmrate,1.44e-11,0.5,0); /* UMIST   */
	newreact("O,NH2=>OH,NH",hmrate,1.39e-11,0,40); /* UMIST   */
	newreact("O,NH2=>NO,H2",hmrate,8.3e-12,0,0); /* UMIST   */
	newreact("O,NH3=>OH,NH2",hmrate,1.89e-11,0,4003); /* UMIST   */
	newreact("O,CN=>NO,C",hmrate,3.81e-11,0.5,14545); /* UMIST   */
	newreact("O,HCN=>CN,OH",hmrate,6.21e-10,0,12439); /* UMIST   */
	newreact("O,HCN=>OCN,H",hmrate,1.36e-12,1.38,3693); /* UMIST   */
	newreact("O,N2=>NO,N",hmrate,2.51e-10,0,38602); /* UMIST   */
	newreact("O,NO=>O2,N",hmrate,1.18e-11,0,20413); /* UMIST   */
	newreact("O,HNO=>NO,OH",hmrate,6.e-11,0,0); /* UMIST   */
	newreact("O,HNO=>O2,NH",hmrate,2.94e-12,0.5,3500); /* UMIST   */
	newreact("O,HNO=>NO2,H",hmrate,1.44e-12,0.5,0); /* UMIST   */
	newreact("O,HS=>S,OH",hmrate,1.74e-11,0.67,956); /* UMIST   */
	newreact("O,HS=>SO,H",hmrate,2.32e-10,0,0); /* UMIST   */
	newreact("O,OCN=>O2,CN",hmrate,4.02e-10,-1.43,3501); /* UMIST   */
	newreact("O,SiN=>NO,Si",hmrate,2.5e-11,0.5,0); /* UMIST   */
	newreact("O,SiN=>SiO,N",hmrate,5.75e-11,0.1,200); /* UMIST   */
	newreact("O,N2O=>NO,NO",hmrate,1.15e-10,0,13400); /* UMIST   */
	newreact("O,N2O=>O2,N2",hmrate,1.66e-10,0,14100); /* UMIST   */
	newreact("O,CS=>SO,C",hmrate,4.68e-11,0.5,28940); /* UMIST   */
	newreact("O,NS=>S,NO",hmrate,5.e-11,0.5,0); /* UMIST   */
	newreact("O,NS=>SO,N",hmrate,1.73e-11,0.5,4000); /* UMIST   */
	newreact("O,SO=>S,O2",hmrate,0.00000000000066,0,2760); /* UMIST   */
	newreact("O,S2=>SO,S",hmrate,1.73e-11,0.5,0); /* UMIST   */
	newreact("NH2,OH=>NH3,O",hmrate,0.000000000000208,0.76,262); /* UMIST   */
	newreact("NH2,OH=>H2O,NH",hmrate,1.5e-12,0,0); /* UMIST   */
	newreact("NH2,NO=>N2,OH,H",hmrate,1.49e-12,0,0); /* UMIST   */
	newreact("NH2,NO=>N2,H2O",hmrate,4.27e-11,-2.5,331); /* UMIST   */
	newreact("CH4,CN=>HCN,CH3",hmrate,3.14e-12,1.53,504); /* UMIST   */
	newreact("OH,NH3=>H2O,NH2",hmrate,1.47e-13,2.05,7); /* UMIST   */
	newreact("OH,CN=>HCN,O",hmrate,1.0e-11,0,1000); /* UMIST   */
	newreact("OH,CN=>OCN,H",hmrate,7.e-11,0,0); /* UMIST   */
	newreact("OH,HCN=>CN,H2O",hmrate,1.87e-13,1.5,3887); /* UMIST   */
	newreact("OH,NO=>NO2,H",hmrate,5.2e-12,0,15100); /* UMIST   */
	newreact("OH,S=>HS,O",hmrate,6.6e-11,0,0); /* UMIST   */
	newreact("OH,S=>SO,H",hmrate,6.6e-11,0,0); /* UMIST   */
	newreact("OH,N2O=>HNO,NO",hmrate,1.04E-17,4.33,12623); /* UMIST   */
	newreact("OH,CS=>OCS,H",hmrate,0.0000000000000939,1.12,800); /* UMIST   */
	newreact("NH3,CN=>HCN,NH2",hmrate,2.6e-11,-1.1,0); /* UMIST   */
	newreact("CN,NO=>OCN,N",hmrate,1.62e-10,0,21205); /* UMIST   */
	newreact("CN,HNO=>NO,HCN",hmrate,1.5e-11,0.5,0); /* UMIST   */
	newreact("CN,S=>NS,C",hmrate,5.71e-11,0.5,32010); /* UMIST   */
	newreact("CN,S=>CS,N",hmrate,1.73e-11,0.5,0); /* UMIST   */
	newreact("N2,O2=>N2O,O",hmrate,1.0e-10,0,55200); /* UMIST   */
	newreact("NO,NO=>O2,N2",hmrate,2.51e-11,0,30653); /* UMIST   */
	newreact("NO,NO=>N2O,O",hmrate,7.22e-12,0,33155); /* UMIST   */
	newreact("NO,HNO=>N2O,OH",hmrate,1.41e-11,0,14890); /* UMIST   */
	newreact("NO,S=>NS,O",hmrate,2.94e-11,0.5,17465); /* UMIST   */
	newreact("NO,S=>SO,N",hmrate,1.75e-10,0,20200); /* UMIST   */
	newreact("O2,S=>SO,O",hmrate,2.28e-12,0.52,0); /* UMIST   */
	newreact("S,SO=>S2,O",hmrate,1.73e-11,0.5,11500); /* UMIST   */
	newreact("C,NH+=>N,CH+",hmrate,1.6e-9,0,0); /* UMIST   */
	newreact("C+,NH=>CH+,N",hmrate,5e-10,0.,0.); /* Zsargo & Federman 2003   */
	newreact("C+,NH2=>HCN+,H",hmrate,1.1e-9,0,0); /* UMIST   */
	newreact("C,NH2+=>NH,CH+",hmrate,1.2e-9,0,0); /* UMIST   */
	newreact("C+,NH3=>HCN+,H2",hmrate,7.e-11,0,0); /* UMIST   */
	newreact("C,NH3+=>NH,CH2+",hmrate,1.0e-11,0,0); /* UMIST   */
	newreact("C,HCN+=>CN,CH+",hmrate,1.1e-9,0,0); /* UMIST   */
	newreact("C,HNO+=>NO,CH+",hmrate,1.0e-9,0,0); /* UMIST   */
	newreact("C+,HS=>CS+,H",hmrate,1.1e-9,0,0); /* UMIST   */
	newreact("C,HS+=>CS+,H",hmrate,9.9e-10,0,0); /* UMIST   */
	newreact("C+,OCN=>CO+,CN",hmrate,3.8e-9,0,0); /* UMIST   */
	newreact("C+,NS=>CS+,N",hmrate,7.6e-10,0,0); /* UMIST   */
	newreact("C+,SO=>S,CO+",hmrate,2.6e-10,0,0); /* UMIST   */
	newreact("C+,SO=>CS+,O",hmrate,2.6e-10,0,0); /* UMIST   */
	newreact("CH+,N=>CN+,H",hmrate,1.9e-10,0,0); /* UMIST   */
	newreact("CH,N+=>CN+,H",hmrate,3.6e-10,0,0); /* UMIST   */
	newreact("CH,NH+=>CH2+,N",hmrate,9.9e-10,0,0); /* UMIST   */
	newreact("CH+,NH=>CN+,H2",hmrate,7.6e-10,0,0); /* UMIST   */
	newreact("CH+,NH2=>HCN+,H2",hmrate,1.1e-9,0,0); /* UMIST   */
	newreact("CH,NH2+=>NH,CH2+",hmrate,3.5e-10,0,0); /* UMIST   */
	newreact("CH+,NH3=>NH4+,C",hmrate,4.05e-10,0,0); /* UMIST   */
	newreact("CH,NH3+=>NH4+,C",hmrate,6.9e-10,0,0); /* UMIST   */
	newreact("CH,HCN+=>CN,CH2+",hmrate,6.3e-10,0,0); /* UMIST   */
	newreact("CH,HNO+=>NO,CH2+",hmrate,6.2e-10,0,0); /* UMIST   */
	newreact("CH+,S=>HS+,C",hmrate,4.7e-10,0,0); /* UMIST   */
	newreact("CH+,S=>CS+,H",hmrate,4.7e-10,0,0); /* UMIST   */
	newreact("CH,S+=>CS+,H",hmrate,6.2e-10,0,0); /* UMIST   */
	newreact("CH,HS+=>S,CH2+",hmrate,5.8e-10,0,0); /* UMIST   */
	newreact("N,CH2+=>HCN+,H",hmrate,2.2e-10,0,0); /* UMIST   */
	newreact("N+,NH=>N2+,H",hmrate,3.7e-10,0,0); /* UMIST   */
	newreact("N,NH+=>N2+,H",hmrate,1.3e-9,0,0); /* UMIST   */
	newreact("N+,CH4=>HCN+,H2,H",hmrate,5.6e-11,0,0); /* UMIST   */
	newreact("N,OH+=>NO+,H",hmrate,8.9e-10,0,0); /* UMIST   */
	newreact("N+,NH3=>NH2+,NH",hmrate,2.16e-10,0,0); /* UMIST   */
	newreact("N,H2O+=>HNO+,H",hmrate,1.9e-10,0,0); /* UMIST   */
	newreact("N+,NO=>N2+,O",hmrate,7.9e-11,0,0); /* UMIST   */
	newreact("N+,O2=>NO,O+",hmrate,3.66e-11,0,0); /* UMIST   */
	newreact("N+,O2=>NO+,O",hmrate,2.63e-10,0,0); /* UMIST   */
	newreact("N,O2+=>NO+,O",hmrate,1.8e-10,0,0); /* UMIST   */
	newreact("N,HS+=>NS+,H",hmrate,7.4e-10,0,0); /* UMIST   */
	newreact("N,SiO+=>NO+,Si",hmrate,9.e-11,0,0); /* UMIST   */
	newreact("N,SiO+=>NO,Si+",hmrate,2.1e-10,0,0); /* UMIST   */
	newreact("N,SO+=>NS+,O",hmrate,5.e-11,0,0); /* UMIST   */
	newreact("N+,OCS=>CS+,NO",hmrate,7.e-11,0,0); /* UMIST   */
	newreact("CH2,NH+=>CH3+,N",hmrate,1.4e-9,0,0); /* UMIST   */
	newreact("CH2,NH2+=>CH3+,NH",hmrate,4.9e-10,0,0); /* UMIST   */
	newreact("CH2+,NH3=>NH4+,CH",hmrate,1.26e-9,0,0); /* UMIST   */
	newreact("CH2,NH3+=>NH2,CH3+",hmrate,9.6e-10,0,0); /* UMIST   */
	newreact("CH2,HCN+=>CN,CH3+",hmrate,8.7e-10,0,0); /* UMIST   */
	newreact("CH2,HNO+=>NO,CH3+",hmrate,8.6e-10,0,0); /* UMIST   */
	newreact("CH2+,S=>HCS+,H",hmrate,1.4e-9,0,0); /* UMIST   */
	newreact("CH2,S+=>HCS+,H",hmrate,1.0e-11,0,0); /* UMIST   */
	newreact("NH+,NH=>NH2+,N",hmrate,1.0e-9,0,0); /* UMIST   */
	newreact("NH+,O=>OH+,N",hmrate,1.0e-9,0,0); /* UMIST   */
	newreact("NH,O+=>NO+,H",hmrate,3.6e-10,0,0); /* UMIST   */
	newreact("NH+,NH2=>NH3+,N",hmrate,1.5e-9,0,0); /* UMIST   */
	newreact("NH,NH2+=>NH3+,N",hmrate,7.3e-10,0,0); /* UMIST   */
	newreact("NH+,OH=>H2O+,N",hmrate,1.0e-9,0,0); /* UMIST   */
	newreact("NH,OH+=>NH2+,O",hmrate,3.6e-10,0,0); /* UMIST   */
	newreact("NH+,NH3=>NH4+,N",hmrate,6.e-10,0,0); /* UMIST   */
	newreact("NH,NH3+=>NH4+,N",hmrate,7.1e-10,0,0); /* UMIST   */
	newreact("NH,CH5+=>CH4,NH2+",hmrate,7.1e-10,0,0); /* UMIST   */
	newreact("NH+,H2O=>NH3+,O",hmrate,1.75e-10,0,0); /* UMIST   */
	newreact("NH+,H2O=>H3O+,N",hmrate,1.05e-9,0,0); /* UMIST   */
	newreact("NH+,H2O=>HNO+,H2",hmrate,3.5e-10,0,0); /* UMIST   */
	newreact("NH,H2O+=>H3O+,N",hmrate,7.1e-10,0,0); /* UMIST   */
	newreact("NH+,H2O=>OH,NH2+",hmrate,8.75e-10,0,0); /* UMIST   */
	newreact("NH+,CN=>HCN+,N",hmrate,1.6e-9,0,0); /* UMIST   */
	newreact("NH,HCN+=>CN,NH2+",hmrate,6.5e-10,0,0); /* UMIST   */
	newreact("NH,CO+=>HCO+,N",hmrate,3.2e-10,0,0); /* UMIST   */
	newreact("NH,Si+=>SiN+,H",hmrate,1.0e-9,0,0); /* UMIST   */
	newreact("NH,HNO+=>NO,NH2+",hmrate,6.3e-10,0,0); /* UMIST   */
	newreact("NH,O2+=>HNO+,O",hmrate,3.2e-10,0,0); /* UMIST   */
	newreact("NH+,O2=>NO+,OH",hmrate,2.05e-10,0,0); /* UMIST   */
	newreact("NH,O2+=>NO2+,H",hmrate,3.2e-10,0,0); /* UMIST   */
	newreact("NH+,S=>HS+,N",hmrate,6.9e-10,0,0); /* UMIST   */
	newreact("NH+,S=>NS+,H",hmrate,6.9e-10,0,0); /* UMIST   */
	newreact("NH,S+=>NS+,H",hmrate,6.3e-10,0,0); /* UMIST   */
	newreact("CH3+,NH3=>NH4+,CH2",hmrate,3.4e-10,0,0); /* UMIST   */
	newreact("CH3+,S=>HCS+,H2",hmrate,1.4e-9,0,0); /* UMIST   */
	newreact("O,NH2+=>HNO+,H",hmrate,7.2e-11,0,0); /* UMIST   */
	newreact("O,NH3+=>HNO+,H2",hmrate,1.0e-11,0,0); /* UMIST   */
	newreact("O+,CN=>NO+,C",hmrate,1.0e-9,0,0); /* UMIST   */
	newreact("O+,HCN=>CO+,NH",hmrate,1.2e-9,0,0); /* UMIST   */
	newreact("O+,HCN=>NO+,CH",hmrate,1.2e-9,0,0); /* UMIST   */
	newreact("O+,HCN=>HCO+,N",hmrate,1.2e-9,0,0); /* UMIST   */
	newreact("O+,N2=>NO+,N",hmrate,1.2e-12,0,0); /* UMIST   */
	newreact("O,N2+=>NO+,N",hmrate,1.3e-10,0,0); /* UMIST   */
	newreact("O,HNO+=>NO2+,H",hmrate,1.e-12,0,0); /* UMIST   */
	newreact("O,HS+=>S+,OH",hmrate,2.9e-10,0,0); /* UMIST   */
	newreact("O,HS+=>SO+,H",hmrate,2.9e-10,0,0); /* UMIST   */
	newreact("O,SiN+=>SiO+,N",hmrate,1.0e-9,0,0); /* UMIST   */
	newreact("O+,N2O=>NO+,NO",hmrate,6.3e-10,0,0); /* UMIST   */
	newreact("O,CS+=>S,CO+",hmrate,6.e-11,0,0); /* UMIST   */
	newreact("O,HCS+=>S,HCO+",hmrate,5.e-12,0,0); /* UMIST   */
	newreact("O,HCS+=>OCS+,H",hmrate,5.e-12,0,0); /* UMIST   */
	newreact("O+,NO2=>O2,NO+",hmrate,8.3e-10,0,0); /* UMIST   */
	newreact("O,NS+=>S,NO+",hmrate,6.1e-10,0,0); /* UMIST   */

	/* Create linear list of species and populate it... */
	coreactions.list =  (struct COmole_rate_s **)MALLOC((size_t)coreactions.n*
																		 sizeof(struct COmole_rate_s *));

	/* ...first active species */
	i = makeplist(mole_priv.reactab,(void **)coreactions.list,
								coreactions.n,NULL); 
	ASSERT (i == coreactions.n); 

}
STATIC void newreact(const char label[], 
										 double (*fun)(struct COmole_rate_s *rate), double a, double b, double c)
{
	struct COmole_rate_s *rate;
	struct molecule    *sp;
	data_u *p;
	int i,j,prod,exists;
	char buf[7];

	DEBUG_ENTRY("newreact()");

	coreactions.n++;
	rate = (struct COmole_rate_s *) MALLOC (sizeof(struct COmole_rate_s));
	p = addentry(label,0,mole_priv.reactab,&exists);
	p->p = (void *) rate;
	rate->label = (char *) p->key;
	rate->fun = fun;
	if(fun == hmrate && b == 0. && c == 0.) /* Reaction rate is actually simpler */
	{
		rate->fun = constrate;
	}
	rate->a = a;
	rate->b = b;
	rate->c = c;
	rate->rk = 0.0; /* Sane initial value, in case used early */

	rate->index = coreactions.n-1;

	rate->nreactants = rate->nrates = rate->nproducts = rate->photon = 0;
	j = prod = 0;
	for(i=0;!i || label[i-1]!='\0';i++) 
	{
		if(label[i] == ',' || label[i] == '=' || label[i] == '\0') 
		{
			buf[j] = '\0';
			j = 0;
			sp = findspecies(buf);
			if(sp != &null_mole) 
			{
				if(prod == 0) 
				{
					fixit(); /* Bodge, should really test if active when filling matrix not when defining reaction */
					if(sp->active) 
					{
						if(rate->nreactants >= MAXREACTANTS) 
						{
							fprintf(stderr,"Mole_co_etc: Too many reactants in %s, only %d allowed\n",label,MAXREACTANTS);
							cdEXIT( EXIT_FAILURE );
						}
						rate->reactants[rate->nreactants] = sp;
						rate->nreactants++;
					}
					if(rate->nrates >= MAXREACTANTS) 
					{
						fprintf(stderr,"Mole_co_etc: Too many rate species in %s, only %d allowed\n",label,MAXREACTANTS);
						cdEXIT( EXIT_FAILURE );
					}
					rate->rate_species[rate->nrates] = sp;
					rate->nrates++;
				} 
				else 
				{
					fixit(); /* Bodge, should really test if active when filling matrix not when defining reaction */
					if(sp->active) 
					{
						if(rate->nproducts >= MAXPRODUCTS) 
						{
						  fprintf(stderr,"Mole_co_etc: Too many products in %s, only %d allowed\n",label,MAXPRODUCTS);
							cdEXIT( EXIT_FAILURE );
						}
						rate->products[rate->nproducts] = sp;
						rate->nproducts++;
					}
				}
			}
			else 
			{
				if(0) 
					fprintf(stderr,"Could not find %s\n",buf);
				if(strncmp(buf,"PHOTON",6) == 0) {
					if(prod == 0)
						rate->photon--;
					else
						rate->photon++;
				}
				fixit(); /* Should do something proper about non-network species */
			}
			if(label[i] == '=') 
			{
				i++;
				prod = 1;
			}
		} 
		else 
		{
			buf[j] = label[i];
			j++;
		}
	}

	/* >> chng 06 Oct 10 rjrw: use 1/(1/m1+1/m2) for reduced mass to prevent underflow */
	if(rate->nrates == 2)
	{
		rate->reduced_mass = 1./(1./rate->rate_species[0]->mole_mass+1./rate->rate_species[1]->mole_mass);
	}
	else
	{
		rate->reduced_mass = 0.;
	}
}


/*
 * Functions to specify chemical rates -- note that the rate->a overall scale
 * factor is applied in CO_update_rks
 *
 */

#include "phycon.h"
#include "physconst.h"
#include "doppvel.h"

STATIC double noneq_offset(struct COmole_rate_s *rate);

STATIC double hmrate(struct COmole_rate_s *rate) 
{
	double te;

	DEBUG_ENTRY("hmrate()");

	te = phycon.te+noneq_offset(rate);

	return pow(te/300.,rate->b)*exp(-rate->c/te);
}

/* Add in non-equilibrium chemistry.  This is done by assuming
 * that turbulence reduces the temperature barrier, thereby
 * enhancing reaction rates for molecules such as CH+.  The 
 * "effective temperature is defined according to 
 * >>refer Federman et al. 1996, MNRAS. L41-46  */

/* The effective temperature is defined as:
 * T(effective) = T + (1/3)*reduced_mass*turbulent_velocity^2/BOLTZMANN_CONSTANT
 */
STATIC double noneq_offset(struct COmole_rate_s *rate)
{	
  /* This logic could be cached by using distinct rate functions in newreact */
	int nreact, n;
	bool lgFact;

	DEBUG_ENTRY("noneq_offset()");

	lgFact = false;
	if(co.lgNonEquilChem)
	{ 
		if(co.lgNeutrals) 
		{
			lgFact = true;
		}
		else
		{
			nreact = rate->nreactants;
			for(n=0;n<nreact;n++) 
			{
				if(rate->reactants[0]->nElec != 0)
				{
					lgFact = true;
					break;
				}
			}
		}
	}

	if( lgFact ) 
		return 0.333f*POW2(DoppVel.TurbVel)/BOLTZMANN*rate->reduced_mass;
	else
		return 0.;
}
STATIC double constrate(struct COmole_rate_s *) 
{
	return 1.;
}
	/* hmi.UV_Cont_rel2_Habing_TH85_depth is field relative to Habing background, dimensionless */
	/* >>chng 04 apr 01, move from TH85 to DB96, also correct leading coef to use
	 * UMIST database value */
	/* CO_photo_dissoc_rate = 1.67e-10f*hmi.UV_Cont_rel2_Habing_TH85_depth;*/

	/* TRY MOVING PHOTORATES OUT OF LOOP */

	/* >>chng 02 jul 04 -- The following are dissociation rates for various molecular species
	For right now we will calculate this rate by the standard form:
		(alpha)*Go*exp(-Beta*AV)
	when the command "set Leiden hack UMIST rates" is used.  Otherwise we
	will just let cloudy calculate the value of the UV radiation field */
#include "rfield.h"
STATIC double th85rate(struct COmole_rate_s *rate) 
{
	double rk;

	DEBUG_ENTRY("th85rate()");

	if(co.lgUMISTrates || rate->b == 0.0)
	{
		rk = hmi.UV_Cont_rel2_Habing_TH85_depth/1.66;
	}
	else 
	{
		rk = hmi.UV_Cont_rel2_Habing_TH85_face/1.66*exp(-(rate->b*rfield.extin_mag_V_point));
	}

	return rk;
}
#include "secondaries.h"
	/* >> chng aug 24, 05 NPA This is the cosmic ray ionization rate used in the molecular network.  
	 * TH85 and the e-mail from Amiel Sternberg has each cosmic ray ionization rate as a 
	 * leading coefficient multiplied by a scale factor.
	 * The leading coefficient is defined as the cosmic ray rate for H + CR = > H+
	 * + e- .  For molecules in the heavy element molecular network, this scale
	 * factor is derived by taking the rate for:

	X + CRPHOT => Y + Z 
	and dividing it by the rate:
	H + CRP => H+ + e-

	This scale factor is 2.17 for all cosmic ray reactions in this network 
		crnu_rate = secondaries.csupra[ipHYDROGEN][0];*/
STATIC double crnurate(struct COmole_rate_s *) 
{
	return 2.17*secondaries.csupra[ipHYDROGEN][0];
}
#include "ionbal.h"
STATIC double co_lnu_c_o_lnu(struct COmole_rate_s *)
{
	double val = 0;
	int ns, ion;
	/* inner shell photoionization of CO, assume rates are same as K 1s and 2s
	 * shell of C and O */
	/* >>chng 04 may 26, upper limit should be ns<2, had been ns<2 so picked up
	 * valence shell, which was incorrect */

	DEBUG_ENTRY("co_lnu_c_o_lnu()");

	for( ns=0; ns<2; ++ns )
	{
		ion = 0;
		val += ionbal.PhotoRate_Shell[ipCARBON][ion][ns][0];
		val += ionbal.PhotoRate_Shell[ipOXYGEN][ion][ns][0];
	}

	return val;
}

STATIC double ele_ion_ladder(struct COmole_rate_s *rate)
{
	long int ipElem;

	ipElem = rate->reactants[0]->nelem_hevmol;

	if(rate->reactants[0]->nElec == 1)
		return ionbal.RateRecomTot[ipElem][0]+gv.GrainChTrRate[ipElem][1][0];
	else
		return ionbal.RateIonizTot(ipElem,0)+gv.GrainChTrRate[ipElem][0][1];
}
	/******************************** Gas-Grain Chemistry**********************************/

	/*  The Gas-Grain Chemistry rates are taken from: 
	>>refer	Hasegawa, T. I. & Herbst, E. 1993, MNRAS, 261, 83

	So far only CO depletion onto grains is considered, however, this code can be generalized 
	if desired to other molecules, using the data in the above paper.  There are three important reactions 
	to determine the abundance of a molecule on grain surfaces deep in molecular clouds.  The
	rate of accretion of molecule X onto grains is

	R(accretion) = PI*[grain_radius]^2*[thermal velocity]*[density_of_grains]

	Two processes remove molecule X from the grain surface, and that is thermal evaporation, due
	to the heat of the grain, and cosmic ray deabsorption.  The first of these rates come from the 
	above paper, and depends primarily on the dust temperature.  The cosmic ray rate is a constant,
	calculated in Hasegawa and Herbst.  

	For each molecule desired, I have created a new species which is the density of that molecule
	on the grain surface */


	/* evaporation rate of molecule on grain is:

	k(evap) = [vibrational absorption frequency]*exp[-binding_energy/dust_temperature]

	The binding energies come from Hasegawa and Herbst, Table 4.  The vibrational frequency comes from
	equation 3 of >>refer Hasegawa, T. I., Herbst, E., & Leung, C. M. 1992, ApJSS, 82, 167

	[vibrational absorption frequency] = 
	SQRT[[2*number_of_sites_for_grain_absorption*binding_energy]/[PI^2*mass_of_molecule]]

	**********************************************************************************************/

STATIC double vib_evap(struct COmole_rate_s *rate)
{
	double binding_energy,  exponent, vib_freq, number_of_sites /* on grain */;

	DEBUG_ENTRY("vib_evap()");

	exponent = 0.0;

	binding_energy = rate->b;
	/*>>chng 06 nov 28 only include source from molecules if we have an estimated first 
	 * solution - first test is that we have called mole at least twice,
	 * second test is that we are on a later iteration */
	if( conv.nTotalIoniz > 1 || iteration > 1 )
	{
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			exponent += exp(-binding_energy/gv.bin[nd]->tedust);
		}	
	}
	number_of_sites = 1.5e15;

	vib_freq = sqrt(2*number_of_sites*binding_energy/(PI*PI*rate->reactants[0]->mole_mass));

	/*>>chng 06 jan 11, NPA - In some H+ regions the grain temperatures are so low
	that molecular freeze out occurs.  This should not happen, because the ices
	should sublimate in such a harsh environment.  Therefore, we introduce an
	artificial sublimation rate to destroy ices.  THIS IS NOT A PHYSICAL RATE!!!!
	only a rate that gets the desired, realistic result */
	/*>>chng 06 sep 03 rjrw -- include this in standard evaporation rate coeff (the artificial part is the sexp term) */
	/** \todo	0	find physical theory for this process */
	/* Rate comes from Table curve and assumes that rate is high (~1) in H+
	 * region and negligible ( < 1e-30) in molecular cloud - designed to stop
	 * freeze-out above 100 K */

	return vib_freq*exponent+sexp( 555.89/phycon.sqrte - 5.55 );
}
STATIC double grn_abs(struct COmole_rate_s *rate)
{
	double den_times_area;

	DEBUG_ENTRY("grn_abs()");

	den_times_area = 0.0;

	fixit(); /* Should cache value */
	/* calculate the rates that are dependent on grain physics.  This includes grain density, 
	cross sectional area, and dust temperature of each constituent.  Note that 

	gv.bin[nd]->IntArea/4.*gv.bin[nd]->cnv_H_pCM3

	is the integrated projected grain surface area per cm^3 of gas for each grain size bin */

	/* >>chng 06 feb 28, turn off this rate when no grain molecules */
	/* >>chng 06 dec 05 rjrw: do this in newreact rather than rate */
	for( size_t nd=0; nd < gv.bin.size(); nd++ )
	{
		/* >>chng 06 mar 04, update expression for projected grain surface area, PvH */
		den_times_area += gv.bin[nd]->IntArea/4.*gv.bin[nd]->cnv_H_pCM3;
	}

	return den_times_area*sqrt(8.*BOLTZMANN*phycon.te/(PI*rate->reactants[0]->mole_mass));
}

#include "rt.h"
STATIC double th85rate_co(struct COmole_rate_s *rate)
{
	double esc_co;
	/******************************************************************************************
	*	   First define the rate, which is of the form:
	*	   
	*		R = (Ro)*(Go*exp(-3.2Av))*Beta(tau(CO))
	*
	*	   where:
	*
	*	   Ro = 1.67*e-10
	*	   (Go*exp(-3.2Av)) = hmi.UV_Cont_rel2_Habing_TH85_depth
	*	   tauCO = 4.4e-15 * findspecies("CO")->hevcol / (DopplerWidth/1e5) /
	*		(1. + phycon.sqrte*0.6019); 
	*       tauC = 1.6*e17*N(C)
	*	   Beta(tau(CO)) = esca0k2(esc_co) 
	********************************************************************************************/   
	/* eqn 12 of 
	 * >>refer	CO	dissoc	Hollenbach, D.J., Takahashi, T., & Tielens, A. 1991, ApJ, 377, 192
	 * based on
	 * >>refer	CO	dissoc	Black, J.H., & van Dishoeck, E.F. 1988, ApJ, 334, 771 */
	esc_co = 4.4e-15 * rate->reactants[0]->hevcol / 
			 /* the line width in km/s */
			 ( GetDopplerWidth( rate->reactants[0]->mole_mass/(realnum)ATOMIC_MASS_UNIT )/1e5) /
			 /* this term accounts for populations within ground elec state */
			 (1. + phycon.sqrte*0.6019);
	return esca0k2(esc_co)*th85rate(rate);
}
STATIC double oh_c2h2_co_ch3(struct COmole_rate_s *rate)
{
	/* This rate will blow up if the temperature gets too low, therefore 
	 * set this rate for T < 500 equal to the rate at 500 K */
	if(phycon.te > 500)
	{		
		return hmrate(rate);
	}
	else
	{
		return 6.3E-18;
	}
}
STATIC double h_hnc_hcn_h(struct COmole_rate_s *rate)
{
	if(phycon.te > 100)
	{		
		return hmrate(rate);
	}
	else
	{
		return 1e-15;
	}
}
