/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
* others.  For conditions of distribution and use see copyright notice in license.txt */
/*H2_CollidRateRead read collision rates */
/*H2_CollidRateEvalAll - set H2 collision rates */
#include "cddefines.h" 
#include "phycon.h"
#include "dense.h"
#include "taulines.h"
#include "input.h"
#include "h2.h"
#include "h2_priv.h"
#include "mole.h"

/* set true to print all collision rates then quit */
#define PRT_COLL	false

/* following are related to H2 - He collision data set 
* H2_He_coll_init, H2_He_coll */
#define N_H2_HE_FIT_PAR 8
static realnum ***H2_He_coll_fit_par;
static bool **lgDefn_H2He_coll;
/* following are related to H2 - H2 collision data set 
 * H2_ORH2_coll_init, H2_ORH2_coll */
/*>>chng 08 feb 27, GS*/
#define N_H2_ORH2_FIT_PAR 8
static realnum ***H2_ORH2_coll_fit_par;
static bool **lgDefn_H2ORH2_coll;

/* following are related to H2 - H2 collision data set 
 * H2_H2para_coll_init, H2_H2para_coll */
/*>>chng 08 feb 27, GS*/
#define N_H2_PAH2_FIT_PAR 8
static realnum ***H2_PAH2_coll_fit_par;
static bool **lgDefn_H2PAH2_coll;

/* compute rate coefficient for a single quenching collision */
STATIC realnum H2_CollidRateEvalOne( 
	/*returns collision rate coefficient, cm-3 s-1 for quenching collision
	 * from this upper state */
	 long iVibHi, long iRotHi,long iVibLo,
	 /* to this lower state */
	long iRotLo, long ipHi , long ipLo , 
	/* colliders are H, He, H2(ortho), H2(para), and H+ */
	long  nColl )
{
	double fitted;
	realnum rate;
	double t3Plus1 = phycon.te/1000. + 1.;
	double t3Plus1Squared = POW2(t3Plus1);
	/* these are fits to the existing collision data 
	 * used to create g-bar rates */
	double gbarcoll[N_X_COLLIDER][3] = 
	{
		{-9.9265 , -0.1048 , 0.456  },
		{-8.281  , -0.1303 , 0.4931 },
		{-10.0357, -0.0243 , 0.67   },
		{-8.6213 , -0.1004 , 0.5291 },
		{-9.2719 , -0.0001 , 1.0391 }
	};

	DEBUG_ENTRY( "H2_CollidRateEvalOne()" );

	/* first do special cases 
	 * ORNL He collision data */
	if( nColl == 1 && mole.lgH2_He_ORNL )
	{
		/* ORNL in use */

		/* H2 - He collisions 
		* The H2 - He collisional data set is controlled by the 
		* atom H2 He collisions command
		* either mole.lgH2_He_Meudon or mole.lgH2_He_ORNL must be true
		* (this is asserted).  Meudon is collider 1 and Stancil is 5. */
		/*>>chng 07 apr 04, by GS - bugfix, had used 1 for collider for g-bar */
		if( (fitted=H2_He_coll(ipHi , ipLo , phycon.te ))>0. )
		{
			/* >>chng 05 oct 02, add this collider 
			* CHECK - this is deexcitation - is that correct 
			* this is new H2 - He collision data - use it if positive 
			* not all states have collision data */
			rate = (realnum)fitted*mole.lgColl_deexec_Calc;
			/*if( PRT_COLL )
				fprintf(ioQQQ,"col fit\t%li\t%li\t%.4e\t%li\t%li\t%li\t%li\t%li\t%.4e\t%.4e\n",
				ipLo,ipHi,phycon.te,nColl,
				iVibHi,iRotHi,iVibLo,iRotLo,
				energy_wn[0][iVibHi][iRotHi] - energy_wn[0][iVibLo][iRotLo],
				H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][nColl] );*/
		}

		/* use g-bar g bar guess of rate coefficient for 
		 * collision with He, this does not change ortho & para 
		 * turn mole.lgColl_gbar on/off with atom h2 gbar on off */
		else if( mole.lgColl_gbar  && 
			(H2_lgOrtho[0][iVibHi][iRotHi]-H2_lgOrtho[0][iVibLo][iRotLo]==0) )
		{
			/* the fit is log(K)=y_0+a*((x)^b), where K is the rate coefficient,
			* and x is the energy in wavenumbers */
			double ediff = energy_wn[0][iVibHi][iRotHi] - energy_wn[0][iVibLo][iRotLo];

			/* do not let energy difference be smaller than 100 wn, the smallest
			* difference for which we fit the rate coefficients */
			ediff = MAX2(100., ediff );
			rate = (realnum)pow(10. ,
				gbarcoll[nColl][0] + gbarcoll[nColl][1] * 
				pow(ediff,gbarcoll[nColl][2]) )*mole.lgColl_deexec_Calc;

			/*if( PRT_COLL )
				fprintf(ioQQQ,"col gbr\t%li\t%li\t%li\t%li\t%li\t%.4e\t%.4e\n",
				nColl+10,
				iVibHi,iRotHi,iVibLo,iRotLo,
				energy_wn[0][iVibHi][iRotHi] - energy_wn[0][iVibLo][iRotLo],
				H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][nColl] );*/
		}
		else
			rate = 0;
	}
	else if( nColl==4 )
	{
		/* collisions of H2 with protons - of this group, these are only
		 * that cause ortho - para conversion */
		/* >>refer	H2	coll Hp	Gerlich, D., 1990, J. Chem. Phys., 92, 2377-2388 */
		if( CollRateFit[iVibHi][iRotHi][iVibLo][iRotLo][nColl][1] != 0 )
		{
			rate = CollRateFit[iVibHi][iRotHi][iVibLo][iRotLo][nColl][0] * 1e-10f *
				/* sec fit coef was dE in milli eV */
				(realnum)sexp( CollRateFit[iVibHi][iRotHi][iVibLo][iRotLo][nColl][1]/1000./phycon.te_eV)*mole.lgColl_deexec_Calc;
			/*if( PRT_COLL )
				fprintf(ioQQQ,"col fit\t%li\t%li\t%li\t%li\t%li\t%.4e\t%.4e\n",
				nColl,
				iVibHi,iRotHi,iVibLo,iRotLo,
				energy_wn[0][iVibHi][iRotHi] - energy_wn[0][iVibLo][iRotLo],
				H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][nColl] );*/
		}
		/* this is option to use guess of rate coefficient for ortho-para
		 * conversion by collision with protons */
		/* turn mole.lgColl_gbar on/off with atom h2 gbar on off */
		else if( mole.lgColl_gbar )
		{
			/* the fit is log(K)=y_0+a*((x)^b), where K is the rate coefficient,
			* and x is the energy in wavenumbers */
			double ediff = energy_wn[0][iVibHi][iRotHi] - energy_wn[0][iVibLo][iRotLo];
			ediff = MAX2(100., ediff );
			rate = (realnum)pow(10. ,
				gbarcoll[nColl][0] + gbarcoll[nColl][1] * 
				pow(ediff ,gbarcoll[nColl][2])	)*mole.lgColl_deexec_Calc;

			/*if( PRT_COLL )
				fprintf(ioQQQ,"col gbr\t%li\t%li\t%li\t%li\t%li\t%.4e\t%.4e\n",
				nColl+10,
				iVibHi,iRotHi,iVibLo,iRotLo,
				energy_wn[0][iVibHi][iRotHi] - energy_wn[0][iVibLo][iRotLo],
				H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][nColl] );*/
		}
		else
			rate = 0;
	}
	/*>>chng 08 feb 27, GS*/
	else if( nColl == 2 && mole.lgH2_ORH2_ORNL )
	{
		/* ORNL in use */
		/* H2 - H2ortho collisions 
		* The H2 - H2ortho collisional data set is controlled by the 
		* atom H2 H2ortho collisions command
		* either mole.lgH2_ORH2_Meudon or mole.lgH2_PAH2_ORNL must be true*/
		if( (fitted=H2_ORH2_coll(ipHi , ipLo , phycon.te ))>0. )
		{
			rate = (realnum)fitted*mole.lgColl_deexec_Calc;
			if( PRT_COLL )
				fprintf(ioQQQ,"col fit\t%li\t%li\t%.4e\t%li\t%li\t%li\t%li\t%li\t%.4e\t%.4e\n",
				ipLo,ipHi,phycon.te,nColl,
				iVibHi,iRotHi,iVibLo,iRotLo,
				energy_wn[0][iVibHi][iRotHi] - energy_wn[0][iVibLo][iRotLo],
				H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][nColl] );
		}
		else
			rate = 0;
	}
	
	/*>>chng 08 feb 27, GS*/
	else if( nColl == 3 && mole.lgH2_PAH2_ORNL )
	{
		/* ORNL in use */
		/* H2 - H2ortho collisions 
		* The H2 - H2ortho collisional data set is controlled by the 
		* atom H2 H2ortho collisions command
		* either mole.lgH2_H2ortho_Meudon or mole.lgH2_H2ortho_ORNL must be true*/
		if( (fitted=H2_PAH2_coll(ipHi , ipLo , phycon.te ))>0. )
		{
			rate = (realnum)fitted*mole.lgColl_deexec_Calc;
			if( PRT_COLL )
				fprintf(ioQQQ,"col fit\t%li\t%li\t%.4e\t%li\t%li\t%li\t%li\t%li\t%.4e\t%.4e\n",
				ipLo,ipHi,phycon.te,nColl,
				iVibHi,iRotHi,iVibLo,iRotLo,
				energy_wn[0][iVibHi][iRotHi] - energy_wn[0][iVibLo][iRotLo],
				H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][nColl] );
		}
		else
			rate = 0;
	}
	
	else if( CollRateFit[iVibHi][iRotHi][iVibLo][iRotLo][nColl][0]!= 0 )
	{
		/* these are the fits from 
		 *>>refer	H2	coll	Le Bourlot, J., Pineau des Forets, 
		 *>>refercon	G., & Flower, D.R. 1999, MNRAS, 305, 802
		 * evaluate collision rates for those with real collision data */

		double r = CollRateFit[iVibHi][iRotHi][iVibLo][iRotLo][nColl][0] + 
			CollRateFit[iVibHi][iRotHi][iVibLo][iRotLo][nColl][1]/t3Plus1 + 
			CollRateFit[iVibHi][iRotHi][iVibLo][iRotLo][nColl][2]/t3Plus1Squared;

		rate = (realnum)pow(10.,r)*mole.lgColl_deexec_Calc;

		/*if( PRT_COLL )
			fprintf(ioQQQ,"col fit\t%li\t%li\t%li\t%li\t%li\t%.4e\t%.4e\n",
			nColl,
			iVibHi,iRotHi,iVibLo,iRotLo,
			energy_wn[0][iVibHi][iRotHi] - energy_wn[0][iVibLo][iRotLo],
			H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][nColl] );*/
	}
	/* this is option to use guess of collision rate coefficient - but only if this is 
	* a downward transition that does not mix ortho and para */
	/* turn mole.lgColl_gbar on/off with atom h2 gbar on off */
	else if( mole.lgColl_gbar  && 
		(H2_lgOrtho[0][iVibHi][iRotHi]-H2_lgOrtho[0][iVibLo][iRotLo]==0) )
	{
		/* the fit is log(K)=y_0+a*((x)^b), where K is the rate coefficient,
		 * and x is the energy in wavenumbers */
		double ediff = energy_wn[0][iVibHi][iRotHi] - energy_wn[0][iVibLo][iRotLo];
		/* do not let energy difference be smaller than 100 wn, the smallest
		 * difference for which we fit the rate coefficients */
		ediff = MAX2(100., ediff );
		rate = (realnum)pow(10. ,
			gbarcoll[nColl][0] + gbarcoll[nColl][1] * 
			pow(ediff,gbarcoll[nColl][2]) )*mole.lgColl_deexec_Calc;
		/* this is hack to change H2 H collision rates */
		if( nColl == 0 && h2.lgH2_H_coll_07 )
			H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][nColl] *= 100.;

		/*if( PRT_COLL )
			fprintf(ioQQQ,"col gbr\t%li\t%li\t%li\t%li\t%li\t%.4e\t%.4e\n",
			nColl+10,
			iVibHi,iRotHi,iVibLo,iRotLo,
			energy_wn[0][iVibHi][iRotHi] - energy_wn[0][iVibLo][iRotLo],
			H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][nColl] );*/
	}
	else
		rate = 0;


	/* >>chng 05 feb 09, add option to kill ortho - para collisions */
	if( !mole.lgH2_ortho_para_coll_on && 
		(H2_lgOrtho[0][iVibHi][iRotHi]-H2_lgOrtho[0][iVibLo][iRotLo]) )
		rate = 0.;

	{
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			fprintf(ioQQQ,"bugcoll\tiVibHi\t%li\tiRotHi\t%li\tiVibLo\t%li\tiRotLo\t%li\tcoll\t%.2e\n",
				iVibHi,iRotHi,iVibLo,iRotLo,
				H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][nColl] );
		}
	}
	return rate;
}

/*H2_CollidRateEvalAll - set H2 collision rate coefficients */
void H2_CollidRateEvalAll( void )
{
	long int numb_coll_trans = 0;
	double excit;
	long int iElecHi , iElecLo , ipHi , iVibHi , iRotHi , 
		ipLo , iVibLo , iRotLo , nColl;

	DEBUG_ENTRY( "H2_CollidRateEvalAll()" );

	if( PRT_COLL )
		fprintf(ioQQQ,"H2 coll deex rate coef\n"
		"VibHi\tRotHi\tVibLo\tRotLo\trate\n");

	iElecHi = 0;
	iElecLo = 0;
	if(mole.nH2_TRACE >= mole.nH2_trace_full) 
		fprintf(ioQQQ,"H2 set collision rates\n");
	/* loop over all possible collisional changes within X 
	* and set collision rates, which only depend on Te
	* will go through array in energy order since coll trans do not
	* correspond to a line 
	* collisional dissociation rate coefficient, units cm3 s-1 */
	H2_coll_dissoc_rate_coef[0][0] = 0.;
	H2_coll_dissoc_rate_coef_H2[0][0] = 0.;
	for( ipHi=0; ipHi<nLevels_per_elec[0]; ++ipHi )
	{
		double energy;

		/* obtain the proper indices for the upper level */
		long int ip = H2_ipX_ener_sort[ipHi];
		iVibHi = ipVib_H2_energy_sort[ip];
		iRotHi = ipRot_H2_energy_sort[ip];

		/* this is a guess of the collisional dissociation rate coefficient -
		 * will be multiplied by the sum of densities of all colliders 
		 * except H2*/
		energy = H2_DissocEnergies[0] - energy_wn[0][iVibHi][iRotHi];
		ASSERT( energy > 0. );
		/* we made this up - Boltzmann factor times rough coefficient */
		H2_coll_dissoc_rate_coef[iVibHi][iRotHi] = 
			1e-14f * (realnum)sexp(energy/phycon.te_wn) * mole.lgColl_dissoc_coll;

		/* collisions with H2 - pre coefficient changed from 1e-8 
		 * (from umist) to 1e-11 as per extensive discussion with Phillip Stancil */
		H2_coll_dissoc_rate_coef_H2[iVibHi][iRotHi] = 
			1e-11f * (realnum)sexp(energy/phycon.te_wn) * mole.lgColl_dissoc_coll;

		/*fprintf(ioQQQ,"DEBUG coll_dissoc_rateee\t%li\t%li\t%.3e\t%.3e\n",
			iVibHi,iRotHi,
			H2_coll_dissoc_rate_coef[iVibHi][iRotHi],
			H2_coll_dissoc_rate_coef_H2[iVibHi][iRotHi]);*/

		for( ipLo=0; ipLo<ipHi; ++ipLo )
		{
			ip = H2_ipX_ener_sort[ipLo];
			iVibLo = ipVib_H2_energy_sort[ip];
			iRotLo = ipRot_H2_energy_sort[ip];

			ASSERT( energy_wn[0][iVibHi][iRotHi] - energy_wn[0][iVibLo][iRotLo] > 0.);

			/* in following the colliders are H, He, H2(ortho), H2(para), and H+ */
			/* fits were read in from the following files: "H2_coll_H.dat" ,
			* "H2_coll_He.dat" , "H2_coll_H2ortho.dat" ,"H2_coll_H2para.dat",
			* "H2_coll_Hp.dat" */

			/* keep track of number of different collision routes */
			++numb_coll_trans;
			/* this is sum over all different colliders, except last two which are special,
			* linear rather than log formula for that one for second to last,
			* and special fitting formula for last */
			if( PRT_COLL  && iVibHi == 1 && iRotHi==3)
				fprintf(ioQQQ,"%li\t%li\t%li\t%li",
				iVibHi,iRotHi,iVibLo,iRotLo);
			for( nColl=0; nColl<N_X_COLLIDER; ++nColl )
			{
				H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][nColl] = 
					H2_CollidRateEvalOne( iVibHi,iRotHi,iVibLo,iRotLo,
					ipHi , ipLo , nColl );
				if( PRT_COLL   && iVibHi == 1 && iRotHi==3)
					fprintf(ioQQQ,"\t%.2e",H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][nColl] );
			} 
			if( PRT_COLL   && iVibHi == 1 && iRotHi==3)
				fprintf(ioQQQ,"\n");
		}
	}
	if( PRT_COLL )
		cdEXIT( EXIT_FAILURE );

	/* at this stage the collision rates that came in from the large data files
	* have been entered into the H2_CollRate array.  Now add on three extra collision
	* terms, the ortho para atomic H collision rates from
	* >>>refer	H2	collision	Sun, Y., & Dalgarno, A., 1994, ApJ, 427, 1053-1056 
	*/
	nColl = 0;
	iElecHi = 0;
	iElecLo = 0;
	iVibHi = 0;
	iVibLo = 0;

	/* >>chng 02 nov 13, the Sun and Dalgarno rates diverge to + inf below this temp */
	/* >>chng 05 feb 09, do not return zero when T < 100 - instead, don't let T fall below 100 */
	double excit1;
	double te_used = MAX2( 100. , phycon.te );
	/* this is the J=1-0 downward collision rate */
	iRotLo = 0;
	iRotHi = 1;
	excit1 = sexp( H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].EnergyK/te_used);
	excit = sexp( -(POW2(5.30-460./te_used)-21.2) )*1e-13;

	H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][0] = (realnum)(
		excit*H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Lo->g/
		H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Hi->g / 
		/* >>chng 02 nov 13, from 2nd to first */
		SDIV(excit1) )*mole.lgColl_deexec_Calc *
		/* option to disable ortho-para conversion by coll with grains */
		mole.lgH2_ortho_para_coll_on;

	/*if( PRT_COLL )
		fprintf(ioQQQ,"col o-p\t%li\t%li\t%li\t%li\t%li\t%.4e\t%.4e\n",
		nColl,
		iVibHi,iRotHi,iVibLo,iRotLo,
		energy_wn[0][iVibHi][iRotHi] - energy_wn[0][iVibLo][iRotLo],
		H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][nColl] );*/

	/* this is the J=3-0 downward collision rate */
	iRotLo = 0;
	iRotHi = 3;
	excit = sexp( -(POW2(6.36-373./te_used)-34.5) )*1e-13;
	excit1 = sexp( H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].EnergyK/te_used);
	H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][0] = (realnum)(
		excit*H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Lo->g/
		H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Hi->g / 
		SDIV(excit1) )*mole.lgColl_deexec_Calc *
		/* option to disable ortho-para conversion by coll with grains */
		mole.lgH2_ortho_para_coll_on;

	/*if( PRT_COLL )
		fprintf(ioQQQ,"col o-p\t%li\t%li\t%li\t%li\t%li\t%.4e\t%.4e\n",
		nColl,
		iVibHi,iRotHi,iVibLo,iRotLo,
		energy_wn[0][iVibHi][iRotHi] - energy_wn[0][iVibLo][iRotLo],
		H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][nColl] );*/

	/* this is the downward J=2-1 collision rate */
	iRotLo = 1;
	iRotHi = 2;
	excit = sexp( -(POW2(5.35-454./te_used)-23.1 ) )*1e-13;
	excit1 = sexp( H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].EnergyK/te_used);
	H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][0] = (realnum)(
		excit*H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Lo->g/
		H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Hi->g / 
		SDIV(excit1) )*mole.lgColl_deexec_Calc *
		/* option to disable ortho-para conversion by coll with grains */
		mole.lgH2_ortho_para_coll_on;

	/* >>chng 05 nov 30, GS, rates decreases exponentially for low temperature, see Le Bourlot et al. 1999  */
	/* Phillips mail--Apparently, the SD fit is only valid over the range of their  calculations, 100-1000K. 
	* The rate should continue to fall exponentially with decreasing T, something like exp(-3900/T) for 0->1 and  
	* exp[-(3900-170.5)/T] for 1->0. It is definitely, not constant for T  lower than 100 K, as far as we know. 
	* There may actually be a quantum tunneling effect which causes the rate to increase at lower T, but no  
	* one has calculated it (as far as I know) and it might happen below 1K or  so.???*/
	if( phycon.te < 100. )
	{
		/* first term in exp is suggested by Phillip, second temps in paren is to ensure continuity
		* across 100K */
		H2_CollRate[0][1][0][0][0] = (realnum)(H2_CollRate[0][0][1][0][0]*exp(-(3900-170.5)*(1./phycon.te - 1./100.)));
		H2_CollRate[0][3][0][0][0] = (realnum)(H2_CollRate[0][0][3][0][0]*exp(-(3900-1015.1)*(1./phycon.te - 1./100.)));
		H2_CollRate[0][2][0][1][0] = (realnum)(H2_CollRate[0][0][2][0][1]*exp(-(3900-339.3)*(1./phycon.te - 1./100.)));
	}

	/*if( PRT_COLL )
		fprintf(ioQQQ,"col o-p\t%li\t%li\t%li\t%li\t%li\t%.4e\t%.4e\n",
		nColl,
		iVibHi,iRotHi,iVibLo,iRotLo,
		energy_wn[0][iVibHi][iRotHi] - energy_wn[0][iVibLo][iRotLo],
		H2_CollRate[iVibHi][iRotHi][iVibLo][iRotLo][nColl] );*/

	if( mole.nH2_TRACE >= mole.nH2_trace_full )
		fprintf(ioQQQ,
		" collision rates updated for new temp, number of trans is %li\n",
		numb_coll_trans);

	/* quit it we are only printing - but do this after printing coll rates
	if( PRT_COLL )
		cdEXIT( EXIT_FAILURE ); */
	return;
}

/*H2_CollidRateRead read collision rates */
void H2_CollidRateRead( long int nColl )
{
	/* the colliders are H, He, H2 ortho, H2 para, H+ 
	 * these are the default file names.  they can be overridden with
	 * the SET ATOMIC DATA MOLECULE H2 command */

	/*NB thesee must be kept parallel with labels in chH2ColliderLabels */

	const char* cdDATAFILE[N_X_COLLIDER] = 
	{
		/* 0 */"H2_coll_H_07.dat",
		/* 1 */"H2_coll_He_LeBourlot.dat", 
		/* 2 */"H2_coll_H2ortho.dat",
		/* 3 */"H2_coll_H2para.dat",
		/* 4 */"H2_coll_Hp.dat"
	};

	FILE *ioDATA;
	char chLine[FILENAME_PATH_LENGTH_2];
	const char* chFilename;
	long int i, n1;
	long int iVibHi , iVibLo , iRotHi , iRotLo;
	long int magic_expect = -1;

	DEBUG_ENTRY( "H2_CollidRateRead()" );

	if( nColl == 0 )
	{
		/* which H2 - H data set? */
		if( h2.lgH2_H_coll_07 )
		{
			magic_expect = 71106;
			chFilename = "H2_coll_H_07.dat";
		}
		else
		{
			/* the 1999 data set */
			magic_expect = 20429;
			chFilename = "H2_coll_H_99.dat";
		}
	}
	else if( nColl == 1 )
	{
		/* which H2 - He data set? */
		if( mole.lgH2_He_ORNL )
		{
			/* >>chng 07 may 12, update magic number when one transition in H2 - H2 file was
			* replaced - the fitting coefficients had terms of order -8e138 which fpe in float */
			long int magic_found,
				magic_expect = 70513;
			/* special case, new data file from Oak Ridge project -
			* call init routine and return - data is always read in when large H2 is included -
			* but data are only used (for now, mid 2005) when command 
			* atom H2 He OLD (Meudon) NEW (Stancil) and OFF given */
			/*H2_He_coll_init receives the name of the file that contains the fitting coefficients 
			* of all transitions and read into 3d vectors. It outputs 'test.out' to test the arrays*/
			chFilename = "H2_coll_He_ORNL.dat";
			if( (magic_found = H2_He_coll_init( chFilename )) != magic_expect )
			{
				fprintf(ioQQQ,"The H2 - He collision data file H2_coll_He_ORNL.dat does not have the correct magic number.\n");
				fprintf(ioQQQ,"I found %li but expected %li\n", magic_found , magic_expect );

				cdEXIT(EXIT_FAILURE);
			}
			return;
		}
		else
		{
			magic_expect = 20429;
			/* the Le Bourlot et al data */
			chFilename = "H2_coll_He_LeBourlot.dat";
		}
	}
	/*>>chng 08 feb 27, GS*/
	else if( nColl == 2 )
	{
		/* which H2 - He data set? */
		if( mole.lgH2_ORH2_ORNL )
		{
			long int magic_found,
				magic_expect = 80227;
			/* special case, new data file from Oak Ridge project -
			* call init routine and return - data is always read in when large H2 is included -
			* but data are only used (for now, mid 2005) when command 
			* atom H2 H2 OLD (Meudon) NEW (Stancil) and OFF given */
			
			char chPath[FILENAME_PATH_LENGTH_2];
			strcpy( chPath, "h2" );
			strcat( chPath, input.chDelimiter );
			strcat( chPath, "H2_coll_H2ortho_ORNL.dat" );
			if( (magic_found = H2_ORH2_coll_init( chPath )) != magic_expect )
			{
				fprintf(ioQQQ,"The H2 - H2 collision data file H2_coll_H2ortho_ORNL.dat does not have the correct magic number.\n");
				fprintf(ioQQQ,"I found %li but expected %li\n", magic_found , magic_expect );

				cdEXIT(EXIT_FAILURE);
			}
			return;
		}
		else
		{
			magic_expect = 20429;
			/* the Le Bourlot et al data */
			chFilename = "H2_coll_H2ortho_LeBourlot.dat";
		}
	}

	/*>>chng 08 feb 27, GS*/
	else if( nColl == 3 )
	{
		/* which H2 - H2 data set? */
		if( mole.lgH2_PAH2_ORNL )
		{
			long int magic_found,
				magic_expect = 80227;
			
			char chPath[FILENAME_PATH_LENGTH_2];
			strcpy( chPath, "h2" );
			strcat( chPath, input.chDelimiter );
			strcat( chPath, "H2_coll_H2para_ORNL.dat" );
			if( (magic_found = H2_PAH2_coll_init( chPath )) != magic_expect )
			{
				fprintf(ioQQQ,"The H2 - H2 collision data file H2_coll_H2para_ORNL.dat does not have the correct magic number.\n");
				fprintf(ioQQQ,"I found %li but expected %li\n", magic_found , magic_expect );

				cdEXIT(EXIT_FAILURE);
			}
			return;
		}
		else
		{
			magic_expect = 20429;
			/* the Le Bourlot et al data */
			chFilename = "H2_coll_H2para_LeBourlot.dat";
		}
	}
	else
	{
		/* files with only one version */
		magic_expect = 20429;
		chFilename = cdDATAFILE[nColl];
	}

	char chPath[FILENAME_PATH_LENGTH_2];
	strcpy( chPath, "h2" );
	strcat( chPath, input.chDelimiter );
	strcat( chPath, chFilename );
	ioDATA = open_data( chPath, "r" );

	/* read the first line and check that magic number is ok */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " H2_CollidRateRead could not read first line of %s\n", chFilename );
		cdEXIT(EXIT_FAILURE);
	}

	/* magic number */
	n1 = atoi( chLine );

	/* magic number
	* the following is the set of numbers that appear at the start of level1.dat 01 08 10 */
	if( n1 != magic_expect )
	{
		fprintf( ioQQQ, 
			 " H2_CollidRateRead: the version of %s is not the current version.\n", chFilename );
		fprintf( ioQQQ, 
			" I expected to find the number %li and got %li instead.\n" ,
			magic_expect , n1 );
		fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
		cdEXIT(EXIT_FAILURE);
	}

	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
	{
		if( chLine[0]=='#'  )
			continue;
		double a[3];
		sscanf(chLine,"%li\t%li\t%li\t%li\t%le\t%le\t%le", 
			&iVibHi ,&iRotHi , &iVibLo , &iRotLo , &a[0],&a[1],&a[2] );
		/*fprintf(ioQQQ,"DEBUG %li %li %li %li \n",
		iVibHi , iRotHi , iVibLo , iRotLo);*/
		ASSERT( iRotHi >= 0 && iVibHi >= 0 && iRotLo >= 0 && iVibLo >=0 );

		/* check that we actually included the levels in the model representation
		* depending on the potential surface, the collision date set
		* may not agree with our adopted model - skip those */
		if( iVibHi > VIB_COLLID || 
			iVibLo > VIB_COLLID ||
			iRotHi > h2.nRot_hi[0][iVibHi] ||
			iRotLo > h2.nRot_hi[0][iVibLo])
		{
			iVibHi = -1;
			continue;
		}

		/* some collision rates have the same upper and lower indices - skip them */
		if( !( (iVibHi == iVibLo) && (iRotHi == iRotLo  )) )
		{
			/* this is downward transition - make sure that the energy difference is positive */
			ASSERT( (energy_wn[0][iVibHi][iRotHi] - energy_wn[0][iVibLo][iRotLo] ) > 0. );
			for( i=0; i<3; ++i )
			{
				CollRateFit[iVibHi][iRotHi][iVibLo][iRotLo][nColl][i] = (realnum)a[i];
			}

			/* this prints all levels with rates 
			fprintf(ioQQQ,"no\t%li\t%li\t%li\t%li\t%.2e\t%.2e\t%.2e\n", 
			iVibHi,iRotHi,iVibLo,iRotLo,a[0],a[1],a[2]);*/
		}
	}
	fclose( ioDATA );

	return;
}
/* end of H2_CollidRateRead */

/* This code was written by Terry Yun, 2005 */ 
/* H2_He_coll_init receives the name of the file that contains the fitting     
* coefficients of all transitions and read into 3d vectors. 
* It returns magic number*/
long int H2_He_coll_init(const char FILE_NAME_IN[])
{

	int i, j;   
	long int magic;
	double par[N_H2_HE_FIT_PAR];
	char line[INPUT_LINE_LENGTH];
	int h2_i, h2_f, he_i, he_f;/* target, projectile initial and final indices */
	char quality, space;
	/*'space' variable is for reading spaces between characters */
	/* scanf can skip spaces(or tap) when it reads numbers, but not for characters. */
	double error;

	FILE *ifp;

	DEBUG_ENTRY( "H2_He_coll_init()" );

	/* create space for two arrays - H2_He_coll_fit_par will contain the
	* fit parameters in form [init state][final state][param]
	* lgDefn_H2He_coll is flag saying whether data exists */
	H2_He_coll_fit_par = (realnum ***)MALLOC(sizeof(realnum **)*(unsigned)nLevels_per_elec[0] );
	lgDefn_H2He_coll = (bool**)MALLOC(sizeof(bool *)*(unsigned)nLevels_per_elec[0] );
	for( i=1; i<nLevels_per_elec[0]; ++i )
	{
		lgDefn_H2He_coll[i] = (bool*)MALLOC(sizeof(bool)*(unsigned)i );
		H2_He_coll_fit_par[i] = (realnum **)MALLOC(sizeof(realnum *)*(unsigned)i );
		for( j=0; j<i; ++j)
			H2_He_coll_fit_par[i][j] = (realnum *)MALLOC(sizeof(realnum)*(unsigned)N_H2_HE_FIT_PAR );
	}

	/* set a flag saying whether data exits: initially everything is '0' */
	for( i=1; i<nLevels_per_elec[0]; i++ )
	{
		for( j=0; j<i; j++ )
		{
			lgDefn_H2He_coll[i][j] = false;
		}
	}

	/*open the input file and put into 3d arrays*/
	char chPath[FILENAME_PATH_LENGTH_2];
	strcpy( chPath, "h2" );
	strcat( chPath, input.chDelimiter );
	strcat( chPath, FILE_NAME_IN );
	ifp = open_data( chPath, "r" );

	/*read the magic number*/
	if( read_whole_line(line, (int)sizeof(line), ifp)==NULL )
	{
		printf("DISASTER H2_He_coll_init can't read first line of %s\n", FILE_NAME_IN);
		cdEXIT(EXIT_FAILURE);
	}
	sscanf(line, "%li", &magic); 

	/* read the file until the end of line */
	while( read_whole_line(line, (int)sizeof(line), ifp) != NULL )
	{
		/* skip any line that starts with '#' */
		if( line[0]!='#' )
		{
			sscanf(line, "%i%i%i%i%c%c%c%c%lf%lf%lf%lf%lf%lf%lf%lf%lf", &h2_i, 
				&h2_f, &he_i, &he_f,&space, &space, &space, &quality, 
				&error, &par[0], &par[1], &par[2], &par[3], &par[4], 
				&par[5], &par[6], &par[7]); 

			/* >>chng 07 may 28, extensive testing showed rate > 1e38 for
			* temperatures where fitting equation hit a pole.  Exchange
			* with Phillip Stancil:
			"You are correct. The problem is the parameter d which is 
			negative in all the problem cases you found. A pole occurs 
			when d*T/1000 + 1 =0.

			I found that there are a total of 11 transitions with d<1. 
			A quick solution is to take abs(d). I checked one case and the 
			resulting function went smoothly through the pole. On either 
			side of the pole the two functions looked the same on a log-log plot.

			LeBourlot used a similar function, but had d=1. So, they never 
			got a pole. I guess I tried to make the function a little too 
			flexible.

			There is a similar term with g*T/1000 + 1. There are no fits 
			with a negative g, but you might take abs(g) to be on the space 
			side for the future."
			*/
			if( par[3] < 0. )
			{
				/*fprintf(ioQQQ,"DEBUG negative [3]=%.2e\n", par[3]);*/
				par[3] = -par[3];
			}
			if( par[6] < 0. )
			{
				/*fprintf(ioQQQ,"DEBUG negative [6]=%.2e\n", par[6]);*/
				par[6] = -par[6];
			}
			/* input file is on physics or Fortran scale so lowest energy
			* index is 1.  Store on C scale */
			--h2_f;
			--h2_i;
			/* set true when the indices are defined */
			lgDefn_H2He_coll[h2_i][h2_f] = true;
			for( i=0; i<N_H2_HE_FIT_PAR; i++ )
				/* assigning the parameters to 3d array */
				H2_He_coll_fit_par[h2_i][h2_f][i] = (realnum)par[i];
		}
	}
	fclose(ifp);  
	return magic;
}

/*>>chng 08 feb 27, GS*/
long int H2_ORH2_coll_init(const char FILE_NAME_IN[] )
{

	int i, j;   
	long int magic;
	double par[N_H2_ORH2_FIT_PAR];
	char line[INPUT_LINE_LENGTH];
	int h2_i, h2_f, h2ortho_i, h2ortho_f;/* target, projectile initial and final indices */
	char quality, space;
	/*'space' variable is for reading spaces between characters */
	/* scanf can skip spaces(or tap) when it reads numbers, but not for characters. */
	double error;

	FILE *ifp;

	DEBUG_ENTRY( "H2_ORH2_coll_init()" );

	/* create space for two arrays - H2_H2_coll_fit_par will contain the
	* fit parameters in form [init state][final state][param]
	* lgDefn_H2H2_coll is flag saying whether data exists */
	H2_ORH2_coll_fit_par = (realnum ***)MALLOC(sizeof(realnum **)*(unsigned)nLevels_per_elec[0] );
	lgDefn_H2ORH2_coll = (bool**)MALLOC(sizeof(bool *)*(unsigned)nLevels_per_elec[0] );
	for( i=1; i<nLevels_per_elec[0]; ++i )
	{
		lgDefn_H2ORH2_coll[i] = (bool*)MALLOC(sizeof(bool)*(unsigned)i );
		H2_ORH2_coll_fit_par[i] = (realnum **)MALLOC(sizeof(realnum *)*(unsigned)i );
		for( j=0; j<i; ++j)
			H2_ORH2_coll_fit_par[i][j] = (realnum *)MALLOC(sizeof(realnum)*(unsigned)N_H2_ORH2_FIT_PAR );
	}

	/* set a flag saying whether data exits: initially everything is '0' */
	for( i=1; i<nLevels_per_elec[0]; i++ )
	{
		for( j=0; j<i; j++ )
		{
			lgDefn_H2ORH2_coll[i][j] = false;
		}
	}

	/*open the input file and put into 3d arrays*/
	ifp = open_data( FILE_NAME_IN, "r" );

	/*read the magic number*/
	if( read_whole_line(line, (int)sizeof(line), ifp)==NULL )
	{
		printf("DISASTER H2_ORH2_coll_init can't read first line of %s\n", FILE_NAME_IN);
		cdEXIT(EXIT_FAILURE);
	}
	sscanf(line, "%li", &magic); 
	
	/* read the file until the end of line */
	while( read_whole_line(line, (int)sizeof(line), ifp) != NULL )
	{
		/* skip any line that starts with '#' */
		if( line[0]!='#' )
		{
			sscanf(line, "%i%i%i%i%c%c%c%c%lf%lf%lf%lf%lf%lf%lf%lf%lf", &h2_i, 
				&h2_f, &h2ortho_i, &h2ortho_f,&space, &space, &space, &quality, 
				&error, &par[0], &par[1], &par[2], &par[3], &par[4], 
				&par[5], &par[6], &par[7]); 

			
			if( par[3] < 0. )
			{
				/*fprintf(ioQQQ,"DEBUG negative [3]=%.2e\n", par[3]);*/
				par[3] = -par[3];
			}
			if( par[6] < 0. )
			{
				/*fprintf(ioQQQ,"DEBUG negative [6]=%.2e\n", par[6]);*/
				par[6] = -par[6];
			}
			/* input file is on physics or Fortran scale so lowest energy
			* index is 1.  Store on C scale */
			--h2_f;
			--h2_i;
			/* set true when the indices are defined */
			lgDefn_H2ORH2_coll[h2_i][h2_f] = true;
			for( i=0; i<N_H2_ORH2_FIT_PAR; i++ )
				/* assigning the parameters to 3d array */
				H2_ORH2_coll_fit_par[h2_i][h2_f][i] = (realnum)par[i];
		}
	}
	fclose(ifp);  
	return magic;
}

/*>>chng 08 feb 27, GS*/
long int H2_PAH2_coll_init(const char FILE_NAME_IN[] )
{

	int i, j;   
	long int magic;
	double par[N_H2_PAH2_FIT_PAR];
	char line[INPUT_LINE_LENGTH];
	int h2_i, h2_f, h2para_i, h2para_f;/* target, projectile initial and final indices */
	char quality, space;
	/*'space' variable is for reading spaces between characters */
	/* scanf can skip spaces(or tap) when it reads numbers, but not for characters. */
	double error;

	FILE *ifp;

	DEBUG_ENTRY( "H2_PAH2_coll_init()" );

	/* create space for two arrays - H2_H2_coll_fit_par will contain the
	* fit parameters in form [init state][final state][param]
	* lgDefn_H2H2_coll is flag saying whether data exists */
	H2_PAH2_coll_fit_par = (realnum ***)MALLOC(sizeof(realnum **)*(unsigned)nLevels_per_elec[0] );
	lgDefn_H2PAH2_coll = (bool**)MALLOC(sizeof(bool *)*(unsigned)nLevels_per_elec[0] );
	for( i=1; i<nLevels_per_elec[0]; ++i )
	{
		lgDefn_H2PAH2_coll[i] = (bool*)MALLOC(sizeof(bool)*(unsigned)i );
		H2_PAH2_coll_fit_par[i] = (realnum **)MALLOC(sizeof(realnum *)*(unsigned)i );
		for( j=0; j<i; ++j)
			H2_PAH2_coll_fit_par[i][j] = (realnum *)MALLOC(sizeof(realnum)*(unsigned)N_H2_PAH2_FIT_PAR );
	}

	/* set a flag saying whether data exits: initially everything is '0' */
	for( i=1; i<nLevels_per_elec[0]; i++ )
	{
		for( j=0; j<i; j++ )
		{
			lgDefn_H2PAH2_coll[i][j] = false;
		}
	}

	/*open the input file and put into 3d arrays*/
	ifp = open_data( FILE_NAME_IN, "r" );
	/*read the magic number*/
	if( read_whole_line(line, (int)sizeof(line), ifp)==NULL )
	{
		printf("DISASTER H2_PAH2_coll_init can't read first line of %s\n", FILE_NAME_IN);
		cdEXIT(EXIT_FAILURE);
	}
	sscanf(line, "%li", &magic); 

	/* read the file until the end of line */
	while( read_whole_line(line, (int)sizeof(line), ifp) != NULL )
	{
		/* skip any line that starts with '#' */
		if( line[0]!='#' )
		{
			sscanf(line, "%i%i%i%i%c%c%c%c%lf%lf%lf%lf%lf%lf%lf%lf%lf", &h2_i, 
				&h2_f, &h2para_i, &h2para_f,&space, &space, &space, &quality, 
				&error, &par[0], &par[1], &par[2], &par[3], &par[4], 
				&par[5], &par[6], &par[7]); 

			if( par[3] < 0. )
			{
				/*fprintf(ioQQQ,"DEBUG negative [3]=%.2e\n", par[3]);*/
				par[3] = -par[3];
			}
			if( par[6] < 0. )
			{
				/*fprintf(ioQQQ,"DEBUG negative [6]=%.2e\n", par[6]);*/
				par[6] = -par[6];
			}
			/* input file is on physics or Fortran scale so lowest energy
			* index is 1.  Store on C scale */
			--h2_f;
			--h2_i;
			/* set true when the indices are defined */
			lgDefn_H2PAH2_coll[h2_i][h2_f] = true;
			for( i=0; i<N_H2_PAH2_FIT_PAR; i++ )
				/* assigning the parameters to 3d array */
				H2_PAH2_coll_fit_par[h2_i][h2_f][i] = (realnum)par[i];
		}
	}
	fclose(ifp);  
	return magic;
}

/* This code was written by Terry Yun, 2005 */ 
/* H2_He_coll Interpolate the rate coefficients 
* It receives initial and final transition and temperature
* The range of the temperature is between 2K - 1e8K */
double H2_He_coll(int init, int final, double temp)
{

	double k, t3Plus1, b2, c2, f2, t3;

	DEBUG_ENTRY( "H2_He_coll()" );

	/* Invalid entries returns '-1':the initial indices are smaller than the final indices */
	if( temp<2 || temp > 1e8 )
	{
		k = -1;
	}
	else if( init <= final )
	{
		k = -1;
	}
	/* Invalid returns '-1': the indices are greater than 302 or smaller than 0 */
	else if( init < 0 || init >302 || final < 0 || final > 302 )
	{
		k = -1;
	}
	/* Undefined indices returns '0' */
	else if( !lgDefn_H2He_coll[init][final] )
	{
		k = -1;
	}
	/* defined indices */
	else if( lgDefn_H2He_coll[init][final] )
	{
		/* the fitting equation we used:
		* k_(vj,v'j') = 10^(a + b/(d*T/10^3+1) + c/t^2)
		* + 10^(e + f/(g*T/10^3+1)**h)
		* T in K, k in cm3/s. */
		double sum1 , sum2;

		/* this kludge is to bypass errors in fitting formula - there
		* is a pole possible at around 1e6 K and rates can diverge
		* and high temperatures */
		/** \todo	1	fix this hack - Phillip Stancil is refitting */
		temp = MIN2( 1e4 , temp );

		t3 = temp/1e3; 
		/* t3Plus1 = T/1000 + 1*/
		t3Plus1 = t3+1; 
		b2 = H2_He_coll_fit_par[init][final][1]/(H2_He_coll_fit_par[init][final][3]*t3+1);
		c2 = H2_He_coll_fit_par[init][final][2]/(t3Plus1*t3Plus1);
		{
			enum {DEBUG_LOC=false};
			if( DEBUG_LOC )
			{
				fprintf(ioQQQ,"bug H2 He coll\t%i %i %.3e %.3e %.3e \n",
					init,final,
					H2_He_coll_fit_par[init][final][5],
					H2_He_coll_fit_par[init][final][6]*t3+1, 
					H2_He_coll_fit_par[init][final][7]
				/*pow(H2_He_coll_fit_par[init][final][6]*t3+1,H2_He_coll_fit_par[init][final][7])*/
				);
			}
		}
		/* this is log of f2 - see whether it is within bounds */
		sum1 = H2_He_coll_fit_par[init][final][7] * log10( H2_He_coll_fit_par[init][final][6]*t3+1. );
		/* this protects against overflow */
		if( fabs(sum1)< 38. )
		{
			/* >>chng 06 jun 07, Mitchell Martin, from following to equivalent below.  This was basically a
			* bug in the pgcc compiler that caused h2_pdr_leiden_f1 to throw an fpe
			* the changed code is equivalent */
			/*f2 = H2_He_coll_fit_par[init][final][5]/pow(H2_He_coll_fit_par[init][final][6]*t3+1.,H2_He_coll_fit_par[init][final][7]);*/
			f2 = H2_He_coll_fit_par[init][final][5]/pow(10. , sum1 );
		}
		else
			f2= 0.;
		{
			enum {DEBUG_LOC=false };
			if( DEBUG_LOC )
			{
				fprintf(ioQQQ,"bug H2 He coll\t%i %i %.3e %.3e %.3e %.3e %.3e sum %.3e %.3e \n",
					init,final,
					H2_He_coll_fit_par[init][final][0],
					b2, 
					c2,
					H2_He_coll_fit_par[init][final][4],
					f2  ,
					H2_He_coll_fit_par[init][final][0]+b2+c2,
					H2_He_coll_fit_par[init][final][4]+f2);
			}
		}

		sum1 = H2_He_coll_fit_par[init][final][0]+b2+c2;
		sum2 = H2_He_coll_fit_par[init][final][4]+f2;
		k = 0.;
		if( fabs(sum1) < 38. )
			k += pow(10., sum1 );
		if( fabs(sum2) < 38. )
			k += pow(10., sum2 );

		if( k>1e-6 )
		{
			/* rate of 1e-6 is largest possible rate according to Phillip
			* Stancil email 07 may 27 */
			fprintf(ioQQQ,"PROBLEM H2-He rate coefficient hit cap, upper index of %i"
				" lower index of %i rate was %.2e resetting to 1e-6, Te=%e\n",
				init , final , k , phycon.te );
			k = 1e-6;
		}
		{
			enum {DEBUG_LOC=false};
			if( DEBUG_LOC )
			{
				fprintf(ioQQQ,"DEBUG H2 He rate %.3e \n",
					k);
			}
		}
	}
	/*unknown invalid entries returns '99'*/
	else
	{
		k = 99;
	}
	return k;
}
/*chng 08 feb 27, GS*/
double H2_ORH2_coll(int init, int final, double temp)
{

	double k, t, b2, c2, f2, t2;

	DEBUG_ENTRY( "H2_ORH2_coll()" );

	/* Invalid entries returns '-1':the initial indices are smaller than the final indices */
	if( temp<1 || temp > 1e4 )
	{
		k = -1;
	}
	else if( init <= final )
	{
		k = -1;
	}
	/* Invalid returns '-1': the indices are greater than 302 or smaller than 0 */
	else if( init < 0 || init >302 || final < 0 || final > 302 )
	{
		k = -1;
	}
	/* Undefined indices returns '0' */
	else if( !lgDefn_H2ORH2_coll[init][final] )
	{
		k = -1;
	}
	/* defined indices */
	else if( lgDefn_H2ORH2_coll[init][final] )
	{
		/* the fitting equation we used:
		* k_(vj,v'j') = 10^(a + b/(d*T/10^3+1) + c/t^2)
		* + 10^(e + f/(g*T/10^3+1)**h)
		* T in K, k in cm3/s. */
		double sum1 , sum2;

		
		t2 = temp/1e3; 
		/* t = T*10^-3 + 1*/
		t = t2+1; 
		b2 = H2_ORH2_coll_fit_par[init][final][1]/(H2_ORH2_coll_fit_par[init][final][3]*t2+1);
		c2 = H2_ORH2_coll_fit_par[init][final][2]/(t*t);
		{
			enum {DEBUG_LOC=false};
			if( DEBUG_LOC )
			{
				fprintf(ioQQQ,"bug H2 H2ortho coll\t%i %i %.3e %.3e %.3e \n",
					init,final,
					H2_ORH2_coll_fit_par[init][final][5],
					H2_ORH2_coll_fit_par[init][final][6]*t2+1, 
					H2_ORH2_coll_fit_par[init][final][7]
				/*pow(H2_ORH2_coll_fit_par[init][final][6]*t2+1,H2_ORH2_coll_fit_par[init][final][7])*/
				);
			}
		}
		/* this is log of f2 - see whether it is within bounds */
		sum1 = H2_ORH2_coll_fit_par[init][final][7] * log10( H2_ORH2_coll_fit_par[init][final][6]*t2+1. );
		/* this protects against overflow */
		if( fabs(sum1)< 38. )
		{
			
			f2 = H2_ORH2_coll_fit_par[init][final][5]/pow(10. , sum1 );
		}
		else
			f2= 0.;
		{
			enum {DEBUG_LOC=false };
			if( DEBUG_LOC )
			{
				fprintf(ioQQQ,"bug H2 H2ortho coll\t%i %i %.3e %.3e %.3e %.3e %.3e sum %.3e %.3e \n",
					init,final,
					H2_ORH2_coll_fit_par[init][final][0],
					b2, 
					c2,
					H2_ORH2_coll_fit_par[init][final][4],
					f2  ,
					H2_ORH2_coll_fit_par[init][final][0]+b2+c2,
					H2_ORH2_coll_fit_par[init][final][4]+f2);
			}
		}

		sum1 = H2_ORH2_coll_fit_par[init][final][0]+b2+c2;
		sum2 = H2_ORH2_coll_fit_par[init][final][4]+f2;
		k = 0.;
		if( fabs(sum1) < 38. )
			k += pow(10., sum1 );
		if( fabs(sum2) < 38. )
			k += pow(10., sum2 );

		
		{
			enum {DEBUG_LOC=false};
			if( DEBUG_LOC )
			{
				fprintf(ioQQQ,"DEBUG H2 H2ortho rate %.3e \n",
					k);
			}
		}
	}
	/*unknown invalid entries returns '99'*/
	else
	{
		k = 99;
	}
	return k;
}

/*>>chng 08 feb 27, GS*/
double H2_PAH2_coll(int init, int final, double temp)
{

	double k, t, b2, c2, f2, t2;

	DEBUG_ENTRY( "H2_PAH2_coll()" );

	/* Invalid entries returns '-1':the initial indices are smaller than the final indices */
	if( temp<1 || temp > 1e4 )
	{
		k = -1;
	}
	else if( init <= final )
	{
		k = -1;
	}
	/* Invalid returns '-1': the indices are greater than 302 or smaller than 0 */
	else if( init < 0 || init >302 || final < 0 || final > 302 )
	{
		k = -1;
	}
	/* Undefined indices returns '0' */
	else if( !lgDefn_H2PAH2_coll[init][final] )
	{
		k = -1;
	}
	/* defined indices */
	else if( lgDefn_H2PAH2_coll[init][final] )
	{
		/* the fitting equation we used:
		* k_(vj,v'j') = 10^(a + b/(d*T/10^3+1) + c/t^2)
		* + 10^(e + f/(g*T/10^3+1)**h)
		* T in K, k in cm3/s. */
		double sum1 , sum2;

		
		t2 = temp/1e3; 
		/* t = T*10^-3 + 1*/
		t = t2+1; 
		b2 = H2_PAH2_coll_fit_par[init][final][1]/(H2_PAH2_coll_fit_par[init][final][3]*t2+1);
		c2 = H2_PAH2_coll_fit_par[init][final][2]/(t*t);
		{
			enum {DEBUG_LOC=false};
			if( DEBUG_LOC )
			{
				fprintf(ioQQQ,"bug H2 H2para coll\t%i %i %.3e %.3e %.3e \n",
					init,final,
					H2_PAH2_coll_fit_par[init][final][5],
					H2_PAH2_coll_fit_par[init][final][6]*t2+1, 
					H2_PAH2_coll_fit_par[init][final][7]
				/*pow(H2_PAH2_coll_fit_par[init][final][6]*t2+1,H2_PAH2_coll_fit_par[init][final][7])*/
				);
			}
		}
		/* this is log of f2 - see whether it is within bounds */
		sum1 = H2_PAH2_coll_fit_par[init][final][7] * log10( H2_PAH2_coll_fit_par[init][final][6]*t2+1. );
		/* this protects against overflow */
		if( fabs(sum1)< 38. )
		{
			
			f2 = H2_PAH2_coll_fit_par[init][final][5]/pow(10. , sum1 );
		}
		else
			f2= 0.;
		{
			enum {DEBUG_LOC=false };
			if( DEBUG_LOC )
			{
				fprintf(ioQQQ,"bug H2 H2para coll\t%i %i %.3e %.3e %.3e %.3e %.3e sum %.3e %.3e \n",
					init,final,
					H2_PAH2_coll_fit_par[init][final][0],
					b2, 
					c2,
					H2_PAH2_coll_fit_par[init][final][4],
					f2  ,
					H2_PAH2_coll_fit_par[init][final][0]+b2+c2,
					H2_PAH2_coll_fit_par[init][final][4]+f2);
			}
		}

		sum1 = H2_PAH2_coll_fit_par[init][final][0]+b2+c2;
		sum2 = H2_PAH2_coll_fit_par[init][final][4]+f2;
		k = 0.;
		if( fabs(sum1) < 38. )
			k += pow(10., sum1 );
		if( fabs(sum2) < 38. )
			k += pow(10., sum2 );

		
		{
			enum {DEBUG_LOC=false};
			if( DEBUG_LOC )
			{
				fprintf(ioQQQ,"DEBUG H2 H2para rate %.3e \n",
					k);
			}
		}
	}
	/*unknown invalid entries returns '99'*/
	else
	{
		k = 99;
	}
	return k;
}

