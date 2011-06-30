/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*H2_Create create variables for the H2 molecule, called by ContCreatePointers after continuum
 * mesh has been set up */
#include "cddefines.h" 
#include "physconst.h" 
#include "mole.h"
#include "taulines.h"
#include "lines_service.h"
#include "opacity.h" 
#include "hmi.h" 
#include "ipoint.h"
#include "grainvar.h"
#include "h2.h"
#include "h2_priv.h"

/* if this is set true then code will print energies and stop */
/*@-redef@*/
enum {DEBUG_ENER=false};
/*@+redef@*/

/* this is equation 8 of Takahashi 2001, clearer definition is given in
 * equation 5 and following discussion of
 * >>refer	H2	formation	Takahashi, J., & Uehara, H., 2001, ApJ, 561, 843-857
 * 0.27eV, convert into wavenumbers */
static double XVIB[H2_TOP] = { 0.70 , 0.60 , 0.20 };
static double Xdust[H2_TOP] = { 0.04 , 0.10 , 0.40 };

/* this is energy difference between bottom of potential well and 0,0
 * the Takahashi energy scale is from the bottom,
 * 2201.9 wavenumbers  */
static const double energy_off = 0.273*FREQ_1EV/SPEEDLIGHT;

STATIC double EH2_eval(  long int iVib , int ipH2 )
{
	double EH2_here;
	double Evm = H2_DissocEnergies[0]* XVIB[ipH2] + energy_off;

	double Ev = (energy_wn[0][iVib][0]+energy_off);
	/* equation 9 of Takahashi 2001 which is only an approximation
	double EH2 = H2_DissocEnergies[0] * (1. - Xdust[ipH2] ); */
	/* equation 1, 2 of 
	 * Takahashi, Junko, & Uehara, Hideya, 2001, ApJ, 561, 843-857,
	 * this is heat deposited on grain by H2 formation in this state */
	double Edust = H2_DissocEnergies[0] * Xdust[ipH2] *
		( 1. - ( (Ev - Evm) / (H2_DissocEnergies[0]+energy_off-Evm)) *
		( (1.-Xdust[ipH2])/2.) );
	ASSERT( Edust >= 0. );

	/* energy is total binding energy less energy lost on grain surface 
	 * and energy offset */
	EH2_here = H2_DissocEnergies[0]+energy_off - Edust;
	ASSERT( EH2_here >= 0.);

	return EH2_here;
}

/*H2_vib_dist evaluates the vibration distribution for H2 formed on grains */
STATIC double H2_vib_dist( long int iVib , int ipH2 , double EH2)
{
	double G1[H2_TOP] = { 0.3 , 0.4 , 0.9 };
	double G2[H2_TOP] = { 0.6 , 0.6 , 0.4 };
	double Evm = H2_DissocEnergies[0]* XVIB[ipH2] + energy_off;
	double Fv;
	if( (energy_wn[0][iVib][0]+energy_off) <= Evm )
	{
		/* equation 4 of Takahashi 2001 */
		Fv = sexp( POW2( (energy_wn[0][iVib][0]+energy_off - Evm)/(G1[ipH2]* Evm ) ) );
	}
	else
	{
		/* equation 5 of Takahashi 2001 */
		Fv = sexp( POW2( (energy_wn[0][iVib][0]+energy_off - Evm)/(G2[ipH2]*(EH2 - Evm ) ) ) );
	}
	return Fv;
}


/*H2_Create create variables for the H2 molecule, called by
 * ContCreatePointers after continuum mesh has been set up */
void H2_Create(void)
{
	long int i , iElecHi , iElecLo;
	long int iVibHi , iVibLo;
	long int iRotHi , iRotLo;
	long int iElec, iVib , iRot;
	long int nColl,
		nlines;
	int ier;
	int nEner;
	/* >>chng 03 nov 26, from enum H2_type to int */
	int ipH2;
	realnum sum , sumj , sumv , sumo , sump;

	/* this is flag set above - when true h2 code is not executed - this is way to
	 * avoid this code when it is not working */
	/* only malloc vectors one time per core load */
	if( lgH2_READ_DATA || !h2.lgH2ON )
		return;

	DEBUG_ENTRY( "H2_Create()" );

	/* print string if H2 debugging is enabled */
	if( mole.nH2_TRACE )
		fprintf(ioQQQ," H2_Create called in DEBUG mode.\n");

	/* this was option to print all electronic states in the main printout - but
	 * number of electronic states was not known at initialization so set to -1,
	 * will set properly now */
	if( h2.nElecLevelOutput < 1 )
		 h2.nElecLevelOutput = mole.n_h2_elec_states;

	/* this var is in h2.h and prevents h2 from being change once committed here */
	lgH2_READ_DATA = true;

	/* create special vector that saves collision rates within ground */
	/* this will contain a vector for collisions within the X ground electronic state,
	 * CollRateFit[vib_up][rot_up][vib_lo][rot_lo][coll_type][3] */
	/* N_X_COLLIDER is number of different species that collide */
	iElecHi = 0;
	/* the current data set is limited to vib hi <= 3 */
	/* will define collision rates for all possible transitions within X */
	CollRateFit.reserve(h2.nVib_hi[iElecHi]+1);
	H2_CollRate.reserve(h2.nVib_hi[iElecHi]+1);
	for( iVibHi = 0; iVibHi <= h2.nVib_hi[iElecHi]; ++iVibHi )
	{
		CollRateFit.reserve(iVibHi,h2.nRot_hi[iElecHi][iVibHi]+1);
		H2_CollRate.reserve(iVibHi,h2.nRot_hi[iElecHi][iVibHi]+1);
		for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
		{
			CollRateFit.reserve(iVibHi,iRotHi,h2.nVib_hi[iElecHi]+1);
			H2_CollRate.reserve(iVibHi,iRotHi,h2.nVib_hi[iElecHi]+1);
			for( iVibLo=0; iVibLo<(h2.nVib_hi[iElecHi]+1); ++iVibLo )
			{
				CollRateFit.reserve(iVibHi,iRotHi,iVibLo,h2.nRot_hi[iElecHi][iVibLo]+1);
				H2_CollRate.reserve(iVibHi,iRotHi,iVibLo,h2.nRot_hi[iElecHi][iVibLo]+1);
				for( iRotLo=0; iRotLo<=h2.nRot_hi[iElecHi][iVibLo]; ++iRotLo )
				{
					H2_CollRate.reserve(iVibHi,iRotHi,iVibLo,iRotLo,N_X_COLLIDER);
					CollRateFit.reserve(iVibHi,iRotHi,iVibLo,iRotLo,N_X_COLLIDER);
					for( nColl=0; nColl<N_X_COLLIDER; ++nColl )
					{
						/* the last one - the three coefficients */
						CollRateFit.reserve(iVibHi,iRotHi,iVibLo,iRotLo,nColl,3);
					}
				}
			}
		}
	}

	CollRateFit.alloc();
	H2_CollRate.alloc();

	/* zero out the collisional rates since only a minority of them are known*/
	CollRateFit.zero();
	H2_CollRate.zero();

	/* create space for the electronic levels */
	H2_populations.reserve(mole.n_h2_elec_states);
	pops_per_vib.reserve(mole.n_h2_elec_states);
	H2_dissprob.reserve(mole.n_h2_elec_states);

	for( iElec = 0; iElec<mole.n_h2_elec_states; ++iElec )
	{

		if( mole.nH2_TRACE  >= mole.nH2_trace_full)
			fprintf(ioQQQ,"elec %li highest vib= %li\n", iElec , h2.nVib_hi[iElec] );

		ASSERT( h2.nVib_hi[iElec] > 0 );

		/* h2.nVib_hi is now the highest vibrational level before dissociation,
		 * now allocate space to hold the number of rotation levels */
		H2_populations.reserve(iElec,h2.nVib_hi[iElec]+1);
		pops_per_vib.reserve(iElec,h2.nVib_hi[iElec]+1);
		if( iElec > 0 )
			H2_dissprob.reserve(iElec,h2.nVib_hi[iElec]+1);

		/* now loop over all vibrational levels, and find out how many rotation levels there are */
		/* ground is special since use tabulated data - there are 14 vib states,
		 * ivib=14 is highest */
		for( iVib = 0; iVib <= h2.nVib_hi[iElec]; ++iVib )
		{
			/* lastly create the space for the rotation quantum number */
			H2_populations.reserve(iElec,iVib,h2.nRot_hi[iElec][iVib]+1);
			if( iElec > 0 )
				H2_dissprob.reserve(iElec,iVib,h2.nRot_hi[iElec][iVib]+1);
		}
	}

	H2_populations.alloc();
	H2_populations_LTE.alloc( H2_populations.clone() );
	H2_old_populations.alloc( H2_populations.clone() );
	H2_Boltzmann.alloc( H2_populations.clone() );
	H2_stat.alloc( H2_populations.clone() );
	energy_wn.alloc( H2_populations.clone() );
	H2_rad_rate_out.alloc( H2_populations.clone() );
	H2_lgOrtho.alloc( H2_populations.clone() );

	pops_per_vib.alloc();

	H2_dissprob.alloc();
	H2_disske.alloc( H2_dissprob.clone() );

	/* set this one time, will never be set again, but might be printed */
	H2_rad_rate_out.zero();

	/* these do not have electronic levels - all within X */
	H2_ipPhoto.reserve(h2.nVib_hi[0]+1);

	/* space for the vibration levels */
	for( iVib = 0; iVib <= h2.nVib_hi[0]; ++iVib )
	{
		/* space for the rotation quantum number */
		H2_ipPhoto.reserve(iVib,h2.nRot_hi[0][iVib]+1);
	}

	H2_ipPhoto.alloc();
	H2_col_rate_in.alloc( H2_ipPhoto.clone() );
	H2_col_rate_out.alloc( H2_ipPhoto.clone() );
	H2_rad_rate_in.alloc( H2_ipPhoto.clone() );
	H2_coll_dissoc_rate_coef.alloc( H2_ipPhoto.clone() );
	H2_coll_dissoc_rate_coef_H2.alloc( H2_ipPhoto.clone() );
	H2_X_colden.alloc( H2_ipPhoto.clone() );
	H2_X_rate_from_elec_excited.alloc( H2_ipPhoto.clone() );
	H2_X_rate_to_elec_excited.alloc( H2_ipPhoto.clone() );
	H2_X_colden_LTE.alloc( H2_ipPhoto.clone() );
	H2_X_formation.alloc( H2_ipPhoto.clone() );
	H2_X_Hmin_back.alloc( H2_ipPhoto.clone() );

	for( iVib = 0; iVib <= h2.nVib_hi[0]; ++iVib )
	{
		for( iRot=h2.Jlowest[0]; iRot<=h2.nRot_hi[0][iVib]; ++iRot )
		{
			/* >>chng 04 jun 14, set these to bad numbers */
			H2_rad_rate_in[iVib][iRot] = -BIGFLOAT;
			H2_coll_dissoc_rate_coef[iVib][iRot] = -BIGFLOAT;
			H2_coll_dissoc_rate_coef_H2[iVib][iRot] = -BIGFLOAT;
		}
	}
	/* zero out the matrices */
	H2_X_colden.zero();
	H2_X_colden_LTE.zero();
	H2_X_formation.zero();
	H2_X_Hmin_back.zero();
	/* rates [cm-3 s-1] from elec excited states into X only vib and rot */
	H2_X_rate_from_elec_excited.zero();
	/* rates [s-1] to elec excited states from X only vib and rot */
	H2_X_rate_to_elec_excited.zero();

	/* distribution function for populations following formation from H minus H- */
	H2_X_hminus_formation_distribution.reserve(nTE_HMINUS);
	for( i=0; i<nTE_HMINUS; ++i )
	{
		H2_X_hminus_formation_distribution.reserve(i,h2.nVib_hi[0]+1);
		/* space for the vibration levels */
		for( iVib = 0; iVib <= h2.nVib_hi[0]; ++iVib )
		{
			H2_X_hminus_formation_distribution.reserve(i,iVib,h2.nRot_hi[0][iVib]+1);
		}
	}
	H2_X_hminus_formation_distribution.alloc();
	H2_X_hminus_formation_distribution.zero();
	H2_Read_hminus_distribution();

	/* >>chng 05 jun 20, do not use this, which is highly processed - use ab initio
	 * rates of excitation to electronic levels instead */
	/* read in cosmic ray distribution information
	H2_Read_Cosmicray_distribution(); */

	/* grain formation matrix */
	H2_X_grain_formation_distribution.reserve(H2_TOP);
	for( ipH2=0; ipH2<(int)H2_TOP; ++ipH2 )
	{
		H2_X_grain_formation_distribution.reserve(ipH2,h2.nVib_hi[0]+1);

		/* space for the vibration levels */
		for( iVib = 0; iVib <= h2.nVib_hi[0]; ++iVib )
		{
			H2_X_grain_formation_distribution.reserve(ipH2,iVib,h2.nRot_hi[0][iVib]+1);
		}
	}
	H2_X_grain_formation_distribution.alloc();
	H2_X_grain_formation_distribution.zero();

	/* space for the energy vector is now malloced, must fill it in,
	 * defines array energy_wn[nelec][iVib][iRot] */
	for( iElec=0; iElec<mole.n_h2_elec_states; ++iElec )
	{
		/* get energies out of files into array energy_wn[nelec][iVib][iRot] */
		H2_ReadEnergies(iElec);

		/* get dissociation probabilities and energies - ground state is stable */
		if( iElec > 0 )
			H2_ReadDissprob(iElec);
	}

	/* >>02 oct 18, add photodissociation, H2 + hnu => 2H + KE */
	/* we now have ro-vib energies, now set up threshold array offsets
	 * for photodissociation */
	for( iVib = 0; iVib <= h2.nVib_hi[0]; ++iVib )
	{
		for( iRot=h2.Jlowest[0]; iRot<=h2.nRot_hi[0][iVib]; ++iRot )
		{
			/* this is energy needed to get up to n=3 electronic continuum 
			 * H2 cannot dissociate following absorption of a continuum photon into the
			 * continuum above X, which would require little energy, and would be given
			 * by H2_DissocEnergies[0], because that process violates momentum conservation 
			 * these would be the triplet states - permitted are into singlets
			 * the effective full wavelength range of this process is from Lya to 
			 * Lyman limit in shielded regions 
			 * tests show limits are between 850A and 1220A - so Lya is included */
			/** \todo	1	add this as a Lya excitation process */
			/*>>KEYWORD	Allison & Dalgarno; continuum dissociation; */
			double thresh = (H2_DissocEnergies[1] - energy_wn[0][iVib][iRot])*WAVNRYD;
			/*fprintf(ioQQQ,"DEBUG\t%.2f\t%f\n", RYDLAM/thresh , thresh);*/
			/* in theory we should be able to assert that thesh just barely reaches
			 * lya, but actual numbers reach down to 0.749 ryd */
			ASSERT( thresh > 0.74 );
			H2_ipPhoto[iVib][iRot] = ipoint(thresh);
		}
	}

	nH2_energies = 0;
	for( iElec=0; iElec<mole.n_h2_elec_states; ++iElec)
	{
		/* the number of levels within the molecule */
		nH2_energies += nLevels_per_elec[iElec];
	}

	if( mole.nH2_TRACE >= mole.nH2_trace_full ) 
	{
		for( iElec=0; iElec<mole.n_h2_elec_states; ++iElec)
		{
			/* print the number of levels within iElec */
			fprintf(ioQQQ,"\t(%li %li)", iElec ,  nLevels_per_elec[iElec] );
		}
		fprintf(ioQQQ,
			" H2_Create: there are %li electronic levels, in each level there are",
			mole.n_h2_elec_states);
		fprintf(ioQQQ,
			" for a total of %li levels.\n", nH2_energies );
	}

	/* now read in the various sets of collision data */
	for( nColl=0; nColl<N_X_COLLIDER; ++nColl )
	{
		/* ground state has tabulated data */
		H2_CollidRateRead(nColl);
	}
	/* option to add gaussian random mole */
	if( mole.lgH2_NOISE )
	{
		iElecHi = 0;
		/* loop over all transitions */
		for( iVibHi = 0; iVibHi <= VIB_COLLID; ++iVibHi )
		{
			for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
			{
				for( iVibLo=0; iVibLo<(VIB_COLLID+1); ++iVibLo )
				{
					for( iRotLo=0; iRotLo<=h2.nRot_hi[iElecHi][iVibLo]; ++iRotLo )
					{
						/* first set of expressions are series that adds to log of rate,
						 * so we will add the gaussian noise */
						/* >>chng 05 dec 13, GS, last two fits are different,
						 * loop had been to N_X_COLLIDER-1 and so included Stancil data
						 * with noise became negative */
						for( nColl=0; nColl<N_X_COLLIDER-2; ++nColl )
						{
							/* the gaussian random number, many possible collision rates
							 * have no data, and CollRateFit[][][][][][0] is zero - do not
							 * scramble these, only scramble the non-zero rates */
							if( CollRateFit[iVibHi][iRotHi][iVibLo][iRotLo][nColl][0] != 0. )
							{
								/* this returns the log of the random noise */
								realnum r = (realnum)RandGauss( mole.xMeanNoise , mole.xSTDNoise );
								/* check that coefficient 0 is the one we want to hit with the mole */
								/* these are used at line 2990 below, */
								CollRateFit[iVibHi][iRotHi][iVibLo][iRotLo][nColl][0] += r;
							}
						}
						/* >>chng 04 feb 19, break out last one which is linear and must be treated separately */
						/* for late one is linear so use pow */
						if( CollRateFit[iVibHi][iRotHi][iVibLo][iRotLo][N_X_COLLIDER-2][0] != 0. )
						{
							/* this returns the log of the random noise */
							realnum r = (realnum)RandGauss( mole.xMeanNoise , mole.xSTDNoise );
							/* check that coefficient 0 is the one we want to hit with the mole */
							/* these are used at line 2990 below, */
							CollRateFit[iVibHi][iRotHi][iVibLo][iRotLo][N_X_COLLIDER-2][0] *= pow((realnum)10.f,r);
						}
					}
				}
			}
		}
	}


	/* create arrays for energy sorted referencing of e, v, J */
	H2_energies = (realnum*)MALLOC(sizeof(realnum)*(unsigned)nH2_energies );
	H2_ipX_ener_sort = (long int*)MALLOC(sizeof(long int)*(unsigned)nH2_energies );
	ipElec_H2_energy_sort = (long int*)MALLOC(sizeof(long int)*(unsigned)nH2_energies );
	ipVib_H2_energy_sort = (long int*)MALLOC(sizeof(long int)*(unsigned)nH2_energies );
	ipRot_H2_energy_sort = (long int*)MALLOC(sizeof(long int)*(unsigned)nH2_energies );

	/* this will be total collision rate from an upper to a lower level within X */
	H2_X_source = (realnum*)MALLOC(sizeof(realnum)*(unsigned)nLevels_per_elec[0] );
	H2_X_sink = (realnum*)MALLOC(sizeof(realnum)*(unsigned)nLevels_per_elec[0] );

	H2_X_coll_rate.reserve(nLevels_per_elec[0]);
	/* now expand out to include all lower levels as lower state */
	for( i=1; i<nLevels_per_elec[0]; ++i )
	{
		H2_X_coll_rate.reserve(i,i);
	}
	H2_X_coll_rate.alloc();

	/* create a vector of sorted energies for X */
	ipEnergySort.reserve(mole.n_h2_elec_states);
	for( iElecHi=0; iElecHi<mole.n_h2_elec_states; ++iElecHi )
	{
		ipEnergySort.reserve(iElecHi,h2.nVib_hi[iElecHi]+1);
		for( iVibHi = 0; iVibHi <= h2.nVib_hi[iElecHi]; ++iVibHi )
		{
			ipEnergySort.reserve(iElecHi,iVibHi,h2.nRot_hi[iElecHi][iVibHi]+1);
		}
	}
	ipEnergySort.alloc();

	nEner = 0;
	for( iElecHi=0; iElecHi<mole.n_h2_elec_states; ++iElecHi )
	{
		/* get set of energies */
		for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
		{
			for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
			{
				H2_energies[nEner] = (realnum)energy_wn[iElecHi][iVibHi][iRotHi];
				ipElec_H2_energy_sort[nEner] = iElecHi;
				ipVib_H2_energy_sort[nEner] = iVibHi;
				ipRot_H2_energy_sort[nEner] = iRotHi;
				ipEnergySort[iElecHi][iVibHi][iRotHi] = -1;
				++nEner;
			}
		}
	}

	ASSERT( nH2_energies == nEner );

	/* sort the energy levels so that we can do top-down trickle of states */
	/*spsort netlib routine to sort array returning sorted indices */
	spsort(
		/* input array to be sorted */
		H2_energies, 
		/* number of values in the molecule */
		nH2_energies, 
		/* permutation output array */
		H2_ipX_ener_sort, 
		/* flag saying what to do - 1 sorts into increasing order, not changing
		* the original vector, -1 sorts into decreasing order. 2, -2 change vector */
		1, 
		/* error condition, should be 0 */
		&ier);

	/* now loop over the energies confirming the order */
	for( nEner=0; nEner<nH2_energies; ++nEner )
	{
		if( nEner+1 < nLevels_per_elec[0] )
			ASSERT( H2_energies[H2_ipX_ener_sort[nEner]] < 
			H2_energies[H2_ipX_ener_sort[nEner+1]] );
		/* following will print quantum indices and energies */
		/*fprintf(ioQQQ,"%li\t%li\t%.3e\n",
			ipVib_H2_energy_sort[H2_ipX_ener_sort[nEner]],
			ipRot_H2_energy_sort[H2_ipX_ener_sort[nEner]],
			H2_energies[H2_ipX_ener_sort[nEner]]);*/
		i = H2_ipX_ener_sort[nEner];
		iElec = ipElec_H2_energy_sort[i];
		iRot = ipRot_H2_energy_sort[i];
		iVib = ipVib_H2_energy_sort[i];
		/* this allows v,J to map into energy sorted array */
		ipEnergySort[iElec][iVib][iRot] = nEner;
	}

	/* now find number of levels in H2g */
	for( nEner=0; nEner<nLevels_per_elec[0]; ++nEner )
	{
		i = H2_ipX_ener_sort[nEner];
		iRot = ipRot_H2_energy_sort[i];
		iVib = ipVib_H2_energy_sort[i];
		if( energy_wn[0][iVib][iRot] > ENERGY_H2_STAR	)
			break;
		nEner_H2_ground = nEner;
	}
	/* need to increment it so that this is the number of levels, not the index
	 * of the highest level */
	++nEner_H2_ground;

	/* this is the number of levels to do with the matrix - set with the
	 * atom h2 matrix command, keyword ALL means to do all of X in the matrix
	 * but number of levels within X was not known when the command was parsed,
	 * so this was set to -1 to defer setting to all until now */
	if( nXLevelsMatrix<0 )
	{
		nXLevelsMatrix = nLevels_per_elec[0];
	}
	else if( nXLevelsMatrix > nLevels_per_elec[0] )
	{
		fprintf( ioQQQ, 
			" The total number of levels used in the matrix solver was set to %li but there are only %li levels in X.\n Sorry.\n",
			nXLevelsMatrix ,
			nLevels_per_elec[0]);
		cdEXIT(EXIT_FAILURE);
	}

	/* create the main array of lines */
	H2Lines.reserve(mole.n_h2_elec_states);

	nlines = 0;
	for( iElecHi=0; iElecHi<mole.n_h2_elec_states; ++iElecHi )
	{
		H2Lines.reserve(iElecHi,h2.nVib_hi[iElecHi]+1);
		for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
		{
			H2Lines.reserve(iElecHi,iVibHi,h2.nRot_hi[iElecHi][iVibHi]+1);
			for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
			{
				/* now the lower levels */
				/* NB - X is the only lower level considered here, since we are only 
				 * concerned with excited electronic levels as a photodissociation process
				 * code exists to relax this assumption - simply change following to iElecHi */
				long int lim_elec_lo = 0;
				H2Lines.reserve(iElecHi,iVibHi,iRotHi,1);
				for( iElecLo=0; iElecLo<=lim_elec_lo; ++iElecLo )
				{
					/* want to include all vib states in lower level if 
					 * different elec level, but only lower vib levels if 
					 * same elec level */
					long int nv = h2.nVib_hi[iElecLo];
					/* within X, no transitions v_hi < v_lo transitions exist */ 
					if( iElecLo==iElecHi )
						nv = iVibHi;
					H2Lines.reserve(iElecHi,iVibHi,iRotHi,iElecLo,nv+1);
					for( iVibLo=0; iVibLo<=nv; ++iVibLo )
					{
						long nr = h2.nRot_hi[iElecLo][iVibLo];
						if( iElecLo==iElecHi && iVibHi==iVibLo )
							/* max because cannot malloc 0 bytes */
							nr = MAX2(1,iRotHi-1);
						H2Lines.reserve(iElecHi,iVibHi,iRotHi,iElecLo,iVibLo,nr+1);
						nlines += nr+1;
					}
				}
			}
		}
	}

	H2Lines.alloc();
	H2_SaveLine.alloc( H2Lines.clone() );
	lgH2_line_exists.alloc( H2Lines.clone() );

	/* zero out array used to save emission line intensities */
	H2_SaveLine.zero();
	/* set flag saying that space exists */
	lgH2_line_exists.zero();

	if( mole.nH2_TRACE >= mole.nH2_trace_full )
		fprintf(ioQQQ," There are a total of %li lines in the entire H2 molecule.\n", nlines );

	/* junk the transitions */
	for( iElecHi=0; iElecHi<mole.n_h2_elec_states; ++iElecHi )
	{
		for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
		{
			for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
			{
				/* NB - X is the only lower level considered here, since we are only 
				 * concerned with excited electronic levels as a photodissociation process
				 * code exists to relax this assumption - simply change following to iElecHi */
				long int lim_elec_lo = 0;
				for( iElecLo=0; iElecLo<=lim_elec_lo; ++iElecLo )
				{
					/* want to include all vib states in lower level if different elec level,
					 * but only lower vib levels if same elec level */
					long int nv = h2.nVib_hi[iElecLo];
					if( iElecLo==iElecHi )
						nv = iVibHi;
					for( iVibLo=0; iVibLo<=nv; ++iVibLo )
					{
						long nr = h2.nRot_hi[iElecLo][iVibLo];
						if( iElecLo==iElecHi && iVibHi==iVibLo )
							nr = iRotHi-1;
						for( iRotLo=h2.Jlowest[iElecLo]; iRotLo<=nr; ++iRotLo )
						{
							H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Junk();
						}
					}
				}
			}
		}
	}

	/* now set up state pointers and zero out the transitions */
	for( iElecHi=0; iElecHi<mole.n_h2_elec_states; ++iElecHi )
	{
		for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
		{
			for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
			{
				H2Lines[iElecHi][iVibHi][iRotHi][0][0][0].Hi = AddState2Stack();

				/* NB - X is the only lower level considered here, since we are only 
				 * concerned with excited electronic levels as a photodissociation process
				 * code exists to relax this assumption - simply change following to iElecHi */
				long int lim_elec_lo = 0;
				for( iElecLo=0; iElecLo<=lim_elec_lo; ++iElecLo )
				{
					/* want to include all vib states in lower level if different elec level,
					 * but only lower vib levels if same elec level */
					long int nv = h2.nVib_hi[iElecLo];
					if( iElecLo==iElecHi )
						nv = iVibHi;
					for( iVibLo=0; iVibLo<=nv; ++iVibLo )
					{
						long nr = h2.nRot_hi[iElecLo][iVibLo];
						if( iElecLo==iElecHi && iVibHi==iVibLo )
							nr = iRotHi-1;
						for( iRotLo=h2.Jlowest[iElecLo]; iRotLo<=nr; ++iRotLo )
						{
							H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Hi = 
								H2Lines[iElecHi][iVibHi][iRotHi][0][0][0].Hi;
							/* lower level is higher level of a previous transition. */
							H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Lo = 
								H2Lines[iElecLo][iVibLo][iRotLo][0][0][0].Hi;

							/* set initial values for each line structure */
							H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Zero();
						}
					}
				}
			}
		}
	}

	/* space for the energy vector is now malloced, must read trans probs from table */
	for( iElec=0; iElec<mole.n_h2_elec_states; ++iElec )
	{
		/* ground state has tabulated data */
		H2_ReadTransprob(iElec);
	}

	/* set all statistical weights - ours is total statistical weight - 
	 * including nuclear spin */
	for( iElecHi=0; iElecHi<mole.n_h2_elec_states; ++iElecHi )
	{
		for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
		{
			for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
			{
				/* unlike atoms, for H2 nuclear spin is taken into account - so the
				 * statistical weight of even and odd J states differ by factor of 3 - see page 166, sec par
				 * >>>refer	H2	H2_stat wght	Shull, J.M., & Beckwith, S., 1982, ARAA, 20, 163-188 */
				if( is_odd(iRotHi+H2_nRot_add_ortho_para[iElecHi]) )
				{
					/* ortho */
					H2_lgOrtho[iElecHi][iVibHi][iRotHi] = true;
					H2_stat[iElecHi][iVibHi][iRotHi] = 3.f*(2.f*iRotHi+1.f);
				}
				else
				{
					/* para */
					H2_lgOrtho[iElecHi][iVibHi][iRotHi] = false;
					H2_stat[iElecHi][iVibHi][iRotHi] = (2.f*iRotHi+1.f);
				}
			}
		}
	}

	/* set up transition parameters */
	for( iElecHi=0; iElecHi<mole.n_h2_elec_states; ++iElecHi )
	{
		for( iVibHi=0; iVibHi<=h2.nVib_hi[iElecHi]; ++iVibHi )
		{
			for( iRotHi=h2.Jlowest[iElecHi]; iRotHi<=h2.nRot_hi[iElecHi][iVibHi]; ++iRotHi )
			{
				/* now the lower levels */
				/* NB - X is the only lower level considered here, since we are only 
				 * concerned with excited electronic levels as a photodissociation process
				 * code exists to relax this assumption - simply change following to iElecHi */
				long int lim_elec_lo = 0;
				for( iElecLo=0; iElecLo<=lim_elec_lo; ++iElecLo )
				{
					/* want to include all vib states in lower level if different elec level,
					 * but only lower vib levels if same elec level */
					long int nv = h2.nVib_hi[iElecLo];
					if( iElecLo==iElecHi )
						nv = iVibHi;
					for( iVibLo=0; iVibLo<=nv; ++iVibLo )
					{
						long nr = h2.nRot_hi[iElecLo][iVibLo];

						if( iElecLo==iElecHi && iVibHi==iVibLo )
							nr = iRotHi-1;
						for( iRotLo=h2.Jlowest[iElecLo]; iRotLo<=nr; ++iRotLo )
						{
							/* NB this must be kept parallel with nelem and ionstag in H2Lines transition struc,
							* since that struc expects to find the abundances here - abund set in hmole.c */
							H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Hi->nelem = LIMELM+3;
							/* this does not mean anything for a molecule */
							H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Hi->IonStg = 1;

							/* statistical weights of lower and upper levels */
							H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Lo->g = H2_stat[iElecLo][iVibLo][iRotLo];
							H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Hi->g = H2_stat[iElecHi][iVibHi][iRotHi];

							/* energy of the transition in wavenumbers */
							H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].EnergyWN = 
								(realnum)(energy_wn[iElecHi][iVibHi][iRotHi] - energy_wn[iElecLo][iVibLo][iRotLo]);

							/*wavelength of transition in Angstroms */
							if( H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].EnergyWN > SMALLFLOAT)
								H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].WLAng = 
								(realnum)(1.e8f/H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].EnergyWN /
								RefIndex( H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].EnergyWN));

							H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].EnergyK = 
								(realnum)(T1CM)*H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].EnergyWN;

							/* energy of photon in ergs */
							H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].EnergyErg = 
								(realnum)(ERG1CM)*H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].EnergyWN;

							/* only do this if radiative transition exists */
							if( lgH2_line_exists[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo] )
							{
								/* line redistribution function - will use complete redistribution */
								/* >>chng 04 mar 26, should include damping wings, especially for electronic
								 * transitions, had used doppler core only */
								H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->iRedisFun = ipCRDW;

								/* line optical depths in direction towards source of ionizing radiation */
								H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->TauIn = opac.taumin;
								H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->TauCon = opac.taumin;
								/* outward optical depth */
								H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->TauTot = 1e20f;


								H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->dampXvel = 
									(realnum)(H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->Aul/
									H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].EnergyWN/PI4);

								H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->gf = 
									(realnum)(GetGF(H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->Aul,
									H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].EnergyWN,
									H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Hi->g ) );

								/* derive the absorption coefficient, call to function is gf, wl (A), g_low */
								H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->opacity = (realnum)(
									abscf(H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->gf,
									H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].EnergyWN,
									H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Lo->g));

								if( iElecHi > 0 )
								{
									/* cosmic ray and non-thermal suprathermal excitation 
									 * to singlet state of H2 (B,C,B',D)
									 * cross section is equation 5 of 
									 *>>refer	H2	cs	Liu, W. & Dalgarno, A. 1994, ApJ, 428, 769
									 * relative to H I Lya cross section 
									 * this is used in mole_h2.cpp to derive H2 electronic excitations
									 * from the H I Lya rate
									 * the following is dimensionless scale factor for excitation 
									 * relative to H I Lya
									 */
									H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Coll.col_str = (realnum)(
										pow3(H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].WLAng*1e-8)*
										(H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Hi->g/
										H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Lo->g)*
										H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Emis->Aul*
										log(3.)*HPLANCK/(160.f*pow3(PI)*0.5*1e-8*EN1EV)/6.e-17);
									ASSERT(H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Coll.col_str>0.);
								}
								else
								{
									/* excitation within X - not treated this way */
									H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Coll.col_str = 0.;
								}
							}
							else
							{
								/* Aul is zero but cosmic ray collisions are not
								 * zero in this case - this is the Aul = 0 branch */
								H2Lines[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo].Coll.col_str = 0.;
							}
						}
					}
				}
			}
		}
	}

	/* define branching ratios for deposition of H2 formed on grain surfaces,
	 * set true to use Takahashi distribution, false to use Draine & Bertoldi */

	/* loop over all types of grain surfaces */
	/* >>chng 02 oct 08, resolved grain types */
	/* number of different grain types H2_TOP is set in grainvar.h,
	 * types are ice, silicate, graphite */
	for( ipH2=0; ipH2<(int)H2_TOP; ++ipH2 )
	{
		sum = 0.;
		sumj = 0.;
		sumv = 0.;
		sumo = 0.;
		sump = 0.;
		iElec = 0;
		/* first is Draine distribution */
		if( hmi.chGrainFormPump == 'D' )
		{
			/* H2 formation temperature, for equation 19, page 271, of
			* >>refer	H2	formation distribution	Draine, B.T., & Bertoldi, F., 1996, ApJ, 468, 269-289
			*/
			double T_H2_FORM = 50000.;
			for( iVib = 0; iVib <= h2.nVib_hi[0]; ++iVib )
			{
				for( iRot=h2.Jlowest[0]; iRot<=h2.nRot_hi[0][iVib]; ++iRot )
				{
					/* no distinction between grain surface composition */
					H2_X_grain_formation_distribution[ipH2][iVib][iRot] = 
						/* first term is nuclear H2_stat weight */
						(1.f+2.f*H2_lgOrtho[iElec][iVib][iRot]) * (1.f+iVib) *
						(realnum)sexp( energy_wn[iElec][iVib][iRot]*T1CM/T_H2_FORM );
					sum += H2_X_grain_formation_distribution[ipH2][iVib][iRot];
					sumj += iRot * H2_X_grain_formation_distribution[ipH2][iVib][iRot];
					sumv += iVib * H2_X_grain_formation_distribution[ipH2][iVib][iRot];
					if( H2_lgOrtho[iElec][iVib][iRot] )
					{
						sumo += H2_X_grain_formation_distribution[ipH2][iVib][iRot];
					}
					else
					{
						/* >>chng 02 nov 14, [0][iVib][iRot] -> [ipH2][iVib][iRot], PvH */
						sump += H2_X_grain_formation_distribution[ipH2][iVib][iRot];
					}
				}
			}
		}
		else if( hmi.chGrainFormPump == 'T' )
		{
			/* Takahashi 2001 distribution */
			double Xrot[H2_TOP] = { 0.14 , 0.15 , 0.15 };
			double Xtrans[H2_TOP] = { 0.12 , 0.15 , 0.25 };
			/* first normalize the vibration distribution function */
			double sumvib = 0.;
			double EH2;

			for( iVib = 0; iVib <= h2.nVib_hi[0]; ++iVib )
			{
				double vibdist;
				EH2 = EH2_eval( iVib , ipH2 );
				vibdist = H2_vib_dist( iVib , ipH2 , EH2);
				sumvib += vibdist;
			}
			/* this branch, use distribution function from
			* >>refer	grain	physics	Takahashi, Junko, 2001, ApJ, 561, 254-263 */
			for( iVib = 0; iVib <= h2.nVib_hi[0]; ++iVib )
			{
				double Ev = (energy_wn[iElec][iVib][0]+energy_off);
				double Fv;
				/* equation 10 of Takahashi 2001, extra term is energy offset between bottom of potential
				 * the 0,0 level */
				double Erot;
				/*fprintf(ioQQQ," Evvvv\t%i\t%li\t%.3e\n", ipH2 ,iVib , Ev*WAVNRYD*EVRYD);*/

				EH2 = EH2_eval( iVib , ipH2 );

				/* equation 3 of Taktahashi & Uehara */
				Erot = (EH2 - Ev) * Xrot[ipH2] / (Xrot[ipH2] + Xtrans[ipH2]);

				/* email exchange with Junko Takahashi - 
				Thank you for your E-mail.
				I did not intend to generate negative Erot.
				I cut off the populations if their energy levels are negative, and made the total
				population be unity by using normalization factors (see, e.g., Eq. 12).

				I hope that my answer is of help to you and your work is going well.
				With best wishes,
				Junko

				>Thanks for the reply.  By cutting off the population, should we set the
				>population to zero when Erot becomes negative, or should we set Erot to
				>a small positive number? 

				I just set the population to zero when Erot becomes negative.
				Our model is still a rough one for the vibration-rotation distribution function
				of H2 newly formed on dust, because we have not yet had any exact
				experimental or theoretical data about it.
				With best wishes,
				Junko

				 */

				if( Erot > 0. )
				{
					/* the vibrational distribution */
					Fv = H2_vib_dist( iVib , ipH2 , EH2) / sumvib;
					/*fprintf(ioQQQ," vibbb\t%li\t%.3e\n", iVib , Fv );*/

					for( iRot=h2.Jlowest[0]; iRot<=h2.nRot_hi[0][iVib]; ++iRot )
					{
						/* equation 6 of Takahashi 2001 */
						double gaussian = 
							sexp( POW2( (energy_wn[iElec][iVib][iRot] - energy_wn[iElec][iVib][0] - Erot) /
							(0.5 * Erot) ) );
						/* equation 7 of Takahashi 2001 */
						double thermal_dist = 
							sexp( (energy_wn[iElec][iVib][iRot] - energy_wn[iElec][iVib][0]) /
							Erot );

						/* take the mean of the two */
						double aver = ( gaussian + thermal_dist ) / 2.;
						/*fprintf(ioQQQ,"rottt\t%i\t%li\t%li\t%.3e\t%.3e\t%.3e\t%.3e\n",
							ipH2,iVib,iRot,
							(energy_wn[iElec][iVib][iRot]+energy_off)*WAVNRYD*EVRYD,
							gaussian, thermal_dist , aver );*/

						/* thermal_dist does become > 1 since Erot can become negative */
						ASSERT( gaussian <= 1. /*&& thermal_dist <= 10.*/ );

						H2_X_grain_formation_distribution[ipH2][iVib][iRot] = (realnum)(
							/* first term is nuclear H2_stat weight */
							(1.f+2.f*H2_lgOrtho[iElec][iVib][iRot]) * Fv * (2.*iRot+1.) * aver );

						sum += H2_X_grain_formation_distribution[ipH2][iVib][iRot];
						sumj += iRot * H2_X_grain_formation_distribution[ipH2][iVib][iRot];
						sumv += iVib * H2_X_grain_formation_distribution[ipH2][iVib][iRot];
						if( H2_lgOrtho[iElec][iVib][iRot] )
						{
							sumo += H2_X_grain_formation_distribution[ipH2][iVib][iRot];
						}
						else
						{
							sump += H2_X_grain_formation_distribution[ipH2][iVib][iRot];
						}

					}
				}
				else
				{
					/* this branch Erot is non-positive, so no distribution */
					for( iRot=h2.Jlowest[0]; iRot<=h2.nRot_hi[0][iVib]; ++iRot )
					{
						H2_X_grain_formation_distribution[ipH2][iVib][iRot] = 0.;
					}
				}
			}
		}
		else if( hmi.chGrainFormPump == 't' )
		{ 
			/* thermal distribution at 1.5 eV, as suggested by Amiel & Jaques */
			/* thermal distribution, upper right column of page 239 of
			 *>>refer	H2	formation	Le Bourlot, J, 1991, A&A, 242, 235 
			 * set with command
			 * set h2 grain formation pumping thermal */
			double T_H2_FORM = 17329.;
			for( iVib = 0; iVib <= h2.nVib_hi[0]; ++iVib )
			{
				for( iRot=h2.Jlowest[0]; iRot<=h2.nRot_hi[0][iVib]; ++iRot )
				{
					/* no distinction between grain surface composition */
					H2_X_grain_formation_distribution[ipH2][iVib][iRot] = 
						/* first term is nuclear H2_stat weight */
						H2_stat[0][iVib][iRot] *
						(realnum)sexp( energy_wn[0][iVib][iRot]*T1CM/T_H2_FORM );
					sum += H2_X_grain_formation_distribution[ipH2][iVib][iRot];
					sumj += iRot * H2_X_grain_formation_distribution[ipH2][iVib][iRot];
					sumv += iVib * H2_X_grain_formation_distribution[ipH2][iVib][iRot];
					if( H2_lgOrtho[iElec][iVib][iRot] )
					{
						sumo += H2_X_grain_formation_distribution[ipH2][iVib][iRot];
					}
					else
					{
						/* >>chng 02 nov 14, [0][iVib][iRot] -> [ipH2][iVib][iRot], PvH */
						sump += H2_X_grain_formation_distribution[ipH2][iVib][iRot];
					}
				}
			}
		}
		else
			TotalInsanity();

		if( mole.nH2_TRACE >= mole.nH2_trace_full )
			fprintf(ioQQQ, "H2 form grains mean J= %.3f mean v = %.3f ortho/para= %.3f\n", 
				sumj/sum , sumv/sum , sumo/sump );

		iElec = 0;
		/* now rescale so that integral is unity */
		for( iVib = 0; iVib <= h2.nVib_hi[0]; ++iVib )
		{
			for( iRot=h2.Jlowest[0]; iRot<=h2.nRot_hi[0][iVib]; ++iRot )
			{
				H2_X_grain_formation_distribution[ipH2][iVib][iRot] /= sum;
				/* print the distribution function */
				/*if( energy_wn[iElec][iVib][iRot] < 5200. )
				fprintf(ioQQQ,"disttt\t%i\t%li\t%li\t%li\t%.4e\t%.4e\t%.4e\t%.4e\n",
					ipH2, iVib , iRot, (long)H2_stat[0][iVib][iRot] ,
					energy_wn[iElec][iVib][iRot] , 
					energy_wn[iElec][iVib][iRot]*T1CM , 
					H2_X_grain_formation_distribution[ipH2][iVib][iRot],
					H2_X_grain_formation_distribution[ipH2][iVib][iRot]/H2_stat[0][iVib][iRot]
					);*/
			}
		}
	}

	/* at this stage the full electronic, vibration, and rotation energies have been defined,
	 * this is an option to print the energies */
	{
		/* set following true to get printout, false to not print energies */
		if( DEBUG_ENER )
		{
			/* print title for quantum numbers and energies */
			/*fprintf(ioQQQ,"elec\tvib\trot\tenergy\n");*/
			for( iElec=0; iElec<mole.n_h2_elec_states; ++iElec )
			{
				/* now must specify the number of rotation levels within the vib levels */
				for( iVib=0; iVib<=h2.nVib_hi[iElec]; ++iVib )
				{
					for( iRot=0; iRot<=h2.nRot_hi[iElec][iVib]; ++iRot )
					{
						fprintf(ioQQQ,"%li\t%li\t%li\t%.5e\n",
							iElec, iVib, iRot ,
							energy_wn[iElec][iVib][iRot]);
					}
				}
			}
			/* this will exit the program after printing the level energies */
			cdEXIT(EXIT_SUCCESS);
		}
	}

	return;
}
