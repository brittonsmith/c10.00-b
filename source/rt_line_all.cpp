/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*RT_line_all do escape and destruction probabilities for all lines in code. */
#include "cddefines.h"
#include "taulines.h"
#include "atomfeii.h"
#include "dense.h"
#include "conv.h"
#include "atoms.h"
#include "rfield.h"
#include "wind.h"
#include "iso.h"
#include "h2.h"
#include "opacity.h"
#include "trace.h"
#include "lines_service.h"
#include "atmdat.h"
#include "hydrogenic.h"
#include "rt.h"
#include "cosmology.h"
#include "physconst.h"
#include "doppvel.h"
#include "mole.h"

/*RT_line_all do escape and destruction probabilities for all lines in code. */
void RT_line_all( void )
{
	long int i,
		ion,
		ipISO,
		nelem;
	long ipHi , ipLo;

	/* this flag says whether to update the line escape probabilities */
	bool lgPescUpdate = conv.lgFirstSweepThisZone || conv.lgIonStageTrimed;

	DEBUG_ENTRY( "RT_line_all()" );

	if( trace.lgTrace )
		fprintf( ioQQQ, "     RT_line_all called\n" );

	/* skip line transfer if requested, and not first call in this zone */
	if( !rfield.lgDoLineTrans && conv.nPres2Ioniz )
	{
		return;
	}

	/* this array is huge and takes significant time to zero out or update, 
	 * only do so when needed, */
	if( conv.lgLastSweepThisZone )
	{
		/* zero out fine opacity array */
		memset(rfield.fine_opac_zone , 0 , (unsigned long)rfield.nfine*sizeof(realnum) );

		if( 0 && cosmology.lgDo )
		{
			realnum dVel = (realnum)SPEEDLIGHT * ( 1.f - 1.f/POW2(1.f+cosmology.redshift_start-cosmology.redshift_current) );
			rfield.ipFineConVelShift = -(long int)( dVel/ rfield.fine_opac_velocity_width + 0.5 );
		}
		else
		{
			/* this is offset within fine opacity array caused by Doppler shift
			 * between first zone, the standard of rest, and current position 
			 * in case of acceleration away from star in outflowing wind 
			 * dWind is positive, this means that the rest frame original
			 * velocity is red shifted to lower energy */
			realnum dWind = wind.windv - wind.windv0;

			/* will add ipVelShift to index of original energy so redshift has 
			 * to be negative */
			rfield.ipFineConVelShift = -(long int)(dWind / rfield.fine_opac_velocity_width+0.5);
		}
	}

#	if 0
	/* this code is a copy of the code within line_one that does cloaking
	 * within this zone.  it can be used to see how a particular line
	 * is being treated. */
	{
#include "doppvel.h"
		double dTau , aa;
		ipISO = 0; nelem = 0;ipLo = 0;
		ipHi = 23;
		dTau =  Transitions[ipISO][nelem][ipHi][ipLo].Emis->PopOpc * 
			Transitions[ipISO][nelem][ipHi][ipLo].Emis->opacity / 
			GetDopplerWidth(dense.AtomicWeight[nelem]) + opac.opacity_abs[Transitions[ipISO][nelem][ipHi][ipLo].ipCont-1];
		dTau *= radius.drad_x_fillfac_mean;
		aa = log(1. + dTau ) / SDIV(dTau);
		fprintf(ioQQQ,"DEBUG dTau\t%.2f\t%.5e\t%.5e\t%.5e\n",fnzone,dTau,
			radius.drad_x_fillfac_mean,
			 aa);
	}
#	endif

	/* find Stark escape probabilities for hydrogen itself */
	if( lgPescUpdate )
		RT_stark();

	/*if( lgUpdateFineOpac )
		fprintf(ioQQQ,"debuggg\tlgUpdateFineOpac in rt_line_all\n");*/
	for( ipISO=ipH_LIKE; ipISO < NISO; ++ipISO )
	{
		/* loop over all iso-electronic sequences */
		for( nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			/* parent ion stage, for H is 1, for He is 1 for He-like and 
			 * 2 for H-like */
			ion = nelem+1-ipISO;

			/* element turned off */
			if( !dense.lgElmtOn[nelem] )
				continue;
			/* need we consider this ion? */
			if( ion <= dense.IonHigh[nelem] )
			{
				/* loop over all lines */
				for( ipHi=1; ipHi < iso.numLevels_local[ipISO][nelem]; ++ipHi )
				{
					for( ipLo=0; ipLo < ipHi; ++ipLo )
					{
						/* negative ipCont means this is not a real line, so do not
						 * transfer it */
						if( Transitions[ipISO][nelem][ipHi][ipLo].ipCont < 0 ) 
							continue;

						/* generate escape prob, pumping rate, destruction prob, 
						 * inward outward fracs */
						RT_line_one( &Transitions[ipISO][nelem][ipHi][ipLo],
							     true,(realnum)iso.pestrk[ipISO][nelem][ipLo][ipHi],
								 GetDopplerWidth(dense.AtomicWeight[nelem]));

						/* set true to print pump rates*/
						enum {DEBUG_LOC=false};
						if( DEBUG_LOC && nelem==1&& ipLo==0  /*&& iteration==2*/ )
						{
							fprintf(ioQQQ,"DEBUG pdest\t%3li\t%.2f\t%.3e\n",
								ipHi ,
								fnzone,
								Transitions[ipISO][nelem][ipHi][ipLo].Emis->Pdest);
						}
					}
				}

				/* now update two photon induced rates */
				atmdat_2phot_rate(ipISO , nelem);

				/* this is option to not do destruction probabilities
				 * for case b no pdest option */
				if( opac.lgCaseB_no_pdest )
				{
					ipLo = 0;
					/* okay to let this go to numLevels_max. */
					for( ipHi=ipLo+1; ipHi < iso.numLevels_max[ipISO][nelem]; ++ipHi )
					{
						if( Transitions[ipISO][nelem][ipHi][ipLo].ipCont <= 0 )
							continue;

						/* hose the previously computed destruction probability */
						Transitions[ipISO][nelem][ipHi][ipLo].Emis->Pdest = SMALLFLOAT; 
					}
				}

				ipLo = 0;
				/* these are the extra Lyman lines for the iso sequences */
				/*for( ipHi=2; ipHi < iso.nLyman[ipISO]; ipHi++ )*/
				/* only update if significant abundance and need to update fine opac */
				if( dense.xIonDense[nelem][ion] > 1e-30 && conv.lgLastSweepThisZone )
				{
					for( ipHi=StatesElemNEW[nelem][nelem-ipISO][iso.numLevels_max[ipISO][nelem]-1].n; ipHi < iso.nLyman[ipISO]; ipHi++ )
					{
						/* we just want the population of the ground state */
						ExtraLymanLines[ipISO][nelem][ipHi].Emis->PopOpc = StatesElemNEW[nelem][nelem-ipISO][0].Pop;
						ExtraLymanLines[ipISO][nelem][ipHi].Lo->Pop =
							StatesElemNEW[nelem][nelem-ipISO][ipLo].Pop;

						/* actually do the work */
						RT_line_one( &ExtraLymanLines[ipISO][nelem][ipHi], true, 0.f,
							GetDopplerWidth(dense.AtomicWeight[nelem]));
					}
				}

			}/* if nelem if ion <=dense.IonHigh */
		}/* loop over nelem */
	}/* loop over ipISO */

	/* note that pesc and dest are updated no matter what escprob logic we
	 * specify, and that these are not updated if we have overrun the
	 * optical depth scale.  only update here in that case */
	if( lgTauGood( &Transitions[ipH_LIKE][ipHYDROGEN][iso.nLyaLevel[ipH_LIKE]][0] ) )
	{
		/* add on destruction of hydrogen Lya by FeII
		 * now add in FeII deexcitation via overlap,
		 * but only as loss, not photoionization, source
		 * dstfe2lya is Ly-alpha deexcit by FeII overlap - conversion into FeII em */
		/* find FeII overlap destruction rate, 
		 * this does NOTHING when large FeII atom is turned on, 
		 * in this case evaluation is done in call to FeIILevelPops */
		t_fe2ovr_la::Inst().atoms_fe2ovr();
		/*fprintf(ioQQQ,"DEBUG fe2 %.2e %.2e\n", hydro.dstfe2lya ,
			hydro.HLineWidth);*/

		/* >>chng 00 jan 06, let dest be large than 1 to desaturate the atom */
		/* >>chng 01 apr 01, add test for tout >= 0., 
		 * must not add to Pdest when it has not been refreshed here */
		Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->Pdest += hydro.dstfe2lya;
	}

	/* is continuum pumping of H Lyman lines included?  yes, but turned off
	 * with atom h-like Lyman pumping off command */
	if( !hydro.lgLymanPumping )
	{
		ipISO = ipH_LIKE;
		nelem = ipHYDROGEN;
		ipLo = 0;
		for( ipHi=ipLo+1; ipHi < iso.numLevels_max[ipISO][nelem]; ++ipHi )
		{
			Transitions[ipISO][nelem][ipHi][ipLo].Emis->pump = 0.;
		}
	}

	{
		/* following should be set true to print ots contributors for he-like lines*/
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC && nzone>433 /*&& iteration==2*/ )
		{
			/* option to dump a line  */
			DumpLine(&Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][0] );
		}
	}

	/* level 1 lines */
	for( i=1; i <= nLevel1; i++ )
	{
		RT_line_one( &TauLines[i], true,0.f, GetDopplerWidth(dense.AtomicWeight[TauLines[i].Hi->nelem-1]) );
	}

	{
		// 09 mar 28, this transition, between the highest two levels in the
		// very simple version of the model Fe II atom, is troublesome due 
		// to Lya pumping, which can invert the highest two level populations.  
		// The electron scattering escape probability shows noise as a result.  
		// The sim nova_photos will throw an assert if this is removed
		// This model is unphysical since highest two levels are close together
		// and Lya is populating the highest one.  The problem does not occur
		// when the large FeII atom is turned on - that has realistic level
		// energies.  
		// evaluate electron scattering escape probability one time per zone
		static realnum P_elec_esc_ipTFe56;
		static long nSave;
		if( nzone <= 1 || nzone!=nSave )
		{
			P_elec_esc_ipTFe56 = TauLines[ipTFe56].Emis->Pelec_esc;
			nSave = nzone;
		}
		else
			TauLines[ipTFe56].Emis->Pelec_esc = P_elec_esc_ipTFe56;
	}

	/* external database lines */
	for( long ipSpecies=0; ipSpecies<nSpecies; ipSpecies++ )
	{
		if( Species[ipSpecies].lgActive )
		{
			for( long ipHi=1; ipHi<Species[ipSpecies].numLevels_local; ipHi++ )
			{
				transition *t = &dBaseTrans[ipSpecies][ipHi][0];
				for(long ipLo = 0; ipLo < ipHi; ipLo++ )
				{
					if( t->ipCont > 0 )
					{
						RT_line_one( t, true, 0.f, GetDopplerWidth( t->Hi->sp->fmolweight ) );
					}
					++t;
				}
			}
		}
	}

	/* this is a major time sink for this routine - only evaluate on last
	 * sweep when fine opacities are updated since only effect of UTAs is
	 * to pump inner shell lines and add to total opacity */
	if( conv.lgFirstSweepThisZone || conv.lgLastSweepThisZone )
	{
		for( i=0; i < nUTA; i++ )
		{
			/* these are not defined in cooling routines so must be set here */
			UTALines[i].Emis->PopOpc = dense.xIonDense[UTALines[i].Hi->nelem-1][UTALines[i].Hi->IonStg-1];
			UTALines[i].Lo->Pop = dense.xIonDense[UTALines[i].Hi->nelem-1][UTALines[i].Hi->IonStg-1];
			UTALines[i].Hi->Pop = 0.;
			RT_line_one( &UTALines[i], true,0.f, 
				GetDopplerWidth(dense.AtomicWeight[UTALines[i].Hi->nelem-1]) );
		}
	}

	/* do satellite lines for iso sequences gt H-like
	 * H-like ions only have one electron, no satellite lines. */
	for( ipISO=ipHE_LIKE; ipISO<NISO; ++ipISO )
	{
		/* loop over all iso-electronic sequences */
		for( nelem=ipISO; nelem<LIMELM; ++nelem )
		{
			if( dense.lgElmtOn[nelem] && iso.lgDielRecom[ipISO] )
			{
				for( ipLo=0; ipLo < iso.numLevels_local[ipISO][nelem]; ++ipLo )
				{
					RT_line_one( &SatelliteLines[ipISO][nelem][ipLo], true,0.f, 
						GetDopplerWidth(dense.AtomicWeight[nelem]) );
				}
			}
		}
	}

	/* the large H2 molecule */
	H2_RTMake( );
	return;
}
