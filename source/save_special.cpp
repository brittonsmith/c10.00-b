/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*SaveSpecial generate output for the save special command */
#include "cddefines.h"
#include "wind.h"
#include "opacity.h"
#include "dense.h"
#include "taulines.h"
#include "radius.h"
#include "phycon.h"
#include "save.h"

/*SaveSpecial generate output for the save special command */
void SaveSpecial(FILE* ioPUN , 
  const char *chTime)
{
	/*long int i;*/

	DEBUG_ENTRY( "SaveSpecial()" );

	if( strncmp(chTime,"LAST",4) == 0 )
	{
		/* code to execute only after last zone */
#		if 0
		long ipISO , nelem , limit , i;
		double EdenAbund , fach;
#		include "physconst.h"
#		include "hydrogenic.h"
		PunFeII( ioPUN );*/
		ipISO = ipHYDROGEN;
		nelem = ipHYDROGEN;

		/* in all following the factor of two is because a single
		 * decay produces two photons */
		EdenAbund = StatesElemNEW[nelem][nelem-ipH_LIKE][ipH2s].Pop*8.226*pow(1.+nelem,6);
		fprintf(ioPUN," 2s = %.3e\n", EdenAbund);

		/* upper limit to H-like 2-phot is energy of La, which is in ipCont-1 cell */
		limit = Transitions[ipH_LIKE][nelem][ipH2p][ipH1s].ipCont-1;
		/* remember sum of rates, this will add up to twice the real rate since
		 * each transition makes two photons */
		for( i=0; i < limit; i++ )
		{
			/*>>chng 01 jan 23, previous change had doubled cross section for H two-photon,
			 * so here we divide by 2 to get old answer */
			/** \todo	2	this most likely needs to be changed in light of new 2nu treatment	*/
			fach = iso.As2nu[ipISO][nelem][i]/2.f;
			fach *= rfield.anu2[i]/rfield.widflx[i]*EN1RYD;
			fprintf(ioPUN,"%.3e\t%.3e\t%.3e\n", 
				RYDLAM/1e4/rfield.anu[i] , fach , fach*(realnum)EdenAbund );
		}
#		endif

	}
	else
	{
		/* code to do for every zone */
		fprintf(ioPUN,"%.5e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n",
			radius.Radius ,
			wind.AccelCont ,
			wind.fmul ,
			opac.opacity_sct[1000],
			dense.eden , 
			dense.xMassDensity,
			dense.gas_phase[ipHYDROGEN] );

#		if 0
		long int iElecHi=1 , iVibHi=0, iRotHi=0;
		/* code to execute after every zone */
		if( h2.lgH2ON )
		{
			ASSERT( H2Lines[iElecHi][iVibHi][iRotHi][0][0][1].ipCont > 0 );

			fprintf(ioPUN,"DEBUG Oion\t%li\t%.2f\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
				iteration,fnzone,
				radius.depth,
				H2Lines[iElecHi][iVibHi][iRotHi][0][0][1].Lo->Pop,
				H2Lines[iElecHi][iVibHi][iRotHi][0][0][1].Emis->PopOpc,
				H2Lines[iElecHi][iVibHi][iRotHi][0][0][1].Emis->opacity,
				H2Lines[iElecHi][iVibHi][iRotHi][0][0][1].Emis->TauCon,
				H2Lines[iElecHi][iVibHi][iRotHi][0][0][1].Emis->TauIn,
				H2Lines[iElecHi][iVibHi][iRotHi][0][0][1].Emis->pump,
				radius.drad
				);
		}
#		endif
		/*DumpLine(&Transitions[ipHE_LIKE][ipHELIUM][ipHe2p1P][0] );*/

	}
	return;
}
