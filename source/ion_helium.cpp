/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*IonHelium solve ionization balance for helium */
#include "cddefines.h"
#include "dense.h"
#include "trace.h"
#include "conv.h"
#include "iso.h"
#include "ionbal.h"
#include "taulines.h"

void IonHelium( void )
{
	bool lgDebugPrint=false;

	DEBUG_ENTRY( "IonHelium()" );

	/* option to "turn off" helium */
	if( !dense.lgElmtOn[ipHELIUM] )
	{
		dense.xIonDense[ipHELIUM][0] = 0.;
		dense.xIonDense[ipHELIUM][1] = 0.;
		dense.xIonDense[ipHELIUM][2] = 0.;
		return;
	}

	/* populations */
	/* >>chng 01 may 09, add option to set ionization with element name ioniz cmnd */
	if( dense.lgSetIoniz[ipHELIUM] )
	{
		dense.xIonDense[ipHELIUM][2] = dense.SetIoniz[ipHELIUM][2]*dense.gas_phase[ipHELIUM];
		dense.xIonDense[ipHELIUM][1] = dense.SetIoniz[ipHELIUM][1]*dense.gas_phase[ipHELIUM];
		dense.xIonDense[ipHELIUM][0] = dense.SetIoniz[ipHELIUM][0]*dense.gas_phase[ipHELIUM];
		return;
	}

	lgDebugPrint = false;
#	if 0
	if( nzone > 197 )
		lgDebugPrint = true;
#	endif

	/* find ionization balance */
	ion_solver( ipHELIUM , lgDebugPrint );

	if( trace.lgHeBug )
	{
		fprintf( ioQQQ, "  %li IonHelium returns;  nzone %ld He0:%.4e He+:%.4e He+2:%.4e\n",
			conv.nTotalIoniz ,
			nzone,
			dense.xIonDense[ipHELIUM][0]/dense.gas_phase[ipHELIUM] ,
			dense.xIonDense[ipHELIUM][1]/dense.gas_phase[ipHELIUM] ,
			dense.xIonDense[ipHELIUM][2]/dense.gas_phase[ipHELIUM] );

		fprintf( ioQQQ, "     He+ /He0:%s smp:%.4e rec:%.4e ion:%.4e rad rec:%.4e 1s Pop:%.4e\n",
			iso.chTypeAtomUsed[ipHE_LIKE][ipHELIUM],
			iso.xIonSimple[ipHE_LIKE][ipHELIUM],
			ionbal.RateRecomTot[ipHELIUM][0],
			ionbal.RateIonizTot(ipHELIUM,0),
			ionbal.RR_rate_coef_used[ipHELIUM][0],
			StatesElemNEW[ipHELIUM][ipHELIUM-ipHE_LIKE][0].Pop );

		fprintf( ioQQQ, "     He+2/He+:%s smp:%.4e rec:%.4e ion:%.4e rad rec:%.4e\n",
			iso.chTypeAtomUsed[ipH_LIKE][ipHELIUM],
			iso.xIonSimple[ipH_LIKE][ipHELIUM] ,
			ionbal.RateRecomTot[ipHELIUM][1],
			ionbal.RateIonizTot(ipHELIUM,1),
			ionbal.RR_rate_coef_used[ipHELIUM][1] );

		fprintf( ioQQQ, "\n" );
	}
	return;
}
