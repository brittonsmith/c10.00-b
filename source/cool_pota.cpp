/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolPota compute potassium cooling */
#include "cddefines.h"
#include "taulines.h"
#include "phycon.h"
#include "lines_service.h"
#include "atoms.h"
#include "cooling.h"

void CoolPota(void)
{
	double cs;

	DEBUG_ENTRY( "CoolPota()" );

	/* potasium lines
	 * KI 7745 */
	cs = 7.231e-4*phycon.te*phycon.te03*phycon.te02;
	PutCS(cs,&TauLines[ipKI7745]);
	atom_level2(&TauLines[ipKI7745]);

	/* [K III] 4.62 microns
	 * Y(ik) from 
	 * >>refer	k3	cs	Pelan, J., & Berrington, K.A. 1995, A&A Suppl, 110, 209 */
	PutCS(2.2,&TauLines[ipxK03462]);
	atom_level2(&TauLines[ipxK03462]);

	/* [KIV] 5.983, 15.39 mic, cs from 
	 * >>refer	k4	cs	Galavis, M.E., Mendoza, C., & Zeippen, C.J. 1995, A&AS, 111, 347 */
	PutCS(4.3,&TauLines[ipxK04598]);
	PutCS(1.13,&TauLines[ipxK04154]);
	PutCS(1.3,&TauDummy);
	/* atom_level3(  t10,t21,t20) */
	atom_level3(&TauLines[ipxK04598],&TauLines[ipxK04154],&TauDummy);

	/* [KVI] 8.823, 5.575 mic, cs from 
	 * >>refer	k6	cs	Galavis, M.E., Mendoza, C., & Zeippen, C.J. 1995, A&AS, 111, 347 */
	cs = MIN2(1.505,0.274*phycon.te10*phycon.te05/phycon.te001/
	  phycon.te001);
	PutCS(cs,&TauLines[ipxK06882]);

	cs = MIN2(4.632,1.909*phycon.te10/phycon.te003);
	cs = MAX2(4.0,cs);
	PutCS(cs,&TauLines[ipxK06557]);
	PutCS(1.2,&TauDummy);

	atom_level3(&TauLines[ipxK06882],&TauLines[ipxK06557],&TauDummy);

	/* [K VII] 3.189 microns cs from 
	 * >>refer	k7	cs	Saraph, H.E., & Storey, P.J. A&AS, 115, 151 */
	PutCS(4.5,&TauLines[ipxK07319]);
	atom_level2(&TauLines[ipxK07319]);

	/* K 11 4249.99A, cs from 
	 * >>refer	k11	cs	Saraph, H.E. & Tully, J.A. 1994, A&AS, 107, 29 */
	cs = MIN2(0.172,0.0109*phycon.te20*phycon.te02/
	  phycon.te001/phycon.te001);
	cs = MAX2(0.111,cs);
	PutCS(0.115,&TauLines[ipxK11425]);

	atom_level2(&TauLines[ipxK11425]);

	return;
}
