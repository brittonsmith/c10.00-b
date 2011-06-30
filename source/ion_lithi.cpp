/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*IonLithi compute ionization balance for lithium */
#include "cddefines.h"
#include "dense.h"
#include "ionbal.h"

void IonLithi(void)
{
	const int NDIM = ipLITHIUM+1;

	static const double dicoef[2][NDIM] = { {2.54e-3,6.15e-3,0.}, {4.42e-2,5.88e-2,0.} };
	static const double dite[2][NDIM] = { {1.57e5,1.41e5,0.}, {3.74e5,1.41e5,0.} };
	static const double ditcrt[NDIM] = {1.2e4,1.2e4,1e20};
	static const double aa[NDIM] = {0.,0.,0.};
	static const double bb[NDIM] = {0.,0.,0.};
	static const double cc[NDIM] = {0.,0.,0.};
	static const double dd[NDIM] = {0.,0.,0.};
	static const double ff[NDIM] = {0.1,0.1,0.};

	DEBUG_ENTRY( "IonLithi()" );

	/* lithium nelem=3
	 * data are for carbon
	 *
	 * fluorescent fields, first dim is stage of ionization, sec is shell
	 *
	 * rates from Shull and van Steenberg, Ap.J. Sup 48, 95. */
	/* DATA GRDEFF/0.10,0.10,0.10/
	 * GRDEFF is fraction of recombinations to ground state, used for
	 * outward diffuse fields
	 *
	 * rec from +3, +4 from Arnaud et al Ast Ap Sup 60 425. (1985)
	 * rec from fully ionized uses Seaton '79 in ionrat */
	/* Pequignot and Aldrovandi Ast Ap 161, 169. */

	if( !dense.lgElmtOn[ipLITHIUM] )
	{
		return;
	}

	/* zero out ionization balance arrays */
	ion_zero(ipLITHIUM);

	ion_photo(ipLITHIUM,false);

	/* find collisional ionization rates */
	ion_collis(ipLITHIUM);

	/* get recombination coefficients */
	ion_recomb(false,(const double*)dicoef,(const double*)dite,ditcrt,aa,bb,cc,dd,ff,ipLITHIUM);

	/* solve for ionization balance */
	ion_solver(ipLITHIUM,false);
	return;
}
