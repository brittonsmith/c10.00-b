/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef CA_H_
#define CA_H_

/** ca2mly.h */
EXTERN struct t_ca {

	/** amount of Lya removed by photoionization of excited CaII levels */
	realnum Ca2RmLya;

	/** cooling due to CaII */
	realnum Cak, 
	  Cah, 
	  Cax, 
	  Cay, 
	  Caz, 
	  Caf1, 
	  Caf2;
	double
	  Cair, 
	  c7306, 
	  Cakh, 
	  Ca3688, 
	  Ca2112, 
	  Ca5620, 
	  Ca4941, 
	  c3997, 
	  c2414, 
	  Ca6087,
	  c5311;

	/** parameters dealing with dest of CaII by Lya */
	realnum dstCala, 
	  dCakh, 
	  dCaf12, 
	  Ca3d, 
	  Ca4p;

	/** summed pop of CaII excited states */
	realnum popca2ex;

	}	ca;


#endif /* CA_H_ */
