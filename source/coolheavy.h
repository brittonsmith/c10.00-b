/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef COOLHEAVY_H_
#define COOLHEAVY_H_

/** coolheavy.h */
EXTERN struct t_CoolHeavy {
	double c1170, 
	  c2428, 
	  c2125, 
	  c7136, 
	  c7007, 
	  c4626,
	  c2691,
	  fs4413, 
	  c409, 
	  c5192, 
	  c3109, 
	  Ar4740, 
	  Ar4711, 
	  Ar2868, 
	  Ar2854, 
	  Ar7263, 
	  Ar7171, 
	  Ar7331, 
	  Ar7237, 
	  c8579, 
	  c6164, 
	  c3679, 
	  c5525, 
	  c3350, 
	  c8494, 
	  Cl5539, 
	  Cl5519, 
	  Cl3354, 
	  Cl3344, 
	  Cl8504, 
	  Cl8436, 
	  Cl8552, 
	  Cl8483, 
	  c8047, 
	  c3119, 
	  c5324, 
	  Cr4l31, 
	  Cr4l21, 
	  Cr5l31, 
	  Cr5l21, 
	  Cr5l32, 
	  Cr3l21, 
	  c3892, 
	  c5177, 
	  c2568, 
	  c1242, 
	  c2170, 
	  c5303, 
	  c7082, 
	  c1118, 
	  c1328, 
	  c592, 
	  c263, 
	  c242, 
	  Fe231, 
	  Fe232, 
	  Fe221, 
	  c1806, 
	  c2569, 
	  c1357, 
	  c2972, 
	  c1365, 
	  c2067, 
	  c4017, 
	  c3343, 
	  c3869, 
	  c4720, 
	  c2424, 
	  c2975, 
	  c1565, 
	  c3426, 
	  c1134, 
	  eebrm, 
	  h2line, 
	  HD,
	  colmet, 
	  tccool, 
	  expans, 
	  cextxx, 
	  c5577, 
	  c6300, 
	  c6300_frac_emit,
	  c5577_frac_emit,
	  c3727, 
	  c7325, 
	  c4363, 
	  c5007, 
	  coolOi, 
	  c6363, 
	  c6731, 
	  c10330, 
	  c6312, 
	  c9532, 
	  c1715, 
	  c1213, 
	  c1826, 
	  c488, 
	  Sc22p08m, 
	  Sc24p1m, 
	  Sc24p2m, 
	  Sc33936, 
	  Sc42100, 
	  Sc45058, 
	  Sc43595, 
	  Ti3l31, 
	  Ti3l32, 
	  Ti3l21, 
	  V38830, 
	  V47741, 
	  V49496, 
	  V44p2m, 
	  V38507, 
	  cyntrn, 
	  heavfb,
	  p2_21,
	  p2_32,
	  p2_31;
	realnum S6733, 
	  S6718, 
	  S4070, 
	  S4078, 
	  S10323, 
	  S10289, 
	  S10373, 
	  S10339;

	  /* Fe12 lines added 01 Aug
	  double fe13_1216 , fe13_3000 , fe13_2302; */

	  /** flag set false if no free-free heating command entered */
	  bool lgFreeOn;
	  double brems_cool_h, 
		 brems_cool_hminus,
		 brems_cool_he, 
		 brems_cool_metals, 
		 /** total brems heating, all opacity sources */
		 brems_heat_total,
		 /** net brems cooling, sum of cool minus heat */
		 brems_cool_net;

	/** the [OII] lines */
	realnum O3730, 
	  O3726, 
	  O2471, 
	  O7323, 
	  O7332;

	/** these are ratios of radiative to total decays from n=2 and 3
	 * of excited OII, used for recombination contribution */
	double O2_A3_tot ,
		O2_A2_tot;

	/** cooling due to H + H+ => H2+ */
	realnum H2PlsCool;

	}	CoolHeavy;

#endif /* COOLHEAVY_H_ */
