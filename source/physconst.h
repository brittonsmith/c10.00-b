/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef PHYSCONST_H_
#define PHYSCONST_H_

/**\file physconst.h \verbatim
 * physical constants used by Cloudy, mostly taken from
 * >>refer	phys	const	Mohr P.J., Taylor B.N., & Newell D.B., Codata 2006, http://www.physics.nist.gov/constants
 * <BR><BR>
 * NB - these are all printed with the "print constants" command, 
 * which is in parse_print.cpp, so any new constants
 * added here must also be added to the prt_constants routine 
 * this is in the func_test case in the auto test suite \endverbatim
 */
/*#include "physconst.h"*/

/*********************************************************************
 * first come math constants                                         *
 *********************************************************************/

/** the number e */
const double EE = 2.718281828459045235360287;

/** the Euler constant (aka Euler-Mascheroni constant or gamma) */
const double EULER = 0.577215664901532860606512090082;

/** pi */
const double PI = 3.141592653589793238462643;

/** 2*pi */
const double PI2 = 6.283185307179586476925287;

/** 4*pi */
const double PI4 = 12.56637061435917295385057;

/** 8*pi */
const double PI8 = 25.13274122871834590770115;

/** sqrt(2) */
const double SQRT2 = 1.414213562373095048801689;

/** sqrt(pi) */
const double SQRTPI = 1.772453850905516027298167;

/** sqrt(pi/2) */
const double SQRTPIBY2 = 1.253314137315500251207883;

/** ln(2) */
const double LN_TWO = 0.6931471805599453094172321;

/** ln(10) */
const double LN_TEN = 2.302585092994045684017991;

/** log(e) */
const double LOG10_E = 0.4342944819032518276511289;

/** factor that converts optical depth into extinction in mags,
 * 2.5 log e */
const double OPTDEP2EXTIN = 1.085736204758129569127822;

/** 180/pi */
const double RADIAN = 57.29577951308232087679815;

/*********************************************************************
 * astronomical constants go here                                    *
 *********************************************************************/

/** solar mass in gram
 * >>refer	phys	const	http://pdg.lbl.gov/2010/reviews/rpp2010-rev-astrophysical-constants.pdf */
const double SOLAR_MASS = 1.9884e33;

/** solar luminosity erg s-1
 * >>refer	phys	const	http://pdg.lbl.gov/2010/reviews/rpp2010-rev-astrophysical-constants.pdf */
const double SOLAR_LUMINOSITY = 3.8427e33;

/** astronomical unit, cm, nearly the length of the semimajor
 * axis of the Earth's elliptical orbit around the sun */
/* >>refer	phys	const	http://pdg.lbl.gov/2010/reviews/rpp2010-rev-astrophysical-constants.pdf */
const double AU = 1.49597870700e13;

/*********************************************************************
 * fundamental constants go next, eventually rest should be defined  *
 * in terms of these, these are Codata 2006 values.                  *
 *********************************************************************/

/** atomic mass unit, gram */
const double ATOMIC_MASS_UNIT = 1.660538782e-24;

/** electron mass, gram */
const double ELECTRON_MASS = 9.10938215e-28;

/** proton mass, gram */
const double PROTON_MASS = 1.672621637e-24;

/** this is the Boltzmann factor, erg/K */
const double BOLTZMANN = 1.3806504e-16;

/** speed of light, cm/s */
const double SPEEDLIGHT = 2.99792458e10;

/** Planck's constant */
const double HPLANCK = 6.62606896e-27;

/** Gravitational constant, cm^3/g/s^2 */
const double GRAV_CONST = 6.67428e-8;

/** elementary charge, in C in SI units, to use this must convert to cgs */
const double ELEM_CHARGE = 1.602176487e-19;

/** infinite mass rydberg constant, in cm^-1 */
const double RYD_INF = 1.0973731568527e5;

/** ionization potential of real hydrogen atom, in inf mass ryd, based on Codata 98,
 * uncertainty 24e-12, calculated by Peter van Hoof */
const double HIONPOT = 0.999466508219;

/*********************************************************************
 * below here should be derived constants                            *
 *                                                                   *
 * NB - explicit values in comments are approximate                  *
 *      and are not maintained !                                     *
 *********************************************************************/

/** number of arcsec in 1 radian, 206264.806 */
const double AS1RAD = RADIAN*3600.;

/** number of square arcsec in 1 steradian, 4.254517e10 */
const double SQAS1SR = pow2(AS1RAD);

/** number of square arcsec in the whole sky, 5.3463838e11 */
const double SQAS_SKY = PI4*SQAS1SR;

/** parsec in cm, 3.085678e18 */
const double PARSEC = AU*AS1RAD;

/** megaparsec in cm, 3.085678e24 */
const double MEGAPARSEC = 1.e6*PARSEC;

/** h/2pi = 1.05457e-27 */
const double H_BAR = HPLANCK/(2.*PI);

/** elementary charge, in ESU, 4.8032e-10 */
const double ELEM_CHARGE_ESU = ELEM_CHARGE*SPEEDLIGHT/10.;

/** electric constant, in F/m, 8.854e-12 */
const double ELECTRIC_CONST = 1.e11/(PI4*pow2(SPEEDLIGHT));

/** this is the factor that appears in front of Boltzmann factor to get
 * LTE level populations for hydrogenic ions. It is given in the
 * first parts of section 5 of part 2 of hazy, and is actually
 * ( planck^2 / (2 pi m_e k ) )^3/2, but cannot evaluate powers here,
 * must raise this to 3/2 when used, HION_LTE_POP = 5.556e-11 cm^2 K */
const double HION_LTE_POP = pow2(HPLANCK)/(PI2*BOLTZMANN*ELECTRON_MASS);

/** SAHA is ( h^2/2/pi/m/k )^3/2, is correct constant for free electron
 * SAHA = 4.14132e-16 cm^3 K^(3/2) */
const double SAHA = sqrt(pow3(HION_LTE_POP));

/** number of ergs per wavenumber, 1.9864e-16 */
const double ERG1CM = HPLANCK*SPEEDLIGHT;

/** degrees kelvin per unit wavenumber, 1.4388 */
const double T1CM = HPLANCK*SPEEDLIGHT/BOLTZMANN;

/** number of Ryd per wavenumber, 9.11267e-6 */
const double WAVNRYD = 1./RYD_INF;

/** Angstrom per infinite mass Ryd, 911.2671 */
const double RYDLAM = 1.e8/RYD_INF;

/** ergs per inf mass Ryd, 2.180e-11 */
const double EN1RYD = HPLANCK*SPEEDLIGHT*RYD_INF;

/** the temperature of 1 Rydberg
 te1ryd is h/k is temp of 1 Rydberg, 1.579e5 */
const double TE1RYD = HPLANCK*SPEEDLIGHT*RYD_INF/BOLTZMANN;

/** Kelvins per eV, 1.1604e4 */
const double EVDEGK = ELEM_CHARGE*1.e7/BOLTZMANN;

/** eV per inf mass Ryd, 13.606 */
const double EVRYD = HPLANCK*SPEEDLIGHT*RYD_INF/ELEM_CHARGE*1.e-7;

/** ergs per eV, 1.602176e-012 */
const double EN1EV = EN1RYD/EVRYD;

/** frequency of one Ryd for infinite mass nuclei, 3.289842e15 */
const double FR1RYD = SPEEDLIGHT*RYD_INF;

/**2 h FR1RYD^3 / c^2 for infinite mass nucleus, 0.5250 */
const double HNU3C2 = 2.*HPLANCK*SPEEDLIGHT*pow3(RYD_INF);

/** frequency of ionization potential of H (not inf mass), 3.288087e15 - never used */
const double FR1RYDHYD = SPEEDLIGHT*RYD_INF*HIONPOT;

/** H_BAR in eV sec, 6.582e-16 */
const double HBAReV = H_BAR/EN1EV;  

/** wavelength (A) of ionization potential of Hydrogen, 911.7535 - never used */
const double RYDLAMHYD = RYDLAM/HIONPOT;

/** Stefan-Boltzmann constant, 5.6704e-5 */
const double STEFAN_BOLTZ = pow2(PI*pow2(BOLTZMANN))/(60.*pow3(H_BAR)*pow2(SPEEDLIGHT));

/** the frequency of one eV, 2.418e14 */
const double FREQ_1EV = SPEEDLIGHT*RYD_INF/EVRYD;

/** the fine-structure constant a= 2pi e^2/hc 7.297 352 533 x 10-3 */
const double FINE_STRUCTURE = pow2(ELEM_CHARGE_ESU)/SPEEDLIGHT/H_BAR;

/** the square of the fine-structure constant */
const double FINE_STRUCTURE2 = pow2(FINE_STRUCTURE);

/** Bohr radius in cm, 5.29177249e-9 */
const double BOHR_RADIUS_CM = FINE_STRUCTURE/(PI4*RYD_INF);

/** the two photon constant as defined by Breit & Teller, as in equation 4 of Spitzer & Greenstein 51, 2.18313 */
const double TWO_PHOT_CONST = 9.*pow3(FINE_STRUCTURE2)*FR1RYD/2048.;

/** this is the square of the value roughly equal to 8.629e-6 that appears in converting 
 * collision strengths to rates. The constant is h^2/((2PI*me)^3/2 * k^1/2). */
const double COLL_CONST = SAHA*BOLTZMANN/HPLANCK;

/** this is the square of the value roughly equal to 4.123e11 that appears in the integration
 * of photoionization cross-sections to obtain recombination coefficients. */
const double MILNE_CONST = SPEEDLIGHT*sqrt(pow3(FINE_STRUCTURE2)*pow3(TE1RYD)/PI);

/** This is the constant used in converting oscillator strengths to As. The formula is
 * Aul = TRANS_PROB_CONST * f(u,l) * wavenumber^2. TRANS_PROB_CONST is 0.667025 */ 
const double TRANS_PROB_CONST = PI4*HPLANCK*FINE_STRUCTURE/ELECTRON_MASS;

#endif /* PHYSCONST_H_ */
