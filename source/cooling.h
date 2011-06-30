/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef COOLING_H_
#define COOLING_H_


/**CoolZero set cooling and heating stack to zero */
void CoolZero(void);

/**CoolAdd add coolants to the cooling stack, called in evaluation of cooling function 
\param *chLabel
\param xlambda
\param cool
*/
void CoolAdd(
  const char *chLabel, 
  realnum xlambda, 
  double cool);

/**CoolSum  total cooling from all entries into cooling stack */
void CoolSum(double *total);

/**CoolEvaluate main routine to call others, to evaluate total cooling 
\param tot total cooling */
void CoolEvaluate(double *tot);

/**coolpr stores coolants before block printed, when printing cooling agents 
\param *io the label for the coolant
\param *chLabel
\param lambda the wavelength
\param ratio the ratio of this coolant, to total cooling, may be negative
\param *chJOB which job, either ZERO, DOIT, or DONE 
*/
void coolpr(
	FILE *io,
	const char *chLabel ,
	realnum lambda,
	double ratio,
	const char *chJOB );

/**HeatSum evaluate all heating agents to determine total heating for this zone,
 * called at end of ionize */
void HeatSum(void);

/**HeatZero zeroes out the heating array, called at start of ionize*/
void HeatZero(void);

/* cooling functions for the heavy elements */
void CoolAlum(void);

void CoolArgo(void);

void CoolCalc(void);

void CoolCarb(void);

void CoolChlo(void);

void CoolChro(void);

void CoolCoba(void);

void CoolDima(void);

void CoolFluo(void);

void CoolIron(void);

void CoolMagn(void);

void CoolMang(void);

void CoolNeon(void);

void CoolNick(void);

void CoolNitr(void);

void CoolOxyg(void);

void CoolPhos(void);

void CoolPota(void);

void CoolScan(void);

void CoolSili(void);

void CoolSodi(void);

void CoolSulf(void);

void CoolTita(void);

void CoolVana(void);

void CoolZinc(void);

void oi_cs(double& cs3P23P1,
	double& cs3P23P0,
	double& cs3P13P0,
	double& cs3P1D2,
	double& cs1D21S0,
	double& cs3P1S0);

void oi_othercs(double& csh01,
	double& cshe01,
	double& csh201,
	double& csh12,
	double& cshe12,
	double& csh212,
	double& csh02,
	double& cshe02,
	double& csh202,
	double& csh2o01,
	double& csh2o02,
	double& csh2o12,
	double& csh2p01,
	double& csh2p02,
	double& csh2p12,
	double& csp01,
	double& csp02,
	double& csp12);

void oii_cs(double& oii_cs4S32D5,
	double& oii_cs4S32D3,
	double& oii_cs2D52D3,
	double& oii_cs4S32P3,
	double& oii_cs2D52P3,
	double& oii_cs2D32P3,
	double& oii_cs4S32P1,
	double& oii_cs2D52P1,
	double& oii_cs2D32P1,
	double& oii_cs2P32P1,
	double& oii_cs4S34P);

void oiii_cs(double& oiii_cs3P25S2,
	double& oiii_cs3P15S2,
	double& oiii_cs3P05S2,
	double& oiii_cs3P1D2,
	double& oiii_cs1D21S0,
	double& oiii_cs3P1S0,
	double& oiii_cs3P03P1,
	double& oiii_cs3P13P2,
	double& oiii_cs3P03P2,
	double& oiii_cs3P3D);

void oiv_cs(double& oiv_cs2P2D,double& oiv_cs2P12P3);

void ov_cs(double& ov_cs1S01P1,double& ov_cs1S03P);

void sii_cs(double& sii_cs4S32D3,
	double& sii_cs4S32D5,
	double& sii_cs4S32P1,
	double& sii_cs4S32P3,
	double& sii_cs2D32D5,
	double& sii_cs2D32P1,
	double& sii_cs2D52P1,
	double& sii_cs2D32P3,
	double& sii_cs2D52P3,
	double& sii_cs2P12P3,
	double& sii_cs4S34P);

void siii_cs(double& siii_cs3P03P1,
	double& siii_cs3P03P2,
	double& siii_cs3P01D2,
	double& siii_cs3P01S0,
	double& siii_cs3P13P2,
	double& siii_cs3P11D2,
	double& siii_cs3P11S0,
	double& siii_cs3P21D2,
	double& siii_cs3P21S0,
	double& siii_cs1D21S0,
	double& siii_cs3P3D,
	double& siii_cs3P5S2);

void siv_cs(double& siv_cs2P12P3);

void sviii_cs(double& sviii_cs2P32P1);

void neiii_cs(double& neiii_cs3P13P0,
		double& neiii_cs3P23P1,
		double& neiii_cs3P23P0,
		double& neiii_cs3P1D2,
		double& neiii_cs3P1S0);

#endif /* COOLING_H_ */
