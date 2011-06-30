/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef RFIELD_H_
#define RFIELD_H_

/* rfield.h */
#include "energy.h"

/** wavelength of V filter in Angstroms */
const double WL_V_FILT = 5500.;

/** wavelength of B filter in Angstroms */
const double WL_B_FILT = 4400.;

/** parameters to do with incident continuum */
/** limit to number of spectra that can be entered */
const int LIMSPC = 100;

/** the limit to the size of the coarse continuum array */
const int NCELL = 130000;

/** zero out rfield arrays between certain limits, code in zero.c */
void rfield_opac_zero( long lo , long ihi );

/** set true when malloced, init to false */
extern bool lgRfieldMalloced;

namespace Illuminate {
	typedef enum { FORWARD , REVERSE , SYMMETRIC } IlluminationType ;
}

EXTERN struct t_rfield {

	/** ================================================================================= */
	/** the following define the continuum energy scale and its limits */

	/** nflux is number of continuum points needed to get to high energy
	 * end of this continuum.  this is initially set to nupper in ContCreatePointers,
	 * and then trimmed down in ContSetIntensity.  It is finally reset in ConvInitSolution,
	 * to make sure continuum includes all possible line and continuum emission
	 * produced by the code */
	long int nflux;

	/** number of frequency cells needed to get full energy resolution of code */
	long int nupper;

	/** the energy of the lower limit low-energy limit of the continuum */
	realnum emm;

	/** the energy of the upper limit high-energy limit of the continuum */
	realnum egamry;

	/**faintest high energy flux to consider, set with set flxfnt command */
	realnum FluxFaint;

	/** energy in Ryd of center of cell */
	realnum *anu;

	/** original set of energy units, always the same no matter what the
	 * shape of the continuum - actual energy units change with shape */
	realnum *AnuOrg;

	/** width of cells in Rydberg */
	realnum *widflx;

	/** used to keep track of number of lines per freq interval */
	long int *line_count;

	/** bremsstrahlung occupation number - only used for induced two photon */
	realnum *OccNumbBremsCont;

	/** outward emitted continuum */
	realnum *OccNumbContEmitOut;

	/** these are log, sqrt, square, and cube of anu array */
	realnum *anulog, 
	  *anusqr, 
	  *anu2, 
	  *anu3;

	/** ================================================================================= */
	/** the following are the arrays containing the local radiation field */

	/**flux is photons per cell N.B. width of cells vary with energy, given by widflx */
	realnum **flux;

	/** this is the isotropic part of the constant continuum */
	realnum *flux_isotropic;

	/** this is the variable and constant parts of the above */
	realnum *flux_beam_time , *flux_beam_const;

	/** the accumulated flux, sum from this energy to infinity */
	realnum *flux_accum;

	/** extinction factor set with extinguish command */
	realnum *ExtinguishFactor;
	realnum ExtinguishLeakage, 
	 ExtinguishColumnDensity, 
	 ExtinguishLowEnergyLimit,
	 /** the constant that multiplies the column density to get optical depth */
	 ExtinguishConvertColDen2OptDepth ,
	 /** the power on the energy */
	 ExtinguishEnergyPowerLow;

	/** this is set true is one of incident continua is expected to have
	 * all H-ionizing radiation blocked.  That is done separately with
	 * the extinguish command.  */
	bool lgMustBlockHIon;

	/* this is set true if H-ionizing radiation is blocked with extinguish
	 * command */
	bool lgBlockHIon;

	/** option to not do line transfer, set false with no line transfer command */
	bool lgDoLineTrans;

	/** says whether to constantly reevaluate opacities, normally true,
	 * set false with no opacity reevaluate command */
	bool lgOpacityReevaluate;

	/** this flag says that CMB has been set */
	bool lgCMB_set;

	/** says whether to constantly reevaluate ionization, normally true,
	 * set false with no ionization reevaluate command */
	bool lgIonizReevaluate;

	/** convoc is the conversion factor from rfield to OccNumbIncidCont */
	realnum *convoc;

	/** OccNumbIncidCont is the continuum occupation number for the 
	 * attenuated incident ONLY */
	realnum *OccNumbIncidCont;

	/** OccNumbDiffCont is the continuum occupation number, for local diffuse continuum */
	realnum *OccNumbDiffCont;

	/** array of Boltzmann factors for the continuum energy grid and current 
	 * temperature */
	double *ContBoltz;

	/** ConEmitLocal is local diffuse continuum, units photons cm-3 s-1 cell-1, 
	 * for each zone, evaluated in RT_diffuse */
	realnum **ConEmitLocal/* [depth][energy]*/;

	/** the local source function - diffuse emission, photons cell-1 cm-2 s-1 */
	realnum **ConSourceFcnLocal/* [depth][energy]*/;

	/** reflected diffuse emission continuum */
	realnum **ConEmitReflec;

	/** outward diffuse emission continuum (not the interactive one),
	 * this is incremented in radinc because of interplay between absorption
	 * and emission - get the outward bremsstrahlung right 
	 * photons cell-1 cm-2 s-1 */
	realnum **ConEmitOut;

	/** this is set in RT_diffuse and carries interactive continua */
	realnum *ConInterOut;

	/** ConRefIncid is reflected portion of incident continuum */
	realnum **ConRefIncid;

	/** these are energy-by-energy sums of various arrays, used to save time in
	 * evaluating rate integrals */
	realnum *SummedCon;
	realnum *SummedDif;
	realnum *SummedOcc;
	realnum *SummedDifSave;

	/** this will control array of locally destroyed continuum photons,
	 * zeroed and evaluated in RT_diffuse, currently only two-photon */
	realnum *ConOTS_local_photons,
		/** the local photoionization rate corresponding to above photons */
		*ConOTS_local_OTS_rate;

	/** computed in rt_diffuse.cpp, escaping continuum emission 
	 * added to beam in radius_increment */
	realnum *DiffuseEscape;

	/** saves total two photon continuum for debugging, set in RT_diffuse */
	realnum *TotDiff2Pht;

	/** otsline and otscon - local ots fields for line and continua
	 * outlin outward line fields */
	realnum 
		/** the local ots line rates */
		*otslin/*[NC_ELL]*/, 
		/** the local ots continuum rates */
		*otscon/*[NC_ELL]*/, 
		**otssav/*[NC_ELL][2]*/;

	/** outward directed line emission photons cm-2 s-1 */
	realnum **outlin,
		*outlin_noplot;

	/** local diffuse line emission, photons cm-3 s-1 */
	realnum *DiffuseLineEmission;

	/** reflected line */
	realnum **reflin/*[NC_ELL]*/;

	/** save incident continuum for later iterations */
	realnum **flux_total_incident;
	realnum *flux_beam_const_save , *flux_time_beam_save , *flux_isotropic_save;

	/** ==1 for time steady, when continuum varies with time, is scale factor */
	realnum time_continuum_scale;

	/** this is zero or one depending whether pumping by diffuse fields
	 * is turned on (1) or turned off (0) */
	realnum DiffPumpOn;

	/** string identifying the first line that occurred at this energy */
	char **chLineLabel/*[NC_ELL][5]*/;

	/** string identifying the first continuum edge that occurred at this energy */
	char **chContLabel/*[NC_ELL][5]*/;

	/** free free gaunt factor for all charges */
	/** First index is EFFECTIVE CHARGE!	*/
	realnum **gff/*[LIMELM][NC_ELL]*/;

	/** flag which, if set to false, causes gauntff.dat, if it exists, to be read in,
	 * or if true, causes the file gauntff.dat to be created. rfield.gff[][] is
	 * filled by interpolation on the values in this file	*/
	bool lgCompileGauntFF;

	/** method for transferring diffuse continuum
	 * either 'ots' or 'oux' where x is n for which */
	char chDffTrns[4];

	/** another flag, true if outward only - 
	 * used to multiply the ConInterOut continuum for creating the interactive
	 * continuum - when false, no outward only, this is not added */
	bool lgOutOnly;

	/** ipEnergyBremsThin is index for lowest energy thin to ff abs and plasma frequency
	 * EnergyBremsThin is energy there, Ryd */
	long int ipEnergyBremsThin;
	realnum EnergyBremsThin;

	/** index of highest cell with positive Boltzmann factor */
	long int ipMaxBolt;

	/** turn off continuum pumping, set with 'no induced processes' command */
	bool lgInducProcess;

	/** saves for the upward and downward Compton shifts */
	double *comup, 
	  *comdn;

	/** array indices for centers of B and V filters */
	long int ipB_filter , ipV_filter;

	/** these are the lower and upper bounds for the G0 radiation field
	 * used by Tielens & Hollenbach in their PDR work */
	long int ipG0_TH85_lo , ipG0_TH85_hi;

	/** these are the lower and upper bounds for the G0 radiation field
	 * used by Tielens & Hollenbach in their PDR work */
	long int ipG0_DB96_lo , ipG0_DB96_hi;

	/** these are the lower and upper bounds for the special G0 radiation field */
	long int ipG0_spec_lo , ipG0_spec_hi;

	/** this is the wavelength where Bertoldi & Draine estimate the Habing field */
	long int ip1000A;

	/** extinction in magnitudes at B and V filters for a point source,
	 * this does not discount forward scattering */
	double extin_mag_B_point , extin_mag_V_point;

	/** extinction in magnitudes at B and V filters for a resolved source,
	 * this does discount forward scattering, so is appropriate for a resolved source */
	double extin_mag_B_extended , extin_mag_V_extended;

	/** these are total opacities at these wavelengths, used to stop at exact Av */
	double opac_mag_B_point, opac_mag_V_point, opac_mag_B_extended , opac_mag_V_extended;

	/** coefficients for fitting Tarter expressions for Compton 
	 * heating and cooling, over full energy range of code */
	realnum *csigh, 
	  *csigc;

	double comtot, 
	  cmheat, 
	  cmcool, 
	  comoff, 
	  cinrat;

	/** set true if Compton cooling underflows */
	bool lgComUndr;

	double totpow[LIMSPC], 
	  slope[LIMSPC], 
	  cutoff[LIMSPC][3], 
	  spfac[LIMSPC];

	/** option to have continuum intensity be time dependent */
	bool lgTimeVary[LIMSPC];

	/* beamed or isotropic continuum?  if isotropic then does not vary
	 * with time */
	bool lgBeamed[LIMSPC];

	/** 1 / cos( illumination angle, angle measured from normal,
	 * default is angle=zero, normal illumination, DirectCos = 1 */
	realnum OpticalDepthScaleFactor[LIMSPC];

	Illuminate::IlluminationType Illumination[LIMSPC];

	/** nShape is SED shape index number, this must equal the number
	 * of field intensities that are specified
	 * ipSpec is radiation field source number 
	 */
	long int nShape, 
	  ipSpec;

	/** these are used for interpolate and table commands, 
	 * all have two indices, continuum and frequency
	 * tNuRyd is the linear energy Rydberg of continuum point */
	/** must be explicit arrays again so that
	 * table commands will work before continuum defined.*/
	vector<Energy> tNu[LIMSPC];
	vector<realnum> tslop[LIMSPC];
	/**< this is the log of f_nu of continuum point */
	vector<realnum> tFluxLog[LIMSPC];
	
	long ncont[LIMSPC];

	/** used to store a check on the continuum mesh resolution scale factor
	 * this is used for the output of the SAVE TRANSMITTED CONTINUUM command */
	double RSFCheck[LIMSPC];

	/** this flag is set true if we malloced out the three previous arrays -
	 * these can be returned once the continuum is generated */
	bool lgContMalloc[LIMSPC];

	/** energy range over which the intensity is integrated for normalizing
	 * each continuum source that contributes to the total source */
	double range[LIMSPC][2];

	/** chSpNorm says how spectrum was normalized - intensity or luminosity case
	 * chRSpec says whether per unit area or tot 4 pi  */
	char chSpNorm[LIMSPC][5], 
	  chRSpec[LIMSPC][5], 
	  chSpType[LIMSPC][6];

	/** these are total numbers of photons over various energy ranges */
	realnum qhtot, 
	  qhe, 
	  qheii, 
	  qbal, 
	  qrad, 
	  qtot;

	/** hydrogen ionization parameter */
	realnum uh;

	/** helium ion ionization parameter */
	realnum uheii;

	/** lgUSphON flag set when we hit Stromgren radius in spherical geometry	*/
	bool lgUSphON;
	/** the Stromgren radius var set to get u spherical */
	realnum rstrom;

	/** flag set if incident radiation field less than
	 * 10x the Habing ISM field */
	bool lgHabing;

	/** heaviest element to be considered - default is iron - the fine opacity
	 * array resolution depends on heaviest element and lowest temperature */
	long int fine_opac_nelem;

	/** number of resolution elements over width of this element, default is 4 */
	long int fine_opac_nresolv;

	/** number of cm/s of each cell in fine mesh */
	realnum fine_opac_velocity_width;

	/** these are transmission coefficients for fine opacity degraded to coarse continuum 
	 * for the last zone */
	realnum *trans_coef_zone;
	/** these are TOTAL transmission coefficients for fine opacity degraded to coarse continuum */
	realnum *trans_coef_total;

	/** and array indices for lower and upper bounds of each coarse continuum mapped onto the
	 * fine continuum 0 (false) if fine continuum does not extend to these cells */
	long int *ipnt_coarse_2_fine;

	/** low and high bounds of fine continuum - set by need to include all resonance lines */
	realnum fine_ener_lo, fine_ener_hi;
	/** the number of fine continuum cells actually used in this continuum */
	long nfine;
	/** the number of fine continuum cells malloced - will be several million */
	long nfine_malloc;
	/** the dimensionless resolution of the fine continuum - dE/E */
	double fine_resol;
	/** the fine continuum opacity array  */
	realnum *fine_opac_zone;
	/** total optical depth array for fine continuum */
	realnum *fine_opt_depth;
	/** energies at center of each bin for fine continuum */
	realnum *fine_anu;

	/** shift in fine continuum rest frame scale due to velocity gradient, 
	 * rest frame is velocity of first zone, positive means that first zone 
	 * is blue shifted relative to current zone.  This is a decelerating
	 * flow.  Evaluated in RT_line_all */
	long int ipFineConVelShift;

	/** option to turn off fine opacities with no fine opacity command */
	bool lgOpacityFine;
	/** says that fine optical depths will be saved, so save them */
	bool lgSaveOpacityFine;

	/** flag saying whether to include H1 Lya ots in the radiation field - 
	 * usually true, set false with no Lya ots command */
	bool lgLyaOTS;

	/** flag saying whether to include HeII Lya and rec cont ots in the 
	 * radiation field - - usually true, set false with no HeII ots command */
	bool lgHeIIOTS;

	/** flag saying whether to kill outward only lines */
	bool lgKillOutLine;

	/** flag saying whether to kill outward only continuum */
	bool lgKillOutCont;

	/** flag saying whether to kill ots lines */
	bool lgKillOTSLine;

	/** following deal with plasma frequency, which enters the continuum
	 * array for even moderate densities due to very low frequencies considered */
	/** set true if plasma freq enters energy array */
	bool lgPlasNu;
	/** plasma frequency for current position in slab */
	realnum plsfrq, 
		/** store highest energy plasma frequency encountered*/
		plsfrqmax;
	// the zone where the plasma frequency is evaluated 
	long int nZonePlsFrqEval;
	/** pointer to current plasma frequency */
	long int ipPlasma, 
		/** pointer to largest plasma freq encountered */
	  ipPlasmax;

	/** these are series of flags that say whether different parts of the
	 * continuum where entered ok */
	bool lgMMok, 
	  lgHPhtOK, 
	  lgXRayOK, 
	  lgGamrOK;

	/** lowest energy counted as gamma rays, 100 keV	*/
	realnum EnerGammaRay;
	long int ipEnerGammaRay;

	/** lgHionRad set to .true. if no hydrogen ionizing radiation */
	bool lgHionRad;

	/** these store photon occupation numbers at various energies in the continuum */
	realnum occmax, 
	  occmnu, 
	  tbrmax, 
	  tbrmnu, 
	  tbr4nu, 
	  occ1nu;

	/** flag saying that occupation number at 1 Ryd > 1 */
	bool lgOcc1Hi;

	/** intensity [erg cm-2 s-1] in incident and diffuse continua in
	 * current zone */
	realnum EnergyIncidCont ,
		EnergyDiffCont;

	// constructor
	t_rfield( )
	{
		nZonePlsFrqEval = -1;

		// the constant that multiplies the column density to get optical depth at 1 Ryd
		ExtinguishConvertColDen2OptDepth = (realnum)6.22e-18;
		// the power on the energy for the extinction
		ExtinguishEnergyPowerLow = (realnum)-2.43;

		// these are the low and high energy bounds of the continuum
		emm = 1.001e-8f;
		egamry = 7.354e6f;

		// space not created yet
		for( long i=0; i < LIMSPC; i++ )
			lgContMalloc[i] = false;
	}

}	rfield;


#endif /* RFIELD_H_ */
