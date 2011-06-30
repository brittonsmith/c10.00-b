/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "physconst.h"
#include "optimize.h"
#include "continuum.h"
#include "called.h"
#include "rfield.h"
#include "thirdparty.h"
#include "stars.h"
/*lint -e785 too few initializers */
/*lint -e801 use of go to depreciated */

/** this is the initial assumed size of the Starburst grid, may be increased during execution if needed */
static const int NSB99 = 1250;
/** maximum number of separate time steps in a Starburst99 model */
static const int MNTS = 200;

/** this is the number of points in each of the stellar continua */
static const int NRAUCH = 19951;
/** The number of models in the original Rauch H-Ca set (first version May 1998, current May 2001) */
static const int NMODS_HCA = 66;
/** The number of models in the new Rauch H-Ni set, Nov 2002 */
static const int NMODS_HNI = 51;
/** The number of models in the new Rauch PG1159 set, Jan 2006 */
static const int NMODS_PG1159 = 71;
/** The number of models in the Rauch Hydrogen only set, Feb 2003 */
static const int NMODS_HYDR = 100;
/** The number of models in the Rauch Helium only set, Jun 2004 */
static const int NMODS_HELIUM = 81;
/** The number of models in the Rauch H+He set, Aug 2004 */
static const int NMODS_HpHE = 117;

/* set to 1 to turn on debug print statements in these routines */
#define DEBUGPRT 0

#define FREE_CHECK(PTR) { ASSERT( PTR != NULL ); free( PTR ); PTR = NULL; }
#define FREE_SAFE(PTR) { if( PTR != NULL ) free( PTR ); PTR = NULL; }

static const bool lgSILENT = false;
static const bool lgVERBOSE = true;

static const bool lgLINEAR = false;
static const bool lgTAKELOG = true;

typedef enum {
	IS_UNDEFINED, IS_FIRST, IS_SECOND
} IntStage;

/** store the parameters of a single atmosphere model */
typedef struct
{
	double par[MDIM];
	int modid;
	char chGrid;
} mpp;

/** \todo - check rebinning of Tlusty models
 ** \todo - why was it necessary to change stars_tlusty.in? (change from r43 to r50?)
 ** \todo - check all interpolation modes of CoStar
 ** \todo - compare models with original code, dump atmospheres!
 ** \todo - check all Edges arrays...
 ** \todo - update Doxygen documentation
 */

/* this is the structure of the binary atmosphere file (VERSION 20100902[01]):
 *
 *               ============================
 *               * int32 VERSION            *
 *               * int32 MDIM               *
 *               * int32 MNAM               *
 *               * int32 ndim               *
 *               * int32 npar               *
 *               * int32 nmods              *
 *               * int32 ngrid              *
 *               * uint32 nOffset           *
 *               * uint32 nBlocksize        *
 *               * double mesh_elo          *
 *               * double mesh_ehi          *
 *               * double mesh_res_factor   *
 *               * char md5sum[NMD5]        *
 *               * char names[MDIM][MNAM+1] *
 *               * mpp telg[nmods]          *
 *               * realnum anu[ngrid]       *
 *               * realnum mod1[ngrid]      *
 *               *    ...                   *
 *               * realnum modn[ngrid]      *
 *               ============================
 *
 * nOffset == 7*sizeof(int32) + 2*sizeof(uint32) + 3*sizeof(double) +
 *            (NMD5 + MDIM*(MNAM+1))*sizeof(char) + nmods*sizeof(mpp)
 * nBlocksize == ngrid*size(realnum) */

/** store all the relevant information on a binary atmosphere file */
typedef struct
{
	/** the name of the binary atmosphere file */
	string name;
	/** if true, more relaxed rules for matching log(g) will be used */
	bool lgIsTeffLoggGrid;
	/** where should we search for the binary atmosphere file */
	access_scheme scheme;
	/** the file handle for this file */
	FILE *ioIN;
	/** the identifier for this grid used in the Cloudy output,
	 * this *must* be exactly 12 characters long */
	const char *ident;
	/** the Cloudy command to recompile the binary atmosphere file */
	const char *command;
	/** which interpolation mode is requested */
	IntMode imode;
	/** the number of dimensions in the grid */
	int32 ndim;
	/** the number of parameters for each model; npar >= ndim */
	int32 npar;
	/** the number of stellar atmosphere models in this file */
	int32 nmods;
	/** the number of grid points per model, should equal rfield.nupper */
	int32 ngrid;
	/** the offset to the first data block (the anu grid) */
	uint32 nOffset;
	/** the size of each model block in bytes */
	uint32 nBlocksize;
	/** these are the model parameters in the same
	 * sequence they are stored in the binary file */
	mpp *telg;    /* telg[nmods] */
	/** these are the unique values for each of the model parameters */
	double **val; /* val[ndim][nval[n]] */
	/** nval[n] is the number of unique values in val[n][*] */
	long *nval;   /* nval[ndim] */
	/** jlo/jhi will hold indices into the binary model file: jlo/jhi(i,...,n)
	 * will point to the model with parameters val[0][i],...,val[ndim-1][n],
	 * or its closest approximation in log(g) in case the model doesn't exist
	 * and lgIsTeffLoggGrid is true.
	 * jlo will hold the model with the highest log(g) <= than requested
	 * jhi will hold the model with the lowest log(g) >= than requested
	 * in case no suitable model could be found either array will hold -2 */
	long *jlo;   /* jlo(nval[0],...,nval[ndim-1]) */
	long *jhi;   /* jhi(nval[0],...,nval[ndim-1]) */
	/** this array will hold the designations for each dimension of the grid */
	char names[MDIM][MNAM+1];
	/** this array holds the length of each CoStar track */
	long *trackLen; /* trackLen[nTracks] */
	/** this is the number of CoStar tracks */
	long nTracks;
	/** jval will hold indices into the CoStar grid: jval(nModels,nTracks) */
	long *jval;
} stellar_grid;

/* internal routines */
STATIC bool lgCompileAtmosphereCoStar(const char[],const char[],const realnum[],long,process_counter&);
STATIC void InterpolateGridCoStar(const stellar_grid*,const double[],double*,double*);
STATIC void FindHCoStar(const stellar_grid*,long,double,long,realnum*,long*,long*);
STATIC void FindVCoStar(const stellar_grid*,double,realnum*,long[]);
STATIC void CoStarListModels(const stellar_grid*);
STATIC int RauchInitializeSub(const char[],const char[],const vector<mpp>&,long,long,
			       long,const double[],int);
STATIC void RauchReadMPP(vector<mpp>&,vector<mpp>&,vector<mpp>&,vector<mpp>&,vector<mpp>&,vector<mpp>&);
inline void getdataline(fstream&,string&);
STATIC bool lgCompileAtmosphere(const char[],const char[],const realnum[],long,process_counter&);
STATIC void InitGrid(stellar_grid*,bool);
STATIC bool lgValidBinFile(const char*,process_counter&,access_scheme);
STATIC bool lgValidAsciiFile(const char*,access_scheme);
STATIC void InitGridCoStar(stellar_grid*);
STATIC void CheckVal(const stellar_grid*,double[],long*,long*);
STATIC void InterpolateRectGrid(const stellar_grid*,const double[],double*,double*);
STATIC void FreeGrid(stellar_grid*);
STATIC void InterpolateModel(const stellar_grid*,const double[],double[],const long[],
			     const long[],long[],long,vector<realnum>&,IntStage);
STATIC void InterpolateModelCoStar(const stellar_grid*,const double[],double[],const long[],
				   const long[],long[],long,long,vector<realnum>&);
STATIC void GetBins(const stellar_grid*,vector<Energy>&);
STATIC void GetModel(const stellar_grid*,long,vector<realnum>&,bool,bool);
STATIC void SetLimits(const stellar_grid*,double,const long[],const long[],const long[],
		      const realnum[],double*,double*);
STATIC void SetLimitsSub(const stellar_grid*,double,const long[],const long[],long[],long,
			 double*,double*);
STATIC void InitIndexArrays(stellar_grid*,bool);
STATIC void FillJ(const stellar_grid*,long[],double[],long,bool);
STATIC long JIndex(const stellar_grid*,const long[]);
STATIC void SearchModel(const mpp[],bool,long,const double[],long,long*,long*);
STATIC void FindIndex(const double[],long,double,long*,long*,bool*);
STATIC bool lgFileReadable(const char*, process_counter&,access_scheme);
STATIC void ValidateGrid(const stellar_grid*,double);
STATIC bool lgValidModel(const vector<Energy>&,const vector<realnum>&,double,double);
STATIC void RebinAtmosphere(long,const realnum[],const realnum[],realnum[],long,const realnum[]);
STATIC realnum RebinSingleCell(realnum,realnum,const realnum[],const realnum[],const realnum[],long);
STATIC long RebinFind(const realnum[],long,realnum);


/* the version number for the ascii/binary atmosphere files */
static const long int VERSION_ASCII = 20060612L;
/* binary files are incompatible when floats are converted to doubles */
#ifdef FLT_IS_DBL
static const long int VERSION_BIN = 201009020L;
#else
static const long int VERSION_BIN = 201009021L;
#endif
static const long int VERSION_RAUCH_MPP = 20090324;

/** List all the available TABLE STAR <grid> commands by checking installed *.mod files */
void AtmospheresAvail( void )
{
	DEBUG_ENTRY( "AtmospheresAvail()" );

	/* This routine makes a list of all the stellar atmosphere grids that are valid,
	 * giving the parameters for use in the input script as well. It is simply a long
	 * list of if-statements, so if any grid is added to Cloudy, it should be added in
	 * this routine as well.
	 *
	 * NB NB NB -- test this routine regularly to see if the list is still complete! */

	fprintf( ioQQQ, "\n I will now list all stellar atmosphere grids that are ready to be used (if any).\n" );
	fprintf( ioQQQ, " User-defined stellar atmosphere grids will not be included in this list.\n\n" );

	process_counter dum;

	/* we always look in the data directory regardless of where we are,
	 * it would be very confusing to the user if we did otherwise... */
	access_scheme as = AS_DATA_ONLY_TRY;

	if( lgValidBinFile( "atlas_fp10k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z+1.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fp05k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z+0.5 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fp03k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z+0.3 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fp02k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z+0.2 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fp01k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z+0.1 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fp00k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z+0.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm01k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-0.1 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm02k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-0.2 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm03k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-0.3 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm05k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-0.5 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm10k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-1.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm15k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-1.5 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm20k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-2.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm25k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-2.5 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm30k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-3.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm35k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-3.5 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm40k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-4.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm45k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-4.5 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm50k2.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas Z-5.0 <Teff> [ <log(g)> ]\n" );

	if( lgValidBinFile( "atlas_fp05k2_odfnew.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas odfnew Z+0.5 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fp02k2_odfnew.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas odfnew Z+0.2 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fp00k2_odfnew.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas odfnew Z+0.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm05k2_odfnew.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas odfnew Z-0.5 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm10k2_odfnew.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas odfnew Z-1.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm15k2_odfnew.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas odfnew Z-1.5 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm20k2_odfnew.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas odfnew Z-2.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "atlas_fm25k2_odfnew.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas odfnew Z-2.5 <Teff> [ <log(g)> ]\n" );

	if( lgValidBinFile( "atlas_3d.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas 3-dim <Teff> <log(g)> <log(Z)>\n" );

	if( lgValidBinFile( "atlas_3d_odfnew.mod", dum, as ) )
		fprintf( ioQQQ, "   table star atlas odfnew 3-dim <Teff> <log(g)> <log(Z)>\n" );

	if( lgValidBinFile( "Sc1_costar_solar.mod", dum, as ) )
		fprintf( ioQQQ, "   table star costar solar (see Hazy for parameters)\n" );
	if( lgValidBinFile( "Sc1_costar_halo.mod", dum, as ) )
		fprintf( ioQQQ, "   table star costar halo (see Hazy for parameters)\n" );

	if( lgValidBinFile( "kurucz79.mod", dum, as ) )
		fprintf( ioQQQ, "   table star kurucz79 <Teff>\n" );

	if( lgValidBinFile( "mihalas.mod", dum, as ) )
		fprintf( ioQQQ, "   table star mihalas <Teff>\n" );

	if( lgValidBinFile( "rauch_h-ca_solar.mod", dum, as ) )
		fprintf( ioQQQ, "   table star rauch H-Ca solar <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "rauch_h-ca_halo.mod", dum, as ) )
		fprintf( ioQQQ, "   table star rauch H-Ca halo <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "rauch_h-ca_3d.mod", dum, as ) )
		fprintf( ioQQQ, "   table star rauch H-Ca 3-dim <Teff> <log(g)> <log(Z)>\n" );

	if( lgValidBinFile( "rauch_h-ni_solar.mod", dum, as ) )
		fprintf( ioQQQ, "   table star rauch H-Ni solar <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "rauch_h-ni_halo.mod", dum, as ) )
		fprintf( ioQQQ, "   table star rauch H-Ni halo <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "rauch_h-ni_3d.mod", dum, as ) )
		fprintf( ioQQQ, "   table star rauch H-Ni 3-dim <Teff> <log(g)> <log(Z)>\n" );

	if( lgValidBinFile( "rauch_pg1159.mod", dum, as ) )
		fprintf( ioQQQ, "   table star rauch pg1159 <Teff> [ <log(g)> ]\n" );

	if( lgValidBinFile( "rauch_hydr.mod", dum, as ) )
		fprintf( ioQQQ, "   table star rauch hydrogen <Teff> [ <log(g)> ]\n" );

	if( lgValidBinFile( "rauch_helium.mod", dum, as ) )
		fprintf( ioQQQ, "   table star rauch helium <Teff> [ <log(g)> ]\n" );

	if( lgValidBinFile( "rauch_h+he_3d.mod", dum, as ) )
		fprintf( ioQQQ, "   table star rauch H+He <Teff> <log(g)> <frac(He)>\n" );

	if( lgValidBinFile( "starburst99.mod", dum, as ) )
		fprintf( ioQQQ, "   table star \"starburst99.mod\" <age>\n" );
	if( lgValidBinFile( "starburst99_2d.mod", dum, as ) )
		fprintf( ioQQQ, "   table star \"starburst99_2d.mod\" <age> <Z>\n" );

	if( lgValidBinFile( "obstar_merged_p03.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty OBstar Z+0.3 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "obstar_merged_p00.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty OBstar Z+0.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "obstar_merged_m03.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty OBstar Z-0.3 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "obstar_merged_m07.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty OBstar Z-0.7 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "obstar_merged_m10.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty OBstar Z-1.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "obstar_merged_m99.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty OBstar Z-inf <Teff> [ <log(g)> ]\n" );

	if( lgValidBinFile( "obstar_merged_3d.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty OBstar 3-dim <Teff> <log(g)> <log(Z)>\n" );

	if( lgValidBinFile( "bstar2006_p03.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Bstar Z+0.3 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "bstar2006_p00.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Bstar Z+0.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "bstar2006_m03.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Bstar Z-0.3 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "bstar2006_m07.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Bstar Z-0.7 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "bstar2006_m10.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Bstar Z-1.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "bstar2006_m99.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Bstar Z-inf <Teff> [ <log(g)> ]\n" );

	if( lgValidBinFile( "bstar2006_3d.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Bstar 3-dim <Teff> <log(g)> <log(Z)>\n" );

	if( lgValidBinFile( "ostar2002_p03.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Ostar Z+0.3 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "ostar2002_p00.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Ostar Z+0.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "ostar2002_m03.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Ostar Z-0.3 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "ostar2002_m07.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Ostar Z-0.7 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "ostar2002_m10.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Ostar Z-1.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "ostar2002_m15.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Ostar Z-1.5 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "ostar2002_m17.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Ostar Z-1.7 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "ostar2002_m20.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Ostar Z-2.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "ostar2002_m30.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Ostar Z-3.0 <Teff> [ <log(g)> ]\n" );
	if( lgValidBinFile( "ostar2002_m99.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Ostar Z-inf <Teff> [ <log(g)> ]\n" );

	if( lgValidBinFile( "ostar2002_3d.mod", dum, as ) )
		fprintf( ioQQQ, "   table star tlusty Ostar 3-dim <Teff> <log(g)> <log(Z)>\n" );

	if( lgValidBinFile( "kwerner.mod", dum, as ) )
		fprintf( ioQQQ, "   table star werner <Teff> [ <log(g)> ]\n" );

	if( lgValidBinFile( "wmbasic.mod", dum, as ) )
		fprintf( ioQQQ, "   table star wmbasic <Teff> <log(g)> <log(Z)>\n" );
	return;
}

/* AtlasCompile rebin Kurucz stellar models to match energy grid of code */
/* >>chng 05 nov 16, added return value to indicate success (0) or failure (1) */
int AtlasCompile(process_counter& pc)
{
	/* these contain frequencies for the major absorption edges */
	realnum Edges[3];

	bool lgFail = false;

	DEBUG_ENTRY( "AtlasCompile()" );

	/* This is a program to re-bin the Kurucz stellar models spectrum to match the 
	 * CLOUDY grid.  For wavelengths shorter than supplied in the Kurucz files,
	 * the flux will be set to zero.  At long wavelengths a Rayleigh-Jeans
	 * extrapolation will be used. */

	/* This version uses power-law interpolation between the points of the stellar
	 * model.*/

	fprintf( ioQQQ, " AtlasCompile on the job.\n" );

	/* define the major absorption edges that require special attention during rebinning
	 *
	 * NB the frequencies should be chosen here such that they are somewhere in between
	 * the two frequency points that straddle the edge in the atmosphere model, the
	 * software in RebinAtmosphere will seek out the exact values of those two points
	 * e.g.: in the CoStar models the H I edge is straddled by wavelength points at
	 * 911.67 and 911.85 A, so Edges[0] should be chosen somewhere in between (e.g. at 911.76A).
	 *
	 * NB beware not to choose edges too close to one another (i.e. on the order of the
	 * resolution of the Cloudy frequency grid). E.g. the He II Balmer edge nearly coincides
	 * with the H I Ly edge, they should be treated as one edge. Trying to separate them will
	 * almost certainly lead to erroneous behaviour in RebinAtmosphere */
	Edges[0] = (realnum)(RYDLAM/911.76);
	Edges[1] = (realnum)(RYDLAM/504.26);
	Edges[2] = (realnum)(RYDLAM/227.84);

	access_scheme as = AS_LOCAL_ONLY_TRY;

	/* >>chng 05 nov 19, add support for non-solar metalicities as well as odfnew models, PvH */
	if( lgFileReadable( "atlas_fp10k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fp10k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fp10k2.ascii", "atlas_fp10k2.mod", Edges, 3L, pc );
	if( lgFileReadable( "atlas_fp05k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fp05k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fp05k2.ascii", "atlas_fp05k2.mod", Edges, 3L, pc );
	if( lgFileReadable( "atlas_fp03k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fp03k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fp03k2.ascii", "atlas_fp03k2.mod", Edges, 3L, pc );
	if( lgFileReadable( "atlas_fp02k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fp02k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fp02k2.ascii", "atlas_fp02k2.mod", Edges, 3L, pc );
	if( lgFileReadable( "atlas_fp01k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fp01k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fp01k2.ascii", "atlas_fp01k2.mod", Edges, 3L, pc );
	if( lgFileReadable( "atlas_fp00k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fp00k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fp00k2.ascii", "atlas_fp00k2.mod", Edges, 3L, pc );
	if( lgFileReadable( "atlas_fm01k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm01k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm01k2.ascii", "atlas_fm01k2.mod", Edges, 3L, pc );
	if( lgFileReadable( "atlas_fm02k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm02k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm02k2.ascii", "atlas_fm02k2.mod", Edges, 3L, pc );
	if( lgFileReadable( "atlas_fm03k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm03k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm03k2.ascii", "atlas_fm03k2.mod", Edges, 3L, pc );
	if( lgFileReadable( "atlas_fm05k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm05k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm05k2.ascii", "atlas_fm05k2.mod", Edges, 3L, pc );
	if( lgFileReadable( "atlas_fm10k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm10k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm10k2.ascii", "atlas_fm10k2.mod", Edges, 3L, pc );
	if( lgFileReadable( "atlas_fm15k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm15k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm15k2.ascii", "atlas_fm15k2.mod", Edges, 3L, pc );
	if( lgFileReadable( "atlas_fm20k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm20k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm20k2.ascii", "atlas_fm20k2.mod", Edges, 3L, pc );
	if( lgFileReadable( "atlas_fm25k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm25k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm25k2.ascii", "atlas_fm25k2.mod", Edges, 3L, pc );
	if( lgFileReadable( "atlas_fm30k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm30k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm30k2.ascii", "atlas_fm30k2.mod", Edges, 3L, pc );
	if( lgFileReadable( "atlas_fm35k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm35k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm35k2.ascii", "atlas_fm35k2.mod", Edges, 3L, pc );
	if( lgFileReadable( "atlas_fm40k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm40k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm40k2.ascii", "atlas_fm40k2.mod", Edges, 3L, pc );
	if( lgFileReadable( "atlas_fm45k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm45k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm45k2.ascii", "atlas_fm45k2.mod", Edges, 3L, pc );
	if( lgFileReadable( "atlas_fm50k2.ascii", pc, as ) && !lgValidBinFile( "atlas_fm50k2.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm50k2.ascii", "atlas_fm50k2.mod", Edges, 3L, pc );

	if( lgFileReadable( "atlas_fp05k2_odfnew.ascii", pc, as ) &&
	    !lgValidBinFile( "atlas_fp05k2_odfnew.mod", pc, as ) )
	    
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fp05k2_odfnew.ascii", "atlas_fp05k2_odfnew.mod",
							Edges, 3L, pc );
	if( lgFileReadable( "atlas_fp02k2_odfnew.ascii", pc, as ) &&
	    !lgValidBinFile( "atlas_fp02k2_odfnew.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fp02k2_odfnew.ascii", "atlas_fp02k2_odfnew.mod",
							Edges, 3L, pc );
	if( lgFileReadable( "atlas_fp00k2_odfnew.ascii", pc, as ) &&
	    !lgValidBinFile( "atlas_fp00k2_odfnew.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fp00k2_odfnew.ascii", "atlas_fp00k2_odfnew.mod",
							Edges, 3L, pc );
	if( lgFileReadable( "atlas_fm05k2_odfnew.ascii", pc, as ) &&
	    !lgValidBinFile( "atlas_fm05k2_odfnew.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm05k2_odfnew.ascii", "atlas_fm05k2_odfnew.mod",
							Edges, 3L, pc );
	if( lgFileReadable( "atlas_fm10k2_odfnew.ascii", pc, as ) &&
	    !lgValidBinFile( "atlas_fm10k2_odfnew.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm10k2_odfnew.ascii", "atlas_fm10k2_odfnew.mod",
							Edges, 3L, pc );
	if( lgFileReadable( "atlas_fm15k2_odfnew.ascii", pc, as ) &&
	    !lgValidBinFile( "atlas_fm15k2_odfnew.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm15k2_odfnew.ascii", "atlas_fm15k2_odfnew.mod",
							Edges, 3L, pc );
	if( lgFileReadable( "atlas_fm20k2_odfnew.ascii", pc, as ) &&
	    !lgValidBinFile( "atlas_fm20k2_odfnew.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm20k2_odfnew.ascii", "atlas_fm20k2_odfnew.mod",
							Edges, 3L, pc );
	if( lgFileReadable( "atlas_fm25k2_odfnew.ascii", pc, as ) &&
	    !lgValidBinFile( "atlas_fm25k2_odfnew.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_fm25k2_odfnew.ascii", "atlas_fm25k2_odfnew.mod",
							Edges, 3L, pc );

	if( lgFileReadable( "atlas_3d.ascii", pc, as ) && !lgValidBinFile( "atlas_3d.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_3d.ascii", "atlas_3d.mod", Edges, 3L, pc );

	if( lgFileReadable( "atlas_3d_odfnew.ascii", pc, as ) &&
	    !lgValidBinFile( "atlas_3d_odfnew.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "atlas_3d_odfnew.ascii", "atlas_3d_odfnew.mod", Edges, 3L, pc );
	return lgFail;
}

/* AtlasInterpolate read in and interpolate on Kurucz grid of atmospheres, originally by K Volk */
long AtlasInterpolate(double val[], /* val[nval] */
		      long *nval,
		      long *ndim,
		      const char *chMetalicity,
		      const char *chODFNew,
		      bool lgList,
		      double *Tlow,
		      double *Thigh)
{
	char chIdent[13];
	stellar_grid grid;

	DEBUG_ENTRY( "AtlasInterpolate()" );

	grid.name = "atlas_";
	if( *ndim == 3 )
		grid.name += "3d";
	else
	{
		grid.name += "f";
		grid.name += chMetalicity;
		grid.name += "k2";
	}
	grid.name += chODFNew;
	grid.name += ".mod";
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	if( *ndim == 3 )
	{
		strcpy( chIdent, "3-dim" );
	}
	else
	{
		strcpy( chIdent, "Z " );
		strcat( chIdent, chMetalicity );
	}
	strcat( chIdent, ( strlen(chODFNew) == 0 ? " Kurucz" : " ODFNew" ) );
	grid.ident = chIdent;
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	/* Note on the interpolation (solar abundance grid): 26 October 2000 (Peter van Hoof)
	 *
	 * I computed the effective temperature for a random sample of interpolated
	 * atmospheres by integrating the flux as shown above and compared the results
	 * with the expected effective temperature using DELTA = (COMP-EXPEC)/EXPEC.
	 *
	 * I found that the average discrepancy was:
	 *
	 *     DELTA = -0.10% +/- 0.06% (sample size 5000)
	 *
	 * The most extreme discrepancies were
	 *     -0.30% <= DELTA <= 0.21%
	 *
	 * The most negative discrepancies were for Teff =  36 -  39 kK, log(g) = 4.5 - 5
	 * The most positive discrepancies were for Teff = 3.5 - 4.0 kK, log(g) = 0 - 1
	 *
	 * The interpolation in the ATLAS grid is clearly very accurate */

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	FreeGrid( &grid );
	return rfield.nupper;
}

/* CoStarCompile rebin costar stellar models to match energy grid of code*/
int CoStarCompile(process_counter& pc)
{
	realnum Edges[3];
	bool lgFail = false;

	DEBUG_ENTRY( "CoStarCompile()" );

	fprintf( ioQQQ, " CoStarCompile on the job.\n" );

	/* define the major absorption edges that require special attention during rebinning
	 *
	 * NB the frequencies should be chosen here such that they are somewhere in between
	 * the two frequency points that straddle the edge in the atmosphere model, the
	 * software in RebinAtmosphere will seek out the exact values of those two points
	 * e.g.: in the CoStar models the H I edge is straddled by wavelength points at
	 * 911.67 and 911.85 A, so Edges[0] should be chosen somewhere in between (e.g. at 911.76A).
	 *
	 * NB beware not to choose edges too close to one another (i.e. on the order of the
	 * resolution of the Cloudy frequency grid). E.g. the He II Balmer edge nearly coincides
	 * with the H I Ly edge, they should be treated as one edge. Trying to separate them will
	 * almost certainly lead to erroneous behaviour in RebinAtmosphere */
	Edges[0] = (realnum)(RYDLAM/911.76);
	Edges[1] = (realnum)(RYDLAM/504.26);
	Edges[2] = (realnum)(RYDLAM/227.84);

	access_scheme as = AS_LOCAL_ONLY_TRY;

	if( lgFileReadable( "Sc1_costar_z020_lb.fluxes", pc, as ) && !lgValidBinFile( "Sc1_costar_solar.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphereCoStar( "Sc1_costar_z020_lb.fluxes", "Sc1_costar_solar.mod",
							      Edges, 3L, pc );
	if( lgFileReadable( "Sc1_costar_z004_lb.fluxes", pc, as ) && !lgValidBinFile( "Sc1_costar_halo.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphereCoStar( "Sc1_costar_z004_lb.fluxes", "Sc1_costar_halo.mod",
							      Edges, 3L, pc );
	return lgFail;
}

/* CoStarInterpolate read in and interpolate on CoStar grid of atmospheres */
long CoStarInterpolate(double val[], /* requested model parameters */
		       long *nval,
		       long *ndim,
		       IntMode imode, /* which interpolation mode is requested */
		       bool lgHalo,  /* flag indicating whether solar (==0) or halo (==1) abundances */
		       bool lgList,
		       double *val0_lo,
		       double *val0_hi)
{
	stellar_grid grid;

	DEBUG_ENTRY( "CoStarInterpolate()" );

	grid.name = ( lgHalo ? "Sc1_costar_halo.mod" : "Sc1_costar_solar.mod" );
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	grid.ident = "      costar";
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	/* listing the models in the grid is implemented in CoStarListModels() */
	InitGrid( &grid, false );
	/* now sort the models according to track */
	InitGridCoStar( &grid );
	/* override default interpolation mode */
	grid.imode = imode;

	if( lgList )
	{
		CoStarListModels( &grid );
		cdEXIT(EXIT_SUCCESS);
	}

	CheckVal( &grid, val, nval, ndim );

	/* Note on the interpolation: 26 October 2000 (Peter van Hoof)
	 *
	 * I computed the effective temperature for a random sample of interpolated
	 * atmospheres by integrating the flux as shown above and compared the results
	 * with the expected effective temperature using DELTA = (COMP-EXPEC)/EXPEC.
	 *
	 * I found that the average discrepancy was:
	 *
	 *     DELTA = -1.16% +/- 0.69% (SOLAR models, sample size 4590)
	 *     DELTA = -1.17% +/- 0.70% (HALO models, sample size 4828)
	 *
	 * The most extreme discrepancies for the SOLAR models were
	 *     -3.18% <= DELTA <= -0.16%
	 *
	 * The most negative discrepancies were for  Teff = 35 kK, log(g) = 3.5
	 * The least negative discrepancies were for Teff = 50 kK, log(g) = 4.1
	 *
	 * The most extreme discrepancies for the HALO models were
	 *     -2.90% <= DELTA <= -0.13%
	 *
	 * The most negative discrepancies were for  Teff = 35 kK, log(g) = 3.5
	 * The least negative discrepancies were for Teff = 50 kK, log(g) = 4.1
	 *
	 * Since Cloudy checks the scaling elsewhere there is no need to re-scale 
	 * things here, but this inaccuracy should be kept in mind since it could
	 * indicate problems with the flux distribution */

	InterpolateGridCoStar( &grid, val, val0_lo, val0_hi );

	FreeGrid( &grid );
	return rfield.nupper;
}

/* GridCompile rebin user supplied stellar models to match energy grid of code */
bool GridCompile(const char *InName)
{
	bool lgFail = false;
	realnum Edges[1];
	string OutName( InName );

	DEBUG_ENTRY( "GridCompile()" );

	fprintf( ioQQQ, " GridCompile on the job.\n" );

	// replace filename extension with ".mod"
	string::size_type ptr = OutName.find( '.' );
	ASSERT( ptr != string::npos );
	OutName.replace( ptr, string::npos, ".mod" );

	process_counter dum;
	lgFail = lgCompileAtmosphere( InName, OutName.c_str(), Edges, 0L, dum );

	if( !lgFail )
	{
		stellar_grid grid;

		/* the file must be local */
		grid.name = OutName;
		grid.scheme = AS_LOCAL_ONLY;
		grid.ident = "bogus ident.";
		grid.command = "bogus command.";

		InitGrid( &grid, false );

		/* check whether the models in the grid have the correct effective temperature */

		if( strcmp( grid.names[0], "Teff" ) == 0 )
		{
			fprintf( ioQQQ, " GridCompile: checking effective temperatures...\n" );
			ValidateGrid( &grid, 0.02 );
		}

		FreeGrid( &grid );
	}

	return lgFail;
}

/* GridInterpolate read in and interpolate on user supplied grid of atmospheres */
long GridInterpolate(double val[], /* val[nval] */
		     long *nval,
		     long *ndim,
		     const char *FileName,
		     bool lgList,
		     double *Tlow,
		     double *Thigh)
{
	char chIdent[13];
	stellar_grid grid;

	DEBUG_ENTRY( "GridInterpolate()" );

	// make filename without extension
	string chTruncName( FileName );
	string::size_type ptr = chTruncName.find( '.' );
	if( ptr != string::npos )
		chTruncName.replace( ptr, string::npos, "" );

	grid.name = FileName;
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	sprintf( chIdent, "%12.12s", chTruncName.c_str() );
	grid.ident = chIdent;
	/* the Cloudy command needed to recompile the binary model file */
	string chString( "COMPILE STARS \"" + chTruncName + ".ascii\"" );
	grid.command = chString.c_str();

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	FreeGrid( &grid );
	return rfield.nupper;
}

/* Kurucz79Compile rebin Kurucz 1979 stellar models to match energy grid of code */
int Kurucz79Compile(process_counter& pc)
{
	realnum Edges[1];

	bool lgFail = false;

	DEBUG_ENTRY( "Kurucz79Compile()" );

	fprintf( ioQQQ, " Kurucz79Compile on the job.\n" );

	/* following atmospheres LTE from Kurucz 1979, Ap.J. Sup 40, 1. and
	 * Kurucz (1989) private communication, newer opacities */

	access_scheme as = AS_LOCAL_ONLY_TRY;

	if( lgFileReadable( "kurucz79.ascii", pc, as ) && !lgValidBinFile( "kurucz79.mod", pc, as ) )
		lgFail = lgCompileAtmosphere( "kurucz79.ascii", "kurucz79.mod", Edges, 0L, pc );
	return lgFail;
}

/* Kurucz79Interpolate read in and interpolate on Kurucz79 grid of atmospheres */
long Kurucz79Interpolate(double val[], /* val[nval] */
			 long *nval,
			 long *ndim,
			 bool lgList,
			 double *Tlow,
			 double *Thigh)
{
	stellar_grid grid;

	DEBUG_ENTRY( "Kurucz79Interpolate()" );

	grid.name = "kurucz79.mod";
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	grid.ident = " Kurucz 1979";
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	FreeGrid( &grid );
	return rfield.nupper;
}

/* MihalasCompile rebin Mihalas stellar models to match energy grid of code */
int MihalasCompile(process_counter& pc)
{
	realnum Edges[1];

	bool lgFail = false;

	DEBUG_ENTRY( "MihalasCompile()" );

	fprintf( ioQQQ, " MihalasCompile on the job.\n" );

	/* following atmospheres NLTE from Mihalas, NCAR-TN/STR-76 */

	access_scheme as = AS_LOCAL_ONLY_TRY;

	if( lgFileReadable( "mihalas.ascii", pc, as ) && !lgValidBinFile( "mihalas.mod", pc, as ) )
		lgFail = lgCompileAtmosphere( "mihalas.ascii", "mihalas.mod", Edges, 0L, pc );
	return lgFail;
}

/* MihalasInterpolate read in and interpolate on Mihalas grid of atmospheres */
long MihalasInterpolate(double val[], /* val[nval] */
			long *nval,
			long *ndim,
			bool lgList,
			double *Tlow,
			double *Thigh)
{
	stellar_grid grid;

	DEBUG_ENTRY( "MihalasInterpolate()" );

	grid.name = "mihalas.mod";
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	grid.ident = "     Mihalas";
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	FreeGrid( &grid );
	return rfield.nupper;
}

/* RauchCompile create ascii and mod files for Rauch atmospheres */
int RauchCompile(process_counter& pc)
{
	bool lgFail = false;

	/* these contain frequencies for the major absorption edges */
	realnum Edges[3];

	/* Before running this program issue the following command where the Rauch
	 * model atmosphere files are kept (0050000_50_solar_bin_0.1 and so on)
	 *
	 *   ls *solar_bin_0.1 > rauchmods.list
	 *
	 * and check to see that there are 66 lines in the file.
	 */

	vector<mpp> telg1(NMODS_HCA);
	vector<mpp> telg2(NMODS_HNI);
	vector<mpp> telg3(NMODS_PG1159);
	vector<mpp> telg4(NMODS_HYDR);
	vector<mpp> telg5(NMODS_HELIUM);
	vector<mpp> telg6(NMODS_HpHE);

	/* metalicities of the solar and halo grid */
	static const double par2[2] = { 0., -1. };

	/* Helium fraction by mass */
	static const double par3[11] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };

	DEBUG_ENTRY( "RauchCompile()" );

	fprintf( ioQQQ, " RauchCompile on the job.\n" );

	RauchReadMPP( telg1, telg2, telg3, telg4, telg5, telg6 );

	process_counter dum;
	access_scheme as = AS_LOCAL_ONLY_TRY;

	/* this is the H-Ca grid */
	if( lgFileReadable( "0050000_50_solar_bin_0.1", dum, as ) && !lgValidAsciiFile( "rauch_h-ca_solar.ascii", as ) )
	{
		fprintf( ioQQQ, " Creating rauch_h-ca_solar.ascii....\n" );
		lgFail = lgFail || RauchInitializeSub( "rauch_h-ca_solar.ascii", "_solar_bin_0.1", 
						       telg1, NMODS_HCA, 1, 1, par2, 1 );
	}

	if( lgFileReadable( "0050000_50_halo__bin_0.1", dum, as ) && !lgValidAsciiFile( "rauch_h-ca_halo.ascii", as ) )
	{
		fprintf( ioQQQ, " Creating rauch_h-ca_halo.ascii....\n" );
		lgFail = lgFail || RauchInitializeSub( "rauch_h-ca_halo.ascii", "_halo__bin_0.1",
						       telg1, NMODS_HCA, 1, 1, par2, 1 );
	}

	if( lgFileReadable( "0050000_50_solar_bin_0.1", dum, as ) &&
	    lgFileReadable( "0050000_50_halo__bin_0.1", dum, as ) &&
	    !lgValidAsciiFile( "rauch_h-ca_3d.ascii", as ) )
	{
		fprintf( ioQQQ, " Creating rauch_h-ca_3d.ascii....\n" );
		lgFail = lgFail || RauchInitializeSub( "rauch_h-ca_3d.ascii", "_solar_bin_0.1",
						       telg1, NMODS_HCA, 1, 2, par2, 1 );
		lgFail = lgFail || RauchInitializeSub( "rauch_h-ca_3d.ascii", "_halo__bin_0.1",
						       telg1, NMODS_HCA, 2, 2, par2, 1 );
	}

	/* this is the H-Ni grid */
	if( lgFileReadable( "0050000_50_solar_iron.bin_0.1", dum, as ) &&
	    !lgValidAsciiFile( "rauch_h-ni_solar.ascii", as ) )
	{
		fprintf( ioQQQ, " Creating rauch_h-ni_solar.ascii....\n" );
		lgFail = lgFail || RauchInitializeSub( "rauch_h-ni_solar.ascii", "_solar_iron.bin_0.1",
						       telg2, NMODS_HNI, 1, 1, par2, 1 );
	}

	if( lgFileReadable( "0050000_50_halo__iron.bin_0.1", dum, as ) &&
	    !lgValidAsciiFile( "rauch_h-ni_halo.ascii", as ) )
	{
		fprintf( ioQQQ, " Creating rauch_h-ni_halo.ascii....\n" );
		lgFail = lgFail || RauchInitializeSub( "rauch_h-ni_halo.ascii", "_halo__iron.bin_0.1",
						       telg2, NMODS_HNI, 1, 1, par2, 1 );
	}

	if( lgFileReadable( "0050000_50_solar_iron.bin_0.1", dum, as ) &&
	    lgFileReadable( "0050000_50_halo__iron.bin_0.1", dum, as ) &&
	    !lgValidAsciiFile( "rauch_h-ni_3d.ascii", as ) )
	{
		fprintf( ioQQQ, " Creating rauch_h-ni_3d.ascii....\n" );
		lgFail = lgFail || RauchInitializeSub( "rauch_h-ni_3d.ascii", "_solar_iron.bin_0.1",
						       telg2, NMODS_HNI, 1, 2, par2, 1 );
		lgFail = lgFail || RauchInitializeSub( "rauch_h-ni_3d.ascii", "_halo__iron.bin_0.1",
						       telg2, NMODS_HNI, 2, 2, par2, 1 );
	}

	/* this is the hydrogen deficient PG1159 grid */
	if( lgFileReadable( "0040000_5.00_33_50_02_15.bin_0.1", dum, as ) &&
	    !lgValidAsciiFile( "rauch_pg1159.ascii", as ) )
	{
		fprintf( ioQQQ, " Creating rauch_pg1159.ascii....\n" );
		lgFail = lgFail || RauchInitializeSub( "rauch_pg1159.ascii", "_33_50_02_15.bin_0.1",
						       telg3, NMODS_PG1159, 1, 1, par2, 2 );
	}

	/* this is the pure hydrogen grid */
	if( lgFileReadable( "0020000_4.00_H_00005-02000A.bin_0.1", dum, as ) &&
	    !lgValidAsciiFile( "rauch_hydr.ascii", as ) )
	{
		fprintf( ioQQQ, " Creating rauch_hydr.ascii....\n" );
		lgFail = lgFail || RauchInitializeSub( "rauch_hydr.ascii", "_H_00005-02000A.bin_0.1",
						       telg4, NMODS_HYDR, 1, 1, par2, 2 );
	}

	/* this is the pure helium grid */
	if( lgFileReadable( "0050000_5.00_He_00005-02000A.bin_0.1", dum, as ) &&
	    !lgValidAsciiFile( "rauch_helium.ascii", as ) )
	{
		fprintf( ioQQQ, " Creating rauch_helium.ascii....\n" );
		lgFail = lgFail || RauchInitializeSub( "rauch_helium.ascii", "_He_00005-02000A.bin_0.1",
						       telg5, NMODS_HELIUM, 1, 1, par2, 2 );
	}

	/* this is the 3D grid for arbitrary H+He mixtures */
	if( lgFileReadable( "0050000_5.00_H+He_1.000_0.000_00005-02000A.bin_0.1", dum, as ) &&
	    !lgValidAsciiFile( "rauch_h+he_3d.ascii", as ) )
	{
		fprintf( ioQQQ, " Creating rauch_h+he_3d.ascii....\n" );
		lgFail = lgFail || RauchInitializeSub( "rauch_h+he_3d.ascii", "_H+He_1.000_0.000_00005-02000A.bin_0.1",
						       telg6, NMODS_HpHE,  1, 11, par3, 2 );
		lgFail = lgFail || RauchInitializeSub( "rauch_h+he_3d.ascii", "_H+He_0.900_0.100_00005-02000A.bin_0.1",
						       telg6, NMODS_HpHE,  2, 11, par3, 2 );
		lgFail = lgFail || RauchInitializeSub( "rauch_h+he_3d.ascii", "_H+He_0.800_0.200_00005-02000A.bin_0.1",
						       telg6, NMODS_HpHE,  3, 11, par3, 2 );
		lgFail = lgFail || RauchInitializeSub( "rauch_h+he_3d.ascii", "_H+He_0.700_0.300_00005-02000A.bin_0.1",
						       telg6, NMODS_HpHE,  4, 11, par3, 2 );
		lgFail = lgFail || RauchInitializeSub( "rauch_h+he_3d.ascii", "_H+He_0.600_0.400_00005-02000A.bin_0.1",
						       telg6, NMODS_HpHE,  5, 11, par3, 2 );
		lgFail = lgFail || RauchInitializeSub( "rauch_h+he_3d.ascii", "_H+He_0.500_0.500_00005-02000A.bin_0.1",
						       telg6, NMODS_HpHE,  6, 11, par3, 2 );
		lgFail = lgFail || RauchInitializeSub( "rauch_h+he_3d.ascii", "_H+He_0.400_0.600_00005-02000A.bin_0.1",
						       telg6, NMODS_HpHE,  7, 11, par3, 2 );
		lgFail = lgFail || RauchInitializeSub( "rauch_h+he_3d.ascii", "_H+He_0.300_0.700_00005-02000A.bin_0.1",
						       telg6, NMODS_HpHE,  8, 11, par3, 2 );
		lgFail = lgFail || RauchInitializeSub( "rauch_h+he_3d.ascii", "_H+He_0.200_0.800_00005-02000A.bin_0.1",
						       telg6, NMODS_HpHE,  9, 11, par3, 2 );
		lgFail = lgFail || RauchInitializeSub( "rauch_h+he_3d.ascii", "_H+He_0.100_0.900_00005-02000A.bin_0.1",
						       telg6, NMODS_HpHE, 10, 11, par3, 2 );
		lgFail = lgFail || RauchInitializeSub( "rauch_h+he_3d.ascii", "_H+He_0.000_1.000_00005-02000A.bin_0.1",
						       telg6, NMODS_HpHE, 11, 11, par3, 2 );
	}

	/* define the major absorption edges that require special attention during rebinning
	 *
	 * NB the frequencies should be chosen here such that they are somewhere in between
	 * the two frequency points that straddle the edge in the atmosphere model, the
	 * software in RebinAtmosphere will seek out the exact values of those two points
	 * e.g.: in the CoStar models the H I edge is straddled by wavelength points at
	 * 911.67 and 911.85 A, so Edges[0] should be chosen somewhere in between (e.g. at 911.76A).
	 *
	 * NB beware not to choose edges too close to one another (i.e. on the order of the
	 * resolution of the Cloudy frequency grid). E.g. the He II Balmer edge nearly coincides
	 * with the H I Ly edge, they should be treated as one edge. Trying to separate them will
	 * almost certainly lead to erroneous behaviour in RebinAtmosphere */
	Edges[0] = 0.99946789f;
	Edges[1] = 1.8071406f;
	Edges[2] = 3.9996377f;

	if( lgFileReadable( "rauch_h-ca_solar.ascii", pc, as ) && !lgValidBinFile( "rauch_h-ca_solar.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "rauch_h-ca_solar.ascii", "rauch_h-ca_solar.mod",Edges,3L, pc );
	if( lgFileReadable( "rauch_h-ca_halo.ascii", pc, as ) && !lgValidBinFile( "rauch_h-ca_halo.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "rauch_h-ca_halo.ascii", "rauch_h-ca_halo.mod", Edges, 3L, pc );
	if( lgFileReadable( "rauch_h-ca_3d.ascii", pc, as ) && !lgValidBinFile( "rauch_h-ca_3d.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "rauch_h-ca_3d.ascii", "rauch_h-ca_3d.mod", Edges, 3L, pc );

	if( lgFileReadable( "rauch_h-ni_solar.ascii", pc, as ) && !lgValidBinFile( "rauch_h-ni_solar.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "rauch_h-ni_solar.ascii", "rauch_h-ni_solar.mod",Edges,3L, pc );
	if( lgFileReadable( "rauch_h-ni_halo.ascii", pc, as ) && !lgValidBinFile( "rauch_h-ni_halo.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "rauch_h-ni_halo.ascii", "rauch_h-ni_halo.mod", Edges, 3L, pc );
	if( lgFileReadable( "rauch_h-ni_3d.ascii", pc, as ) && !lgValidBinFile( "rauch_h-ni_3d.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "rauch_h-ni_3d.ascii", "rauch_h-ni_3d.mod", Edges, 3L, pc );

	if( lgFileReadable( "rauch_pg1159.ascii", pc, as ) && !lgValidBinFile( "rauch_pg1159.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "rauch_pg1159.ascii", "rauch_pg1159.mod", Edges, 3L, pc );

	if( lgFileReadable( "rauch_hydr.ascii", pc, as ) && !lgValidBinFile( "rauch_hydr.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "rauch_hydr.ascii", "rauch_hydr.mod", Edges, 3L, pc );

	if( lgFileReadable( "rauch_helium.ascii", pc, as ) && !lgValidBinFile( "rauch_helium.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "rauch_helium.ascii", "rauch_helium.mod", Edges, 3L, pc );

	if( lgFileReadable( "rauch_h+he_3d.ascii", pc, as ) && !lgValidBinFile( "rauch_h+he_3d.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "rauch_h+he_3d.ascii", "rauch_h+he_3d.mod", Edges, 3L, pc );
	return lgFail;
}

/* RauchInterpolateHCa get one of the Rauch H-Ca model atmospheres, originally by K. Volk */
long RauchInterpolateHCa(double val[], /* val[nval] */
			 long *nval,
			 long *ndim,
			 bool lgHalo,
			 bool lgList,
			 double *Tlow,
			 double *Thigh)
{
	stellar_grid grid;

	DEBUG_ENTRY( "RauchInterpolateHCa()" );

	if( *ndim == 3 )
		grid.name = "rauch_h-ca_3d.mod";
	else
		grid.name = ( lgHalo ? "rauch_h-ca_halo.mod" : "rauch_h-ca_solar.mod" );
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	grid.ident = "  H-Ca Rauch";
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	FreeGrid( &grid );
	return rfield.nupper;
}

/* RauchInterpolateHNi get one of the Rauch H-Ni model atmospheres */
long RauchInterpolateHNi(double val[], /* val[nval] */
			 long *nval,
			 long *ndim,
			 bool lgHalo,
			 bool lgList,
			 double *Tlow,
			 double *Thigh)
{
	stellar_grid grid;

	DEBUG_ENTRY( "RauchInterpolateHNi()" );

	if( *ndim == 3 )
		grid.name = "rauch_h-ni_3d.mod";
	else
		grid.name = ( lgHalo ? "rauch_h-ni_halo.mod" : "rauch_h-ni_solar.mod" );
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	grid.ident = "  H-Ni Rauch";
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	FreeGrid( &grid );
	return rfield.nupper;
}

/* RauchInterpolatePG1159 get one of the Rauch PG1159 model atmospheres */
long RauchInterpolatePG1159(double val[], /* val[nval] */
			    long *nval,
			    long *ndim,
			    bool lgList,
			    double *Tlow,
			    double *Thigh)
{
	stellar_grid grid;

	DEBUG_ENTRY( "RauchInterpolatePG1159()" );

	grid.name = "rauch_pg1159.mod";
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	grid.ident = "PG1159 Rauch";
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	FreeGrid( &grid );
	return rfield.nupper;
}

/* RauchInterpolateHydr get one of the Rauch pure hydrogen model atmospheres */
long RauchInterpolateHydr(double val[], /* val[nval] */
			  long *nval,
			  long *ndim,
			  bool lgList,
			  double *Tlow,
			  double *Thigh)
{
	stellar_grid grid;

	DEBUG_ENTRY( "RauchInterpolateHydr()" );

	grid.name = "rauch_hydr.mod";
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	grid.ident = "  Hydr Rauch";
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	FreeGrid( &grid );
	return rfield.nupper;
}

/* RauchInterpolateHelium get one of the Rauch pure helium model atmospheres */
long RauchInterpolateHelium(double val[], /* val[nval] */
			    long *nval,
			    long *ndim,
			    bool lgList,
			    double *Tlow,
			    double *Thigh)
{
	stellar_grid grid;

	DEBUG_ENTRY( "RauchInterpolateHelium()" );

	grid.name = "rauch_helium.mod";
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	grid.ident = "Helium Rauch";
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	FreeGrid( &grid );
	return rfield.nupper;
}

/* RauchInterpolateHpHe get one of the Rauch hydrogen plus helium model atmospheres */
long RauchInterpolateHpHe(double val[], /* val[nval] */
			  long *nval,
			  long *ndim,
			  bool lgList,
			  double *Tlow,
			  double *Thigh)
{
	stellar_grid grid;

	DEBUG_ENTRY( "RauchInterpolateHpHe()" );

	grid.name = "rauch_h+he_3d.mod";
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	grid.ident = "  H+He Rauch";
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	FreeGrid( &grid );
	return rfield.nupper;
}

/* StarburstInitialize does the actual work of preparing the ascii file */
bool StarburstInitialize(const char chInName[],
			 const char chOutName[])
{
	char chLine[INPUT_LINE_LENGTH];          /* used for getting input lines */

	bool lgHeader = true;
	long int i, j, nmods, ngp;

	size_t nsb_sz = (size_t)NSB99;

	double *wavl, *fluxes[MNTS], Age[MNTS], lwavl;

	FILE *ioOut,  /* pointer to output file we came here to create*/
	  *ioIn;      /* pointer to input files we will read */

	DEBUG_ENTRY( "StarburstInitialize()" );

	for( i=0; i < MNTS; i++ )
		fluxes[i] = NULL;

	/* grab some space for the wavelengths and fluxes */
	wavl = (double *)MALLOC( sizeof(double)*nsb_sz);

	ioIn = open_data( chInName, "r", AS_LOCAL_ONLY );

	lwavl = 0.;
	nmods = 0;
	ngp = 0;

	while( read_whole_line( chLine, INPUT_LINE_LENGTH, ioIn ) != NULL )
	{
		if( !lgHeader )
		{
			double cage, cwavl, cfl1;

			/* format: age/yr wavl/Angstrom log10(flux_total) log10(flux_stellar) log10(flux_neb) */
			/* we are only interested in the total flux, so we ignore the remaining numbers */
			if( sscanf( chLine, " %le %le %le", &cage, &cwavl, &cfl1 ) != 3 )
			{
				fprintf( ioQQQ, "syntax error in data of Starburst grid.\n" );
				goto error;
			}

			if( cwavl < lwavl )
			{
				++nmods;
				ngp = 0;

				if( nmods >= MNTS )
				{
					fprintf( ioQQQ, "too many time steps in Starburst grid.\n" );
					fprintf( ioQQQ, "please increase MNTS and recompile.\n" );
					goto error;
				}
			}

			if( ngp == 0 )
			{
				fluxes[nmods] = (double *)MALLOC( sizeof(double)*nsb_sz);
				Age[nmods] = cage;
			}

			if( ngp >= (long)nsb_sz )
			{
				/* this should only be needed when nmods == 0 */
				ASSERT( nmods == 0 );

				nsb_sz *= 2;
				fluxes[0] = (double *)REALLOC(fluxes[0],sizeof(double)*nsb_sz);
				wavl = (double *)REALLOC(wavl,sizeof(double)*nsb_sz);
			}

			if( !fp_equal(Age[nmods],cage,10) )
			{
				fprintf( ioQQQ, "age error in Starburst grid.\n" );
				goto error;
			}

			if( nmods == 0 )
				wavl[ngp] = cwavl;
			else
			{
				if( !fp_equal(wavl[ngp],cwavl,10) )
				{
					fprintf( ioQQQ, "wavelength error in Starburst grid.\n" );
					goto error;
				}
			}

			/* arbitrarily renormalize to flux in erg/cm^2/s/A at 1kpc */
			/* constant is log10( 4*pi*(kpc/cm)^2 ) */
			fluxes[nmods][ngp] = pow( 10., cfl1 - 44.077911 );

			lwavl = cwavl;
			++ngp;
		}

		if( lgHeader && strncmp( &chLine[1], "TIME [YR]", 9 ) == 0 )
			lgHeader = false;
	}

	if( lgHeader )
	{
		/* this happens when the "TIME [YR]" string was not found in column 1 of the file */
		fprintf( ioQQQ, "syntax error in header of Starburst grid.\n" );
		goto error;
	}

	++nmods;

	/* finished - close the unit */
	fclose(ioIn);

	ioOut = open_data( chOutName, "w", AS_LOCAL_ONLY );

	fprintf( ioOut, "  %ld\n", VERSION_ASCII );
	fprintf( ioOut, "  %d\n", 1 );
	fprintf( ioOut, "  %d\n", 1 );
	fprintf( ioOut, "  Age\n" );
	fprintf( ioOut, "  %ld\n", nmods );
	fprintf( ioOut, "  %ld\n", ngp );
	/* Starburst99 models give the wavelength in Angstrom */
	fprintf( ioOut, "  lambda\n" );
	/* conversion factor to Angstrom */
	fprintf( ioOut, "  %.8e\n", 1. );
	/* Starburst99 models give the total flux F_lambda in erg/s/A, will be renormalized at 1 kpc */
	fprintf( ioOut, "  F_lambda\n" );
	/* this factor is irrelevant since Teff check will not be carried out */
	fprintf( ioOut, "  %.8e\n", 1. );
	/* write out the Ages */
	for( i=0; i < nmods; i++ )
	{
		fprintf( ioOut, "  %.3e", Age[i] );
		if( ((i+1)%4) == 0 )
			fprintf( ioOut, "\n" );
	}
	if( (i%4) != 0 )
		fprintf( ioOut, "\n" );

	fprintf( ioQQQ, " Writing: " );

	/* write out the wavelength grid */
	for( j=0; j < ngp; j++ )
	{
		fprintf( ioOut, "  %.4e", wavl[j] );
		/* want to have 5 numbers per line */
		if( ((j+1)%5) == 0 )
			fprintf( ioOut, "\n" );
	}
	/* need to throw a newline if we did not end on an exact line */
	if( (j%5) != 0 )
		fprintf( ioOut, "\n" );

	/* print to screen to show progress */
	fprintf( ioQQQ, "." );
	fflush( ioQQQ );

	for( i=0; i < nmods; i++ )
	{
		for( j=0; j < ngp; j++ )
		{
			fprintf( ioOut, "  %.4e", fluxes[i][j] );
			/* want to have 5 numbers per line */
			if( ((j+1)%5) == 0 )
				fprintf( ioOut, "\n" );
		}
		/* need to throw a newline if we did not end on an exact line */
		if( (j%5) != 0 )
			fprintf( ioOut, "\n" );

		/* print to screen to show progress */
		fprintf( ioQQQ, "." );
		fflush( ioQQQ );
	}

	fprintf( ioQQQ, " done.\n" );

	fclose(ioOut);

	/* free the space grabbed above */
	for( i=0; i < MNTS; i++ )
		FREE_SAFE( fluxes[i] );
	FREE_CHECK( wavl );
	return false;

error:
	for( i=0; i < MNTS; i++ )
		FREE_SAFE( fluxes[i] );
	FREE_CHECK( wavl );
	return true;
}

/* StarburstCompile, rebin Starburst99 model output to match energy grid of code */
bool StarburstCompile(process_counter& pc)
{
	realnum Edges[1];

	bool lgFail = false;

	DEBUG_ENTRY( "StarburstCompile()" );

	fprintf( ioQQQ, " StarburstCompile on the job.\n" );

	process_counter dum;
	access_scheme as = AS_LOCAL_ONLY_TRY;

	if( lgFileReadable( "starburst99.stb99", dum, as ) && !lgValidAsciiFile( "starburst99.ascii", as ) )
		lgFail = lgFail || StarburstInitialize( "starburst99.stb99", "starburst99.ascii" );
	if( lgFileReadable( "starburst99.ascii", pc, as ) && !lgValidBinFile( "starburst99.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "starburst99.ascii", "starburst99.mod", Edges, 0L, pc );

	if( lgFileReadable( "starburst99_2d.ascii", pc, as ) && !lgValidBinFile( "starburst99_2d.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "starburst99_2d.ascii", "starburst99_2d.mod", Edges, 0L, pc );
	return lgFail;
}

/* TlustyCompile rebin Tlusty BSTAR2006/OSTAR2002 stellar models to match energy grid of code */
int TlustyCompile(process_counter& pc)
{
	/* these contain frequencies for the major absorption edges */
	realnum Edges[1];

	bool lgFail = false;

	DEBUG_ENTRY( "TlustyCompile()" );

	fprintf( ioQQQ, " TlustyCompile on the job.\n" );

	access_scheme as = AS_LOCAL_ONLY_TRY;

	if( lgFileReadable( "obstar_merged_p03.ascii", pc, as ) && !lgValidBinFile( "obstar_merged_p03.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "obstar_merged_p03.ascii", "obstar_merged_p03.mod", Edges, 0L, pc );
	if( lgFileReadable( "obstar_merged_p00.ascii", pc, as ) && !lgValidBinFile( "obstar_merged_p00.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "obstar_merged_p00.ascii", "obstar_merged_p00.mod", Edges, 0L, pc );
	if( lgFileReadable( "obstar_merged_m03.ascii", pc, as ) && !lgValidBinFile( "obstar_merged_m03.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "obstar_merged_m03.ascii", "obstar_merged_m03.mod", Edges, 0L, pc );
	if( lgFileReadable( "obstar_merged_m07.ascii", pc, as ) && !lgValidBinFile( "obstar_merged_m07.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "obstar_merged_m07.ascii", "obstar_merged_m07.mod", Edges, 0L, pc );
	if( lgFileReadable( "obstar_merged_m10.ascii", pc, as ) && !lgValidBinFile( "obstar_merged_m10.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "obstar_merged_m10.ascii", "obstar_merged_m10.mod", Edges, 0L, pc );
	if( lgFileReadable( "obstar_merged_m99.ascii", pc, as ) && !lgValidBinFile( "obstar_merged_m99.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "obstar_merged_m99.ascii", "obstar_merged_m99.mod", Edges, 0L, pc );

	if( lgFileReadable( "obstar_merged_3d.ascii", pc, as ) && !lgValidBinFile( "obstar_merged_3d.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "obstar_merged_3d.ascii", "obstar_merged_3d.mod", Edges, 0L, pc );

	if( lgFileReadable( "bstar2006_p03.ascii", pc, as ) && !lgValidBinFile( "bstar2006_p03.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "bstar2006_p03.ascii", "bstar2006_p03.mod", Edges, 0L, pc );
	if( lgFileReadable( "bstar2006_p00.ascii", pc, as ) && !lgValidBinFile( "bstar2006_p00.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "bstar2006_p00.ascii", "bstar2006_p00.mod", Edges, 0L, pc );
	if( lgFileReadable( "bstar2006_m03.ascii", pc, as ) && !lgValidBinFile( "bstar2006_m03.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "bstar2006_m03.ascii", "bstar2006_m03.mod", Edges, 0L, pc );
	if( lgFileReadable( "bstar2006_m07.ascii", pc, as ) && !lgValidBinFile( "bstar2006_m07.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "bstar2006_m07.ascii", "bstar2006_m07.mod", Edges, 0L, pc );
	if( lgFileReadable( "bstar2006_m10.ascii", pc, as ) && !lgValidBinFile( "bstar2006_m10.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "bstar2006_m10.ascii", "bstar2006_m10.mod", Edges, 0L, pc );
	if( lgFileReadable( "bstar2006_m99.ascii", pc, as ) && !lgValidBinFile( "bstar2006_m99.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "bstar2006_m99.ascii", "bstar2006_m99.mod", Edges, 0L, pc );

	if( lgFileReadable( "bstar2006_3d.ascii", pc, as ) && !lgValidBinFile( "bstar2006_3d.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "bstar2006_3d.ascii", "bstar2006_3d.mod", Edges, 0L, pc );

	if( lgFileReadable( "ostar2002_p03.ascii", pc, as ) && !lgValidBinFile( "ostar2002_p03.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "ostar2002_p03.ascii", "ostar2002_p03.mod", Edges, 0L, pc );
	if( lgFileReadable( "ostar2002_p00.ascii", pc, as ) && !lgValidBinFile( "ostar2002_p00.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "ostar2002_p00.ascii", "ostar2002_p00.mod", Edges, 0L, pc );
	if( lgFileReadable( "ostar2002_m03.ascii", pc, as ) && !lgValidBinFile( "ostar2002_m03.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "ostar2002_m03.ascii", "ostar2002_m03.mod", Edges, 0L, pc );
	if( lgFileReadable( "ostar2002_m07.ascii", pc, as ) && !lgValidBinFile( "ostar2002_m07.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "ostar2002_m07.ascii", "ostar2002_m07.mod", Edges, 0L, pc );
	if( lgFileReadable( "ostar2002_m10.ascii", pc, as ) && !lgValidBinFile( "ostar2002_m10.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "ostar2002_m10.ascii", "ostar2002_m10.mod", Edges, 0L, pc );
	if( lgFileReadable( "ostar2002_m15.ascii", pc, as ) && !lgValidBinFile( "ostar2002_m15.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "ostar2002_m15.ascii", "ostar2002_m15.mod", Edges, 0L, pc );
	if( lgFileReadable( "ostar2002_m17.ascii", pc, as ) && !lgValidBinFile( "ostar2002_m17.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "ostar2002_m17.ascii", "ostar2002_m17.mod", Edges, 0L, pc );
	if( lgFileReadable( "ostar2002_m20.ascii", pc, as ) && !lgValidBinFile( "ostar2002_m20.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "ostar2002_m20.ascii", "ostar2002_m20.mod", Edges, 0L, pc );
	if( lgFileReadable( "ostar2002_m30.ascii", pc, as ) && !lgValidBinFile( "ostar2002_m30.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "ostar2002_m30.ascii", "ostar2002_m30.mod", Edges, 0L, pc );
	if( lgFileReadable( "ostar2002_m99.ascii", pc, as ) && !lgValidBinFile( "ostar2002_m99.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "ostar2002_m99.ascii", "ostar2002_m99.mod", Edges, 0L, pc );

	if( lgFileReadable( "ostar2002_3d.ascii", pc, as ) && !lgValidBinFile( "ostar2002_3d.mod", pc, as ) )
		lgFail = lgFail || lgCompileAtmosphere( "ostar2002_3d.ascii", "ostar2002_3d.mod", Edges, 0L, pc );
	return lgFail;
}

/* TlustyInterpolate get one of the Tlusty OBSTAR_MERGED/BSTAR2006/OSTAR2002 model atmospheres */
long TlustyInterpolate(double val[], /* val[nval] */
		       long *nval,
		       long *ndim,
		       tl_grid tlg,
		       const char *chMetalicity,
		       bool lgList,
		       double *Tlow,
		       double *Thigh)
{
	char chIdent[13];
	stellar_grid grid;

	DEBUG_ENTRY( "TlustyInterpolate()" );

	if( tlg == TL_OBSTAR )
		grid.name = "obstar_merged_";
	else if( tlg == TL_BSTAR )
		grid.name = "bstar2006_";
	else if( tlg == TL_OSTAR )
		grid.name = "ostar2002_";
	else
		TotalInsanity();
	if( *ndim == 3 )
		grid.name += "3d";
	else
		grid.name += chMetalicity;
	grid.name += ".mod";
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	if( *ndim == 3 )
	{
		strcpy( chIdent, "3-dim" );
	}
	else
	{
		strcpy( chIdent, "Z " );
		strcat( chIdent, chMetalicity );
	}
	if( tlg == TL_OBSTAR )
		strcat( chIdent, " OBstar" );
	else if( tlg == TL_BSTAR )
		strcat( chIdent, " Bstr06" );
	else if( tlg == TL_OSTAR )
		strcat( chIdent, " Ostr02" );
	else
		TotalInsanity();
	grid.ident = chIdent;
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	FreeGrid( &grid );
	return rfield.nupper;
}

/* WernerCompile rebin Werner stellar models to match energy grid of code */
/* >>chng 05 nov 16, added return value to indicate success (0) or failure (1) */
int WernerCompile(process_counter& pc)
{
	/* these contain frequencies for the major absorption edges */
	realnum Edges[3];

	bool lgFail = false;

	DEBUG_ENTRY( "WernerCompile()" );

	fprintf( ioQQQ, " WernerCompile on the job.\n" );

	/* define the major absorption edges that require special attention during rebinning
	 *
	 * NB the frequencies should be chosen here such that they are somewhere in between
	 * the two frequency points that straddle the edge in the atmosphere model, the
	 * software in RebinAtmosphere will seek out the exact values of those two points
	 * e.g.: in the CoStar models the H I edge is straddled by wavelength points at
	 * 911.67 and 911.85 A, so Edges[0] should be chosen somewhere in between (e.g. at 911.76A).
	 *
	 * NB beware not to choose edges too close to one another (i.e. on the order of the
	 * resolution of the Cloudy frequency grid). E.g. the He II Balmer edge nearly coincides
	 * with the H I Ly edge, they should be treated as one edge. Trying to separate them will
	 * almost certainly lead to erroneous behaviour in RebinAtmosphere */
	Edges[0] = 0.99946789f;
	Edges[1] = 1.8071406f;
	Edges[2] = 3.9996377f;

	/* The "kwerner.ascii" file is a modified ascii dump of the Klaus Werner 
	 * stellar model files which he gave to me in 1992.  The first set of values 
	 * is the frequency grid (in Ryd) followed by the atmosphere models in order
	 * of increasing temperature and log(g). The following comments are already
	 * incorporated in the modified kwerner.ascii file that is supplied with Cloudy.
	 *
	 * >>chng 00 oct 18, The frequency grid was slightly tweaked compared to the
	 * original values supplied by Klaus Werner to make it monotonically increasing;
	 * this is due to there being fluxes above and below certain wavelengths where
	 * the opacity changes (i.e. the Lyman and Balmer limits for example) which are 
	 * assigned the same wavelength in the original Klaus Werner files. PvH
	 *
	 * >>chng 00 oct 20, StarEner[172] is out of sequence. As per the Klaus Werner comment,
	 * it should be omitted. The energy grid is very dense in this region and was most likely
	 * intended to sample an absorption line which was not included in this particular grid.
	 * StarFlux[172] is therefore always equal to the flux in neighbouring points (at least
	 * those with slightly smaller energies). It is therefore safe to ignore this point. PvH
	 *
	 * >>chng 00 oct 20, As per the same comment, StarFlux[172] is also deleted. PvH */

	access_scheme as = AS_LOCAL_ONLY_TRY;

	if( lgFileReadable( "kwerner.ascii", pc, as ) && !lgValidBinFile( "kwerner.mod", pc, as ) )
		lgFail = lgCompileAtmosphere( "kwerner.ascii", "kwerner.mod", Edges, 3L, pc );
	return lgFail;
}

/* WernerInterpolate read in and interpolate on Werner grid of PN atmospheres, originally by K Volk */
long WernerInterpolate(double val[], /* val[nval] */
		       long *nval,
		       long *ndim,
		       bool lgList,
		       double *Tlow,
		       double *Thigh)
{
	stellar_grid grid;

	DEBUG_ENTRY( "WernerInterpolate()" );

	/* This subroutine was added (28 dec 1992) to read from the set of
	 * hot white dwarf model atmospheres from Klaus Werner at Kiel. The 
	 * values are read in (energy in Rydberg units, f_nu in cgs units)
	 * for any of the 20 models. Each model had 513 points before rebinning.
	 * The Rayleigh-Jeans tail was extrapolated. */

	grid.name = "kwerner.mod";
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	grid.ident = "Klaus Werner";
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	/* Note on the interpolation: 26 October 2000 (Peter van Hoof)
	 *
	 * I computed the effective temperature for a random sample of interpolated
	 * atmospheres by integrating the flux as shown above and compared the results
	 * with the expected effective temperature using DELTA = (COMP-EXPEC)/EXPEC.
	 *
	 * I found that the average discrepancy was:
	 *
	 *     DELTA = -0.71% +/- 0.71% (sample size 5000)
	 *
	 * The most extreme discrepancies were
	 *     -4.37% <= DELTA <= 0.24%
	 *
	 * The most negative discrepancies were for Teff =  95 kK, log(g) = 5
	 * The most positive discrepancies were for Teff = 160 kK, log(g) = 8
	 *
	 * Since Cloudy checks the scaling elsewhere there is no need to re-scale 
	 * things here, but this inaccuracy should be kept in mind since it could
	 * indicate problems with the flux distribution */

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	FreeGrid( &grid );
	return rfield.nupper;
}

/* WMBASICCompile rebin WMBASIC stellar models to match energy grid of code */
int WMBASICCompile(process_counter& pc)
{
	/* these contain frequencies for the major absorption edges */
	realnum Edges[3];

	bool lgFail = false;

	DEBUG_ENTRY( "WMBASICCompile()" );

	fprintf( ioQQQ, " WMBASICCompile on the job.\n" );

	/* define the major absorption edges that require special attention during rebinning */
	Edges[0] = 0.99946789f;
	Edges[1] = 1.8071406f;
	Edges[2] = 3.9996377f;

	access_scheme as = AS_LOCAL_ONLY_TRY;

	if( lgFileReadable( "wmbasic.ascii", pc, as ) && !lgValidBinFile( "wmbasic.mod", pc, as ) )
		lgFail = lgCompileAtmosphere( "wmbasic.ascii", "wmbasic.mod", Edges, 3L, pc );
	return lgFail;
}

/* WMBASICInterpolate read in and interpolate on WMBASIC grid of hot star atmospheres */
long WMBASICInterpolate(double val[], /* val[nval] */
			long *nval,
			long *ndim,
			bool lgList,
			double *Tlow,
			double *Thigh)
{
	stellar_grid grid;

	DEBUG_ENTRY( "WMBASICInterpolate()" );

	grid.name = "wmbasic.mod";
	grid.scheme = AS_DATA_OPTIONAL;
	/* identification of this atmosphere set, used in
	 * the Cloudy output, *must* be 12 characters long */
	grid.ident = "     WMBASIC";
	/* the Cloudy command needed to recompile the binary model file */
	grid.command = "COMPILE STARS";

	InitGrid( &grid, lgList );

	CheckVal( &grid, val, nval, ndim );

	InterpolateRectGrid( &grid, val, Tlow, Thigh );

	FreeGrid( &grid );
	return rfield.nupper;
}

/* CompileAtmosphereCoStar rebin costar stellar atmospheres to match cloudy energy grid, 
 * called by the compile stars command */
STATIC bool lgCompileAtmosphereCoStar(const char chFNameIn[],
				      const char chFNameOut[],
				      const realnum Edges[], /* Edges[nedges] */
				      long nedges,
				      process_counter& pc)
{
	char chLine[INPUT_LINE_LENGTH];
	char names[MDIM][MNAM+1];
	int32 val[7];
	uint32 uval[2];
	double dval[3];
	char md5sum[NMD5];
	long int i, j, nskip, nModels, nWL;

	/* these will be malloced into large work arrays*/
	realnum *StarEner = NULL, *StarFlux = NULL, *CloudyFlux = NULL;
	/* this will hold all the model parameters */
	mpp *telg = NULL;

	FILE *ioIN;  /* used for input */
	FILE *ioOUT; /* used for output */

	DEBUG_ENTRY( "CompileAtmosphereCoStar()" );

	/* This is a program to re-bin the costar stellar model spectra to match the 
	 * Cloudy grid.  For short wavelengths I will use a power law extrapolation
	 * of the model values (which should be falling rapidly) if needed.  At long
	 * wavelengths I will assume Rayleigh-Jeans from the last stellar model point
	 * to extrapolate to 1 cm wavelength. */

	/* This version uses power-law interpolation between the points of the stellar model. */

	/* read the original data file obtained off the web, 
	 * open as read only */
	try
	{
		ioIN = open_data( chFNameIn, "r", AS_LOCAL_ONLY );
	}
	catch( cloudy_exit )
	{
		goto error;
	}
	fprintf( ioQQQ, " CompileAtmosphereCoStar got %s.\n", chFNameIn );

	/* get first line and see how many more to skip */
	if( read_whole_line( chLine, (int)sizeof(chLine), ioIN ) == NULL )
	{
		fprintf( ioQQQ, " CompileAtmosphereCoStar fails reading nskip.\n" );
		goto error;
	}
	sscanf( chLine, "%li", &nskip );

	/* now skip the header information */
	for( i=0; i < nskip; ++i )
	{
		if( read_whole_line( chLine, (int)sizeof(chLine), ioIN ) == NULL )
		{
			fprintf( ioQQQ, " CompileAtmosphereCoStar fails skipping header.\n" );
			goto error;
		}
	}

	/* now get number of models and number of wavelengths */
	if( read_whole_line( chLine, (int)sizeof(chLine), ioIN ) == NULL )
	{
		fprintf( ioQQQ, " CompileAtmosphereCoStar fails reading nModels, nWL.\n" );
		goto error;
	}
	sscanf( chLine, "%li%li", &nModels, &nWL );

	if( nModels <= 0 || nWL <= 0 )
	{
		fprintf( ioQQQ, " CompileAtmosphereCoStar scanned off impossible values for nModels=%li or nWL=%li\n",
			 nModels, nWL );
		goto error;
	}

	/* allocate space for the stellar parameters */
	telg = (mpp *)CALLOC( (size_t)nModels, sizeof(mpp) );

	/* get all model parameters for the atmospheres */
	for( i=0; i < nModels; ++i )
	{
		if( read_whole_line( chLine, (int)sizeof(chLine), ioIN ) == NULL )
		{
			fprintf( ioQQQ, " CompileAtmosphereCoStar fails reading model parameters.\n" );
			goto error;
		}
		/* first letter on line is indicator of grid */
		telg[i].chGrid = chLine[0];
		/* get the model id number */
		sscanf( chLine+1, "%i", &telg[i].modid );
		/* get the temperature */
		sscanf( chLine+23, "%lg", &telg[i].par[0] );
		/* get the surface gravity */
		sscanf( chLine+31, "%lg", &telg[i].par[1] );
		/* get the ZAMS mass */
		sscanf( chLine+7, "%lg", &telg[i].par[2] );
		/* get the model age */
		sscanf( chLine+15, "%lg", &telg[i].par[3] );

		/* the code in parse_table.cpp implicitly depends on this! */
		ASSERT( telg[i].par[2] > 10. );
		ASSERT( telg[i].par[3] > 10. );

		/* convert ZAMS masses to logarithms */
		telg[i].par[2] = log10(telg[i].par[2]);
	}

	/* this will be the file we create, that will be read to compute models, 
	 * open to write binary */
	try
	{
		ioOUT = open_data( chFNameOut, "wb", AS_LOCAL_ONLY );
	}
	catch( cloudy_exit )
	{
		goto error;
	}

	val[0] = (int32)VERSION_BIN;
	val[1] = (int32)MDIM;
	val[2] = (int32)MNAM;
	val[3] = (int32)2;
	val[4] = (int32)4;
	val[5] = (int32)nModels;
	val[6] = (int32)rfield.nupper;
	uval[0] = sizeof(val) + sizeof(uval) + sizeof(dval) + sizeof(md5sum) +
		sizeof(names) + nModels*sizeof(mpp); /* nOffset */
	uval[1] = rfield.nupper*sizeof(realnum); /* nBlocksize */
	dval[0] = double(rfield.emm);
	dval[1] = double(rfield.egamry);
	dval[2] = double(continuum.ResolutionScaleFactor);

	strncpy( md5sum, continuum.mesh_md5sum.c_str(), sizeof(md5sum) );

	strncpy( names[0], "Teff\0\0",  MNAM+1 );
	strncpy( names[1], "log(g)",    MNAM+1 );
	strncpy( names[2], "log(M)",    MNAM+1 );
	strncpy( names[3], "Age\0\0\0", MNAM+1 );

	if( fwrite( val, sizeof(val), 1, ioOUT ) != 1 ||
	    fwrite( uval, sizeof(uval), 1, ioOUT ) != 1 ||
	    /* write out the lower, upper bound of the energy mesh, and the res scale factor */
	    fwrite( dval, sizeof(dval), 1, ioOUT ) != 1 ||
	    /* write out the (modified) md5 checksum of continuum_mesh.ini */
	    fwrite( md5sum, sizeof(md5sum), 1, ioOUT ) != 1 ||
	    fwrite( names, sizeof(names), 1, ioOUT ) != 1 ||
	    /* write out the array of {Teff,log(g)} pairs */
	    fwrite( telg, sizeof(mpp), (size_t)nModels, ioOUT ) != (size_t)nModels ||
	    /* write out the cloudy energy grid for later sanity checks */
	    fwrite( rfield.AnuOrg, (size_t)uval[1], 1, ioOUT ) != 1 )
	{
		fprintf( ioQQQ, " CompileAtmosphereCoStar failed writing header of output file.\n" );
		goto error;
	}

	/* MALLOC some workspace */
	StarEner = (realnum *)MALLOC( sizeof(realnum)*nWL );
	StarFlux = (realnum *)MALLOC( sizeof(realnum)*nWL );
	CloudyFlux = (realnum *)MALLOC( (size_t)uval[1] );

	fprintf( ioQQQ, " Compiling: " );

	/* get the star data */
	for( i=0; i < nModels; ++i )
	{
		/* get number to skip */
		if( read_whole_line( chLine, (int)sizeof(chLine), ioIN ) == NULL )
		{
			fprintf( ioQQQ, " CompileAtmosphereCoStar fails reading the skip to next spectrum.\n" );
			goto error;
		}
		sscanf( chLine, "%li", &nskip );

		for( j=0; j < nskip; ++j )
		{
			if( read_whole_line( chLine, (int)sizeof(chLine), ioIN ) == NULL )
			{
				fprintf( ioQQQ, " CompileAtmosphereCoStar fails doing the skip.\n" );
				goto error;
			}
		}

		/* now read in the wavelength and flux for this star, read in 
		 * backwards since we want to be in increasing energy order rather
		 * than wavelength */
		for( j=nWL-1; j >= 0; --j )
		{
			if( read_whole_line( chLine, (int)sizeof(chLine), ioIN ) == NULL )
			{
				fprintf( ioQQQ, " CompileAtmosphereCoStar fails reading the spectral data.\n" );
				goto error;
			}
			double help1, help2;
			sscanf( chLine, "%lg %lg", &help1, &help2 );

			/* continuum flux was log, convert to linear, also do
			 * conversion from "astrophysical" flux to F_nu in cgs units */
			StarFlux[j] = (realnum)(PI*pow(10.,help2));
			/* StarEner was in Angstroms, convert to Ryd */
			StarEner[j] = (realnum)(RYDLAM/help1);

			/* sanity check */
			if( j < nWL-1 )
				ASSERT( StarEner[j] < StarEner[j+1] );
		}

		/* this will do the heavy lifting, and define arrays used below,
		 * NB the lowest energy point in these grids appears to be bogus.
		 * tell rebin about nWL-1 */
		RebinAtmosphere(nWL-1, StarEner+1, StarFlux+1, CloudyFlux, nedges, Edges );

		/* write the continuum out as a binary file */
		if( fwrite( CloudyFlux, (size_t)uval[1], 1, ioOUT ) != 1 )
		{
			fprintf( ioQQQ, " CompileAtmosphereCoStar failed writing star flux.\n" );
			goto error;
		}

		fprintf( ioQQQ, "." );
		fflush( ioQQQ );
	}

	fprintf( ioQQQ, " done.\n" );

	fclose( ioIN );
	fclose( ioOUT );

	FREE_CHECK( telg );
	FREE_CHECK( StarEner );
	FREE_CHECK( StarFlux );
	FREE_CHECK( CloudyFlux );

	fprintf( ioQQQ, "\n CompileAtmosphereCoStar completed ok.\n\n" );
	++pc.nOK;
	return false;

error:
	FREE_SAFE( telg );
	FREE_SAFE( StarEner );
	FREE_SAFE( StarFlux );
	FREE_SAFE( CloudyFlux );
	++pc.nFail;
	return true;
}

/* InterpolateGridCoStar read in and interpolate on costar grid of windy O atmospheres */
STATIC void InterpolateGridCoStar(const stellar_grid *grid, /* struct with all the grid parameters */
				  const double val[], /* val[0]: Teff for imode = 1,2; M_ZAMS for imode = 3;
						       * age for imode = 4 */
				                      /* val[1]: nmodid for imode = 1; log(g) for imode = 2;
						       * age for imode = 3; M_ZAMS for imode = 4 */
				  double *val0_lo,
				  double *val0_hi)
{
	long i, j, k, nmodid, off, ptr;
	long *indloTr, *indhiTr, useTr[2];
	long indlo[2], indhi[2], index[2];
	realnum *ValTr;
	double lval[2], aval[4];

	DEBUG_ENTRY( "InterpolateGridCoStar()" );

	switch( grid->imode )
	{
	case IM_COSTAR_TEFF_MODID:
	case IM_COSTAR_TEFF_LOGG:
		lval[0] = val[0];
		lval[1] = val[1];
		off = 0;
		break;
	case IM_COSTAR_MZAMS_AGE:
		lval[0] = log10(val[0]); /* use log10(M_ZAMS) internally */
		lval[1] = val[1];
		off = 2;
		break;
	case IM_COSTAR_AGE_MZAMS:
		/* swap parameters, hence mimic IM_COSTAR_MZAMS_AGE */
		lval[0] = log10(val[1]); /* use log10(M_ZAMS) internally */
		lval[1] = val[0];
		off = 2;
		break;
	default:
		fprintf( ioQQQ, " InterpolateGridCoStar called with insane value for imode: %d.\n", grid->imode );
		cdEXIT(EXIT_FAILURE);
	}

	nmodid = (long)(lval[1]+0.5);

	ASSERT( rfield.lgContMalloc[rfield.nShape] );

	/* read in the saved cloudy energy scale so we can confirm this is a good image */
	GetBins( grid, rfield.tNu[rfield.nShape] );

#	if DEBUGPRT	
	/* check whether the models in the grid have the correct effective temperature */
	ValidateGrid( grid, 0.005 );
#	endif

	/* now allocate some temp workspace */
	ValTr = (realnum *)MALLOC( sizeof(realnum)*grid->nTracks );
	indloTr = (long *)MALLOC( sizeof(long)*grid->nTracks );
	indhiTr = (long *)MALLOC( sizeof(long)*grid->nTracks );

	/* first do horizontal search, i.e. search along individual tracks */
	for( j=0; j < grid->nTracks; j++ )
	{
		if( grid->imode == IM_COSTAR_TEFF_MODID )
		{
			if( grid->trackLen[j] >= nmodid ) {
				index[0] = nmodid - 1;
				index[1] = j;
				ptr = grid->jval[JIndex(grid,index)];
				indloTr[j] = ptr;
				indhiTr[j] = ptr;
				ValTr[j] = (realnum)grid->telg[ptr].par[off];
			}
			else
			{
				indloTr[j] = -2;
				indhiTr[j] = -2;
				ValTr[j] = -FLT_MAX;
			}
		}
		else
		{
			FindHCoStar( grid, j, lval[1], off, ValTr, indloTr, indhiTr );
		}
	}

#	if DEBUGPRT
	for( j=0; j < grid->nTracks; j++ ) 
	{
		if( indloTr[j] >= 0 ) 
			printf( "track %c: models %c%d, %c%d, val %g\n",
				(char)('A'+j), grid->telg[indloTr[j]].chGrid, grid->telg[indloTr[j]].modid,
				grid->telg[indhiTr[j]].chGrid, grid->telg[indhiTr[j]].modid, ValTr[j]);
	}
#	endif

	/* now do vertical search, i.e. interpolate between tracks */
	FindVCoStar( grid, lval[0], ValTr, useTr );

	/* This should only happen when InterpolateGridCoStar is called in non-optimizing mode,
	 * when optimizing InterpolateGridCoStar should report back to optimize_func()...
	 * The fact that FindVCoStar allows interpolation between non-adjoining tracks
	 * should guarantee that this will not happen. */
	if( useTr[0] < 0 )
	{
		fprintf( ioQQQ, " The parameters for the requested CoStar model are out of range.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	ASSERT( useTr[0] >= 0 && useTr[0] < grid->nTracks );
	ASSERT( useTr[1] >= 0 && useTr[1] < grid->nTracks );
	ASSERT( indloTr[useTr[0]] >= 0 && indloTr[useTr[0]] < (int)grid->nmods );
	ASSERT( indhiTr[useTr[0]] >= 0 && indhiTr[useTr[0]] < (int)grid->nmods );
	ASSERT( indloTr[useTr[1]] >= 0 && indloTr[useTr[1]] < (int)grid->nmods );
	ASSERT( indhiTr[useTr[1]] >= 0 && indhiTr[useTr[1]] < (int)grid->nmods );

#	if DEBUGPRT
	printf( "interpolate between tracks %c and %c\n", (char)('A'+useTr[0]), (char)('A'+useTr[1]) );
#	endif

	indlo[0] = indloTr[useTr[0]];
	indhi[0] = indhiTr[useTr[0]];
	indlo[1] = indloTr[useTr[1]];
	indhi[1] = indhiTr[useTr[1]];

	InterpolateModelCoStar( grid, lval, aval, indlo, indhi, index, 0, off, rfield.tslop[rfield.nShape] );

	for( i=0; i < rfield.nupper; i++ )
	{
		rfield.tslop[rfield.nShape][i] = (realnum)pow((realnum)10.f,rfield.tslop[rfield.nShape][i]);
		if( rfield.tslop[rfield.nShape][i] < 1e-37 )
			rfield.tslop[rfield.nShape][i] = 0.;
	}

	if( false )
	{
		FILE *ioBUG = fopen( "interpolated.txt", "w" );
		for( k=0; k < rfield.nupper; ++k )
			fprintf( ioBUG, "%e %e\n", rfield.tNu[rfield.nShape][k].Ryd(), rfield.tslop[rfield.nShape][k] );
		fclose( ioBUG );
	}

	/* sanity check: see whether this model has the correct effective temperature */
	if( ! lgValidModel( rfield.tNu[rfield.nShape], rfield.tslop[rfield.nShape], aval[0], 0.05 ) )
		TotalInsanity();

	/* set limits for optimizer */
	SetLimits( grid, lval[0], NULL, NULL, useTr, ValTr, val0_lo, val0_hi );

	/* now write some final info */
	if( called.lgTalk )
	{
		fprintf( ioQQQ, "                       * c<< FINAL: T_eff = %7.1f, ", aval[0] );
		fprintf( ioQQQ, "log(g) = %4.2f, M(ZAMS) = %5.1f, age = ", aval[1], pow(10.,aval[2]) );
		fprintf( ioQQQ, PrintEfmt("%8.2e",aval[3]) );
		fprintf( ioQQQ, "  >>> *\n" );
	}

	FREE_CHECK( indhiTr );
	FREE_CHECK( indloTr );
	FREE_CHECK( ValTr );
	return;
}

/* find which models to use for interpolation along a given evolutionary track */
STATIC void FindHCoStar(const stellar_grid *grid,
			long track,
			double par2,   /* requested log(g) or age */
			long off,      /* determines which parameter to match 0 -> log(g), 2 -> age */
			realnum *ValTr,/* ValTr[track]: Teff/log(M) value for interpolated model along track */
			long *indloTr, /* indloTr[track]: model number for first model used in interpolation */
			long *indhiTr) /* indhiTr[track]: model number for second model used in interpolation */
{
	long index[2], j, mod1, mod2;

	DEBUG_ENTRY( "FindHCoStar()" );

	indloTr[track] = -2;
	indhiTr[track] = -2;
	ValTr[track] = -FLT_MAX;

	index[1] = track;

	for( j=0; j < grid->trackLen[track]; j++ )
	{
		index[0] = j;
		mod1 = grid->jval[JIndex(grid,index)];

		/* do we have an exact match ? */
		if( fabs(par2-grid->telg[mod1].par[off+1]) <= 10.*FLT_EPSILON*fabs(grid->telg[mod1].par[off+1]) )
		{
			indloTr[track] = mod1;
			indhiTr[track] = mod1;
			ValTr[track] = (realnum)grid->telg[mod1].par[off];
			return;
		}
	}

	for( j=0; j < grid->trackLen[track]-1; j++ )
	{
		index[0] = j;
		mod1 = grid->jval[JIndex(grid,index)];
		index[0] = j+1;
		mod2 = grid->jval[JIndex(grid,index)];

		/* do we interpolate ? */
		if( (par2 - grid->telg[mod1].par[off+1])*(par2 - grid->telg[mod2].par[off+1]) < 0. )
		{
			double frac;

			indloTr[track] = mod1;
			indhiTr[track] = mod2;
			frac = (par2 - grid->telg[mod2].par[off+1])/
				(grid->telg[mod1].par[off+1] - grid->telg[mod2].par[off+1]);
			ValTr[track] = (realnum)(frac*grid->telg[mod1].par[off] + 
				(1.-frac)*grid->telg[mod2].par[off] );
			break;
		}
	}
	return;
}

/* find which tracks to use for interpolation in between tracks */
STATIC void FindVCoStar(const stellar_grid *grid,
			double par1,  /* requested Teff or ZAMS mass */
			realnum *ValTr, /* internal workspace */
			long useTr[]) /* useTr[0]: track number for first track to be used in interpolation
				       *            (i.e., 0 means 'A', etc.)
				       * useTr[1]: track number for second track to be used in interpolation
				       * NOTE: FindVCoStar raises a flag when interpolating between non-adjoining
				       *       tracks, i.e. when (useTr[1]-useTr[0]) > 1 */
{
	long j;

	DEBUG_ENTRY( "FindVCoStar()" );

	useTr[0] = -1;
	useTr[1] = -1;

	for( j=0; j < grid->nTracks; j++ )
	{
		/* do we have an exact match ? */
		if( ValTr[j] != -FLT_MAX && fabs(par1-(double)ValTr[j]) <= 10.*FLT_EPSILON*fabs(ValTr[j]) )
		{
			useTr[0] = j;
			useTr[1] = j;
			break;
		}
	}

	if( useTr[0] >= 0 )
	{
		return;
	}

	for( j=0; j < grid->nTracks-1; j++ )
	{
		if( ValTr[j] != -FLT_MAX )
		{
			long int i,j2;

			/* find next valid track */
			j2 = 0;
			for( i = j+1; i < grid->nTracks; i++ )
			{
				if( ValTr[i] != -FLT_MAX )
				{
					j2 = i;
					break;
				}
			}

			/* do we interpolate ? */
			if( j2 > 0 && ((realnum)par1-ValTr[j])*((realnum)par1-ValTr[j2]) < 0.f )
			{
				useTr[0] = j;
				useTr[1] = j2;
				break;
			}
		}
	}

	/* raise caution when we interpolate between non-adjoining tracks */
	continuum.lgCoStarInterpolationCaution = ( useTr[1]-useTr[0] > 1 );
	return;
}

/* Make a listing of all the models in the CoStar grid */
STATIC void CoStarListModels(const stellar_grid *grid)
{
	long index[2], maxlen, n;

	DEBUG_ENTRY( "CoStarListModels()" );

	maxlen = 0;
	for( n=0; n < grid->nTracks; n++ )
		maxlen = MAX2( maxlen, grid->trackLen[n] );

	fprintf( ioQQQ, "\n" );
	fprintf( ioQQQ, "  Track\\Index |" );
	for( n = 0; n < maxlen; n++ )
		fprintf( ioQQQ, "     %5ld      ", n+1 );
	fprintf( ioQQQ, "\n" );
	fprintf( ioQQQ, "--------------|" );
	for( n = 0; n < maxlen; n++ )
		fprintf( ioQQQ, "----------------" );
	fprintf( ioQQQ, "\n" );

	for( index[1]=0; index[1] < grid->nTracks; ++index[1] )
	{
		long ptr;
		double Teff, alogg, Mass;

		fprintf( ioQQQ, " %c", (char)('A'+index[1]) );
		index[0] = 0;
		ptr = grid->jval[JIndex(grid,index)];
		Mass = pow(10.,grid->telg[ptr].par[2]);
		fprintf( ioQQQ, " (%3.0f Msol) |", Mass );

		for( index[0]=0; index[0] < grid->trackLen[index[1]]; ++index[0] )
		{
			ptr = grid->jval[JIndex(grid,index)];
			Teff = grid->telg[ptr].par[0];
			alogg = grid->telg[ptr].par[1];
			fprintf( ioQQQ, "  (%6.1f,%4.2f)", Teff, alogg );
		}
		fprintf( ioQQQ, "\n" );
	}
	return;
}

/*  RauchInitializeSub does the actual work of preparing the ascii file */
STATIC int RauchInitializeSub(const char chFName[],
			      const char chSuff[],
			      const vector<mpp>& telg,
			      long nmods,
			      long n,
			      long ngrids,
			      const double par2[], /* par2[ngrids] */
			      int format)
{
	char chLine[INPUT_LINE_LENGTH]; /* used for getting input lines */

	FILE *ioOut,  /* pointer to output file we came here to create*/
	  *ioIn;      /* pointer to input files we will read */

	long int i, j;

	double *wavl, *fluxes;

	DEBUG_ENTRY( "RauchInitializeSub()" );

	/* grab some space for the wavelengths and fluxes */
	wavl = (double *)MALLOC( sizeof(double)*NRAUCH);
	fluxes = (double *)MALLOC( sizeof(double)*NRAUCH);

	try
	{
		if( n == 1 )
			ioOut = open_data( chFName, "w", AS_LOCAL_ONLY );
		else
			ioOut = open_data( chFName, "a", AS_LOCAL_ONLY );
	}
	catch( cloudy_exit )
	{
		goto error;
	}

	if( n == 1 )
	{
		fprintf( ioOut, "  %ld\n", VERSION_ASCII );
		fprintf( ioOut, "  %d\n", ( ngrids == 1 ? 2 : 3 ) );
		fprintf( ioOut, "  %d\n", ( ngrids == 1 ? 2 : 3 ) );
		fprintf( ioOut, "  Teff\n" );
		fprintf( ioOut, "  log(g)\n" );
		if( ngrids == 2 )
			fprintf( ioOut, "  log(Z)\n" );
		else if( ngrids == 11 )
			fprintf( ioOut, "  f(He)\n" );
		/* NB - this is based on the assumption that each of the planes in the cubic grid is the same */
		fprintf( ioOut, "  %ld\n", nmods*ngrids );
		fprintf( ioOut, "  %d\n", NRAUCH );
		/* Rauch models give the wavelength in Angstrom */
		fprintf( ioOut, "  lambda\n" );
		/* conversion factor to Angstrom */
		fprintf( ioOut, "  %.8e\n", 1. );
		/* Rauch models give the "Astrophysical" flux F_lambda in erg/cm^2/s/cm */
		fprintf( ioOut, "  F_lambda\n" );
		/* the factor PI*1e-8 is needed to convert to "regular" flux in erg/cm^2/s/Angstrom */
		fprintf( ioOut, "  %.8e\n", PI*1.e-8 );
		/* NB - this is based on the assumption that each of the planes in the cubic grid is the same */
		for( j=0; j < ngrids; j++ )
		{
			/* write out the {Teff,log(g)} grid */
			for( i=0; i < nmods; i++ )
			{
				if( ngrids == 1 )
					fprintf( ioOut, "  %.0f %.1f", telg[i].par[0], telg[i].par[1] );
				else
					fprintf( ioOut, "  %.0f %.1f %.1f", telg[i].par[0], telg[i].par[1], par2[j] );
				if( ((i+1)%4) == 0 )
					fprintf( ioOut, "\n" );
			}
			if( (i%4) != 0 )
				fprintf( ioOut, "\n" );
		}

		fprintf( ioQQQ, " Writing: " );
	}

	for( i=0; i < nmods; i++ )
	{
		/* must create name of next stellar atmosphere */
		if( format == 1 )
			sprintf( chLine, "%7.7ld_%2ld", (long)(telg[i].par[0]+0.5), (long)(10.*telg[i].par[1]+0.5) );
		else if( format == 2 )
			sprintf( chLine, "%7.7ld_%.2f", (long)(telg[i].par[0]+0.5), telg[i].par[1] );
		else
		{
			fprintf( ioQQQ, " insanity in RauchInitializeSub\n" );
			ShowMe();
			cdEXIT(EXIT_FAILURE);
		}
		string chFileName( chLine );
		chFileName += chSuff;
		/* now open next stellar atmosphere for reading*/
		try
		{
			ioIn = open_data( chFileName.c_str(), "r", AS_LOCAL_ONLY );
		}
		catch( cloudy_exit )
		{
			goto error;
		}

		/* get first line */
		j = 0;
		if( read_whole_line( chLine, (int)sizeof(chLine), ioIn ) == NULL )
		{
			fprintf( ioQQQ, " RauchInitializeSub error in atmosphere file %4ld%4ld\n", 
				 i, j );
			goto error;
		}
		/* >>chng 02 nov 20, now keep reading them until don't hit the *
		 * since number of comments may change */
		while( chLine[0] == '*' )
		{
			if( read_whole_line( chLine, (int)sizeof(chLine), ioIn ) == NULL )
			{
				fprintf( ioQQQ, " RauchInitializeSub error in atmosphere file %4ld%4ld\n", 
					 i, j );
				goto error;
			}
			++j;
		}

		for( j=0; j < NRAUCH; j++ )
		{
			double ttemp, wl;
			/* get the input line */
			/* >>chng 02 nov 20, don't reread very first line image since we got it above */
			if( j > 0 )
			{
				if(read_whole_line( chLine, (int)sizeof(chLine), ioIn )==NULL )
				{
					fprintf( ioQQQ, " RauchInitializeSub error in atmosphere file %4ld%4ld\n", 
						 i, j );
					goto error;
				}
			}

			/* scan off wavelength and flux)*/
			if( sscanf( chLine, "%lf %le", &wl, &ttemp ) != 2 )
			{
				fprintf( ioQQQ, " RauchInitializeSub error in atmosphere file %4ld%4ld\n", 
					 i, j );
				goto error;
			}

			if( i == 0 )
				wavl[j] = wl;
			else
			{
				/* check if this model is on the same wavelength grid as the first */
				if( !fp_equal(wavl[j],wl,10) )
				{
					fprintf( ioQQQ, " RauchInitializeSub error in atmosphere file %4ld%4ld\n", 
						 i, j );
					goto error;
				}
			}
			fluxes[j] = ttemp; 
		}

		/* finished - close the unit */
		fclose(ioIn);

		/* now write to output file */
		if( i == 0 && n == 1 )
		{
			/* wavelength grid is the same for all models, so write only once */
			for( j=0; j < NRAUCH; j++ )
			{
				fprintf( ioOut, "  %.4e", wavl[j] );
				/* want to have 5 numbers per line */
				if( ((j+1)%5) == 0 )
					fprintf( ioOut, "\n" );
			}
			/* need to throw a newline if we did not end on an exact line */
			if( (j%5) != 0 )
				fprintf( ioOut, "\n" );
		}

		for( j=0; j < NRAUCH; j++ )
		{
			fprintf( ioOut, "  %.4e", fluxes[j] );
			/* want to have 5 numbers per line */
			if( ((j+1)%5) == 0 )
				fprintf( ioOut, "\n" );
		}
		/* need to throw a newline if we did not end on an exact line */
		if( (j%5) != 0 )
			fprintf( ioOut, "\n" );

		/* print to screen to show progress */
		fprintf( ioQQQ, "." );
		fflush( ioQQQ );
	}

	if( n == ngrids )
		fprintf( ioQQQ, " done.\n" );

	fclose(ioOut);

	/* free the space grabbed above */
	FREE_CHECK( fluxes );
	FREE_CHECK( wavl );
	return 0;

error:
	FREE_CHECK( fluxes );
	FREE_CHECK( wavl );
	return 1;
}

STATIC void RauchReadMPP(vector<mpp>& telg1,
			 vector<mpp>& telg2,
			 vector<mpp>& telg3,
			 vector<mpp>& telg4,
			 vector<mpp>& telg5,
			 vector<mpp>& telg6)
{
	DEBUG_ENTRY( "RauchReadMPP()" );

	const char fnam[] = "rauch_models.dat";
	fstream ioDATA;
	open_data( ioDATA, fnam, mode_r );

	string line;
	getdataline( ioDATA, line );
	long version;
	istringstream iss( line );
	iss >> version;
	if( version != VERSION_RAUCH_MPP )
	{
		fprintf( ioQQQ, " RauchReadMPP: the version of %s is not the current version.\n", fnam );
		fprintf( ioQQQ, " Please obtain the current version from the Cloudy web site.\n" );
		fprintf( ioQQQ, " I expected to find version %ld and got %ld instead.\n",
			 VERSION_RAUCH_MPP, version );
		cdEXIT(EXIT_FAILURE);
	}

	getdataline( ioDATA, line );
	unsigned long ndata;
	istringstream iss2( line );
	iss2 >> ndata;
	ASSERT( ndata == telg1.size() );
	// this implicitly assumes there is exactly one comment line between
	// the number of data points and the start of the data
	getline( ioDATA, line );
	// read data for H-Ca grid
	for( unsigned long i=0; i < ndata; ++i )
		ioDATA >> telg1[i].par[0] >> telg1[i].par[1];
	getline( ioDATA, line );
		
	getdataline( ioDATA, line );
	istringstream iss3( line );
	iss3 >> ndata;
	ASSERT( ndata == telg2.size() );
	getline( ioDATA, line );
	// read data for H-Ni grid
	for( unsigned long i=0; i < ndata; ++i )
		ioDATA >> telg2[i].par[0] >> telg2[i].par[1];
	getline( ioDATA, line );
		
	getdataline( ioDATA, line );
	istringstream iss4( line );
	iss4 >> ndata;
	ASSERT( ndata == telg3.size() );
	getline( ioDATA, line );
	// read data for PG1159 grid
	for( unsigned long i=0; i < ndata; ++i )
		ioDATA >> telg3[i].par[0] >> telg3[i].par[1];
	getline( ioDATA, line );
		
	getdataline( ioDATA, line );
	istringstream iss5( line );
	iss5 >> ndata;
	ASSERT( ndata == telg4.size() );
	getline( ioDATA, line );
	// read data for pure H grid
	for( unsigned long i=0; i < ndata; ++i )
		ioDATA >> telg4[i].par[0] >> telg4[i].par[1];
	getline( ioDATA, line );
		
	getdataline( ioDATA, line );
	istringstream iss6( line );
	iss6 >> ndata;
	ASSERT( ndata == telg5.size() );
	getline( ioDATA, line );
	// read data for pure He grid
	for( unsigned long i=0; i < ndata; ++i )
		ioDATA >> telg5[i].par[0] >> telg5[i].par[1];
	getline( ioDATA, line );
		
	getdataline( ioDATA, line );
	istringstream iss7( line );
	iss7 >> ndata;
	ASSERT( ndata == telg6.size() );
	getline( ioDATA, line );
	// read data for pure H+He grid
	for( unsigned long i=0; i < ndata; ++i )
		ioDATA >> telg6[i].par[0] >> telg6[i].par[1];
	getline( ioDATA, line );
		
	getdataline( ioDATA, line );
	istringstream iss8( line );
	iss8 >> version;
	ASSERT( version == VERSION_RAUCH_MPP );

	return;
}

inline void getdataline(fstream& ioDATA,
			string& line)
{
	do
	{
		getline( ioDATA, line );
	}
	while( line[0] == '#' );
	return;
}

/* lgCompileAtmosphere does the actual rebinning onto the Cloudy grid and writes the binary file */
/* >>chng 01 feb 12, added return value to indicate success (0) or failure (1) */
STATIC bool lgCompileAtmosphere(const char chFNameIn[],
				const char chFNameOut[],
				const realnum Edges[], /* Edges[nedges] */
				long nedges,
				process_counter& pc)
{
	FILE *ioIN;  /* used for input */
	FILE *ioOUT; /* used for output */

	char chDataType[11];
	char names[MDIM][MNAM+1];

	bool lgFreqX, lgFreqY, lgFlip;
	int32 val[7];
	uint32 uval[2];
	double dval[3];
	char md5sum[NMD5];
	long int i, imod, version, nd, ndim, npar, nmods, ngrid;

	/* these will be malloced into large work arrays */
	realnum *StarEner = NULL, *StarFlux = NULL, *CloudyFlux = NULL, *scratch = NULL;

	double convert_wavl, convert_flux;

	mpp *telg = NULL;

	DEBUG_ENTRY( "lgCompileAtmosphere()" );

	try
	{
		ioIN = open_data( chFNameIn, "r", AS_LOCAL_ONLY );
	}
	catch( cloudy_exit )
	{
		goto error;
	}
	fprintf( ioQQQ, " lgCompileAtmosphere got %s.\n", chFNameIn );

	/* read version number */
	if( fscanf( ioIN, "%ld", &version ) != 1 )
	{
		fprintf( ioQQQ, " lgCompileAtmosphere failed reading VERSION.\n" );
		goto error;
	}

	if( version != VERSION_ASCII )
	{
		fprintf( ioQQQ, " lgCompileAtmosphere: there is a version number mismatch in"
			 " the ascii atmosphere file: %s.\n", chFNameIn );
		fprintf( ioQQQ, " lgCompileAtmosphere: Please recreate this file or download the"
			 " latest version following the instructions on the Cloudy website.\n" );
		goto error;
	}

	/* >>chng 06 jun 10, read the dimension of the grid, PvH */
	if( fscanf( ioIN, "%ld", &ndim ) != 1 || ndim <= 0 || ndim > MDIM )
	{
		fprintf( ioQQQ, " lgCompileAtmosphere failed reading valid dimension of grid.\n" );
		goto error;
	}

	/* >>chng 06 jun 12, read the number of model parameters, PvH */
	if( fscanf( ioIN, "%ld", &npar ) != 1 || npar <= 0 || npar < ndim || npar > MDIM )
	{
		fprintf( ioQQQ, " lgCompileAtmosphere failed reading valid no. of model parameters.\n" );
		goto error;
	}

	/* make sure valgrind doesn't trip on the binary write of this array */
	memset( names, '\0', MDIM*(MNAM+1) );

	for( nd=0; nd < npar; nd++ )
	{
		if( fscanf( ioIN, "%6s", names[nd] ) != 1 )
		{
			fprintf( ioQQQ, " lgCompileAtmosphere failed reading parameter label.\n" );
			goto error;
		}
	}

	/* >>chng 05 nov 18, read the following extra parameters from the ascii file, PvH */
	if( fscanf( ioIN, "%ld", &nmods ) != 1 || nmods <= 0 )
	{
		fprintf( ioQQQ, " lgCompileAtmosphere failed reading valid number of models.\n" );
		goto error;
	}

	if( fscanf( ioIN, "%ld", &ngrid ) != 1 || ngrid <= 1 )
	{
		fprintf( ioQQQ, " lgCompileAtmosphere failed reading valid number of grid points.\n" );
		goto error;
	}

	/* read data type for wavelengths, allowed values are lambda, nu */
	if( fscanf( ioIN, "%10s", chDataType ) != 1 )
	{
		fprintf( ioQQQ, " lgCompileAtmosphere failed reading wavl DataType string.\n" );
		goto error;
	}

	if( strcmp( chDataType, "lambda" ) == 0 )
		lgFreqX = false;
	else if( strcmp( chDataType, "nu" ) == 0 )
		lgFreqX = true;
	else {
		fprintf( ioQQQ, " lgCompileAtmosphere found illegal wavl DataType: %s.\n", chDataType );
		goto error;
	}

	if( fscanf( ioIN, "%le", &convert_wavl ) != 1 || convert_wavl <= 0. )
	{
		fprintf( ioQQQ, " lgCompileAtmosphere failed reading valid wavl conversion factor.\n" );
		goto error;
	}

	/* read data type for flux, allowed values F_lambda, H_lambda, F_nu, H_nu */
	if( fscanf( ioIN, "%10s", chDataType ) != 1 )
	{
		fprintf( ioQQQ, " lgCompileAtmosphere failed reading flux DataType string.\n" );
		goto error;
	}

	if( strcmp( chDataType, "F_lambda" ) == 0 || strcmp( chDataType, "H_lambda" ) == 0 )
		lgFreqY = false;
	else if( strcmp( chDataType, "F_nu" ) == 0 || strcmp( chDataType, "H_nu" ) == 0 )
		lgFreqY = true;
	else {
		fprintf( ioQQQ, " lgCompileAtmosphere found illegal flux DataType: %s.\n", chDataType );
		goto error;
	}

	if( fscanf( ioIN, "%le", &convert_flux ) != 1 || convert_flux <= 0. )
	{
		fprintf( ioQQQ, " lgCompileAtmosphere failed reading valid flux conversion factor.\n" );
		goto error;
	}

	telg = (mpp *)CALLOC( (size_t)nmods, sizeof(mpp) );

	for( i=0; i < nmods; i++ )
	{
		for( nd=0; nd < npar; nd++ )
		{
			if( fscanf( ioIN, "%le", &telg[i].par[nd] ) != 1 )
			{
				fprintf( ioQQQ, " lgCompileAtmosphere failed reading valid model parameter.\n" );
				goto error;
			}
		}
		if( telg[i].par[0] <= 0. )
		{
			fprintf( ioQQQ, " lgCompileAtmosphere failed reading valid %s.\n", names[0] );
			goto error;
		}
	}

	try
	{
		ioOUT = open_data( chFNameOut, "wb", AS_LOCAL_ONLY );
	}
	catch( cloudy_exit )
	{
		goto error;
	}

	val[0] = (int32)VERSION_BIN;
	val[1] = (int32)MDIM;
	val[2] = (int32)MNAM;
	val[3] = (int32)ndim;
	val[4] = (int32)npar;
	val[5] = (int32)nmods;
	val[6] = (int32)rfield.nupper;
	uval[0] = sizeof(val) + sizeof(uval) + sizeof(dval) + sizeof(md5sum) +
		sizeof(names) + nmods*sizeof(mpp); /* nOffset */
	uval[1] = rfield.nupper*sizeof(realnum); /* nBlocksize */
	dval[0] = double(rfield.emm);
	dval[1] = double(rfield.egamry);
	dval[2] = double(continuum.ResolutionScaleFactor);

	strncpy( md5sum, continuum.mesh_md5sum.c_str(), sizeof(md5sum) );

	if( fwrite( val, sizeof(val), 1, ioOUT ) != 1 ||
	    fwrite( uval, sizeof(uval), 1, ioOUT ) != 1 ||
	    /* write out the lower, upper bound of the energy mesh, and the res scale factor */
	    fwrite( dval, sizeof(dval), 1, ioOUT ) != 1 ||
	    /* write out the (modified) md5 checksum of continuum_mesh.ini */
	    fwrite( md5sum, sizeof(md5sum), 1, ioOUT ) != 1 ||
	    fwrite( names, sizeof(names), 1, ioOUT ) != 1 ||
	    /* write out the array of {Teff,log(g)} pairs */
	    fwrite( telg, sizeof(mpp), (size_t)nmods, ioOUT ) != (size_t)nmods ||
	    /* write out the cloudy energy grid for later sanity checks */
	    fwrite( rfield.AnuOrg, (size_t)uval[1], 1, ioOUT ) != 1 )
	{
		fprintf( ioQQQ, " lgCompileAtmosphere failed writing header of output file.\n" );
		goto error;
	}

	/* MALLOC some workspace */
	StarEner = (realnum *)MALLOC( sizeof(realnum)*ngrid );
	scratch = (realnum *)MALLOC( sizeof(realnum)*ngrid );
	StarFlux = (realnum *)MALLOC( sizeof(realnum)*ngrid );
	CloudyFlux = (realnum *)MALLOC( (size_t)uval[1] );

	/* read wavelength grid */
	for( i=0; i < ngrid; i++ )
	{
		double help;
		if( fscanf( ioIN, "%lg", &help ) != 1 )
		{
			fprintf( ioQQQ, " lgCompileAtmosphere failed reading wavelength.\n" );
			goto error;
		}
		/* this conversion makes sure that scratch[i] is
		 * either wavelength in Angstrom or frequency in Hz */
		scratch[i] = (realnum)(help*convert_wavl);

		ASSERT( scratch[i] > 0.f );
	}

	lgFlip = ( !lgFreqX && scratch[0] < scratch[1] ) || ( lgFreqX && scratch[0] > scratch[1] );

	/* convert continuum over to increasing frequency in Ryd */
	for( i=0; i < ngrid; i++ )
	{
		/* convert scratch[i] to frequency in Ryd */
		if( lgFreqX )
			scratch[i] /= (realnum)FR1RYD;
		else
			scratch[i] = (realnum)(RYDLAM/scratch[i]);

		if( lgFlip )
			StarEner[ngrid-i-1] = scratch[i];
		else
			StarEner[i] = scratch[i];
	}

	ASSERT( StarEner[0] > 0.f );
	/* make sure the array is in ascending order */
	for( i=1; i < ngrid; i++ )
		ASSERT( StarEner[i] > StarEner[i-1] );

	fprintf( ioQQQ, " Compiling: " );

	for( imod=0; imod < nmods; imod++ )
	{
		const realnum CONVERT_FNU = (realnum)(1.e8*SPEEDLIGHT/POW2(FR1RYD));

		/* now read the stellar fluxes */
		for( i=0; i < ngrid; i++ )
		{
			double help;
			if( fscanf( ioIN, "%lg", &help ) != 1 )
			{
				fprintf( ioQQQ, " lgCompileAtmosphere failed reading star flux.\n" );
				goto error;
			}
			/* this conversion makes sure that scratch[i] is either
			 * F_nu in erg/cm^2/s/Hz or F_lambda in erg/cm^2/s/A */
			scratch[i] = (realnum)(help*convert_flux);

			/* this can underflow on the Wien tail */
			ASSERT( scratch[i] >= 0.f );
		}

		for( i=0; i < ngrid; i++ )
		{
			if( lgFlip )
				StarFlux[ngrid-i-1] = scratch[i];
			else
				StarFlux[i] = scratch[i];
		}

		for( i=0; i < ngrid; i++ )
		{
			/* this converts to F_nu in erg/cm^2/s/Hz */
			if( !lgFreqY )
				StarFlux[i] *= CONVERT_FNU/POW2(StarEner[i]);
			ASSERT( StarFlux[i] >= 0.f );
		}

		/* the re-binned values are returned in the "CloudyFlux" array */
		RebinAtmosphere( ngrid, StarEner, StarFlux, CloudyFlux, nedges, Edges );

		/* write the continuum out as a binary file */
		if( fwrite( CloudyFlux, (size_t)uval[1], 1, ioOUT ) != 1 )
		{
			fprintf( ioQQQ, " lgCompileAtmosphere failed writing star flux.\n" );
			goto error;
		}

		fprintf( ioQQQ, "." );
		fflush( ioQQQ );
	}

	fprintf( ioQQQ, " done.\n" );

	fclose(ioIN);
	fclose(ioOUT);

	/* now free up the memory we claimed */
	FREE_CHECK( CloudyFlux );
	FREE_CHECK( StarFlux );
	FREE_CHECK( StarEner );
	FREE_CHECK( scratch );
	FREE_CHECK( telg );

	fprintf( ioQQQ, " lgCompileAtmosphere completed ok.\n\n" );
	++pc.nOK;
	return false;

error:
	FREE_SAFE( CloudyFlux );
	FREE_SAFE( StarFlux );
	FREE_SAFE( StarEner );
	FREE_SAFE( scratch );
	FREE_SAFE( telg );
	++pc.nFail;
	return true;
}

STATIC void InitGrid(stellar_grid *grid,
		     bool lgList)
{
	long nd;
	int32 version, mdim, mnam;
	double mesh_elo, mesh_ehi;
	char md5sum[NMD5];

	DEBUG_ENTRY( "InitGrid()" );

	try
	{
		grid->ioIN = open_data( grid->name.c_str(), "rb", grid->scheme );
	}
	catch( cloudy_exit )
	{
		/* something went wrong */
		/* NB NB - DO NOT CHANGE THE FOLLOWING ERROR MESSAGE! checkall.pl picks it up */
		fprintf( ioQQQ, " Error: stellar atmosphere file not found.\n" );
		fprintf(ioQQQ , "\n\n If the path is set then it is possible that the stellar atmosphere data files do not exist.\n");
		fprintf(ioQQQ , " Have the stellar data files been downloaded and compiled with the COMPILE STARS command?\n");
		fprintf(ioQQQ , " If you are simply running the test suite and do not need the stellar continua then you should simply ignore this failure\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* >>chng 01 oct 17, add version and size to this array */
	if( fread( &version, sizeof(version), 1, grid->ioIN ) != 1 ||
	    fread( &mdim, sizeof(mdim), 1, grid->ioIN ) != 1 ||
	    fread( &mnam, sizeof(mnam), 1, grid->ioIN ) != 1 ||
	    fread( &grid->ndim, sizeof(grid->ndim), 1, grid->ioIN ) != 1 ||
	    fread( &grid->npar, sizeof(grid->npar), 1, grid->ioIN ) != 1 ||
	    fread( &grid->nmods, sizeof(grid->nmods), 1, grid->ioIN ) != 1 ||
	    fread( &grid->ngrid, sizeof(grid->ngrid), 1, grid->ioIN ) != 1 ||
	    fread( &grid->nOffset, sizeof(grid->nOffset), 1, grid->ioIN ) != 1 ||
	    fread( &grid->nBlocksize, sizeof(grid->nBlocksize), 1, grid->ioIN ) != 1 ||
	    fread( &mesh_elo, sizeof(mesh_elo), 1, grid->ioIN ) != 1 ||
	    fread( &mesh_ehi, sizeof(mesh_ehi), 1, grid->ioIN ) != 1 ||
	    fread( &rfield.RSFCheck[rfield.nShape], sizeof(rfield.RSFCheck[rfield.nShape]), 1, grid->ioIN ) != 1 ||
	    fread( md5sum, sizeof(md5sum), 1, grid->ioIN ) != 1 )
	{
		fprintf( ioQQQ, " InitGrid failed reading header.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* do some sanity checks */
	if( version != VERSION_BIN )
	{
		fprintf( ioQQQ, " InitGrid: there is a version mismatch between"
			 " the compiled atmospheres file I expected and the one I found.\n" );
		fprintf( ioQQQ, " InitGrid: Please recompile the stellar"
			 " atmospheres file with the command: %s.\n", grid->command );
		cdEXIT(EXIT_FAILURE);
	}

	if( mdim != MDIM || mnam != MNAM )
	{
		fprintf( ioQQQ, " InitGrid: the compiled atmospheres file is produced"
			 " with an incompatible version of Cloudy.\n" );
		fprintf( ioQQQ, " InitGrid: Please recompile the stellar"
			 " atmospheres file with the command: %s.\n", grid->command );
		cdEXIT(EXIT_FAILURE);
	}

	if( !fp_equal( double(rfield.emm), mesh_elo ) ||
	    !fp_equal( double(rfield.egamry), mesh_ehi ) ||
	    strncmp( continuum.mesh_md5sum.c_str(), md5sum, NMD5 ) != 0 )
	{
		fprintf( ioQQQ, " InitGrid: the compiled atmospheres file is produced"
			 " with an incompatible frequency grid.\n" );
		fprintf( ioQQQ, " InitGrid: Please recompile the stellar"
			 " atmospheres file with the command: %s.\n", grid->command );
		cdEXIT(EXIT_FAILURE);
	}

	ASSERT( grid->ndim > 0 && grid->ndim <= MDIM );
	ASSERT( grid->npar >= grid->ndim && grid->npar <= MDIM );
	ASSERT( grid->nmods > 0 );
	ASSERT( grid->ngrid > 0 );
	ASSERT( grid->nOffset > 0 );
	ASSERT( grid->nBlocksize > 0 );

	rfield.nupper = grid->ngrid;

	if( fread( &grid->names, sizeof(grid->names), 1, grid->ioIN ) != 1 )
	{
		fprintf( ioQQQ, " InitGrid failed reading names array.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	grid->lgIsTeffLoggGrid = ( grid->ndim >= 2 &&
				   strcmp( grid->names[0], "Teff" ) == 0 &&
				   strcmp( grid->names[1], "log(g)" ) == 0 );

	grid->telg = (mpp *)MALLOC( sizeof(mpp)*grid->nmods );
	grid->val = (double **)MALLOC( sizeof(double*)*grid->ndim );
	for( nd=0; nd < grid->ndim; nd++ )
	{
		grid->val[nd] = (double *)MALLOC( sizeof(double)*grid->nmods );
	}
	grid->nval = (long *)MALLOC( sizeof(long)*grid->ndim );

	if( fread( grid->telg, sizeof(mpp), grid->nmods, grid->ioIN ) != (size_t)grid->nmods )
	{
		fprintf( ioQQQ, " InitGrid failed reading model parameter block.\n" );
		cdEXIT(EXIT_FAILURE);
	}

#	ifdef SEEK_END
	/* sanity check: does the file have the correct length ? */
	/* NOTE: this operation is not necessarily supported by all operating systems
	 * but if the preprocessor symbol SEEK_END exists it is assumed to be supported */
	int res = fseek( grid->ioIN, 0, SEEK_END );
	if( res == 0 )
	{
		long End = ftell( grid->ioIN );
		long Expected = grid->nOffset + (grid->nmods+1)*grid->nBlocksize;
		if( End != Expected )
		{
			fprintf( ioQQQ, " InitGrid: Problem performing sanity check for size of binary file.\n" );
			fprintf( ioQQQ, " InitGrid: I expected to find %ld bytes, but actually found %ld bytes.\n",
				 Expected, End );
			fprintf( ioQQQ, " InitGrid: Please recompile the stellar"
				 " atmospheres file with the command: %s.\n", grid->command );
			cdEXIT(EXIT_FAILURE);
		}
	}
#	endif

	InitIndexArrays( grid, lgList );

	/* set default interpolation mode */
	grid->imode = IM_RECT_GRID;
	/* these are only used by CoStar grids */
	grid->trackLen = NULL;
	grid->nTracks = 0;
	grid->jval = NULL;
	return;
}

/* check whether a binary atmosphere exists and is up-to-date */
STATIC bool lgValidBinFile(const char *binName, process_counter& pc, access_scheme scheme)
{
	int32 version, mdim, mnam;
	double mesh_elo, mesh_ehi, mesh_res_factor;
	char md5sum[NMD5];
	stellar_grid grid;

	DEBUG_ENTRY( "lgValidBinFile()" );

	//
	// NB NB NB
	//
	// this routine is called when either of these two commands is issued:
	// 
	// TABLE STAR AVAIL
	// COMPILE STAR [ additional parameters ]
	//
	// both these commands execute as soon as they are parsed and then terminate
	// hence it is safe to assume that no SET CONTINUUM RESOLUTION command will follow!
	//
	// !!! THE TEST BELOW FOR VALIDITY OF THE FILE DEPENDS ON THAT ASSUMPTION !!!
	//

	grid.name = binName;

	if( (grid.ioIN = open_data( grid.name.c_str(), "rb", scheme )) == NULL )
		return false;

	if( fread( &version, sizeof(version), 1, grid.ioIN ) != 1 ||
	    fread( &mdim, sizeof(mdim), 1, grid.ioIN ) != 1 ||
	    fread( &mnam, sizeof(mnam), 1, grid.ioIN ) != 1 ||
	    fread( &grid.ndim, sizeof(grid.ndim), 1, grid.ioIN ) != 1 ||
	    fread( &grid.npar, sizeof(grid.npar), 1, grid.ioIN ) != 1 ||
	    fread( &grid.nmods, sizeof(grid.nmods), 1, grid.ioIN ) != 1 ||
	    fread( &grid.ngrid, sizeof(grid.ngrid), 1, grid.ioIN ) != 1 ||
	    fread( &grid.nOffset, sizeof(grid.nOffset), 1, grid.ioIN ) != 1 ||
	    fread( &grid.nBlocksize, sizeof(grid.nBlocksize), 1, grid.ioIN ) != 1 ||
	    fread( &mesh_elo, sizeof(mesh_elo), 1, grid.ioIN ) != 1 ||
	    fread( &mesh_ehi, sizeof(mesh_ehi), 1, grid.ioIN ) != 1 ||
	    fread( &mesh_res_factor, sizeof(mesh_res_factor), 1, grid.ioIN ) != 1 ||
	    fread( md5sum, sizeof(md5sum), 1, grid.ioIN ) != 1 )
	{
		fclose( grid.ioIN );
		return false;
	}

	/* do some sanity checks */
	if( version != VERSION_BIN || mdim != MDIM || mnam != MNAM ||
	    !fp_equal( double(rfield.emm), mesh_elo ) ||
	    !fp_equal( double(rfield.egamry), mesh_ehi ) ||
	    !fp_equal( double(continuum.ResolutionScaleFactor), mesh_res_factor ) ||
	    strncmp( continuum.mesh_md5sum.c_str(), md5sum, NMD5 ) != 0 )
	{
		fclose( grid.ioIN );
		return false;
	}

#	ifdef SEEK_END
	/* sanity check: does the file have the correct length ? */
	/* NOTE: this operation is not necessarily supported by all operating systems
	 * but if the preprocessor symbol SEEK_END exists it is assumed to be supported */
	int res = fseek( grid.ioIN, 0, SEEK_END );
	if( res == 0 )
	{
		long End = ftell( grid.ioIN );
		long Expected = grid.nOffset + (grid.nmods+1)*grid.nBlocksize;
		if( End != Expected )
		{
			fclose( grid.ioIN );
			return false;
		}
	}
#	endif

	fclose( grid.ioIN );
	++pc.notProcessed; // the file is up-to-date -> no processing 
	return true;
}

/* check whether a ascii atmosphere file exists and is up-to-date */
STATIC bool lgValidAsciiFile(const char *ascName, access_scheme scheme)
{
	long version;
	FILE *ioIN;

	DEBUG_ENTRY( "lgValidAsciiFile()" );

	/* can we read the file? */
	if( (ioIN = open_data( ascName, "r", scheme )) == NULL )
		return false;

	/* check version number */
	if( fscanf( ioIN, "%ld", &version ) != 1 || version != VERSION_ASCII )
	{
		fclose( ioIN );
		return false;
	}

	fclose( ioIN );
	return true;
}

/* sort CoStar models according to track and index number, store indices in grid->jval[] */
STATIC void InitGridCoStar(stellar_grid *grid) /* the grid parameters */
{
	char track;
	bool lgFound;
	long i, index[2];

	DEBUG_ENTRY( "InitGridCoStar()" );

	ASSERT( grid->ndim == 2 );
	ASSERT( grid->jlo != NULL && grid->jhi != NULL );

	grid->jval = grid->jlo;
	FREE_CHECK( grid->jhi );
	grid->jlo = grid->jhi = NULL;

	/* invalidate contents set by InitGrid first */
	memset( grid->jval, 0xff, (size_t)(grid->nval[0]*grid->nval[1]*sizeof(long)) );

	grid->trackLen = (long *)CALLOC( (size_t)grid->nmods, sizeof(long) );

	index[1] = 0;
	while( true )
	{
		index[0] = 0;
		track = (char)('A'+index[1]);
		do
		{
			lgFound = false;
			for( i=0; i < grid->nmods; i++ )
			{
				if( grid->telg[i].chGrid == track && grid->telg[i].modid == index[0]+1 )
				{
					grid->jval[JIndex(grid,index)] = i;
					++index[0];
					lgFound = true;
					break;
				}
			}
		}
		while( lgFound );

		if( index[0] == 0 )
			break;

		grid->trackLen[index[1]] = index[0];
		++index[1];
	}

	grid->nTracks = index[1];
	return;
}

STATIC void CheckVal(const stellar_grid *grid,
		     double val[], /* val[ndim] */
		     long *nval,
		     long *ndim)
{
	DEBUG_ENTRY( "CheckVal()" );

	if( *ndim == 0 )
		*ndim = (long)grid->ndim;
	if( *ndim == 2 && *nval == 1 && grid->lgIsTeffLoggGrid )
	{
		/* default gravity is maximum gravity */
		val[*nval] = grid->val[1][grid->nval[1]-1];
		++(*nval);
	}
	if( *ndim != (long)grid->ndim )
	{
		fprintf( ioQQQ, " A %ld-dim grid was requested, but a %ld-dim grid was found.\n",
			 *ndim, (long)grid->ndim );
		cdEXIT(EXIT_FAILURE);
	}
	if( *nval < *ndim )
	{
		fprintf( ioQQQ, " A %ld-dim grid was requested, but only %ld parameters were entered.\n",
			 *ndim, *nval );
		cdEXIT(EXIT_FAILURE);
	}
}

STATIC void InterpolateRectGrid(const stellar_grid *grid,
				const double val[], /* val[ndim] */
				double *Tlow,
				double *Thigh)
{
	bool lgInvalid;
	long int i,
	  *indlo,
	  *indhi,
	  *index,
	  k,
	  nd;
	double *aval;

	DEBUG_ENTRY( "InterpolateRectGrid()" );

	/* create some space */
	indlo = (long *)MALLOC((size_t)(grid->ndim*sizeof(long)) );
	indhi = (long *)MALLOC((size_t)(grid->ndim*sizeof(long)) );
	index = (long *)MALLOC((size_t)(grid->ndim*sizeof(long)) );
	aval = (double *)MALLOC((size_t)(grid->npar*sizeof(double)) );

	ASSERT( rfield.lgContMalloc[rfield.nShape] );
	ASSERT( grid->nBlocksize == rfield.nupper*sizeof(realnum) );

	/* save energy scale for check against code's in conorm (scale not yet defined when this routine called) */
	GetBins( grid, rfield.tNu[rfield.nShape] );

#	if DEBUGPRT
	/* check whether the models have the correct effective temperature, for debugging only */
	ValidateGrid( grid, 0.02 );
#	endif

	/* now generate pointers for models to use */
	for( nd=0; nd < grid->ndim; nd++ )
	{
		FindIndex( grid->val[nd], grid->nval[nd], val[nd], &indlo[nd], &indhi[nd], &lgInvalid );
		if( lgInvalid )
		{
			fprintf( ioQQQ, 
				 " Requested parameter %s = %.2f is not within the range %.2f to %.2f\n",
				 grid->names[nd], val[nd], grid->val[nd][0], grid->val[nd][grid->nval[nd]-1] );
			cdEXIT(EXIT_FAILURE);
		}
	}

	InterpolateModel( grid, val, aval, indlo, indhi, index, grid->ndim, rfield.tslop[rfield.nShape], IS_UNDEFINED );

	/* print the parameters of the interpolated model */
	if( called.lgTalk )
	{
		if( grid->npar == 1 )
			fprintf( ioQQQ, 
				 "                       * c<< FINAL:  %6s = %13.2f"
				 "                                          >>> *\n", 
				 grid->names[0], aval[0] );
		else if( grid->npar == 2 )
			fprintf( ioQQQ, 
				 "                       * c<< FINAL:  %6s = %10.2f"
				 "   %6s = %8.5f                         >>> *\n", 
				 grid->names[0], aval[0], grid->names[1], aval[1] );
		else if( grid->npar == 3 )
			fprintf( ioQQQ, 
				 "                       * c<< FINAL:  %6s = %7.0f"
				 "   %6s = %5.2f   %6s = %5.2f              >>> *\n", 
				 grid->names[0], aval[0], grid->names[1], aval[1],
				 grid->names[2], aval[2] );
		else if( grid->npar >= 4 )
		{
			fprintf( ioQQQ, 
				 "                       * c<< FINAL:  %4s = %7.0f"
				 " %6s = %4.2f %6s = %5.2f %6s = ",
				 grid->names[0], aval[0], grid->names[1], aval[1],
				 grid->names[2], aval[2], grid->names[3] );
			fprintf( ioQQQ, PrintEfmt( "%9.2e", aval[3] ) );
			fprintf( ioQQQ, "  >>> *\n" );
		}
	}	

	for( i=0; i < rfield.nupper; i++ )
	{
		rfield.tslop[rfield.nShape][i] = (realnum)pow((realnum)10.f,rfield.tslop[rfield.nShape][i]);
		if( rfield.tslop[rfield.nShape][i] < 1e-37 )
			rfield.tslop[rfield.nShape][i] = 0.;
	}

	if( false )
	{
		FILE *ioBUG = fopen( "interpolated.txt", "w" );
		for( k=0; k < rfield.nupper; ++k )
			fprintf( ioBUG, "%e %e\n", rfield.tNu[rfield.nShape][k].Ryd(), rfield.tslop[rfield.nShape][k] );
		fclose( ioBUG );
	}

	if( strcmp( grid->names[0], "Teff" ) == 0 )
	{
		if( ! lgValidModel( rfield.tNu[rfield.nShape], rfield.tslop[rfield.nShape], val[0], 0.10 ) )
			TotalInsanity();
	}

	/* set limits for optimizer */
	SetLimits( grid, val[0], indlo, indhi, NULL, NULL, Tlow, Thigh );

	FREE_CHECK( aval );
	FREE_CHECK( index );
	FREE_CHECK( indhi );
	FREE_CHECK( indlo );
	return;
}

STATIC void FreeGrid(stellar_grid *grid)
{
	long i;

	DEBUG_ENTRY( "FreeGrid()" );

	/* this was opened/allocated in InitGrid and subsidiaries,
	 * this should become a destructor in C++ */
	fclose( grid->ioIN );
	FREE_CHECK( grid->telg );
	for( i = 0; i < grid->ndim; i++ )
		FREE_CHECK( grid->val[i] );
	FREE_CHECK( grid->val );
	FREE_CHECK( grid->nval );
	FREE_SAFE( grid->jlo );
	FREE_SAFE( grid->jhi );
	FREE_SAFE( grid->trackLen )
	FREE_SAFE( grid->jval );
	return;
}

STATIC void InterpolateModel(const stellar_grid *grid,
			     const double val[],
			     double aval[],
			     const long indlo[],
			     const long indhi[],
			     long index[],
			     long nd,
			     vector<realnum>& flux1,
			     IntStage stage)
{
	bool lgDryRun;
	long i, ind, j;

	DEBUG_ENTRY( "InterpolateModel()" );

	--nd;

	lgDryRun = ( flux1.size() == 0 );

	if( nd < 0 )
	{
		long n = JIndex(grid,index);
		if( stage == IS_FIRST )
			ind = ( grid->jlo[n] >= 0 ) ? grid->jlo[n] : grid->jhi[n];
		else if( stage == IS_SECOND ) 
			ind = ( grid->jhi[n] >= 0 ) ? grid->jhi[n] : grid->jlo[n];
		else if( grid->ndim == 1 )
			/* in this case grid->jlo[n] and grid->jhi[n] should be identical */
			ind = grid->jlo[n];
		else
			TotalInsanity();

		if( ind < 0 )
		{
			fprintf( ioQQQ, " The requested interpolation could not be completed, sorry.\n" );
			fprintf( ioQQQ, " No suitable match was found for a model with" );
			for( i=0; i < grid->ndim; i++ )
				fprintf( ioQQQ, " %s=%.6g ", grid->names[i], grid->val[i][index[i]] );
			fprintf( ioQQQ, "\n" );
			cdEXIT(EXIT_FAILURE);
		}

		for( i=0; i < grid->npar; i++ )
			aval[i] = grid->telg[ind].par[i];

		if( !lgDryRun )
		{
			for( i=0; i < grid->ndim && called.lgTalk; i++ )
			{
				if( !fp_equal(grid->val[i][index[i]],aval[i],10) )
				{
					fprintf( ioQQQ, " No exact match was found for a model with" );
					for( j=0; j < grid->ndim; j++ )
						fprintf( ioQQQ, " %s=%.6g ", grid->names[j], grid->val[j][index[j]] );
					fprintf( ioQQQ, "- using the following model instead:\n" );
					break;
				}
			}

			GetModel( grid, ind, flux1, lgVERBOSE, lgTAKELOG );
		}
	}
	else
	{
		vector<realnum> flux2(rfield.nupper);
		double *aval2;

#		if !defined NDEBUG
		const realnum SECURE = 10.f*FLT_EPSILON;
#		endif

		aval2 = (double*)MALLOC((size_t)(grid->npar*sizeof(double)) );

		/* Interpolation is carried out first in the parameter with nd == 0 (usually
		 * Teff), then the parameter with nd == 1 (usually log(g)), etc. One or two
		 * atmosphere models are read depending on whether the parameter was matched
		 * exactly or not. If needed, logarithmic interpolation is done.
		 */

		if( nd == 1 )
			stage = IS_FIRST;

		index[nd] = indlo[nd];
		InterpolateModel( grid, val, aval, indlo, indhi, index, nd, flux1, stage );

		if( nd == 1 )
			stage = IS_SECOND;

		index[nd] = indhi[nd];
		vector<realnum> empty;
		InterpolateModel( grid, val, aval2, indlo, indhi, index, nd, empty, stage );

		if( !fp_equal(aval2[nd],aval[nd],10) )
		{
			double fr1, fr2, fc1 = 0., fc2 = 0.;

			if( !lgDryRun )
				InterpolateModel( grid, val, aval2, indlo, indhi, index, nd, flux2, stage );

			fr1 = (aval2[nd]-val[nd])/(aval2[nd]-aval[nd]);
			/* when interpolating in log(g) it can happen that fr1 is outside the range 0 .. 1
			 * this can be the case when the requested log(g) was not present in the grid
			 * and it had to be approximated by another model. In this case do not extrapolate */
			if( nd == 1 )
				fr1 = MIN2( MAX2( fr1, 0. ), 1. );
			fr2 = 1. - fr1;

			ASSERT( 0.-SECURE <= fr1 && fr1 <= 1.+SECURE );

			if( !lgDryRun )
			{
#				if DEBUGPRT
				fprintf( ioQQQ, "interpolation nd=%ld fr1=%g\n", nd, fr1 );
#				endif

				/* special treatment for high-temperature Rauch models */
				if( nd == 0 && strcmp( grid->names[nd], "Teff" ) == 0 )
				{
					/* The following is an approximate scaling to use for the range of 
					 * temperatures above 200000 K in the H-Ca Rauch models where the
					 * temperature steps are large and thus the interpolations are over
					 * large ranges.  For the lower temperatures I assume that there is
					 * no need for this.
					 *
					 * It should be remembered that this interpolation is not exact, and 
					 * the possible error at high temperatures might be large enough to 
					 * matter. (Kevin Volk)
					 */
					fc1 = ( val[nd] > 200000. ) ? log10(val[nd]/grid->val[nd][indlo[nd]])*4. : 0.;
					fc2 = ( val[nd] > 200000. ) ? log10(val[nd]/grid->val[nd][indhi[nd]])*4. : 0.;
				}

				for( i=0; i < rfield.nupper; ++i )
					flux1[i] = (realnum)(fr1*(flux1[i]+fc1) + fr2*(flux2[i]+fc2));
			}

			for( i=0; i < grid->npar; i++ )
				aval[i] = fr1*aval[i] + fr2*aval2[i];
		}

		FREE_CHECK( aval2 );
	}
	return;
}

STATIC void InterpolateModelCoStar(const stellar_grid *grid,
				   const double val[],
				   double aval[],
				   const long indlo[],
				   const long indhi[],
				   long index[],
				   long nd,
				   long off,
				   vector<realnum>& flux1)
{
	long i, ind;

	DEBUG_ENTRY( "InterpolateModelCoStar()" );

	if( nd == 2 )
	{
		ind = ( index[1] == 0 ) ? indlo[index[0]] : indhi[index[0]];

		GetModel( grid, ind, flux1, lgVERBOSE, lgTAKELOG );

		for( i=0; i < grid->npar; i++ )
			aval[i] = grid->telg[ind].par[i];
	}
	else
	{
		bool lgSkip;
#		if !defined NDEBUG
		const realnum SECURE = 10.f*FLT_EPSILON;
#		endif

		/* Interpolation is carried out first along evolutionary tracks, then
		 * in between evolutionary tracks. Between 1 and 4 atmosphere models are read
		 * depending on whether the parameter/track was matched exactly or not.
		 */

		index[nd] = 0;
		InterpolateModelCoStar( grid, val, aval, indlo, indhi, index, nd+1, off, flux1 );

		lgSkip = ( nd == 1 ) ?  ( indhi[index[0]] == indlo[index[0]] ) :
			( indlo[0] == indlo[1] && indhi[0] == indhi[1] );

		if( ! lgSkip )
		{
			vector<realnum> flux2(rfield.nupper);
			double fr1, fr2, *aval2;

			aval2 = (double*)MALLOC((size_t)(grid->npar*sizeof(double)) );

			index[nd] = 1;
			InterpolateModelCoStar( grid, val, aval2, indlo, indhi, index, nd+1, off, flux2 );

			fr1 = (aval2[nd+off]-val[nd])/(aval2[nd+off]-aval[nd+off]);
			fr2 = 1. - fr1;

#			if DEBUGPRT
			fprintf( ioQQQ, "interpolation nd=%ld fr1=%g\n", nd, fr1 );
#			endif

			ASSERT( 0.-SECURE <= fr1 && fr1 <= 1.+SECURE );

			for( i=0; i < rfield.nupper; ++i )
				flux1[i] = (realnum)(fr1*flux1[i] + fr2*flux2[i]);

			for( i=0; i < grid->npar; i++ )
				aval[i] = fr1*aval[i] + fr2*aval2[i];

			FREE_CHECK( aval2 );
		}
	}
	return;
}

STATIC void GetBins(const stellar_grid *grid,
		    vector<Energy>& ener)
{
	DEBUG_ENTRY( "GetBins()" );

	/* make sure ident is exactly 12 characters long, otherwise output won't fit */
	ASSERT( strlen(grid->ident) == 12 );
	
	ASSERT( grid->nBlocksize == rfield.nupper*sizeof(realnum) );

	/* skip over ind stars */
	/* >>chng 01 oct 18, add nOffset */
	if( fseek( grid->ioIN, (long)(grid->nOffset), SEEK_SET ) != 0 )
	{
		fprintf( ioQQQ, " Error finding atmosphere frequency bins\n");
		cdEXIT(EXIT_FAILURE);
	}

	vector<realnum> data(rfield.nupper);
	if( fread( get_ptr(data), 1, grid->nBlocksize, grid->ioIN ) != grid->nBlocksize )
	{
		fprintf( ioQQQ, " Error reading atmosphere frequency bins\n" );
		cdEXIT(EXIT_FAILURE);
	}

	for( long i=0; i < rfield.nupper; ++i )
		ener[i].set(data[i]);
	
	return;
}

STATIC void GetModel(const stellar_grid *grid,
		     long ind,
		     vector<realnum>& flux,
		     bool lgTalk,
		     bool lgTakeLog)
{
	long i;

	DEBUG_ENTRY( "GetModel()" );

	/* add 1 to account for frequency grid that is stored in front of all the atmospheres */
	ind++;

	/* make sure ident is exactly 12 characters long, otherwise output won't fit */
	ASSERT( strlen(grid->ident) == 12 );
	/* ind == 0 is the frequency grid, ind == 1 .. nmods are the atmosphere models */
	ASSERT( ind >= 0 && ind <= grid->nmods );

	/* skip over ind stars */
	/* >>chng 01 oct 18, add nOffset */
	if( fseek( grid->ioIN, (long)(ind*grid->nBlocksize+grid->nOffset), SEEK_SET ) != 0 )
	{
		fprintf( ioQQQ, " Error seeking atmosphere %ld\n", ind );
		cdEXIT(EXIT_FAILURE);
	}

	if( fread( get_ptr(flux), 1, grid->nBlocksize, grid->ioIN ) != grid->nBlocksize )
	{
		fprintf( ioQQQ, " Error trying to read atmosphere %ld\n", ind );
		cdEXIT(EXIT_FAILURE);
	}

	/* print the parameters of the atmosphere model */
	if( called.lgTalk && lgTalk )
	{
		/* ind-1 below since telg doesn't have the entry for the frequency grid */
		if( grid->npar == 1 )
			fprintf( ioQQQ, 
				 "                       * c<< %s model%5ld read.  "
				 "  %6s = %13.2f                 >>> *\n", 
				 grid->ident, ind, grid->names[0], grid->telg[ind-1].par[0] );
		else if( grid->npar == 2 )
			fprintf( ioQQQ, 
				 "                       * c<< %s model%5ld read.  "
				 "  %6s = %10.2f %6s = %8.5f  >>> *\n", 
				 grid->ident, ind, grid->names[0], grid->telg[ind-1].par[0],
				 grid->names[1], grid->telg[ind-1].par[1] );
		else if( grid->npar == 3 )
			fprintf( ioQQQ, 
				 "                       * c<< %s model%5ld read. "
				 " %6s=%7.0f %6s=%5.2f %6s=%5.2f >>> *\n", 
				 grid->ident, ind, grid->names[0], grid->telg[ind-1].par[0],
				 grid->names[1], grid->telg[ind-1].par[1],
				 grid->names[2], grid->telg[ind-1].par[2] );
		else if( grid->npar >= 4 )
		{
			fprintf( ioQQQ, 
				 "                       * c< %s mdl%4ld"
				 " %4s=%5.0f %6s=%4.2f %6s=%5.2f %6s=",
				 grid->ident, ind, grid->names[0], grid->telg[ind-1].par[0],
				 grid->names[1], grid->telg[ind-1].par[1],
				 grid->names[2], grid->telg[ind-1].par[2], grid->names[3] );
			fprintf( ioQQQ, PrintEfmt( "%9.2e", grid->telg[ind-1].par[3] ) );
			fprintf( ioQQQ, " >> *\n" ); 
		}
	}	

	if( lgTakeLog )
	{
		/* convert to logs since we will interpolate in log flux */
		for( i=0; i < rfield.nupper; ++i )
			flux[i] = (realnum)log10( MAX2( 1e-37, (double)flux[i] ) );
	}
	return;
}

STATIC void SetLimits(const stellar_grid *grid,
		      double val,
		      const long indlo[],
		      const long indhi[],
		      const long useTr[],
		      const realnum ValTr[],
		      double *loLim,
		      double *hiLim)
{
	DEBUG_ENTRY( "SetLimits()" );

	if( optimize.lgVarOn )
	{
		int ptr0,ptr1;
		long index[MDIM], j;
		const double SECURE = (1. + 20.*(double)FLT_EPSILON);

		*loLim = +DBL_MAX;
		*hiLim = -DBL_MAX;

		switch( grid->imode )
		{
		case IM_RECT_GRID:
			*loLim = -DBL_MAX;
			*hiLim = +DBL_MAX;
			SetLimitsSub( grid, val, indlo, indhi, index, grid->ndim, loLim, hiLim );
			break;
		case IM_COSTAR_TEFF_MODID:
		case IM_COSTAR_TEFF_LOGG:
		case IM_COSTAR_MZAMS_AGE:
			for( j=0; j < grid->nTracks; j++ )
			{
				if( ValTr[j] != -FLT_MAX )
				{
					/* M_ZAMS is already logarithm, Teff is linear */
					double temp = ( grid->imode == IM_COSTAR_MZAMS_AGE ) ?
						pow(10.,(double)ValTr[j]) : ValTr[j];
					*loLim = MIN2(*loLim,temp);
					*hiLim = MAX2(*hiLim,temp);
				}
			}
			break;
		case IM_COSTAR_AGE_MZAMS:
			index[0] = 0;
			index[1] = useTr[0];
			ptr0 = grid->jval[JIndex(grid,index)];
			index[1] = useTr[1];
			ptr1 = grid->jval[JIndex(grid,index)];
			*loLim = MAX2(grid->telg[ptr0].par[3],grid->telg[ptr1].par[3]);
#			if DEBUGPRT
			printf( "set limit 0: (models %d, %d) %f %f\n",
			       ptr0+1, ptr1+1, grid->telg[ptr0].par[3], grid->telg[ptr1].par[3] );
#			endif
			index[0] = grid->trackLen[useTr[0]]-1;
			index[1] = useTr[0];
			ptr0 = grid->jval[JIndex(grid,index)];
			index[0] = grid->trackLen[useTr[1]]-1;
			index[1] = useTr[1];
			ptr1 = grid->jval[JIndex(grid,index)];
			*hiLim = MIN2(grid->telg[ptr0].par[3],grid->telg[ptr1].par[3]);
#			if DEBUGPRT
			printf( "set limit 1: (models %d, %d) %f %f\n",
			       ptr0+1, ptr1+1, grid->telg[ptr0].par[3], grid->telg[ptr1].par[3] );
#			endif
			break;
		default:
			fprintf( ioQQQ, " SetLimits called with insane value for imode: %d.\n", grid->imode );
			cdEXIT(EXIT_FAILURE);
		}

		ASSERT( fabs(*loLim) < DBL_MAX && fabs(*hiLim) < DBL_MAX );

		/* check sanity of optimization limits */
		if( *hiLim <= *loLim )
		{
			fprintf( ioQQQ, " no room to optimize: lower limit %.4f, upper limit %.4f.\n",
				 *loLim,*hiLim );
			cdEXIT(EXIT_FAILURE);
		}

		/* make a bit of room for round-off errors */
		*loLim *= SECURE;
		*hiLim /= SECURE;

#		if DEBUGPRT
		printf("set limits: %g %g\n",*loLim,*hiLim);
#		endif
	}
	else
	{
		*loLim = 0.;
		*hiLim = 0.;
	}
	return;
}

STATIC void SetLimitsSub(const stellar_grid *grid,
			 double val,
			 const long indlo[],
			 const long indhi[],
			 long index[],
			 long nd,
			 double *loLim,
			 double *hiLim)
{
	long n;

	DEBUG_ENTRY( "SetLimitsSub()" );

	--nd;

	if( nd < 1 )
	{
		double loLoc = +DBL_MAX;
		double hiLoc = -DBL_MAX;

		for( index[0]=0; index[0] < grid->nval[0]; ++index[0] )
		{
			/* grid->val[0][i] is the array of Par0 values (Teff/Age/...) in the
			 * grid, which it is sorted in strict monotonically increasing order.
			 * This routine searches for the largest range [loLoc,hiLoc] in Par0
			 * such that loLoc <= val <= hiLoc, and at least one model exists for
			 * each Par0 value in this range. This assures that interpolation is
			 * safe and the optimizer will not trip... */
			n = JIndex(grid,index);
			if( grid->jhi[n] < 0 )
			{
				/* there are no models with this value of Par0 */
				/* this value of Par0 should be outside of allowed range */
				if( grid->val[0][index[0]] < val )
					loLoc = DBL_MAX;
				/* this is beyond the legal range, so terminate the search */
				if( grid->val[0][index[0]] > val )
					break;
			}
			else
			{
				/* there are models with this value of Par0 */
				/* update range to include this value of Par0 */
				if( grid->val[0][index[0]] <= val )
				{
					/* remember lowest legal value of loLoc
					 * -> only update if previous value was illegal */
					if( loLoc == DBL_MAX )
						loLoc = grid->val[0][index[0]];
				}
				if( grid->val[0][index[0]] >= val )
				{
					/* remember highest legal value of hiLoc
					 * -> always update */
					hiLoc = grid->val[0][index[0]];
				}
			}
		}

		ASSERT( fabs(loLoc) < DBL_MAX && fabs(hiLoc) < DBL_MAX && loLoc <= hiLoc );

		*loLim = MAX2(*loLim,loLoc);
		*hiLim = MIN2(*hiLim,hiLoc);
	}
	else
	{
		index[nd] = indlo[nd];
		SetLimitsSub( grid, val, indlo, indhi, index, nd, loLim, hiLim );

		if( indhi[nd] != indlo[nd] )
		{
			index[nd] = indhi[nd];
			SetLimitsSub( grid, val, indlo, indhi, index, nd, loLim, hiLim );
		}
	}
	return;
}

STATIC void InitIndexArrays(stellar_grid *grid,
			    bool lgList)
{
	long i, *index, j, jsize, nd;
	double *val;

	DEBUG_ENTRY( "InitIndexArrays()" );

	ASSERT( grid->telg != NULL );
	ASSERT( grid->nmods > 0 );

	jsize = 1;

	/* this loop creates a list of all unique model parameter values in increasing order */
	for( nd=0; nd < grid->ndim; nd++ )
	{
		double pval = grid->telg[0].par[nd];
		grid->val[nd][0] = pval;
		grid->nval[nd] = 1;

		for( i=1; i < grid->nmods; i++ )
		{
			bool lgOutOfRange;
			long i1, i2;

			pval = grid->telg[i].par[nd];
			FindIndex( grid->val[nd], grid->nval[nd], pval, &i1, &i2, &lgOutOfRange );
			/* if i1 < i2, the new parameter value was not present yet and
			 * it needs to be inserted in between i1 and i2 --> first move
			 * all entries from i2 to grid->nval[nd]-1 one slot upward and
			 * then insert the new value at i2; this also works correctly
			 * if lgOutOfRange is set, hence no special check is needed */ 
			if( i1 < i2 )
			{
				/* val[nd] has grid->nmods entries, so cannot overflow */
				for( j = grid->nval[nd]-1; j >= i2; j-- )
					grid->val[nd][j+1] = grid->val[nd][j];
				grid->val[nd][i2] = pval;
				grid->nval[nd]++;
			}
		}

		jsize *= grid->nval[nd];

#		if DEBUGPRT
		printf( "%s[%ld]:", grid->names[nd], grid->nval[nd] );
		for( i=0; i < grid->nval[nd]; i++ )
			printf( " %g", grid->val[nd][i] );
		printf( "\n" );
#		endif
	}

	index = (long *)MALLOC( sizeof(long)*grid->ndim );
	val = (double *)MALLOC( sizeof(double)*grid->ndim );

	/* this memory will be freed in the calling function */
	grid->jlo = (long *)MALLOC( sizeof(long)*jsize );
	grid->jhi = (long *)MALLOC( sizeof(long)*jsize );

	/* set up square array of model indices; this will be used to
	 * choose the correct models for the interpolation process */
	FillJ( grid, index, val, grid->ndim, lgList );

	FREE_CHECK( val );
	FREE_CHECK( index );

	if( lgList )
		cdEXIT(EXIT_SUCCESS);
	return;
}

STATIC void FillJ(const stellar_grid *grid,
		  long index[], /* index[grid->ndim] */
		  double val[], /* val[grid->ndim] */
		  long nd,
		  bool lgList)
{
	DEBUG_ENTRY( "FillJ()" );

	--nd;

	if( nd < 0 )
	{
		long n = JIndex(grid,index);
		SearchModel( grid->telg, grid->lgIsTeffLoggGrid, grid->nmods, val, grid->ndim,
			     &grid->jlo[n], &grid->jhi[n] );
	}
	else
	{
		for( index[nd]=0; index[nd] < grid->nval[nd]; index[nd]++ )
		{
			val[nd] = grid->val[nd][index[nd]];
			FillJ( grid, index, val, nd, lgList );
		}
	}

	if( lgList && nd == MIN2(grid->ndim-1,1) )
	{
		fprintf( ioQQQ, "\n" );
		if( grid->ndim > 2 )
		{
			fprintf( ioQQQ, "subgrid for" );
			for( long n = nd+1; n < grid->ndim; n++ )
				fprintf( ioQQQ, " %s=%g", grid->names[n], val[n] );
			fprintf( ioQQQ, ":\n\n" );
		}
		if( grid->ndim > 1 )
		{
			fprintf( ioQQQ, "%6.6s\\%6.6s |", grid->names[0], grid->names[1] );
			for( long n = 0; n < grid->nval[1]; n++ )
				fprintf( ioQQQ, " %9.3g", grid->val[1][n] );
			fprintf( ioQQQ, "\n" );
			fprintf( ioQQQ, "--------------|" );
			for( long n = 0; n < grid->nval[1]; n++ )
				fprintf( ioQQQ, "----------" );
		}
		else
		{
			fprintf( ioQQQ, "%13.13s |\n", grid->names[0] );
			fprintf( ioQQQ, "--------------|----------" );
		}
		fprintf( ioQQQ, "\n" );
		for( index[0]=0; index[0] < grid->nval[0]; index[0]++ )
		{
			fprintf( ioQQQ, "%13.7g |", grid->val[0][index[0]] );
			if( grid->ndim > 1 )
			{
				for( index[1]=0; index[1] < grid->nval[1]; index[1]++ )
					if( grid->jlo[JIndex(grid,index)] == grid->jhi[JIndex(grid,index)] &&
					    grid->jlo[JIndex(grid,index)] >= 0 )
						fprintf( ioQQQ, " %9ld", grid->jlo[JIndex(grid,index)]+1 );
					else
						fprintf( ioQQQ, "        --" );
			}
			else
			{
				fprintf( ioQQQ, " %9ld", grid->jlo[JIndex(grid,index)]+1 );
			}
			fprintf( ioQQQ, "\n" );
		}
		fprintf( ioQQQ, "\n" );
	}
	return;
}

STATIC long JIndex(const stellar_grid *grid,
		   const long index[]) /* index[grid->ndim] */
{
	long i, ind, mul;

	DEBUG_ENTRY( "JIndex()" );

	ind = 0;
	mul = 1;
	for( i=0; i < grid->ndim; i++ )
	{
		ind += index[i]*mul;
		mul *= grid->nval[i];
	}
	return ind;
}

STATIC void SearchModel(const mpp telg[], /* telg[nmods] */
			bool lgIsTeffLoggGrid,
			long nmods,
			const double val[], /* val[ndim] */
			long ndim,
			long *index_low,
			long *index_high)
{
	long i, nd;
	double alogg_low = -DBL_MAX, alogg_high = DBL_MAX;

	DEBUG_ENTRY( "SearchModel()" );

	/* given values for the model parameters, this routine searches for the atmosphere
	 * model that is the best match. If all parameters can be matched simultaneously the
	 * choice is obvious, but this cannot always be achieved (typically for high Teff, the
	 * low log(g) models will be missing). If lgIsTeffLoggGrid is true, the rule is that
	 * all parameters except log(g) must always be matched (such a model is not always
	 * guaranteed to exist). If all requested parameters can be matched exactly, both
	 * index_low and index_high will point to that model. If all parameters except log(g)
	 * can be matched exactly, it will return the model with the lowest log(g) value larger
	 * than the requested value in index_high, and the model with the highest log(g) value
	 * lower than the requested value in index_low. If either requirement cannot be
	 * fulfilled, -2 will be returned. When lgIsTeffLoggGrid is false, all parameters must
	 * be matched and both index_low and index_high will point to that model. If no such
	 * model can be found, -2 will be returned. */

	*index_low = *index_high = -2;
	for( i=0; i < nmods; i++ )
	{
		bool lgNext = false;
		/* ignore models with different parameters */
		for( nd=0; nd < ndim; nd++ )
		{
			if( nd != 1 && !fp_equal(telg[i].par[nd],val[nd],10) )
			{
				lgNext = true;
				break;
			}
		}
		if( lgNext )
			continue;

		/* an exact match is found */
		if( ndim == 1 || fp_equal(telg[i].par[1],val[1],10) )
		{
			*index_low = i;
			*index_high = i;
			return;
		}
		if( lgIsTeffLoggGrid )
		{
			/* keep a record of the highest log(g) model smaller than alogg */
			if( telg[i].par[1] < val[1] && telg[i].par[1] > alogg_low )
			{
				*index_low = i;
				alogg_low = telg[i].par[1];
			}
			/* also keep a record of the lowest log(g) model greater than alogg */
			if( telg[i].par[1] > val[1] && telg[i].par[1] < alogg_high )
			{
				*index_high = i;
				alogg_high = telg[i].par[1];
			}
		}
	}
	return;
}

STATIC void FindIndex(const double xval[], /* xval[NVAL] */
		      long NVAL,
		      double x,
		      long *ind1,
		      long *ind2,
		      bool *lgInvalid)
{
	bool lgOutLo, lgOutHi;
	long i;

	DEBUG_ENTRY( "FindIndex()" );

	/* this routine searches for indices ind1, ind2 such that
	 *   xval[ind1] < x < xval[ind2]
	 * if x is equal to one of the values in xval, then
	 *   ind1 == ind2  and  xval[ind1] == x
	 *
	 * if x is outside the range xval[0] ... xval[NVAL-1]
	 * then lgInvalid will be set to true
	 *
	 * NB NB -- this routine implicitly assumes that xval is
	 *          strictly monotonically increasing!
	 */

	ASSERT( NVAL > 0 );

	/* is x outside of range xval[0] ... xval[NVAL-1]? */
	lgOutLo = ( x-xval[0] < -10.*DBL_EPSILON*fabs(xval[0]) );
	lgOutHi = ( x-xval[NVAL-1] > 10.*DBL_EPSILON*fabs(xval[NVAL-1]) );

	if( lgOutLo || lgOutHi )
	{
		/* pretend there are two fictitious array elements
		 *   xval[-1] = -Inf  and  xval[NVAL] = +Inf,
		 * and return ind1 and ind2 accordingly. This behavior
		 * is needed for InitIndexArrays() to work correctly */
		*ind1 = lgOutLo ? -1 : NVAL-1;
		*ind2 = lgOutLo ?  0 : NVAL;
		*lgInvalid = true;
		return;
	}

	*lgInvalid = false;

	/* there are more efficient ways of doing this, e.g. a binary search.
	 * However, the xval arrays typically only have 1 or 2 dozen elements,
	 * so the overhead is negligible and the clarity of this code is preferred */

	/* first look for an "exact" match */
	for( i=0; i < NVAL; i++ )
	{
		if( fp_equal(xval[i],x,10) )
		{
			*ind1 = i;
			*ind2 = i;
			return;
		}
	}

	/* no match was found -> bracket the x value */
	for( i=0; i < NVAL-1; i++ )
	{
		if( xval[i] < x && x < xval[i+1] )
		{
			*ind1 = i;
			*ind2 = i+1;
			return;
		}
	}

	/* this should never be reached ! */
	fprintf( ioQQQ, " insanity in FindIndex\n" );
	ShowMe();
	cdEXIT(EXIT_FAILURE);
}

STATIC bool lgFileReadable(const char *chFnam, process_counter& pc, access_scheme scheme)
{
	DEBUG_ENTRY( "lgFileReadable()" );

	FILE *ioIN;

	ioIN = open_data( chFnam, "r", scheme );
	if( ioIN != NULL )
	{
		fclose( ioIN );
		++pc.nFound;
		return true;
	}
	else
	{
		return false;
	}
}

/*ValidateGrid: check each model in the grid to see if it has the correct Teff */
STATIC void ValidateGrid(const stellar_grid *grid,
			 double toler)
{
	long i, k, nd;
	vector<Energy> anu(rfield.nupper);
	vector<realnum> flux(rfield.nupper);

	DEBUG_ENTRY( "ValidateGrid()" );

	if( strcmp( grid->names[0], "Teff" ) != 0 )
	{
		return;
	}

	GetBins( grid, anu );

	for( i=0; i < grid->nmods; i++ ) 
	{
		fprintf( ioQQQ, "testing model %ld ", i+1 );
		for( nd=0; nd < grid->npar; nd++ )
			fprintf( ioQQQ, " %s %g", grid->names[nd], grid->telg[i].par[nd] );

		GetModel( grid, i, flux, lgSILENT, lgLINEAR );

		if( lgValidModel( anu, flux, grid->telg[i].par[0], toler ) )
			fprintf( ioQQQ, "   OK\n" );

		if( false )
		{
			FILE *ioBUG = fopen( "atmosphere_dump.txt", ( i == 0 ) ? "w" : "a" );

			fprintf( ioBUG, "######## MODEL %ld", i+1 );
			for( nd=0; nd < grid->npar; nd++ )
				fprintf( ioBUG, " %s %g", grid->names[nd], grid->telg[i].par[nd] );
			fprintf( ioBUG, "####################\n" );

			for( k=0; k < rfield.nupper; ++k )
				fprintf( ioBUG, "%e %e\n", anu[k].Ryd(), flux[k] );

			fclose( ioBUG );
		}
	}
	return;
}

STATIC bool lgValidModel(const vector<Energy>& anu,
			 const vector<realnum>& flux,
			 double Teff,
			 double toler)
{
	bool lgPassed = true;
	long k;
	double chk, lumi;

	DEBUG_ENTRY( "lgValidModel()" );

	ASSERT( Teff > 0. );

	lumi = 0.;
	/* rebinned models are in cgs F_nu units */
	for( k=1; k < rfield.nupper; k++ )
		lumi += (anu[k].Ryd() - anu[k-1].Ryd())*(flux[k] + flux[k-1])/2.;

	/* now convert luminosity to effective temperature */
	chk = pow(lumi*FR1RYD/STEFAN_BOLTZ,0.25);
	/* the allowed tolerance is set by the caller in toler */
	if( fabs(Teff - chk) > toler*Teff ) {
		fprintf( ioQQQ, "\n*** WARNING, Teff discrepancy for this model, expected Teff %.2f, ", Teff);
		fprintf( ioQQQ, "integration yielded Teff %.2f, delta %.2f%%\n", chk, (chk/Teff-1.)*100. );
		lgPassed = false;
	}
	return lgPassed;
}

/*RebinAtmosphere: generic routine for rebinning atmospheres onto Cloudy grid */
STATIC void RebinAtmosphere(long nCont, /* the number of points in the incident continuum*/
			    const realnum StarEner[], /* StarEner[nCont], the freq grid for the model, in Ryd*/
			    const realnum StarFlux[], /* StarFlux[nCont], the original model flux */
			    realnum CloudyFlux[],     /* CloudyFlux[NC_ELL], the model flux on the cloudy grid */
			    long nEdge, /* the number of bound-free continuum edges in AbsorbEdge */
			    const realnum AbsorbEdge[]) /* AbsorbEdge[nEdge], energies of the edges */
{
	bool lgDone;
	long int ind,
	  j,
	  k;
	/* >>chng 00 jun 02, demoted next two to realnum, PvH */
	realnum BinHigh, 
	  BinLow,
	  BinMid,
	  BinNext,
	  *EdgeLow=NULL,
	  *EdgeHigh=NULL,
	  *StarPower;

	DEBUG_ENTRY( "RebinAtmosphere()" );

	if( nEdge > 0 )
	{
		EdgeLow = (realnum*)MALLOC( sizeof(realnum)*(unsigned)nEdge );
		EdgeHigh = (realnum*)MALLOC( sizeof(realnum)*(unsigned)nEdge );
	}

	/* this loop should be before the next loop, otherwise models with a
	 * very strong He II edge (i.e. no flux beyond that edge) will fail */
	for( j=0; j < nEdge; j++ )
	{
		ind = RebinFind(StarEner,nCont,AbsorbEdge[j]);

		/* sanity check */
		ASSERT( ind >= 0 && ind+1 < nCont );

		EdgeLow[j] = StarEner[ind];
		EdgeHigh[j] = StarEner[ind+1];
	}

	/* cut off that part of the Wien tail that evaluated to zero */
	/* >> chng 05 nov 22, inverted loop, slightly faster PvH */
	/*for( j=nCont-1; j >= 0; j-- )*/
	for( j=0; j < nCont; j++ )
	{
		if( StarFlux[j] == 0.f )
		{
			nCont = j;
			break;
		}
	}
	ASSERT( nCont > 0 );

	StarPower = (realnum *)MALLOC( sizeof(realnum)*(unsigned)(nCont-1) );

	for( j=0; j < nCont-1; j++ )
	{
		double ratio_x, ratio_y;

		/* >>chng 05 nov 22, add sanity check to prevent invalid fp operations */
		ASSERT( StarEner[j+1] > StarEner[j] );

		/* >>chng 06 aug 11, on some systems (e.g., macbook pro) y/x can get evaluated as y*(1/x);
		 * this causes overflows if x is a denormalized number, hence we force a cast to double, PvH */
		ratio_x = (double)StarEner[j+1]/(double)StarEner[j];
		ratio_y = (double)StarFlux[j+1]/(double)StarFlux[j];
		StarPower[j] = (realnum)(log(ratio_y)/log(ratio_x));
	}

	for( j=0; j < rfield.nupper; j++ )
	{
		/* >>chng 05 nov 22, modified BinLow, BinHigh, BinNext to make boundaries match exactly, PvH */
		/* BinLow is lower bound of this continuum cell */
		BinLow = ( j > 0 ) ?
			(realnum)sqrt(rfield.anu[j-1]*rfield.anu[j]) : (realnum)sqrt(POW3(rfield.anu[0])/rfield.anu[1]);

		/* BinHigh is upper bound of this continuum cell */
		BinHigh = ( j+1 < rfield.nupper ) ?
			(realnum)sqrt(rfield.anu[j]*rfield.anu[j+1]) : rfield.anu[rfield.nupper-1];

		/* BinNext is upper bound of next continuum cell */
		BinNext = ( j+2 < rfield.nupper ) ?
			(realnum)sqrt(rfield.anu[j+1]*rfield.anu[j+2]) : rfield.anu[rfield.nupper-1];

		lgDone = false;

		/* >>chng 00 aug 14, take special care not to interpolate over major edges,
		 * the region in between EdgeLow and EdgeHigh should be avoided,
		 * the spectrum is extremely steep there, leading to significant roundoff error, PvH */
		for( k=0; k < nEdge; k++ )
		{
			if( BinLow < EdgeLow[k] && BinNext > EdgeHigh[k] )
			{
				BinMid = 0.99999f*EdgeLow[k];
				CloudyFlux[j] = RebinSingleCell(BinLow,BinMid,StarEner,StarFlux,StarPower,nCont);
				j++;

				/* sanity check */
				ASSERT( j < rfield.nupper );

				BinMid = 1.00001f*EdgeHigh[k];
				CloudyFlux[j] = RebinSingleCell(BinMid,BinNext,StarEner,StarFlux,StarPower,nCont);
				lgDone = true;
				break;
			}
		}

		/* default case when we are not close to an edge */
		if( !lgDone )
		{
			CloudyFlux[j] = RebinSingleCell(BinLow,BinHigh,StarEner,StarFlux,StarPower,nCont);
		}
	}

	FREE_CHECK( StarPower );
	FREE_SAFE( EdgeHigh );
	FREE_SAFE( EdgeLow );
	return;
}

STATIC realnum RebinSingleCell(realnum BinLow,
			       realnum BinHigh,
			       const realnum StarEner[],  /* StarEner[nCont] */
			       const realnum StarFlux[],  /* StarFlux[nCont] */
			       const realnum StarPower[], /* StarPower[nCont-1] */
			       long nCont)
{
	long int i, 
	  ipHi, 
	  ipLo;
	double anu,
	  retval,
	  widflx;
	double sum,
	  v1, 
	  val, 
	  x1, 
	  x2;

	DEBUG_ENTRY( "RebinSingleCell()" );

	/* >>chng 05 nov 22, use geometric mean instead of arithmetic mean, PvH */
	anu = sqrt(BinLow*BinHigh);
	/* >>chng 05 nov 22, reduce widflx if cell sticks out above highest frequency in model, PvH */
	widflx = MIN2(BinHigh,StarEner[nCont-1])-BinLow;

	if( BinLow < StarEner[0] )
	{
		/* this is case where Cloudy's continuum is below stellar continuum,
		 * (at least for part of the cell), so we do Rayleigh Jeans extrapolation */
		retval = (realnum)(StarFlux[0]*pow(anu/StarEner[0],2.));
	}
	else if( BinLow > StarEner[nCont-1] )
	{
		/* case where cloudy continuum is entirely above highest stellar point */
		retval = 0.0e00;
	}
	else
	{
		/* now go through stellar continuum to find bins corresponding to
		 * this cloudy cell, stellar continuum defined through nCont cells */
		ipLo = RebinFind(StarEner,nCont,BinLow);
		ipHi = RebinFind(StarEner,nCont,BinHigh);
		/* sanity check */
		ASSERT( ipLo >= 0 && ipLo < nCont-1 && ipHi >= ipLo );

		if( ipHi == ipLo )
		{
			/* Do the case where the cloudy cell and its edges are between
			 * two adjacent stellar model points: do power-law interpolation  */
			retval = (realnum)(StarFlux[ipLo]*pow(anu/StarEner[ipLo],(double)StarPower[ipLo]));
		}
		else
		{
			/* Do the case where the cloudy cell and its edges span two or more
			 * stellar model cells:  add segments with power-law interpolation up to
			 * do the averaging.*/

			sum = 0.;

			/* ipHi points to stellar point at high end of cloudy continuum cell,
			 * if the Cloudy cell extends beyond the stellar grid, ipHi == nCont-1
			 * and the MIN2(ipHi,nCont-2) prevents access beyond allocated memory
			 * ipLo points to low end, above we asserted that 0 <= ipLo < nCont-1 */
			for( i=ipLo; i <= MIN2(ipHi,nCont-2); i++ )
			{
				double pp1 = StarPower[i] + 1.;

				if( i == ipLo )
				{
					x1 = BinLow;
					x2 = StarEner[i+1];
					v1 = StarFlux[i]*pow(x1/StarEner[i],(double)StarPower[i]);
					/*v2 = StarFlux[i+1];*/
				}

				else if( i == ipHi )
				{
					x2 = BinHigh;
					x1 = StarEner[i];
					/*v2 = StarFlux[i]*pow(x2/StarEner[i],StarPower[i]);*/
					v1 = StarFlux[i];
				}

				/*if( i > ipLo && i < ipHi )*/
				else
				{
					x1 = StarEner[i];
					x2 = StarEner[i+1];
					v1 = StarFlux[i];
					/*v2 = StarFlux[i+1];*/
				}

				if( fabs(pp1) < 0.001 )
				{
					val = x1*v1*log(x2/x1);
				}
				else
				{
					val = pow(x2/x1,pp1) - 1.;
					val = val*x1*v1/pp1;
				}
				sum += val;
			}

			retval = sum/widflx;
		}
	}
	return (realnum)retval;
}

STATIC long RebinFind(const realnum array[], /* array[nArr] */
		      long nArr,
		      realnum val)
{
	long i1,
	  i2,
	  i3,
	  ind = -2,
	  sgn;

	DEBUG_ENTRY( "RebinFind()" );

	/* sanity check */
	ASSERT( nArr > 1 );

	/* return ind(val) such that array[ind] <= val <= array[ind+1],
	 *
	 * NB NB: this routine assumes that array[] increases monotonically !
	 *
	 * the first two clauses indicate out-of-bounds conditions and
	 * guarantee that when val1 <= val2, also ind(val1) <= ind(val2) */

	if( val < array[0] )
	{
		ind = -1;
	}
	else if( val > array[nArr-1] )
	{
		ind = nArr-1;
	}
	else
	{
		/* do a binary search for ind */
		i1 = 0;
		i3 = nArr-1;
		while( i3-i1 > 1 )
		{
			i2 = (i1+i3)/2;
			sgn = sign3(val-array[i2]);

			switch(sgn)
			{
			case -1:
				i3 = i2;
				break;
			case 0:
				ind = i2;
				return( ind );
			case 1:
				i1 = i2;
				break;
			}
		}
		ind = i1;
	}

	/* sanity check */
	ASSERT( ind > -2 );
	return ind;
}
/*lint +e785 too few initializers */
/*lint +e801 use of go to depreciated */
