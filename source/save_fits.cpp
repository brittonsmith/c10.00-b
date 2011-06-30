/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "cddrive.h"
#include "optimize.h"
#include "grid.h"
#include "save.h"
#include "rfield.h"
#include "prt.h"
#include "input.h"
#include "version.h"
#include "physconst.h"

#define RECORDSIZE	2880
#define LINESIZE	80

#if defined(_BIG_ENDIAN) 
	/* the value of A will not be manipulated */
#	define HtoNL(A) (A)	
/*
#	define HtoNS(A) (A)
#	define NtoHS(A) (A)
#	define NtoHL(A) (A)
*/
#else
/* defined(_LITTLE_ENDIAN) */
/* the value of A will be byte swapped */
#	define HtoNL(A) ((((A) & 0xff000000) >> 24) | \
		(((A) & 0x00ff0000) >> 8) | \
		(((A) & 0x0000ff00) << 8) | \
		(((A) & 0x000000ff) << 24))
/*
#	define HtoNS(A) ((((A) & 0xff00) >> 8) | (((A) & 0x00ff) << 8))
#	define NtoHS HtoNS
#	define NtoHL HtoNL
*/
/*#else
error One of BIG_ENDIAN or LITTLE_ENDIAN must be #defined.*/
#endif

#define ByteSwap5(x) ByteSwap((unsigned char *) &x,sizeof(x))

#if !defined(_BIG_ENDIAN) 
STATIC void ByteSwap(unsigned char * b, int n)
{
	register int i = 0;
	register int j = n-1;
	while (i<j)
	{
		char temp = b[i];
		b[i] = b[j];
		b[j] = temp;
		/* std::swap(b[i], b[j]); */
		i++, j--;
	}
	return;
}
#endif

static FILE *ioFITS_OUTPUT;
static long bytesAdded = 0;
static long bitpix = 8;
static long pcount = 0;
static long gcount = 1;
static long maxParamValues = 0;
const char ModelUnits[2][17] = {"'dimensionless '", "'photons/cm^2/s'" };

STATIC void punchFITS_PrimaryHeader( bool lgAddModel );
STATIC void punchFITS_ParamHeader( /* long *numParamValues, */ long nintparm, long naddparm );
STATIC void punchFITS_ParamData( char **paramNames, long *paramMethods, realnum **paramRange,
								realnum **paramData, long nintparm, long naddparm,
								long *numParamValues );
STATIC void punchFITS_EnergyHeader( long numEnergies );
STATIC void punchFITS_EnergyData( const vector<realnum>& Energies, long EnergyOffset );
STATIC void punchFITS_SpectraHeader( bool lgAdditiveModel, long nintparm, long naddparm, long totNumModels, long numEnergies );
STATIC void punchFITS_SpectraData( realnum **interpParameters, multi_arr<realnum,3>& theSpectrum, int option,
								  long totNumModels, long numEnergies, long nintparm, long naddparm  );
STATIC void punchFITS_GenericHeader( long numEnergies );
STATIC void punchFITS_GenericData( long numEnergies, long ipLoEnergy, long ipHiEnergy );
STATIC void writeCloudyDetails( void );
STATIC long addComment( const char *CommentToAdd );
STATIC long addKeyword_txt( const char *theKeyword, const void *theValue, const char *theComment, long Str_Or_Log );
STATIC long addKeyword_num( const char *theKeyword, long theValue, const char *theComment);

void saveFITSfile( FILE* ioPUN, int option )
{

	long i;

	DEBUG_ENTRY( "saveFITSfile()" );

	if( !grid.lgGridDone && option != NUM_OUTPUT_TYPES )
	{
		/* don't do any of this until flag set true. */
		return;
	}

	ioFITS_OUTPUT = ioPUN;

#if	0
	{

		long i,j;
		FILE *asciiDump;

		if( (asciiDump = fopen(  "gridspectra.con","w" ) ) == NULL )
		{
			fprintf( ioQQQ, "saveFITSfile could not open save file for writing.\nSorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		for( i=0; i<grid.numEnergies-1; i++ )
		{
			fprintf( asciiDump, "%.12e\t", grid.Energies[i] );
			for( j=0; j<grid.totNumModels; j++ )
			{
				fprintf( asciiDump, "%.12e\t", grid.Spectra[4][j][i] );
			}
			fprintf( asciiDump, "\n" );
		}
		fclose( asciiDump );
	}
#endif

	ASSERT( option >= 0 );

	/* This is generic FITS option */
	if( option == NUM_OUTPUT_TYPES )
	{
		punchFITS_PrimaryHeader( false );
		punchFITS_GenericHeader( rfield.nflux - 1 );
		punchFITS_GenericData( rfield.nflux -1, 0, rfield.nflux -2 );
	}
	/* These are specially designed XSPEC outputs. */
	else if( option < NUM_OUTPUT_TYPES )
	{
		bool lgAdditiveModel;

		/* option 10 is exp(-tau). */
		if( option == 10 )
		{
			/* false says not an additive model */
			lgAdditiveModel = false;
		}
		else
		{
			lgAdditiveModel = true;			
		}

		punchFITS_PrimaryHeader( lgAdditiveModel );

		for( i=0; i<grid.nintparm+grid.naddparm; i++ )
		{
			maxParamValues = MAX2( maxParamValues, grid.numParamValues[i] );
		}

		ASSERT( maxParamValues >= 2 );

		punchFITS_ParamHeader( /* grid.numParamValues, */ grid.nintparm, grid.naddparm );
		punchFITS_ParamData( grid.paramNames, grid.paramMethods, grid.paramRange, grid.paramData,
			grid.nintparm, grid.naddparm, grid.numParamValues );
		punchFITS_EnergyHeader( grid.numEnergies );
		punchFITS_EnergyData( grid.Energies, grid.ipLoEnergy );
		punchFITS_SpectraHeader( lgAdditiveModel, grid.nintparm, grid.naddparm, grid.totNumModels, grid.numEnergies);
		punchFITS_SpectraData( grid.interpParameters, grid.Spectra, option,
			grid.totNumModels, grid.numEnergies, grid.nintparm, grid.naddparm );
	}
	else 
	{
		fprintf( ioQQQ, "PROBLEM - undefined option encountered in saveFITSfile. \n" );
		cdEXIT( EXIT_FAILURE );
	}
	return;
}

STATIC void punchFITS_PrimaryHeader( bool lgAddModel )
{
	const char *ModelName = "'CLOUDY'";

	DEBUG_ENTRY( "punchFITS_PrimaryHeader()" );

	bytesAdded = 0;

	fixit(); /* bitpix is wrong when realnum is double? */

	bytesAdded += addKeyword_txt( "SIMPLE"	, "T",					"file does conform to FITS standard", 1 );
	bytesAdded += addKeyword_num( "BITPIX"	, bitpix,				"number of bits per data pixel" );
	bytesAdded += addKeyword_num( "NAXIS"	, 0,					"number of data axes" );
	bytesAdded += addKeyword_txt( "EXTEND"	, "T",					"FITS dataset may contain extensions", 1 );
	bytesAdded += addKeyword_txt( "CONTENT" , "'MODEL   '",			"spectrum file contains time intervals and event", 0 );
	bytesAdded += addKeyword_txt( "MODLNAME", ModelName,			"Model name", 0 );
	bytesAdded += addKeyword_txt( "MODLUNIT", ModelUnits[lgAddModel],	"Model units", 0 );
	bytesAdded += addKeyword_txt( "REDSHIFT", "T",				"If true then redshift will be included as a par", 1 );
	if( lgAddModel == true )
	{
		bytesAdded += addKeyword_txt( "ADDMODEL", "T",				"If true then this is an additive table model", 1 );
	}
	else
	{
		bytesAdded += addKeyword_txt( "ADDMODEL", "F",				"If true then this is an additive table model", 1 );
	}

	/* bytes are added here as well */
	writeCloudyDetails();

	bytesAdded += addKeyword_txt( "HDUCLASS", "'OGIP    '",			"Format conforms to OGIP/GSFC conventions", 0 );
	bytesAdded += addKeyword_txt( "HDUCLAS1", "'XSPEC TABLE MODEL'","Extension contains an image", 0 );
	bytesAdded += addKeyword_txt( "HDUVERS"	, "'1.0.0   '",			"Version of format (OGIP memo OGIP-92-001)", 0 );
		/* After everything else */
	bytesAdded += fprintf(ioFITS_OUTPUT, "%-80s", "END" );

	ASSERT( bytesAdded%LINESIZE == 0 );

	/* Now add blanks */
	while( bytesAdded%RECORDSIZE > 0 )
	{
		bytesAdded += fprintf(ioFITS_OUTPUT, "%-1s", " " );
	}
	return;
}

STATIC void punchFITS_ParamHeader( /* long *numParamValues, */ long nintparm, long naddparm )
{
	long numFields = 10;
	long naxis, naxis1, naxis2;
	char theValue[20];
	char theValue_temp[20];

	DEBUG_ENTRY( "punchFITS_ParamHeader()" );

	ASSERT( nintparm+naddparm <= LIMPAR );

	/* Make sure the previous blocks are the right size */
	ASSERT( bytesAdded%RECORDSIZE == 0 );

	naxis = 2;
	/* >>chng 06 aug 23, change to maximum number of parameter values. */
	naxis1 = 44+4*maxParamValues;
	naxis2 = nintparm+naddparm;

	bytesAdded += addKeyword_txt( "XTENSION", "'BINTABLE'",			"binary table extension", 0  );
	bytesAdded += addKeyword_num( "BITPIX"	, bitpix,					"8-bit bytes" );
	bytesAdded += addKeyword_num( "NAXIS"	, naxis,					"2-dimensional binary table" );
	bytesAdded += addKeyword_num( "NAXIS1"	, naxis1,				"width of table in bytes" );
	bytesAdded += addKeyword_num( "NAXIS2"	, naxis2,					"number of rows in table" );
	bytesAdded += addKeyword_num( "PCOUNT"	, pcount,					"size of special data area" );
	bytesAdded += addKeyword_num( "GCOUNT"	, gcount,					"one data group (required keyword)" );
	bytesAdded += addKeyword_num( "TFIELDS"	, numFields,			"number of fields in each row" );
	bytesAdded += addKeyword_txt( "TTYPE1"	, "'NAME    '",			"label for field   1", 0  );
	bytesAdded += addKeyword_txt( "TFORM1"	, "'12A     '",			"data format of the field: ASCII Character", 0  );
	bytesAdded += addKeyword_txt( "TTYPE2"	, "'METHOD  '",			"label for field   2", 0  );
	bytesAdded += addKeyword_txt( "TFORM2"	, "'J       '",			"data format of the field: 4-byte INTEGER", 0  );
	bytesAdded += addKeyword_txt( "TTYPE3"	, "'INITIAL '",			"label for field   3", 0  );
	bytesAdded += addKeyword_txt( "TFORM3"	, "'E       '",			"data format of the field: 4-byte REAL", 0  );
	bytesAdded += addKeyword_txt( "TTYPE4"	, "'DELTA   '",			"label for field   4", 0  );
	bytesAdded += addKeyword_txt( "TFORM4"	, "'E       '",			"data format of the field: 4-byte REAL", 0  );
	bytesAdded += addKeyword_txt( "TTYPE5"	, "'MINIMUM '",			"label for field   5", 0  );
	bytesAdded += addKeyword_txt( "TFORM5"	, "'E       '",			"data format of the field: 4-byte REAL", 0  );
	bytesAdded += addKeyword_txt( "TTYPE6"	, "'BOTTOM  '",			"label for field   6", 0  );
	bytesAdded += addKeyword_txt( "TFORM6"	, "'E       '",			"data format of the field: 4-byte REAL", 0  );
	bytesAdded += addKeyword_txt( "TTYPE7"	, "'TOP     '",			"label for field   7", 0  );
	bytesAdded += addKeyword_txt( "TFORM7"	, "'E       '",			"data format of the field: 4-byte REAL", 0  );
	bytesAdded += addKeyword_txt( "TTYPE8"	, "'MAXIMUM '",			"label for field   8", 0  );
	bytesAdded += addKeyword_txt( "TFORM8"	, "'E       '",			"data format of the field: 4-byte REAL", 0  );
	bytesAdded += addKeyword_txt( "TTYPE9"	, "'NUMBVALS'",			"label for field   9", 0  );
	bytesAdded += addKeyword_txt( "TFORM9"	, "'J       '",			"data format of the field: 4-byte INTEGER", 0  );
	bytesAdded += addKeyword_txt( "TTYPE10"	, "'VALUE   '",			"label for field  10", 0  );

	/* >>chng 06 aug 23, use maxParamValues instead of numParamValues */
	/* The size of this array is dynamic, set to size of the maximum of the numParamValues vector */
	sprintf( theValue_temp,		"%ld%s", maxParamValues, "E" );
	sprintf( theValue,		"%s%-8s%s", "'",theValue_temp,"'" );
	bytesAdded += addKeyword_txt( "TFORM10"	, theValue,			"data format of the field: 4-byte REAL", 0  );

	bytesAdded += addKeyword_txt( "EXTNAME"	, "'PARAMETERS'",		"name of this binary table extension", 0  );
	bytesAdded += addKeyword_txt( "HDUCLASS", "'OGIP    '",			"Format conforms to OGIP/GSFC conventions", 0  );
	bytesAdded += addKeyword_txt( "HDUCLAS1", "'XSPEC TABLE MODEL'","model spectra for XSPEC", 0  );
	bytesAdded += addKeyword_txt( "HDUCLAS2", "'PARAMETERS'",		"Extension containing paramter info", 0  );
	bytesAdded += addKeyword_txt( "HDUVERS"	, "'1.0.0   '",			"Version of format (OGIP memo OGIP-92-001)", 0  );
	bytesAdded += addKeyword_num( "NINTPARM", nintparm,				"Number of interpolation parameters" );
	bytesAdded += addKeyword_num( "NADDPARM", naddparm,				"Number of additional parameters" );
	/* After everything else */
	bytesAdded += fprintf(ioFITS_OUTPUT, "%-80s", "END" );

	ASSERT( bytesAdded%LINESIZE == 0 );

	/* Now add blanks */
	while( bytesAdded%RECORDSIZE > 0 )
	{
		bytesAdded += fprintf(ioFITS_OUTPUT, "%-1s", " " );
	}
	return;
}

STATIC void punchFITS_ParamData( char **paramNames,
					   long *paramMethods,
					   realnum **paramRange,
					   realnum **paramData,
					   long nintparm,
					   long naddparm,
					   long *numParamValues )
{
	long i, j;

	DEBUG_ENTRY( "punchFITS_ParamData()" );

	ASSERT( nintparm+naddparm <= LIMPAR );

	/* Now add the parameters data */
	for( i=0; i<nintparm+naddparm; i++ )
	{
		int32 numTemp;

#define	LOG2LINEAR 0

		paramMethods[i] = HtoNL(paramMethods[i]);
		/* >>chng 06 aug 23, numParamValues is now an array.  */
		numTemp = HtoNL(numParamValues[i]);

#if LOG2LINEAR
		/* change to linear */
		paramRange[i][0] = (realnum)pow( 10., (double)paramRange[i][0] );
		paramRange[i][1] = (realnum)pow( 10., (double)paramRange[i][1] );
		paramRange[i][2] = (realnum)pow( 10., (double)paramRange[i][2] );
		paramRange[i][3] = (realnum)pow( 10., (double)paramRange[i][3] );
		paramRange[i][4] = (realnum)pow( 10., (double)paramRange[i][4] );
		paramRange[i][5] = (realnum)pow( 10., (double)paramRange[i][5] );
#endif

#if !defined(_BIG_ENDIAN) 
		ByteSwap5( paramRange[i][0] );
		ByteSwap5( paramRange[i][1] );
		ByteSwap5( paramRange[i][2] );
		ByteSwap5( paramRange[i][3] );
		ByteSwap5( paramRange[i][4] );
		ByteSwap5( paramRange[i][5] );
#endif

		/* >>chng 06 aug 23, numParamValues is now an array.  */
		for( j=0; j<numParamValues[i]; j++ )
		{

#if LOG2LINEAR
			paramData[i][j] = (realnum)pow( 10., (double)paramData[i][j] );
#endif

#if !defined(_BIG_ENDIAN) 
			ByteSwap5( paramData[i][j] );
#endif
		}

		bytesAdded += fprintf(ioFITS_OUTPUT, "%-12s", paramNames[i] );
		bytesAdded += (long)fwrite( &paramMethods[i],	1,				  sizeof(int32),   ioFITS_OUTPUT );
		bytesAdded += (long)fwrite( paramRange[i],		1,				6*sizeof(realnum),   ioFITS_OUTPUT );
		bytesAdded += (long)fwrite( &numTemp,			1,				  sizeof(int32),   ioFITS_OUTPUT );
		/* >>chng 06 aug 23, numParamValues is now an array.  */
		bytesAdded += (long)fwrite( paramData[i],		1, (unsigned)numParamValues[i]*sizeof(realnum),   ioFITS_OUTPUT );

		for( j=numParamValues[i]+1; j<=maxParamValues; j++ )
		{
			realnum filler = -10.f;
			bytesAdded += (long)fwrite( &filler,		1, sizeof(realnum),   ioFITS_OUTPUT );
		}
	}

	/* Switch the endianness again */
	for( i=0; i<nintparm+naddparm; i++ )
	{
		paramMethods[i] = HtoNL(paramMethods[i]);

#if !defined(_BIG_ENDIAN) 
		ByteSwap5( paramRange[i][0] );
		ByteSwap5( paramRange[i][1] );
		ByteSwap5( paramRange[i][2] );
		ByteSwap5( paramRange[i][3] );
		ByteSwap5( paramRange[i][4] );
		ByteSwap5( paramRange[i][5] );
#endif

		/* >>chng 06 aug 23, numParamValues is now an array.  */
		for( j=0; j<numParamValues[i]; j++ )
		{
#if !defined(_BIG_ENDIAN) 
			ByteSwap5( paramData[i][j] );
#endif
		}
	}

	while( bytesAdded%RECORDSIZE > 0 )
	{
		int	tempInt = 0;
		bytesAdded += (long)fwrite( &tempInt, 1, 1,   ioFITS_OUTPUT );
	}
	return;
}

STATIC void punchFITS_EnergyHeader( long numEnergies )
{
	long numFields = 2;
	long naxis, naxis1, naxis2;

	DEBUG_ENTRY( "punchFITS_EnergyHeader()" );

	/* Make sure the previous blocks are the right size */
	ASSERT( bytesAdded%RECORDSIZE == 0 );

	naxis = 2;
	naxis1 = 2*sizeof(realnum);
	naxis2 = numEnergies;

	bytesAdded += addKeyword_txt( "XTENSION", "'BINTABLE'",			"binary table extension", 0 );
	bytesAdded += addKeyword_num( "BITPIX"	, bitpix,				"8-bit bytes" );
	bytesAdded += addKeyword_num( "NAXIS"	, naxis,				"2-dimensional binary table" );
	bytesAdded += addKeyword_num( "NAXIS1"	, naxis1,				"width of table in bytes" );
	bytesAdded += addKeyword_num( "NAXIS2"	, naxis2,				"number of rows in table" );
	bytesAdded += addKeyword_num( "PCOUNT"	, pcount,				"size of special data area" );
	bytesAdded += addKeyword_num( "GCOUNT"	, gcount,				"one data group (required keyword)" );
	bytesAdded += addKeyword_num( "TFIELDS"	, numFields,			"number of fields in each row" );
	bytesAdded += addKeyword_txt( "TTYPE1"	, "'ENERG_LO'",			"label for field   1", 0  );
	bytesAdded += addKeyword_txt( "TFORM1"	, "'E       '",			"data format of the field: 4-byte REAL", 0  );
	bytesAdded += addKeyword_txt( "TTYPE2"	, "'ENERG_HI'",			"label for field   2", 0  );
	bytesAdded += addKeyword_txt( "TFORM2"	, "'E       '",			"data format of the field: 4-byte REAL", 0  );
	bytesAdded += addKeyword_txt( "EXTNAME"	, "'ENERGIES'",			"name of this binary table extension", 0  );
	bytesAdded += addKeyword_txt( "HDUCLASS", "'OGIP    '",			"Format conforms to OGIP/GSFC conventions", 0  );
	bytesAdded += addKeyword_txt( "HDUCLAS1", "'XSPEC TABLE MODEL'","model spectra for XSPEC", 0  );
	bytesAdded += addKeyword_txt( "HDUCLAS2", "'ENERGIES'",			"Extension containing energy bin info", 0  );
	bytesAdded += addKeyword_txt( "HDUVERS"	, "'1.0.0   '",			"Version of format (OGIP memo OGIP-92-001)", 0  );
	/* After everything else */
	bytesAdded += fprintf(ioFITS_OUTPUT, "%-80s", "END" );

	ASSERT( bytesAdded%LINESIZE == 0 );

	while( bytesAdded%RECORDSIZE > 0 )
	{
		bytesAdded += fprintf(ioFITS_OUTPUT, "%-1s", " " );
	}
	return;
}

STATIC void punchFITS_EnergyData( const vector<realnum>& Energies, long EnergyOffset )
{
	DEBUG_ENTRY( "punchFITS_EnergyData()" );

	/* Now add the energies data */
	for( unsigned long i=0; i < Energies.size(); i++ )
	{
		realnum EnergyLow, EnergyHi;
		ASSERT( i+EnergyOffset < (unsigned)rfield.nupper );
		/* Convert all of these to kev */
		EnergyLow = 0.001f*(realnum)EVRYD*(Energies[i]-rfield.widflx[i+EnergyOffset]/2.f);

		if( i == Energies.size()-1 )
		{
			EnergyHi = 0.001f*(realnum)EVRYD*(Energies[i] + rfield.widflx[i+EnergyOffset]/2.f);
		}
		else
		{
			EnergyHi = 0.001f*(realnum)EVRYD*(Energies[i+1] - rfield.widflx[i+EnergyOffset+1]/2.f);
		}

#if !defined(_BIG_ENDIAN) 
		ByteSwap5(EnergyLow);
		ByteSwap5(EnergyHi);
#endif

		bytesAdded += (long)fwrite( &EnergyLow	, 1, sizeof(realnum), ioFITS_OUTPUT );
		bytesAdded += (long)fwrite( &EnergyHi	, 1, sizeof(realnum), ioFITS_OUTPUT );
	}

	while( bytesAdded%RECORDSIZE > 0 )
	{
		int	tempInt = 0;
		bytesAdded += (long)fwrite( &tempInt, 1, 1,   ioFITS_OUTPUT );
	}
	return;
}

STATIC void punchFITS_SpectraHeader( bool lgAddModel, long nintparm, long naddparm, long totNumModels, long numEnergies )
{
	long i, numFields = 2+naddparm;
	long naxis, naxis1, naxis2;
	char theKeyword1[8];
	char theKeyword2[8];
	char theKeyword3[8];
	char theValue1[20];
	char theValue2[20];
	char theValue2temp[20];
	char theValue[20];
	char theValue_temp[20];
	char theComment1[47];

	DEBUG_ENTRY( "punchFITS_SpectraHeader()" );

	ASSERT( nintparm + naddparm <= LIMPAR );

	/* Make sure the previous blocks are the right size */
	ASSERT( bytesAdded%RECORDSIZE == 0 );

	naxis = 2;
	naxis1 = ( numEnergies*(naddparm+1) + nintparm ) * (long)sizeof(realnum);
	naxis2 = totNumModels; 

	bytesAdded += addKeyword_txt( "XTENSION", "'BINTABLE'",			"binary table extension", 0  );
	bytesAdded += addKeyword_num( "BITPIX"	, bitpix,				"8-bit bytes" );
	bytesAdded += addKeyword_num( "NAXIS"	, naxis,				"2-dimensional binary table" );
	bytesAdded += addKeyword_num( "NAXIS1"	, naxis1,				"width of table in bytes" );
	bytesAdded += addKeyword_num( "NAXIS2"	, naxis2,				"number of rows in table" );
	bytesAdded += addKeyword_num( "PCOUNT"	, pcount,				"size of special data area" );
	bytesAdded += addKeyword_num( "GCOUNT"	, gcount,				"one data group (required keyword)" );
	bytesAdded += addKeyword_num( "TFIELDS"	, numFields,			"number of fields in each row" );

	/******************************************/
	/* These are the interpolation parameters */
	/******************************************/
	bytesAdded += addKeyword_txt( "TTYPE1"	, "'PARAMVAL'",			"label for field   1", 0 );
	/* The size of this array is dynamic, set to size of nintparm */
	sprintf( theValue2temp,		"%ld%s", nintparm, "E" );
	sprintf( theValue2,		"%s%-8s%s", "'",theValue2temp,"'" );
	bytesAdded += addKeyword_txt( "TFORM1"	, theValue2,			"data format of the field: 4-byte REAL", 0  );

	/******************************************/
	/* This is the interpolated spectrum      */	
	/******************************************/
	bytesAdded += addKeyword_txt( "TTYPE2"	, "'INTPSPEC'",	"label for field 2", 0  );
	/* The size of this array is dynamic, set to size of numEnergies */
	sprintf( theValue_temp,		"%ld%s", numEnergies, "E" );
	sprintf( theValue,		"%s%-8s%s", "'",theValue_temp,"'" );
	bytesAdded += addKeyword_txt( "TFORM2"	, theValue,	"data format of the field: 4-byte REAL", 0  );
	bytesAdded += addKeyword_txt( "TUNIT2"	, ModelUnits[lgAddModel],	"physical unit of field", 0  );

	/******************************************/
	/* These are the additional parameters    */
	/******************************************/
	for( i=1; i<=naddparm; i++ )
	{
		sprintf( theKeyword1,	"%s%ld", "TTYPE", i+2 );
		sprintf( theKeyword2,	"%s%ld", "TFORM", i+2 );
		sprintf( theKeyword3,	"%s%ld", "TUNIT", i+2 );

		sprintf( theValue1,		"%s%2.2ld%s", "'ADDSP", i, "'" );
		/* The size of this array is dynamic, set to size of numEnergies */
		sprintf( theValue2temp,		"%ld%s", numEnergies, "E" );
		sprintf( theValue2,		"%s%-8s%s", "'",theValue2temp,"'" );

		sprintf( theComment1,	"%s%ld", "label for field ", i+2 );

		bytesAdded += addKeyword_txt( theKeyword1	, theValue1,		theComment1, 0  );
		bytesAdded += addKeyword_txt( theKeyword2	, theValue2,		"data format of the field: 4-byte REAL", 0  );
		bytesAdded += addKeyword_txt( theKeyword3	, ModelUnits[lgAddModel],	"physical unit of field", 0  );
	}

	bytesAdded += addKeyword_txt( "EXTNAME"	, "'SPECTRA '",			"name of this binary table extension", 0  );
	bytesAdded += addKeyword_txt( "HDUCLASS", "'OGIP    '",			"Format conforms to OGIP/GSFC conventions", 0  );
	bytesAdded += addKeyword_txt( "HDUCLAS1", "'XSPEC TABLE MODEL'","model spectra for XSPEC", 0  );
	bytesAdded += addKeyword_txt( "HDUCLAS2", "'MODEL SPECTRA'",	"Extension containing model spectra", 0  );
	bytesAdded += addKeyword_txt( "HDUVERS"	, "'1.0.0   '",			"Version of format (OGIP memo OGIP-92-001)", 0  );
	/* After everything else */
	bytesAdded += fprintf(ioFITS_OUTPUT, "%-80s", "END" );

	ASSERT( bytesAdded%LINESIZE == 0 );

	while( bytesAdded%RECORDSIZE > 0 )
	{
		bytesAdded += fprintf(ioFITS_OUTPUT, "%-1s", " " );
	}
	return;
}

STATIC void punchFITS_SpectraData( realnum **interpParameters, multi_arr<realnum,3>& theSpectrum, int option,
						   long totNumModels, long numEnergies, long nintparm, long naddparm )
{
	long i;
	long naxis2 = totNumModels;

	DEBUG_ENTRY( "punchFITS_SpectraData()" );

	ASSERT( nintparm + naddparm <= LIMPAR );

	/* Now add the spectra data */
	for( i=0; i<naxis2; i++ )
	{

#if !defined(_BIG_ENDIAN)
		for( long j = 0; j<numEnergies; j++ )
		{
			ByteSwap5( theSpectrum[option][i][j] );
		}

		for( long j = 0; j<nintparm; j++ )
		{
			ByteSwap5( interpParameters[i][j] );
		}
#endif

		/* The interpolation parameters vector */
		bytesAdded += (long)fwrite( interpParameters[i],	1,	 (unsigned)nintparm*sizeof(realnum), ioFITS_OUTPUT );
		/* The interpolated spectrum */
		bytesAdded += (long)fwrite( &theSpectrum[option][i][0],			1, (unsigned)numEnergies*sizeof(realnum), ioFITS_OUTPUT );

#if !defined(_BIG_ENDIAN)
		/* Switch the endianness back to native. */
		for( long j = 0; j<numEnergies; j++ )
		{
			ByteSwap5( theSpectrum[option][i][j] );
		}

		for( long j = 0; j<nintparm; j++ )
		{
			ByteSwap5( interpParameters[i][j] );
		}
#endif

		/* >>chng 06 aug 23, disable additional parameters for now */
		if( naddparm > 0 )
		{
			/* The additional parameters */

			/* bytesAdded += (long)fwrite( &theSpectrum[option][i][0],		1, (unsigned)numEnergies*sizeof(realnum), ioFITS_OUTPUT ); */
			/* \todo 2	must create another array if we are to save additional parameter information. */
			fprintf( ioQQQ, " Additional parameters not currently supported.\n" );
			cdEXIT( EXIT_FAILURE );
		}
	}

	while( bytesAdded%RECORDSIZE > 0 )
	{
		int	tempInt = 0;
		bytesAdded += (long)fwrite( &tempInt, 1, 1,   ioFITS_OUTPUT );
	}
	return;
}

STATIC void punchFITS_GenericHeader( long numEnergies )
{
	long numFields = 2;
	long naxis, naxis1, naxis2;

	DEBUG_ENTRY( "punchFITS_GenericHeader()" );

	/* Make sure the previous blocks are the right size */
	ASSERT( bytesAdded%RECORDSIZE == 0 );

	naxis = 2;
	naxis1 = numFields*(long)sizeof(realnum);
	naxis2 = numEnergies;

	bytesAdded += addKeyword_txt( "XTENSION", "'BINTABLE'",			"binary table extension", 0 );
	bytesAdded += addKeyword_num( "BITPIX"	, bitpix,				"8-bit bytes" );
	bytesAdded += addKeyword_num( "NAXIS"	, naxis,				"2-dimensional binary table" );
	bytesAdded += addKeyword_num( "NAXIS1"	, naxis1,				"width of table in bytes" );
	bytesAdded += addKeyword_num( "NAXIS2"	, naxis2,				"number of rows in table" );
	bytesAdded += addKeyword_num( "PCOUNT"	, pcount,				"size of special data area" );
	bytesAdded += addKeyword_num( "GCOUNT"	, gcount,				"one data group (required keyword)" );
	bytesAdded += addKeyword_num( "TFIELDS"	, numFields,			"number of fields in each row" );
	bytesAdded += addKeyword_txt( "TTYPE1"	, "'ENERGY  '",			"label for field   1", 0  );
	bytesAdded += addKeyword_txt( "TFORM1"	, "'E       '",			"data format of the field: 4-byte REAL", 0  );
	bytesAdded += addKeyword_txt( "TTYPE2"	, "'TRN_SPEC'",			"label for field   2", 0  );
	bytesAdded += addKeyword_txt( "TFORM2"	, "'E       '",			"data format of the field: 4-byte REAL", 0  );
	bytesAdded += addKeyword_txt( "EXTNAME"	, "'SPECTRA '",			"name of this binary table extension", 0  );
	bytesAdded += addKeyword_txt( "HDUCLASS", "'OGIP    '",			"Format conforms to OGIP/GSFC conventions", 0  );
	bytesAdded += addKeyword_txt( "HDUCLAS1", "'XSPEC TABLE MODEL'","model spectra for XSPEC", 0  );
	bytesAdded += addKeyword_txt( "HDUCLAS2", "'ENERGIES'",			"Extension containing energy bin info", 0  );
	bytesAdded += addKeyword_txt( "HDUVERS"	, "'1.0.0   '",			"Version of format (OGIP memo OGIP-92-001)", 0  );
	/* After everything else */
	bytesAdded += fprintf(ioFITS_OUTPUT, "%-80s", "END" );

	ASSERT( bytesAdded%LINESIZE == 0 );

	while( bytesAdded%RECORDSIZE > 0 )
	{
		bytesAdded += fprintf(ioFITS_OUTPUT, "%-1s", " " );
	}
	return;
}

STATIC void punchFITS_GenericData( long numEnergies, long ipLoEnergy, long ipHiEnergy )
{
	long i;

	DEBUG_ENTRY( "punchFITS_GenericData()" );

	realnum *TransmittedSpectrum;

	TransmittedSpectrum = (realnum*)MALLOC(sizeof(realnum)*(unsigned)(numEnergies) );

	cdSPEC2( 8, numEnergies, ipLoEnergy, ipHiEnergy, TransmittedSpectrum );

	/* Now add the energies data */
	for( i=0; i<numEnergies; i++ )
	{
		realnum Energy;
		Energy = rfield.AnuOrg[i];

#if !defined(_BIG_ENDIAN) 
		ByteSwap5(Energy);
		ByteSwap5(TransmittedSpectrum[i]);
#endif

		bytesAdded += (long)fwrite( &Energy	,				1, sizeof(realnum), ioFITS_OUTPUT );
		bytesAdded += (long)fwrite( &TransmittedSpectrum[i],1, sizeof(realnum), ioFITS_OUTPUT );
	}

	while( bytesAdded%RECORDSIZE > 0 )
	{
		int	tempInt = 0;
		bytesAdded += (long)fwrite( &tempInt, 1, 1,   ioFITS_OUTPUT );
	}

	free( TransmittedSpectrum );
	return;
}

STATIC void writeCloudyDetails( void )
{
	char timeString[30]="";
	char tempString[70];
	time_t now;
	long i, j, k;

	/* usually print date and time info - do not if "no times" command entered, 
	 * which set this flag false */
	now = time(NULL);
	if( prt.lgPrintTime ) 
	{
		/* now add date of this run */
		/* now print this time at the end of the string.  the system put cr at the end of the string */
		strcpy( timeString , ctime(&now) );
	}
	/* ctime puts a carriage return at the end, but we can't have that in a fits file.
	 * remove the carriage return here. */
	for( i=0; i<30; i++ )
	{
		if( timeString[i] == '\n' )
		{
			timeString[i] = ' ';
		}
	}

	strcpy( tempString, "Generated by Cloudy " );
	// strncat guarantees that terminating 0 byte will always be written...
	strncat( tempString, t_version::Inst().chVersion, sizeof(tempString)-strlen(tempString) );
	bytesAdded += addComment( tempString );
	bytesAdded += addComment( t_version::Inst().chInfo );
	strcpy( tempString, "--- " );
	strcat( tempString, timeString );
	bytesAdded += addComment( tempString );
	bytesAdded += addComment( "Input string was as follows: " );
	/* >>chng 05 nov 24, from <nSave to <=nSave bug caught by PvH */
	for( i=0; i<=input.nSave; i++ )
	{
		char firstLine[70], extraLine[65];

		for( j=0; j<INPUT_LINE_LENGTH; j++ )
		{
			if( input.chCardSav[i][j] =='\0' )
				break;
		}

		ASSERT( j < 200 );
		for( k=0; k< MIN2(69, j); k++ )
		{
			firstLine[k] = input.chCardSav[i][k];
		}
		firstLine[k] = '\0';
		bytesAdded += addComment( firstLine );
		if( j >= 69 )
		{
			for( k=69; k< 133; k++ )
			{
				extraLine[k-69] = input.chCardSav[i][k];
			}
			/* >> chng 06 jan 05, this was exceeding array bounds. */
			extraLine[64] = '\0'; 
			strcpy( tempString, "more " );
			strcat( tempString, extraLine );
			bytesAdded += addComment( tempString );
			if( j >= 133 )
			{
				for( k=133; k< 197; k++ )
				{
					extraLine[k-133] = input.chCardSav[i][k];
				}
				extraLine[64] = '\0';
				strcpy( tempString, "more " );
				strcat( tempString, extraLine );
				bytesAdded += addComment( tempString );
			}
		}
	}

	return;
}

STATIC long addKeyword_txt( const char *theKeyword, const void *theValue, const char *theComment, long Str_Or_Log )
{
	long numberOfBytesWritten = 0;

	DEBUG_ENTRY( "addKeyword_txt()" );

	/* False means string, true means logical */
	if( Str_Or_Log == 0 )
	{
		numberOfBytesWritten = fprintf(ioFITS_OUTPUT, "%-8s%-2s%-20s%3s%-47s",
			theKeyword,
			"= ",
			(char *)theValue,
			" / ",
			theComment );
	}
	else
	{
		ASSERT( Str_Or_Log == 1 );
		numberOfBytesWritten = fprintf(ioFITS_OUTPUT, "%-8s%-2s%20s%3s%-47s",
			theKeyword,
			"= ",
			(char *)theValue,
			" / ",
			theComment );
	}

	ASSERT( numberOfBytesWritten%LINESIZE == 0 );
	return numberOfBytesWritten;
}

STATIC long addKeyword_num( const char *theKeyword, long theValue, const char *theComment)
{
	long numberOfBytesWritten = 0;

	DEBUG_ENTRY( "addKeyword_num()" );

	numberOfBytesWritten = fprintf(ioFITS_OUTPUT, "%-8s%-2s%20ld%3s%-47s",
		theKeyword,
		"= ",
		theValue,
		" / ",
		theComment );

	ASSERT( numberOfBytesWritten%LINESIZE == 0 );
	return numberOfBytesWritten;
}

long addComment( const char *CommentToAdd )
{
	long i, numberOfBytesWritten = 0;
	char tempString[70] = "                                                                     ";

	DEBUG_ENTRY( "addComment()" );

	strncpy( &tempString[0], CommentToAdd, 69 );
	ASSERT( (int)strlen( tempString ) <= 70 );	

	/* tabs violate FITS standard, replace them with spaces. */
	for( i=0; i<69; i++ )
	{
		if( tempString[i] == '\t' )
		{
			tempString[i] = ' ';
		}
	}

	numberOfBytesWritten = fprintf(ioFITS_OUTPUT, "COMMENT   %-70s", tempString );

	ASSERT( numberOfBytesWritten%LINESIZE == 0 );
	return numberOfBytesWritten;
}
