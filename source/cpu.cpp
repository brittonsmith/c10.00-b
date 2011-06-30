/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/**\file cpu.cpp implement hardware dependent definitions */
#if defined(__HP_aCC)
/* this is for the HP compiler on the sdx */
extern "C" unsigned long fegettrapenable();
extern "C" void fesettrapenable(unsigned long); 
#endif

#if defined(__ia64) && defined(__INTEL_COMPILER)
extern "C" unsigned long fpgetmask();
extern "C" void fpsetmask(unsigned long);
#endif	

#if defined(__sun) || defined(__sgi)
#include <ieeefp.h>
#if defined(HAVE_SUNMATH) || defined(FLUSH_DENORM_TO_ZERO)
#include <sunmath.h>
#endif
#endif

#if defined(__alpha) && defined(__linux) && defined(__GNUC__)
#define __USE_GNU
#include <fenv.h>
#endif

#ifdef __unix
#include <unistd.h>
#endif

/* the redefinition of float in cddefines.h can cause problems in system headers
 * hence these includes MUST come after the system header includes above */
#include "cddefines.h"
#include "path.h"
#include "trace.h"

STATIC NORETURN void AbortErrorMessage( const char* fname, vector<string>& PathList, access_scheme scheme );

/* NB NB - this constructor needs to be called before any of the user code is executed !! */
t_cpu::t_cpu()
{
	DEBUG_ENTRY( "t_cpu()" );

#ifdef CATCH_SIGNAL
	// set up signal handlers so that we can properly terminate MPI runs...

#	ifdef __unix
	p_action.sa_handler = &signal_handler;
	sigemptyset( &p_action.sa_mask );
	p_action.sa_flags = SA_NODEFER;

	p_default.sa_handler = SIG_DFL;
	sigemptyset( &p_default.sa_mask );
	p_default.sa_flags = SA_NODEFER;

	for( int sig=1; sig <= 31; sig++ )
	{
		// is the signal valid?
		if( sigaction( sig, NULL, NULL ) == 0 )
			// these two are for suspending and resuming a job
			if( sig != SIGSTOP && sig != SIGCONT )
				sigaction( sig, action(), NULL );
	}
#	endif

#	ifdef _MSC_VER
	signal( SIGABRT, &signal_handler );
	signal( SIGFPE, &signal_handler );
	signal( SIGILL, &signal_handler );
	signal( SIGINT, &signal_handler );
	signal( SIGSEGV, &signal_handler );
	signal( SIGTERM, &signal_handler );
#	endif
#endif

	/* >>chng 05 dec 14, add test of endianness of the CPU, PvH */
	endian.c[0] = 0x12;
	endian.c[1] = 0x34;
	endian.c[2] = 0x56;
	endian.c[3] = 0x78;

	/* >>chng 05 dec 15, add signaling NaN for float and double to cpu struct, PvH */
	/* in C++ this should be replaced by numeric_limits<TYPE>::signaling_NaN() */
	if( sizeof(sys_float) == 4 )
	{
#		ifdef __mips
		/* definition of signaling and quiet NaN is reversed on MIPS */
		Float_SNaN_Value = 0xffffffff;
#		else
		if( big_endian() || little_endian() )
		{
			/* this should work on most modern CPU's */
			Float_SNaN_Value = 0xffbfffff;
		}
		else
		{
			/* this is an unusual CPU -> bit pattern for SNaN is unknown */
			Float_SNaN_Value = -1;
		}
#		endif
	}
	else
	{
		/* this is an unusual CPU -> bit pattern for SNaN is unknown */
		Float_SNaN_Value = -1;
	}

#	ifdef HAVE_INT64

	if( sizeof(double) == 8 )
	{
#		ifdef __mips
		/* definition of signaling and quiet NaN is reversed on MIPS */
		Double_SNaN_Value = 0xffffffffffffffff;
#		else
		/* this should work on most modern CPU's */
		Double_SNaN_Value = 0xfff7ffffffbfffff;
#		endif
	}
	else
	{
		/* this is an unusual CPU -> bit pattern for SNaN is unknown */
		Double_SNaN_Value = -1;
	}

#	else

	if( sizeof(double) == 8 )
	{
#		ifdef __mips
		/* definition of signaling and quiet NaN is reversed on MIPS */
		Double_SNaN_Value[0] = 0xffffffff;
		Double_SNaN_Value[1] = 0xffffffff;
#		else
		if( big_endian() )
		{
			/* this should work on most modern CPU's */
			Double_SNaN_Value[0] = 0xfff7ffff;
			Double_SNaN_Value[1] = 0xffbfffff;
		}
		else if( little_endian() )
		{
			/* this should work on most modern CPU's */
			Double_SNaN_Value[0] = 0xffbfffff;
			Double_SNaN_Value[1] = 0xfff7ffff;
		}
		else
		{
			/* this is an unusual CPU -> bit pattern for SNaN is unknown */
			Double_SNaN_Value[0] = -1;
			Double_SNaN_Value[1] = -1;
		}
#		endif
	}
	else
	{
		/* this is an unusual CPU -> bit pattern for SNaN is unknown */
		Double_SNaN_Value[0] = -1;
		Double_SNaN_Value[1] = -1;
	}

#	endif

	/* set FP environment to trap FP exceptions */
	enable_traps();

	ioStdin = stdin;
	ioQQQ = stdout;
	ioPrnErr = stderr;
	lgPrnErr = false;

	test_float = FLT_MIN;
	test_double = DBL_MIN;

	/* default is for failed asserts not to abort */
	p_lgAssertAbort = false;

	const char *str;

	/* determine the no. of CPUs on this machine; used by PHYMIR, grid command, .... */
#	if defined(_SC_NPROCESSORS_ONLN)  /* Linux, Sun Sparc, DEC Alpha */
	n_avail_CPU = sysconf(_SC_NPROCESSORS_ONLN);
#	elif defined(_SC_NPROC_ONLN)      /* SGI Iris */
	n_avail_CPU = sysconf(_SC_NPROC_ONLN);
#	elif defined(_SC_CRAY_NCPU)       /* Cray */
	n_avail_CPU = sysconf(_SC_CRAY_NCPU);
#	elif defined(_WIN32)
	str = getenv( "NUMBER_OF_PROCESSORS" );
	if( str != NULL )
	{
		int found = sscanf( str, "%ld", &n_avail_CPU );
		if( found != 1 )
			n_avail_CPU = 1;
	}
	else
	{
		n_avail_CPU = 1;
	}
#	else                              
	/* Other systems, supply no. of CPUs on OPTIMIZE PHYMIR command line */
	n_avail_CPU = 1;
#	endif
	/* the constructor is run before MPI starts, so this flag and the rank are not available yet */
	p_lgMPI = false;
	/* the default is for all ranks to cooperate on the same sim */
	p_lgMPISingleRankMode = false;
	n_rank = 0;

#	ifdef _WIN32	
	str = getenv( "COMPUTERNAME" );
#	else
	str = getenv( "HOSTNAME" );
#	endif

	if( str != NULL )
		strncpy( HostName, str, STDLEN );
	else
		strncpy( HostName, "unknown", STDLEN );
	HostName[STDLEN-1] = '\0';

	/* pick up the path from the environment, if set by user */
	const char *path = getenv( "CLOUDY_DATA_PATH" );

	/* if the environment variable was not set, the default set in path.h takes effect */
	string chSearchPathRaw = ( path != NULL ) ? string( path ) : string( CLOUDY_DATA_PATH );

#	ifdef _WIN32
	string separator( ";" );
	char chExpectedEnd = '\\';
#	else
	string separator( ":" );
	char chExpectedEnd = '/';
#	endif

	chSearchPath.push_back( "" ); // the current working directory should be first and last
	Split( chSearchPathRaw, separator, chSearchPath, SPM_RELAX );
	chSearchPath.push_back( "" );

	for( vector<string>::size_type i=0; i < chSearchPath.size(); ++i )
	{
		if( chSearchPath[i].length() > 0 )
		{
			/* get last valid char */
			char chEnd = *chSearchPath[i].rbegin();

			/* make sure path ends with directory separator */
			if( chEnd != chExpectedEnd )
				chSearchPath[i] += chExpectedEnd;
		}
	}

	nFileDone = 0;
}

void t_cpu::enable_traps() const
{
	/* >>chng 01 aug 07, added code to circumvent math library bug with g++ on
	 * alpha-linux machines, see bug report 51072 on http://bugzilla.redhat.com, PvH */
	/* >>chng 01 apr 17, added code for Solaris and SGI operating systems, PvH */
	/* this routine contains no code for alphas or crays, they do not need
	 * special code to enable FP exceptions since they are enabled by default */

	/* there is no command line option on MS Visual Studio to force crash */
#	if defined(_MSC_VER)
	volatile unsigned int NewMask;

	/* | is a bitwise inclusive or, turns on bits
	 * 0|0 = 0
	 * 0|1 = 1|0 = 1|1 = 1 */
	NewMask = _EM_ZERODIVIDE | _EM_OVERFLOW | _EM_INVALID;
	/* ~ is the unary bitwise complement - all bits flip */ 
	NewMask = ~NewMask;
	_controlfp( NewMask , _MCW_EM );

	/* this is the code for Linux PC (but not Linux alpha) to force crash */
	/* >>chng 04 apr 26, added support for AMD64, enable FPE traps for SSE/SSE2, PvH */
	/* >>chng 06 aug 12, added support for Apple MacOSX, and hopefully also Solaris x86, PvH */
#	elif defined(__GNUC__) && ( defined(__i386) || defined(__amd64) )
	volatile unsigned int Old_Mask, New_Mask;
#	if defined(__SSE__) || defined(__SSE2__)
	volatile unsigned int SSE_Mask;
#	endif

#	define _FPU_MASK_IM  0x01  /* Invalid          */
#	define _FPU_MASK_DM  0x02  /* Denormalized     */
#	define _FPU_MASK_ZM  0x04  /* Division-by-zero */
#	define _FPU_MASK_OM  0x08  /* Overflow         */
#	define _FPU_MASK_UM  0x10  /* Underflow        */
#	define _FPU_MASK_PM  0x20  /* Inexact          */

	/*  | is a bitwise inclusive or, turns on bits     */
	/* 0|0 = 0                                         */
	/* 0|1 = 1|0 = 1|1 = 1                             */

	/*  ~ is the unary bitwise complement - all bits flip */

	/* this enables FPE traps for regular i387 FP instructions */

	volatile unsigned int UnMask = ~((unsigned int)( _FPU_MASK_ZM | _FPU_MASK_IM  | _FPU_MASK_OM ));

	__asm__ volatile("fnstcw %0" : "=m" (*&Old_Mask));

	New_Mask = Old_Mask & UnMask;

	__asm__ volatile("fldcw %0" : : "m" (*&New_Mask));

#	if defined(__SSE__) || defined(__SSE2__)

#	if defined(FLUSH_DENORM_TO_ZERO)
	/* using this causes denormalized numbers to be flushed to zero,
	 * which will speed up the code on Pentium 4 processors */
	SSE_Mask = 0x9900;
#	else
	/* this version allows denormalized numbers to be retained */
	SSE_Mask = 0x1900;
#	endif

	/* this enables FPE traps for SSE/SSE2 instructions */

	__asm__ volatile( "ldmxcsr %0" : : "m" (*&SSE_Mask) );

#	endif

	/* this is for IA64 systems running g++ or icc (e.g. SGI, HP, ...) */
#	elif defined(__ia64)

#	define FPSR_TRAP_VD     (1 << 0)        /* invalid op trap disabled */
#	define FPSR_TRAP_DD     (1 << 1)        /* denormal trap disabled */
#	define FPSR_TRAP_ZD     (1 << 2)        /* zero-divide trap disabled */
#	define FPSR_TRAP_OD     (1 << 3)        /* overflow trap disabled */
#	define FPSR_TRAP_UD     (1 << 4)        /* underflow trap disabled */
#	define FPSR_TRAP_ID     (1 << 5)        /* inexact trap disabled */

#	define FPSR_SF0_FTZ     (1 << 6)        /* flush denormalized numbers to zero */

#	if defined(__GNUC_EXCL__)
	/* __asm__ instructions are not supported by icc as of v9.0 */
#	define _IA64_REG_AR_FPSR    40

#	define ia64_getreg( regnum )       __asm__ volatile( "mov %0=ar%1" : "=r" (fpsr) : "i"(regnum) )
#	define ia64_setreg( regnum, val )  __asm__ volatile( "mov ar%0=%1" :: "i" (regnum), "r"(val): "memory" )
#	define ia64_serialize              __asm__ volatile( "srlz.i" );

	volatile unsigned long fpsr, flags = FPSR_TRAP_VD | FPSR_TRAP_ZD | FPSR_TRAP_OD;

	ia64_getreg( _IA64_REG_AR_FPSR );
	fpsr &= ~flags;
#	if defined(FLUSH_DENORM_TO_ZERO)
	fpsr |= FPSR_SF0_FTZ;
#	endif
	ia64_setreg( _IA64_REG_AR_FPSR, fpsr );
	/* this prevents RAW and WAW dependency violations in case this ever gets inlined... */
	ia64_serialize;

#	elif defined(__INTEL_COMPILER)
	/* this is for icc on IA64 SGI machines */
	unsigned long fpsr = fpgetmask();
	fpsr |= FPSR_TRAP_VD | FPSR_TRAP_ZD | FPSR_TRAP_OD;
	fpsetmask( fpsr );
#	elif defined(__HP_aCC)
	/* this is for the HP compiler on the sdx */
	unsigned long fpsr = fegettrapenable(); 
	fpsr |= FPSR_TRAP_VD | FPSR_TRAP_ZD | FPSR_TRAP_OD;
	fesettrapenable( fpsr ); 
#	endif /* defined(__GNUC_EXCL__) */

	/* this is for Solaris and SGI to force crash */
#	elif defined(__sun) || defined(__sgi)

	fp_except mask;

	/* >>chng 05 dec 30, accept FLUSH_DENORM_TO_ZERO as a synonym for HAVE_SUNMATH, PvH */
#	if defined(HAVE_SUNMATH) || defined(FLUSH_DENORM_TO_ZERO)

	/* >>chng 01 oct 09, disable gradual underflow on ultrasparc whith g++
	 * (only needed for versions < 3.1 or >= 4.3.0, see Note 1).
	 *
	 * compile this module with:
	 *     g++ [ other options... ] -I<include-dir> -DHAVE_SUNMATH -c cpu.cpp
	 * link the program with:
	 *     g++ -L<library-dir> -o cloudy.exe *.o -lsunmath
	 *
	 * you probably need to use -I<include-dir> and -L<library-dir> to point the
	 * compiler/linker to the location of the sunmath.h header file and libsunmath.so
	 * library (e.g., -I/opt/SUNWspro/prod/include/cc -L/opt/SUNWspro/lib; note that
	 * the actual location may vary from one installation to another).
	 * See also bug report 4487 on http://gcc.gnu.org/bugzilla/
	 *
	 * Note 1: Starting with g++ 3.1, bug 4487 has been solved: -funsafe-math-optimizations
	 * will automatically disable gradual underflow. Hence using nonstandard_arithmetic()
	 * is no longer necessary. The option -funsafe-math-optimizations should be included
	 * both when compiling and linking:
	 *
	 * g++ [ other options... ] -funsafe-math-optimizations -c *.c
	 * g++ [ other options... ] -funsafe-math-optimizations -o cloudy.exe *.o
	 *
	 * Starting with g++ 4.3.0 the -funsafe-math-optimizations option can no longer be
	 * used as it implicitly enables -fno-trapping-math, which is unsafe for Cloudy
	 * because we do trap floating point exceptions.
	 *
	 * Note 2: Don't use nonstandard_arithmetic() with CC (the SunWorks/Forte compiler);
	 * use the -fast commandline option instead to disable gradual underflow (or use
	 * -fnonstd if you don't want all the other options enabled by -fast). The option
	 * -fast (or -fnonstd) should be included both when compiling and linking:
	 *
	 * CC [ other options... ] -fast -c *.c
	 * CC -fast -o cloudy.exe *.o
	 *
	 * PvH */
	nonstandard_arithmetic();
#	endif

	/* enable floating point exceptions on sun and sgi */
	mask = fpgetmask();
	mask = mask | FP_X_INV | FP_X_OFL | FP_X_DZ;
	fpsetmask(mask);

#	elif defined(__alpha) && defined(__linux) && defined(__GNUC__)

	/* the following is not supported on all hardware platforms, but certainly for EV56
	 * and later. earlier versions may work as well, but that has not been tested.
	 * for details see https://bugzilla.redhat.com/bugzilla/show_bug.cgi?id=51072 */
#	ifdef FE_NONIEEE_ENV
	/* this prevents the infamous math library bug when compiling with gcc on alpha-linux
	 * machines. if this doesn't work on your system, the only alternative is to link
	 * against the Compaq math library: gcc *.o -lcpml -lm, or use ccc itself, PvH */
	fesetenv(FE_NONIEEE_ENV);
#	endif

#	endif
	return;
}

void t_cpu::signal_handler(int sig)
{
	// when an FPE is caught, the mask is reset...
	cpu.enable_traps();
#	ifdef _MSC_VER
	// at this point the signal handler has reverted to the default handler
	signal( sig, &signal_handler );
#	endif
	throw bad_signal( sig );
}


void t_cpu::printDataPath() const
{
	fprintf(ioQQQ, "The path is:\n");
	for( vector<string>::size_type i=1; i < cpu.chSearchPath.size()-1; ++i )
		fprintf( ioQQQ, "   ==%s==\n", cpu.chSearchPath[i].c_str() );
}

// this routine generates a list of all full paths to the locations where we should look for the file
void t_cpu::getPathList( const char* fname, vector<string>& PathList, access_scheme scheme ) const
{
	DEBUG_ENTRY( "getPathList()" );

	vector<string>::size_type begin, end;

	switch( scheme )
	{
	case AS_DATA_ONLY:
	case AS_DATA_ONLY_TRY:
	case AS_DATA_OPTIONAL:
		begin = 1;
		end = cpu.chSearchPath.size()-1;
		break;
	case AS_DATA_LOCAL:
	case AS_DATA_LOCAL_TRY:
		begin = 1;
		end = cpu.chSearchPath.size();
		break;
	case AS_LOCAL_DATA:
	case AS_LOCAL_DATA_TRY:
		begin = 0;
		end = cpu.chSearchPath.size()-1;
		break;
	case AS_LOCAL_ONLY:
	case AS_LOCAL_ONLY_TRY:
		begin = 0;
		end = 1;
		break;
	default:
		TotalInsanity();
	}

	PathList.clear();
	string FileName( fname );
	for( vector<string>::size_type i=begin; i < end; ++i )
		PathList.push_back( cpu.chSearchPath[i] + FileName );
}

STATIC NORETURN void AbortErrorMessage( const char* fname, vector<string>& PathList, access_scheme scheme )
{
	DEBUG_ENTRY( "AbortErrorMessage()" );

	if( scheme == AS_DATA_OPTIONAL )
		// presence is optional -> make warning less scary...
		fprintf( ioQQQ, "\nI could not open the data file %s\n\n", fname );
	else
		fprintf( ioQQQ, "\nPROBLEM DISASTER I could not open the data file %s\n\n", fname );
	if( cpu.firstOpen() || scheme == AS_DATA_ONLY )
	{
		// failed on very first open -> most likely path is not correct
		// failed on AS_DATA_ONLY -> CLOUDY_DATA_PATH may point to obsolete data dir
		fprintf( ioQQQ, "Although there may be other reasons you have received this error,\n");
		fprintf( ioQQQ, "the most likely are that the path has not been properly set\n");
		fprintf( ioQQQ, "or that the path points to an old version of the data.\n\n");
		fprintf( ioQQQ, "Please have a look at the file path.h in the source directory\n");
		fprintf( ioQQQ, "to check how the variable CLOUDY_DATA_PATH is set - \n");
		fprintf( ioQQQ, "it should give the location of the data files I need.\n");
		fprintf( ioQQQ, "These are the files in the data download from the web site.\n\n");
		fprintf( ioQQQ, "Recompile the code with the correct data path set in path.h\n");
		fprintf( ioQQQ, "or use the shell command \nexport CLOUDY_DATA_PATH=\"/path/to/data\"\n to set the\n");
		fprintf( ioQQQ, "path from a bash command prompt.\n\n");
		cpu.printDataPath();
	}
	else
	{
		// failed on search including local directory -> most likely the file name
		// was mistyped on a compile command, or Cloudy is run in the wrong directory
		// if scheme == AS_DATA_OPTIONAL, this most likely is a stellar grid that is not installed.
		fprintf( ioQQQ, "These are all the paths I tried:\n" );
		for( vector<string>::const_iterator ptr=PathList.begin(); ptr != PathList.end(); ++ptr )
			fprintf( ioQQQ, "   ==%s==\n", ptr->c_str() );
		// AS_DATA_OPTIONAL files should provide their own message (currently only stellar grids)
		if( scheme != AS_DATA_OPTIONAL )
		{
			fprintf( ioQQQ, "\nAlthough there may be other reasons you have received this error,\n");
			fprintf( ioQQQ, "the most likely are that you mistyped the file name, or that you\n");
			fprintf( ioQQQ, "are running Cloudy in the wrong directory. If you are running a\n");
			fprintf( ioQQQ, "COMPILE command, this needs to be done in the data directory.\n\n");
			fprintf( ioQQQ, "Otherwise, please have a look at the file path.h in the source\n");
			fprintf( ioQQQ, "directory to check how the variable CLOUDY_DATA_PATH is set - \n");
			fprintf( ioQQQ, "it should give the location of the data files I need.\n");
			fprintf( ioQQQ, "These are the files in the data download from the web site.\n\n");
			fprintf( ioQQQ, "Recompile the code with the correct data path set in path.h\n");
			fprintf( ioQQQ, "or use the shell command \nexport CLOUDY_DATA_PATH=\"/path/to/data\"\n to set the\n");
			fprintf( ioQQQ, "path from a bash command prompt.\n\n");
		}
	}
	fprintf(ioQQQ, "Sorry.\n\n\n");
	cdEXIT(EXIT_FAILURE);
}

FILE* open_data( const char* fname, const char* mode, access_scheme scheme )
{
	DEBUG_ENTRY( "open_data()" );

	bool lgAbort = ( scheme == AS_DATA_ONLY || scheme == AS_DATA_OPTIONAL || scheme == AS_DATA_LOCAL ||
			 scheme == AS_LOCAL_DATA || scheme == AS_LOCAL_ONLY );

	vector<string> PathList;
	cpu.getPathList( fname, PathList, scheme );

	FILE* handle = NULL;
	vector<string>::const_iterator ptr;
	for( ptr=PathList.begin(); ptr != PathList.end() && handle == NULL; ++ptr )
	{
		handle = fopen( ptr->c_str(), mode );
		if( trace.lgTrace )
			fprintf( ioQQQ, " open_data trying %s mode %s handle %p\n", ptr->c_str(), mode, handle );
	}

	if( handle == NULL && lgAbort )
		AbortErrorMessage( fname, PathList, scheme );

	++cpu.nFileDone;

	return handle;
}

void open_data( fstream& stream, const char* fname, ios_base::openmode mode, access_scheme scheme )
{
	DEBUG_ENTRY( "open_data()" );

	bool lgAbort = ( scheme == AS_DATA_ONLY || scheme == AS_DATA_OPTIONAL || scheme == AS_DATA_LOCAL ||
			 scheme == AS_LOCAL_DATA || scheme == AS_LOCAL_ONLY );

	vector<string> PathList;
	cpu.getPathList( fname, PathList, scheme );

	ASSERT( !stream.is_open() );
	vector<string>::const_iterator ptr;
	for( ptr=PathList.begin(); ptr != PathList.end() && !stream.is_open(); ++ptr )
	{
		stream.open( ptr->c_str(), mode );
		if( trace.lgTrace )
			fprintf( ioQQQ, " open_data trying %s succes? %c\n", ptr->c_str(), TorF(stream.is_open()) );
	}

	if( !stream.is_open() && lgAbort )
		AbortErrorMessage( fname, PathList, scheme );

	++cpu.nFileDone;
}

/** define routines for setting single and double precision signaling NaN
 * The bit pattern for an SNaN is implementation defined, but this should
 * work on most modern CPU's. The system definition is preferred, so in
 * C++ this should be replaced by numeric_limits<TYPE>::signaling_NaN() */

void set_NaN(sys_float &x)
{
	if( sizeof(sys_float) == 4 )
		*reinterpret_cast<int32*>(&x) = cpu.Float_SNaN_Value;
	else
		x = -FLT_MAX;
}

void set_NaN(sys_float x[], /* x[n] */
	     long n)
{
	long i;

	if( sizeof(sys_float) == 4 )
	{
		int32 *y = reinterpret_cast<int32*>(x);
		for( i=0; i < n; i++ )
			*y++ = cpu.Float_SNaN_Value;
	}
	else
	{
		for( i=0; i < n; i++ )
			x[i] = -FLT_MAX;
	}
}

void set_NaN(double &x)
{
	if( sizeof(double) == 8 )
	{
#		ifdef HAVE_INT64
		*reinterpret_cast<int64*>(&x) = cpu.Double_SNaN_Value;
#		else
		int32 *y = reinterpret_cast<int32*>(&x);
		*y++ = cpu.Double_SNaN_Value[0];
		*y = cpu.Double_SNaN_Value[1];
#		endif
	}
	else
		x = -DBL_MAX;
}

/* set_NaN - set NaN */
void set_NaN(double x[], /* x[n] */
	     long n)
{
	long i;

	if( sizeof(double) == 8 )
	{
#		ifdef HAVE_INT64
		int64 *y = reinterpret_cast<int64*>(x);
		for( i=0; i < n; i++ )
			*y++ = cpu.Double_SNaN_Value;
#		else
		int32 *y = reinterpret_cast<int32*>(x);
		for( i=0; i < n; i++ )
		{
			*y++ = cpu.Double_SNaN_Value[0];
			*y++ = cpu.Double_SNaN_Value[1];
		}
#		endif
	}
	else
	{
		for( i=0; i < n; i++ )
			x[i] = -DBL_MAX;
	}
}

/** detect quiet and signaling NaNs in single precision FP */
bool MyIsnan(const sys_float &x)
{
	if( sizeof(sys_float) == 4 && FLT_MAX_EXP-FLT_MIN_EXP+3 == 256 )
	{
		const int32 *p = reinterpret_cast<const int32*>(&x);
		int32 r = *p & 0x7f800000; r ^= 0x7f800000;
		int32 s = *p & 0x007fffff;
		return ( r == 0 && s != 0 );
	}
	else
		/* we don't understand this CPU */
		return false;
}

/** detect quiet and signaling NaNs in double precision FP */
bool MyIsnan(const double &x)
{
	if( sizeof(double) == 8 && DBL_MAX_EXP-DBL_MIN_EXP+3 == 2048 )
	{
#		ifdef HAVE_INT64
		const int64 *p = reinterpret_cast<const int64*>(&x);
		int64 r = *p & 0x7ff0000000000000; r ^= 0x7ff0000000000000;
		int64 s = *p & 0x000fffffffffffff;
		return ( r == 0 && s != 0 );
#		else
		const int32 *p = reinterpret_cast<const int32*>(&x);
		if( cpu.little_endian() )
		{
			int32 r = p[1] & 0x7ff00000; r ^= 0x7ff00000;
			int32 s = p[1] & 0x000fffff; s |= p[0];
			return ( r == 0 && s != 0 );
		}
		else if( cpu.big_endian() )
		{
			int32 r = p[0] & 0x7ff00000; r ^= 0x7ff00000;
			int32 s = p[0] & 0x000fffff; s |= p[1];
			return ( r == 0 && s != 0 );
		}
		else
			/* we don't understand this CPU */
			return false;
#		endif
	}
	else
		/* we don't understand this CPU */
		return false;
}
