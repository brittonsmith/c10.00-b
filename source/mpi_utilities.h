/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef MPI_UTILITIES_H_
#define MPI_UTILITIES_H_

#ifdef MPI_ENABLED

#include <mpi.h>

namespace MPI
{
	//!
	//! Create a set of overloaded functions that allows us to pull in MPI
	//! data types in a safer and more C++ way.
	//!
	//! Some typical ways to use it are:
	//!   double x;
	//!   MPI::COMM_WORLD.Bcast( &x, 1, MPI::type(x), 0 );
	//!
	//! But also:
	//!   typedef unsigned int myType;
	//!   myType t;
	//!   MPI::COMM_WORLD.Bcast( &t, 1, MPI::type(t), 0 );
	//!
	//! or:
	//!   realnum v[20];
	//!   MPI::COMM_WORLD.Bcast( v, 20, MPI::type(v), 0 );
	//!
	//! This is very useful for realnum: MPI::type(v) will do the right thing,
	//! no matter what the setting of FLT_IS_DBL is. Also note that MPI::type
	//! can take both variables and arrays as its argument.
	//!
	//! The list below should cover all POD types that are currently used in
	//! Cloudy.
	//!
	//! PS - the "t" in MPI::type is deliberately chosen to be lower case
	//!      so that it cannot possibly clash with existing MPI symbols.
	//!

	inline const Datatype& type(bool) { return BOOL; }
	inline const Datatype& type(const bool*) { return BOOL; }
	inline const Datatype& type(char) { return CHAR; }
	inline const Datatype& type(const char*) { return CHAR; }
	inline const Datatype& type(unsigned char) { return UNSIGNED_CHAR; }
	inline const Datatype& type(const unsigned char*) { return UNSIGNED_CHAR; }
	inline const Datatype& type(short int) { return SHORT; }
	inline const Datatype& type(const short int*) { return SHORT; }
	inline const Datatype& type(unsigned short int) { return UNSIGNED_SHORT; }
	inline const Datatype& type(const unsigned short int*) { return UNSIGNED_SHORT; }
	inline const Datatype& type(int) { return INT; }
	inline const Datatype& type(const int*) { return INT; }
	inline const Datatype& type(unsigned int) { return UNSIGNED; }
	inline const Datatype& type(const unsigned int*) { return UNSIGNED; }
	inline const Datatype& type(long) { return LONG_INT; }
	inline const Datatype& type(const long*) { return LONG_INT; }
	inline const Datatype& type(unsigned long) { return UNSIGNED_LONG; }
	inline const Datatype& type(const unsigned long*) { return UNSIGNED_LONG; }
	inline const Datatype& type(sys_float) { return FLOAT; }
	inline const Datatype& type(const sys_float*) { return FLOAT; }
	inline const Datatype& type(double) { return DOUBLE; }
	inline const Datatype& type(const double*) { return DOUBLE; }
	inline const Datatype& type(complex<sys_float>) { return COMPLEX; }
	inline const Datatype& type(const complex<sys_float>*) { return COMPLEX; }
	inline const Datatype& type(complex<double>) { return DOUBLE_COMPLEX; }
	inline const Datatype& type(const complex<double>*) { return DOUBLE_COMPLEX; }

	class load_balance
	{
		vector<int> p_jobs;
		unsigned int p_ptr;
		void p_clear0()
		{
			p_jobs.clear();
		}
		void p_clear1()
		{
			p_ptr = 0;
		}
	public:
		load_balance()
		{
			p_clear1();
		}
		explicit load_balance( int nJobs )
		{
			p_clear1();
			init( nJobs );
		}
		~load_balance()
		{
			p_clear0();
		}
		void clear()
		{
			p_clear0();
			p_clear1();
		}
		void init( int nJobs );
		int next_job()
		{
			if( p_ptr < p_jobs.size() )
			{
				int res = p_jobs[p_ptr];
				p_ptr += COMM_WORLD.Get_size();
				return res;
			}
			else
				return -1;
		}
		void finalize()
		{
			// wait for all jobs to finish
			COMM_WORLD.Barrier();
		}
	};
}

#else /* MPI_ENABLED */

namespace MPI
{
	// This global struct is needed so that we can #define away the arguments of
	// calls to MPI routines, which allows us to reduce the number of stubs needed.
	// Since it contains no real data and only an inline function, the fact that
	// it is global creates no problems (it stores nothing in memory). Some compilers
	// (like g++) don't even require this struct to be allocated.
	EXTERN struct t_MPI
	{
		int total_insanity() { return TotalInsanityAsStub<int>(); }
	} COMM_WORLD;
}

// define MPI stubs here, so that we don't get endless #ifdef MPI_ENBLED in the code...
// this way we can use if( cpu.lgMPI() ) { .... } instead
#define Barrier() total_insanity()
#define Bcast(W,X,Y,Z) total_insanity()
#define Reduce(U,V,W,X,Y,Z) total_insanity()

#endif /* MPI_UTILITIES_H_ */

/** GridPointPrefix: generate filename prefix for any files associated with a single point in a grid */
inline string GridPointPrefix(int n)
{
	ostringstream oss;
	oss << "grid" << setfill( '0' ) << setw(9) << n << "_";
	return oss.str();
}

/** process_output: concatenate output files produced in MPI grid run */
void process_output();

/** append_file: append output produced on file <source> to open file descriptor <dest> */
void append_file( FILE*, const char* );

#endif /* _MPI_UTILITIES_H_ */
