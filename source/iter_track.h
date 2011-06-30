/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef ITER_TRACK_H_
#define ITER_TRACK_H_

////////////////////////////////////////////////////////////////////////////////
// The class iter_track was derived from this original source:                //
//          http://www.mymathlib.webtrellis.net/roots/amsterdam.html          //
//                                                                            //
// The original code was heavily modified by: Peter van Hoof, ROB             //
////////////////////////////////////////////////////////////////////////////////

double Amsterdam_Method( double (*f)(double), double a, double fa, double c, double fc,
			 double tol, int max_iter, int& err );

class iter_track
{
	vector< pair<double,double> > p_history;
	double p_result;
	double p_tol;
	int p_a;
	int p_b;
	int p_c;
	bool p_lgRootFound;

	void p_clear0()
	{
		p_history.clear();
	}
	void p_clear1()
	{
		p_history.reserve( 10 );
		set_NaN( p_result );
		p_tol = numeric_limits<double>::max();
		p_a = -1;
		p_b = -1;
		p_c = -1;
		p_lgRootFound = false;
	}
	void p_set_root(double x)
	{
		p_result = x;
		p_lgRootFound = true;
	}
	double p_x(int ip) const
	{
		return p_history[ip].first;
	}
	double p_y(int ip) const
	{
		return p_history[ip].second;
	}
	double p_midpoint() const
	{
		return 0.5*(p_x(p_a) + p_x(p_c));
	}
	double p_numerator(double dab, double dcb, double fa, double fb, double fc)
	{
		return fb*(dab*fc*(fc-fb)-fa*dcb*(fa-fb));
	}
	double p_denominator(double fa, double fb, double fc)
	{
		return (fc-fb)*(fa-fb)*(fa-fc);
	}

public:
	iter_track()
	{
		p_clear1();
	}
	~iter_track()
	{
		p_clear0();
	}
	void clear()
	{
		p_clear0();
		p_clear1();
	}
	void set_tol(double tol)
	{
		p_tol = tol;
	}
	double bracket_width() const
	{
		return p_x(p_c) - p_x(p_a);
	}
	bool lgConverged()
	{
		if( p_lgRootFound )
			return true;
		if( bracket_width() < p_tol )
		{
			p_result = p_midpoint();
			return true;
		}
		return false;
	}
	double root() const
	{
		return p_result;
	}
	int init_bracket( double x1, double fx1, double x2, double fx2 )
	{
		// fx1 and fx2 must have opposite sign, or be zero
		int s1 = sign3(fx1);
		int test = s1*sign3(fx2);
		if( test > 0 )
			return -1;
		if( test == 0 )
			p_set_root( ( s1 == 0 ) ? x1 : x2 );

		p_history.push_back( pair<double,double>(x1,fx1) );
		p_history.push_back( pair<double,double>(x2,fx2) );
		p_a = ( x1 < x2 ) ? 0 : 1;
		p_c = ( x1 < x2 ) ? 1 : 0;
		return 0;
	}
	void add( double x, double fx )
	{
		p_history.push_back( pair<double,double>(x,fx) );
		p_b = p_history.size()-1;
		if( fx == 0. )
			p_set_root( x );
	}
	double next_val();
	double next_val(double max_rel_step)
	{
		double next = next_val();
		double last = p_history.back().first;
		double rel_step = safe_div( next, last ) - 1.;
		rel_step = sign( min(abs(rel_step),abs(max_rel_step)), rel_step );
		return (1.+rel_step)*last;
	}
	// these routines return a numerical estimate of the derivative
	// by making a linear least squares fit of y(x) to the last n steps
	double deriv(int n, double& sigma) const;
	double deriv(double& sigma) const
	{
		return deriv( p_history.size(), sigma );
	}
	double deriv(int n) const
	{
		double sigma;
		return deriv( n, sigma );
	}
	double deriv() const
	{
		double sigma;
		return deriv( p_history.size(), sigma );
	}
	// these routines return a numerical estimate of the root
	// by making a linear least squares fit of x(y) to the last n steps
	double zero_fit(int n, double& sigma) const;
	double zero_fit(double& sigma) const
	{
		return zero_fit( p_history.size(), sigma );
	}
	double zero_fit(int n) const
	{
		double sigma;
		return zero_fit( n, sigma );
	}
	double zero_fit() const
	{
		double sigma;
		return zero_fit( p_history.size(), sigma );
	}
	// the following two methods are debugging aids
	void print_status() const
	{
		dprintf( ioQQQ, "a %i %.15e %.15e\n", p_a, p_x(p_a), p_y(p_a) );
		dprintf( ioQQQ, "b %i %.15e %.15e\n", p_b, p_x(p_b), p_y(p_b) );
		dprintf( ioQQQ, "c %i %.15e %.15e\n", p_c, p_x(p_c), p_y(p_c) );
	}
	void print_history() const
	{
		fprintf(ioQQQ," Temperature      C-H     history\n");
		for( unsigned int i=0; i < p_history.size(); ++i )
			fprintf( ioQQQ, "%.15e %.15e\n", p_x(i), p_y(i) );
	}
};

template<class T>
class iter_track_basic
{
	T p_lo_bound;
	T p_hi_bound;
	void p_clear1()
	{
		p_lo_bound = numeric_limits<T>::max();
		p_hi_bound = numeric_limits<T>::min();
	}
public:
	iter_track_basic()
	{
		p_clear1();
	}
	void clear()
	{
		p_clear1();
	}
	T next_val( T current, T next_est )
	{
		if( next_est < current )
			p_hi_bound = current;
		else
			p_lo_bound = current;
		if( p_lo_bound < p_hi_bound )
			return (T)0.5*(p_lo_bound + p_hi_bound);
		else
			return next_est;
	}
};

#endif
