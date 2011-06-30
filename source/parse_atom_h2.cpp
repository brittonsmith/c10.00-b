/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseAtomH2 parse information from the atom command line */
#include "cddefines.h" 
#include "hmi.h" 
#include "mole.h" 
#include "h2.h" 
#include "h2_priv.h" 
#include "parser.h" 
#include "thirdparty.h"

/*ParseAtomH2 parse information from the atom command line */
void ParseAtomH2(Parser &p )
{
	long int j;

	DEBUG_ENTRY( "ParseAtomH2()" );

	/* this command has a 2 in the H2 label - must not parse the two by
	 * accident.  Get the first number off the line image, and confirm that
	 * it is a 2 */
	j = (long int)p.FFmtRead();
	if( j != 2 )
	{
		fprintf( ioQQQ, " Something is wrong with the order of the numbers on this line.\n" );
		fprintf( ioQQQ, " The first number I encounter should be a 2.\n Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* the mere calling of this routine turns the large H2 molecule on */
	h2.lgH2ON = true;

	if( p.nMatch("LEVE") )
	{
		/* number of electronic levels */

		/* lgH2_READ_DATA is false at start of calculation, set true when 
		 * space allocated for the H lines.  Once done we must ignore all 
		 * future changes in the number of levels */
		if( !lgH2_READ_DATA )
		{
			mole.n_h2_elec_states = (long int)p.FFmtRead();
			if( p.lgEOL() )
			{
				if( p.nMatch("LARG") )
				{
					/* LARGE is option to use the most number of electronic levels */
					mole.n_h2_elec_states = N_H2_ELEC;
				}
				else
				{
					p.NoNumb("number of electronic levels");
				}
			}

			/* do not allow fewer than 3 - that includes Lyman & Werner bands */
			if( mole.n_h2_elec_states < 3 )
			{
				fprintf( ioQQQ, " This would be too few electronic levels - resetting to 3.\n" );
				mole.n_h2_elec_states = 3;
			}
			/* N_H2_ELEC is in h2.h and is the greatest number of elec lev possible */
			else if( mole.n_h2_elec_states > N_H2_ELEC )
			{
				fprintf( ioQQQ, 
					" This would be too many levels, the limit is %i.\n" , 
					N_H2_ELEC);
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	else if( p.nMatch("LIMI") )
	{
		/* the limit to the H2 / Htot ratio - 
		 * if smaller than this, do not compute large H2 mole */
		mole.H2_to_H_limit = p.FFmtRead();
		if( p.lgEOL() )
		{
			/* did not find a number, either mistake or key " off" */
			if( p.nMatch( " OFF"  ) )
			{
				/* turn off limit */
				mole.H2_to_H_limit = -1.;
			}
			else
			{
				p.NoNumb( "limit to the H2 / Htot ratio" );
			}
		}
		else
		{
			/* got a number, check if negative and so a log */
			/* a number <= 0 is the log of the ratio */
			if( mole.H2_to_H_limit <= 0. )
				mole.H2_to_H_limit = pow(10., mole.H2_to_H_limit);
		}
	}
	else if( p.nMatch("GBAR" ) )
	{
		/* option to either use, or not use, gbar approximation for low X 
		 * levels with no collision data - by default it is on */
		if( p.nMatch(" OFF" ) )
		{
			mole.lgColl_gbar = false;
		}
		else if( p.nMatch(" ON " ) )
		{
			mole.lgColl_gbar = true;
		}
		else
		{
			fprintf( ioQQQ, 
				" The gbar approximation must be off (\" OFF\") or on (\" ON \").\n");
			cdEXIT(EXIT_FAILURE);
		}
	}
	/* option to turn collisional effects off or on */
	else if( p.nMatch("COLL" ) )
	{
		/* option to turn collisional dissociation off or on */
		if( p.nMatch("DISS" ) )
		{
			/* option to turn collisions off */
			if( p.nMatch(" ON " ) )
			{
				/* this is the default, leave collisions off */
				mole.lgColl_dissoc_coll = true;
			}
			else
			{
				/* default (and only reason for this command) is to turn off collisions */
				mole.lgColl_dissoc_coll = false;
			}
		}
		/* option to turn collisional dissociation off or on 
		 * >>chng 06 mar 01, had been simply if - so all collisions were turned off
		 * when dissociation collisions turned off - 
		 * due to bucket else at end */
		else if( p.nMatch("ORTH" ) && p.nMatch("PARA" ) )
		{
			/* option to turn ortho - para collisions with particles off */
			if( p.nMatch(" ON " ) )
			{
				/* this is the default, leave collisions off */
				mole.lgH2_ortho_para_coll_on = true;
			}
			else
			{
				/* default (and only reason for this command) is to turn off 
				 * ortho-para collisions */
				mole.lgH2_ortho_para_coll_on = false;
			}
		}

		/* option to turn collisional effects off or on */
		else if( p.nMatch("GRAI" ) )
		{
			/* option to turn collisions off */
			if( p.nMatch(" ON" ) )
			{
				/* this is the default, leave collisions off */
				mole.lgH2_grain_deexcitation = true;
			}
			else
			{
				/* default (and only reason for this command) is to turn off collisions */
				mole.lgH2_grain_deexcitation = false;
			}
		}
		else if( p.nMatch(" HE " ) )
		{
			/* atom H2 He collisions ORNL (the default), Le BOURlot, and OFF
			 * which data set for He collisions,
			 * Teck Lee et al. ApJ to be submitted */
			if( p.nMatch(" NEW" ) || p.nMatch("ORNL" ) )
			{
				/* use the new coefficients */
				mole.lgH2_He_ORNL = true;
			}
			else if( p.nMatch(" OLD" ) || p.nMatch("BOUR" ) )
			{
				/* use the coefficients from
				 *>>refer	H2	collision	Le Bourlot, J., Pineau des Forets, 
				 *>>refercon	G., & Flower, D.R. 1999, MNRAS, 305, 802*/
				mole.lgH2_He_ORNL = false;
			}
			else
			{
				fprintf( ioQQQ, 
					" I did not find a keyword on this ATOM H2 HE command - I know about the keys ORNL and Le BOURlot\n");
				cdEXIT(EXIT_FAILURE);
			}
		}

		/*>>chng 08 feb 27, GS*/
		else if( p.nMatch("ORH2" ) )
		{
			/* atom H2 H2ortho collisions ORNL (the default), Le BOURlot, and OFF
			 * which data set for H2 collisions,
			 * Teck Lee et al. ApJ to be submitted */
			if( p.nMatch("ORNL" ) )
			{
				/* use the new coefficients */
				mole.lgH2_ORH2_ORNL = true;
			}
			else if( p.nMatch("BOUR" ) )
			{
				/* use the coefficients from
				 *>>refer	H2	collision	Le Bourlot, J., Pineau des Forets, 
				 *>>refercon	G., & Flower, D.R. 1999, MNRAS, 305, 802*/
				mole.lgH2_ORH2_ORNL = false;
			}
			else
			{
				fprintf( ioQQQ, 
					" I did not find a keyword on this ATOM H2 ohH2 command - I know about the keys ORNL and Le BOURlot\n");
				cdEXIT(EXIT_FAILURE);
			}
		}

		else if( p.nMatch("PAH2" ) )
		{
			/* atom H2 H2ortho collisions ORNL (the default), Le BOURlot, and OFF
			 * which data set for He collisions,
			 * Teck Lee et al. ApJ to be submitted */
			if(  p.nMatch("ORNL" ) )
			{
				/* use the new coefficients */
				mole.lgH2_PAH2_ORNL = true;
			}
			else if( p.nMatch("BOUR" ) )
			{
				/* use the coefficients from
				 *>>refer	H2	collision	Le Bourlot, J., Pineau des Forets, 
				 *>>refercon	G., & Flower, D.R. 1999, MNRAS, 305, 802*/
				mole.lgH2_PAH2_ORNL = false;
			}
			else
			{
				fprintf( ioQQQ, 
					" I did not find a keyword on this ATOM H2 paH2 command - I know about the keys ORNL and Le BOURlot\n");
				cdEXIT(EXIT_FAILURE);
			}
		}

		else
		{
			/* option to turn all collisions off */
			if( p.nMatch(" ON " ) )
			{
				/* this is the default, leave collisions on */
				mole.lgColl_deexec_Calc = true;
			}
			else
			{
				/* default (and only reason for this command) is to turn off collisions */
				mole.lgColl_deexec_Calc = false;
			}
		}
	}

	/* set number of levels in matrix, but not trace matrix option */
	else if( p.nMatch("MATR" ) && !p.nMatch("TRAC" ) )
	{
		/* matrix option sets the number of levels that will
		 * be included in the matrix solution */
		nXLevelsMatrix = (long)p.FFmtRead();
		if( p.nMatch(" ALL") )
		{
			/* " all" means do all of X, but space has not yet been allocated,
			 * so we do know know how many levels are within X - set special
			 * flag that will be used then this is known */
			/*nXLevelsMatrix = nLevels_per_elec[0];*/
			nXLevelsMatrix = -1;
		}
		else if( p.lgEOL() && !(p.nMatch(" OFF") || p.nMatch("NONE") ) )
		{
			/* this branch hit eol but OFF or NONE is not on line - this is a mistake */
			fprintf( ioQQQ, 
				" The total number of levels used in the matrix solver must be entered, or keywords ALL or NONE entered.\n Sorry.\n");
			cdEXIT(EXIT_FAILURE);
		}
		/* cannot check less than total number of levels within x since not yet set 
		 * We do not certify that matrix limits are greater than 1 -
		 * zero or <0 limits just turns if off, as did the off option */
	}
	else if( p.nMatch(" LTE" ) )
	{
		/* LTE option causes code to assume LTE for level populations  */
		mole.lgH2_LTE = true;
	}

	else if( p.nMatch("TRAC" ) )
	{
		/* these are used to set trace levels of output 
		mole.nH2_trace_final = 1;
		mole.nH2_trace_iterations = 2;
		mole.nH2_trace_full = 3;
		mole.nH2_trace_matrix = 4*/

		/* turns on trace printout - there are multiple levels */
		if( p.nMatch("FINA" ) )
		{
			/* FINAL gives only final information when solver exits */
			mole.nH2_TRACE = mole.nH2_trace_final;
		}
		else if( p.nMatch("ITER" ) )
		{
			/* follow iterations within each call */
			mole.nH2_TRACE = mole.nH2_trace_iterations;
		}
		else if( p.nMatch("FULL" ) )
		{
			/* full details of solution - this is also the default*/
			mole.nH2_TRACE = mole.nH2_trace_full;
		}
		else if( p.nMatch("MATR" ) )
		{
			/* print the matrices used for X */
			mole.nH2_TRACE = mole.nH2_trace_matrix;
		}
		else
		{
			/* full details of solution is also the default*/
			mole.nH2_TRACE = mole.nH2_trace_full;
		}
	}
	else if( p.nMatch("NOIS" ) )
	{
		unsigned int iseed;
		/* check on effects of uncertainties in collision rates */
		mole.lgH2_NOISE = true;
		mole.lgH2_NOISECOSMIC = true;

		/* optional mean - default is 0 */
		mole.xMeanNoise = p.FFmtRead();
		if( p.lgEOL() )
			mole.xMeanNoise = 0.;

		/* this is the standard deviation for the mole, with default */
		mole.xSTDNoise = p.FFmtRead();
		if( p.lgEOL() )
			mole.xSTDNoise = 0.5;

		/* this may be a seed for the random number generator.  if no seed is
		 * set then use system time, and always get different sequence */
		iseed = (unsigned int)p.FFmtRead();
		/* returned 0 if eol hit */
		if( iseed > 0 )
		{
			/* user set seed */
			init_genrand( iseed );
		}
		else
		{
			init_genrand( (unsigned)time( NULL ) );
		}
	}

	else if( p.nMatch("THER" ) )
	{
		/* change the treatment of the heating - cooling effects of H2,
		 * options are simple (use TH85 expressions) and full (use large molecule)*/
		if( p.nMatch("SIMP" ) )
		{
			hmi.lgH2_Thermal_BigH2 = false;
		}
		else if( p.nMatch("FULL" ) )
		{
			/* this is the default - use big atom */
			hmi.lgH2_Thermal_BigH2 = true;
		}
	}

	else if( p.nMatch("CHEM" ) )
	{
		/* atom h2 chemistry simple command
		 * change the treatment of the chemistry - formation and destruction,
		 * options are simple (use TH85 expressions) and full (use large molecule)*/
		if( p.nMatch("SIMP" ) )
		{
			hmi.lgH2_Chemistry_BigH2 = false;
		}
		else if( p.nMatch("FULL" ) )
		{
			/* this is the default - use big atom */
			hmi.lgH2_Chemistry_BigH2 = true;
		}
	}

	/* there is no final branch - if we do not find a keyword, simply
	 * turn on the H2 molecule */
	return;
}
