/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ion_recom_calculate calculate radiative and dielectronic recombination rate coefficients */
/*Badnell_rec_init This code is written by Terry Yun, 2005 *
 * It reads rate coefficient fits into 3D arrays and output array.out for testing *
 * The testing can be commented out */
/*Badnell_DR_rate_eval This code is written by Terry Yun, 2005 *
 * It interpolates the rate coefficients in a given temperature.*
   It receives ATOMIC_NUM_BIG, NELECTRONS values, temperature and returns the rate coefficient*
   It returns
        '-2': initial <= final
              init < 0 or init >302 or final < 0 or final > 302
        '-1': the transition is not defined
        '99': unknown invalid entries */ 
#include "cddefines.h"
#include "phycon.h"
#include "elementnames.h"
#include "atmdat.h"
#include "iso.h"
#include "ionbal.h"
#include "dense.h"
#include "taulines.h"

static const int MAX_FIT_PAR_DR = 9;
static double ***DRFitParPart1; 
static double ***DRFitParPart2;
static int **nDRFitPar;

static const int MAX_FIT_PAR_RR = 6;
static double ***RRFitPar; 
static long int *nsumrec;

/* flags to recall that we have read the fits from the main data files */
static bool **lgDRBadnellDefined ,
	**lgDRBadnellDefinedPart2,
	**lgRRBadnellDefined;
static bool lgMustMallocRec=true;

/* these enable certain debugging print statements */
/* #define PRINT_DR */
/* #define PRINT_RR */

#if defined(PRINT_DR) || defined(PRINT_RR)
static const char FILE_NAME_OUT[] = "array.out";
#endif

/**\verbatim Badnell_DR_rate_eval This code is written by Terry Yun, 2005 
It interpolates the rate coefficients in a given temperature.
It receives atomic number on Physics scale, with H = 1, 
and the number of core electrons before recombination, and returns the rate coefficient*
It returns
'-2': initial <= final
init < 0 or init >302 or final < 0 or final > 302
'-1': the transition is not defined
'99': unknown invalid entries                         
\endverbatim
\param z_val atomic number on C scale - He is 1
\param n_val number of core electrons before capture of free electron 
*/ 
STATIC double Badnell_DR_rate_eval(
	/* atomic number on C scale - He - 1 */
	int nAtomicNumberCScale, 
	/* number of core electrons before capture of free electron */
	int n_core_e_before_recomb )
{

	double RateCoefficient, sum;
	int i;

	DEBUG_ENTRY( "Badnell_DR_rate_eval()" );
	ASSERT( nAtomicNumberCScale>=0 && nAtomicNumberCScale<LIMELM );

	if( nAtomicNumberCScale==ipIRON && n_core_e_before_recomb>=12 && 
		n_core_e_before_recomb<=18 )
	{
		/* these data are from 
		*>>refer	Fe	DR	Badnell, N., 2006, ApJ, 651, L73
		* Fe 8+ to Fe 12+, but also include Fe13+ and Fe 14+,
		* so these are 26-8=18 to 26-14=12 
		* increasing number of bound electrons, 0 is 14 elec, 1 is 15 elec 
		* Fe 3p^q, q=2-6
		* this is DR fit coefficients given in table 1 of Badnell 06 */
		double cFe_q[7][8] =
		{
			{5.636e-4, 7.390e-3, 3.635e-2, 1.693e-1, 3.315e-2, 2.288e-1, 7.316e-2, 0.},
			{1.090e-3, 7.801e-3, 1.132e-2, 4.740e-2, 1.990e-1, 3.379e-2, 1.140e-1, 1.250e-1},
			{3.266e-3, 7.637e-3, 1.005e-2, 2.527e-2, 6.389e-2, 1.564e-1, 0.,       0.},
			{1.074e-3, 6.080e-3, 1.887e-2, 2.540e-2, 7.580e-2, 2.773e-1, 0.,       0.},
			{9.073e-4, 3.777e-3, 1.027e-2, 3.321e-2, 8.529e-2, 2.778e-1, 0.,       0.},
			{5.335e-4, 1.827e-3, 4.851e-3, 2.710e-2, 8.226e-2, 3.147e-1, 0.,       0.},
			{7.421e-4, 2.526e-3, 4.605e-3, 1.489e-2, 5.891e-2, 2.318e-1, 0.,       0.}
		};

		/* Table 2 of Badnell 06 */
		double EFe_q[7][8] =
		{
			{3.628e3, 2.432e4, 1.226e5, 4.351e5, 1.411e6, 6.589e6, 1.030e7, 0},
			{1.246e3, 1.063e4, 4.719e4, 1.952e5, 5.637e5, 2.248e6, 7.202e6, 3.999e9},
			{1.242e3, 1.001e4, 4.466e4, 1.497e5, 3.919e5, 6.853e5, 0.     , 0.},
			{1.387e3, 1.048e4, 3.955e4, 1.461e5, 4.010e5, 7.208e5, 0.     , 0.},
			{1.525e3, 1.071e4, 4.033e4, 1.564e5, 4.196e5, 7.580e5, 0.     , 0.},
			{2.032e3, 1.018e4, 4.638e4, 1.698e5, 4.499e5, 7.880e5, 0.     , 0.},
			{3.468e3, 1.353e4, 3.690e4, 1.957e5, 4.630e5, 8.202e5, 0.     , 0.}
		};
		/* nion is for the above block of numbers */
		long int nion = n_core_e_before_recomb - 12;
		ASSERT( nion>=0 && nion <=6 );

		sum = 0;
		i = 0;
		/* loop over all non-zero terms */
		for(i=0; i<8; ++i  )
		{
			sum += (cFe_q[nion][i] * sexp( EFe_q[nion][i]/phycon.te));
		}

		/*RateCoefficient = pow(phycon.te, -1.5) * sum;*/
		RateCoefficient = sum / phycon.te32;

		return RateCoefficient;
	}

	/*Invalid entries returns '-2':more electrons than protons */
	else if( nAtomicNumberCScale < n_core_e_before_recomb )
	{
		RateCoefficient = -2;
	}
	/*Invalid entries returns '-2' if nAtomicNumberCScale and n_core_e_before_recomb are out of the range*/
	else if( nAtomicNumberCScale >= LIMELM )
	{
		RateCoefficient = -2;
	}
	/*undefined z and n returns '-1'*/
	else if( !lgDRBadnellDefined[nAtomicNumberCScale][n_core_e_before_recomb] )
	{
		RateCoefficient = -1;
	}
	else if( lgDRBadnellDefined[nAtomicNumberCScale][n_core_e_before_recomb] )
	{
		/* this branch, recombination coefficient has been defined */
		sum = 0;
		i = 0;
		/* loop over all non-zero terms */
		for(i=0; i<nDRFitPar[nAtomicNumberCScale][n_core_e_before_recomb]; ++i  )
		{
			sum += (DRFitParPart1[nAtomicNumberCScale][n_core_e_before_recomb][i] *
				sexp( DRFitParPart2[nAtomicNumberCScale][n_core_e_before_recomb][i]/phycon.te));
		}

		/*RateCoefficient = pow(phycon.te, -1.5) * sum;*/
		RateCoefficient = sum / phycon.te32;
	}
	/*unknown invalid entries returns '-99'*/
	else
	{
		RateCoefficient = -99;
	}

	return RateCoefficient;
}

/**Badnell_RR_rate_eval
\param z_val atomic number on C scale - He - 1
\param n_val number of core electrons before capture of free electron
*/
STATIC double Badnell_RR_rate_eval(
			/* atomic number on C scale - He - 1 */
			int nAtomicNumberCScale, 
			/* number of core electrons before capture of free electron */
			int n_core_e_before_recomb )
{
	double RateCoefficient;
	double B, D, F;

	DEBUG_ENTRY( "Badnell_RR_rate_eval()" );

	ASSERT( nAtomicNumberCScale>=0 && nAtomicNumberCScale<LIMELM );

	if( nAtomicNumberCScale==ipIRON && 
		n_core_e_before_recomb>=12 && n_core_e_before_recomb<=18 )
	{
		/* RR rate coefficients from Table 3 of
		*>>refer	Fe	RR	Badnell, N. 2006, ApJ, 651, L73 
		* Fe 8+ to Fe 12+, but also include Fe13+ and Fe 14+,
		* so these are 26-8=18 to 26-14=12 
		* increasing number of bound electrons, 0 is 14 elec, 1 is 15 elec 
		* Fe 3p^q, q=2-6
		* this is DR fit coefficients given in table 1 of Badnell 06 */
		double parFeq[7][6] ={
			{1.179e-9 , 0.7096, 4.508e2, 3.393e7, 0.0154, 3.977e6},
			{1.050e-9 , 0.6939, 4.568e2, 3.987e7, 0.0066, 5.451e5},
			{9.832e-10, 0.7146, 3.597e2, 3.808e7, 0.0045, 3.952e5},
			{8.303e-10, 0.7156, 3.531e2, 3.554e7, 0.0132, 2.951e5},
			{1.052e-9 , 0.7370, 1.639e2, 2.924e7, 0.0224, 4.291e5},
			{1.338e-9 , 0.7495, 7.242e1, 2.453e7, 0.0404, 4.199e5},
			{1.263e-9 , 0.7532, 5.209e1, 2.169e7, 0.0421, 2.917e5}
		};

		double temp;
		/* nion is for the above block of numbers */
		long int nion = n_core_e_before_recomb - 12;
		ASSERT( nion>=0 && nion <=6 );

		temp = -parFeq[nion][5]/phycon.te; /* temp = (-T2/T) */
		B = parFeq[nion][1] + parFeq[nion][4]*exp(temp);
		D = sqrt(phycon.te/parFeq[nion][2]); /* D = (T/T0)^1/2 */
		F = sqrt(phycon.te/parFeq[nion][3]); /* F = (T/T1)^1/2 */
		RateCoefficient = parFeq[nion][0]/(D*pow((1.+D),(1.-B))*pow((1.+F),(1.+B)));

		return RateCoefficient;
	}

	/*Invalid entries returns '-2':if the z_values are smaller than equal to the n_values */
	else if( nAtomicNumberCScale < n_core_e_before_recomb )
	{
		RateCoefficient = -2;
	}
	/*Invalid entries returns '-2' if nAtomicNumberCScale and n_core_e_before_recomb are out of the range*/
	else if( nAtomicNumberCScale >= LIMELM )
	{
		RateCoefficient = -2;
	}
	/*undefined z and n returns '-1'*/
	else if( !lgRRBadnellDefined[nAtomicNumberCScale][n_core_e_before_recomb] )
	{
		RateCoefficient = -1;
	}
	/* coefficients:A=RRFitPar[0], B=RRFitPar[1], T0=RRFitPar[2], T1=RRFitPar[3], DRFitParPart1=RRFitPar[4], T2=RRFitPar[5] */
	else if( lgRRBadnellDefined[nAtomicNumberCScale][n_core_e_before_recomb] )
	{

		/* RateCoefficient=A*[(T/T0)^1/2*(1+(T/T0)^1/2)^1-B*(1+(T/T1)^1/2)^1+B]^-1 
		where B = B + DRFitParPart1*exp(-T2/T) */
		double temp;

		temp = -RRFitPar[nAtomicNumberCScale][n_core_e_before_recomb][5]/phycon.te; /* temp = (-T2/T) */
		B = RRFitPar[nAtomicNumberCScale][n_core_e_before_recomb][1] + 
			RRFitPar[nAtomicNumberCScale][n_core_e_before_recomb][4]*exp(temp);
		D = sqrt(phycon.te/RRFitPar[nAtomicNumberCScale][n_core_e_before_recomb][2]); /* D = (T/T0)^1/2 */
		F = sqrt(phycon.te/RRFitPar[nAtomicNumberCScale][n_core_e_before_recomb][3]); /* F = (T/T1)^1/2 */
		RateCoefficient = RRFitPar[nAtomicNumberCScale][n_core_e_before_recomb][0]/(D*pow((1.+D),(1.-B))*pow((1.+F),(1.+B)));
	}

	/*unknown invalid entries returns '-99'*/
	else
		RateCoefficient = -99;

	return RateCoefficient;
}


/*Badnell_rec_init This code is written by Terry Yun, 2005 *
 * It reads rate coefficient fits into 3D arrays and output array.out for testing *
 * The testing can be commented out */
void Badnell_rec_init( void )
{

	double par_C[MAX_FIT_PAR_DR];
	double par_E[MAX_FIT_PAR_DR];
	char chLine[INPUT_LINE_LENGTH];
	int NuclearCharge=-1, NumberElectrons=-1;
	int i, k;
	int count, number;
	double temp_par[MAX_FIT_PAR_RR];
	int M_state, W_state,
		nelem;

	const int NBLOCK = 2;
	int data_begin_line[NBLOCK];/*it tells you where the data set begins(begins with 'Z')*/
	int length_of_line;	/*this variable checks for a blank line*/
	FILE *ioDATA;
	const char* chFilename;
	int yr, mo, dy;
	char *chs;

	const int BIGGEST_INDEX_TO_USE = 103;

	/* Declaration of data file name array - done by Kausalya */
	long TheirIndexToOurIndex[BIGGEST_INDEX_TO_USE];
	char string[120];
	double value;
	bool lgEOL;
	long int i1, ipISO = ipHE_LIKE;
	long INDX=0,INDP=0,N=0,S=0,L=0,J=0,maxINDX=0,loopindex=0,max_N_of_data=-1;
	bool lgFlag = true;

	static int nCalled = 0;

	const char* cdDATAFILE[] = 
	{
		/* the list of filenames for Badnell DR, one to two electron */
		"",
		"nrb00#h_he1ic12.dat", 
		"nrb00#h_li2ic12.dat",
		"nrb00#h_be3ic12.dat", 
		"nrb00#h_b4ic12.dat", 
		"nrb00#h_c5ic12.dat", 
		"nrb00#h_n6ic12.dat", 
		"nrb00#h_o7ic12.dat", 
		"nrb00#h_f8ic12.dat", 
		"nrb00#h_ne9ic12.dat", 
		"nrb00#h_na10ic12.dat", 
		"nrb00#h_mg11ic12.dat", 
		"nrb00#h_al12ic12.dat", 
		"nrb00#h_si13ic12.dat", 
		"nrb00#h_p14ic12.dat", 
		"nrb00#h_s15ic12.dat", 
		"nrb00#h_cl16ic12.dat", 
		"nrb00#h_ar17ic12.dat", 
		"nrb00#h_k18ic12.dat", 
		"nrb00#h_ca19ic12.dat", 
		"nrb00#h_sc20ic12.dat", 
		"nrb00#h_ti21ic12.dat", 
		"nrb00#h_v22ic12.dat", 
		"nrb00#h_cr23ic12.dat", 
		"nrb00#h_mn24ic12.dat", 
		"nrb00#h_fe25ic12.dat",
		"nrb00#h_co26ic12.dat", 
		"nrb00#h_ni27ic12.dat", 
		"nrb00#h_cu28ic12.dat", 
		"nrb00#h_zn29ic12.dat"
	};
	//End of modification

	DEBUG_ENTRY( "Badnell_rec_init()" );

	/* must only do this once */
	if( nCalled > 0 )
	{
		return;
	}
	++nCalled;

#	if defined(PRINT_DR) || defined(PRINT_RR)
	FILE *ofp = open_data( FILE_NAME_OUT, "w", AS_LOCAL_ONLY );
#	endif

	/* Modification  done by Kausalya 
	 * Start - Try to open all the 29 data files.*/
	for(nelem=ipHELIUM;nelem<LIMELM;nelem++)
	{
		if( nelem < 2 || dense.lgElmtOn[nelem] )
		{
			ioDATA= open_data( cdDATAFILE[nelem], "r" );

			lgFlag = true;
			ASSERT(ioDATA);

			for( i=0; i<BIGGEST_INDEX_TO_USE; i++ )
				TheirIndexToOurIndex[i] = -1;

			/* Reading lines */
			while(lgFlag)
			{
				if(read_whole_line(string,sizeof(string),ioDATA)!=NULL)
				{
					if( nMatch("INDX  INDP ",string) )
					{
						/* ignore next line of data */
						if( read_whole_line( string , (int)sizeof(string) , ioDATA ) == NULL )
						{
							fprintf( ioQQQ, " Badnell data file appears to be corrupted.\n");
							cdEXIT(EXIT_FAILURE);
						}

						/* This one should be real data */
						while( read_whole_line(string, (int)sizeof(string), ioDATA) != NULL )
						{
							if( strcmp(string,"\n")==0 )
							{
								lgFlag = false;
								break;
							}

							i1=3;
							INDX=(long)FFmtRead(string,&i1,INPUT_LINE_LENGTH,&lgEOL);
							if( INDX >= BIGGEST_INDEX_TO_USE )
							{
								INDX--;
								lgFlag = false;
								break;
							}

							ASSERT( INDX < BIGGEST_INDEX_TO_USE );									 

							INDP=(long)FFmtRead(string,&i1,INPUT_LINE_LENGTH,&lgEOL);
							ASSERT( INDP >= 1 );									 

							if(INDP==1)
							{
								if( (i1=nMatch("1S1 ",string)) > 0 )
								{
									i1 += 4;
									N=(long)FFmtRead(string,&i1,INPUT_LINE_LENGTH,&lgEOL);
									ASSERT( N>=1 );
								}
								else
								{
									TotalInsanity();
								}

								if( (i1=nMatch("     (",string)) > 0 )
								{
									i1 += 6;
									S=(long)FFmtRead(string,&i1,INPUT_LINE_LENGTH,&lgEOL);
									/* S in file is 3 or 1, we need 1 or 0 */
									ASSERT( S==1 || S==3 );
								}
								else
								{
									TotalInsanity();
								}

								/* move i1 one further to get L */
								i1++;
								L=(long)FFmtRead(string,&i1,INPUT_LINE_LENGTH,&lgEOL);
								ASSERT( L >= 0 && L < N );

								/* move i1 two further to get J */
								i1 += 2;
								J=(long)FFmtRead(string,&i1,INPUT_LINE_LENGTH,&lgEOL);
								ASSERT( J <= ( L + (int)((S+1)/2) ) && 
									J >= ( L - (int)((S+1)/2) ) && J >= 0 );

								/* if line in data file is higher N than highest considered, stop reading.  */
								if( N<= iso.n_HighestResolved_max[ipHE_LIKE][nelem] + iso.nCollapsed_max[ipHE_LIKE][nelem] )
									TheirIndexToOurIndex[INDX] = iso.QuantumNumbers2Index[ipHE_LIKE][nelem][N][L][S];
								else
								{
									/* Current line is not being used, 
									 * decrement INDX so maxINDX is set correctly below. */
									INDX--;
									lgFlag = false;
									break;
								}

								/* Must adjust index if in 2^3Pj term */
								if( N==2 && L==1 && S==3 )
								{
									if( J==0 )
										TheirIndexToOurIndex[INDX] = 3;
									else if( J==1 )
										TheirIndexToOurIndex[INDX] = 4;
									else
									{
										ASSERT( J==2 );
										ASSERT( TheirIndexToOurIndex[INDX] == 5 );
									}
								}
								max_N_of_data = MAX2( max_N_of_data, N );
							}    
							else                                                                                      
							{
								// Stop parsing the tuple since INDP!=1
								lgFlag = false;   
							}     
						}                              
					}                                                                           
				}
				else
				{    
					// End of file is reached.
					lgFlag = false;
				}   
			}

			maxINDX =INDX;
			ASSERT( maxINDX > 0 );
			ASSERT( maxINDX < BIGGEST_INDEX_TO_USE );
			/* reset INDX */
			INDX = 0;
			lgFlag = true;
			while(lgFlag)
			{
				if(read_whole_line(string,sizeof(string),ioDATA)!=NULL)
				{
					/* to access the first table whose columns are INDX ,INDP */
					if( nMatch("INDX TE= ",string) )
					{
						lgFlag = false;
						/* we found the beginning of the data array */
						/* ignore next line of data */
						if( read_whole_line( string , (int)sizeof(string) , ioDATA ) == NULL )
						{
							fprintf( ioQQQ, " Badnell data file appears to be corrupted.\n");
							cdEXIT(EXIT_FAILURE);
						}

						/* This one should be real data */
						while( read_whole_line(string, (int)sizeof(string), ioDATA) != NULL )
						{
							/* If we find this string, we have reached the end of the table. */
							if( nMatch("PRTF",string) || INDX >= maxINDX || INDX<0 )
								break;

							i1=3;
							INDX=(long)FFmtRead(string,&i1,INPUT_LINE_LENGTH,&lgEOL);
							if( INDX>maxINDX )
								break;

							for(loopindex=0;loopindex<10;loopindex++)
							{
								value=(double)FFmtRead(string,&i1,INPUT_LINE_LENGTH,&lgEOL);
								if( TheirIndexToOurIndex[INDX] < iso.numLevels_max[ipHE_LIKE][nelem] && 
									TheirIndexToOurIndex[INDX] > 0 )
								{
									iso.DielecRecombVsTemp[ipHE_LIKE][nelem][TheirIndexToOurIndex[INDX]][loopindex] += value;
								}
							}

							/* data are broken into two lines, read second line here */
							if( read_whole_line( string , (int)sizeof(string) , ioDATA ) == NULL )
							{
								fprintf( ioQQQ, " Badnell data file appears to be corrupted.\n");
								cdEXIT(EXIT_FAILURE);
							}

							/* start of data for second line */
							i1 = 13;
							for(loopindex=10;loopindex<19;loopindex++)
							{
								value=(double)FFmtRead(string,&i1,INPUT_LINE_LENGTH,&lgEOL);
								if( TheirIndexToOurIndex[INDX] < iso.numLevels_max[ipHE_LIKE][nelem] && 
									TheirIndexToOurIndex[INDX] > 0 )
								{
									iso.DielecRecombVsTemp[ipHE_LIKE][nelem][TheirIndexToOurIndex[INDX]][loopindex] += value;
								}
							}
						}
					}
				}
				else
					lgFlag = false;
			}  
			fclose(ioDATA);
			ASSERT( maxINDX > 0 );
			ASSERT( maxINDX < BIGGEST_INDEX_TO_USE );
			ASSERT( max_N_of_data > 0 );

			if( max_N_of_data < iso.n_HighestResolved_max[ipHE_LIKE][nelem] + iso.nCollapsed_max[ipHE_LIKE][nelem] )
			{
				long indexOfMaxN;
				L = -1;
				S = -1; 

				/* This loop extrapolates nLS data to nLS states */
				for( i=TheirIndexToOurIndex[maxINDX]+1; 
					i<iso.numLevels_max[ipHE_LIKE][nelem]-iso.nCollapsed_max[ipHE_LIKE][nelem]; i++ )
				{
					L = L_(i);
					S = S_(i);

					if( L > 4 )
						continue;

					indexOfMaxN = iso.QuantumNumbers2Index[ipHE_LIKE][nelem][max_N_of_data][L][S];
					for(loopindex=0;loopindex<19;loopindex++)
					{
						iso.DielecRecombVsTemp[ipHE_LIKE][nelem][i][loopindex] = 
							iso.DielecRecombVsTemp[ipHE_LIKE][nelem][indexOfMaxN][loopindex] *
							pow( (double)max_N_of_data/(double)StatesElemNEW[nelem][nelem-ipHE_LIKE][i].n, 3.);
					}
				}

				/* Get the N of the highest resolved singlet P (in the model, not the data) */
				indexOfMaxN = 
					iso.QuantumNumbers2Index[ipHE_LIKE][nelem][iso.n_HighestResolved_max[ipHE_LIKE][nelem]][1][1];

				/* This loop extrapolates nLS data to collapsed n levels, just use highest singlet P data */
				for( i=iso.numLevels_max[ipHE_LIKE][nelem]-iso.nCollapsed_max[ipHE_LIKE][nelem];
					i<iso.numLevels_max[ipHE_LIKE][nelem]; i++ )
				{
					for(loopindex=0;loopindex<19;loopindex++)
					{
						iso.DielecRecombVsTemp[ipHE_LIKE][nelem][i][loopindex] = 
							iso.DielecRecombVsTemp[ipHE_LIKE][nelem][indexOfMaxN][loopindex] *
							pow( (double)iso.n_HighestResolved_max[ipHE_LIKE][nelem]/
							(double)StatesElemNEW[nelem][nelem-ipHE_LIKE][i].n, 3.);
					}
				}
			}
		}
	}

	for( i=0; i<NBLOCK; ++i )
	{
		/* set to really large negative number - crash if used before being redefined */
		data_begin_line[i] = INT_MIN;
	}

	chFilename = "badnell_dr.dat";
	ioDATA = open_data( chFilename, "r" );

	count = 0;
	number = 0;

	/*Find out the number line where the data starts 
	 * there are two main blocks of data and each starts with a Z in column 2 */
	while( read_whole_line(chLine, (int)sizeof(chLine), ioDATA) != NULL )
	{
		count++;

		if( chLine[2]=='Z' )
		{
			/* number has to be 0 or 1, and indicates the first or second block of data
			 * count is the line number for the start of that block */
			data_begin_line[number] = count;
			ASSERT( number < NBLOCK );
			number++;
		}
	}

	/*set a flag for a undefined data*/
	if( lgMustMallocRec )
	{
		nsumrec = (long *)MALLOC(LIMELM*sizeof(long) );
		nDRFitPar = (int**)MALLOC( LIMELM*sizeof( int*) );
		lgDRBadnellDefined = (bool **)MALLOC( LIMELM*sizeof(bool*) );
		lgDRBadnellDefinedPart2 = (bool **)MALLOC( LIMELM*sizeof(bool*) );
		lgRRBadnellDefined = (bool **)MALLOC( LIMELM*sizeof(bool*) );

		DRFitParPart1 = (double ***)MALLOC( LIMELM*sizeof(double**) );
		DRFitParPart2 = (double ***)MALLOC( LIMELM*sizeof(double**) );
		RRFitPar = (double ***)MALLOC( LIMELM*sizeof(double**) );
	}

	for( nelem=0; nelem<LIMELM; nelem++ )
	{
		if( lgMustMallocRec )
		{
			nDRFitPar[nelem] = (int*)MALLOC( (nelem+1)*sizeof( int) );
			lgDRBadnellDefined[nelem] = (bool *)MALLOC( (nelem+1)*sizeof(bool) );
			lgDRBadnellDefinedPart2[nelem] = (bool *)MALLOC( (nelem+1)*sizeof(bool) );
			lgRRBadnellDefined[nelem] = (bool *)MALLOC( (nelem+1)*sizeof(bool) );

			DRFitParPart1[nelem] = (double **)MALLOC( (nelem+1)*sizeof(double*) );
			DRFitParPart2[nelem] = (double **)MALLOC( (nelem+1)*sizeof(double*) );
			RRFitPar[nelem] = (double **)MALLOC( (nelem+1)*sizeof(double*) );
		}
		for( long ion=0; ion<nelem+1; ++ion )
		{
			if( lgMustMallocRec )
			{
				DRFitParPart1[nelem][ion] = (double *)MALLOC( MAX_FIT_PAR_DR*sizeof(double) );
				DRFitParPart2[nelem][ion] = (double *)MALLOC( MAX_FIT_PAR_DR*sizeof(double) );
				RRFitPar[nelem][ion] = (double *)MALLOC( MAX_FIT_PAR_RR*sizeof(double) );
			}
			lgDRBadnellDefined[nelem][ion] = false;
			lgDRBadnellDefinedPart2[nelem][ion] = false;   
			lgRRBadnellDefined[nelem][ion] = false;

			/*set fitting coefficients to zero initially*/
			for( k=0; k<MAX_FIT_PAR_DR; k++ )
			{
				DRFitParPart1[nelem][ion][k] = 0;
				DRFitParPart2[nelem][ion][k] = 0;  
			}
			for( k=0; k<MAX_FIT_PAR_RR; k++ )
			{
				RRFitPar[nelem][ion][k] = 0;
			}
		}
	}
	lgMustMallocRec = false;

	count = 0;
	/*Start from beginning to read in again*/  
	fseek(ioDATA, 0, SEEK_SET);
	/* read magic number for DR data */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " DISASTER PROBLEM Badnell_rec_init could not read first line of badnell_dr.dat.\n");
		cdEXIT(EXIT_FAILURE);
	}
	count++;

	/* look for ')' on the line, magic number comes after it */
	if( (chs = strchr_s(chLine, ')'))==NULL )
	{
		/* format is incorrect */
		fprintf( ioQQQ, " DISASTER PROBLEM Badnell_rec_init data file incorrect format.\n");
		cdEXIT(EXIT_FAILURE);
	}

	++chs;
	sscanf(chs, "%4i%2i%2i",&yr, &mo, &dy);
	/* check magic number - the date on the line */
	int dr_yr = 2007, dr_mo = 10, dr_dy = 27;
	if((yr != dr_yr) || (mo != dr_mo) || (dy != dr_dy))
	{
		fprintf(ioQQQ,
			"DISASTER PROBLEM Badnell_rec_init The version of %s I found (%i %i %i) is not the current version (%i %i %i).\n",
			chFilename, yr, mo, dy, dr_yr, dr_mo, dr_dy);
		fprintf(ioQQQ," The first line of the file is the following\n %s\n", chLine );
		cdEXIT(EXIT_FAILURE);
	}

	while( read_whole_line(chLine, (int)sizeof(chLine), ioDATA) != NULL )
	{
		count++;
		length_of_line = (int)strlen(chLine);

		/*reading in coefficients DRFitParPart1 */
		if( count > data_begin_line[0] && count < data_begin_line[1] && length_of_line >3 )
		{
			/*set array par_C to zero */
			for( i=0; i<MAX_FIT_PAR_DR; i++ )
			{
				par_C[i] = 0;
			}
			sscanf(chLine, "%i%i%i%i%lf%lf%lf%lf%lf%lf%lf%lf%lf",
				&NuclearCharge, &NumberElectrons, &M_state, &W_state, &par_C[0], &par_C[1], &par_C[2],
				&par_C[3], &par_C[4], &par_C[5], &par_C[6], &par_C[7], &par_C[8]);
			/* data files have atomic number on physics scale, convert to C scale
			 * for following code */
			long int NuclearChargeM1 = NuclearCharge-1;

			if(M_state == 1 && NuclearChargeM1 < LIMELM )
			{
				/*Set a flag to '1' when the indices are defined */
				ASSERT( NumberElectrons < LIMELM );
				ASSERT( NuclearChargeM1 < LIMELM );
				lgDRBadnellDefined[NuclearChargeM1][NumberElectrons] = true;

				/*counting the number of coefficients */
				nDRFitPar[NuclearChargeM1][NumberElectrons] = 9;
				for( i=8; i>=0; i-- )
				{
					if( par_C[i] == 0 )
						--nDRFitPar[NuclearChargeM1][NumberElectrons];
					else
						break;
				}

				/*assign the values into array */
				for( i=0; i<9; i++ )
					DRFitParPart1[NuclearChargeM1][NumberElectrons][i] = par_C[i];
			}
		}
	}

	/*starting again to read in E values */
	fseek(ioDATA, 0, SEEK_SET);
	count = 0; 
	while( read_whole_line(chLine, (int)sizeof(chLine), ioDATA) != NULL )
	{
		count++;
		length_of_line = (int)strlen(chLine);
		if( count > data_begin_line[1] && length_of_line > 3 )
		{

			/*set array par_E to zero*/
			for( i=0; i<MAX_FIT_PAR_DR; i++ )
			{
				par_E[i] = 0;
			}
			sscanf(chLine, "%i%i%i%i%lf%lf%lf%lf%lf%lf%lf%lf%lf",
				&NuclearCharge, &NumberElectrons, &M_state, &W_state, &par_E[0], &par_E[1], &par_E[2],
				&par_E[3], &par_E[4], &par_E[5], &par_E[6], &par_E[7], &par_E[8]);
			/* data file is on physics scale but we use C scale */
			long int NuclearChargeM1 = NuclearCharge-1;

			if(M_state == 1 && NuclearChargeM1<LIMELM)
			{
				ASSERT( NumberElectrons < LIMELM );
				ASSERT( NuclearChargeM1 < LIMELM );
				lgDRBadnellDefinedPart2[NuclearChargeM1][NumberElectrons] = true;

				/*counting the number of coefficients */
				nDRFitPar[NuclearChargeM1][NumberElectrons] = 9;
				for( i=8; i>=0; i-- )
				{
					if( par_E[i] == 0 )
						--nDRFitPar[NuclearChargeM1][NumberElectrons];
					else
						break;
				}

				/*assign the values into array*/
				for( i=0; i<nDRFitPar[NuclearChargeM1][NumberElectrons]; i++ )
					DRFitParPart2[NuclearChargeM1][NumberElectrons][i] = par_E[i];
			}
		}
	}

	fclose( ioDATA );

	/*output coefficients for defined values for testing */
#	ifdef PRINT_DR
	for( nelem=0; nelem<LIMELM; nelem++ )
	{
		for( int ion=0; ion<nelem+1;++ion )
		{
			if( lgDRBadnellDefined[nelem][ion] )
			{
				fprintf(ofp, "%i %i %e %e %e %e %e %e %e %e %e\n",
					nelem, ion, DRFitParPart1[nelem][ion][0], 
					DRFitParPart1[nelem][ion][1], DRFitParPart1[nelem][ion][2], 
					DRFitParPart1[nelem][ion][3], DRFitParPart1[nelem][ion][4],
					DRFitParPart1[nelem][ion][5], DRFitParPart1[nelem][ion][6], 
					DRFitParPart1[nelem][ion][7], DRFitParPart1[nelem][ion][8]);
			}
		}
	}
	for( nelem=0; nelem<LIMELM; nelem++ )
	{
		for( int ion=0; ion<nelem+1; ion++ )
		{
			if( lgDRBadnellDefinedPart2[nelem][ion] )
			{
				fprintf(ofp, "%i %i %e %e %e %e %e %e %e %e %e\n",
					nelem, ion, DRFitParPart2[nelem][ion][0], 
					DRFitParPart2[nelem][ion][1], DRFitParPart2[nelem][ion][2], 
					DRFitParPart2[nelem][ion][3], DRFitParPart2[nelem][ion][4],
					DRFitParPart2[nelem][ion][5], DRFitParPart2[nelem][ion][6], 
					DRFitParPart2[nelem][ion][7], DRFitParPart2[nelem][ion][8]);
			}
		}
	}
	fclose(ofp);
#	endif

	/*checking for the match of lgDRBadnellDefined and lgDRBadnellDefinedPart2 - 
	 * Both have to be defined*/
	bool lgDRBadnellBothDefined = true;
	for( nelem=0; nelem<LIMELM; nelem++ )
	{
		for( int ion=0; ion<nelem+1; ion++ )
		{
			/* check that first and second half of DR fitting coefficients 
			 * are both defined */
			if( lgDRBadnellDefined[nelem][ion] != lgDRBadnellDefinedPart2[nelem][ion] )
			{
				fprintf( ioQQQ, "DR %i, RR %i: %c %c\n", nelem, ion, 
					 TorF(lgDRBadnellDefined[nelem][ion]), 
					 TorF(lgDRBadnellDefinedPart2[nelem][ion]));
				fprintf( ioQQQ, "PROBLEM ion_recomb_Badnell first and second half of Badnell DR not consistent.\n");
				lgDRBadnellBothDefined = false;
			}
		}
	}

	if( !lgDRBadnellBothDefined )
	{
		/* disaster - DR files are not consistent */
		fprintf(ioQQQ,
			"DISASTER PROBLEM The DR data files are corrupted - part 1 and 2 do not agree.\n");
		fprintf(ioQQQ," Start again with a fresh copy of the data directory\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* now do radiative recombination */
	chFilename = "badnell_rr.dat";
	ioDATA = open_data( chFilename, "r" );

	/* read magic number for RR data */
	{
		if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
		{
			fprintf( ioQQQ, " DISASTER PROBLEM Badnell_rec_init could not read first line of badnell_rr.dat.\n");
			cdEXIT(EXIT_FAILURE);
		}
		/* this is just before date, which we use for magic number */
		if( (chs = strchr_s(chLine, ')'))==NULL )
		{
			/* format is incorrect */
			fprintf( ioQQQ, " DISASTER PROBLEM Badnell_rec_init data file incorrect format.\n");
			cdEXIT(EXIT_FAILURE);
		}
		++chs;
		sscanf(chs, "%4i%2i%2i", &yr, &mo, &dy);
		int rr_yr = 2007, rr_mo = 10, rr_dy = 27;
		if((yr != rr_yr)||(mo != rr_mo)||(dy != rr_dy))
		{
			fprintf(ioQQQ,"DISASTER PROBLEM The version of %s I found (%i %i %i) is not the current version (%i %i %i).\n",
				chFilename, yr, mo, dy, rr_yr, rr_mo, rr_dy);
			fprintf(ioQQQ," The line was as follows:\n %s\n", chLine );
			cdEXIT(EXIT_FAILURE);
		}
	}

	while( read_whole_line(chLine, (int)sizeof(chLine), ioDATA) != NULL )
	{
		/*read in coefficients - first set array par to zero */
		for( i=0; i<MAX_FIT_PAR_RR; i++ )
		{
			temp_par[i] = 0;
		}
		if(chLine[0] != '#')
		{
			sscanf(chLine, "%i%i%i%i%lf%lf%lf%lf%lf%lf",
				&NuclearCharge, &NumberElectrons, &M_state, &W_state, &temp_par[0], &temp_par[1],
				&temp_par[2], &temp_par[3], &temp_par[4], &temp_par[5]);
			long NuclearChargeM1 = NuclearCharge-1;

			if(M_state == 1 && NuclearChargeM1<LIMELM)
			{
				ASSERT( NuclearChargeM1 < LIMELM );
				ASSERT( NumberElectrons <= LIMELM );
				/*Set a flag to '1' when the indices are defined */  
				lgRRBadnellDefined[NuclearChargeM1][NumberElectrons] = true;
				/*assign the values into array */
				for( i=0; i<MAX_FIT_PAR_RR; i++ )
					RRFitPar[NuclearChargeM1][NumberElectrons][i] = temp_par[i];
			}
		}
	}

	/*output coefficients for defined values for testing */
#	ifdef PRINT_RR
	count = 0;
	for( nelem=0; nelem<LIMELM; nelem++ )
	{
		for( ion=0; ion<nelem+1; ion++ )
		{
			if( lgRRBadnellDefined[nelem][ion] )
			{
				fprintf(ofp, "%i %i %e %e %e %e %e %e\n",
					nelem, ion, RRFitPar[nelem][ion][0], 
					RRFitPar[nelem][ion][1], RRFitPar[nelem][ion][2], 
					RRFitPar[nelem][ion][3],
					RRFitPar[nelem][ion][4], RRFitPar[nelem][ion][5]);
				count++;
			}
		}
	}
	fprintf(ofp, "total lines are %i ", count);

	fclose(ofp);
#	endif

	fclose(ioDATA);

	{
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
			{
				int ion;
				fprintf(ioQQQ,"\nDEBUG rr rec\t%i",nelem);
				for( ion=0; ion<=nelem; ++ion )
				{
					fprintf(ioQQQ,"\t%.2e", Badnell_RR_rate_eval(nelem+1 , nelem-ion ) );
				}
				fprintf(ioQQQ,"\n");
				fprintf(ioQQQ,"DEBUG dr rec\t%i",nelem);
				for( ion=0; ion<=nelem; ++ion )
				{
					fprintf(ioQQQ,"\t%.2e", Badnell_DR_rate_eval(nelem+1 , nelem-ion ) );
				}
				fprintf(ioQQQ,"\n");
			}
			cdEXIT(EXIT_FAILURE);
		}
	}
	return;
}

/*ion_recom_calculate calculate radiative and dielectronic recombination rate coefficients */
void ion_recom_calculate( void )
{
	static double TeUsed = -1;
	long int ion , 
		nelem ,
		i;

	DEBUG_ENTRY( "ion_recom_calculate()" );

	/* do not reevaluate if change in temperature is small */
	if( fp_equal(phycon.te,TeUsed) )
		return;

	/* save mean rec coefficients - used to derive rates for those ions with none */
	for( ion=0; ion < LIMELM; ++ion )
		nsumrec[ion] = 0;

	TeUsed = phycon.te;
	/* NB - for following loop important to go over all elements and all
	* ions, not just active ones, since mean DR is used as the guess for
	* DR rates for unknown species. */
	multi_arr<double,2> DR_Badnell_rate_stack( LIMELM, LIMELM );
	for( nelem=ipHYDROGEN; nelem < LIMELM; ++nelem )
	{
		for( ion=0; ion < nelem+1; ++ion )
		{
			long int n_bnd_elec_before_recom ,
				n_bnd_elec_after_recom;

			n_bnd_elec_before_recom = nelem-ion;
			n_bnd_elec_after_recom = nelem-ion+1;
			/* Badnell dielectronic recombination rate coefficients */
			if( (ionbal.DR_Badnell_rate_coef[nelem][ion] = 
				Badnell_DR_rate_eval(
				/* atomic number on C scale */
				nelem, 
				/* number of core electrons before capture of free electron,
				 * for bare ion this is zero */
				n_bnd_elec_before_recom )) >= 0. )
			{
				ionbal.lgDR_Badnell_rate_coef_exist[nelem][ion] = true;
				/* keep track of mean DR rate for this ion */
				DR_Badnell_rate_stack[ion][nsumrec[ion]++] =
					ionbal.DR_Badnell_rate_coef[nelem][ion];
			}
			else
			{
				/* real rate does not exist, will use mean later */
				ionbal.lgDR_Badnell_rate_coef_exist[nelem][ion] = false;
				ionbal.DR_Badnell_rate_coef[nelem][ion] = 0.;
			}

			/* Badnell radiative recombination rate coefficients */
			if( (ionbal.RR_Badnell_rate_coef[nelem][ion] = Badnell_RR_rate_eval(
				/* atomic number on C scale */
				nelem, 
				/* number of core electrons before capture of free electron */
				n_bnd_elec_before_recom )) >= 0. )
			{
				ionbal.lgRR_Badnell_rate_coef_exist[nelem][ion] = true;
			}
			else
			{
				/* real rate does not exist, will use mean later */
				ionbal.lgRR_Badnell_rate_coef_exist[nelem][ion] = false;
				ionbal.RR_Badnell_rate_coef[nelem][ion] = 0.;
			}

			/* save D. Verner's radiative recombination rate coefficient 
			* needed for rec cooling, cm3 s-1 */
			ionbal.RR_Verner_rate_coef[nelem][ion] = 
				t_ADfA::Inst().rad_rec( 
				/* number number of physics scale */
				nelem+1 , 
				/* number of protons on physics scale */
				n_bnd_elec_after_recom , 
				phycon.te );
		}
	}

	/* this branch starts idiosyncratic single ions */
	double Fe_Gu_c[9][6] = {
		{ 2.50507e-11, 5.60226e-11, 1.85001e-10, 3.57495e-9, 1.66321e-7, 0. },/*fit params for Fe+6*/
		{ 9.19610e-11, 2.92460e-10, 1.02120e-9, 1.14852e-8, 3.25418e-7, 0. }, /* fitting params for Fe+7 */
		{ 9.02625e-11, 6.22962e-10, 5.77545e-9, 1.78847e-8, 3.40610e-7, 0. }, /* fitting params for Fe+8 */
		{ 9.04286e-12, 9.68148e-10, 4.83636e-9, 2.48159e-8, 3.96815e-7, 0. }, /* fitting params for Fe+9 */
		{ 6.77873e-10, 1.47252e-9, 5.31202e-9, 2.54793e-8, 3.47407e-7, 0. }, /* fitting params for Fe+10 */
		{ 1.29742e-9, 4.10172e-9, 1.23605e-8, 2.33615e-8, 2.97261e-7, 0. }, /* fitting params for Fe+11 */
		{ 8.78027e-10, 2.31680e-9, 3.49333e-9, 1.16927e-8, 8.18537e-8, 1.54740e-7 },/*fit params for Fe+12*/
		{ 2.23178e-10, 1.87313e-9, 2.86171e-9, 1.38575e-8, 1.17803e-7, 1.06251e-7 },/*fit params for Fe+13*/
		{ 2.17263e-10, 7.35929e-10, 2.81276e-9, 1.32411e-8, 1.15761e-7, 4.80389e-8 }/*fit params for Fe+14*/
	},

		Fe_Gu_E[9][6] = {
			{ 8.30501e-2, 8.52897e-1, 3.40225e0, 2.23053e1, 6.80367e1, 0. }, /* fitting params for Fe+6 */
			{ 1.44392e-1, 9.23999e-1, 5.45498e0, 2.04301e1, 7.06112e1, 0. }, /* fitting params for Fe+7 */ 
			{ 5.79132e-2, 1.27852e0, 3.22439e0, 1.79602e1, 6.96277e1, 0. }, /* fitting params for Fe+8 */
			{ 1.02421e-1, 1.79393e0, 4.83226e0, 1.91117e1, 6.80858e1, 0. }, /* fitting params for Fe+9 */
			{ 1.24630e-1, 6.86045e-1, 3.09611e0, 1.44023e1, 6.42820e1, 0. }, /* fitting params for Fe+10 */
			{ 1.34459e-1, 6.63028e-1, 2.61753e0, 1.30392e1, 6.10222e1, 0. }, /* fitting params for Fe+11 */
			{ 7.79748e-2, 5.35522e-1, 1.88407e0, 8.38459e0, 3.38613e1, 7.89706e1 }, /*fitting params for Fe+12*/
			{ 8.83019e-2, 6.12756e-1, 2.36035e0, 9.61736e0, 3.64467e1, 8.72406e1 }, /*fitting params for Fe+13*/
			{ 1.51322e-1, 5.63155e-1, 2.57013e0, 9.08166e0, 3.69528e1, 1.08067e2 } /* fitting params for Fe+14*/
	};

	double s_c[5][2] = {
		{0.1565e-3, 0.1617e-2},	/* fitting parameters for S+1 */
		{0.3026e-3, 0.2076e-1},	/* fitting parameters for S+2 */
		{0.3177e-1, 0.6309e-3}, /* fitting parameters for S+3 */
		{0.2464e-1, 0.5553e-3}, /* fitting parameters for S+4 */
		{0.1924e-1, 0.5111e-3}	/* fitting parameters for s+5 */
	},

		s_E[5][2] = {
			{0.1157e6, 0.1868e6},	/* fitting parameters for S+1 */
			{0.9662e5, 0.1998e6},	/* fitting parameters for S+2 */
			{0.1928e6, 0.6126e5},	/* fitting parameters for S+3 */
			{0.1806e6, 0.3142e5},	/* fitting parameters for S+4 */
			{0.1519e6, 0.4978e5}	/* fitting parameters for s+5 */
	};

	/* do a series of special cases for Fe DR  */
	nelem = ipIRON;

	/*the sum of the fitting parameter calculations*/
	double fitSum = 0;

	/* Fe+14 - Fe+13 - ion = 0 is A^+ -> A^0 so off by one*/
	ion = 13;
	/* >>chng Fe+14 -> Fe+13, experimental DR from 
	* >>refer	Fe+14	DR	Lukic, D.V. et al 2007, astroph 0704-0905 
	* following is the MCBP entry from Table 3,
	* units of Fe14_c are cm^3 s^-1 K^1.5 */
	double fe14_c[6]={7.07E-04,7.18E-03,2.67E-02,3.15E-02,1.62E-01,5.37E-04};
	/* units of fe14_E are eV THIS IS DIFFERENT FROM OTHER PAPERS BY
	* THESE AUTHORS - THEY CHANGE TEMPERATURE UNITS WITHIN THE SAME
	* PARAGRAPH!!! */
	double fe14_E[6]={4.12E-01,2.06E+00,1.03E+01,2.20E+01,4.22E+01,3.41E+03};
	/* update DR rate - always do this since this reference is more
	* recent that previous paper */
	for( i=0; i<6; i++ )
	{
		fitSum += fe14_c[i] * sexp( fe14_E[i]/phycon.te_eV );
	}

	ionbal.lgDR_Badnell_rate_coef_exist[nelem][ion] = true;
	ionbal.DR_Badnell_rate_coef[nelem][ion] = fitSum / phycon.te32;
	DR_Badnell_rate_stack[ion][nsumrec[ion]++] =
		ionbal.DR_Badnell_rate_coef[nelem][ion];

	/* >>chng 06 jun 21 by Mitchell Martin, added new DR rate coefficient
	* for Fe XIV (Fe+13) to Fe XIII (Fe+12) recombination calculation 
	*>>refer	Fe12	DR	Schmidt E.W. et al. 2006, ApJ, 641, 157L */
	/*fitting parameters used to calculate the DR rate for Fe+13 -> Fe+12*/
	double fe13_c[10]={ 3.55e-4, 2.40e-3, 7.83e-3, 1.10e-2, 3.30e-2, 1.45e-1, 8.50e-2, 2.59e-2, 8.93e-3, 9.80e-3 },
		fe13_E[10]={ 2.19e-2, 1.79e-1, 7.53e-1, 2.21e0, 9.57e0, 3.09e1, 6.37e1, 2.19e2, 1.50e3, 7.86e3 };
	/* do Fe+13 -> Fe+12 from Savin et al. if not already done */
	ion = 12;
	if( !ionbal.lgDR_Badnell_rate_coef_exist[nelem][ion] )
	{
		for( i=0; i<10; i++ )
		{
			fitSum += fe13_c[i] * sexp( fe13_E[i]/phycon.te_eV );
		}

		/* update DR rate is not already done */
		ionbal.lgDR_Badnell_rate_coef_exist[nelem][ion] = true;
		ionbal.DR_Badnell_rate_coef[nelem][ion] = fitSum / phycon.te32;
		DR_Badnell_rate_stack[ion][nsumrec[ion]++] =
			ionbal.DR_Badnell_rate_coef[nelem][ion];
	}

	/* this is the temperature in eV evaluated to the 3/2 power */
	double te_eV32 = pow( phycon.te_eV, 1.5 );

	/* >>chng 06 jul 07 by Mitchell Martin, added DR rate coefficient 
	* calculations for Fe+6->Fe+5 through Fe+14->Fe+13
	* this is still for nelem = ipIRON from the previous calculation 
	* starts with Fe+6 -> Fe+5 and does the next ion with each iteration */
	for( ion=0; ion<9; ion++ )
	{
		/* only do this rate if not already done by a previous approximation */
		if( !ionbal.lgDR_Badnell_rate_coef_exist[nelem][ion+5] )
		{
			fitSum = 0; /* resets the fitting parameter calculation */
			for( i=0; i<6; i++ )
			{
				fitSum += Fe_Gu_c[ion][i] * sexp( Fe_Gu_E[ion][i]/phycon.te_eV );
			}
			ionbal.lgDR_Badnell_rate_coef_exist[nelem][ion+5] = true;
			ionbal.DR_Badnell_rate_coef[nelem][ion+5] = fitSum / te_eV32;
			DR_Badnell_rate_stack[ion+5][nsumrec[ion+5]++] =
				ionbal.DR_Badnell_rate_coef[nelem][ion+5];
		}
	}
	/* this is end of Fe DR rates */

	/* now get mean of the DR rates - may be used for ions with no DR rates */
	for( ion=0; ion < LIMELM; ++ion )
	{
		ionbal.DR_Badnell_rate_coef_mean_ion[ion] = 0.;
		if( nsumrec[ion] > 0 )
		{
			// we have collected all the DR rates above, now we calculate the logarithmic average
			// of the top 3/4 of the entries. this procedure has been checked not to result in
			// discontinuities in the average rate as a function of temperature.
			// NB NB if this code is changed, check again that the resulting average rates
			// are continuous, otherwise discontiuities in the ionization balance and heating-
			// cooling function may result.

			// first sort the stack so that we can easily determine the largest rates
			sort( DR_Badnell_rate_stack.ptr( ion, 0 ), DR_Badnell_rate_stack.ptr( ion, nsumrec[ion] ) );

			double maxval = DR_Badnell_rate_stack[ion][nsumrec[ion]-1];
			if( maxval > 0. )
			{
				// now calculate logarithmic average of top three quarters of entries
				// we discard the lower quarter because many of these are very low or 0
				// NB NB we MUST use a fixed number of values at any temperature to
				// calculate the average since discontinuities would result otherwise!
				double sum = 0.;
				for( int j = nsumrec[ion]/4; j < nsumrec[ion]; ++j )
				{
					double rate_one = log(max(DR_Badnell_rate_stack[ion][j],1e-5*maxval));
					sum += rate_one;
				}
				sum /= (double)(nsumrec[ion] - nsumrec[ion]/4);
				ionbal.DR_Badnell_rate_coef_mean_ion[ion] = exp(sum);
			}
		}
		/*fprintf(ioQQQ,"DEBUG %li %.2e\n", ion , ionbal.DR_Badnell_rate_coef_mean_ion[ion] );*/
	}

	/* >>chng 06 jun 28 by Mitchell Martin, added DR rate coefficient 
	* calculations for SII, SIII, SIV, SV, and SVI*/
	nelem = ipSULPHUR;
	/* starts with S+1 -> S0 and does the next ion with each iteration*/
	for( ion=0; ion<5; ion++ )
	{
		/* only do this rate if not already done by a previous approximation
		* so following used for ion <= 3 */
		if( !ionbal.lgDR_Badnell_rate_coef_exist[nelem][ion] )
		{
			/* >>chng 06 jun 28 by Mitchell Martin, added DR rate 
			* coefficient  for SII, SIII, SIV, SV, and SVI*/
			fitSum = 0;
			for( i=0; i<2; i++ )
			{
				fitSum += s_c[ion][i] * sexp( s_E[ion][i]/phycon.te );
			}
			ionbal.lgDR_Badnell_rate_coef_exist[nelem][ion] = true;
			ionbal.DR_Badnell_rate_coef[nelem][ion] = fitSum / phycon.te32;
			/*>>chng 07 apr 28, use max of this or mean */
			/* change DR data set for S 
			* three cases for S DR - ionbal.nDR_S_guess
			* 0, default larger of guess and Badnell
			* 1, pure Badnell
			* 3, scaled oxygen */
			if( ionbal.nDR_S_guess==0 )
			{
				/* default larger of guess or Badnell */
				ionbal.DR_Badnell_rate_coef[nelem][ion] = 
					MAX2( ionbal.DR_Badnell_rate_coef[nelem][ion] ,
					ionbal.DR_Badnell_rate_coef_mean_ion[ion]);
			}
			else if( ionbal.nDR_S_guess==1 )
			{
				/* pure Badnell - compiler will remove this non op */
				ionbal.DR_Badnell_rate_coef[nelem][ion] = 
					ionbal.DR_Badnell_rate_coef[nelem][ion];
			}
			else if( ionbal.nDR_S_guess==2 )
			{
				/* scaled oxygen */
				ionbal.DR_Badnell_rate_coef[nelem][ion] = 
					ionbal.DR_Badnell_rate_coef[ipOXYGEN][ion]*
					ionbal.DR_S_scale[ion];
			}
			else
				TotalInsanity();
		}
	}

	/* this set true with PRINT on any of the Badnell set recombination commands */
	if( ionbal.lgRecom_Badnell_print )
	{
		fprintf(ioQQQ,"DEBUG Badnell recombination RR, then DR, T=%.3e\n", phycon.te );
		for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
		{
			fprintf(ioQQQ,"nelem=%li %s, RR then DR\n",
				nelem , elementnames.chElementNameShort[nelem] );
			for( ion=0; ion<nelem+1; ++ion )
			{
				if( ionbal.lgRR_Badnell_rate_coef_exist[nelem][ion] )
				{
					fprintf(ioQQQ," %.2e", ionbal.RR_Badnell_rate_coef[nelem][ion] );
				}
				else
				{
					fprintf(ioQQQ," %.2e", -1. );
				}
			}
			fprintf(ioQQQ,"\n" );
			for( ion=0; ion<nelem+1; ++ion )
			{
				if( ionbal.lgDR_Badnell_rate_coef_exist[nelem][ion] )
				{
					fprintf(ioQQQ," %.2e", ionbal.DR_Badnell_rate_coef[nelem][ion] );
				}
				else
				{
					fprintf(ioQQQ," %.2e", -1. );
				}
			}
			fprintf(ioQQQ,"\n\n" );
		}
		/* now print mean recombination and standard deviation */
		fprintf(ioQQQ,"mean DR recombination ion mean stddev stddev/mean \n" );
		for( ion=0; ion<LIMELM; ++ion )
		{
			double stddev;
			stddev = 0.;
			for( nelem=ion; nelem<LIMELM; ++nelem )
			{
				stddev += POW2( 
					ionbal.DR_Badnell_rate_coef[nelem][ion]- 
					ionbal.DR_Badnell_rate_coef_mean_ion[ion] );
			}
			stddev = sqrt( stddev / MAX2( 1 , nsumrec[ion]-1 ) );
			fprintf(ioQQQ," %2li %.2e %.2e %.2e\n", 
				ion , 
				ionbal.DR_Badnell_rate_coef_mean_ion[ion] , 
				stddev ,
				stddev / SDIV(ionbal.DR_Badnell_rate_coef_mean_ion[ion]) );
		}
		cdEXIT( EXIT_SUCCESS );
	}
	return;
}
