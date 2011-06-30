/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
* others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "lines_service.h"
#include "physconst.h"
#include "taulines.h"
#include "trace.h"
#include "string.h"
#include "thirdparty.h"
#include <fstream>

static const bool DEBUGSTATE = false;

void atmdat_CHIANTI_readin( long intNS, char *chPrefix )
{
	DEBUG_ENTRY( "atmdat_CHIANTI_readin()" );

	int intsplinepts,intTranType,intxs;
	long int nMolLevs;
	FILE *ioElecCollData=NULL, *ioProtCollData=NULL;
	realnum  fstatwt,fenergyK,fenergyWN,fWLAng,fenergy,feinsteina;
	double fScalingParam,fEnergyDiff,fGF,*xs,*spl,*spl2;

	char chLine[FILENAME_PATH_LENGTH_2] ,
		chEnFilename[FILENAME_PATH_LENGTH_2],
		chTraFilename[FILENAME_PATH_LENGTH_2],
		chEleColFilename[FILENAME_PATH_LENGTH_2],
		chProColFilename[FILENAME_PATH_LENGTH_2];

	realnum *dBaseStatesEnergy;
	long int *intNewIndex=NULL,*intOldIndex=NULL;
	int interror;
	bool lgProtonData=false,lgEneLevOrder;

	// this is the largest number of levels allowed by the chianti format, I3
	static const int MAX_NUM_LEVELS = 999;

	Species[intNS].lgMolecular = false;

	strcpy( chEnFilename , chPrefix );
	strcpy( chTraFilename , chPrefix );	
	strcpy( chEleColFilename , chPrefix );		
	strcpy( chProColFilename , chPrefix );			

	/*For the CHIANTI DATABASE*/
	/*Open the energy levels file*/
	strcat( chEnFilename , ".elvlc");
	uncaps( chEnFilename );
	if(DEBUGSTATE)
		printf("The data file is %s \n",chEnFilename);

	/*Open the files*/
	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_CHIANTI_readin opening %s:",chEnFilename);

	fstream elvlcstream;
	open_data(elvlcstream, chEnFilename,mode_r);

	/*Open the transition probabilities file*/
	strcat( chTraFilename , ".wgfa");
	uncaps( chTraFilename );
	if(DEBUGSTATE)
		printf("The data file is %s \n",chTraFilename);

	/*Open the files*/
	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_CHIANTI_readin opening %s:",chTraFilename);

	fstream wgfastream;
	open_data(wgfastream, chTraFilename,mode_r);

	/*Open the electron collision data*/
	strcat( chEleColFilename , ".splups");
	uncaps( chEleColFilename );
	if(DEBUGSTATE)
		printf("The data file is %s \n",chEleColFilename);
	/*Open the files*/
	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_CHIANTI_readin opening %s:",chEleColFilename);

	ioElecCollData = open_data( chEleColFilename, "r" );

	/*Open the proton collision data*/
	strcat( chProColFilename , ".psplups");
	uncaps( chProColFilename );
	if(DEBUGSTATE)
		printf("The data file is %s \n",chProColFilename);
	/*Open the files*/
	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_CHIANTI_readin opening %s:",chProColFilename);

	/*We will set a flag here to indicate if the proton collision strengths are available */
	if( ( ioProtCollData = fopen( chProColFilename , "r" ) ) != NULL )
	{
		lgProtonData = true;
		//fclose( ioProtCollData );
		//ioProtCollData = NULL;
	}
	else
	{
		lgProtonData = false;
	}

	//nMolLevs set to -1 since while loop counts 1 extra
	//eof_col is used get the first 4 charcters per line to find end of file
	const int eof_col = 5;
	nMolLevs = -1;
	lgEneLevOrder = true;
	if (elvlcstream.is_open())
	{
		int nj = 0;
		char otemp[eof_col];
		//Look for -1 to determine number of levels in elvlc file
		while(nj != -1)
		{
			elvlcstream.get(otemp,eof_col);
			elvlcstream.ignore(INT_MAX,'\n');
			nj = atoi(otemp);
			nMolLevs++;
		}
		elvlcstream.seekg(0,ios::beg);
	}

	long HighestIndexInFile = nMolLevs;
	nMolLevs = MIN2( nMolLevs, MAX_NUM_LEVELS );

	if( nMolLevs <= 0 )
	{
		fprintf( ioQQQ, "The number of energy levels is non-positive in datafile %s.\n", chEnFilename );
		fprintf( ioQQQ, "The file must be corrupted.\n" );
		cdEXIT( EXIT_FAILURE );
	}

	Species[intNS].numLevels_max = nMolLevs;
	Species[intNS].numLevels_local = Species[intNS].numLevels_max;

	/*Malloc the energy array to check if the energy levels are in order*/
	dBaseStatesEnergy = (realnum*)MALLOC((unsigned long)(nMolLevs)*sizeof(realnum));
	/*malloc the States array*/
	dBaseStates[intNS] = (quantumState *)MALLOC((size_t)(nMolLevs)*sizeof(quantumState));
	/*malloc the Transition array*/
	dBaseTrans[intNS] = (transition **)MALLOC(sizeof(transition*)*(unsigned long)nMolLevs);
	dBaseTrans[intNS][0] = NULL;
	for( long ipHi = 1; ipHi < nMolLevs; ipHi++)
	{
		dBaseTrans[intNS][ipHi] = (transition *)MALLOC(sizeof(transition)*(unsigned long)ipHi);
		for( long ipLo = 0; ipLo < ipHi; ipLo++)
		{
			dBaseTrans[intNS][ipHi][ipLo].Junk();
			dBaseTrans[intNS][ipHi][ipLo].Lo = &dBaseStates[intNS][ipLo];
			dBaseTrans[intNS][ipHi][ipLo].Hi = &dBaseStates[intNS][ipHi];
		}
	}

	for( long ipLev = 0; ipLev<NUM_COLLIDERS; ipLev++ )
	{
		CollRatesArray[intNS][ipLev] = NULL;
	}
	/*Allocating space just for the electron*/
	CollRatesArray[intNS][0] = (double**)MALLOC((nMolLevs)*sizeof(double*));
	for( long ipHi = 0; ipHi<nMolLevs; ipHi++ )
	{
		CollRatesArray[intNS][0][ipHi] = (double*)MALLOC((unsigned long)(nMolLevs)*sizeof(double));
		for( long ipLo = 0; ipLo<nMolLevs; ipLo++)
		{
			CollRatesArray[intNS][0][ipHi][ipLo] = 0.0;
		}
	}
	/*Allocating space for the proton*/
	if(lgProtonData)
	{
		CollRatesArray[intNS][1] = (double**)MALLOC((nMolLevs)*sizeof(double*));
		for( long ipHi = 0; ipHi<nMolLevs; ipHi++ )
		{
			CollRatesArray[intNS][1][ipHi] = (double*)MALLOC((unsigned long)(nMolLevs)*sizeof(double));
			for( long ipLo = 0; ipLo<nMolLevs; ipLo++)
			{
				CollRatesArray[intNS][1][ipHi][ipLo] = 0.0;
			}
		}
	}

	/*Reading in the energy and checking that they are in order*/
	// C++ io works in terms of cursor movement rather than position on line
	//length (+1) of the nrg in the elvlc file
	const int lvl_nrg_col=16;
	//# of columns skipped from the left to get to nrg start
	const int lvl_skipto_nrg = 40;
	//# of columns to skip over the rydberg energy, we don't use it
	const int lvl_skip_ryd = 15;

	//Read in nrg levels to see if they are in order
	for( long ipLev=0; ipLev<nMolLevs; ipLev++ )
	{
		if(elvlcstream.is_open())
		{
			char thtemp[lvl_nrg_col],obtemp[lvl_nrg_col];
			elvlcstream.seekg(lvl_skipto_nrg,ios::cur);
			elvlcstream.get(thtemp,lvl_nrg_col);
			elvlcstream.seekg(lvl_skip_ryd,ios::cur);
			fenergy = (realnum) atof(thtemp);
			//This if was used originally used to read exp nrgs when the
			//theo nrgs were not present. Now it only looks for theo nrg
			//if the exp nrg is zero.
			if(fenergy == 0. && elvlcstream.peek() !='\n')
			{
				elvlcstream.get(obtemp,lvl_nrg_col);
				if(atof(obtemp) != 0.)
				{
					fenergy = (realnum) atof(obtemp);
				}
				//else
					//fprintf( ioQQQ," WARNING: Switched to theoretical energies for species %s, level %3li\n", Species[intNS].chLabel, ipLev );
			}

			elvlcstream.ignore(INT_MAX,'\n');
			dBaseStatesEnergy[ipLev] = fenergy;
		}
		else
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEnFilename);
			cdEXIT(EXIT_FAILURE);
		}

		/*To check if energy levels are in order*/
		/*If the energy levels are not in order a flag is set*/
		if(ipLev>0)
		{
			if (dBaseStatesEnergy[ipLev] < dBaseStatesEnergy[ipLev-1])
			{
				lgEneLevOrder = false;	
			}
		}
	}

	// malloc vector for new energy-sorted indices
	intNewIndex = (long int*)MALLOC((unsigned long)(nMolLevs)* sizeof(long int));

	if(lgEneLevOrder == false)
	{
		/*If the energy levels are not in order and the flag is set(lgEneLevOrder=FALSE)
		then we sort the energy levels to find the relation matrix between the old and new indices*/
		intOldIndex = (long int*)MALLOC((unsigned long)(nMolLevs)* sizeof(long int));
		/*We now do a modified quick sort*/
		spsort(dBaseStatesEnergy,nMolLevs,intOldIndex,2,&interror);
		/*intNewIndex has the key to map*/
		for( long i=0; i<nMolLevs; i++ )
		{
			ASSERT( intOldIndex[i] < nMolLevs );
			intNewIndex[intOldIndex[i]] = i;	
		}
		if(DEBUGSTATE)
		{
			for( long i=0; i<nMolLevs; i++ )
			{
				printf("The %ld value of the relation matrix is %ld \n",i,intNewIndex[i]);
			}
			for( long i=0; i<nMolLevs; i++ )
			{
				printf("The %ld value of the sorted energy array is %f \n",i,dBaseStatesEnergy[i]);
			}
		}
		free( intOldIndex );
	}
	else
	{
		/* if energies already in correct order, new index is same as original. */
		for( long i=0; i<nMolLevs; i++ )
		{
			intNewIndex[i] = i;	
		}
	}

	/* insist that there the intNewIndex values are unique */
	for( long i=0; i<nMolLevs; i++ )
	{
		for( long j=i+1; j<nMolLevs; j++ )
		{
			ASSERT( intNewIndex[i] != intNewIndex[j] );	
		}
	}


	//After we read in the energies we rewind the energy levels file
	elvlcstream.seekg(0,ios::beg);

	//lvl_skipto_statwt is the # of columns to skip to statwt from left
	const int lvl_skipto_statwt = 37;
	//lvl_statwt_col is the length (+1) of statwt
	const int lvl_statwt_col = 4;
	//Read in stat weight and energy
	for( long ipLevOld=0; ipLevOld<nMolLevs; ipLevOld++ )
	{
		long ipLevNew = intNewIndex[ipLevOld];
		char gtemp[lvl_statwt_col],thtemp[lvl_nrg_col],obtemp[lvl_nrg_col];
	
		/*information needed for label*/
		strcpy(dBaseStates[intNS][ipLevNew].chLabel, "    ");
		strncpy(dBaseStates[intNS][ipLevNew].chLabel,Species[intNS].chLabel, 4);
		dBaseStates[intNS][ipLevNew].chLabel[4] = '\0';
		
		//associate with species
		dBaseStates[intNS][ipLevNew].sp = &Species[intNS];

		//Move cursor to and read statwt
		elvlcstream.seekg(lvl_skipto_statwt,ios::cur);
		elvlcstream.get(gtemp,lvl_statwt_col);
		fstatwt = (realnum)atof(gtemp);

		if(fstatwt <= 0.)
		{
			fprintf( ioQQQ, " WARNING: A positive non zero value is expected for the statistical weight but was not found in %s\n", chEnFilename);
			fstatwt = 1.f;
			//cdEXIT(EXIT_FAILURE);
		}
		dBaseStates[intNS][ipLevNew].g = fstatwt;

		//Read nrg
		elvlcstream.get(thtemp,lvl_nrg_col);
		elvlcstream.seekg(lvl_skip_ryd,ios::cur);
		fenergy = (realnum) atof(thtemp);
		//This if was used originally used to read exp nrgs when the
		//theo nrgs were not present. Now it only looks for theo nrg
		//if the exp nrg is zero.
		if(fenergy == 0. && elvlcstream.peek() !='\n')
		{
			elvlcstream.get(obtemp,lvl_nrg_col);
			if(atof(obtemp) != 0.)
			{
			fenergy = (realnum) atof(obtemp);
			}
		}		
		elvlcstream.ignore(INT_MAX,'\n');
		dBaseStates[intNS][ipLevNew].energy.set(fenergy,"cm^-1");
	}
	//Close the level stream
	elvlcstream.close();

	/* fill in all transition energies, can later overwrite for specific radiative transitions */
	for( long ipHi=1; ipHi<nMolLevs; ipHi++ )
	{
		for( long ipLo=0; ipLo<ipHi; ipLo++ )
		{
			fenergyWN = (realnum)MAX2( 1E-10, dBaseStates[intNS][ipHi].energy.WN() - dBaseStates[intNS][ipLo].energy.WN() );
			fenergyK = (realnum)(fenergyWN*T1CM);

			dBaseTrans[intNS][ipHi][ipLo].EnergyK = fenergyK;
			dBaseTrans[intNS][ipHi][ipLo].EnergyWN = fenergyWN;
			dBaseTrans[intNS][ipHi][ipLo].EnergyErg = (realnum)ERG1CM *fenergyWN;
			dBaseTrans[intNS][ipHi][ipLo].WLAng = (realnum)(1e+8/fenergyWN/RefIndex(fenergyWN));
		}
	}

	/************************************************************************/
	/*** Read in the level numbers, Einstein A and transition wavelength  ***/
	/************************************************************************/

	//Count the number of rows first
	long wgfarows = -1;
	if (wgfastream.is_open())
	{
		int nj = 0;
		char otemp[eof_col];
		while(nj != -1)
		{
			wgfastream.get(otemp,eof_col);
			wgfastream.ignore(INT_MAX,'\n');
			nj = atoi(otemp);
			wgfarows++;
		}
		wgfastream.seekg(0,ios::beg);
	}
	else 
		fprintf( ioQQQ, " The data file %s is corrupted .\n",chTraFilename);

	//line_index_col is the length(+1) of the level indexes in the WGFA file
	const int line_index_col = 6;
	//line_nrg_to_eina is the # of columns to skip from wavelength to eina in WGFA file
	const int line_nrg_to_eina =15;
	//line_eina_col is the length(+1) of einsteinA in WGFA
	const int line_eina_col = 16;
	char lvltemp[line_index_col];
	//Start reading WGFA file
	if (wgfastream.is_open())
	{
		for (long i = 0;i<wgfarows;i++)
		{
			wgfastream.get(lvltemp,line_index_col);
			long ipLo = atoi(lvltemp)-1;
			wgfastream.get(lvltemp,line_index_col);
			long ipHi = atoi(lvltemp)-1;

			if( ipLo >= nMolLevs || ipHi >= nMolLevs )
			{
				// skip these lines
				wgfastream.ignore(INT_MAX,'\n');
				continue;
			}
	
			if( ipHi == ipLo )
			{
				fprintf( ioQQQ," WARNING: Upper level = lower for a radiative transition in %s. Ignoring.\n", chTraFilename );
				wgfastream.ignore(INT_MAX,'\n');
				continue;
			}
	
			/* translate to properly sorted indices */
			ipLo = intNewIndex[ipLo];
			ipHi = intNewIndex[ipHi];
	
			ASSERT( ipHi != ipLo );
			ASSERT(ipHi >= 0);
			ASSERT(ipLo >= 0);

			// sometimes the CHIANTI datafiles list the highest index first as in the middle of these five lines in ne_10.wgfa:
			//    ...
			//    8   10       187.5299      0.000e+00      4.127e+05                 3d 2D1.5 -                  4s 2S0.5           E2
			//    9   10       187.6573      0.000e+00      6.197e+05                 3d 2D2.5 -                  4s 2S0.5           E2
			//   11   10   4842624.0000      1.499e-05      9.423e-06                 4p 2P0.5 -                  4s 2S0.5           E1
			//    1   11         9.7085      1.892e-02      6.695e+11                 1s 2S0.5 -                  4p 2P0.5           E1
			//    2   11        48.5157      6.787e-02      9.618e+10                 2s 2S0.5 -                  4p 2P0.5           E1
			//    ...
			// so, just set ipHi (ipLo) equal to the max (min) of the two indices.
			// NB NB NB it looks like this may depend upon whether one uses observed or theoretical energies.

			if( ipHi < ipLo )
			{
				long ipLoTemp = ipLo;
				long ipHiTemp = ipHi;
				ipHi = max( ipHiTemp, ipLoTemp );
				ipLo = min( ipHiTemp, ipLoTemp );
				//fprintf( ioQQQ," WARNING: Swapped level indices for species %6s, indices %3li %3li.\n", Species[intNS].chLabel, ipLoTemp, ipHiTemp );
			}

			//Read in wavelengths
			char trantemp[lvl_nrg_col];
			wgfastream.get(trantemp,lvl_nrg_col);
			fWLAng = (realnum)atof(trantemp);

			/* \todo 2 CHIANTI labels the H 1 2-photon transition as z wavelength of zero.
			 * Should we just ignore all of the wavelengths in this file and use the
			 * difference of level energies instead. */
			if( fWLAng <= 0. )
			{
				//if( fWLAng < 0.)
					//fprintf( ioQQQ," WARNING: Negative wavelength for species %6s, indices %3li %3li \n", Species[intNS].chLabel, ipLo, ipHi);
				fWLAng = (realnum)(1e8/
					(dBaseStates[intNS][ipHi].energy.WN() - dBaseStates[intNS][ipLo].energy.WN()));
			}
			//Skip from end of Wavelength to Einstein A and read in
			wgfastream.seekg(line_nrg_to_eina,ios::cur);
			wgfastream.get(trantemp,line_eina_col);
			feinsteina = (realnum)atof(trantemp);
			if( feinsteina < 1e-20 )
			{
				static bool notPrintedYet = true;
				if( notPrintedYet )
				{
					fprintf( ioQQQ," WARNING: Radiative rate(s) equal to zero in %s.\n", chTraFilename );
					notPrintedYet = false;
				}
				wgfastream.ignore(INT_MAX,'\n');
				continue;
			}

			fixit(); // may need to do something with these rates
			//Read in the rest of the line and look for auto
			wgfastream.getline(chLine,INT_MAX);
			if( nMatch("auto", chLine) )
			{
				if( dBaseTrans[intNS][ipHi][ipLo].Emis != NULL )
				{
					dBaseTrans[intNS][ipHi][ipLo].Emis->AutoIonizFrac =
						feinsteina/(dBaseTrans[intNS][ipHi][ipLo].Emis->Aul + feinsteina);
					ASSERT( dBaseTrans[intNS][ipHi][ipLo].Emis->AutoIonizFrac >= 0. );
					ASSERT( dBaseTrans[intNS][ipHi][ipLo].Emis->AutoIonizFrac <= 1. );
				}
				continue;
			}

			dBaseTrans[intNS][ipHi][ipLo].Emis
				= AddLine2Stack(feinsteina, &dBaseTrans[intNS][ipHi][ipLo]);
			dBaseTrans[intNS][ipHi][ipLo].WLAng = fWLAng;
		 
			fenergyWN = (realnum)(1e+8/fWLAng);
			fenergyK = (realnum)(fenergyWN*T1CM);
			dBaseTrans[intNS][ipHi][ipLo].EnergyK = fenergyK;
			dBaseTrans[intNS][ipHi][ipLo].EnergyWN = fenergyWN;
			dBaseTrans[intNS][ipHi][ipLo].EnergyErg = (realnum)ERG1CM *fenergyWN;
			dBaseTrans[intNS][ipHi][ipLo].WLAng = (realnum)(1e+8/fenergyWN/RefIndex(fenergyWN));
			dBaseTrans[intNS][ipHi][ipLo].Emis->gf = (realnum)GetGF(dBaseTrans[intNS][ipHi][ipLo].Emis->Aul,
				dBaseTrans[intNS][ipHi][ipLo].EnergyWN, dBaseTrans[intNS][ipHi][ipLo].Hi->g);

			if(DEBUGSTATE)
			{
				fprintf( ioQQQ, "The lower level is %ld \n",ipLo);
				fprintf( ioQQQ, "The upper level is %ld \n",ipHi);
				fprintf( ioQQQ, "The Einstein A is %f \n",dBaseTrans[intNS][ipHi][ipLo].Emis->Aul);
				fprintf( ioQQQ, "The wavelength of the transition is %f \n",dBaseTrans[intNS][ipHi][ipLo].WLAng );
			}
		}
	}
	else fprintf( ioQQQ, " The data file %s is corrupted .\n",chTraFilename);
	wgfastream.close();
	
	/* Malloc space for splines */
	AtmolCollSplines[intNS] = (CollSplinesArray***)MALLOC(nMolLevs *sizeof(CollSplinesArray**));
	for( long ipHi=0; ipHi<nMolLevs; ipHi++ )
	{
		AtmolCollSplines[intNS][ipHi] = (CollSplinesArray **)MALLOC((unsigned long)(nMolLevs)*sizeof(CollSplinesArray *));
		for( long ipLo=0; ipLo<nMolLevs; ipLo++ )
		{
			AtmolCollSplines[intNS][ipHi][ipLo] = 
				(CollSplinesArray *)MALLOC((unsigned long)(NUM_COLLIDERS)*sizeof(CollSplinesArray ));

			for( long k=0; k<NUM_COLLIDERS; k++ )
			{
				/* initialize all spline variables */
				AtmolCollSplines[intNS][ipHi][ipLo][k].collspline = NULL;
				AtmolCollSplines[intNS][ipHi][ipLo][k].SplineSecDer = NULL;
				AtmolCollSplines[intNS][ipHi][ipLo][k].nSplinePts = -1; 
				AtmolCollSplines[intNS][ipHi][ipLo][k].intTranType = -1;
				AtmolCollSplines[intNS][ipHi][ipLo][k].gf = BIGDOUBLE;
				AtmolCollSplines[intNS][ipHi][ipLo][k].EnergyDiff = BIGDOUBLE;
				AtmolCollSplines[intNS][ipHi][ipLo][k].ScalingParam = BIGDOUBLE;
			}
		}
	}

	/************************************/
	/*** Read in the collisional data ***/
	/************************************/

	// ipCollider 0 is electrons, 1 is protons
	for( long ipCollider=0; ipCollider<=1; ipCollider++ )
	{
		char chFilename[FILENAME_PATH_LENGTH_2];

		if( ipCollider==0 )
		{
			strcpy( chFilename, chEleColFilename );
		}
		else if( ipCollider==1 )
		{
			if( !lgProtonData )
				break;
			fprintf( ioQQQ," WARNING: Chianti proton collision data have different format than electron data files!\n" );
			strcpy( chFilename, chProColFilename );
		}
		else
			TotalInsanity();

		/*Dummy string used for convenience*/
		strcpy(chLine,"A");

		fstream splupsstream;
		open_data(splupsstream, chFilename,mode_r);

		//cs_eof_col is the length(+1) of the first column used for finding the end of file
		const int cs_eof_col = 4;
		//cs_index_col is the length(+1) of the indexes in the CS file
		const int cs_index_col = 4;
		//cs_trantype_col is the length(+1) of the transition type in the CS file
		const int cs_trantype_col = 4;
		//cs_values_col is the length(+1) of the other values in the CS file
		//including: GF, nrg diff, scaling parameter, and spline points
		const int cs_values_col = 11;
		//Determine the number of rows in the CS file
		if (splupsstream.is_open())
		{
			int nj = 0;
			//splupslines is -1 since the loop runs 1 extra time
			long splupslines = -1;
			char otemp[cs_eof_col];
			while(nj != -1)
			{
				splupsstream.get(otemp,cs_eof_col);
				splupsstream.ignore(INT_MAX,'\n');
				nj = atoi(otemp);
				splupslines++;
			}
			splupsstream.seekg(0,ios::beg);
	
			for (int m = 0;m<splupslines;m++)
			{
				if( ipCollider == 0 )
				{
					splupsstream.seekg(6,ios::cur);
				}

				/* level indices */
				splupsstream.get(otemp,cs_index_col);
				long ipLo = atoi(otemp)-1;
				splupsstream.get(otemp,cs_index_col);
				long ipHi = atoi(otemp)-1;

				if( ipLo >= nMolLevs || ipHi >= nMolLevs )
				{
					// skip these transitions
					splupsstream.ignore(INT_MAX,'\n');
					continue;
				}

				if( ipLo >= HighestIndexInFile || ipHi >= HighestIndexInFile )
				{
					fprintf( ioQQQ," WARNING: Problem with indices in CHIANTI file %s.  Often due to indices > 100 in fixed format. Disabling model.\n", chFilename );
					fixit(); // think about easy way to disable and move on.  Zero numLevels_max?  Point something to NULL?
					break;
				}

				ASSERT( ipLo < nMolLevs );
				ASSERT( ipHi < nMolLevs );
				/*Transition Type*/
				splupsstream.get(otemp,cs_trantype_col);
				intTranType = atoi(otemp);
				/*gf*/
				char qtemp[cs_values_col];
				splupsstream.get(qtemp,cs_values_col);
				fGF = atof(qtemp);
				/*Energy Difference*/
				splupsstream.get(qtemp,cs_values_col);
				fEnergyDiff = atof(qtemp);
				/*Scaling Parameter*/
				splupsstream.get(qtemp,cs_values_col);
				fScalingParam = atof(qtemp);

				/* translate to properly sorted indices */
				ipLo = intNewIndex[ipLo];
				ipHi = intNewIndex[ipHi];

				ASSERT( ipLo < nMolLevs );
				ASSERT( ipHi < nMolLevs );
				ASSERT( AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].collspline == NULL );
				ASSERT( AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].SplineSecDer == NULL );

				const int CHIANTI_SPLINE_MAX=9, CHIANTI_SPLINE_MIN=5;
				STATIC_ASSERT(CHIANTI_SPLINE_MAX > CHIANTI_SPLINE_MIN);

				/*We malloc the space here*/
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].collspline =
					(double *)MALLOC((unsigned long)(CHIANTI_SPLINE_MAX)*sizeof(double));
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].SplineSecDer =
					(double *)MALLOC((unsigned long)(CHIANTI_SPLINE_MAX)*sizeof(double));

				/* always read at least CHIANTI_SPLINE_MIN */
				for( intsplinepts=0; intsplinepts<=CHIANTI_SPLINE_MAX; intsplinepts++ )
				{
					//Look at the next character to see if it is the end of line.
					char p = splupsstream.peek();
					if( p == '\n' )
					{
						/* there should be 5 or 9 spline points.  If we got EOL,
						 * insist that we were trying to read the 6th or 10th. */
						if( (intsplinepts != CHIANTI_SPLINE_MAX) && (intsplinepts != CHIANTI_SPLINE_MIN) )
						{
							static bool notPrintedYet = true;
							if( notPrintedYet )
							{
								fprintf( ioQQQ, " WARNING: Wrong number of spline points in %s.\n", chFilename);
								notPrintedYet = false;
							}
							for( long j=intsplinepts-1; j<CHIANTI_SPLINE_MAX; j++ )
								AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].collspline[j] = 0.;
						}
						else
							break;
					}
					else
					{
						if( intsplinepts >= CHIANTI_SPLINE_MAX )
						{
							fprintf( ioQQQ, " WARNING: More spline points than expected in %s, indices %3li %3li.  Ignoring extras.\n", chFilename, ipHi, ipLo );
							break;
						}
						ASSERT( intsplinepts < CHIANTI_SPLINE_MAX );
						double temp;
						//Store a single spline point then look for more
						splupsstream.get(qtemp,cs_values_col);
						temp = atof(qtemp);
						if(temp < 0)
						{
							temp = 0.;
						}
						AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].collspline[intsplinepts] = temp;
					}

				}

				ASSERT( (intsplinepts == CHIANTI_SPLINE_MAX) || (intsplinepts == CHIANTI_SPLINE_MIN) );

				/*The zeroth element contains the number of spline points*/
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].nSplinePts = intsplinepts;
				/*Transition type*/
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].intTranType = intTranType;
				/*gf value*/
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].gf = fGF;
				/*Energy difference between two levels in Rydbergs*/
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].EnergyDiff = fEnergyDiff;
				/*Scaling parameter C*/
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].ScalingParam = fScalingParam;

				/*Once the spline points have been filled,fill the second derivatives structure*/
				/*Creating spline points array*/
				xs = (double *)MALLOC((unsigned long)intsplinepts*sizeof(double));
				spl = (double *)MALLOC((unsigned long)intsplinepts*sizeof(double));
				spl2 = (double *)MALLOC((unsigned long)intsplinepts*sizeof(double));

				// should be able to just loop to intsplinepts -- ASSERT above protects
				if(intsplinepts == CHIANTI_SPLINE_MIN)
				{
					for(intxs=0;intxs<CHIANTI_SPLINE_MIN;intxs++)
					{
						xs[intxs] = 0.25*intxs;
						spl[intxs] = AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].collspline[intxs];
					}
				}
				else if(intsplinepts == CHIANTI_SPLINE_MAX)
				{
					for(intxs=0;intxs<CHIANTI_SPLINE_MAX;intxs++)
					{
						xs[intxs] = 0.125*intxs;
						spl[intxs] = AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].collspline[intxs];
					}
				}
				else
					TotalInsanity();

				spline(xs, spl,intsplinepts,2e31,2e31,spl2);

				/*Filling the second derivative structure*/
				for( long i=0; i<intsplinepts; i++)
				{
					AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].SplineSecDer[i] = spl2[i];
				}

				free( xs );
				free( spl );
				free( spl2 );
				splupsstream.ignore(INT_MAX,'\n');
			}
			splupsstream.close();
		}
		
	}
	
	// free some memory
	free( dBaseStatesEnergy );
	free( intNewIndex );

	// close open file handles
	fclose( ioElecCollData );
	if( lgProtonData )
		fclose( ioProtCollData );

	return;
}
