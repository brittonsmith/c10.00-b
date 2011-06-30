/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
* others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "lines_service.h"
#include "physconst.h"
#include "taulines.h"
#include "trace.h"
#include "string.h"
#include "thirdparty.h"
#include "mole.h"

#define DEBUGSTATE false

/*Separate Routine to read in the molecules*/
void atmdat_LAMDA_readin( long intNS, char *chEFilename )
{
	DEBUG_ENTRY( "atmdat_LAMDA_readin()" );

	int nMolLevs = -10000,intCollIndex = -10000,intLColliderIndex = -10000;
	/* type is set to 0 for non CHIANTI and 1 for CHIANTI*/
	long int intlnct,intrtct,intgrtct,intCollPart,
		intDCollPart,intCollTran,nCollTrans,intCollTemp,intcollindex;
	/*intrtct refers to radiative transitions count*/
	FILE *ioLevData;
	realnum  fstatwt,fenergyK,fenergyWN,fenergy,feinsteina;

	char chLine[FILENAME_PATH_LENGTH_2] ,
		//chEFilename[FILENAME_PATH_LENGTH_2],
		*chColltemp,*chCollName;

	ASSERT( intNS >= 0 );

	const int MAX_NUM_LEVELS = 70;

	Species[intNS].lgMolecular = true;

	/*Load the species name in the Species array structure*/
	if(DEBUGSTATE)
		printf("The name of the %li species is %s \n",intNS+1,Species[intNS].chLabel);

	/*Open the files*/
	if( trace.lgTrace )
		fprintf( ioQQQ," moldat_readin opening %s:",chEFilename);

	ioLevData = open_data( chEFilename, "r" );

	nMolLevs = 0;
	intrtct = 0;
	intgrtct = 0;
	intlnct = 0;
	while(intlnct < 3)
	{
		intlnct++;
		if(read_whole_line( chLine , (int)sizeof(chLine) , ioLevData ) == NULL )
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
			cdEXIT(EXIT_FAILURE);
		}
	}
	/*Extracting out the molecular weight*/
	if(read_whole_line( chLine , (int)sizeof(chLine) , ioLevData ) == NULL )
	{
		fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
		cdEXIT(EXIT_FAILURE);
	}
	Species[intNS].fmolweight = (realnum)atof(chLine);

	/*Discard this line*/
	if(read_whole_line( chLine , (int)sizeof(chLine) , ioLevData ) == NULL )
	{
		fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
		cdEXIT(EXIT_FAILURE);
	}

	// sections of the file are separated by line that begin with "!"
	ASSERT( chLine[0] == '!' );

	/*Reading in the number of energy levels*/
	if(read_whole_line( chLine , (int)sizeof(chLine) , ioLevData ) == NULL )
	{
		fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
		cdEXIT(EXIT_FAILURE);
	}
	nMolLevs = atoi(chLine);

	long HighestIndexInFile = nMolLevs;
	/* truncate model to preset max */
	nMolLevs = MIN2( nMolLevs, MAX_NUM_LEVELS );

	if(nMolLevs <= 0)
	{
		fprintf( ioQQQ, "The number of energy levels is non-positive in datafile %s.\n", chEFilename );
		fprintf( ioQQQ, "The file must be corrupted.\n" );
		cdEXIT( EXIT_FAILURE );
	}

	Species[intNS].numLevels_max = nMolLevs;
	Species[intNS].numLevels_local = Species[intNS].numLevels_max;

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

	for( intcollindex = 0; intcollindex <NUM_COLLIDERS; intcollindex++ )
	{
		CollRatesArray[intNS][intcollindex] = NULL;
	}


	if(read_whole_line( chLine , (int)sizeof(chLine) , ioLevData ) == NULL )
	{
		fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
		cdEXIT(EXIT_FAILURE);
	}

	// sections of the file are separated by line that begin with "!"
	ASSERT( chLine[0] == '!' );

	for( long ipLev=0; ipLev<HighestIndexInFile; ipLev++)
	{	
		if(read_whole_line( chLine , (int)sizeof(chLine) , ioLevData ) == NULL )
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
			cdEXIT(EXIT_FAILURE);
		}

		// skip these levels
		if( ipLev >= nMolLevs )
			continue;

		/*information needed for label*/
		strcpy(dBaseStates[intNS][ipLev].chLabel, "    ");
		strncpy(dBaseStates[intNS][ipLev].chLabel,Species[intNS].chLabel, 4);
		dBaseStates[intNS][ipLev].chLabel[4] = '\0';
		// pad label to exactly four characters.
		if( dBaseStates[intNS][ipLev].chLabel[2]=='\0' )
		{
			dBaseStates[intNS][ipLev].chLabel[2]=' ';
			dBaseStates[intNS][ipLev].chLabel[3]=' ';
		}
		else if( dBaseStates[intNS][ipLev].chLabel[3]=='\0' )
		{
			dBaseStates[intNS][ipLev].chLabel[3]=' ';
		}

		//associate with species
		dBaseStates[intNS][ipLev].sp = &Species[intNS];

		long i = 1;
		bool lgEOL;
		long index;

		index = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
		fenergy = (realnum)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
		fstatwt = (realnum)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );

		ASSERT( index == ipLev + 1 );
		dBaseStates[intNS][ipLev].energy.set(fenergy,"cm^-1");
		dBaseStates[intNS][ipLev].g = fstatwt;

		if (ipLev > 0)
		{
			if (dBaseStates[intNS][ipLev].energy.WN() < dBaseStates[intNS][ipLev-1].energy.WN())
			{
				fprintf( ioQQQ, " The energy levels are not in order in species %s at index %li.\n",
					Species[intNS].chLabel, ipLev );
				cdEXIT(EXIT_FAILURE);
			}
		}
		if(DEBUGSTATE)
		{
			printf("The converted energy is %f \n",dBaseStates[intNS][ipLev].energy.WN());
			printf("The value of g is %f \n",dBaseStates[intNS][ipLev].g);
		}
	}

	/* fill in all transition energies, can later overwrite for specific radiative transitions */
	for( long ipHi=1; ipHi<nMolLevs; ipHi++ )
	{
		for( long ipLo=0; ipLo<ipHi; ipLo++ )
		{
			fenergyWN = (realnum)(dBaseStates[intNS][ipHi].energy.WN() - dBaseStates[intNS][ipLo].energy.WN());
			fenergyK = (realnum)(fenergyWN*T1CM);

			dBaseTrans[intNS][ipHi][ipLo].EnergyK = fenergyK;
			dBaseTrans[intNS][ipHi][ipLo].EnergyWN = fenergyWN;
			dBaseTrans[intNS][ipHi][ipLo].EnergyErg = (realnum)ERG1CM *fenergyWN;

			/* there are OH hyperfine levels where i+1 and i have exactly
			 * the same energy.  The refractive index routine will FPE with
			 * an energy of zero - so we do this test */
			if( fenergyWN>SMALLFLOAT )
				dBaseTrans[intNS][ipHi][ipLo].WLAng = (realnum)(1e+8f/fenergyWN/
					RefIndex(fenergyWN));
			else
				dBaseTrans[intNS][ipHi][ipLo].WLAng = 1e30f;
		}
	}

	if(read_whole_line( chLine , (int)sizeof(chLine) , ioLevData ) == NULL )
	{
		fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
		cdEXIT(EXIT_FAILURE);
	}
	if(chLine[0]!='!')
	{
		fprintf( ioQQQ, " The number of energy levels in file %s is not correct.\n",chEFilename);
		cdEXIT(EXIT_FAILURE);
	}
	if(read_whole_line( chLine , (int)sizeof(chLine) , ioLevData ) == NULL )
	{
		fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
		cdEXIT(EXIT_FAILURE);
	}
	intgrtct = atoi(chLine);
	/*The above gives the number of radiative transitions*/
	if(intgrtct <= 0)
	{
		fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
		cdEXIT(EXIT_FAILURE);
	}
	if(DEBUGSTATE)
	{
		printf("The number of radiative transitions is %li \n",intgrtct);
	}
	if(read_whole_line( chLine , (int)sizeof(chLine) , ioLevData ) == NULL )
	{
		fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
		cdEXIT(EXIT_FAILURE);
	}

	for( intrtct=0; intrtct<intgrtct; intrtct++)
	{
		if(read_whole_line( chLine , (int)sizeof(chLine) , ioLevData ) == NULL )
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
			cdEXIT(EXIT_FAILURE);
		}

		long i = 1;
		bool lgEOL;
		long index, ipHi, ipLo;

		index = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
		ipHi = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL ) - 1;
		ipLo = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL ) - 1;

		ASSERT( ipHi >= 0 );
		ASSERT( ipLo >= 0 );

		// skip these lines
		if( ipLo >= nMolLevs || ipHi >= nMolLevs )
			continue;

		feinsteina = (realnum)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
		/* don't need the energy in GHz, so throw it away. */
		FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
		fenergyK = (realnum)((dBaseStates[intNS][ipHi].energy.WN() -dBaseStates[intNS][ipLo].energy.WN())*T1CM);
		ASSERT( index == intrtct + 1 );

		dBaseTrans[intNS][ipHi][ipLo].Emis = AddLine2Stack(feinsteina, &dBaseTrans[intNS][ipHi][ipLo]);
		dBaseTrans[intNS][ipHi][ipLo].EnergyK = fenergyK;
		ASSERT( !isnan( dBaseTrans[intNS][ipHi][ipLo].EnergyK ) );
		fenergyWN = (realnum)((fenergyK)/T1CM);
		dBaseTrans[intNS][ipHi][ipLo].EnergyWN = fenergyWN;
		dBaseTrans[intNS][ipHi][ipLo].EnergyErg = (realnum)ERG1CM *fenergyWN;
		dBaseTrans[intNS][ipHi][ipLo].WLAng = (realnum)(1e+8/fenergyWN/RefIndex(fenergyWN));

		dBaseTrans[intNS][ipHi][ipLo].Emis->gf = (realnum)GetGF(dBaseTrans[intNS][ipHi][ipLo].Emis->Aul,
			dBaseTrans[intNS][ipHi][ipLo].EnergyWN, dBaseTrans[intNS][ipHi][ipLo].Hi->g);

		if(DEBUGSTATE)
		{
			printf("The upper level is %ld \n",ipHi+1);
			printf("The lower level is %ld \n",ipLo+1);
			printf("The Einstein A  is %E \n",dBaseTrans[intNS][ipHi][ipLo].Emis->Aul);
			printf("The Energy of the transition is %E \n",dBaseTrans[intNS][ipHi][ipLo].EnergyK);
		}
	}

	if(read_whole_line( chLine , (int)sizeof(chLine)  , ioLevData ) == NULL )
	{
		fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
		cdEXIT(EXIT_FAILURE);
	}

	if(chLine[0]!='!')
	{
		fprintf( ioQQQ, " The number of radiative transitions in file %s is not correct.\n",chEFilename);
		cdEXIT(EXIT_FAILURE);
	}

	/*Getting the number of collisional partners*/
	if(read_whole_line( chLine , (int)sizeof(chLine)  , ioLevData ) == NULL )
	{
		fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
		cdEXIT(EXIT_FAILURE);
	}
	else
	{
		intCollPart = atoi(chLine);
	}
	/*Checking the number of collisional partners does not exceed 9*/
	if(intCollPart > NUM_COLLIDERS-1)
	{
		fprintf( ioQQQ, " The number of colliders is greater than what is expected in file %s.\n", chEFilename );
		cdEXIT(EXIT_FAILURE);
	}
	/*Creating the duplicate of the number of collisional partners which is reduced*/
	intDCollPart = intCollPart;
	while(intDCollPart > 0)
	{
		if(read_whole_line( chLine , (int)sizeof(chLine)  , ioLevData ) == NULL )
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
			cdEXIT(EXIT_FAILURE);
		}

		ASSERT( chLine[0] == '!' );

		if(read_whole_line( chLine , (int)sizeof(chLine)  , ioLevData ) == NULL )
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
			cdEXIT(EXIT_FAILURE);
		}
		/*Extract out the name of the collider*/
		/*The following are the rules expected in the datafiles to extract the names of the collider*/
		/*The line which displays the species and the collider starts with a number*/
		/*This refers to the collider in the Leiden database*/
		/*In the Leiden database 1 referes to H2,2 to para-H2,3 to ortho-H2
		4 to electrons,5 to H and 6 to He*/
		chCollName = strtok(chLine," ");
		/*Leiden Collider Index*/
		intLColliderIndex = atoi(chCollName);
		/*In Cloudy,We assign the following indices for the colliders*/
		/*electron=0,proton=1,atomic hydrogen=2,He=3,He+=4,He++=5,oH2=6,pH2=7,H2=8*/

		if(intLColliderIndex == 1)
		{
			intCollIndex = 8;
		}
		else if(intLColliderIndex == 2)
		{
			intCollIndex = 7;
		}
		else if(intLColliderIndex == 3)
		{
			intCollIndex = 6;
		}
		else if(intLColliderIndex == 4)
		{
			intCollIndex = 0;
		}
		else if(intLColliderIndex == 5)
		{
			intCollIndex = 2;
		}
		else if(intLColliderIndex == 6)
		{
			intCollIndex = 3;
		}
		else
		{
			TotalInsanity();
		}
		/*This is where we allocate memory if the collider exists*/
		/*Needed to take care of the he collisions*/
		CollRatesArray[intNS][intCollIndex] = (double**)MALLOC((unsigned long)(nMolLevs) * sizeof(double*));
		for( long ipHi = 0; ipHi<nMolLevs; ipHi++ )
		{
			CollRatesArray[intNS][intCollIndex][ipHi] = (double*)MALLOC((unsigned long)(nMolLevs) * sizeof(double));
			for( long ipLo = 0; ipLo<nMolLevs; ipLo++ )
			{
				CollRatesArray[intNS][intCollIndex][ipHi][ipLo] = 0.0;
			}
		}
		if(read_whole_line( chLine , (int)sizeof(chLine)  , ioLevData ) == NULL )
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
			cdEXIT(EXIT_FAILURE);
		}

		ASSERT( chLine[0] == '!' );

		if(read_whole_line( chLine , (int)sizeof(chLine)  , ioLevData ) == NULL )
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
			cdEXIT(EXIT_FAILURE);
		}
		/*Number of coll trans*/
		intCollTran = atoi(chLine);

		if(read_whole_line( chLine , (int)sizeof(chLine)  , ioLevData ) == NULL )
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
			cdEXIT(EXIT_FAILURE);
		}

		ASSERT( chLine[0] == '!' );

		if(read_whole_line( chLine , (int)sizeof(chLine)  , ioLevData ) == NULL )
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
			cdEXIT(EXIT_FAILURE);
		}
		/*Number of coll temps*/
		intCollTemp = atoi(chLine);
		/*Storing the number of collisional temperatures*/
		AtmolCollRateCoeff[intNS][intCollIndex].ntemps = intCollTemp;
		/*Mallocing*/
		AtmolCollRateCoeff[intNS][intCollIndex].temps = 
			(double *)MALLOC((unsigned long)intCollTemp*sizeof(double));
		AtmolCollRateCoeff[intNS][intCollIndex].collrates = 
			(double***)MALLOC((unsigned long)(nMolLevs)*sizeof(double**));
		for( long ipHi=0; ipHi<nMolLevs; ipHi++ )
		{
			AtmolCollRateCoeff[intNS][intCollIndex].collrates[ipHi] = 
				(double **)MALLOC((unsigned long)(nMolLevs)*sizeof(double*));
			for( long ipLo=0; ipLo<nMolLevs; ipLo++ )
			{
				AtmolCollRateCoeff[intNS][intCollIndex].collrates[ipHi][ipLo] = 
					(double *)MALLOC((unsigned long)(intCollTemp)*sizeof(double));
				//initialize to zero
				memset( AtmolCollRateCoeff[intNS][intCollIndex].collrates[ipHi][ipLo], 0, (unsigned long)(intCollTemp)*sizeof(double) );
			}
		}

		/*Discard the header line*/
		if(read_whole_line( chLine , (int)sizeof(chLine)  , ioLevData ) == NULL )
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
			cdEXIT(EXIT_FAILURE);
		}

		ASSERT( chLine[0] == '!' );

		/*Getting the collisional Temperatures*/
		if(read_whole_line( chLine , (int)sizeof(chLine)  , ioLevData ) == NULL )
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
			cdEXIT(EXIT_FAILURE);
		}
		/*Filling the collisional temps array*/
		chColltemp =strtok(chLine," ");
		AtmolCollRateCoeff[intNS][intCollIndex].temps[0] =(realnum) atof(chColltemp);
		for( long ipTe=1; ipTe<intCollTemp; ipTe++ )
		{
			chColltemp =strtok(NULL," ");
			AtmolCollRateCoeff[intNS][intCollIndex].temps[ipTe] = atof(chColltemp);
		}
		/*Discard the header line*/
		if(read_whole_line( chLine , (int)sizeof(chLine)  , ioLevData ) == NULL )
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
			cdEXIT(EXIT_FAILURE);
		}
		/*Getting all the collisional transitions data*/
		for( nCollTrans=0; nCollTrans<intCollTran; nCollTrans++ )
		{
			if(read_whole_line( chLine , (int)sizeof(chLine)  , ioLevData ) == NULL )
			{
				fprintf( ioQQQ, " The data file %s is corrupted .\n",chEFilename);
				cdEXIT(EXIT_FAILURE);
			}

			long i = 1;
			bool lgEOL;
			long index, ipHi, ipLo;

			index = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
			ipHi = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL ) - 1;
			ipLo = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL ) - 1;

			// skip these lines
			if( ipLo >= nMolLevs || ipHi >= nMolLevs )
				continue;

			/* Indices between the very highest levels seem to be reversed */
			if( ipHi < ipLo )
			{
				ASSERT( ipLo == nMolLevs - 1);
				long temp = ipHi;
				ipHi = ipLo;
				ipLo = temp;
			}

			for( long j=0; j<intCollTemp; j++ )
			{
				AtmolCollRateCoeff[intNS][intCollIndex].collrates[ipHi][ipLo][j] = 
					(double)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
			}

			ASSERT( index == nCollTrans + 1 );

			if(DEBUGSTATE)
			{
				printf("The values of up and lo are %ld & %ld \n",ipHi,ipLo);
				printf("The collision rates are");
				for(i=0;i<intCollTemp;i++)
				{
					printf("\n %e",AtmolCollRateCoeff[intNS][intCollIndex].collrates[ipHi][ipLo][i]);
				}
				printf("\n");
			}
		}

		intDCollPart = intDCollPart -1;
	}

	fclose( ioLevData );

	return;
}
