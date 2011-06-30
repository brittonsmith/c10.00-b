/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef ELEMENTNAMES_H_
#define ELEMENTNAMES_H_

/**\file elementnames.h */
/** set of names of the chemical elements, long and short forms */
EXTERN struct t_elementnames {

	/** following used for prints in each zone, full name.
	 * the LIMELEM element is 12CO, 
	 * +1 is 13CO, +2 is H2 */
	char chElementName[LIMELM+3][11];

	/** labels for match on element name
	 * this must be caps for present logic in matches */
	char chElementNameShort[LIMELM+3][5];

	/** two letter very short form of element name, used to make
	 * emission line labels */
	char chElementSym[LIMELM+3][3];

	/** this is series of two char numbers, beginning with " 1" and
	 * ending with "31" */
	char chIonStage[LIMELM+4][3];

	/** string giving ionization stage as roman numerals */
	char chIonRoman[LIMELM+4][7];

	/* here lies the first C++ code, written by PvH, 2006 Nov 20, after
	 * the conversion from C to C++.  It is a constructor. */
	t_elementnames() {
		strncpy( chElementName[0],  "Hydrogen  ", 11 );
		strncpy( chElementName[1],  "Helium    ", 11 );
		strncpy( chElementName[2],  "Lithium   ", 11 );
		strncpy( chElementName[3],  "Beryllium ", 11 );
		strncpy( chElementName[4],  "Boron     ", 11 );
		strncpy( chElementName[5],  "Carbon    ", 11 );
		strncpy( chElementName[6],  "Nitrogen  ", 11 );
		strncpy( chElementName[7],  "Oxygen    ", 11 );
		strncpy( chElementName[8],  "Fluorine  ", 11 );
		strncpy( chElementName[9],  "Neon      ", 11 );
		strncpy( chElementName[10], "Sodium    ", 11 );
		strncpy( chElementName[11], "Magnesium ", 11 );
		strncpy( chElementName[12], "Aluminium ", 11 );
		strncpy( chElementName[13], "Silicon   ", 11 );
		strncpy( chElementName[14], "Phosphorus", 11 );
		strncpy( chElementName[15], "Sulphur   ", 11 );
		strncpy( chElementName[16], "Chlorine  ", 11 );
		strncpy( chElementName[17], "Argon     ", 11 );
		strncpy( chElementName[18], "Potassium ", 11 );
		strncpy( chElementName[19], "Calcium   ", 11 );
		strncpy( chElementName[20], "Scandium  ", 11 );
		strncpy( chElementName[21], "Titanium  ", 11 );
		strncpy( chElementName[22], "Vanadium  ", 11 );
		strncpy( chElementName[23], "Chromium  ", 11 );
		strncpy( chElementName[24], "Manganese ", 11 );
		strncpy( chElementName[25], "Iron      ", 11 );
		strncpy( chElementName[26], "Cobalt    ", 11 );
		strncpy( chElementName[27], "Nickel    ", 11 );
		strncpy( chElementName[28], "Copper    ", 11 );
		strncpy( chElementName[29], "Zinc      ", 11 );
		strncpy( chElementName[30], "C12O18    ", 11 );
		strncpy( chElementName[31], "C13O18    ", 11 );
		strncpy( chElementName[32], "H2        ", 11 );

		strncpy( chElementNameShort[0],  "HYDR", 5 );
		strncpy( chElementNameShort[1],  "HELI", 5 );
		strncpy( chElementNameShort[2],  "LITH", 5 );
		strncpy( chElementNameShort[3],  "BERY", 5 );
		strncpy( chElementNameShort[4],  "BORO", 5 );
		strncpy( chElementNameShort[5],  "CARB", 5 );
		strncpy( chElementNameShort[6],  "NITR", 5 );
		strncpy( chElementNameShort[7],  "OXYG", 5 );
		strncpy( chElementNameShort[8],  "FLUO", 5 );
		strncpy( chElementNameShort[9],  "NEON", 5 );
		strncpy( chElementNameShort[10], "SODI", 5 );
		strncpy( chElementNameShort[11], "MAGN", 5 );
		strncpy( chElementNameShort[12], "ALUM", 5 );
		strncpy( chElementNameShort[13], "SILI", 5 );
		strncpy( chElementNameShort[14], "PHOS", 5 );
		strncpy( chElementNameShort[15], "SULP", 5 );
		strncpy( chElementNameShort[16], "CHLO", 5 );
		strncpy( chElementNameShort[17], "ARGO", 5 );
		strncpy( chElementNameShort[18], "POTA", 5 );
		strncpy( chElementNameShort[19], "CALC", 5 );
		strncpy( chElementNameShort[20], "SCAN", 5 );
		strncpy( chElementNameShort[21], "TITA", 5 );
		strncpy( chElementNameShort[22], "VANA", 5 );
		strncpy( chElementNameShort[23], "CHRO", 5 );
		strncpy( chElementNameShort[24], "MANG", 5 );
		strncpy( chElementNameShort[25], "IRON", 5 );
		strncpy( chElementNameShort[26], "COBA", 5 );
		strncpy( chElementNameShort[27], "NICK", 5 );
		strncpy( chElementNameShort[28], "COPP", 5 );
		strncpy( chElementNameShort[29], "ZINC", 5 );
		strncpy( chElementNameShort[30], "12CO", 5 );
		strncpy( chElementNameShort[31], "13CO", 5 );
		strncpy( chElementNameShort[32], "H2  ", 5 );

		strncpy( chElementSym[0],  "H ", 3 );
		strncpy( chElementSym[1],  "He", 3 );
		strncpy( chElementSym[2],  "Li", 3 );
		strncpy( chElementSym[3],  "Be", 3 );
		strncpy( chElementSym[4],  "B ", 3 );
		strncpy( chElementSym[5],  "C ", 3 );
		strncpy( chElementSym[6],  "N ", 3 );
		strncpy( chElementSym[7],  "O ", 3 );
		strncpy( chElementSym[8],  "F ", 3 );
		strncpy( chElementSym[9],  "Ne", 3 );
		strncpy( chElementSym[10], "Na", 3 );
		strncpy( chElementSym[11], "Mg", 3 );
		strncpy( chElementSym[12], "Al", 3 );
		strncpy( chElementSym[13], "Si", 3 );
		strncpy( chElementSym[14], "P ", 3 );
		strncpy( chElementSym[15], "S ", 3 );
		strncpy( chElementSym[16], "Cl", 3 );
		strncpy( chElementSym[17], "Ar", 3 );
		strncpy( chElementSym[18], "K ", 3 );
		strncpy( chElementSym[19], "Ca", 3 );
		strncpy( chElementSym[20], "Sc", 3 );
		strncpy( chElementSym[21], "Ti", 3 );
		strncpy( chElementSym[22], "V ", 3 );
		strncpy( chElementSym[23], "Cr", 3 );
		strncpy( chElementSym[24], "Mn", 3 );
		strncpy( chElementSym[25], "Fe", 3 );
		strncpy( chElementSym[26], "Co", 3 );
		strncpy( chElementSym[27], "Ni", 3 );
		strncpy( chElementSym[28], "Cu", 3 );
		strncpy( chElementSym[29], "Zn", 3 );
		strncpy( chElementSym[30], "12", 3 );
		strncpy( chElementSym[31], "13", 3 );
		strncpy( chElementSym[32], "H2", 3 );

		strncpy( chIonStage[0],  " 1", 3 );
		strncpy( chIonStage[1],  " 2", 3 );
		strncpy( chIonStage[2],  " 3", 3 );
		strncpy( chIonStage[3],  " 4", 3 );
		strncpy( chIonStage[4],  " 5", 3 );
		strncpy( chIonStage[5],  " 6", 3 );
		strncpy( chIonStage[6],  " 7", 3 );
		strncpy( chIonStage[7],  " 8", 3 );
		strncpy( chIonStage[8],  " 9", 3 );
		strncpy( chIonStage[9],  "10", 3 );
		strncpy( chIonStage[10], "11", 3 );
		strncpy( chIonStage[11], "12", 3 );
		strncpy( chIonStage[12], "13", 3 );
		strncpy( chIonStage[13], "14", 3 );
		strncpy( chIonStage[14], "15", 3 );
		strncpy( chIonStage[15], "16", 3 );
		strncpy( chIonStage[16], "17", 3 );
		strncpy( chIonStage[17], "18", 3 );
		strncpy( chIonStage[18], "19", 3 );
		strncpy( chIonStage[19], "20", 3 );
		strncpy( chIonStage[20], "21", 3 );
		strncpy( chIonStage[21], "22", 3 );
		strncpy( chIonStage[22], "23", 3 );
		strncpy( chIonStage[23], "24", 3 );
		strncpy( chIonStage[24], "25", 3 );
		strncpy( chIonStage[25], "26", 3 );
		strncpy( chIonStage[26], "27", 3 );
		strncpy( chIonStage[27], "28", 3 );
		strncpy( chIonStage[28], "29", 3 );
		strncpy( chIonStage[29], "30", 3 );
		strncpy( chIonStage[30], "31", 3 );
		/* this is special for molecule */
		strncpy( chIonStage[31], "CO", 3 );
		strncpy( chIonStage[32], "CO", 3 );
		strncpy( chIonStage[33], "H2", 3 );

		strncpy( chIonRoman[0],  "I", 7 );
		strncpy( chIonRoman[1],  "II", 7 );
		strncpy( chIonRoman[2],  "III", 7 );
		strncpy( chIonRoman[3],  "IV", 7 );
		strncpy( chIonRoman[4],  "V", 7 );
		strncpy( chIonRoman[5],  "VI", 7 );
		strncpy( chIonRoman[6],  "VII", 7 );
		strncpy( chIonRoman[7],  "VIII", 7 );
		strncpy( chIonRoman[8],  "IX", 7 );
		strncpy( chIonRoman[9],  "X", 7 );
		strncpy( chIonRoman[10], "XI", 7 );
		strncpy( chIonRoman[11], "XII", 7 );
		strncpy( chIonRoman[12], "XIII", 7 );
		strncpy( chIonRoman[13], "XIV", 7 );
		strncpy( chIonRoman[14], "XV", 7 );
		strncpy( chIonRoman[15], "XVI", 7 );
		strncpy( chIonRoman[16], "XVII", 7 );
		strncpy( chIonRoman[17], "XVIII", 7 );
		strncpy( chIonRoman[18], "XIX", 7 );
		strncpy( chIonRoman[19], "XX", 7 );
		strncpy( chIonRoman[20], "XXI", 7 );
		strncpy( chIonRoman[21], "XXII", 7 );
		strncpy( chIonRoman[22], "XXIII", 7 );
		strncpy( chIonRoman[23], "XXIV", 7 );
		strncpy( chIonRoman[24], "XXV", 7 );
		strncpy( chIonRoman[25], "XXVI", 7 );
		strncpy( chIonRoman[26], "XXVII", 7 );
		strncpy( chIonRoman[27], "XXVIII", 7 );
		strncpy( chIonRoman[28], "XXIX", 7 );
		strncpy( chIonRoman[29], "XXX", 7 );
		strncpy( chIonRoman[30], "XXXI", 7 );
		/* this is special for molecule */
		strncpy( chIonRoman[31], "  ", 7 );
		strncpy( chIonRoman[32], "  ", 7 );
		strncpy( chIonRoman[33], "  ", 7 );
	};

} elementnames;

#endif /* ELEMENTNAMES_H_ */
