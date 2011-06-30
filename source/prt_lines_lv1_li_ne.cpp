/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*lines_lv1_li_ne place lines of elements lithium through neon into lines storage stack */
/*GetLineRec return recombination coefficient*hnu*eden*n_ion for C, N, or O recombination lines from Dima's list,
 * also zero's line in master stack so not entered second time in later dump of all rec lines */
#include "cddefines.h"
#include "carb.h"
#include "nitro.h"
#include "oxy.h"
#include "coolheavy.h"
#include "atmdat.h"
#include "doppvel.h"
#include "ionbal.h"
#include "dense.h"
#include "phycon.h"
#include "physconst.h"
#include "atoms.h"
#include "mole.h"
#include "embesq.h"
#include "taulines.h"
#include "trace.h"
#include "lines_service.h"
#include "lines.h"
STATIC double GetLineRec(
	/* this is the number of the emission line in the stack of lines, on the C scale */
	long int ip, 
	/* the multiplet wavelength */
  long int lWl);

void lines_lv1_li_ne(void)
{
	long int ipnt;
	double 
	  chem ,
	  corr, 
	  ct4363, 
	  ctRate, 
	  efac, 
	  effec, 
	  efficn2, 
	  fac, 
	  HBeta ,
	  p386, 
	  pump,
	  r4363, 
	  r6584, 
	  raten3, 
	  rb, 
	  rec, 
	  rn3mor, 
	  rn3tot, 
	  rnii, 
	  rp300, 
	  rp386, 
	  r12 ,
	  r13 ,
	  sum,
		rate_OH_dissoc;
	double rec7323 , rec7332, rec3730 , rec3726 , rec2471,
		reco23tot , reco22tot;

	DEBUG_ENTRY( "lines_lv1_li_ne()" );

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "   lines_lv1_li_ne called\n" );
	}

	/* level 1 ines */
	ipnt = StuffComment( "level 1 lines" );
	linadd( 0., (realnum)ipnt , "####", 'i',
		" start level 1 ines" );

	linadd(CoolHeavy.colmet,0,"Mion",'c',
		" cooling due to collisional ionization of heavy elements" );

	/* lithium */

	/* beryllium */

	/* boron*/

	/* carbon*/

	/* >>chng 97 may 02, added better rec coefficient
	 * C I 1656 recombination, all agents */
	rec = GetLineRec(3,1657)*(1.-TauLines[ipT1656].Emis->ColOvTot);
	rec = GetLineRec(3,1657)*emit_frac(&TauLines[ipT1656]);
	PutExtra(rec);

	PutLine(&TauLines[ipT1656],
		" C 1 1656, collision strength from van Regemoter");

	linadd(rec,1656,"REC ",'i',
		" C 1 1656 recomb; n.b. coll deexcitation not in" );

	PntForLine(9850.,"C  1",&ipnt);
	lindst(carb.c9850,9850,"C Ic",ipnt,'c',true,
		" C 1 9850, coll excit" );

	/* >>chng 97 may 02, added better rec coefficient
	 * C 1 9850, recombination contribution rec coefficient from 
	 * >>refer	C1	rec	Escalante, Vladimir, & Victor, G.A., 1990, ApJS 73, 513.
	 * r9850 is correction for collisional deexcitation as in carb cool
	 * >>chng 97 aug 2, had factor of rec, changed to r9850, this
	 * was a big mistake */
	rec = (GetLineRec(1,9088) + GetLineRec(2,9658))*carb.r9850;
	lindst(rec,9850,"C Ir",ipnt,'r',true,
		" C I 9850 recombination contribution" );

	// total intensity added as information line but in outward
	// beam since cooling done above
	linadd(rec+carb.c9850,9850,"TOTL",'i',
		" total intensity, all processes, C I 9850");


	PntForLine(8727.,"C  1",&ipnt);
	lindst(carb.c8727,8727,"C  1",ipnt,'c',true,
		"C 1 8727; equivalent to 4363" );

	linadd(carb.c8727*1.22e-6,4621,"C  1",'c',
		" 1S - 3P" );

	PutLine(&TauLines[ipT610],
		"  C 1 610 micron ");

	PutLine(&TauLines[ipT370],
		"  C 1 370 micron ");

	PutLine(&TauLines[ipT157],
		"  C 2 158 micron, both e- and H0, H2 in excitation ");

	linadd(
		TauLines[ipC2_2325].Emis->xIntensity+
		TauLines[ipC2_2324].Emis->xIntensity+
		TauLines[ipC2_2329].Emis->xIntensity+
		TauLines[ipC2_2328].Emis->xIntensity+
		TauLines[ipC2_2327].Emis->xIntensity+carb.p2326*8.6e-12,
		2326,"TOTL",'i',
		" total intensity of C II] 2326, all lines in the multiplet " );
	PutLine(&TauLines[ipC2_2325],
		" ");
	PutLine(&TauLines[ipC2_2324],
		" ");
	PutLine(&TauLines[ipC2_2329],
		" ");
	PutLine(&TauLines[ipC2_2328],
		" ");
	PutLine(&TauLines[ipC2_2327],
		" ");

	linadd(carb.p2326*8.6e-12,2326,"Phot",'i' ,
		" photoproduction, Helfand and Trefftz");

	/* >>chng 97 may 02, better rec coef */
	/*rec = GetLineRec(11,1335)*(1.-TauLines[ipT1335].ColOvTot);*/
	/* >>chng 02 jul 01, add function to return emission probability */
	rec = GetLineRec(11,1335)*emit_frac(&TauLines[ipT1335]);

	/* total intensity of C 2 1335 */
	PutExtra(MAX2(0.,rec));
	PutLine(&TauLines[ipT1335],
		" total intensity of C 2 1335");

	linadd(rec,1335,"REC ",'i',
		" C 2 1335 recombination," );

	/* the CII 3918.98/3920.68 and 6578.05/6582.88 multiplets,
	 * contributions by both continuum pumping through XUV line
	 * and recombination */
	/* this is the driving line, pump is photons cm^-3 s^-1 */
	if( nWindLine > 0 )
	{
		pump = TauLine2[186].Emis->pump*TauLine2[186].Emis->PopOpc;
	}
	else
	{
		pump = 0.;
	}


	// pumped 3920, count as recomb since does remove energy
	PntForLine(3920.,"C  2",&ipnt);
	lindst(pump*0.387 * 5.08e-12/(1.+dense.eden/1e12) ,3920,"pump",ipnt,'r',true ,
		"  CII 3918.98/3920.68 is only pumped, no recombination part");

	/* recombination and specific pump for CII 5684 */
	rec = GetLineRec(8, 6580 );
	/* convert UV pump rate to intensity with branching ratio and hnu */
	pump *= 0.305 * 0.387 * 3.02e-12;
	linadd(rec/(1.+dense.eden/1e12),6580,"C 2r",'i',
		" recombination part of C II 6580 line " );
	linadd(pump/(1.+dense.eden/1e12),6580,"C 2p",'i',
		" pumped part of line C II 6580" );

	// pumped 3920, count as recomb since does remove energy
	PntForLine(6580.,"C  2",&ipnt);
	lindst((rec+pump)/(1.+dense.eden/1e12),6580,"TOTL",ipnt,'r',true ,
		" total intensity, all processes, C II 6580");

	/* C 3 977
	 * recombination contribution from nussbaumer and story 84 */
	/*rec = GetLineRec(179,977)*(1.-TauLines[ipT977].ColOvTot);*/
	/* >>chng 02 jul 01, add function to compute emission fraction */
	rec = GetLineRec(179,977)*emit_frac(&TauLines[ipT977]);

	/* continuum pumped C 3 977 by continuum near 386A */
	rp386 = TauLines[ipT386].Emis->pump*TauLines[ipT386].Emis->PopOpc;

	/* higher lines, same process */
	rp300 = TauLines[ipT310].Emis->pump*TauLines[ipT310].Emis->PopOpc + 
		TauLines[ipT291].Emis->pump*TauLines[ipT291].Emis->PopOpc + 
		TauLines[ipT280].Emis->pump*TauLines[ipT280].Emis->PopOpc + 
		TauLines[ipT274].Emis->pump*TauLines[ipT274].Emis->PopOpc + 
		TauLines[ipT270].Emis->pump*TauLines[ipT270].Emis->PopOpc;

	/* total line intensity due to pumping */
	/*p386 = (rp386 + rp300)*2.03e-11*(1.-TauLines[ipT977].ColOvTot);*/
	/* >>chng 02 jul 01, add function to compute emission fraction */
	p386 = (rp386 + rp300)*2.03e-11*emit_frac(&TauLines[ipT977]);

	/* total C 3 977 including recombination and pumping */
	PutExtra(p386+MAX2(0.,rec));

	PutLine(&TauLines[ipT977],
		"  total C 3] 977, recombination + collisional + pumped excitation ");

	linadd(rec,977,"C3 R",'i',
		" dielectronic recombination contribution to C 3 977 " );

	linadd(p386,977,"P386",'r',
		" C 3 977 pumped by continuum near 386A" );

	/* C 3  1909 collision, both lines together */
	fac = embesq.em1908 + TauLines[ipT1909].Emis->xIntensity;
	lindst(fac,1909,"TOTL",TauLines[ipT1909].ipCont,'i',false,
		"C 3  1909 collision, both lines together");

	// the faster line, which might become optically thick
	PutLine(&TauLines[ipT1909],
		"C III 19091");

	// the very slow transition
	PntForLine(1907.,"C  3",&ipnt);
	lindst(embesq.em1908,1907,"C  3",ipnt,'c',true,
		" C 3 1908 j-2 to ground" );

	lindst(embesq.em13C1910,1910,"13C3",ipnt,'c',true,
		" the 13C forbidden line of C III " );

	/*corr = 1.-TauLines[ipT1909].ColOvTot;*/
	/* >>chng 02 jul 01, add function to compute emission fraction */
	corr = emit_frac(&TauLines[ipT1909]);
	fac = dense.eden*dense.xIonDense[ipCARBON][3]/(phycon.te/phycon.te10)*
	  ionbal.DielSupprs[1][2];

	linadd(3.1e-19*fac*corr,1909,"C3 R",'i',
		" C 3 1909 recombination from Storey" );

	linadd(carb.p1909*corr*1.05e-11,1909,"Phot",'i',
		" C 3 1909 following relax following inner shell photoionization" );

	/* >>chng 97 may 02, better rec ocef */
	/*rec = GetLineRec(178,1176)*(1.-TauLines[ipc31175].ColOvTot);*/
	/* >>chng 02 jul 01, add function to compute emission fraction */
	rec = GetLineRec(178,1176)*emit_frac(&TauLines[ipc31175]);
	PutExtra(MAX2(0.,rec));

	PutLine(&TauLines[ipc31175],
		"  C 3* 1175, excited state line, above 1909 ");

	linadd(MAX2(0.,rec),1175,"Rec ",'i',
		" dielectronic recombination contribution to C 3 1175 " );

	/* recombination C 4 1549 from C 5
	 * >>chng 97 may 02, better rec coef */
	/*rec = GetLineRec(25,1549)*(1.-TauLines[ipT1550].ColOvTot);*/
	/* >>chng 02 jul 01, add function to compute emission fraction */
	rec = GetLineRec(25,1549)*emit_frac(&TauLines[ipT1550]);

	linadd(
		TauLines[ipT1550].Emis->xIntensity+
		TauLines[ipT1548].Emis->xIntensity+
		rec,1549,"TOTL",'i',"total intensity of C 4 1549, all processes " );

	sum = 
		TauLines[ipT1550].Emis->xIntensity*TauLines[ipT1550].Emis->FracInwd + 
		TauLines[ipT1548].Emis->xIntensity*TauLines[ipT1548].Emis->FracInwd;

	linadd(sum+rec*TauLines[ipT1550].Emis->FracInwd,1549,"Inwd",'i',
		"inward part of C 4 " );

	PutExtra(rec*.3333);
	PutLine(&TauLines[ipT1550],
		" ");

	PutExtra(rec*.6666);
	PutLine(&TauLines[ipT1548],
		" ");

	linadd((TauLines[ipT1550].Emis->ots+
		TauLines[ipT1548].Emis->ots)*TauLines[ipT1548].EnergyErg,
	  1549,"DEST",'i',
	  " part of line destroyed by photoionization of Balmer continuum " );

	linadd(rec,1549,"C4 r",'i',
		" recombination C 4 1549 from CV" );

	PutLine(&TauLines[ipT312],
		"  Li seq 2s 3p Li seq transition");

	/*************************/

	/* nitrogen */

	// these are the FUV lines that pump [N I]
	PutLine(&TauLines[ipNI_pumpIndirect],"FUV indirect excitation");
	for( int i=0; i < NI_NDP; ++i )
		PutLine(&TauLines[ipNI_pumpDirect[i]],"FUV direct excitation");

	PutLine(&TauLines[ipT1200]," N I 1200, full multiplet, all processes");

	PntForLine(3466.52,"N  1",&ipnt);
	lindst( nitro.xN3466, 3466, "N  1", ipnt, 'c', true, " [N I] 3466.50 only" );

	lindst( nitro.xN3467, 3468, "N  1", ipnt, 'c', true, " [N I] 3466.54 only" );

	// doublet together
	lindst( nitro.xN3467+nitro.xN3466, 3467, "TOTL", ipnt, 'i', false,
		" [N I] 3466.54, 3466.50 together" );

	PntForLine(5199.,"N  1",&ipnt);
	/**** Terry's addition **************/
	lindst( nitro.rec5199, 5199, "TOTr", ipnt, 'r', true,
		" estimate of contribution to [N I] 5199 by recombination" );

	/* this is upper limit to production of 5200 by chemistry - assume every photo dissociation
	 * populates upper level 
	 * co.nitro_dissoc_rate is the total N photo dissociation rate, cm-3 s-1 */
	// count as recombination since removes energy but is not coolng
	chem = co.nitro_dissoc_rate * 3.82e-12 * nitro.quench_5200;
	lindst( chem, 5199, "chem", ipnt, 'r', true,
		" upper limit to [N I] 5199 produced by chemistry" );

	/* this is upper limit to production of 5200 by charge transfer -  
	 * atmdat.HCharExcRecTo_N0_2D is the rate coefficient (cm3 s-1) for N+(3P) + H0 -> H+ + N0(2D) */
	ctRate = atmdat.HCharExcRecTo_N0_2D*dense.xIonDense[ipHYDROGEN][0]*dense.xIonDense[ipNITROGEN][1] * 
		3.82e-12 * nitro.quench_5200;
	lindst( ctRate, 5199, "H CT", ipnt, 'r', true,
		" upper limit to [N I] 5199 produced by charge transfer" );

	// estimate of pumping contribution to 5200 doublet
	lindst( nitro.pump5199, 5199, "pump", ipnt, 'r', true,
		" estimate of contribution to [N I] 5199 by FUV pumping" );

	/*  >>chng 06 jul 08, add chem process to total assuming in proportion to stat weight  */
	lindst( nitro.xN5200 + (ctRate+chem)*0.6, 5200, "N  1", ipnt, 'c', true,
		" [N I] 5200 - all processes - stat weight is 6 - total in term is 10" );

	lindst( nitro.xN5198 + (ctRate+chem)*0.4, 5198, "N  1", ipnt, 'c', true,
		" [N I] 5198 - all processes - stat weight is 4 - total in term is 10" );

	// doublet together
	lindst( nitro.xN5200+nitro.xN5198 + (ctRate+chem), 5199, "TOTL", ipnt, 'i', false,
		" [N I] 5200 + 5198 together" );

	PntForLine(10403.,"N  1",&ipnt);
	lindst( nitro.xN10397, 10397, "N  1", ipnt, 'c', true, " [N I] 10397.7 only");
	lindst( nitro.xN10398, 10398, "N  1", ipnt, 'c', true, " [N I] 10398.2 only");
	lindst( nitro.xN10407, 10407, "N  1", ipnt, 'c', true, " [N I] 10407.2 only");
	lindst( nitro.xN10408, 10408, "N  1", ipnt, 'c', true, " [N I] 10407.6 only");

	// entire multiplet together
	lindst( nitro.xN10398+nitro.xN10397+nitro.xN10408+nitro.xN10407, 10403, "TOTL", ipnt, 'i', false,
		" [N I] 10398.2, 10397.7, 10407.6, 10407.2 together");

	/*************************/

	PutLine(&TauLines[ipT671] ,"  6 lines with fake collision strengths" );
	PutLine(&TauLines[ipT315] ,"  6 lines with fake collision strengths" );
	PutLine(&TauLines[ipT333] ,"  6 lines with fake collision strengths" );
	PutLine(&TauLines[ipT324] ,"  6 lines with fake collision strengths" );
	PutLine(&TauLines[ipT374g],"  6 lines with fake collision strengths" );
	PutLine(&TauLines[ipT374x],"  6 lines with fake collision strengths" );

	PntForLine(6584.,"N  2",&ipnt);
	lindst(nitro.c6584/(1.+1./2.951),6584,"N  2",ipnt,'c',true ,
		"  N 2 6584 alone " );

	PntForLine(6548.,"N  2",&ipnt);
	lindst(nitro.c6584/(1.+2.951),6548,"N  2",ipnt,'c',true,
		"  N 2 6548 alone "  );

	efficn2 = 4e-3/(4e-3 + 5.18e-6*dense.eden/phycon.sqrte);
	r6584 = 8e-22/(phycon.te70/phycon.te03/phycon.te03)*efficn2;
	linadd(r6584*dense.xIonDense[ipNITROGEN][2]*dense.eden,6584,"REC ",  'i',
		"  N 2 6584 alone, recombination contribution"  );

	/* helium charge transfer from 
	>>refer	n2	CT	Sun Sadeghpour, Kirby Dalgarno and Lafyatis, CfA preprint 4208 */
	ctRate = 1.8e-11*dense.xIonDense[ipHELIUM][0]*dense.xIonDense[ipNITROGEN][2]*1.146/(1.146 + 
	  0.87*dense.cdsqte)*3.46e-12;

	/* >>chng 01 jul 09, add recombination contribution to 5755 */
	/* >>refer	n2	rec	Liu, X.W., Storey, P.J., Barlow, M.J., Danziger, I.J., Cohen, M.,
	 * >>refercon	& Bryce, M., 2000, MNRAS, 312, 585 */
	/* they give intensity in terms of hbeta intensity as their equation 1 */
	if( dense.xIonDense[ipHYDROGEN][1] > SMALLFLOAT )
	{
		/* this test on >0 is necessary because for sims with no H-ionizing radiation
		 * the H+ density is initially zero */
		/* H beta recombination, assuming old case B, needed since HS tables have
		 * only a narrow temperature range - at this point units are ergs cm^3 s-1 */
		HBeta = (pow(10.,-20.89 - 0.10612*POW2(phycon.alogte - 4.4)))/phycon.te;

		/* now convert to ergs cm-3 s-1 
		 * >>chng 05 mar 17, this step was missing, so recombination intensity off by density squared,
		 * bug reported by Marcelo Castellanos */
		HBeta *= dense.eden * dense.xIonDense[ipHYDROGEN][1];

		/* CoolHeavy.xN2_A3_tot is fraction of excitations that produce a photon 
		 * and represents the correction for collisional deexcitation */
		/*>>chng 05 dec 16, Liu et al. (2000) eqn 1 uses t = Te/10^4 K, not Te so phycon.te30 
		 * is too large: (10^4)^0.3 = 16 - div by 15.8489 - bug caught by Kevin Blagrave */
		rec = nitro.xN2_A3_tot * HBeta *
			3.19 * phycon.te30 / 15.84893 * dense.xIonDense[ipNITROGEN][2]/dense.xIonDense[ipHYDROGEN][1];
	}
	else
	{
		HBeta = 0.;
		rec = 0.;
	}

	linadd(nitro.c5755+ctRate+rec ,5755,"N  2",'i',
		" N 2 5755  total, collisions plus charge transfer plus recombination"  );

	PntForLine(5755.,"N  2",&ipnt);
	lindst(nitro.c5755,5755,"Coll",ipnt,'c',true,
		" N 2 5755  collisional contribution" );

	lindst(ctRate,5755,"C T ",ipnt,'r',true,
		" N 2 5755  charge transfer contribution " );

	lindst( rec ,5755,"N 2r",ipnt,'r',true,
		" N 2 5755  recombination contribution" );

	PutLine(&TauLines[ipT122],
		"  N 2 fine structure line ");

	PutLine(&TauLines[ipT205],
		"  N 2 fine structure line " );

	PutLine(&TauLines[ipT2140],
		"  N 2 2140 intercombination line " );

	/* >>chng 97 may 02, better rec contribution */
	/*rec = GetLineRec(201,1085)*(1.-TauLines[ipT1085].ColOvTot);*/
	/* >>chng 02 jul 01, add function to compute emission fraction */
	rec = GetLineRec(201,1085)*emit_frac(&TauLines[ipT1085]);
	PutExtra(MAX2(0.,rec));

	PutLine(&TauLines[ipT1085],
		"  N 2 1084, CS guess from g-bar " );

	linadd(MAX2(0.,rec),1085,"Rec ",'i',
		" dielectronic recombination contribution to N 2 1085" );

	/* continuum pumping of N 2 intersystem transition */
	rnii = TauLines[ipT671].Emis->pump*TauLines[ipT671].Emis->PopOpc;

	linadd(rnii*0.377*0.75*3.02e-12*efficn2,6584,"N2cn",'i',
		" continuum pumped N 2 6584 " );

	efficn2 = 1./(1. + COLL_CONST*dense.eden/phycon.sqrte);
	linadd(rnii*0.0117*3.46e-12*efficn2,5755,"N2cn",'i',
		" continuum pumped N 2 5755" );

	/* pumping of the NII 509A line excites 2p4s ^3P^o,
	 * which decays through the 3328, 5679, and 671 multiplets
	 * the NII 3311 (6 lines), 3840 (6 lines), and 3600 (3 lines) multiplets,
	 * contributions by both continuum pumping through XUV line
	 * and recombination */
	/* this is the driving line, pump is photons cm^-3 s^-1 */
	if( nWindLine > 0 )
	{
		pump = TauLine2[265].Emis->pump*TauLine2[265].Emis->PopOpc;
	}
	else
	{
		pump = 0.;
	}

	PntForLine(3311.,"N  2",&ipnt);
	lindst(pump*0.236 * 6.01e-12/(1.+dense.eden/1e12) ,3311,"pump",ipnt,'r',true,
		" NII 3311.42 - 3331.31 (6 lines) are only pumped, no recombination part" );

	PntForLine(3840.,"N  2",&ipnt);
	lindst(pump*0.186 * 5.18e-12/(1.+dense.eden/1e12) ,3840,"pump",ipnt,'r',true,
		" NII 3829.8-3856.06 (6 lines) are only pumped, no recombination part" );

	PntForLine(3609.,"N  2",&ipnt);
	lindst(pump*0.025 * 5.52e-12/(1.+dense.eden/1e12) ,3609,"pump",ipnt,'r',true,
		" NII 3593.60/3609.1/3615.86 (3 lines) are only pumped, no recombination part" );


	PntForLine(4640.,"N  2",&ipnt);
	lindst(pump*0.186*0.595 * 4.31e-12/(1.+dense.eden/1e12) ,4640,"pump",ipnt,'r',true ,
		" NII 4601.5-4643.1 (6 lines) are only pumped, no recombination part");

	PntForLine(5010.,"N  2",&ipnt);
	lindst(pump*0.025*0.442 * 3.97e-12/(1.+dense.eden/1e12) ,5010,"pump",ipnt,'r',true,
		" NII 5002.7/5010.6/5045.1 (3 lines) are only pumped, no recombination part " );

	/* recombination and specific pump for NII 5679 */
	rec = GetLineRec(44, 5679 );
	/* convert UV pump rate to intensity with branching ratio and hnu */
	pump *= 0.236 * 0.626 * 3.50e-12;

	linadd(rec/(1.+dense.eden/1e12),5679,"N 2r",'i',
		" recombination part of N II 5679 line" );

	linadd(pump/(1.+dense.eden/1e12),5679,"N 2p",'i',
		" pumped part of line N II 5679 " );

	PntForLine(5679.,"N  2",&ipnt);
	lindst((rec+pump)/(1.+dense.eden/1e12),5679,"TOTL",ipnt,'r',true,
		" total intensity, all processes, N II 5679 " );

	PutLine(&TauLines[ipT57],
		" [N 3] 57 micron fine structure line");

	linadd(
		TauLines[ipN3_1749].Emis->xIntensity+
		TauLines[ipN3_1747].Emis->xIntensity+
		TauLines[ipN3_1754].Emis->xIntensity+
		TauLines[ipN3_1752].Emis->xIntensity+
		TauLines[ipN3_1751].Emis->xIntensity,
		1750,"TOTL",'i',
		" total intensity of N III] 1750, all lines in the multiplet " );
	PutLine(&TauLines[ipN3_1749],
		" ");
	PutLine(&TauLines[ipN3_1747],
		" ");
	PutLine(&TauLines[ipN3_1754],
		" ");
	PutLine(&TauLines[ipN3_1752],
		" ");
	PutLine(&TauLines[ipN3_1751],
		" ");

	/* continuum pumped "Bowen" N 3
	 * rate system a is populated  */
	raten3 = TauLines[ipT374x].Emis->PopOpc*TauLines[ipT374x].Emis->pump;

	/* rate system b is populated */
	if( DoppVel.TurbVel < 200. )
	{
		rb = TauLines[ipT374x].Emis->PopOpc*TauLines[ipT374x].Emis->pump + 
			TauLines[ipT374g].Emis->PopOpc*TauLines[ipT374g].Emis->pump;
	}
	else
	{
		/* only one line if both fully overlap due to large turb */
		rb = TauLines[ipT374g].Emis->PopOpc*TauLines[ipT374g].Emis->pump;
	}

	rn3mor = 
		TauLines[ipT315].Emis->PopOpc*TauLines[ipT315].Emis->pump*0.448 + 
		TauLines[ipT324].Emis->PopOpc*TauLines[ipT324].Emis->pump*0.78 + 
		TauLines[ipT333].Emis->PopOpc*TauLines[ipT333].Emis->pump*0.434;

	/* pumping of optical N 3 bowen lines */
	rn3tot = (rb + raten3)*0.439 + rn3mor;

	sum = raten3*4.29e-12;
	PntForLine(4640.,"N3cn",&ipnt);
	lindst(sum,4640,"N3cn",ipnt,'r',true ,
		" continuum pumped \"Bowen\" N 3, optically thin excited line ");

	/*  */
	sum = rb*4.29e-12*0.834;
	PntForLine(4634.,"N3cn",&ipnt);
	lindst(sum,4634,"N3cn",ipnt,'r',true ,
		" continuum pumped \"Bowen\" N 3, optically thin excited line");

	sum = rb*4.29e-12*(1. - 0.834);
	PntForLine(4642.,"N3cn",&ipnt);
	lindst(sum,4642,"N3cn",ipnt,'r',true,
		" continuum pumped \"Bowen\" N 3, optically thin excited line" );

	/* total rate for N 3 990
	 * correction factor for collisional deexcitation */
	/*fac = 1.-TauLines[ipT990].ColOvTot;*/
	/* >>chng 02 jul 01, add function to compute emission fraction */
	fac = 1.-emit_frac(&TauLines[ipT990]);

	/* >>chng 97 may 02, better rec coef */
	/*rec = GetLineRec(216,991)*(1.-TauLines[ipT990].ColOvTot);*/
	/* >>chng 02 jul 01, add function to compute emission fraction */
	rec = GetLineRec(216,991)*emit_frac(&TauLines[ipT990] );
	PutExtra(MAX2(0.,rec)+rn3tot*2.01e-11*fac);

	PutLine(&TauLines[ipT990],
		"  N 3 990, all processes ");

	linadd(rec+rn3tot*2.01e-11*fac,990,"extr",'i',
		" total N 3 990, both electron excitation and continuum pumping" );

	linadd(rec,990,"rec ",'i',
		" part of N 3 990 due to recombination " );

	linadd(rn3tot*2.01e-11,990,"N 3p",'r',
		" N 3 989.8, continuum pumped" );

	linadd(embesq.em1486+TauLines[ipT1486].Emis->xIntensity,1486,"TOTL",'i',
		" N 4] 1486, total intensity of both lines" );

	PutLine(&TauLines[ipT1486],
		" ");

	PntForLine(1485.,"N  4",&ipnt);
	lindst(embesq.em1486,1485,"N  4",ipnt,'c',true ,
		" the N IV] slow transition by itself " );

	/* >>chng 97 may 02, better expression for dielectronic recombination */
	/*rec = GetLineRec(287,765)*(1.-TauLines[ipT765].ColOvTot);*/
	/* >>chng 02 jul 01, add function to get emission fraction */
	rec = GetLineRec(287,765)*emit_frac(&TauLines[ipT765] );

	/* dielectronic recombination contribution from Nussbaumer and Storey 1984 */
	PutExtra(rec);

	PutLine(&TauLines[ipT765],
		" N 4 765, collisionally excited");

	linadd(MAX2(0.,rec),765,"rec ",'i',
		" N 4 765 recombination," );

	linadd(TauLines[ipT1243].Emis->xIntensity+TauLines[ipT1239].Emis->xIntensity,1240,"TOTL",'i',
		" continuum pumping of NV 1240, N 5 1240, total emission, collisions plus pumping " );
	sum = TauLines[ipT1243].Emis->xIntensity*TauLines[ipT1243].Emis->FracInwd + TauLines[ipT1239].Emis->xIntensity*
	  TauLines[ipT1239].Emis->FracInwd;

	linadd(sum,1240,"Inwd",'i',
		" inward part of N 5 " );
	PutLine(&TauLines[ipT1243],
		" ");
	PutLine(&TauLines[ipT1239],
		" ");

	PutLine(&TauLines[ipT209],
		"  N 5 209, 2s-3p Li seq ");

	PntForLine(6300.,"O  1",&ipnt);
	lindst(CoolHeavy.c6300,6300,"O  1",ipnt,'c',true ,
		"  oxygen total Oxygen I  6300, including continuum optical depth ");

	/* the intensity of [OI] 6300 line due to OH photodistruction */
	/* last term is fraction that emit rather than collisionally deexcited,
	 * factor of 0.55 is branching ratio from chemistry for producing
	 * OI in correct excited state 
	 * this is used to get OH photo formation of [OI] 6300
	 *>>refer	OI	photoexcitation	Storzer, H., & Hollenbach, D. 2000, ApJ, 539, 751-759 
	 * discussion on bottom left side of page 752 */
	/* rate_OH_dissoc is number of OH destruction events, OH -> O + H, cm-3 s-1,
	 * 0.55 is fraction of OH dissociation that lead to pop of upper level of 6300
	 * r12 is energy emitted per unit vol, erg cm-3 s-1, in 6300, due to OH dest */
	rate_OH_dissoc = CO_findrate("PHOTON,OH=>O,H");
	r12 = rate_OH_dissoc * 0.55 * 3.16e-12 * CoolHeavy.c6300_frac_emit;

	lindst( r12*TauLines[ipT6300].Emis->Aul/(TauLines[ipT6300].Emis->Aul+TauLines[ipT6363].Emis->Aul) , 
		6300., "OH p",ipnt , 'i' , false,
		" the intensity of [OI] 6300 line due to OH photodistruction");

	lindst( r12*TauLines[ipT6363].Emis->Aul/(TauLines[ipT6300].Emis->Aul+TauLines[ipT6363].Emis->Aul) , 
		6363., "OH p",ipnt , 'i' , false,
		" the intensity of [OI] 6363 line due to OH photodistruction ");

	PntForLine(6363.,"O  1",&ipnt);
	lindst(CoolHeavy.c6363,6363,"O  1",ipnt,'c',true,
		" total Oxygen I  6363, including continuum optical depth " );

	PntForLine(5577.,"O  1",&ipnt);
	lindst(CoolHeavy.c5577,5577,"O  1",ipnt,'c',true,
		" auroral OI " );

	r13 = rate_OH_dissoc * 0.05 * 3.57e-12 * 0.94*CoolHeavy.c5577_frac_emit;
	lindst( r13 , 5577., "OH p",ipnt , 'i', false,
		" 94% of excitations to highest level decay via 5577" );

	PutLine(&TauLines[ipT63],
		"  O I fine structure line ");

	PutLine(&TauLines[ipT146],
		"  O I fine structure line ");

	linadd(MAX2(0.,CoolHeavy.coolOi),0,"TOIc",'c',
		" total collisional cooling due to 6-level OI atom" );

	linadd(MAX2(0.,-CoolHeavy.coolOi),0,"TOIh",'h',
		" total collisional heating due to 6-level OI atom " );

	/* OI 8446 from six level atom */
	/* >>chng 04 nov 15, upper level for 8446 was incorrect - was 4 should have been 2
	 * bug caught by Yoshiki Matsuoka *
	 *sum = atoms.popoi[4]*TauLines[ipT8446].Pesc*TauLines[ipT8446].Aul*2.36e-12;*/
	sum = atoms.popoi[2]*TauLines[ipT8446].Emis->Pesc*TauLines[ipT8446].Emis->Aul*2.36e-12;
	PntForLine(8446.,"O  1",&ipnt);

	lindst(sum,8446,"6lev",ipnt,'i',false,
		" * \todo	2	 these are actually in the optical depth arrays, should be moved to call PutLine - but some needed vars are not defined in structure " );

	PntForLine(1304.,"O  1",&ipnt);
	sum = atoms.popoi[1]*TauLines[ipT1304].Emis->Pesc*TauLines[ipT1304].Emis->Aul*1.53e-11;
	lindst(sum,1304,"6lev",ipnt,'i',false,
		" OI 1304 from six level atom " );

	PntForLine(1039.,"O  1",&ipnt);
	sum = atoms.popoi[3]*TauLines[ipT1039].Emis->Pesc*TauLines[ipT1039].Emis->Aul*1.92e-11;
	lindst(sum,1039,"6lev",ipnt,'i',false ,
		" OI 1039 from six level atom");

	PntForLine(4368.,"O  1",&ipnt);
	sum = atoms.popoi[5]*TauLines[ipT4368].Emis->Pesc*TauLines[ipT4368].Emis->Aul*4.55e-12;
	lindst(sum,4368,"6lev",ipnt,'i',false,
		" OI 4368 from six level atom" );

	PntForLine(13100.,"O  1",&ipnt);
	sum = atoms.popoi[3]*TauLines[ipTOI13].Emis->Pesc*TauLines[ipTOI13].Emis->Aul*1.52e-12;
	lindst(sum,13100,"6lev",ipnt,'i',false ,
		"OI 1.3 micron from six level atom");

	PntForLine(11300.,"O  1",&ipnt);
	sum = atoms.popoi[4]*TauLines[ipTOI11].Emis->Pesc*TauLines[ipTOI11].Emis->Aul*1.76e-12;
	lindst(sum,11300,"6lev",ipnt,'i',false ,
		" OI 1.1 micron from six level atom");

	PntForLine(29000.,"O  1",&ipnt);
	sum = atoms.popoi[5]*TauLines[ipTOI29].Emis->Pesc*TauLines[ipTOI29].Emis->Aul*6.86e-13;
	lindst(sum,29000,"6lev",ipnt,'i',false ,
		" OI 2.9 micron from six level atom");

	PntForLine(46000.,"O  1",&ipnt);
	sum = atoms.popoi[5]*TauLines[ipTOI46].Emis->Pesc*TauLines[ipTOI46].Emis->Aul*4.32e-13;
	lindst(sum,46000,"6lev",ipnt,'i',false ,
		" OI 4.6 micron from six level atom");

	/*double rec7323 , rec7332, rec3730 , rec3726 , rec2471 
	 * reco23tot , reco22tot;*/

	/* total recombination to 2P^o, the highest two levels of the 5-level atom,
	 * which produces the 7325 multiplet, last factor accounts for coll deexcitation  
	 * this implements equation 2 of 
	 * refer	o2	rec	Liu, X-W., Storey, P.J., Barlow, M.J., Danziger, I.J.,
	 * refercon	Cohen, M., & Bryce, M., 2000, MNRAS, 312, 585 */
	/* >>chng 05 dec 29, from first eqn, or unknown origin, to second, from indicated
	 * reference.  They agreed within 20% */
	/*reco23tot = 3.484e-11 / ( phycon.sqrte / phycon.te05 * phycon.te003 ) * 
		dense.eden * dense.xIonDense[ipOXYGEN][2] * CoolHeavy.O2_A3_tot *2.72e-12;*/
	if( dense.xIonDense[ipHYDROGEN][1]  > SMALLFLOAT )
	{
		/* this test is necessary because for sims with no H-ionizing radiation
		 * the H+ density is initially zero */
		reco23tot = CoolHeavy.O2_A3_tot * HBeta *
			9.36 * phycon.te40*phycon.te04 / 57.544 * dense.xIonDense[ipOXYGEN][2]/dense.xIonDense[ipHYDROGEN][1];
	}
	else
	{
		reco23tot = 0.;
	}

	sum = CoolHeavy.O2471*2471./7325. + CoolHeavy.O7323 + CoolHeavy.O7332;
	if( sum > SMALLFLOAT )
	{
		/* assume effective branching ratio according to predicted intensities from 5-lev atom*/
		reco23tot /= sum;
	}
	else
	{
		reco23tot = 0.;
	}
	/* these are now ergs per sec unit vol for each transition */
	rec7323 = reco23tot * CoolHeavy.O7323;
	rec7332 = reco23tot * CoolHeavy.O7332;
	rec2471 = reco23tot * CoolHeavy.O2471*2471./7325. * 8.05e-12/2.72e-12;

	/* total recombination to 2D^o, the middle two levels of the 5-level atom,
	 * which produces the 3727 multiplet, last factor accounts for coll deexcit */
	reco22tot = 1.660e-10 / ( phycon.sqrte * phycon.te03 * phycon.te005 ) * 
		dense.eden * dense.xIonDense[ipOXYGEN][2] * CoolHeavy.O2_A2_tot;
	/* assume effective branching ratio according to predicted intensities from 5-lev atom*/
	sum = CoolHeavy.O3726 + CoolHeavy.O3730;
	if( sum > SMALLFLOAT )
	{
		reco22tot /= sum;
	}
	else
	{
		reco22tot = 0.;
	}
	/* these are now ergs per sec unit vol for each transition */
	rec3726 = reco22tot * CoolHeavy.O3726 * 5.34e-12;
	rec3730 = reco22tot * CoolHeavy.O3730 * 5.34e-12;

	/* O II 3727 produced by photoionization OF O0 */
	oxy.s3727 = (realnum)((oxy.s3727 + oxy.s7325*0.5)*5.34e-12*
	  9.7e-5/(9.7e-5 + dense.eden*1.15e-6/phycon.sqrte));

	PntForLine(3727.,"O  2",&ipnt);
	fac = CoolHeavy.c3727+oxy.s3727+rec3726+rec3730;
	lindst(fac ,3727,"TOTL",ipnt,'c',true,
		" O II 3727, all lines of multiplet together " );

	PntForLine(7325.,"O  2",&ipnt);
	fac = CoolHeavy.c7325+rec7323+rec7332;
	lindst( fac ,7325,"TOTL",ipnt,'c',true ,
		" O II 7325, all lines of multiplet together");

	linadd(oxy.s3727,3727,"IONZ",'i',
		" line produced by photoionization of Oo; already in TOTL" );
	oxy.s7325 = (realnum)(oxy.s7325*2.72e-12*0.34/(0.34 + dense.eden*
	  6.04e-6/phycon.sqrte));

	linadd(oxy.s7325,7325,"IONZ",'i',
		" line produced by photoionization of Oo; already in TOTL" );

	linadd(CoolHeavy.c7325,7325,"Coll",'i',
		" collisional contribution to line " );

	linadd(CoolHeavy.c3727,3727,"Coll",'i',
	" collisional contribution to line " );

	linadd(CoolHeavy.O3730,3729,"O II",'i',
		" five level atom calculations; D5/2 - S3/2" );

	linadd(CoolHeavy.O3726,3726,"O II",'i',
		" D3/2 - S3/2 transition" );

	linadd(CoolHeavy.O2471,2471,"O II",'c',
		" both 2P 1/2 and 3/2 to ground " );

	linadd(CoolHeavy.O7323,7323,"O II",'i',
		" P1/2-D5/2 and P3/2-D5/2 together" );

	linadd(CoolHeavy.O7332,7332,"O II",'i',
		" P1/2-D3/2 and P3/2-D3/2 together " );

    linadd( rec3730 ,3729,"O 2r",'i',
		" recombination contribution  refer	o2	rec	Liu, X-W., Storey, P.J., Barlow, M.J., Danziger, I.J.,refercon	Cohen, M., & Bryce, M., 2000, MNRAS, 312, 585  recombination contributions five level atom calculations; D5/2 - S3/2 " );

	linadd( rec3726 ,3726,"O 2r",'i',
		" D3/2 - S3/2 transition" );

	linadd(rec2471,2471,"O 2r",'i',
		" both 2P 1/2 and 3/2 to ground " );
	linadd(rec7323,7323,"O 2r",'i',
		" P1/2-D5/2 and P3/2-D5/2 together " );

	linadd(rec7332,7332,"O 2r",'i',
		" P1/2-D3/2 and P3/2-D3/2 together " );

	PutLine(&TauLines[ipT834],
		"  O II 833.8 coll excit ");

	/* the OII multiplets,
	 * contributions by both continuum pumping through XUV line
	 * and recombination */
	/* this is the driving line, pump is photons cm^-3 s^-1 */
	if( nWindLine > 0 )
	{
		pump = TauLine2[387].Emis->pump*TauLine2[387].Emis->PopOpc;
	}
	else
	{
		pump = 0.;
	}

	PntForLine(3120.,"O  2",&ipnt);
	lindst(pump*0.336 * 6.37e-12/(1.+dense.eden/1e12) ,3120,"pump",ipnt,'r',true,
		" OII 3113.62 - 3139.68 (8 lines) are only pumped, no recombination part" );

	PntForLine(3300.,"O  2",&ipnt);
	lindst(pump*0.147 * 6.03e-12/(1.+dense.eden/1e12) ,3300,"pump",ipnt,'r',true,
		" OII 3277.56 - 3306.45 (6 lines) are only pumped, no recombination part" );

	PntForLine(3762.,"O  2",&ipnt);
	lindst(pump*0.087 * 5.29e-12/(1.+dense.eden/1e12) ,3762,"pump",ipnt,'r',true,
		" OII 3739.76/3762.47/3777.42 (3 lines) are only pumped, no recombination part" );

	/* recombination and specific pump for OII 4638.86-4696.35 (8 lines) */
	rec = GetLineRec(82, 4651 );
	PntForLine(4651.,"O  2",&ipnt);
	lindst(rec,4651,"O 2r",ipnt,'r',true,
		" O II 4651 total recombination, 4638.86-4696.35 (8 lines)  " );

	/* convert UV pump rate to intensity with branching ratio and hnu recombination 
	 * part of O II 4651 line */
	linadd(pump* 0.336 * 0.933 * 4.27e-12/(1.+dense.eden/1e12),4651,"O 2p",'i',
		" pumped part of line O II 4651 " );

	/* recombination and specific pump for OII 4317.14-4366.89 (6 lines) */
	rec = GetLineRec(83, 4341 );

	linadd(rec/(1.+dense.eden/1e12),4341,"O 2r",'i',
		" recombination contribution to O II 4341 line " );

	linadd(pump* 0.147 * 0.661 * 4.58e-12/(1.+dense.eden/1e12),4341,"O 2p",'i',
		" pumped part of line O II 4341 " );

	PntForLine(4341.,"O  2",&ipnt);
	lindst(rec+pump* 0.147 * 0.661 * 4.58e-12/(1.+dense.eden/1e12),4341,"TOTL",ipnt,'r',true,
		" total intensity, all processes, O II 4341" );

	/* recombination and specific pump for OII 3712.74/3727.32/3749.48 (3 lines) */
	rec = GetLineRec(84, 3736 );
	/* convert UV pump rate to intensity with branching ratio and hnu */

	linadd(rec/(1.+dense.eden/1e12),3736,"O 2r",'i',"\n recombination part of O II 3736 line " );
	linadd(pump* 0.087 * 0.763 * 5.33e-12/(1.+dense.eden/1e12),3736,"O 2p",'i',
		" pumped part of line O II 3736" );

	PntForLine(3736.,"O  2",&ipnt);
	lindst((rec+pump* 0.087 * 0.763 * 5.33e-12)/(1.+dense.eden/1e12),3736,"TOTL",ipnt,'r',true,
		" total intensity, all processes, O II 3736" );

	/* O III 1661+1666 */
	/*efac = ((1.-TauLines[ipT1666].ColOvTot) + (1.-TauLines[ipT1661].ColOvTot))*0.5;*/
	efac = (emit_frac(&TauLines[ipT1666]) + emit_frac(&TauLines[ipT1661]))*0.5;

	linadd(TauLines[ipT1666].Emis->xIntensity+TauLines[ipT1661].Emis->xIntensity,1665,"TOTL",'i',
		"total intensity of OIII] 1665, all processes " );
	PutLine(&TauLines[ipT1661]," ");

	PutLine(&TauLines[ipT1666]," ");

	linadd(oxy.p1666*1.20e-11*efac,1665,"Phot",'i',
		" contribution to OIII 1665 due to inner shell (2s^2) ionization " );

	linadd(oxy.AugerO3*1.20e-11*efac*0.27,1665,"Augr",'i',
		" contribution to OIII 1665 due to K-shell ionization " );

	PntForLine(5007.,"O  3",&ipnt);
	lindst(CoolHeavy.c5007/(1.+1./3.01),5007,"O  3",ipnt,'c',true ,
		"  O III  5007 alone, collisions, tot OIII is this times 1.333  fac = c5007/(1.+1./2.887) >>chng 01 may 04, branching ratio had been 2.887, revised to 3 as per refer	o3	as	Storey, P.J., & Zeippen, C.J., 2000, 312, 813-816 ");

	PntForLine(4959.,"O  3",&ipnt);
	lindst(CoolHeavy.c5007/(1.+3.01),4959,"O  3",ipnt,'c',true,
		" O III  4959 alone, collisions, tot OIII is this times 4" );

	PntForLine(4931.,"O  3",&ipnt);
	lindst(CoolHeavy.c5007/(1.+3.01)*4.09e-4 ,4931,"O  3",ipnt,'c',true ,
		"  O III  4931 alone, collisions  >>chng 01 jul 11, added this line  >>refer	o3	as	Nussbaumer, H., & Storey, P., 1981, A&A, 99, 177  >>refer	o3	as	Mathis, J.S., & Liu, X.-W., 1999, ApJ, 521, 212-216 ");

	linadd(oxy.d5007t/1.25,5007,"LOST",'i',
		" O III 5007 lost through excited state photo" );

	/* collisional quenching ratio */
	effec = 1.6/(1.6 + 0.9*dense.cdsqte);

	/* O III 4363 recombination, coefficient from Burgess and Seaton */
	r4363 = 6.3e-21/(phycon.te70*phycon.te10)*dense.eden*dense.xIonDense[ipOXYGEN][3]*
	  effec;

	/* charge exchange, 
	 * >>refer	O3	CT	Dalgarno+Sternberg ApJ Let 257, L87.
	 * scaled to agree with 
	 * >>refer	O3	CT	Gargaud et al AA 208, 251, (1989) */
	ct4363 = phycon.sqrte*1.3e-12*4.561e-12*dense.xIonDense[ipHYDROGEN][0]*dense.xIonDense[ipOXYGEN][3]*
	  effec;

	fac = CoolHeavy.c4363 + r4363 + ct4363;
	linadd(fac,4363,"TOTL",'i',
		" O III 4363, sum of rec, coll, ct excitation" );

	PntForLine(4363.,"O  3",&ipnt);
	lindst(CoolHeavy.c4363,4363,"Coll",ipnt,'c',true,
		" O III 4363,collisions from five level atom " );

	lindst(r4363,4363,"Rec ",ipnt,'r',true,
		" O III 4363 recombination, coefficient from Burgess and Seaton " );

	PntForLine(2321.,"O  3",&ipnt);
	lindst(CoolHeavy.c4363*0.236,2321,"O  3",ipnt,'c',true ,
		" collisional excitation of 2321, 5-level atom");
	linadd(ct4363,4363,"C EX",'i' ,
		" call linadd( c4363*0.236 , 2321 , 'O  3','c') charge exchange, Dalgarno+Sternberg ApJ Let 257, L87. ");

	linadd(dense.xIonDense[ipHYDROGEN][0]*dense.xIonDense[ipOXYGEN][3]*0.225*3.56e-12*1.34e-11*phycon.sqrte,
	  5592,"C EX",'i'," charge exchange rate, D+S " );

	PutLine(&TauLines[ipTO88],
		" O III 88 micron, collisionally excited");

	PutLine(&TauLines[ipT52],
		" O III 52 micron, collisionally excited ");

	/* >>chng 97 may 02, better rec contribution */
	/*rec = GetLineRec(331,835)*(1.-TauLines[ipT835].ColOvTot);*/
	rec = GetLineRec(331,835)*emit_frac(&TauLines[ipT835]);
	PutExtra(MAX2(0.,rec));

	PutLine(&TauLines[ipT835],
		"  O III 834A, collisions and dielectronic recombination ");

	linadd(MAX2(0.,rec),835,"rec ",'i',
		" O III 834A, dielectronic recombination only" );

	PutLine(&TauLines[ipT26],
		"  O IV 26 micron ");

	linadd(
		TauLines[ipO4_1400].Emis->xIntensity+
		TauLines[ipO4_1397].Emis->xIntensity+
		TauLines[ipO4_1407].Emis->xIntensity+
		TauLines[ipO4_1405].Emis->xIntensity+
		TauLines[ipO4_1401].Emis->xIntensity,
		1402,"TOTL",'i',
		" total intensity of O IV] 1402, all lines in the multiplet " );

	PutLine(&TauLines[ipO4_1400],
		" ");
	PutLine(&TauLines[ipO4_1397],
		" ");
	PutLine(&TauLines[ipO4_1407],
		" ");
	PutLine(&TauLines[ipO4_1405],
		" ");
	PutLine(&TauLines[ipO4_1401],
		" ");

	linadd(oxy.p1401*1.42e-11,1401,"InSh",'i',
		" inner shell photoionization, relaxation " );

	/* >>chng 97 may 02, better rec contribution */
	rec = GetLineRec(378,789)*(1.-TauLines[ipT789].Emis->ColOvTot);
	rec = GetLineRec(378,789)*emit_frac(&TauLines[ipT789]);
	PutExtra(MAX2(0.,rec));

	PutLine(&TauLines[ipT789]," O IV 789A");

	linadd(MAX2(0.,rec),789,"rec ",'i',
		" O IV 789A, dielectronic recombination only" );

	/* >>chng 97 may 02, better rec contribution */
	rec = GetLineRec(466,630);
	PutExtra(MAX2(0.,rec));

	PutLine(&TauLines[ipT630],"O V 630, collisional excitation and dielectronic recombination");

	linadd(MAX2(0.,rec),630,"rec ",'i',
		" O V 630A, dielectronic recombination only" );

	linadd(embesq.em1218+TauLines[ipT1214].Emis->xIntensity,1218,"TOTL",'i',
		" O V 1218], total intensity of both lines " );
	PutLine(&TauLines[ipT1214], " ");

	linadd(embesq.em1218,1211,"O  5",'i',
		" the slow transition by itself" );

	linadd(1.4e-21/phycon.te70*dense.eden*dense.xIonDense[ipOXYGEN][5]*
	  /*(1.-TauLines[ipT1214].ColOvTot),5112,"O  5",'i' );*/
	  emit_frac(&TauLines[ipT1214]),5112,"O  5",'i',
	  " BS O V 5112, recombination " );

	linadd(TauLines[ipT1032].Emis->xIntensity+TauLines[ipT1037].Emis->xIntensity,1035,"TOTL",'i',
		" O VI 1035, total of pumping and collisional excitation " );
	sum = TauLines[ipT1032].Emis->xIntensity*TauLines[ipT1032].Emis->FracInwd + 
		TauLines[ipT1037].Emis->xIntensity* TauLines[ipT1037].Emis->FracInwd;

	linadd(sum,1035,"Inwd",'i',
		"  inward part of OVI line" );
	PutLine(&TauLines[ipT1032],
		" ");
	PutLine(&TauLines[ipT1037],
		" ");

	PutLine(&TauLines[ipT150],
		"O VI 150, Li seq 2s 3p ");

	PutLine(&TauLines[ipF0229],
		" fluorine	 * [F II] 29.33 microns");

	PutLine(&TauLines[ipF0267],
		"  [F II] 67.2 microns ");

	PutLine(&TauLines[ipF444],
		"  F IV 44.07 micron ");

	PutLine(&TauLines[ipF425],
		"  F IV 25.83 micron");

	PutLine(&TauLines[ipTNe13],
		"  neon Neon II 12.8 micron ");

	PutLine(&TauLines[ipTNe16],
		"  Ne III fine structure line ");

	PutLine(&TauLines[ipTNe36],
		"  Ne III fine structure line ");

	PntForLine(3869.,"Ne 3",&ipnt);
	lindst(CoolHeavy.c3869/(1.+1./3.318),3869,"Ne 3",ipnt,'c',true,
		"  Ne III  3869, of 3968+3869 doublet" );

	PntForLine(3968.,"Ne 3",&ipnt);
	lindst(CoolHeavy.c3869/(1.+3.318),3968,"Ne 3",ipnt,'c',true,
		" Ne III  3968, of 3968+3869 doublet" );

	PntForLine(3343.,"Ne 3",&ipnt);
	lindst(CoolHeavy.c3343,3343,"Ne 3",ipnt,'c',true,
		" NeIII auroral line " );

	PntForLine(1815.,"Ne 3",&ipnt);
	lindst(CoolHeavy.c3343*1.38,1815,"Ne 3",ipnt,'c',true ,
		"  NeIII auroral line");

	PntForLine(2424.,"Ne 4",&ipnt);
	lindst(CoolHeavy.c2424,2424,"Ne 4",ipnt,'c',true,
		" Ne IV 2424, collisional excitation" );

	PntForLine(4720.,"Ne 4",&ipnt);
	lindst(CoolHeavy.c4720,4720,"Ne 4",ipnt,'c',true,
		"  Ne IV N=3-2 lines, three level atom approx, this is the sum of the 4714.5, 4724.2, 4725.5 lines" );

	PntForLine(1602.,"Ne 4",&ipnt);
	lindst(CoolHeavy.c4720*4.34,1602,"Ne 4",ipnt,'c',true ,
		"  Ne IV N=3 lines, three level atom approx");

	PntForLine(3426.,"Ne 5",&ipnt);
	lindst(CoolHeavy.c3426/(1.+1./2.738),3426,"Ne 5",ipnt,'c',true,
		" Ne V 3426 of 3426, 3346 doublet" );

	PntForLine(3346.,"Ne 5",&ipnt);
	lindst(CoolHeavy.c3426/(1.+2.738),3346,"Ne 5",ipnt,'c',true,
		" Ne V 3346 of 3426, 3346 doublet " );

	PntForLine(2976.,"Ne 5",&ipnt);
	lindst(CoolHeavy.c2975,2976,"Ne 5",ipnt,'c',true,
		" auroral line " );

	PntForLine(1575.,"Ne 5",&ipnt);
	lindst(CoolHeavy.c1565,1575,"Ne 5",ipnt,'c',true,
		   " collisionally excited" );

	PutLine(&TauLines[ipTNe24],"\n  Ne V 24.2, 14.3 micron ");

	PutLine(&TauLines[ipTNe14],"\n  Ne V 24.2, 14.3 micron ");

	PntForLine(1141.,"Ne 5",&ipnt);
	lindst(CoolHeavy.c1134,1141,"Ne 5",ipnt,'c',true," both components of 5S-3P 1146.1, 1137.0 doublet " );

	PutLine(&TauLines[ipxNe0676],"\n  [Ne VI] 7.6 microns ");

	linadd(embesq.em895+TauLines[ipT895].Emis->xIntensity,895,"TOTL",'i',
		" Ne VII 895, collisionally excited, both lines " );

	PutLine(&TauLines[ipT895],
		"  Ne VII 895, only fast transition ");

	linadd(embesq.em895,890,"Ne 7",'i',
		" Ne VII 890, single line " );

	linadd(TauLines[ipT770].Emis->xIntensity+TauLines[ipT780].Emis->xIntensity,774,"TOTL",'i',
		" Ne VIII 774, collisionally excited " );

	sum = TauLines[ipT770].Emis->xIntensity*TauLines[ipT770].Emis->FracInwd + 
		TauLines[ipT780].Emis->xIntensity*TauLines[ipT780].Emis->FracInwd;
	linadd(sum,774,"Inwd",'i',
		" inward part of NeVIII 774 line" );

	PutLine(&TauLines[ipT770],
		" the NeVIII 770 780 doublet ");
	PutLine(&TauLines[ipT780],
		" ");

	PutLine(&TauLines[ipT88],
		"  Ne VIII 88 2s 3p, collisionally excited ");

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "   lines_lv1_li_ne returns\n" );
	}
	return;
}

/*GetLineRec return rec coef*hnu*eden*n_ion for C, N, or O recombination lines from Dima's list,
 * also zero's line in master stack so not entered second time in later dump of all rec lines */
STATIC double GetLineRec(
	/* this is the number of the emission line in the stack of lines, on the C scale */
	long int ip, 
	/* the multiplet wavelength */
  long int lWl)
{
	double GetLineRec_v;

	DEBUG_ENTRY( "GetLineRec()" );

	if( (long)(LineSave.RecCoefCNO[2][ip]+0.5) != lWl )
	{
		fprintf( ioQQQ, " GetLineRec called with incorrect wavelength.\n" );
		fprintf( ioQQQ, " index, call and get wl are %5ld%5ld%5ld\n", 
		  ip, lWl, (long)(LineSave.RecCoefCNO[2][ip]+0.5) );
		cdEXIT(EXIT_FAILURE);
	}

	/* final product is vol emissivity in line */
	GetLineRec_v = LineSave.RecCoefCNO[3][ip]*dense.eden*
		dense.xIonDense[(long)(LineSave.RecCoefCNO[0][ip])-1][(long)(LineSave.RecCoefCNO[0][ip]-LineSave.RecCoefCNO[1][ip]+2)-1]*
	  1.99e-8/LineSave.RecCoefCNO[2][ip];

	/* zero out rec coefficient so that not used again in master dump
	 * this routine cannot be called twice on same line */
	LineSave.RecCoefCNO[3][ip] = 0.;
	return( GetLineRec_v );
}
