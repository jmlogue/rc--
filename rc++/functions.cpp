#include "functions.h"
#include <iomanip> // RAD: so far used only for setprecission() in cmd output

using namespace std;

double g = 9.81;
double pi = 3.141592653;
const int ArraySize = 16;

string strUppercase(string stringvar);
int sgn(double sgnvar);

void cptheta(double CP[4][4], double& windangle, double* Cpwall);

void Flueflow(double& Tin, double& Swrfl, double& dPwind, double& dPtemp, double& P,
	double& h, double& Pint, double& nflue, int& Nrflue, Chimtype* Chim, double& Mflue,
	double& rhoo, double& rhoi, double& dPflue, double& Tout, double& aeq);

void Floorflow3(double& Cfloor, double& Cpfloor, double& dPwind, double& Pint, double& C,
	double& n, double& Mfloor, double& rhoo, double& rhoi, double& dPfloor, double& Hfloor, double& dPtemp);

void Ceilflow(int& ahflag, double& R, double& X, double& Patticint, double& h, double& dPtemp,
	double& dPwind, double& Pint, double& C, double& n, double& Mceil, double& Cattic, double& rhoa,
	double& rhoi, double& dPceil, double& Tattic, double& Tin, double& Tout, double& rhoo,
	double& msupahoff, double& mretahoff, double& supC, double& supn, double& retC, double& retn);

void Neutrallevel2(double& dPtemp, double& dPwind, double* Sw, double& Pint, double* Cpwall, double* Bo, double& h);

void Wallflow3(double& Tin, double& Tout, double& rhoi, double& rhoo, double& Bo, double& Cpwall,
	double& n, double& Cwall, double& h, double& Pint, double& dPtemp, double& dPwind, double& Mwall,
	double& Mwallin, double& Mwallout, double& dPwalltop, double& dPwallbottom, double& Hfloor);

void Fanflow(fantype& fan, double& rhoo, double& rhoi);

void Pipeflow(double& rhoo, double& rhoi, double& CP, double& dPwind, double& dPtemp,
	double& Pint, pipetype& Pipe, double& Tin, double& Tout);

void Windoorflow(double& Tin, double& Tout, double& rhoi, double& rhoo, double& h, double& Bo,
	double& Cpwall, double& n, double& Pint, double& dPtemp, double& dPwind, windoortype& Windoor);

void RoofCptheta(double* Cproof, double& windangle, double* Cppitch, double& pitch);

void Neutrallevel3(double& dPtemp, double& dPwind, double& Patticint, double& Cpr, double& Broofo, double& Hpeak);

void roofflow(double& Tattic, double& Tout, double& rhoa, double& rhoo, double& Broofo, double& Cpr,
	double& Nattic, double& Croof, double& Hpeak, double& Patticint, double& dPtemp, double& dPwind,
	double& Mroof, double& Mroofin, double& Mroofout, double& dProoftop, double& dProofbottom, double& H);

void Aventflow(double& rhoo, double& rhoa, double& CP, double& dPwind, double& dPtemp, double& Patticint,
	AventType& Avent, double& Tattic, double& Tout);

void Soffitflow(double& rhoo, double& rhoa, double& CP, double& dPwind, double& dPtemp, double& Patticint,
	SoffitType& Soffit, double& Soffitfraction, double& Cattic, double& Nattic, double& Tattic, double& Tout);

void AFanflow(fantype& AFan, double& rhoo, double& rhoa);

// ----- MatSEqn forward declarations -----
/*
===================================================================
In theory, a matrix can be exactly singular.  Numerically, exact singularity is a rare occurrence because 
round off error turns exact zeroes into very small numbers. This data corruption makes it necessary to set
a criterion to determine how small a number has to be before we flag it as zero and call the matrix singular.

The following constants set the maximum drop allowed in weighted pivot values without error code -1 (matrix singular) 
being returned.  To increase the singularity sensitivity of the MatSEqn routine, increase the size of feps for 
floating precision calls or deps for double precision calls.  Similarly, decreasing the size of the constants will 
make the routines less sensitive to singularity.
===================================================================
*/

const float feps = .00001f;
const double deps = .00000000001;

int MatSEqn(double A[][ArraySize], double* b, int asize, int asize2, int bsize);
int matlu(double A[][ArraySize], int* rpvt, int* cpvt, int asize, int asize2, int& continuevar);
int matbs(double A[][ArraySize], double* b, double* x, int* rpvt, int* cpvt, int asize);

//**********************************
// Functions definitions...

void heat ( 
	double& Tout, 
	double& rhoref, 
	double& Mceil, 
	double& AL4, 
	double& speed, 
	double& ssolrad, 
	double& nsolrad, 
	double* told, 
	double& ATTICVOL, 
	double& housevol, 
	double& sc, 
	double* b, 
	int& ERRCODE, 
	double& TSKY, 
	double& FLOORAREA, 
	double& pitch, 
	double& location, 
	double& MSUPREG, 
	double& mretreg, 
	double& mretleak, 
	double& MSUPLEAK, 
	double& MAH, 
	double& supRval, 
	double& retRval, 
	double& supdiam, 
	double& retdiam, 
	double& suparea, 
	double& retarea, 
	double& supthick, 
	double& retthick, 
	double& supvol, 
	double& retvol, 
	double& supcp, 
	double& retcp, 
	double& SUPVEL, 
	double& retvel, 
	double& suprho, 
	double& retrho, 
	double& pref, 
	double& HROUT, 
	double& diffuse, 
	double& UA, 
	double& matticenvin, 
	double& matticenvout, 
	double& mhousein, 
	double& mhouseout, 
	double& planarea, 
	double& msupahoff, 
	double& mretahoff, 
	double& solgain, 
	double& windowS, 
	double& windowN, 
	double& windowWE, 
	double& shadcoef, 
	double& mfancyc, 
	double& Hpeak, 
	double& h, 
	double& Whouse,
	double& retlength,
	double& suplength,
	int& rooftype,
	double& M1,
	double& M12,
	double& M15,
	double& M16,
	int& rroof,
	double& rceil,
	int& ahflag, 
	double& dtau,
	double& merv,
	double& mhrv,
	double& SBETA,
	double& CBETA,
	double& L,
	double& dec,
	double& Csol,
	int& idirect,
	double& capacityc,
	double& capacityh,
	double& evapcap,
	double& intgain,
	int bsize
) {

	int RHOSHEATHING;
	int RHOWOOD;
	int RHOSHINGLES;
	int rhodrywall;
	int heatcount;
	int asize;
	int asize2;
	
	double incsolar[4] = {0,0,0,0};
	double A[16][ArraySize];
	double toldcur[16];
	double RHOATTIC, RHOhouse, RHOSUP, RHORET;
	double rhoout;
	double WOODTHICK;
	double SIGMA;
	double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13;
	double pws;
	double PW;
	double A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16;
	double ahouse;
	double mshingles;
	double M2, M3, M4, M5, M6, M7, M8, M9, M10, M11, M13, M14;
	double CPAIR;
	double cp1, cp2, cp3, cp4, cp5, cp6, cp7, cp8, cp9, cp10, cp11, cp12, cp13, cp14, cp15, cp16;
	double KWOOD;
	double kair;
	double muair;
	double Rshingles;
	double Rval2, Rval3, Rval4, Rval5, Rval7, Rval8, Rval9, Rval10, Rval11, Rval14;
	double u;
	double Hnat2, Hnat3, Hnat4, Hnat5, Hnat6, Hnat8, Hnat9, Hnat11, Hnat14;
	double H2, H3, H4, H5, H6, H7, H8, H9, H10, H11, H13, H14;
	double tfilm2, tfilm3, tfilm4, tfilm5, tfilm6, tfilm8, tfilm9, tfilm10, tfilm11, tfilm14;
	double Hforced2, Hforced3, Hforced4, Hforced5, Hforced6, Hforced8, Hforced9, Hforced11, Hforced14;
	double HI11, HI14;
	double F8t2, F8t4, F8t11, F8t14;
	double F2t4, F2t8, F2t11, F2t14;
	double F4t2, F4t8, F4t11, F4t14;
	double F11t2, F11t4;
	double F14t2, F14t4;
	double HR2t4, HR2t8, HR2t11, HR2t14;
	double HR4t2, HR4t8, HR4t11, HR4t14;
	double HR8t4, HR8t2;
	double HR11t2, HR11t4;
	double HR14t2, HR14t4;
	double R2t4, R2t8, R2t11, R2t14;
	double R4t2, R4t8, R4t11, R4t14;
	double R8t2, R8t4;
	double R11t2, R11t4;
	double R14t2, R14t4;
	double EPS1, epsshingles;	
	double hr7;
	double HRG3, HRG5, HRS3, HRS5;
	double FRS, FG3, FG5;
	double R7, RG3, RG5, RS3, RS5;
	double Beta;	
	double TGROUND;
	double alpha3, alpha5;
	double phi, sphi, cphi, cphi2;
	double incsolarS, incsolarW, incsolarN, incsolarE, incsolarvar;
	double S;
	double Gamma;
	double ct;
	double tsolair;


	// Node Identification
	// Node 1 is the Attic Air
	// Node 2 is the Inner North Sheathing
	// Node 3 is the Outer North Sheathing
	// Node 4 is the Inner South Sheathing
	// Node 5 is the Outer South Sheathing
	// Node 6 is all of the Wood (joists, trusses, etc.) lumped together
	// Node 7 is the Ceiling of the House
	// Node 8 is the Floor of the Attic
	// Node 9 is the Inner Gable Wall (both lumped together)
	// Node 10 is the Outer Gable Wall (both lumped together)
	// Node 11 is the Return Duct Outer Surface
	// Node 12 is the Return Duct Air
	// Node 13 is The Mass of the House
	// Node 14 is the Supply Duct Outer Surface
	// Node 15 is the Supply Duct Air
	// Node 16 is the House Air (all one zone)

	// Densities (Reference density is in the input file
	RHOATTIC = rhoref * 293 / told[0];
	RHOhouse = rhoref * 293 / told[15];

	RHOSUP = rhoref * 293 / told[14];
	RHORET = rhoref * 293 / told[11];
	rhoout = rhoref * 293 / Tout;

	for(int i=0; i < 16; i++) {
		for(int j=0; j < 16; j++) {
			A[i][j] = 0;
		}
	}

	WOODTHICK = .015;			// thickness of sheathing material
	SIGMA = 5.669E-08;			// STEPHAN-BOLTZMANN CONST

	// pi = 3.141592653# 		// FF: pi is already decleared as a global variable

	// THE FOLLOWING ARE CONSTANTS TO DETERMINE WATER VAPOUR PRESSURE
	c1 = -5674.5359;
	c2 = 6.3925274;
	c3 = -9.677843E-03;
	c4 = .00000062215701;
	c5 = .000000002074782;
	c6 = 9.484024E-13;
	c7 = 4.1635019;
	c8 = -5800.2206;
	c9 = 1.3914993;
	c10 = -.04860239;
	c11 = .000041764768;
	c12 = -.000000014452093;
	c13 = 6.5459673;

	// THE FOLLOWING IS FROM ASHRAE 1989
	if(Tout > 273) {
		pws = exp(c8 / Tout + c9 + c10 * Tout + c11 * Tout * Tout + c12 * pow(Tout, 3) + c13 * log(Tout));
	} else {
		pws = exp(c1 / Tout + c2 + c3 * Tout + c4 * Tout * Tout + c5 * pow(Tout, 3) + c6 * pow(Tout, 4) + c7 * log(Tout));
	}

	PW = HROUT * pref / (.62198 + HROUT);						// water vapor partial pressure pg 1.9 ASHRAE fundamentals 2009
	PW = PW / 1000 / 3.38;										// CONVERT TO INCHES OF HG
	TSKY = Tout * pow((.55 + .33 * sqrt(PW)), .25);				// TSKY DEPENDS ON PW

	// Surface Area of Nodes
	A2 = planarea / 2 / cos(pitch * pi / 180);					// PITCHED SLOPE AREA
	A3 = A2;													// ALL SHEATHING SURFACES HAVE THE SAME AREA
	A4 = A2;
	A5 = A2;

	// the following are commented out for ConSOl becasue cement tile is flat and does not have increased surface area
	// IF rooftype = 2 OR rooftype = 3 THEN
	//        // tile roof has more surface area for convection heat transfer
	//        A3 = 1.5 * A2
	//        A5 = A3
	// END IF

	A6 = planarea * 1.5;										// ATTIC WOOD SURF AREA
	A7 = planarea;
	A8 = A7;
	A9 = planarea / 2 * tan(pitch * pi / 180);					// total endwall area
	A10 = A9;
	A12 = pi * retlength * retdiam;
	A15 = pi * suplength * supdiam;
	A11 = retarea;
	A14 = suparea;

	// assumes house is square with one 2.5 m and one 3m ceiling
	ahouse = 3 * pow(FLOORAREA, .5) * 2 + 2.5 * pow(FLOORAREA, .5) * 2 + 2 * FLOORAREA;

	A16 = ahouse;
	A13 = 6 * A16;												//  maybe need to really make A13 much bigger  its the surface area of everything in the house
	//  not just the floor  Something like five times the floor area might be appropriate

	RHOSHEATHING = 450;
	RHOWOOD = 500;												// WOOD DENSITY
	RHOSHINGLES = 1100 * 2;										// factor of two because they overlap
	rhodrywall = 800;

	// masses
	mshingles = RHOSHINGLES * .005 * A2;

	if(rooftype == 2 || rooftype == 3) {
		mshingles = 50 * A2;
	}

	M1 = ATTICVOL * RHOATTIC;											// mass of attic air
	M2 = .5 * A2 * RHOSHEATHING * WOODTHICK;							// 1/2 OF TOTAL

	// OTHER 1/2 OUTSIDE SHEATHING
	// WOOD TCOND W/MMC
	M3 = M2 + mshingles;
	M4 = .5 * A4 * RHOSHEATHING * WOODTHICK;
	M5 = M4 + mshingles;
	M6 = 10 * planarea;													// Wild speculation
	M7 = .5 * 4 * planarea;												// MASS OF JOISTS DRYWALL AND INSULATION
	M8 = M7;
	M9 = .5 * RHOWOOD * WOODTHICK * A9;
	M10 = M9;
	M11 = retlength * pi * (retdiam + retthick) * retthick * retrho;	// retarea * retrho * retthick
	M12 = retvol * RHORET;

	//  maybe not - 02/2004 need to increase house mass with furnishings and their area: say 5000kg furnishings
	//  and A13 increased to something like 2.5 times floor area
	M13 = (2.5 * pow(FLOORAREA, .5) * 4 * 2000 * .01 + planarea * .05 * 2000);		// mass of walls (5 cm effctive thickness)+mass of slab (aso 5 cm thick)
	M14 = suplength * pi * (supdiam + supthick) * supthick * suprho;				// suparea * suprho * supthick
	M15 = supvol * RHOSUP;
	M16 = housevol * RHOhouse;

	// specific heats
	cp1 = 1005.7;														// specific heat of air [j/kg C]
	CPAIR = cp1;
	cp2 = 1210;															// CP plywood
	cp3 = 1260;															// CP asphalt shingles
	cp4 = cp2;
	if(rooftype == 2 || rooftype == 3) {
		cp3 = 880;														// CP for tiles roof
		cp5 = cp3;
	}
	cp5 = cp4;
	cp6 = 1630;															// CP wood
	cp7 = 1150;
	cp8 = cp7;
	cp9 = cp2;
	cp10 = cp9;
	cp11 = retcp;														// input
	cp12 = cp1;
	cp13 = 1300;														// combination of wood and drywall
	cp14 = supcp;
	cp15 = cp1;
	cp16 = cp1;

	// conductivities and R-values and the like
	KWOOD = .15;														// check with Iain about this
	kair = .02624;														// MAKE FUNCITON OF TEMP
	muair = .000018462;													// SAME  (these values at 300K)

	if(rooftype == 1) {													// asphalt shingles
		Rshingles = .077;												// ashrae fundamentals 97 pg 24.5
	} else if(rooftype == 2) {											// red clay tile
		Rshingles = .5;
	} else if(rooftype == 3) {											// low coating clay tile
		Rshingles = .5;
	} else if(rooftype == 4) {											// asphalt shingle & white coating
		Rshingles = .077;
	}

	// changed to account for cathedralized attics
	if(rroof == 0) {
		Rval2 = (WOODTHICK / KWOOD) + Rshingles;
	} else {
		Rval2 = rroof + Rshingles;
	}

	Rval3 = Rval2;
	Rval4 = Rval2;
	Rval5 = Rval2;
	Rval7 = rceil;														// EFFECTIVE THERMAL RESISTANCE OF CEILING
	Rval8 = Rval7;

	// Rval9 = 2.3														// rvalue of insulated gable end walls
	Rval9 = .5;															// rvalue of uninsulated gable end walls

	Rval10 = Rval9;
	Rval11 = retRval;
	Rval14 = supRval;

	// most of the surfaces in the attic undergoe both natural and forced convection
	// the overall convection is determined the forced and natural convection coefficients
	// to the THIRD power, adding them, and cube rooting the result.  This is Iain// s idea
	// and it seemed to work for him

	// Characterstic velocity
	u = (matticenvin - matticenvout) / RHOATTIC / AL4 / 4.0;
	if(u == 0)
		u = abs(Mceil) / 2 / RHOATTIC / AL4 * 2 / 4.0;
	if(u == 0)
		u = .1;

	// ITERATION OF TEMPERATURES WITHIN HEAT SUBROUTINE
	// THIS ITERATES BETWEEN ALL TEMPERATURES BEFORE RETURNING TO MAIN PROGRAM
	heatcount = 0;

	while(1) {

		// FF: sets to 0 all elements inside array b and A
		for(int i=0; i < 16; i++) {
			for(int j=0; j < 16; j++) {
				A[i][j] = 0;
			}
			b[i] = 0;
		}

		heatcount = heatcount + 1;

		if(heatcount == 1) {
			for(int i=0; i < 16; i++) {
				toldcur[i] = told[i];
			}
		}

		// inner north sheathing
		Hnat2 = 3.2 * pow(abs(told[1] - told[0]), (1 / 3.0));
		tfilm2 = (told[1] + told[0]) / 2;
		Hforced2 = (18.192 - .0378 * tfilm2) * pow(u, .8);
		H2 = pow((pow(Hnat2, 3) + pow(Hforced2, 3)), .333333);

		// outer north sheathing
		Hnat3 = 3.2 * pow(abs(told[2] - Tout), (1 / 3.0));
		tfilm3 = (told[2] + Tout) / 2;
		Hforced3 = (18.192 - .0378 * tfilm3) * pow(speed, .8);
		H3 = pow((pow(Hnat3, 3) + pow(Hforced3, 3)), .333333);   		// ATTIC INTERNAL CONV COEF

		// inner south sheathing
		Hnat4 = 3.2 * pow(abs(told[3] - told[0]), (1 / 3.0)); 			// Natural convection from Ford
		tfilm4 = (told[3] + told[0]) / 2;                				// film temperature
		Hforced4 = (18.192 - .037 * tfilm4) * pow(u, .8);    			// force convection from Ford
		H4 = pow((pow(Hnat4, 3) + pow(Hforced4, 3)), .333333);       	// ATTIC INTERNAL CONV COEF

		// outer north sheathing
		Hnat5 = 3.2 * pow(abs(told[4] - Tout), (1 / 3.0));
		tfilm5 = (told[4] + Tout) / 2;
		Hforced5 = (18.192 - .0378 * tfilm5) * pow(speed, .8);
		H5 = pow((pow(Hnat5, 3) + pow(Hforced5, 3)), .333333);   		// ATTIC INTERNAL CONV COEF

		// Wood (joists,truss,etc.)
		Hnat6 = 3.2 * pow(abs(told[5] - told[0]), (1 / 3.0));
		tfilm6 = (told[5] + told[0]) / 2;
		Hforced6 = (18.192 - .0378 * (tfilm6)) * pow(u, .8);
		H6 = pow((pow(Hnat6, 3) + pow(Hforced6, 3)), .333333);

		// Underside of Ceiling
		// modified to use fixed numbers from ASHRAE fund ch.3
		//  on 05/18/2000
		H7 = 6;
		if(ahflag != 0) {
			H7 = 9;
		}

		// House Mass
		// uses ceiling heat transfer coefficient as rest fo house heat transfer coefficient	
		H13 = H7;

		// Attic Floor
		Hnat8 = 3.2 * pow(abs(told[7] - told[0]), (1 / 3.0));
		tfilm8 = (told[7] + told[0]) / 2;
		Hforced8 = (18.192 - .0378 * (tfilm8)) * pow(u, .8);
		H8 = pow((pow(Hnat8, 3) + pow(Hforced8, 3)), .333333);

		// Inner side of gable endwalls (lumped together)
		Hnat9 = 3.2 * pow(abs(told[8] - told[0]), (1 / 3.0));
		tfilm9 = (told[8] + told[0]) / 2;
		Hforced9 = (18.192 - .037 * (tfilm9)) * pow(u, .8);
		H9 = pow((pow(Hnat9, 3) + pow(Hforced9, 3)), .333333);

		// Outer side of gable ends
		tfilm10 = (told[9] + Tout) / 2;
		H10 = (18.192 - .0378 * (tfilm10)) * pow(speed, .8);

		// Outer Surface of Return Ducts
		Hnat11 = 3.2 * pow(abs(told[10] - told[0]), (1 / 3.0));
		tfilm11 = (told[10] + told[0]) / 2;
		Hforced11 = (18.192 - .0378 * (tfilm11)) * pow(u, .8);
		H11 = pow((pow(Hnat11, 3) + pow(Hforced11, 3)), .333333);

		// Inner Surface of Return Ducts
		// from Holman   Nu(D) = 0.023*Re(D)^0.8*Pr(D)^0.4
		// Note Use of HI notation
		HI11 = .023 * kair / retdiam * pow((retdiam * RHORET * abs(retvel) / muair), .8) * pow((CPAIR * muair / kair), .4);

		if(HI11 <= 0)
			HI11 = H11;

		// Outer Surface of Supply Ducts
		Hnat14 = 3.2 * pow(abs(told[13] - told[0]), (1 / 3.0));
		tfilm14 = (told[13] + told[0]) / 2;
		Hforced14 = (18.192 - .0378 * (tfilm14)) * pow(u, .8);
		H14 = pow((pow(Hnat14, 3) + pow(Hforced14, 3)), .333333);

		// Inner Surface of Supply Ducts
		// from Holman   Nu(D) = 0.023*Re(D)^0.8*Pr(D)^0.4
		// Note Use of HI notation
		HI14 = .023 * kair / supdiam * pow((supdiam * RHOSUP * SUPVEL / muair), .8) * pow((CPAIR * muair / kair), .4);
		// I think that the above may be an empirical relationship

		if(HI14 <= 0)
			HI14 = H14;
		// Radiation shape factors


		if(location == 1) { //  Ducts in the house
			// convection heat transfer coefficients

			// Outer Surface of Return Ducts
			H11 = 6;
			if(ahflag != 0)
				H11 = 9;

			// Inner Surface of Return Ducts
			// from Holman   Nu(D) = 0.023*Re(D)^0.8*Pr(D)^0.4
			// Note Use of HI notation
			HI11 = .023 * kair / retdiam * pow((retdiam * RHORET * abs(retvel) / muair), .8) * pow((CPAIR * muair / kair), .4);
			// I think that the above may be an imperical relationship
			if(HI11 <= 0)
				HI11 = H11;

			// Outer Surface of Supply Ducts
			H14 = 6;

			if(ahflag != 0)
				H14 = 9;

			// Inner Surface of Supply Ducts
			// from Holman   Nu(D) = 0.023*Re(D)^0.8*Pr(D)^0.4
			// Note Use of HI notation
			HI14 = .023 * kair / supdiam * pow((supdiam * RHOSUP * SUPVEL / muair), .8) * pow((CPAIR * muair / kair), .4);
			// I think that the above may be an empirical relationship

			if(HI14 <= 0)
				HI14 = H14;
			// Radiation shape factors

			// radiation heat transfer coefficients

			// Only 5 nodes (2,4,8,11,14) are involved in radiation transfer in the attic
			// The endwalls have a very small contribution to radiation exchange and are neglected
			// The wood may or may not contribute to radiation exchange, but their geometry is
			// too complex to make any assumptions.  So it is excluded
			// Assumes that the duct is suspended above the floor, completely out of the insulation
			// this will change in the future

			F8t2 = 1 / 2.0;
			F8t4 = F8t2;
			F2t8 = F8t2 * A8 / A2;
			F4t8 = F2t8;
			F2t4 = (1 - F2t8);
			F4t2 = (1 - F4t8);
			F4t11 = 0;
			F4t14 = 0;
			F2t11 = 0;
			F2t14 = 0;
			F14t2 = 0;
			F14t4 = 0;
			F11t4 = 0;
			F11t2 = 0;
			// Radiation Heat Transfer Coefficiencts
			EPS1 = .9;         										// Emissivity of building materials
			epsshingles = .9;  										// this should be a user input

			// make a loop to do this
			// North Sheathing
			R2t4 = (1 - EPS1) / EPS1 + 1 / F2t4 + (1 - EPS1) / EPS1 * (A2 / A4);
			R2t8 = (1 - EPS1) / EPS1 + 1 / F2t8 + (1 - EPS1) / EPS1 * (A2 / A8);
			HR2t4 = SIGMA * (told[1] + told[3]) * (pow(told[1], 2) + pow(told[3], 2)) / R2t4;
			HR2t8 = SIGMA * (told[1] + told[7]) * (pow(told[1], 2) + pow(told[7], 2)) / R2t8;

			// South Sheathing
			R4t2 = (1 - EPS1) / EPS1 + 1 / F4t2 + (1 - EPS1) / EPS1 * (A4 / A2);
			R4t8 = (1 - EPS1) / EPS1 + 1 / F4t8 + (1 - EPS1) / EPS1 * (A4 / A8);
			HR4t2 = SIGMA * (told[3] + told[1]) * (pow(told[3], 2) + pow(told[1], 2)) / R4t2;
			HR4t8 = SIGMA * (told[3] + told[7]) * (pow(told[3], 2) + pow(told[7], 2)) / R4t8;

			// Attic Floor
			R8t4 = (1 - EPS1) / EPS1 + 1 / F8t4 + (1 - EPS1) / EPS1 * (A8 / A4);
			R8t2 = (1 - EPS1) / EPS1 + 1 / F8t2 + (1 - EPS1) / EPS1 * (A8 / A2);
			HR8t4 = SIGMA * (told[7] + told[3]) * (pow(told[7], 2) + pow(told[3], 2)) / R8t4;
			HR8t2 = SIGMA * (told[7] + told[1]) * (pow(told[7], 2) + pow(told[1], 2)) / R8t2;

		} else {
			// ducts in the attic

			// 33.3% of each duct sees each sheathing surface (top third of duct)
			F14t2 = 1 / 2.0;
			F11t2 = F14t2;
			F14t4 = F14t2;
			F11t4 = F14t2;


			// Remaining 50% of each duct surface sees the floor
			// changed, the ducts don't see the floor
			// F11t8 = 0
			// F14t8 = F11t8

			// The ducts don't see each other
			// F11t14 = 0
			// F14t11 = 0

			F8t14 = 0;											// F14t8 * (A14 / 3) / A8
			F8t11 = 0;											// F11t8 * (A11 / 3) / A8

			F2t14 = F14t2 * (A14 / 3) / A2;
			F2t11 = F11t2 * (A11 / 3) / A2;

			F4t14 = F14t4 * (A14 / 3) / A4;
			F4t11 = F11t4 * (A11 / 3) / A4;

			F8t2 = 1 / 2.0;										// (1 - F8t14 - F8t11) / 2
			F8t4 = F8t2;

			F2t8 = F8t2 * A8 / A2;
			F4t8 = F2t8;
			F2t4 = (1 - F2t8 - F2t11 - F2t14);
			F4t2 = (1 - F4t8 - F4t11 - F4t14);


			// Radiation Heat Transfer Coefficiencts
			EPS1 = .9;         									// Emissivity of building materials
			epsshingles = .91;  								// this could be a user input
			if(rooftype == 2 || rooftype == 3) {
				epsshingles = .9;
			}

			// make a loop to do this
			// North Sheathing
			R2t4 = (1 - EPS1) / EPS1 + 1 / F2t4 + (1 - EPS1) / EPS1 * (A2 / A4);
			R2t8 = (1 - EPS1) / EPS1 + 1 / F2t8 + (1 - EPS1) / EPS1 * (A2 / A8);
			R2t11 = (1 - EPS1) / EPS1 + 1 / F2t11 + (1 - EPS1) / EPS1 * (A2 / (A11 / 3));
			R2t14 = (1 - EPS1) / EPS1 + 1 / F2t14 + (1 - EPS1) / EPS1 * (A2 / (A14 / 3));

			HR2t4 = SIGMA * (told[1] + told[3]) * (pow(told[1], 2) + pow(told[3], 2)) / R2t4;
			HR2t8 = SIGMA * (told[1] + told[7]) * (pow(told[1], 2) + pow(told[7], 2)) / R2t8;
			HR2t11 = SIGMA * (told[1] + told[10]) * (pow(told[1], 2) + pow(told[10], 2)) / R2t11;
			HR2t14 = SIGMA * (told[1] + told[13]) * (pow(told[1], 2) + pow(told[13], 2)) / R2t14;

			// South Sheathing
			R4t2 = (1 - EPS1) / EPS1 + 1 / F4t2 + (1 - EPS1) / EPS1 * (A4 / A2);
			R4t8 = (1 - EPS1) / EPS1 + 1 / F4t8 + (1 - EPS1) / EPS1 * (A4 / A8);
			R4t11 = (1 - EPS1) / EPS1 + 1 / F4t11 + (1 - EPS1) / EPS1 * (A4 / (A11 / 3));
			R4t14 = (1 - EPS1) / EPS1 + 1 / F4t14 + (1 - EPS1) / EPS1 * (A4 / (A14 / 3));

			HR4t2 = SIGMA * (told[3] + told[1]) * (pow(told[3], 2) + pow(told[1], 2)) / R4t2;
			HR4t8 = SIGMA * (told[3] + told[7]) * (pow(told[3], 2) + pow(told[7], 2)) / R4t8;
			HR4t11 = SIGMA * (told[3] + told[10]) * (pow(told[3], 2) + pow(told[10], 2)) / R4t11;
			HR4t14 = SIGMA * (told[3] + told[13]) * (pow(told[3], 2) + pow(told[13], 2)) / R4t14;

			// Attic Floor
			R8t4 = (1 - EPS1) / EPS1 + 1 / F8t4 + (1 - EPS1) / EPS1 * (A8 / A4);
			R8t2 = (1 - EPS1) / EPS1 + 1 / F8t2 + (1 - EPS1) / EPS1 * (A8 / A2);
			// R8t11 = (1 - EPS1) / EPS1 + 1 / F8t11 + (1 - EPS1) / EPS1 * (A8 / A11)
			// R8t14 = (1 - EPS1) / EPS1 + 1 / F8t14 + (1 - EPS1) / EPS1 * (A8 / A14)

			HR8t4 = SIGMA * (told[7] + told[3]) * (pow(told[7], 2) + pow(told[3], 2)) / R8t4;
			HR8t2 = SIGMA * (told[7] + told[1]) * (pow(told[7], 2) + pow(told[1], 2)) / R8t2;
			// HR8t11 = SIGMA * (told(8) + told(11)) * (told(8) ^ 2 + told(11) ^ 2) / R8t11
			// HR8t14 = SIGMA * (told(8) + told(14)) * (told(8) ^ 2 + told(14) ^ 2) / R8t14

			// Return Ducts (note, No radiative exchange w/ supply ducts)
			R11t4 = (1 - EPS1) / EPS1 + 1 / F11t4 + (1 - EPS1) / EPS1 * (A11 / A4);
			// R11t8 = (1 - EPS1) / EPS1 + 1 / F11t8 + (1 - EPS1) / EPS1 * (A11 / A8)
			R11t2 = (1 - EPS1) / EPS1 + 1 / F11t2 + (1 - EPS1) / EPS1 * (A11 / A2);

			HR11t4 = SIGMA * (told[10] + told[3]) * (pow(told[10], 2) + pow(told[3], 2)) / R11t4;
			// HR11t8 = SIGMA * (told(11) + told(8)) * (told(11) ^ 2 + told(8) ^ 2) / R11t8
			HR11t2 = SIGMA * (told[10] + told[1]) * (pow(told[10], 2) + pow(told[1], 2)) / R11t2;

			// Supply Ducts (note, No radiative exchange w/ return ducts)
			R14t4 = (1 - EPS1) / EPS1 + 1 / F14t4 + (1 - EPS1) / EPS1 * (A14 / A4);
			// R14t8 = (1 - EPS1) / EPS1 + 1 / F14t8 + (1 - EPS1) / EPS1 * (A14 / A8)
			R14t2 = (1 - EPS1) / EPS1 + 1 / F14t2 + (1 - EPS1) / EPS1 * (A14 / A2);

			HR14t4 = SIGMA * (told[13] + told[3]) * (pow(told[13], 2) + pow(told[3], 2)) / R14t4;
			// HR14t8 = SIGMA * (told(14) + told(8)) * (told(14) ^ 2 + told(8) ^ 2) / R14t8
			HR14t2 = SIGMA * (told[13] + told[1]) * (pow(told[13], 2) + pow(told[1], 2)) / R14t2;		
		}

		// left overs
		// underside of ceiling
		R7 = (1 - EPS1) / EPS1 + 1 + (1 - EPS1) / EPS1 * (A7 / A13);

		// FOR RAD COEF LAST HOUSE TEMP USED AS INITIAL
		// ESTIMATE OF CEIL TEMP
		hr7 = SIGMA * (told[6] + told[12] * (pow(told[6], 2) + pow(told[12], 2)) / R7);
		Beta = pitch;                                				// ROOF PITCH
		FRS = (1 - sc) * (180 - Beta) / 180;      					// ROOF-SKY SHAPE FACTOR

		if(sc < 1) {
			RS5 = (1 - epsshingles) / epsshingles + 1 / FRS;
			HRS5 = SIGMA * (told[4] + TSKY) * (pow(told[4], 2) + pow(TSKY, 2)) / RS5;
		} else {
			HRS5 = 0;
		}

		FG5 = 1 - FRS;                            					// ROOF-GROUND SHAPE FACTOR
		TGROUND = Tout;                           					// ASSUMING GROUND AT AIR TEMP
		RG5 = (1 - epsshingles) / epsshingles + 1 / FG5;
		HRG5 = SIGMA * (told[4] + TGROUND) * (pow(told[4], 2) + pow(TGROUND, 2)) / RG5;

		// asphalt shingles
		if(rooftype == 1) {
			alpha5 = .92;
			alpha3 = .92;
		} else if(rooftype == 2) {
			// red clay tile - edited for ConSol to be light brown concrete
			alpha5 = .58; 											// .67
			alpha3 = .58; 											// .67
		} else if(rooftype == 3) {
			// low coating clay tile
			alpha5 = .5;
			alpha3 = .5;
		} else if(rooftype == 4) {
			// asphalt shingles  & white coating
			alpha5 = .15;
			alpha3 = .15;
		}

		// South Sheathing
		if(sc < 1) {
			RS3 = (1 - epsshingles) / epsshingles + 1 / FRS;
			HRS3 = SIGMA * (told[2] + TSKY) * (pow(told[2], 2) + pow(TSKY, 2)) / RS3;
		} else {
			HRS3 = 0;
		}

		FG3 = 1 - FRS;                            					// ROOF-GROUND SHAPE FACTOR
		RG3 = (1 - epsshingles) / epsshingles + 1 / FG3;
		HRG3 = SIGMA * (told[2] + TGROUND) * (pow(told[2], 2) + pow(TGROUND, 2)) / RG3;
		
		// NODE 1 IS ATTIC AIR
		if(Mceil >= 0) {
			// flow from attic to house
			A[0][0] = M1 * cp1 / dtau + H14 * A14 / 2 + H11 * A11 / 2 + H8 * A8 + H6 * A6 + Mceil * cp1 + msupahoff * cp15 + mretahoff * cp12 + H4 * A4 + H2 * A2 + A9 * H9 - matticenvout * cp1 - mretleak * cp1;
			b[0] = M1 * cp1 * told[0] / dtau + matticenvin * cp1 * Tout + MSUPLEAK * cp1 * toldcur[14];
		} else {
			// flow from house to attic
			A[0][0] = M1 * cp1 / dtau + H14 * A14 / 2 + H11 * A11 / 2 + H8 * A8 + H6 * A6 + H4 * A4 + H2 * A2 + A9 * H9 - matticenvout * cp1 - mretleak * cp1;
			b[0] = M1 * cp1 * told[0] / dtau - Mceil * cp1 * toldcur[15] - msupahoff * cp15 * toldcur[14] - mretahoff * cp12 * toldcur[11] + matticenvin * cp1 * Tout + MSUPLEAK * cp15 * toldcur[14];
		}

		A[0][1] = -H2 * A2;
		A[0][3] = -H4 * A4;
		A[0][5] = -H6 * A6;
		A[0][7] = -H8 * A8;
		A[0][8] = -H9 * A9;
		A[0][10] = -H11 * A11 / 2;
		A[0][13] = -H14 * A14 / 2;

		if(location == 1) {
			// ducts in house
			if(Mceil >= 0) {
				// flow from attic to house
				A[0][0] = M1 * cp1 / dtau + H8 * A8 + H6 * A6 + Mceil * cp1 + msupahoff * cp15 + mretahoff * cp12 + H4 * A4 + H2 * A2 + A9 * H9 - matticenvout * cp1 - mretleak * cp1;
				b[0] = M1 * cp1 * told[0] / dtau + matticenvin * cp1 * Tout + MSUPLEAK * cp1 * toldcur[14];
			} else {
				// flow from house to attic
				A[0][0] = M1 * cp1 / dtau + H8 * A8 + H6 * A6 + H4 * A4 + H2 * A2 + A9 * H9 - matticenvout * cp1 - mretleak * cp1;
				b[0] = M1 * cp1 * told[0] / dtau - Mceil * cp1 * toldcur[15] - msupahoff * cp15 * toldcur[14] - mretahoff * cp12 * toldcur[11] + matticenvin * cp1 * Tout + MSUPLEAK * cp15 * toldcur[14];
			}
			// no duct surface conduciton losses
			A[0][10] = 0;
			A[0][13] = 0;
		}


		// NODE 2 IS INSIDE NORTH SHEATHING
		A[1][0] = -H2 * A2;
		A[1][1] = M2 * cp2 / dtau + H2 * A2 + A2 / Rval2 + HR2t4 * A2 + HR2t8 * A2 + HR2t11 * A2 + HR2t14 * A2;
		b[1] = M2 * cp2 * told[1] / dtau;
		A[1][2] = -A2 / Rval2;
		A[1][3] = -HR2t4 * A2;
		A[1][7] = -HR2t8 * A2;
		A[1][10] = -HR2t11 * A2;
		A[1][13] = -HR2t14 * A2;

		if(location == 1) {
			// ducts in house
			A[1][1] = M2 * cp2 / dtau + H2 * A2 + A2 / Rval2 + HR2t4 * A2 + HR2t8 * A2;			// + HR2t11 * A2 + HR2t14 * A2
			b[1] = M2 * cp2 * told[1] / dtau;
			A[1][10] = 0;													// -HR2t11 * A2
			A[1][13] = 0;													// -HR2t14 * A2
		}

		// NODE 3 IS OUTSIDE NORTH SHEATHING
		A[2][1] = -A2 / Rval3;
		A[2][2] = M3 * cp3 / dtau + H3 * A3 + A2 / Rval3 + HRS3 * A2 + HRG3 * A2;
		b[2] = M3 * cp3 * told[2] / dtau + H3 * A3 * Tout + A2 * nsolrad * alpha3 + HRS3 * A2 * TSKY + HRG3 * A2 * TGROUND;

		// NODE 4 IS INSIDE SOUTH SHEATHING
		A[3][0] = -H4 * A4;
		A[3][1] = -HR4t2 * A4;
		A[3][3] = M4 * cp4 / dtau + H4 * A4 + A4 / Rval4 + HR4t2 * A4 + HR4t8 * A4 + HR4t11 * A4 + HR4t14 * A4;
		b[3] = M4 * cp4 * told[3] / dtau;
		A[3][4] = -A4 / Rval4;
		A[3][7] = -HR4t8 * A4;
		A[3][10] = -HR4t11 * A4;
		A[3][13] = -HR4t14 * A4;

		if(location == 1) {
			A[3][3] = M4 * cp4 / dtau + H4 * A4 + A4 / Rval4 + HR4t2 * A4 + HR4t8 * A4;			// + HR4T11 * A4 + HR4T14 * A4
			b[3] = M4 * cp4 * told[3] / dtau;
			A[3][10] = 0;													// -HR4T11 * A4
			A[3][13] = 0;													// -HR4T14 * A4
		}

		// NODE 5 IS OUTSIDE SOUTH SHEATHING
		A[4][3] = -A4 / Rval5;
		A[4][4] = M5 * cp5 / dtau + H5 * A5 + A4 / Rval5 + HRS5 * A4 + HRG5 * A4;
		b[4] = M5 * cp5 * told[4] / dtau + H5 * A5 * Tout + A4 * ssolrad * alpha5 + HRS5 * A4 * TSKY + HRG5 * A4 * TGROUND;

		// NODE 6 IS MASS OF WOOD IN ATTIC I.E. JOISTS AND TRUSSES
		A[5][0] = -H6 * A6;
		A[5][5] = M6 * cp6 / dtau + H6 * A6;
		b[5] = M6 * cp6 * told[5] / dtau;

		// NODE  7 ON INSIDE OF CEILING
		A[6][6] = M7 * cp7 / dtau + H7 * A7 + hr7 * A7 + A7 / Rval7;
		b[6] = M7 * cp7 / dtau * told[6];
		A[6][7] = -A7 / Rval7;
		A[6][15] = -H7 * A7;
		A[6][12] = -hr7 * A7;

		if(location == 1) {
			// ducts in house
			A[6][6] = M7 * cp7 / dtau + H7 * A7 + hr7 * A7 + A7 / Rval7;
		}

		// NODE 8 ON ATTIC FLOOR
		A[7][0] = -H8 * A8;
		A[7][1] = -HR8t2 * A8;
		A[7][3] = -HR8t4 * A8;
		A[7][6] = -A8 / Rval8;
		A[7][7] = M8 * cp8 / dtau + H8 * A8 + HR8t2 * A8 + HR8t4 * A8 + A8 / Rval8;				// + HR8t11 * A8 + HR8t14 * A8
		b[7] = M8 * cp8 / dtau * told[7];

		// NODE 9 IS INSIDE ENDWALLS THAT ARE BOTH LUMPED TOGETHER
		A[8][0] = -H9 * A9;
		A[8][8] = M9 * cp9 / dtau + H9 * A9 + A9 / Rval9;
		A[8][9] = -A9 / Rval9;
		b[8] = M9 * cp9 * told[8] / dtau;

		// NODE 10 IS OUTSIDE ENDWALLS THAT ARE BOTH LUMPED TOGETHER
		A[9][8] = -A10 / Rval10;
		A[9][9] = M10 * cp10 / dtau + H10 * A10 + A10 / Rval10;
		b[9] = M10 * cp10 * told[9] / dtau + H10 * A10 * Tout;

		// NODE 11 Exterior Return Duct Surface
		// Remember that the fluid properties are evaluated at a constant temperature
		// therefore, the convection on the inside of the ducts is
		A[10][0] = -A11 * H11 / 2;
		A[10][1] = -A11 * HR11t2 / 3;
		A[10][3] = -A11 * HR11t4 / 3;
		A[10][10] = M11 * cp11 / dtau + H11 * A11 / 2 + A12 / (Rval11 + 1 / HI11) + A11 / 3 * HR11t2 + A11 / 3 * HR11t4;
		b[10] = M11 * cp11 * told[10] / dtau;
		A[10][11] = -A12 / (Rval11 + 1 / HI11);

		if(location == 1) {
			// ducts in house
			A[10][0] = 0;
			A[10][1] = 0;
			A[10][3] = 0;
			A[10][10] = M11 * cp11 / dtau + H11 * A11 + A12 / (Rval11 + 1 / HI11);
			b[10] = M11 * cp11 * told[10] / dtau;
			A[10][15] = -A11 * H11;
		}

		// NODE 12 Air in return duct
		A[11][10] = -A12 / (Rval11 + 1 / HI11);

		if(Mceil >= 0) {
			// flow from attic to house
			A[11][11] = M12 * cp12 / dtau + A12 / (Rval11 + 1 / HI11) + MAH * cp12 + mretahoff * cp12;
			b[11] = M12 * cp12 * told[11] / dtau + mretahoff * cp1 * toldcur[0] - mretleak * cp1 * toldcur[0] - mretreg * cp1 * toldcur[15] - mfancyc * cp1 * Tout - mhrv * cp16 * (.3 * Tout + .7 * told[15]) - merv * cp16 * (.4 * Tout + .6 * told[15]);
		} else {
			// flow from house to attic
			A[11][11] = M12 * cp12 / dtau + A12 / (Rval11 + 1 / HI11) + MAH * cp12 - mretahoff * cp12;
			b[11] = M12 * cp12 * told[11] / dtau - mretahoff * cp16 * toldcur[15] - mretleak * cp1 * toldcur[0] - mretreg * cp1 * toldcur[15] - mfancyc * cp1 * Tout - mhrv * cp16 * (.3 * Tout + .7 * told[15]) - merv * cp16 * (.4 * Tout + .6 * told[15]);
		}

		// node 13 is the mass of the structure of the house that interacts

		// with the house air to increase its effective thermal mass
		// August 99, 95% of solar gain goes to house mass, 5% to house air now
		// Solar gain also calculate more carefully

		for(int i=0; i < 4; i++) {
			S = ((i+1) - 1) * pi / 2;
			cphi = (SBETA * sin(L) - sin(dec)) / CBETA / cos(L);
			cphi2 = pow(cphi, 2);

			if(cphi == 1) {
				sphi = sqrt(1 - cphi2);
			} else {
				sphi = 0;
			}

			if(cphi == 0) {
				if(sphi > 0) {
					phi = 3.14159 / 2;
				} else {
					phi = -3.14159 / 2;
				}
			} else if(cphi == 1) {
				phi = 0;
			} else {
				phi = atan(sphi / cphi);
			}

			Gamma = phi - S;
			ct = CBETA * cos(Gamma);

			if(ct <= -.2) {
				incsolar[i] = Csol * idirect * .45;
			} else {
				incsolar[i] = Csol * idirect * (.55 + .437 * CBETA + .313 * pow(CBETA, 2));
			}
		}

		incsolarS = incsolar[0];
		incsolarW = incsolar[1];
		incsolarN = incsolar[2];
		incsolarE = incsolar[3];

		solgain = shadcoef * (windowS * incsolarS + windowWE / 2 * incsolarW + windowN * incsolarN + windowWE / 2 * incsolarE);

		// incident solar radiations averaged for solair temperature
		incsolarvar = (incsolarN + incsolarS + incsolarE + incsolarW) / 4;
		
		A[12][12] = M13 * cp13 / dtau + H13 * A13 + hr7 * A7;
		A[12][15] = -H13 * A13;
		A[12][6] = -hr7 * A7;
		b[12] = M13 * cp13 * told[12] / dtau + .95 * solgain;

		// NODE 14
		A[13][0] = -A14 * H14 / 2;
		A[13][1] = -A14 * HR14t2 / 3;
		A[13][3] = -A14 * HR14t4 / 3;
		A[13][13] = M14 * cp14 / dtau + H14 * A14 / 2 + A15 / (Rval14 + 1 / HI14) + A14 * HR14t2 / 3 + A14 / 3 * HR14t4;
		b[13] = M14 * cp14 * told[13] / dtau;
		A[13][14] = -A15 / (Rval14 + 1 / HI14);
		if(location == 1) {
			// ducts in house
			A[13][0] = 0;					// -A14 * H14 / 2
			A[13][1] = 0;					// -A14 * HR14t2
			A[13][3] = 0;					// -A14 * HR14t4
			A[13][13] = M14 * cp14 / dtau + H14 * A14 + A15 / (Rval14 + 1 / HI14);
			A[13][15] = -A14 * H14;
		}

		// NODE 15 Air in SUPPLY duct
		// capacity is AC unit capcity in Watts
		// this is a sensible heat balance, the moisture is balanced in a separate routine.

		A[14][13] = -A15 / (Rval14 + 1 / HI14);

		if(Mceil >= 0) {
			// flow from attic to house
			A[14][14] = M15 * cp15 / dtau + A15 / (Rval14 + 1 / HI14) + MSUPREG * cp15 + MSUPLEAK * cp15 + msupahoff * cp15;
			b[14] = M15 * cp15 * told[14] / dtau - capacityc + capacityh + evapcap + MAH * cp12 * toldcur[11] + msupahoff * cp1 * toldcur[0];
		} else {
			// flow from house to attic
			A[14][14] = M15 * cp15 / dtau + A15 / (Rval14 + 1 / HI14) + MSUPREG * cp15 + MSUPLEAK * cp15 - msupahoff * cp15;
			b[14] = M15 * cp15 * told[14] / dtau - capacityc + capacityh + evapcap + MAH * cp12 * toldcur[11] - msupahoff * cp16 * toldcur[15];
		}

		// NODE 16 AIR IN HOUSE
		// use solair tmeperature for house UA
		// the .03 is from 1993 AHSRAE Fund. SI 26.5
		tsolair = incsolarvar * .03 + Tout;

		if(Mceil >= 0) {
			// flow from attic to house
			A[15][15] = M16 * cp16 / dtau + H7 * A7 - mretreg * cp16 - mhouseout * cp16 + H13 * A13 + UA;
			b[15] = M16 * cp16 * told[15] / dtau + mhousein * cp16 * Tout + UA * tsolair + .05 * solgain + MSUPREG * cp1 * toldcur[14] + Mceil * cp1 * toldcur[0] + msupahoff * cp15 * toldcur[14] + mretahoff * cp12 * toldcur[11] + intgain;
		} else {
			// flow from house to attic
			A[15][15] = M16 * cp16 / dtau + H7 * A7 - Mceil * cp16 - msupahoff * cp16 - mretahoff * cp16 - mretreg * cp16 - mhouseout * cp16 + H13 * A13 + UA;
			b[15] = M16 * cp16 * told[15] / dtau + mhousein * cp16 * Tout + UA * tsolair + .05 * solgain + MSUPREG * cp1 * toldcur[14] + intgain;
		}

		A[15][6] = -H7 * A7;
		A[15][12] = -H13 * A13;

		if(location == 1) {
			// ducts in house
			if(Mceil >= 0) {
				// flow from attic to house
				A[15][15] = M16 * cp16 / dtau + H7 * A7 + A11 * H11 + A14 * H14 - mretreg * cp16 - mhouseout * cp16 + H13 * A13 + UA;
				b[15] = M16 * cp16 * told[15] / dtau + mhousein * cp16 * Tout + UA * tsolair + .05 * solgain + MSUPREG * cp1 * toldcur[14] + Mceil * cp1 * toldcur[0] + msupahoff * cp15 * toldcur[14] + mretahoff * cp12 * toldcur[11];
			} else {
				// flow from house to attic
				A[15][15] = M16 * cp16 / dtau + H7 * A7 + A11 * H11 + A14 * H14 - Mceil * cp16 - msupahoff * cp16 - mretahoff * cp16 - mretreg * cp16 - mhouseout * cp16 + H13 * A13 + UA;
				b[15] = M16 * cp16 * told[15] / dtau + mhousein * cp16 * Tout + UA * tsolair + .05 * solgain + MSUPREG * cp1 * toldcur[14];
			}
			A[15][10] = -A11 * H11;
			A[15][13] = -A14 * H14;
		}

		asize = sizeof(A)/sizeof(A[0]);
		asize2 = sizeof(A[0])/sizeof(A[0][0]);
		
		ERRCODE = MatSEqn(A, b, asize, asize2, bsize);

		// PRINT HR2T14, HR2t8, HR2t11, HRS3, HRG3, HR4t2, HR4t8, Hr4t11, Hr4t11, HRS5, HRG5, hr7, HR8t2, HR8t4, HR8t11, HR11t2, HR11t4, HR11t8, HR14t2, HR14t4
		
		if(abs(b[0] - toldcur[0]) < .1) {
			break;
		} else {
			for(int i=0; i < 16; i++) {
				toldcur[i] = b[i];
			}
		}
	} // END of DO LOOP
}

void moisture ( 
	double* HR, 
	double* hrold, 
	double& Mw1, 
	double& Mw2, 
	double& Mw3, 
	double& Mw4, 
	double& Mw5, 
	double& dtau, 
	double& matticenvout, 
	double& Mceil, 
	double& msupahoff, 
	double& mretahoff, 
	double& matticenvin, 
	double& HROUT, 
	double& MSUPLEAK, 
	double& MAH, 
	double& mretreg, 
	double& mretleak, 
	double& MSUPREG, 
	double& latcap, 
	double& mhousein, 
	double& mhouseout, 
	double& latload, 
	double& mhrv, 
	double& mfancyc, 
	double& merv, 
	double& MWha 
	) {
		double Q[5];
		double R[5];

		// 5 nodes
		// Node 1 is attic air
		// Node 2 is return air
		// Node 3 is supply air
		// Node 4 is house air
		// Node 5 is house materials that interact with house air only

		// this routine takes humidity ratios and mass flows and calculates W in each zone

		// Attic Air
		if(Mceil >= 0) {
			// flow from attic to house
			Q[0] = Mw1 / dtau - matticenvout - mretleak + Mceil + msupahoff + mretahoff;
			R[0] = Mw1 * hrold[0] / dtau + matticenvin * HROUT + MSUPLEAK * hrold[2];
		} else {
			// flow from house to attic
			Q[0] = Mw1 / dtau - matticenvout - mretleak;
			R[0] = Mw1 * hrold[0] / dtau + matticenvin * HROUT + MSUPLEAK * hrold[2] - Mceil * hrold[3] - msupahoff * hrold[2] - mretahoff * hrold[1];
		}

		// Return Air
		if(Mceil >= 0) { // flow from attic to house
			Q[1] = Mw2 / dtau + MAH + mretahoff;
			R[1] = Mw2 * hrold[1] / dtau - mretleak * hrold[0] + mretahoff * hrold[0] - mretreg * hrold[3] - mfancyc * HROUT - merv * (.64 * HROUT + .36 * hrold[3]) - mhrv * HROUT;
		} else { // flow from house to attic
			Q[1] = Mw2 / dtau - mretahoff + MAH;
			R[1] = Mw2 * hrold[1] / dtau - mretahoff * hrold[3] - mretleak * hrold[0] - mretreg * hrold[3] - mfancyc * HROUT - merv * (.64 * HROUT + .36 * hrold[3]) - mhrv * HROUT;
		}

		// Supply Air
		if(Mceil >= 0) { // flow from attic to house
			Q[2] = Mw3 / dtau + MSUPLEAK + MSUPREG + msupahoff;
			R[2] = Mw3 * hrold[2] / dtau + MAH * hrold[1] + msupahoff * hrold[0] - latcap / 2501000;
		} else { // flow from house to attic
			Q[2] = Mw3 / dtau + MSUPLEAK + MSUPREG - msupahoff;
			R[2] = Mw3 * hrold[2] / dtau + MAH * hrold[1] - msupahoff * hrold[3] - latcap / 2501000;
		}

		// House Air
		if(Mceil >= 0) { // flow from attic to house
			Q[3] = Mw4 / dtau - mhouseout - mretreg + MWha;
			R[3] = Mw4 * hrold[3] / dtau + mhousein * HROUT + MSUPREG * hrold[2] + Mceil * hrold[0] + msupahoff * hrold[2] + mretahoff * hrold[1] + latload + MWha * hrold[4];
		} else { // flow from house to attic
			Q[3] = Mw4 / dtau - mhouseout - mretreg - Mceil - msupahoff - mretahoff + MWha;
			R[3] = Mw4 * hrold[3] / dtau + mhousein * HROUT + MSUPREG * hrold[2] + latload + MWha * hrold[4];
		}

		// House furnishings/storage
		Q[4] = Mw5 / dtau + MWha;
		R[4] = Mw5 * hrold[4] / dtau + MWha * hrold[3];
		
		for(int P=0; P < 5; P++) {
			HR[P] = R[P] / Q[P];
			// PRINT HR(P);
		}
}

void houseleak (
	int& ahflag,
	double& flag, 
	double& U, 
	double& windangle, 
	double& Tin, 
	double& Tattic, 
	double& Tout, 
	double& C, 
	double& n, 
	double& h, 
	double& R, 
	double& X, 
	int& Nrflue, 
	Chimtype* Chim, 
	double* Wallfraction, 
	double* FloorFraction, 
	double* Sw, 
	double& Swrfl, 
	int& Nwindoor, 
	windoortype* Windoor, 
	int& Nfan, 
	fantype* fan, 
	int& Npipe, 
	pipetype* Pipe, 
	double& min, 
	double& mout, 
	double& Pint, 
	double& Mflue, 
	double& Mceil, 
	double* Mfloor, 
	double& Cattic, 
	double& dPflue, 
	double& dPceil, 
	double& dPfloor, 
	int& Crawl, 
	double& Hfloor, 
	string& row, 
	double* Soffitfraction, 
	double& Patticint, 
	double* Cpwall, 
	double& rhoref, 
	double& MSUPREG, 
	double& MAH, 
	double& mretleak, 
	double& MSUPLEAK, 
	double& mretreg, 
	double& mhousein, 
	double& mhouseout, 
	double& supC, 
	double& supn, 
	double& retC, 
	double& retn, 
	double& msupahoff, 
	double& mretahoff, 
	double& aeq 
	) {
		double dtheta = 11.3;
		int nofirst = 0;

		double P = 0.17;							// P is the power law exponent of the wind speed profile at the building site (see LBNL 42361 Walker & Wilson 1998 p.15)
		double Cpwallvar = 0;						// This variable replaces the non-array Cpwall var
		double CPvar = 0;							// This variable replaces the non-array CP var
		double Bovar = 0;							// This variable is to replace Bo[-1] and to account for Bo(0) in BASIC version

		double Mwall[4];
		double Mwallin[4];
		double Mwallout[4];
		double Bo[4];
		double CP[4][4];
		double dPint;
		double rhoi;
		double rhoo;
		double rhoa;
		double nflue;
		double dPwind;
		double dPtemp;
		double dPwalltop = 0;
		double dPwallbottom = 0;
		double Cproof;
		double Cpwalls;
		double Cpattic;		
		double Cpfloor;
		double Cwall;
		double Cfloor;
		
		Mflue = 0;
		Mceil = 0;

		for(int i=0; i < 4; i++) {
			Mfloor[i] = 0;			
			Mwall[i] = 0;
			Mwallin[i] = 0;
			Mwallout[i] = 0;
			Cpwall[i] = 0;
			Bo[i] = 0;
			for(int j=0; j < 4; j++) {
				CP[i][j] = 0;
			}
		}

		rhoi = rhoref * 293 / Tin;
		rhoo = rhoref * 293 / Tout;
		rhoa = rhoref * 293 / Tattic;

		// the following are pressure differences common to all the flow equations
		dPwind = rhoo / 2 * pow(U,2);
		dPtemp = rhoo * g * (Tin - Tout) / Tin;

		// IF dPwind = 0 AND dPtemp = 0 THEN
		//	GOTO nopressure
		//END IF

		// the following are some typical pressure coefficients for rectangular houses
		
		Cproof = -.4;

		for(int i=0; i < 4; i++) {
			CP[i][0] = .6;
			CP[i][1] = -.3;
		}
		// here the variation of each wall Cp with wind angle is accounted for:
		// for row houses:

		if(strUppercase(row) == "R") {
			CP[0][2] = -.2;
			CP[0][3] = -.2;
			CP[1][2] = -.2;
			CP[1][3] = -.2;
			CP[2][2] = -.65;
			CP[2][3] = -.65;
			CP[3][2] = -.65;
			CP[3][3] = -.65;
		} else {
			// for isolated houses
			for(int i=0; i < 4; i++) {
				CP[i][2] = -.65;
				CP[i][3] = -.65;
			}
		}

		cptheta(CP, windangle, Cpwall);

		Cpwalls = 0;
		Cpattic = 0;

		for(int i=0; i < 4; i++) {
			Cpattic = Cpattic + pow(Sw[i], 2) * Cpwall[i] * Soffitfraction[i];
			Cpwalls = Cpwalls + pow(Sw[i], 2) * Cpwall[i] * Wallfraction[i];			// Shielding weighted Cp
		}

		Cpattic = Cpattic + pow(Swrfl, 2) * Cproof * Soffitfraction[4];

		if(flag < 1) {
			Pint = 0;						// a reasonable first guess
			dPint = 200;
		} else {
			dPint = .25;
		}

		do {
			min = 0;
			mout = 0;
			// ORG: IF Cflue THEN
			if(Nrflue) {					//FF: This IF behaves as if(Nrflue != 0)
				nflue = .5;
				
				// ORG: Flueflow Tin, Fluetemp, Swrfl, dPwind, dPtemp, P, H, Pint, Cflue, nflue, Hflue, Mflue, rhoo, rhoi, dPflue, Tout
				Flueflow(Tin, Swrfl, dPwind, dPtemp, P, h, Pint, nflue, Nrflue, Chim, Mflue, rhoo, rhoi, dPflue, Tout, aeq);

				if(Mflue >= 0) {
					min = min + Mflue;
				} else {
					mout = mout + Mflue;
				}
			}
			
			

			if((R - X) / 2) {
				if(Crawl == 1) {
					// for a crawlspace the flow is put into array position 1
					Cpfloor = Cpwalls;
					Cfloor = C * (R - X) / 2;

					Floorflow3(Cfloor, Cpfloor, dPwind, Pint, C, n, Mfloor[0], rhoo, rhoi, dPfloor, Hfloor, dPtemp);
					
					if(Mfloor[0] >= 0) {
						min = min + Mfloor[0];
					} else {
						mout = mout + Mfloor[0];
					}

				} else {
					for(int i=0; i < 4; i++) {
						Cpfloor = pow(Sw[i], 2) * Cpwall[i];
						Cfloor = C * (R - X) / 2 * FloorFraction[i];

						Floorflow3(Cfloor, Cpfloor, dPwind, Pint, C, n, Mfloor[i], rhoo, rhoi, dPfloor, Hfloor, dPtemp);
						
						if(Mfloor[i] >= 0) {							
							min = min + Mfloor[i];
						} else {
							mout = mout + Mfloor[i];
						}
					}					
				}
			}
			
			if((R + X) / 2) {

				Ceilflow(ahflag, R, X, Patticint, h, dPtemp, dPwind, Pint, C, n, Mceil, Cattic, rhoa, rhoi, dPceil, Tattic, Tin, Tout, rhoo, msupahoff, mretahoff, supC, supn, retC, retn);

				if(Mceil >= 0) {
					min = min + Mceil + msupahoff + mretahoff;
				} else {
					mout = mout + Mceil + msupahoff + mretahoff;
				}
			}

			// the neutral level is calculated for each wall:
			Neutrallevel2(dPtemp, dPwind, Sw, Pint, Cpwall, Bo, h);
			
			if(R < 1) {
				for(int i=0; i < 4; i++) {
					Cpwallvar = pow(Sw[i], 2) * Cpwall[i];
					Cwall = C * (1 - R) * Wallfraction[i];
					
					Wallflow3(Tin, Tout, rhoi, rhoo, Bo[i], Cpwallvar, n, Cwall, h, Pint, dPtemp, dPwind, Mwall[i], Mwallin[i], Mwallout[i], dPwalltop, dPwallbottom, Hfloor);
					
					min = min + Mwallin[i];
					mout = mout + Mwallout[i];
				}
			}

			for(int i=0; i < Nfan; i++) {
				if(fan[i].on == 1) {						// for cycling fans hey will somtimes be off and we don;t want ot include them

					Fanflow(fan[i], rhoo, rhoi);					
					
					if(fan[i].m >= 0) {
						min = min + fan[i].m;
					} else {
						mout = mout + fan[i].m;
					}
				}
			}
			
			for(int i=0; i < Npipe; i++) {

				// FF: This if is to counter the 0 index out of bound present in BASIC version
				// example: if pipe[i].wall = 0 then it would look for Sw[0] in BASIC version, which defaults to 0 while in the C++ version
				// it would try to find Sw[-1] and cause error.

				if(Pipe[i].wall - 1 >= 0)
					CPvar = pow(Sw[Pipe[i].wall - 1], 2) * Cpwall[Pipe[i].wall - 1];
				else
					CPvar = 0;

				Pipeflow(rhoo, rhoi, CPvar, dPwind, dPtemp, Pint, Pipe[i], Tin, Tout);
				
				if(Pipe[i].m >= 0) {
					min = min + Pipe[i].m;
				} else {
					mout = mout + Pipe[i].m;
				}
			}
			
			for(int i=0; i < Nwindoor; i++) {
				if(Windoor[i].wall-1 >= 0) {
					Cpwallvar = pow(Sw[Windoor[i].wall-1], 2) * Cpwall[Windoor[i].wall-1];
					Bovar = Bo[Windoor[i].wall-1];
				} else {
					Cpwallvar = 0;
					Bovar = 0;
				}

				Windoorflow(Tin, Tout, rhoi, rhoo, h, Bovar, Cpwallvar, n, Pint, dPtemp, dPwind, Windoor[i]);
				
				min = min + Windoor[i].min;
				mout = mout + Windoor[i].mout;
			}

			// DUCT MASS FLOWS
			// Msup is flow out of supply registers plus leakage to inside
			// Mret is flow into return registers plus leakage to inside

			min = min + MSUPREG;
			mout = mout + mretreg; // Note Mret should be negative

			Pint = Pint - sgn(min + mout) * dPint;
			dPint = dPint / 2;

		} while (dPint > .0001);

		if(Mceil >= 0) { // flow from attic to house
			mhousein = min - Mceil - MSUPREG - msupahoff - mretahoff;
			mhouseout = mout - mretreg;
		} else {
			mhousein = min - MSUPREG;
			mhouseout = mout - Mceil - mretreg - msupahoff - mretahoff;
		}
		// nopressure:
}

void Atticleak ( 
	double& flag, 
	double& U, 
	double& windangle, 
	double& Tin, 
	double& Tout, 
	double& Tattic, 
	double& Cattic, 
	double& Nattic, 
	double& h, 
	double& Hpeak, 
	double& Swrfl, 
	double* Sw, 
	int& Navent, 
	AventType* Avent, 
	SoffitType* Soffit, 
	double& Matticin, 
	double& Matticout, 
	double& Patticint, 
	double& Mceil, 
	string& row, 
	double* Soffitfraction, 
	double& pitch, 
	string& RR, 
	int& Nafans, 
	fantype* AFan, 
	double& rhoref, 
	double& MSUPREG, 
	double& mretleak, 
	double& MSUPLEAK, 
	double& matticenvin, 
	double& matticenvout, 
	double& dtau, 
	double& msupahoff, 
	double& mretahoff 
) {

	double dtheta = 0;	
	double Matticwall[4];
	double Matticwallin[4];
	double Matticwallout[4];
	double Cpwall[4];
	double Batto[4];
	double Cproof[4];
	double Cppitch[4];
	double CP[4][4];

	double rhoi;
	double rhoo;
	double rhoa;
	double matticfloor;
	double dPwind;
	double dPtemp;
	double dPatticint;
	double Croof;
	double Cpr;
	double Broofo = 0;
	double Mroof = 0;
	double dProoftop = 0;
	double dProofbottom = 0;
	double CPvar;

	double Qnet;						// FF: This variable is set to 0 but then never used.

	dtheta = 11.3;

	for(int i=0; i < 4; i++) {
		Matticwall[i] = 0;
		Matticwallin[i] = 0;
		Matticwallout[i] = 0;
		Cpwall[i] = 0;
		Batto[i] = 0;
		Cproof[i] = 0;
		Cppitch[i] = 0;
		for(int j=0; j < 4; j++) {
			CP[i][j] = 0;
		}
	}

	rhoi = rhoref * 293 / Tin;
	rhoo = rhoref * 293 / Tout;
	rhoa = rhoref * 293 / Tattic;

	matticfloor = -Mceil;

	// the following are pressure differences common to all the flow equations
	dPwind = rhoo / 2 * pow(U, 2);
	dPtemp = rhoo * g * (Tattic - Tout) / Tattic;

	if(dPwind == 0 && dPtemp == 0) {
		Qnet = 0;
		return;
	}

	// the following are some typical pressure coefficients for pitched roofs
	if(pitch < 10) {
		Cproof[0] = -.8;
		Cproof[1] = -.4;
	} else if(pitch > 30) {
		Cproof[0] = .3;
		Cproof[1] = -.5;
	} else {
		Cproof[0] = -.4;
		Cproof[1] = -.4;
	}
	if(strUppercase(row) == "R") {
		Cproof[2] = -.2;
		Cproof[3] = -.2;
	} else {
		// for isolated houses
		Cproof[2] = -.6;
		Cproof[3] = -.6;
	}
	if(strUppercase(RR) == "D") {
		Cproof[2] = Cproof[0];
		Cproof[3] = Cproof[1];
		if(strUppercase(row) == "R") {
			Cproof[0] = -.2;
			Cproof[1] = -.2;
		} else {
			// for isolated houses
			Cproof[0] = -.6;
			Cproof[1] = -.6;
		}
	}

	// ***************************
	// the following are some typical pressure coefficients for rectangular houses


	double Cproofvar = -.4;		// FF: this variable isnt used

	for(int i=0; i < 4; i++) {
		CP[i][0] = .6;
		CP[i][1] = -.3;
	}

	// here the variation of each wall Cp with wind angle is accounted for:
	// for row houses:

	if(strUppercase(row) == "R") {
		CP[0][2] = -.2;
		CP[0][3] = -.2;
		CP[1][2] = -.2;
		CP[1][3] = -.2;
		CP[2][2] = -.65;
		CP[2][3] = -.65;
		CP[3][2] = -.65;
		CP[3][3] = -.65;
	} else {
		// for isolated houses
		for(int i=0; i < 4; i++) {
			CP[i][2] = -.65;
			CP[i][3] = -.65;
		}
	}

	// cptheta CP(), windangle, Cpwall()
	cptheta(CP, windangle, Cpwall);

	// for pitched roof leaks
	RoofCptheta(Cproof, windangle, Cppitch, pitch);
	
	if(flag < 2) {
		Patticint = 0;            // a reasonable first guess
		dPatticint = 25;
	} else {
		dPatticint = .25;
	}

	do {
		Matticin = 0;
		Matticout = 0;
		
		// the following section is for the pitched part of the roof where
		// the two pitched faces are assumed to have the same leakage

		Croof = Cattic * Soffitfraction[4] / 2;
		// for first pitched part either front,above wall 1, or side above wall3
		if(strUppercase(RR) == "D") {
			Cpr = Cppitch[2] * pow(Sw[2], 2);
		} else {
			Cpr = Cppitch[0] * pow(Sw[0], 2);
		}

		// the neutral level is calculated separately for each roof pitch
		Neutrallevel3(dPtemp, dPwind, Patticint, Cpr, Broofo, Hpeak);
		
		// developed from wallflow3:
		roofflow(Tattic, Tout, rhoa, rhoo, Broofo, Cpr, Nattic, Croof, Hpeak, Patticint, dPtemp, dPwind, Mroof, Matticwallin[0], Matticwallout[0], dProoftop, dProofbottom, h);

		Matticin = Matticin + Matticwallin[0];
		Matticout = Matticout + Matticwallout[0];

		// for second pitched part either back, above wall 2, or side above wall 4
		if(strUppercase(RR) == "D") {
			Cpr = Cppitch[3] * pow(Sw[3], 2);
		} else {
			Cpr = Cppitch[1] * pow(Sw[1], 2);
		}

		Neutrallevel3(dPtemp, dPwind, Patticint, Cpr, Broofo, Hpeak);

		// developed from wallflow3:
		roofflow(Tattic, Tout, rhoa, rhoo, Broofo, Cpr, Nattic, Croof, Hpeak, Patticint, dPtemp, dPwind, Mroof, Matticwallin[1], Matticwallout[1], dProoftop, dProofbottom, h);

		Matticin = Matticin + Matticwallin[1];
		Matticout = Matticout + Matticwallout[1];

		for(int i=0; i < Navent; i++) {

			// FF: This if is to counter the 0 index out of bound present in BASIC version
			// example: if Avent[i].wall = 0 then it would look for Sw[0] in BASIC version, which defaults to 0 while in the C++ version
			// it would try to find Sw[-1] and cause error. Original version:
			// CPvar = pow(Sw[Avent[i].wall], 2) * Cppitch[Avent[i].wall];
			if(Avent[i].wall - 1 >= 0)
				CPvar = pow(Sw[Avent[i].wall - 1], 2) * Cppitch[Avent[i].wall - 1];
			else
				CPvar = 0;

			Aventflow(rhoo, rhoa, CPvar, dPwind, dPtemp, Patticint, Avent[i], Tattic, Tout);
			
			if(Avent[i].m >= 0) {
				Matticin = Matticin + Avent[i].m;
			} else {
				Matticout = Matticout + Avent[i].m;
			}
		}
		// note that gable vents are the same as soffits
		for(int i=0; i < 4; i++) {
			CPvar = pow(Sw[i], 2) * Cpwall[i];

			Soffitflow(rhoo, rhoa, CPvar, dPwind, dPtemp, Patticint, Soffit[i], Soffitfraction[i], Cattic, Nattic, Tattic, Tout);

			if(Soffit[i].m >= 0)
				Matticin = Matticin + Soffit[i].m;
			else
				Matticout = Matticout + Soffit[i].m;
		}
		
		for(int i=0; i < Nafans; i++) {
			if(AFan[i].on == 1) {        // for cycling fans hey will somtimes be off and we don;t want ot include them

				AFanflow(AFan[i], rhoo, rhoa);

				if(AFan[i].m >= 0)
					Matticin = Matticin + AFan[i].m;
				else
					Matticout = Matticout + AFan[i].m;
			}
		}

		// Attic floor flow is determined first by HOUSELEAK
		// this fixed flow rate will fix the Patticint that is then
		// passed back to HOUSELEAK as the interior pressure of the attic

		// note that mattic floor has a sign change so that inflow to the attic is positive for a negative Mceil
		if(matticfloor >= 0)
			Matticin = Matticin + matticfloor - msupahoff - mretahoff;
		else
			Matticout = Matticout + matticfloor - msupahoff - mretahoff;

		Matticin = Matticin + MSUPLEAK;
		Matticout = Matticout + mretleak;

		Patticint = Patticint - sgn(Matticin + Matticout) * dPatticint;
		dPatticint = dPatticint / 2;
	} while(dPatticint > .0001);

	if(matticfloor >= 0) {
		matticenvin = Matticin - matticfloor - MSUPLEAK + msupahoff + mretahoff;
		matticenvout = Matticout - mretleak;
	} else {
		matticenvin = Matticin - MSUPLEAK;
		matticenvout = Matticout - matticfloor - mretleak + msupahoff + mretahoff;
	}
}

// This funcition substitutes for BASIC UCASE command, it returns a string that is an uppercase version of original
string strUppercase(string stringvar) {
   for(unsigned int i=0; i < stringvar.length(); i++)
	   stringvar[i] = toupper(stringvar[i]);

   return stringvar;
}

// This funcition substitutes for BASIC SGN command, it returns +1 if variable is positive, -1 if negative, 0 if its 0
int sgn(double sgnvar) {

	return (sgnvar > 0) - (sgnvar < 0);
}


void cptheta(double CP[4][4], double& windangle, double* Cpwall) {
	// this function takes Cps from a single wind angle perpendicular to the
	// upwind wall and finds Cps for all the walls for any wind angle
	
	double CPvar = 0;
	double Theta = 0;
	
	for(int i=0; i < 4; i++) {
		Theta = windangle * pi / 180;
		
		if(i == 1) {
			Theta = Theta - pi;
		} else if(i == 2) {
			Theta = Theta - .5 * pi;
		} else if(i == 3) {
			Theta = Theta - 1.5 * pi;
		}
		
		CPvar = (CP[i][0] + CP[i][1]) * pow(pow(cos(Theta),2),.25);
		CPvar = CPvar + (CP[i][0] - CP[i][1]) * cos(Theta) * pow(abs(cos(Theta)),-.25);
		CPvar = CPvar + (CP[i][2] + CP[i][3]) * pow(pow(sin(Theta),2),2);
		CPvar = CPvar + (CP[i][2] - CP[i][3]) * sin(Theta);
		
		Cpwall[i] = CPvar / 2;
	}
}

void Flueflow(double& Tin, double& Swrfl, double& dPwind, double& dPtemp, double& P,
	double& h, double& Pint, double& nflue, int& Nrflue, Chimtype* Chim, double& Mflue,
	double& rhoo, double& rhoi, double& dPflue, double& Tout, double& aeq) {

		// calculates flow through the flue

		// dkm: internal variable to calculate external variable Mflue
		double massflue = 0;
		double Cpflue;

		for(int i=0; i < Nrflue; i++) { 												// dkm
			Cpflue = -.5 * pow(Swrfl,2) * pow((Chim[i].Hflue / h),(2 * P));

			if(Chim[i].Fluetemp == -99) {
				// dPflue = Pint - dPtemp * Chim(i%).Hflue + dPwind * Cpflue
				// IF dPflue <= 0 THEN  	// flow out through flue
				//        massflue = massflue + (-rhoi * Chim(i%).Cflue * (293.15 / Tout) ^ (3 * nflue - 2) * (-dPflue) ^ nflue)
				// ELSE                     // flow in through flue
				//        dPflue = Pint - dPtemp * H + dPwind * Cpflue
				//        IF dPflue <= 0 THEN dPflue = -dPflue	// a quick check for sign of pressure if nearly neutral
				//        massflue = massflue + (rhoo * Chim(i%).Cflue * (293.15 / Tin) ^ (3 * nflue - 2) * (dPflue) ^ nflue)
				// END IF

				dPflue = Pint - dPtemp * Chim[i].Hflue + dPwind * Cpflue;
				if(dPflue >= 0) { 														// flow in through flue
					// dPflue = Pint - dPtemp * H + dPwind * Cpflue
					// IF dPflue <= 0 THEN dPflue = -dPflue		// a quick check for sign of pressure if nearly neutral
					massflue = massflue + (rhoo * Chim[i].Cflue * pow((293.15 / Tin),(3 * nflue - 2)) * pow(dPflue,nflue));
				} else {																// flow out through flue
					massflue = massflue + (-rhoi * Chim[i].Cflue * pow((293.15 / Tout),(3 * nflue - 2)) * pow(-dPflue,nflue));
				}
			} else {
				// for a heated flue:  driving pressure correction:
				dPflue = Pint - dPtemp * Chim[i].Hflue + dPwind * Cpflue - g * rhoi * Chim[i].Hflue * (1 - Tin / Chim[i].Fluetemp);
				// density-viscosity correction:
				if(dPflue >= 0) {
					massflue = massflue + (pow((293.15 / Chim[i].Fluetemp),(3 * nflue - 2)) * rhoo * Chim[i].Cflue * pow(dPflue,nflue));
				} else {
					massflue = massflue + (pow(-(293.15 / Chim[i].Fluetemp),(3 * nflue - 2)) * rhoi * Chim[i].Cflue * pow(-dPflue,nflue));
				}
			}
		}			// dkm

		Mflue = massflue;
}


void Floorflow3(double& Cfloor, double& Cpfloor, double& dPwind, double& Pint, double& C,
	double& n, double& Mfloor, double& rhoo, double& rhoi, double& dPfloor, double& Hfloor, double& dPtemp) {

		// calculates flow through floor level leaks
		dPfloor = Pint + Cpfloor * dPwind - Hfloor * dPtemp;

		if(dPfloor >= 0)
			Mfloor = rhoo * Cfloor * pow(dPfloor,n);
		else
			Mfloor = -rhoi * Cfloor * pow(-dPfloor,n);
}

void Ceilflow(int& ahflag, double& R, double& X, double& Patticint, double& h, double& dPtemp,
	double& dPwind, double& Pint, double& C, double& n, double& Mceil, double& Cattic, double& rhoa,
	double& rhoi, double& dPceil, double& Tattic, double& Tin, double& Tout, double& rhoo,
	double& msupahoff, double& mretahoff, double& supC, double& supn, double& retC, double& retn) {

		double Cceil;

		// PRINT AHflag; msupahoff; retC; supC; retn; supn
		// calculates flow through the ceiling
		Cceil = C * (R + X) / 2;
		dPceil = Pint - Patticint - rhoo * g * ((Tin - Tout) / Tin - (Tattic - Tout) / Tattic) * h;

		if(dPceil >= 0) {
			Mceil = rhoa * Cceil * pow(dPceil,n);
			if(ahflag == 0) {
				msupahoff = rhoa * supC * pow(dPceil,supn);
				mretahoff = rhoa * retC * pow(dPceil,retn);
				if(Cattic == 0) {
					msupahoff = 0;
					mretahoff = 0;
					Mceil = 0;
				}
			} else {
				msupahoff = 0;
				mretahoff = 0;
			}
		} else {
			Mceil = -rhoi * Cceil * pow(-dPceil,n);
			if(ahflag == 0) {
				msupahoff = -rhoi * supC * pow(-dPceil,supn);
				mretahoff = -rhoi * retC * pow(-dPceil,retn);
				if(Cattic == 0) {
					Mceil = 0;
					msupahoff = 0;
					mretahoff = 0;
				}
			} else {
				msupahoff = 0;
				mretahoff = 0;
			}
		}
}

void Neutrallevel2(double& dPtemp, double& dPwind, double* Sw, double& Pint, double* Cpwall, double* Bo, double& h) {
	// calculates the neutral level for each wall
	for(int i=0; i < 4; i++) {
		if(dPtemp)
			Bo[i] = (Pint + pow(Sw[i], 2) * dPwind * Cpwall[i]) / dPtemp / h;
		else
			Bo[i] = 0;
	}
}

void Wallflow3(double& Tin, double& Tout, double& rhoi, double& rhoo, double& Bo, double& Cpwall,
	double& n, double& Cwall, double& h, double& Pint, double& dPtemp, double& dPwind, double& Mwall,
	double& Mwallin, double& Mwallout, double& dPwalltop, double& dPwallbottom, double& Hfloor) {
		
		// calculates the flow through a wall
		double Hwall = h - Hfloor;
		Mwallin = 0;
		Mwallout = 0;
		
		// dummys changed so Hfloor<>0
		double dummy1 = Pint + dPwind * Cpwall - dPtemp * h;
		double dummy2 = Pint + dPwind * Cpwall - dPtemp * Hfloor;

		dPwalltop = dummy1;		
		dPwallbottom = dummy2;

		if(Tin == Tout) {
			if(dummy2 > 0) {
				Mwallin = rhoo * Cwall * pow(dummy2, n);
				Mwallout = 0;
			} else {
				Mwallin = 0;
				Mwallout = -rhoi * Cwall * pow(-dummy2, n);
			}
		} else {
			if(Tin > Tout) {
				if(Bo <= 0) {
					Mwallin = 0;
					Mwallout = rhoi * Cwall / Hwall / dPtemp / (n + 1) * (dummy1 * pow(abs(dummy1), n) - dummy2 * pow(abs(dummy2), n));
				} else if(Bo >= 1) {
					Mwallin = -rhoo * Cwall / Hwall / dPtemp / (n + 1) * (dummy1 * pow(abs(dummy1), n) - dummy2 * pow(abs(dummy2), n));
					Mwallout = 0;
				} else {
					Mwallin = rhoo * Cwall / Hwall / dPtemp / (n + 1) * dummy2 * pow(abs(dummy2), n);
					Mwallout = rhoi * Cwall / Hwall / dPtemp / (n + 1) * dummy1 * pow(abs(dummy1), n);
				}
			} else {
				if(Bo <= 0) {
					Mwallin = -rhoo * Cwall / Hwall / dPtemp / (n + 1) * (dummy1 * pow(abs(dummy1), n) - dummy2 * pow(abs(dummy2), n));
					Mwallout = 0;
				} else if(Bo >= 1) {
					Mwallin = 0;
					Mwallout = rhoi * Cwall / Hwall / dPtemp / (n + 1) * (dummy1 * pow(abs(dummy1), n) - dummy2 * pow(abs(dummy2), n));
				} else {
					Mwallin = -rhoo * Cwall / Hwall / dPtemp / (n + 1) * dummy1 * pow(abs(dummy1), n);
					Mwallout = -rhoi * Cwall / Hwall / dPtemp / (n + 1) * dummy2 * pow(abs(dummy2), n);
				}
			}
		}

		Mwall = Mwallin + Mwallout;
}

void Fanflow(fantype& fan, double& rhoo, double& rhoi) {
	
	// calculates flow through ventilation fans
	if(fan.q < 0)
		fan.m = rhoi * fan.q;
	else
		fan.m = rhoo * fan.q;
}

void Pipeflow(double& rhoo, double& rhoi, double& CP, double& dPwind, double& dPtemp,
	double& Pint, pipetype& Pipe, double& Tin, double& Tout) {
		
		// calculates flow through pipes
		// changed on NOV 7 th 1990 so Pipe.A is Cpipe
		// changed on Nov 9th 1990 for density & viscosity variation of C
		// this is representeed by the temperature ratios

		Pipe.dP = Pint - dPtemp * Pipe.h + dPwind * CP;

		if(Pipe.dP >= 0)
			Pipe.m = rhoo * Pipe.A * pow((293.15 / Tout), (3 * Pipe.n - 2)) * pow(Pipe.dP, Pipe.n);
		else
			Pipe.m = -rhoi * Pipe.A * pow((293.15 / Tin), (3 * Pipe.n - 2)) * pow(-Pipe.dP, Pipe.n);
}

void Windoorflow(double& Tin, double& Tout, double& rhoi, double& rhoo, double& h, double& Bo,
	double& Cpwall, double& n, double& Pint, double& dPtemp, double& dPwind, windoortype& Windoor) {

		// calculates flow through open doors or windows

		double dT;
		double Awindoor;
		double dummy;
		double dummy1;
		double dummy2;
		double Viscosity;
		double Kwindow;
		double Kwindownew;
		double Topcrit;
		double Bottomcrit;
		double Pwindow;
		double ReH;

		dT = Tin - Tout;
		Windoor.dPbottom = 2 * (dPwind * Cpwall + Pint - Windoor.Bottom * dPtemp) / rhoo;
		Windoor.dPtop = 2 * (dPwind * Cpwall + Pint - Windoor.Top * dPtemp) / rhoo;

		if(dT == 0) {
			Awindoor = Windoor.High * Windoor.Wide;
			dummy = Pint + dPwind * Cpwall;
			if(dummy >= 0) {
				Windoor.min = .6 * Awindoor * sqrt(dummy * rhoo * 2);
				Windoor.mout = 0;
			} else {
				Windoor.min = 0;
				Windoor.mout = -.6 * Awindoor * sqrt(-dummy * rhoi * 2);
			}
		} else {
			dummy1 = Windoor.dPbottom;
			dummy2 = Windoor.dPtop;
			if(Bo * h <= Windoor.Bottom) {
				Kwindow = .6;
				dummy = dummy2 * sqrt(abs(dummy2)) - dummy1 * sqrt(abs(dummy1));
				if(dT > 0) {
					Windoor.min = 0;
					Windoor.mout = sqrt(rhoi * rhoo) * Kwindow * Windoor.Wide * Tin / 3 / g / dT * dummy;
				} else {
					Windoor.min = -rhoo * Kwindow * Windoor.Wide * Tin / 3 / g / dT * dummy;
					Windoor.mout = 0;
				}
			} else if(Bo * h > Windoor.Top) {
				Kwindow = .6;
				dummy = dummy1 * sqrt(abs(dummy1)) - dummy2 * sqrt(abs(dummy2));
				if(dT > 0) {
					Windoor.min = rhoo * Kwindow * Windoor.Wide * Tin / 3 / g / dT * dummy;
					Windoor.mout = 0;
				} else {
					Windoor.min = 0;
					Windoor.mout = -sqrt(rhoo * rhoi) * Kwindow * Windoor.Wide * Tin / 3 / g / dT * dummy;
				}
			} else {
				Viscosity = .0000133 + .0000009 * ((Tout + Tin) / 2 - 273);
				Kwindow = .4 + .0045 * dT;
				Topcrit = Windoor.Top - .1 * Windoor.High;
				Bottomcrit = Windoor.Bottom + .1 * Windoor.High;
				if(Bo * h > Topcrit) {
					Pwindow = Topcrit * dPtemp - dPwind * Cpwall;
				} else if(Bo * h < Bottomcrit) {
					Pwindow = Bottomcrit * dPtemp - dPwind * Cpwall;
				} else {
					Pwindow = Pint;
				}
				dummy1 = 2 * (dPwind * Cpwall + Pwindow - Windoor.Bottom * dPtemp) / rhoo;
				dummy2 = 2 * (dPwind * Cpwall + Pwindow - Windoor.Top * dPtemp) / rhoo;
				do {
					Kwindownew = Kwindow;
					if(dT > 0) {
						Windoor.min = rhoo * Kwindownew * Windoor.Wide * Tin / 3 / g / dT * dummy1 * sqrt(abs(dummy1));
						Windoor.mout = sqrt(rhoi * rhoo) * Kwindownew * Windoor.Wide * Tin / 3 / g / dT * dummy2 * sqrt(abs(dummy2));
					} else {
						Windoor.min = -rhoo * Kwindownew * Windoor.Wide * Tin / 3 / g / dT * dummy2 * sqrt(abs(dummy2));
						Windoor.mout = -sqrt(rhoi * rhoo) * Kwindownew * Windoor.Wide * Tin / 3 / g / dT * dummy1 * sqrt(abs(dummy1));
					}
					ReH = 2 * abs(Windoor.min / rhoo + Windoor.mout / rhoi) / Windoor.Wide / Viscosity;
					Kwindow = (.3 + .6 * .00003 * ReH) / (1 + .00003 * ReH);
				} while(!(abs(Kwindownew - Kwindow) < .001));

				if(Bo * h > Topcrit) {
					Kwindow = (.6 - Kwindow) / (.1 * Windoor.High) * (h * Bo - Topcrit) + Kwindow;
				} else if(Bo * h < Bottomcrit) {
					Kwindow = (.6 - Kwindow) / (.1 * Windoor.High) * (Bottomcrit - h * Bo) + Kwindow;
				}
				if(dT > 0) {
					Windoor.min = rhoo * Kwindow * Windoor.Wide * Tin / 3 / g / dT * dummy1 * sqrt(abs(dummy1));
					Windoor.mout = sqrt(rhoi * rhoo) * Kwindow * Windoor.Wide * Tin / 3 / g / dT * dummy2 * sqrt(abs(dummy2));
				} else {
					Windoor.min = -rhoo * Kwindow * Windoor.Wide * Tin / 3 / g / dT * dummy2 * sqrt(abs(dummy2));
					Windoor.mout = -sqrt(rhoi * rhoo) * Kwindow * Windoor.Wide * Tin / 3 / g / dT * dummy1 * sqrt(abs(dummy1));
				}
			}
		}

		Windoor.m = Windoor.min + Windoor.mout;
}

void RoofCptheta(double* Cproof, double& windangle, double* Cppitch, double& pitch) {
	
	double Theta;
	double C2;
	double C;
	double func;
	double S2;
	double S;

	for(int i=0; i < 4; i++) {
		Theta = windangle * pi / 180;
		if(i == 2) {
			Theta = Theta - pi;
		} else if(i == 3) {
			Theta = Theta - .5 * pi;
		} else if(i == 4) {
			Theta = Theta - 1.5 * pi;
		}

		C2 = pow(cos(windangle * 3.14159 / 180), 2);
		C = cos(windangle * 3.14159 / 180);

		if(C != 0) {
			if(pitch != 28) {
				func = (1 - abs(pow(C, 5))) / 2 * ((28 - pitch) / 28) * pow(abs((28 - pitch) / 28), -.99) + (1 + abs(pow(C, 5))) / 2;
			} else {
				func = (1 + abs(pow(C, 5))) / 2;
			}
			C = C * func;
		}

		S2 = pow(sin(windangle * 3.14159 / 180), 2);
		S = sin(windangle * 3.14159 / 180);

		Cppitch[i] = (Cproof[0] + Cproof[1]) * C2;
		Cppitch[i] = Cppitch[i] + (Cproof[0] - Cproof[1]) * C;
		Cppitch[i] = Cppitch[i] + (Cproof[2] + Cproof[3]) * S2;
		Cppitch[i] = Cppitch[i] + (Cproof[2] - Cproof[3]) * S;
		Cppitch[i] = Cppitch[i] / 2;
	}
}

void Neutrallevel3(double& dPtemp, double& dPwind, double& Patticint, double& Cpr, double& Broofo, double& Hpeak) {
	// calculates the neutral level for the attic
	if(dPtemp != 0)
		Broofo = (Patticint + dPwind * Cpr) / dPtemp / Hpeak;
	else
		Broofo = 0;
}

void roofflow(double& Tattic, double& Tout, double& rhoa, double& rhoo, double& Broofo, double& Cpr,
	double& Nattic, double& Croof, double& Hpeak, double& Patticint, double& dPtemp, double& dPwind,
	double& Mroof, double& Mroofin, double& Mroofout, double& dProoftop, double& dProofbottom, double& H) {
		
		double Hroof;
		double dummy1;
		double dummy2;

		// calculates the flow through the pitched section of the roof
		Mroofin = 0;
		Mroofout = 0;
		Hroof = Hpeak - H;
		
		dummy1 = Patticint + dPwind * Cpr - dPtemp * Hpeak;
		dProoftop = dummy1;
		dummy2 = Patticint + dPwind * Cpr - dPtemp * H;
		dProofbottom = dummy2;

		if(Tattic == Tout) {
			if(dummy2 > 0) {
				Mroofin = rhoo * Croof * pow(dummy2, Nattic);
				Mroofout = 0;
			} else {
				Mroofin = 0;
				Mroofout = -rhoa * Croof * pow(-dummy2, Nattic);
			}
		} else {
			if(Tattic > Tout) {
				if(Broofo <= H / Hpeak) {
					Mroofin = 0;
					Mroofout = rhoa * Croof / Hroof / dPtemp / (Nattic + 1) * (dummy1 * pow(abs(dummy1), Nattic) - dummy2 * pow(abs(dummy2), Nattic));
				} else if(Broofo >= 1) {
					Mroofin = -rhoo * Croof / Hroof / dPtemp / (Nattic + 1) * (dummy1 * pow(abs(dummy1), Nattic) - dummy2 * pow(abs(dummy2), Nattic));
					Mroofout = 0;
				} else {
					Mroofin = rhoo * Croof / Hroof / dPtemp / (Nattic + 1) * dummy2 * pow(abs(dummy2), Nattic);
					Mroofout = rhoa * Croof / Hroof / dPtemp / (Nattic + 1) * dummy1 * pow(abs(dummy1), Nattic);
				}
			} else {
				if(Broofo <= H / Hpeak) {
					Mroofin = -rhoo * Croof / Hroof / dPtemp / (Nattic + 1) * (dummy1 * pow(abs(dummy1), Nattic) - dummy2 * pow(abs(dummy2), Nattic));
					Mroofout = 0;
				} else if(Broofo >= 1) {
					Mroofin = 0;
					Mroofout = rhoa * Croof / Hroof / dPtemp / (Nattic + 1) * (dummy1 * pow(abs(dummy1), Nattic) - dummy2 * pow(abs(dummy2), Nattic));
				} else {
					Mroofin = -rhoo * Croof / Hroof / dPtemp / (Nattic + 1) * dummy1 * pow(abs(dummy1), Nattic);
					Mroofout = -rhoa * Croof / Hroof / dPtemp / (Nattic + 1) * dummy2 * pow(abs(dummy2), Nattic);
				}
			}
		}

		Mroof = Mroofin + Mroofout;
}

void Aventflow(double& rhoo, double& rhoa, double& CP, double& dPwind, double& dPtemp, double& Patticint,
	AventType& Avent, double& Tattic, double& Tout) {

		// calculates flow through attic roof vents
		Avent.dP = Patticint - dPtemp * Avent.h + dPwind * CP;

		if(Avent.dP >= 0)
			Avent.m = rhoo * Avent.A * pow((293.15 / Tout), (3 * Avent.n - 2)) * pow(Avent.dP, Avent.n);
		else
			Avent.m = -rhoa * Avent.A * pow((293.15 / Tattic), (3 * Avent.n - 2)) * pow(-Avent.dP, Avent.n);
}

void Soffitflow(double& rhoo, double& rhoa, double& CP, double& dPwind, double& dPtemp, double& Patticint,
	SoffitType& Soffit, double& Soffitfraction, double& Cattic, double& Nattic, double& Tattic, double& Tout) {
		
		// calculates flow through attic soffit vents and Gable end vents	
		Soffit.dP = Patticint - dPtemp * Soffit.h + dPwind * CP;

		if(Soffit.dP >= 0)
			Soffit.m = rhoo * Soffitfraction * Cattic * pow((293.15 / Tout), (3 * Nattic - 2)) * pow(Soffit.dP, Nattic);
		else
			Soffit.m = -rhoa * Soffitfraction * Cattic * pow((293.15 / Tattic), (3 * Nattic - 2)) * pow(-Soffit.dP, Nattic);
}

void AFanflow(fantype& AFan, double& rhoo, double& rhoa) {
	
	// calculates flow through ventilation fans
	if(AFan.q < 0)
		AFan.m = rhoa * AFan.q;
	else
		AFan.m = rhoo * AFan.q;
}

// ----- MatSEqn definitions -----

int MatSEqn(double A[][ArraySize], double* b, int asize, int asize2, int bsize) {
	// Error codes returned:
	//      0  no error                     -1  matrix not invertible
	//     -2  matrix not square            -3  inner dimensions different
	//     -4  matrix dimensions different  -5  result matrix dimensioned incorrectly
	//     any other codes returned are standard BASIC errors
	// 
	// -------------------------------------------------------------------
	
	int errcode = 0;
	int bserrcode = 0;
	int ERR = 0;
		
	int continuevar = 0;
	
	double x[ArraySize];
	int rpvt[ArraySize], cpvt[ArraySize];
	
	for(int i=0; i < ArraySize; i++) {
		x[i] = 0;
		rpvt[i] = 0;
		cpvt[i] = 0;
	}
	
	errcode = matlu(A, rpvt, cpvt, asize, asize2, continuevar);			// Get LU matrix

	if(continuevar != -1) {
		errcode = (errcode + 5) % 200 - 5;
		cout << "\nMatSeqn continue error: " << errcode << endl;
		return errcode;	
	}
	
	// check dimension of b
	if(asize != bsize) {
		ERR = 197;
		errcode = (ERR + 5) % 200 - 5;
		cout << "\nMatSeqn size error: " << errcode << endl;
		return errcode;
	}
	
	bserrcode = matbs(A, b, x, rpvt, cpvt, asize);						// Backsolve system
	
	for(int i=0; i < asize; i++) {
	   b[i] = x[i];														// Put solution in b for return
	}

	if(errcode !=0) {
		errcode = (errcode + 5) % 200 - 5;
		cout << "\nMatSeqn error: " << errcode << endl;
		return errcode;
	}
	
	errcode = (ERR + 5) % 200 - 5;
	return errcode;	
}

int matlu(double A[][ArraySize], int* rpvt, int* cpvt, int asize, int asize2, int& continuevar) {
	int errcode = 0;
	int tempswap;
	int count;
	int r;

	int c;
	int bestrow;
	int bestcol;
	int rp;
	int cp;

	double rownorm[ArraySize];
	double max;
	double temp;
	double oldmax;

	// Checks if A is square, returns error code if not
	if(asize != asize2) {
		errcode = 198;
		continuevar = 0;
		cout << "\nMatlu error: " << errcode << endl;
		return errcode;
	}

	count = 0;														// initialize count, continue
	continuevar = -1;
	
	for(int row = 0; row < asize; row++) {							// initialize rpvt and cpvt
		rpvt[row] = row;
		cpvt[row] = row;
		rownorm[row] = 0;                							// find the row norms of A()

		for(int col = 0; col < asize; col++) {
			rownorm[row] = rownorm[row] + abs(A[row][col]);
		}

		// if any rownorm is zero, the matrix is singular, set error, exit and do not continue
		if(rownorm[row] == 0) {
			errcode = 199;
			continuevar = 0;
			cout << "\nMatlu error: " << errcode << endl;
			return errcode;
		}
	}
	
	for(int pvt = 0; pvt < (asize-1); pvt++) {
		// Find best available pivot
		// checks all values in rows and columns not already used for pivoting
		// and finds the number largest in absolute value relative to its row norm
		max = 0;
		
		for(int row = pvt; row < asize; row++) {
			r = rpvt[row];
			for(int col = pvt; col < asize; col++) {
				c = cpvt[col];
				temp = abs(A[r][c]) / rownorm[r];
				if(temp > max) {
					max = temp;
					bestrow = row;		    						// save the position of new max
					bestcol = col;
				}
			}
		}
		
		// if no nonzero number is found, A is singular, send back error, do not continue
		if(max == 0) {
			errcode = 199;
			continuevar = 0;
			cout << "\nMatlu error: " << errcode << endl;
			return errcode;
		} else if(pvt > 1 && max < (deps * oldmax)) {				// check if drop in pivots is too much
			errcode = 199;
		}
		
		oldmax = max;
		
		// if a row or column pivot is necessary, count it and permute rpvt or cpvt.
		// Note: the rows and columns are not actually switched, only the order in which they are used.		
		if(rpvt[pvt] != rpvt[bestrow]) {
			count = count + 1;
			
			tempswap = rpvt[pvt];
			rpvt[pvt] = rpvt[bestrow];
			rpvt[bestrow] = tempswap;
		}    
		
		if(cpvt[pvt] != cpvt[bestcol]) {
			count = count + 1;
			
			tempswap = cpvt[pvt];
			cpvt[pvt] = cpvt[bestcol];
			cpvt[bestcol] = tempswap;
		}
		
		//Eliminate all values below the pivot
		rp = rpvt[pvt];
		cp = cpvt[pvt];

		for(int row = (pvt+1); row < asize; row++) {
			r = rpvt[row];
			A[r][cp] = -A[r][cp] / A[rp][cp];						// save multipliers
			
			for(int col = (pvt+1); col < asize; col++) {
				c = cpvt[col];										// complete row operations
				A[r][c] = A[r][c] + A[r][cp] * A[rp][c];				
			}
		}
		
	}
	
	// if last pivot is zero or pivot drop is too large, A is singular, send back error
	if(A[rpvt[asize-1]][cpvt[asize-1]] == 0) {
		errcode = 199;
		continuevar = 0;
		cout << "\nMatlu error: " << errcode << endl;
		return errcode;
	} else if(((abs(A[rpvt[asize-1]][cpvt[asize-1]])) / rownorm[rpvt[asize-1]]) < (deps * oldmax)) {
		// if pivot is not identically zero then continue remains TRUE
		errcode = 199;
	}
	
	if(errcode != 0 && errcode < 199) {
		continuevar = 0;
	}
	
	return errcode;
}

int matbs(double A[][ArraySize], double* b, double* x, int* rpvt, int* cpvt, int asize) {

	int c, r;

	// do row operations on b using the multipliers in L to find Lb
	for(int pvt = 0; pvt < (asize-1); pvt++) {
		c = cpvt[pvt];		
		for(int row = pvt+1 ; row < asize; row++) {
			r = rpvt[row];
			b[r] = b[r] + A[r][c] * b[rpvt[pvt]];
		}
	}

	// backsolve Ux=Lb to find x
	for(int row = asize-1; row >= 0; row--) {
		c = cpvt[row];
		r = rpvt[row];
		x[c] = b[r];
		for(int col = (row+1); col < asize; col++) {
			x[c] = x[c] - A[r][cpvt[col]] * x[cpvt[col]];		
		}
		x[c] = x[c] / A[r][c];
	}

	return 0;
}

