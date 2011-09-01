// RAD :: Remove After Debugging

#include <iostream>
#include <fstream>
//#include <string>
#include <iomanip>		// RAD: so far used only for setprecission() in cmd output
#include <ctime>
#include <direct.h>		//used for chdir command
#include "functions.h"

// RAD: system() is evil, example: remove system("pause")
// RAD: do-while loop is orignally day < 22, it is currently at day < 2 for shorter analysis comparison

using namespace std;

// Main function
int main(int argc, char *argv[]) 
{ 	
	// [input filenames]
	// (weather and building filenames are read from the starting batch file)
	// for file locations different to current program directory use the complete file path and name,
	// also remember the escape characters effect. working example:
	// string batchfile_name = "D:\\foldername\\example.csv";
	_chdir("input\\");
	string batchfile_name = "v5batch.csv";
	string shelterfile_name = "bshelter.dat";
	string fanSchedulefile_name1 = "sched1";
	string fanSchedulefile_name2 = "sched2";
	string fanSchedulefile_name3 = "sched3";
	
	// Original test case has totaldays = 22. Reduced to 1 to simplify validation.
	// total days to run the simulation for:
	int totaldays = 365;
	
	// ----------------- [start]--------------------
	char reading[255];
	int Numsim;
	string Simfile[50], Climzone[50], Flnme[50];

	cout << "\nRC++ [begin]" << endl << endl;

	// [start] reading batch file
	ifstream batchfile(batchfile_name); 
	if(!batchfile) { 
		cout << "Cannot open: " << batchfile_name << endl;
		system("pause");
		return 1; 
	} 
	
	batchfile.getline(reading, 255);
	//sscanf_s(reading, "%d", &Numsim);
	Numsim = atoi(reading);

	cout << "Batch file info:" << endl;
	for(int i=0; i < Numsim; i++) {
		getline(batchfile, Simfile[i]);
		getline(batchfile, Climzone[i]);
		getline(batchfile, Flnme[i]);
				
		cout << "Simfile[" << i << "]: " << Simfile[i] << endl;
		cout << "Climzone[" << i << "]: " << Climzone[i] << endl;
		cout << "Flnme[" << i << "]: " << Flnme[i] << endl << endl;;		
	}
	
	batchfile.close();
	// [end] reading batch file

	for(int jen=0; jen < Numsim; jen++) {		
		// The following commented variables are not used within the code but were the declared in the RC - BASIC version
		// REDIM mcinit(7), mtotalinit(7)
		// REDIM Bo(4)
		// REDIM Mwall(4), Mwallin(4), Mwallout(4)
		// REDIM Matticwall(4), Mattwallin(4), Mattwallout(4), Batto(4)

		windoortype Windoor[10] = {0};
		fantype fan[10] = {0};
		fantype AFan[10] = {0};
		pipetype Pipe[10] = {0};		
		AventType Avent[10] = {0};
		SoffitType Soffit[4] = {0};
		Chimtype Chim[6] = {0};
		
		double Sw[4];
		double FloorFraction[4];
		double Wallfraction[4];
		double Swinit[4][361];		
		double Mfloor[4] = {0,0,0,0};
		double Soffitfraction[5];
		double Cpwall[4] = {0,0,0,0};
		//double mvpow[10];
		double mvpow;
		double b[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		double told[16];
		double hrold[5];
		double HR[5];		
		double theat[24];
		double tcool[24];
		
		// Zeroing the variables to create the sums for the .ou2 file
		int hrtot = 0;
		int endrunon = 0;
		int prerunon = 0;
		int numt = 0;
		int bsize = 0;
		double Mcoil = 0;
		double SHR = 0;		
		double toutavg = 0;
		double tatticavg = 0;
		double thouseavg = 0;
		double achavg = 0;
		double gastherm = 0;
		double faneleckWh = 0;
		double compeleckWh = 0;
		double mveleckWh = 0;
		
		// "constants":
		double rhoref = 1.2;	// Reference air density at 20 deg C at sealevel
		int tref = 293;			// Reference room temp 20 deg C
		double pi = 3.1415926;

		double hret = 28;		// Initial number for hret in Btu/lb

		// RIVEC thermostat settings
		for(int i=0; i < 24; i++) {
			if(i < 7 || i > 22)							// Heating
				theat[i] = 273.15 + (65 - 32) * 5.0 / 9.0;
			else
				theat[i] = 273.15 + (68 - 32) * 5.0 / 9.0;
			
			if(i < 7 || i > 16)							// Cooling	
				tcool[i] = 273.15 + (78 - 32) * 5.0 / 9.0;
			else if(i > 6 && i < 13)
				tcool[i] = 273.15 + (83 - 32) * 5.0 / 9.0;
		}

		// Cooling
		tcool[13] = 273.15 + (82 - 32) * 5.0 / 9.0;
		tcool[14] = 273.15 + (81 - 32) * 5.0 / 9.0;
		tcool[15] = 273.15 + (80 - 32) * 5.0 / 9.0;
		tcool[16] = 273.15 + (79 - 32) * 5.0 / 9.0;

		// This reads in the filename and climate zone data from a file instead of input by user
		string f = Simfile[jen];
		string weatherfile = Climzone[jen];
		string Outfile = Flnme[jen];
		int dt = 1;							// Simulation timestep in minutes

		// ****** Opening moisture data file ******
		// OPEN "c:\resave\rcfiles\out\" + Outfile$ + ".hum" FOR OUTPUT AS #7
		_chdir("..\\output\\");
		ofstream moisturefile(Outfile + ".hum"); 
		if(!moisturefile) { 
			cout << "Cannot open: " << Outfile + ".hum" << endl;
			system("pause");
			return 1; 
		}
		_chdir("..\\input\\");
		ifstream buildingfile(f + ".csv"); 
		if(!buildingfile) { 
			cout << "Cannot open: " << f + ".csv" << endl;
			system("pause");
			return 1; 
		}

		//moisturefile << "testing 1\n" << "testing 2" << endl << "testing 3\n";
		
		buildingfile.getline(reading, 255);
		double C = atof(reading);		// Envelope Leakage Coefficient [m3/sPa^n]

		buildingfile.getline(reading, 255);
		double n = atof(reading);		// Envelope Pressure Exponent

		buildingfile.getline(reading, 255);
		double h = atof(reading);		// Eave Height [m]

		buildingfile.getline(reading, 255);
		double R = atof(reading);		// Ceiling Floor Leakage Sum

		buildingfile.getline(reading, 255);
		double X = atof(reading);		// Ceiling Floor Leakage Difference
		
		double lc = (R + X) * 50;		// Percentage of leakage in the ceiling
		double lf = (R - X) * 50;		// Percentage of leakage in the floor
		double lw = 100 - lf - lc;		// Percentage of leakage in the walls

		buildingfile.getline(reading, 255);
		int Nrflue = atoi(reading);
		
		for(int i=0; i < Nrflue; i++) {
			buildingfile.getline(reading, 255);
			Chim[i].Cflue = atof(reading);
			
			buildingfile.getline(reading, 255);
			Chim[i].Hflue = atof(reading);
			
			buildingfile.getline(reading, 255);
			Chim[i].Fluetemp = atof(reading);
		}
		
		// =========================== Leakage Inputs ==============================
		for(int i=0; i < 4; i++) {
			buildingfile.getline(reading, 255);
			Wallfraction[i] = atof(reading);

			buildingfile.getline(reading, 255);
			FloorFraction[i] = atof(reading);
		}

		buildingfile.getline(reading, 255);
		double Swrfl = atof(reading);

		buildingfile.getline(reading, 255);
		int Npipe = atoi(reading);
		
		for(int i=0; i < Npipe; i++) {
			buildingfile.getline(reading, 255);
			Pipe[i].wall = atoi(reading);

			buildingfile.getline(reading, 255);
			Pipe[i].h = atof(reading);

			buildingfile.getline(reading, 255);
			Pipe[i].A = atof(reading);

			buildingfile.getline(reading, 255);
			Pipe[i].n = atof(reading);

			buildingfile.getline(reading, 255);
			Pipe[i].Swf = atof(reading);

			buildingfile.getline(reading, 255);
			Pipe[i].Swoff = atof(reading);
		}

		// =========================== Building Inputs =============================
		buildingfile.getline(reading, 255);
		double Hfloor = atof(reading);

		buildingfile.getline(reading, 255);
		string row = reading;
		
		buildingfile.getline(reading, 255);
		double volume = atof(reading);

		buildingfile.getline(reading, 255);
		double FLOORAREA = atof(reading);		// Conditioned Floor Area (m2)
		
		buildingfile.getline(reading, 255);
		double planarea = atof(reading);		// Footprint of house (m2)
		
		buildingfile.getline(reading, 255);
		double Lhouse = atof(reading);			// Long side of house (m)
		
		buildingfile.getline(reading, 255);
		double Whouse = atof(reading);			// Short side of house (m)
		
		buildingfile.getline(reading, 255);
		double UAh = atof(reading);			// Heating U-Value (Thermal Conductance)
		
		buildingfile.getline(reading, 255);
		double UAc = atof(reading);			// Cooling U-Value (Thermal Conductance)

		double housevol = FLOORAREA * 2.5;

		// ====================== Venting Inputs (Windows/Doors)====================
		buildingfile.getline(reading, 255);
		int Nwindoor = atoi(reading);

		for(int i=0; i < Nwindoor; i++) {
			buildingfile.getline(reading, 255);
			Windoor[i].wall = atoi(reading);

			buildingfile.getline(reading, 255);
			Windoor[i].High = atof(reading);
			
			buildingfile.getline(reading, 255);
			Windoor[i].Wide = atof(reading);
			
			buildingfile.getline(reading, 255);
			Windoor[i].Top = atof(reading);
			
			buildingfile.getline(reading, 255);
			Windoor[i].Bottom = atof(reading);
		}

		// ================== Mechanical Venting Inputs (Fans etc) =================
		buildingfile.getline(reading, 255);
		int Nfan = atoi(reading);

		for(int i=0;  i < Nfan; i++) {
			buildingfile.getline(reading, 255);
			fan[i].power = atof(reading);

			buildingfile.getline(reading, 255);
			fan[i].q = atof(reading);

			buildingfile.getline(reading, 255);
			fan[i].oper = atof(reading);

			fan[i].on = 0;
		}

		buildingfile.getline(reading, 255);
		double windowWE = atof(reading);

		buildingfile.getline(reading, 255);
		double windowN = atof(reading);

		buildingfile.getline(reading, 255);
		double windowS = atof(reading);
		
		buildingfile.getline(reading, 255);
		double shadcoef = atof(reading);

		buildingfile.getline(reading, 255);
		double rceilh = atof(reading);

		buildingfile.getline(reading, 255);
		double rceilc = atof(reading);

		buildingfile.getline(reading, 255);
		double latload = atof(reading);

		buildingfile.getline(reading, 255);
		double intgain1 = atof(reading);

		// =========================== Attic Inputs ================================
		buildingfile.getline(reading, 255);
		double ATTICVOL = atof(reading);

		buildingfile.getline(reading, 255);
		double Cattic = atof(reading);

		buildingfile.getline(reading, 255);
		double Nattic = atof(reading);

		for(int i=0; i < 5; i++) {
			buildingfile.getline(reading, 255);
			Soffitfraction[i] = atof(reading);			
		}
		for(int i=0; i < 4; i++) {
			buildingfile.getline(reading, 255);
			Soffit[i].h = atof(reading);			
		}

		// =========================== Attic Vent Inputs ===========================
		buildingfile.getline(reading, 255);
		int Navent = atoi(reading);

		for(int i=0; i < Navent; i++) {
			buildingfile.getline(reading, 255);
			Avent[i].wall = atoi(reading);

			buildingfile.getline(reading, 255);
			Avent[i].h = atof(reading);

			buildingfile.getline(reading, 255);
			Avent[i].A = atof(reading);
			
			buildingfile.getline(reading, 255);
			Avent[i].n = atof(reading);
		}

		// =========================== Roof Inputs =================================
		buildingfile.getline(reading, 255);
		double pitch = atof(reading);

		buildingfile.getline(reading, 255);
		string RR = reading;

		buildingfile.getline(reading, 255);
		double Hpeak = atof(reading);
		
		buildingfile.getline(reading, 255);
		int Nafans = atoi(reading);

		for(int i = 0; i < Nafans; i++) {
			buildingfile.getline(reading, 255);
			AFan[i].power = atof(reading);

			buildingfile.getline(reading, 255);
			AFan[i].q = atof(reading);
			
			buildingfile.getline(reading, 255);
			AFan[i].oper = atof(reading);
		}

		buildingfile.getline(reading, 255);
		int rroof = atoi(reading);

		buildingfile.getline(reading, 255);
		int rooftype = atoi(reading);

		// =========================== Duct Inputs =================================
		buildingfile.getline(reading, 255);
		double location = atof(reading);

		buildingfile.getline(reading, 255);
		double supthick = atof(reading);
		
		buildingfile.getline(reading, 255);
		double retthick = atof(reading);

		buildingfile.getline(reading, 255);
		double supRval = atof(reading);

		buildingfile.getline(reading, 255);
		double retRval = atof(reading);

		buildingfile.getline(reading, 255);
		double supLF = atof(reading);

		buildingfile.getline(reading, 255);
		double retLF = atof(reading);

		buildingfile.getline(reading, 255);
		double suplength = atof(reading);

		buildingfile.getline(reading, 255);
		double retlength = atof(reading);

		buildingfile.getline(reading, 255);
		double supdiam = atof(reading);

		buildingfile.getline(reading, 255);
		double retdiam = atof(reading);

		buildingfile.getline(reading, 255);
		double cqah = atof(reading);				// Cooling air flow
		
		buildingfile.getline(reading, 255);
		double hqah = atof(reading);				// Heating air flow
		
		buildingfile.getline(reading, 255);
		double supn = atof(reading);
		
		buildingfile.getline(reading, 255);
		double retn = atof(reading);
		
		buildingfile.getline(reading, 255);
		double supC = atof(reading);
		
		buildingfile.getline(reading, 255);
		double retC = atof(reading);
		
		// NOTE: The burried variable is not used but left in code to continue proper file navigation
		buildingfile.getline(reading, 255);
		double burried = atof(reading);
		
		// =========================== Equipment Inputs ============================

		buildingfile.getline(reading, 255);
		double capacityraw = atof(reading);
		
		buildingfile.getline(reading, 255);
		double capacityari = atof(reading);
		
		buildingfile.getline(reading, 255);
		double EERari = atof(reading);
		
		buildingfile.getline(reading, 255);
		double hcapacity = atof(reading);		// Heating capacity kBtu/hour NOW +VE
		
		buildingfile.getline(reading, 255);
		double hfanpower = atof(reading);		// Heating fan power
		
		buildingfile.getline(reading, 255);
		double cfanpower = atof(reading);		// Cooling fan power
		
		buildingfile.getline(reading, 255);
		double charge = atof(reading);
		
		buildingfile.getline(reading, 255);
		double afue = atof(reading);			// AFUE efficiency
		
		buildingfile.getline(reading, 255);
		int bathroomSchedule = atoi(reading);			// Number of bathrooms (basically number of seaprate bathroom fans)
		
		buildingfile.getline(reading, 255);
		int Nbedrooms = atoi(reading);			// Number of occupants (for 62.2 target ventilation calculation)

		buildingfile.close();

		// In case the user enters leakage fraction in % rather than as a fraction
		if(supLF > 1)
			supLF = supLF / 100;
		if(retLF > 1)
			retLF = retLF / 100;
		if(cqah > 100) {				// This assumes that we will not have fan flows greater than 100 m^3/s
			cqah = cqah * .0004719;		// CFM to m^3/s
			hqah = hqah * .0004719;		// CFM to m^3/s
		}

		// m^3/s to cfm
		double qahcfm = cqah / .0004719;				// need air handler flow in cfm for equipment calculations
		// The following is the cooling capacity incorrect air flow correction term
		double qahcorrec = 1.62 - .62 * qahcfm / (400 * capacityraw) + .647 * log(qahcfm / (400 * capacityraw));
		double supcp = 753.624;							// j/kg/K
		double suprho = 16.018 * 2;						// kg/m^3 the factor of two represents the plastic and sprical
		double retcp = 753.624;							// j/kg/K
		double retrho = 16.018 * 2;						// kg/m^3
		double suparea = (supdiam + 2 * supthick) * pi * suplength;
		double retarea = (retdiam + 2 * retthick) * pi * retlength;
		double supvol = (pow(supdiam, 2) * pi / 4) * suplength;
		double retvol = (pow(retdiam, 2) * pi / 4) * retlength;
		hcapacity = hcapacity * .2913 * 1000 * afue;     // AFUE efficiency variable inserted (afue)
		double MWha = .5 * FLOORAREA / 186;

		// MWha is the moisture transport coefficient that scales with floor area (to scale with surface area of moisture
		// Source/sink, it is set to 0.05 for a 2000ft2 house (186 m2)
		
		//---------dkm: new inputs-----------------------------------------------------------------
		int dcv = 0;					// Dose Controlled Ventilation flag. Indicates whether there is a dose/exp controlled fan - changes later by itself
		int qflue = 0;					// NOTE: This variable isnt being used...
		
		//------- Check if Dynamic Fan and Occupancy scheduling should be used ----------------------------------
		int occupancyFlag = 0;
		int dynamicScheduleFlag = 0;

		for(int i=0; i < Nfan; i++) {
			if(fan[i].oper == 23 || fan[i].oper == 24 || fan[i].oper == 25 || fan[i].oper == 26 || fan[i].oper == 27)
				dynamicScheduleFlag = 1;
			if(fan[i].oper == 50)
				occupancyFlag = 1;		// Occupancy flag to determine if occupancy scheduling is used (0 = off, 1 = on)
		}

		// ================ INITIALIZE SIMULATION VARIABLES =======================================
		// Always use urban shelter from data file
		row = "R"; // RAD :  this variable should be read from input file and not modified here, remove after debugging.
		
		// Reading shelter values from a file.  Bshelter.dat contains shelter values
		// Computed for the houses at AHHRF for every degree of wind angle
		
		ifstream shelterfile(shelterfile_name); 
		if(!shelterfile) { 
			cout << "Cannot open: " << shelterfile_name << endl;
			system("pause");
			return 1; 
		}

		// FF: This angle variable is being overwritten to read Swinit in proper location per FOR iteration. Not used in code.
		double angle;
		for(int i=0; i < 361; i++) {
			shelterfile >> angle >> Swinit[0][i] >> Swinit[1][i] >> Swinit[2][i] >> Swinit[3][i];
		}

		shelterfile.close();

		// The speed correction allows wind speeds at meteorological stations to be converted to eaves height wind speed
		double speedcorr = pow((80. / 10), .15) * pow((h / 80), .3);
		double Cceil = C * (R + X) / 2;
		
		// AL4 is used to estimate flow velocities in the attic
		double AL4 = Cattic * sqrt(rhoref / 2) * pow(4, (Nattic - .5));
		
		// New char velocity for unvented attics
		if(AL4 == 0)
			AL4 = Cceil * sqrt(rhoref / 2) * pow(4, (Nattic - .5));
		
		// AL5 is use to estimate flow velocities in the house
		double AL5 = C * sqrt(rhoref / 2) * pow(4, (n - .5));

		// ================= CREATE OUTPUT FILE =================================================
		// OPEN "c:\resave\rcfiles\out\" + Outfile$ + ".out" FOR OUTPUT AS #2
		_chdir("..\\output\\");
		ofstream outputfile(Outfile + ".out"); 
		if(!outputfile) { 
			cout << "Cannot open: " << Outfile + ".out" << endl;
			system("pause");
			return 1; 
		}

		outputfile << "* Input file = " << f << " Weather file = " << weatherfile << endl;
		outputfile << "Time,minutes,tout,tatt,thouse,tsup,tret,ahon,ahpower,hcap, compresspower, mechventpower, hr,qhouse,ach,ccap,SHR,Mcoil,setpoint,housepress,turnover, relexp, reldose, fan1, fan2, fan3, fan4, fan5, fan6, ventsum, rivecOn, nonRivecVentSum, achflue, occupied" << endl;

		// ================== OPEN WEATHER FILE FOR INPUT ========================================
		_chdir("..\\input\\");
		ifstream ws4file(weatherfile + ".ws4"); 
		if(!ws4file) { 
			cout << "Cannot open: " << weatherfile + ".ws4" << endl;
			system("pause");
			return 1; 
		}

		double latitude, altitude;
		ws4file >> latitude >> altitude;

		// set 1 = air handler off, and house air temp below setpoint minus 0.5 degrees
		// set 0 = air handler off, but house air temp above setpoint minus 0.5 degrees
		int set = 0;

		int ahflag = 0;			// Air Handler Flag (0/1 = OFF/ON)
		int ahflagprev = 0;
		int econoflag = 0;		// Economizer Flag (0/1 = OFF/ON)

		// ---------------- Attic and duct system heat transfer nodes --------------------
		for(int k=0; k < 16; k++) {
			told[k] = 293.15;
		}
		told[2] = 278;
		told[4] = 278;

		double tpredattic = told[0];
		double tpredret = told[11];
		double tpredsup = told[14];
		double Tpredhouse = told[15];

		// Initialise time keeping variables
		int daymin = 0;				// Minute number in the current day (1 to 1440)
		int hourmin = 0;			// Minute number in the current hour (1 to 60)

		// Setting initial values of air mass for moisture balance:
		double M1 = ATTICVOL * 1.2;			// Mass of attic air
		double M12 = retvol * 1.2;			// Mass of return air
		double M15 = supvol * 1.2;			// Mass of supply air
		double M16 = housevol * 1.2;		// Mass of house air
		double Mw5 = 60 * FLOORAREA;		// ?

		// Like the mass transport coefficient, the active mass for moisture scales with floor area
		// The following are defined for RIVEC
		double aeq = .001 * (.05 * FLOORAREA + 3.5 * (Nbedrooms + 1)) * 3600 / volume;			// 62.2 target ventilation rate (ACH)
		double reldose = 1;						// Initial value for relative dose
		double turnover = 1 / aeq;					// Initial value for turnover time (hrs)
		double dtau = dt * 60;					// One minute timestep (in seconds)
		double rivecdt = dtau / 3600;				// Rivec timestep is in hours, dtau is simulation timestep in seconds

		// ======================= NEW FAN SCHEDULE INPUTS: WJNT ========================================
		// Read in fan schedule (lists of 1s and 0s, 1 = fan ON, 0 = fan OFF, for every minute of the year)
		// Different schedule file depending on number of bathrooms
		string fanSchedule = "no_fan_schedule";
		ifstream fanschedulefile;
		
		if(dynamicScheduleFlag == 1) {
			switch (bathroomSchedule) {
			case 1:
				fanSchedule = fanSchedulefile_name1;			// One bathroom fan using a dynamic schedule
				break;
			case 2:
				fanSchedule = fanSchedulefile_name2;			// Two bathroom fans using dynamic schedules
				break;
			case 3:
				fanSchedule = fanSchedulefile_name3;			// Three bathroom fans using dynamic schedules
				break;
			}

			fanschedulefile.open(fanSchedule + ".csv"); 
			if(!fanschedulefile) { 
				cout << "Cannot open: " << fanSchedule + ".csv" << endl;
				system("pause");
				return 1; 
			}			
		}
		
		// =========================== END OF NEW INPUTS: WJNT ==========================================
		// ============== THE SIMULATION LOOP FOR MINUTE-BY-MINUTE STARTS HERE: =========================

		time_t starttime, endtime;
		time(&starttime);
		string ts = ctime(&starttime);			// Start wall clock time for simulation
		
		int ahminutes;
		int target;
		int day;
		int idirect;
		int solth;
		int dryerFan = 0;
		int kitchenFan = 0;
		int bathOneFan = 0;
		int bathTwoFan = 0;
		int bathThreeFan = 0;
		int HOUR;
		int ttime;
		int weekendflag;
		int occupied[24] = {0};
		int comptime = 0;
		int comptimecount = 0;
		int rivecOn;
		int itercount;
		int Crawl = 0;
		int ERRCODE = 0;
		int prevday = 0;
		
		double pref;
		double sc;
		double T;
		double Tout;
		double HROUT;
		double speed;
		double direction;		
		double mfancyc;
		double hcap;
		double mhrv;
		double fanheat;
		double ventsumin;
		double ventsumout;
		double turnoverold;
		double reldoseold;
		double nonRivecVentSumIn;
		double nonRivecVentSumOut;
		double nonRivecVentSum;
		double qah;
		double UA;
		double rceil;
		double econodt = 0;
		double SUPVELah;
		double retvelah;
		double QsupREG;
		double QretREG;
		double qretleak;
		double qsupleak;
		double MSUPLEAK1;
		double MRETLEAK1;
		double MSUPREG1;
		double mretreg1;
		double MAH1;
		double MSUPREG;
		double MAH;
		double mretleak;
		double MSUPLEAK;
		double mretreg;
		double SUPVEL;
		double retvel;
		double fanpow;
		double fanpower;
		double msupahoff=0;
		double mretahoff=0;
		double relexp=0;
		double evapcap = 0;
		double latcap = 0;
		double capacity = 0;
		double capacityh = 0;
		double powercon = 0;
		double powercomp = 0;
		double capacityc = 0;
		double Toutf = 0;
		double tretf = 0;
		double dhret = 0;
		double chargecapd = 0;
		double chargeeerd = 0;
		double EER = 0;
		double chargecapw = 0;
		double chargeeerw = 0;
		double Mcoilprevious = 0;
		double Mceilold;
		double limit;
		double matticenvout = 0;
		double Mceil = 0;
		//double mretahaoff;
		double matticenvin = 0;
		double mhousein = 0;
		double mhouseout = 0;
		double RHOATTIC;
		double intgain;
		double Tpredceil;
		double TPREDFLOOR;
		double TPREDSI;
		double TPREDSO;
		double TPREDNI;
		double TPREDNO;
		double TPREDENDI;
		double TPREDENDO;
		double TPREDWOOD;
		double TPREDRETSURF;
		double TPREDSUPSURF;
		double tpredhousemass;
		double flag;
		double min = 0;
		double mout = 0;
		double Pint = 0;
		double Mflue = 0;
		double dPflue = 0;
		double dPceil = 0;
		double dPfloor = 0;
		double Patticint = 0;
		double Matticin = 0;
		double Matticout = 0;
		double TSKY = 0;
		double w = 0;
		double solgain = 0;
		double mhouse = 0;
		double qhouse = 0;
		double achhouse = 0;
		double achflue = 0;
		double ventsum = 0;
		double ventsumD = 0;
		double ventsumoutD = 0;
		double ventsuminD = 0;
		double merv = 0;


		do {
			if(daymin == 1440)					// Resetting number of minutes into day every 24 hours
				daymin = 0;
			
			if(hourmin == 60)
				hourmin = 0;					// Resetting number of minutes into hour

			hourmin++;							// Minute count (of hour)

			if(hourmin == 1) {
				ahminutes = 0;					// Resetting air handler operation minutes for this hour
				mfancyc = 0;					// Fan cycler?
			}
			
			target = hourmin - 40;				// Target is used for fan cycler operation currently set for 20 minutes operation

			if(target < 0)
				target = 0;						// For other time periods, replace the "40" with 60 - operating minutes
			
			hcap = 0;
			mvpow = 0;
			mhrv = 0;
			fanheat = 0;
			econoflag = 0;						// Setting economizer to "off"
			ventsumin = 0;						// Setting sum of supply mechanical ventilation to zero
			ventsumout = 0;						// Setting sum of exhaust mechanical ventilation to zero
			turnoverold = turnover;
			reldoseold = reldose;
			nonRivecVentSumIn = 0;				// Setting sum of non-RIVEC supply mechanical ventilation to zero
			nonRivecVentSumOut = 0;				// Setting sum of non-RIVEC exhaust mechanical ventilation to zero
			
			//********************* get weather data ***************************
			ws4file >> day >> idirect >> solth >> T >> HROUT >> speed >> direction >> pref >> sc;

			// Schedule Inputs: WJNT Apr, 2011
			// Assumes operation of dryer and kitchen fans, then 1 - 3 bathroom fans
			if(dynamicScheduleFlag == 1) {
				switch (bathroomSchedule) {
				case 1:
					fanschedulefile >> dryerFan >> kitchenFan >> bathOneFan;
					break;
				case 2:
					fanschedulefile >> dryerFan >> kitchenFan >> bathOneFan >> bathTwoFan;
					break;
				case 3:
					fanschedulefile >> dryerFan >> kitchenFan >> bathOneFan >> bathTwoFan >> bathThreeFan;
					break;
				}
			}

			sc = sc / 10;					// Converting to decimal fraction
			speed = speed * speedcorr;		// Correct wind speed to building eaves height [m/s]
			
			if(speed < 1)					// Minimum wind velocity allowed is 1 m/s
				speed = 1;

			for(int k=0; k < 4; k++)		// Wind direction as a compass direction?
				Sw[k] = Swinit[k][int (direction + 1) - 1];		// -1 in order to allocate from array (1 to n) to array (0 to n-1)
			
			pref = 1000 * pref;				// Reference Pressure [Pa]
			Tout = 273 + T;					// Outside Air Temperature [K]
			HOUR = int (daymin / 60);		// HOUR is hours into the current day - used to control diurnal cycles for fans
			
			if(HOUR == 24)
				HOUR = 0;
			
			hrtot = hrtot + dt;				// Is number of minutes into simulation - used for beginning and end of cycle fan operation
			ttime = int ((daymin - 1) / 60) + 1;      // Hours into day for thermostat schedule

			if((int (day / 7) == day / 7.) || (int ((day - 1) / 7) == (day - 1) / 7.))
				weekendflag = 1;
			else
				weekendflag = 0;
			
			// if(day == 1) weekendflag = 0;        // Not sure why???
			// ======================= OCCUPANCY SCHEDULING =======================

			if(weekendflag == 1 && occupancyFlag == 1) {		
					for(int i=0; i < 24; i++)
						occupied[i]= 1;
			}
			if(weekendflag == 0 && occupancyFlag == 1) {		
					for(int i=0; i < 24; i++) {
						if(i < 8 || i > 15)
							occupied[i]= 1;
						else
							occupied[i]= 0;
					}
			}

			if(hrtot == 1) {				// Setting initial humidity conditions
				for(int i = 0; i < 5; i++) {
					hrold[i] = HROUT;
					HR[i] = HROUT;
				}
			}

			double rhoout = rhoref * tref / Tout;			// Outside Air Density
			double rhoin = rhoref * tref / Tpredhouse;		// Inside Air Density
			double RHOSUP = rhoref * tref / tpredsup;		// Supply Duct Air Density
			double RHORET = rhoref * tref / tpredret;		// Return Duct Air Density
			
			// Print to screen stage of simulation
			
			// FF: if you want to print a message every house use this code:
			//		cout << "Day = " << day << " HR = " << HOUR << endl;

			// FF: if you want to print a message every day use this code:
			if(prevday != day) {
				cout << "Day = " << day << endl;
				prevday = day;
			}
			
			// Finding the month for solar radiation calculations
			int Month;

			if(day <= 31)
				Month = 1;
			if(day > 31 && day <= 59)
				Month = 2;
			if(day > 60 && day <= 90)
				Month = 3;
			if(day > 91 && day <= 120)
				Month = 4;
			if(day > 121 && day <= 151)
				Month = 5;
			if(day > 152 && day <= 181)
				Month = 6;
			if(day > 182 && day <= 212)
				Month = 7;
			if(day > 213 && day <= 243)
				Month = 8;
			if(day > 244 && day <= 273)
				Month = 9;
			if(day > 274 && day <= 304)
				Month = 10;
			if(day > 305 && day <= 334)
				Month = 11;
			if(day > 335)
				Month = 12;

			int NOONMIN = 720 - daymin;
			double HA = pi * .25 * NOONMIN / 180;			// Hour Angle
			double L = pi * latitude / 180;					// LATITUDE
			
			// SOLAR DECLINATION FROM ASHRAE P.27.2, 27.9IP  BASED ON 21ST OF EACH MONTH
			double dec;
			double Csol;

			switch (Month) {
			case 1:
				dec = -20 * pi / 180;
				Csol = .103;
				break;
			case 2:
				dec = -10.8 * pi / 180;
				Csol = .104;
				break;
			case 3:
				dec = 0;
				Csol = .109;
				break;
			case 4:
				dec = 11.6 * pi / 180;
				Csol = .12;
				break;
			case 5:
				dec = 20 * pi / 180;
				Csol = .13;
				break;
			case 6:
				dec = 23.45 * pi / 180;
				Csol = .137;
				break;
			case 7:
				dec = 20.6 * pi / 180;
				Csol = .138;
				break;
			case 8:
				dec = 12.3 * pi / 180;
				Csol = .134;
				break;
			case 9:
				dec = 0;
				Csol = .121;
				break;
			case 10:
				dec = -10.5 * pi / 180;
				Csol = .111;
				break;
			case 11:
				dec = -19.8 * pi / 180;
				Csol = .106;
				break;
			case 12:
				dec = -23.45 * pi / 180;
				Csol = .103;
				break;
			}

			double SBETA = cos(L) * cos(dec) * cos(HA) + sin(L) * sin(dec);
			double CBETA = sqrt(1 - pow(SBETA, 2));
			double CTHETA = CBETA;
			
			// accounting for slight difference in solar measured data from approximate geometry calcualtions
			if(SBETA < 0)
				SBETA = 0;
			
			double hdirect = SBETA * idirect;
			double diffuse = solth - hdirect;

			// FOR SOUTH ROOF
			double SIGMA = pi * pitch / 180;		//ROOF PITCH ANGLE
			CTHETA = CBETA * 1 * sin(SIGMA) + SBETA * cos(SIGMA);
			double ssolrad = idirect * CTHETA + diffuse;
			
			// FOR NORTH ROOF
			CTHETA = CBETA * -1 * sin(SIGMA) + SBETA * cos(SIGMA);
			
			if(CTHETA < 0)
				CTHETA = 0;
			
			double nsolrad = idirect * CTHETA + diffuse;
			int hcflag;

			if(day == 1)
				hcflag = 1;
			
			// Setting thermostat heat/cool mode:
			if(told[15] < 295.35)			// 295.35K = 22.2C = 71.96F
				hcflag = 1;                 // Turn on HEATING
			else
				hcflag = 2;                 // Turn on COOLING
			
			// ====================== HEATING THERMOSTAT CALCULATIONS ============================
			if(hcflag == 1) {
				qah = hqah;              // Heating Air Flow Rate [m3/s]
				UA = UAh;                // Thermal Conductance [W/K]
				rceil = rceilh;

				if(told[15] > (theat[ttime-1] + .5))
					ahflag = 0;					// Heat off if building air temp above setpoint

				if(ahflag == 0 || ahflag == 100) {
					if(told[15] <= (theat[ttime-1] - .5))
						set = 1;				// Air handler off, and tin below setpoint - 0.5 degrees
					else
						set = 0;				// Air handler off, but tin above setpoint - 0.5 degrees
				}

				if(told[15] < (theat[ttime-1] - .5))
					ahflag = 1;					// turn on air handler/heat

				if(told[15] >= (theat[ttime-1] - .5) && told[15] <= (theat[ttime-1] + .5)) {
					if(set == 1)
						ahflag = 1;
					else
						ahflag = 0;
				}

				if(ahflag == 0 && ahflag != ahflagprev)
					endrunon = hrtot + 1;		// Adding 1 minute runon during which heat from beginning of cycle is put into air stream

				if(ahflag == 1)
					hcap = hcapacity / afue;	// hcap is gas consumed so needs to decide output (hcapacity) by AFUE

				// ====================== COOLING THERMOSTAT CALCULATIONS ========================
			} else {
				hcap = 0;
				qah = cqah;						// Cooling Air Flow Rate
				UA = UAc;
				rceil = rceilc;
				endrunon = 0;
				prerunon = 0;

				if(told[15] < (tcool[ttime-1] - .5))
					ahflag = 0;

				if(ahflag == 0 || ahflag == 100) {
					if(told[15] >= (tcool[ttime-1] + .5))
						set = 1;				// 1 = AH OFF and house temp below thermostat setpoint - 0.5
					else
						set = 0;				// 0 = AH OFF and house temp above thermostat setpoint - 0.5
				}

				if(told[15] >= (tcool[ttime-1] + .5))
					ahflag = 2;

				// if(told[15] >= (tcool[ttime-1] + .5 + tcooladder)) ahflag = 2;  // IN ARTI CODE LINE 645
				if(told[15] < (tcool[ttime-1] + .5) && told[15] >= (tcool[ttime-1] - .5)) {
					// if(told[15] < (tcool[ttime-1] + .5 + tcooladder) && told[15] >= (tcool[ttime-1] - .5 + tcooladder)) { //IN ARTI CODE LINE 646
					if(set == 1)
						ahflag = 2;
					else
						ahflag = 0;
				}

				if(ahflag != ahflagprev && ahflag == 2) {				// First minute of operation
					comptime = 1;
					comptimecount = 1;
				} else if(ahflag == 2 && comptimecount == 1) {			// Second minute of operation
					comptime = 2;
					comptimecount = 0;
				} else if(ahflag == 2 && comptimecount == 0) {
					comptime = 0;
				} else if(ahflag != ahflagprev && ahflag == 0) {		// End of cooling cycle
					comptime = 0;
				}

				econodt = (Tpredhouse - Tout) * 9.0 / 5.0;				// converting indoor-outdoor temp to F for economizer switching

				if(ahflag == 0 && econodt >= 6 && Tpredhouse > 294.15)	// no cooling and 6F dt for economizer operation
					econoflag = 1;
				else
					econoflag = 0;
			}

			ahflagprev = ahflag;
			if(ahflag != 0)
				ahminutes++;			// counting number of air handler operation minutes in the hour

			// the following air flows depend on if we are heating or cooling
			SUPVELah = qah / (pow(supdiam,2) * pi / 4);
			retvelah = qah / (pow(retdiam,2) * pi / 4);
			QsupREG = -qah * supLF + qah;
			QretREG = qah * retLF - qah;
			qretleak = -qah * retLF;
			qsupleak = qah * supLF;
			MSUPLEAK1 = qsupleak * rhoin;
			MRETLEAK1 = qretleak * rhoin;
			MSUPREG1 = QsupREG * rhoin;
			mretreg1 = QretREG * rhoin;
			MAH1 = qah * rhoin;
			mfancyc = 0;
			mvpow = 0;					// Total vent fan power consumption
			
			if(ahflag == 0) {			// AH OFF
				// the 0 at the end of these terms means that the AH is off
				if(hrtot >= endrunon) {  // endrunon is the end of the heating cycle + 1 minutes
					MSUPREG = 0;
					MAH = 0;
					mretleak = 0;
					MSUPLEAK = 0;
					mretreg = 0;
					SUPVEL = abs(msupahoff) / RHOSUP / (pow(supdiam,2) * 3.1415926 / 4);
					retvel = abs(mretahoff) / RHORET / (pow(retdiam,2) * 3.1415926 / 4);
					fanpow = 0;			// Fan power consumption [W]
					fanpower = 0;		// Fan heat into air stream [W]
					hcap = 0;			// Gas burned by furnace NOT heat output (that is hcapacity)
				} else {
					ahflag = 102;		// Note: need to have capacity changed to one quarter of burner capacity during cool down
					MAH = MAH1;
					mretleak = MRETLEAK1;
					MSUPLEAK = MSUPLEAK1;
					mretreg = mretreg1;
					MSUPREG = MSUPREG1;
					SUPVEL = SUPVELah;
					retvel = retvelah;
					fanpower = hfanpower * .85;			// Heating fan power multiplied by an efficiency
					fanpow = hfanpower;
					hcap = hcapacity / afue;

					for(int i = 0; i < Nfan; i++) {     // For outside air into return cycle operation
						if(fan[i].oper == 6 && ahminutes <= 20) {
							mfancyc = fan[i].q * rhoin;
							mretreg = mretreg1 - mfancyc;
						}
						if(fan[i].oper == 10) {			// vent always open during any air handler operation
							mfancyc = fan[i].q * rhoin;
							mretreg = mretreg1 - mfancyc;
						}
						if(fan[i].oper == 14) {			// vent always open during any air handler operation
							mfancyc = fan[i].q * rhoin;
							mretreg = mretreg1 - mfancyc;
						}
						if(fan[i].oper == 18) {			// vent always open during any air handler operation
							mfancyc = fan[i].q * rhoin;
							mretreg = mretreg1 - mfancyc;
						}
						if(fan[i].oper == 15 && ahminutes <= 20) {
							mfancyc = fan[i].q * rhoin;
							mretreg = mretreg1 - mfancyc;
						}
						if(fan[i].oper == 13 && ahminutes <= 20) {
							mfancyc = fan[i].q * rhoin;
							mretreg = mretreg1 - mfancyc;
						}
						if(fan[i].oper == 22 && ahminutes <= 20 && econoflag == 0) {
							mfancyc = fan[i].q * rhoin;
							mretreg = mretreg1 - mfancyc;
						}
					}
				}
			} else {
				if(ahflag == 2) {	//	cooling 
					hcap = 0;
					MAH = MAH1;
					mretleak = MRETLEAK1;
					MSUPLEAK = MSUPLEAK1;
					mretreg = mretreg1;
					MSUPREG = MSUPREG1;
					SUPVEL = SUPVELah;
					retvel = retvelah;
					fanpower = cfanpower * .85;				// Cooling fan power multiplied by an efficiency
					fanpow = cfanpower;
					for(int i = 0; i < Nfan; i++) {			// for outside air into return cycle operation
						if(fan[i].oper == 6 && ahminutes <= 20) {
							mfancyc = fan[i].q * rhoin * cqah / hqah;	// correcting this return intake flow so it remains a constant fraction of air handler flow
							mretreg = mretreg1 - mfancyc;
						}
						if(fan[i].oper == 10) {				// vent always open during any air handler operation
							mfancyc = fan[i].q * rhoin;
							mretreg = mretreg1 - mfancyc;
						}
						if(fan[i].oper == 14) {				// vent always open during any air handler operation
							mfancyc = fan[i].q * rhoin;
							mretreg = mretreg1 - mfancyc;
						}
						if(fan[i].oper == 18) {				// vent always open during any air handler operation
							mfancyc = fan[i].q * rhoin;
							mretreg = mretreg1 - mfancyc;
						}
						if(fan[i].oper == 15 && ahminutes <= 20) {
							mfancyc = fan[i].q * rhoin;
							mretreg = mretreg1 - mfancyc;
						}
						if(fan[i].oper == 13 && ahminutes <= 20) {
							mfancyc = fan[i].q * rhoin;
							mretreg = mretreg1 - mfancyc;
						}
						if(fan[i].oper == 22 && ahminutes <= 20 && econoflag == 0) {
							mfancyc = fan[i].q * rhoin;
							mretreg = mretreg1 - mfancyc;
						}
					}
				} else {
					MAH = MAH1;
					mretleak = MRETLEAK1;
					MSUPLEAK = MSUPLEAK1;
					mretreg = mretreg1;
					MSUPREG = MSUPREG1;
					SUPVEL = SUPVELah;
					retvel = retvelah;
					fanpower = hfanpower * .85;
					fanpow = hfanpower;
					hcap = hcapacity / afue;

					for(int i = 0; i < Nfan; i++) {			// for outside air into return cycle operation
						if(fan[i].oper == 6 && ahminutes <= 20) {
							mfancyc = fan[i].q * rhoin;
							mretreg = mretreg1 - mfancyc;
						}
						if(fan[i].oper == 10) {				// vent always open during any air handler operation
							mfancyc = fan[i].q * rhoin;
							mretreg = mretreg1 - mfancyc;
						}
						if(fan[i].oper == 14) {				// vent always open during any air handler operation
							mfancyc = fan[i].q * rhoin;
							mretreg = mretreg1 - mfancyc;
						}
						if(fan[i].oper == 18) {				// vent always open during any air handler operation
							mfancyc = fan[i].q * rhoin;
							mretreg = mretreg1 - mfancyc;
						}
						if(fan[i].oper == 15 && ahminutes <= 20) {
							mfancyc = fan[i].q * rhoin;
							mretreg = mretreg1 - mfancyc;
						}
						if(fan[i].oper == 13 && ahminutes <= 20) {
							mfancyc = fan[i].q * rhoin;
							mretreg = mretreg1 - mfancyc;
						}
						if(fan[i].oper == 22 && ahminutes <= 20 && econoflag == 0) {
							mfancyc = fan[i].q * rhoin;
							mretreg = mretreg1 - mfancyc;
						}
					}
				}
			}

			for(int i=0; i < Nfan; i++) {
				// ***********************  like fan.oper = 10  but with a 20 minute minimum invoked here:
				if(fan[i].oper == 13) {			// fan cycler operation without exhaust fan
					if(ahminutes < target && ahflag != 1 && ahflag != 2 && ahflag != 102) {		// conditions for turning on air handler only for fan cycler operation
						ahflag = 100;
						ahminutes = ahminutes + 1;
						qah = cqah;
						SUPVELah = qah / (pow(supdiam,2) * pi / 4);
						retvelah = qah / (pow(retdiam,2) * pi / 4);
						QsupREG = -qah * supLF + qah;
						QretREG = qah * retLF - qah;
						qretleak = -qah * retLF;
						qsupleak = qah * supLF;
						MSUPLEAK = qsupleak * rhoin;
						mretleak = qretleak * rhoin;
						// *** fan cycler air is from outside so not simply added to qretleak or mretleak
						// *** but needs to be subtracted from mretreg1

						// FF: Looked at this comment but not sure if outside airflow was already added or not?
						// *************** need to add outside air flow to mretleak Fan[i].Q * rhoin if Fan[i].oper = 1 i.e. continuous exhaust supplemented by occasions "supply"
						//   Fan[i].Q will be an exhaust fan and therfore a negative number
						mfancyc = fan[i].q * rhoin;
						MSUPREG = QsupREG * rhoin;
						mretreg = QretREG * rhoin - mfancyc;
						MAH = qah * rhoin;
						mvpow = mvpow + cfanpower;			// note ah operates at cooling speed
						fanpower = cfanpower * .85;			// for heat
						ventsumin = ventsumin + abs(fan[i].q) * 3600 / volume;
						nonRivecVentSumIn = nonRivecVentSumIn + abs(fan[i].q) * 3600 / volume;
					}
				}

				if(fan[i].oper == 22) {						// Fan cycler operation only when economizer is off
					if(econoflag == 0) {					// Only open CFI if economizer is off
						if(ahminutes < target && ahflag != 1 && ahflag != 2 && ahflag != 102) {		// conditions for turning on air handler only for fan cycler operation
							ahflag = 100;
							ahminutes = ahminutes + 1;
							qah = cqah;
							SUPVELah = qah / (pow(supdiam,2) * pi / 4);
							retvelah = qah / (pow(retdiam,2) * pi / 4);
							QsupREG = -qah * supLF + qah;
							QretREG = qah * retLF - qah;
							qretleak = -qah * retLF;
							qsupleak = qah * supLF;
							MSUPLEAK = qsupleak * rhoin;
							mretleak = qretleak * rhoin;
							// *** fan cycler air is from outside so not simply added to qretleak or mretleak
							// *** but needs to be subtracted from mretreg1
							// *************** need to add outside air flow to mretleak Fan(i).Q * rhoin if Fan(i).oper = 1 i.e. continuous exhaust supplemented by occasions "supply"
							// Fan(i).Q will be an exhaust fan and therfore a negative number
							mfancyc = fan[i].q * rhoin;
							MSUPREG = QsupREG * rhoin;
							mretreg = QretREG * rhoin - mfancyc;
							MAH = qah * rhoin;
							mvpow = mvpow + cfanpower;    // note ah operates at cooling speed
							fanpower = cfanpower * .85;  // for heat
							ventsumin = ventsumin + abs(fan[i].q) * 3600 / volume;
							nonRivecVentSumIn = nonRivecVentSumIn + abs(fan[i].q) * 3600 / volume;
						}
					}
				}
			}

			// =============================== FAN CONTROLS ==============================================

			// *************************************** [RECHECK - START] **********************************************************>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

			for(int i = 0; i < Nfan; i++) {
				if(fan[i].oper == 1) {							// FAN ALWAYS ON
					fan[i].on = 1;
					mvpow = mvpow + fan[i].power;				// vent fan power
					if(fan[i].q > 0) {							// supply fan - its heat needs to be added to the internal gains of the house
						fanheat = fan[i].power * .84;			// 16% efficiency for the particular fan used in this study.
						ventsumin = ventsumin + abs(fan[i].q) * 3600 / volume;
					}
					else											// exhasut fan
						ventsumout = ventsumout + abs(fan[i].q) * 3600 / volume;				
				} else if(fan[i].oper == 2) {					// FIXED SCHEDULE BATHROOM
					if(HOUR > 7.5 && HOUR <= 8) {
						fan[i].on = 1;
						mvpow = mvpow + fan[i].power;
						ventsumout = ventsumout + abs(fan[i].q) * 3600 / volume;
						nonRivecVentSumOut = nonRivecVentSumOut + abs(fan[i].q) * 3600 / volume;
					} else
						fan[i].on = 0;
				} else if(fan[i].oper == 3) {					// FIXED SCHEDULE KITCHEN
					if(HOUR > 17 && HOUR <= 18) {
						fan[i].on = 1;
						mvpow = mvpow + fan[i].power;			// vent fan power
						ventsumout = ventsumout + abs(fan[i].q) * 3600 / volume;
						nonRivecVentSumOut = nonRivecVentSumOut + abs(fan[i].q) * 3600 / volume;
					} else
						fan[i].on = 0;
				} else if(fan[i].oper == 4) {					// FIXED SCHEDULE EXHAUST FAN
					if(hcflag == 1) {
						if(HOUR > 1 && HOUR <= 5)
							fan[i].on = 0;
						else {
							fan[i].on = 1;
							mvpow = mvpow + fan[i].power;		// vent fan power
							ventsumout = ventsumout + abs(fan[i].q) * 3600 / volume;
							nonRivecVentSumOut = nonRivecVentSumOut + abs(fan[i].q) * 3600 / volume;
						}
					}
					if(hcflag == 2) {
						if(HOUR > 15 && HOUR <= 19) {
							fan[i].on = 0;
						} else {
							fan[i].on = 1;
							mvpow = mvpow + fan[i].power;		// vent fan power
							ventsumout = ventsumout + abs(fan[i].q) * 3600 / volume;
							nonRivecVentSumOut = nonRivecVentSumOut + abs(fan[i].q) * 3600 / volume;
						}
					}
				} else if(fan[i].oper == 23) {					// DYNAMIC SCHEDULE DRYER
					if(dryerFan == 1) {
						fan[i].on = 1;
						mvpow = mvpow + fan[i].power;
						ventsumout = ventsumout + abs(fan[i].q) * 3600 / volume;
						nonRivecVentSumOut = nonRivecVentSumOut + abs(fan[i].q) * 3600 / volume;
					} else
						fan[i].on = 0;
				} else if(fan[i].oper == 24) {					// DYNAMIC SCHEDULE KITCHEN FAN
					if(kitchenFan == 1) {
						fan[i].on = 1;
						mvpow = mvpow + fan[i].power;
						ventsumout = ventsumout + abs(fan[i].q) * 3600 / volume;
						nonRivecVentSumOut = nonRivecVentSumOut + abs(fan[i].q) * 3600 / volume;
					} else
						fan[i].on = 0;
				} else if(fan[i].oper == 25) {					// DYNAMIC SCHEDULE BATHROOM 1
					if(bathOneFan == 1) {
						fan[i].on = 1;
						mvpow = mvpow + fan[i].power;
						ventsumout = ventsumout + abs(fan[i].q) * 3600 / volume;
						nonRivecVentSumOut = nonRivecVentSumOut + abs(fan[i].q) * 3600 / volume;
					} else
						fan[i].on = 0;
				} else if(fan[i].oper == 26) {					// DYNAMIC SCHEDULE BATHROOM 2
					if(bathTwoFan == 1) {
						fan[i].on = 1;
						mvpow = mvpow + fan[i].power;
						ventsumout = ventsumout + abs(fan[i].q) * 3600 / volume;
						nonRivecVentSumOut = nonRivecVentSumOut + abs(fan[i].q) * 3600 / volume;
					} else
						fan[i].on = 0;
				} else if(fan[i].oper == 27) {					// DYNAMIC SCHEDULE BATHROOM 3
					if(bathThreeFan == 1) {
						fan[i].on = 1;
						mvpow = mvpow + fan[i].power;
						ventsumout = ventsumout + abs(fan[i].q) * 3600 / volume;
						nonRivecVentSumOut = nonRivecVentSumOut + abs(fan[i].q) * 3600 / volume;
					} else
						fan[i].on = 0;
				} else if(fan[i].oper == 19) {					// DRYER FAN
					if(weekendflag == 1) {
						if(HOUR > 12 && HOUR <= 15) {
							fan[i].on = 1;
							mvpow = mvpow + fan[i].power;
							ventsumout = ventsumout + abs(fan[i].q) * 3600 / volume;
							nonRivecVentSumOut = nonRivecVentSumOut + abs(fan[i].q) * 3600 / volume;
						} else
							fan[i].on = 0;
					} else
						fan[i].on = 0;
				} else if(fan[i].oper == 21) {					// Economizer
					if(econoflag == 1) {
						fan[i].on = 1;
						fanpow = fanpow + fan[i].power;			// vent fan power
						ventsumin = ventsumin + abs(fan[i].q) * 3600 / volume;
						nonRivecVentSumIn = nonRivecVentSumIn + abs(fan[i].q) * 3600 / volume;
					} else
						fan[i].on = 0;
				} else if(fan[i].oper > 30)
					dcv = 1;			
			}
			// *************************************** [RECHECK - END] **********************************************************>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

			// ***************  RIVEC ON/Off decision
			// ***************  RIVEC calculations
			if(nonRivecVentSumIn > nonRivecVentSumOut)						//ventsum based on largestof inflow or outflow
				nonRivecVentSum = nonRivecVentSumIn;
			else
				nonRivecVentSum = nonRivecVentSumOut;
			
			// RIVEC (only done evey ten minutes)
			for(int i=0; i < Nfan; i++) {
				// ---------------------FAN 31---------------------- DEFAULT RIVEC OPERATION BASED ON CONTROL ALGORITHM
				if(dcv == 1 && fan[i].oper == 31) {
					if(hourmin == 1 || hourmin == 11 || hourmin == 21 || hourmin == 31 || hourmin == 41 || hourmin == 51) {
						rivecOn = 1;			// this variable tells us if rivec is calling for the whole house fan to be turned off
						// start with controlled fan on as default					
						if(hcflag == 1) {									// HEATING TIMES
							if(HOUR >= 10 && HOUR < 22) {					// base time period
								if(reldose <= 1 && relexp <= .8)
									rivecOn = 0;
								//IF nonRivecVentSum < aeq THEN rivecOn = 1 (removed from revised algorithm)
							} else if(HOUR >= 22 || HOUR < 2) {				// pre-peak time period
								if(reldose <= 1 && relexp <= 1)
									rivecOn = 0;
								if(nonRivecVentSum > 1.25 * aeq)
									rivecOn = 0;
							} else if(HOUR >= 2 && HOUR < 6) {				// peak time period
									rivecOn = 0;							// always off
							} else {                                        // post-peak time period
								if(reldose <= 1 && relexp <= 1)
									rivecOn = 0;
								if(nonRivecVentSum > 1.25 * aeq)
									rivecOn = 0;
							}
						} else {											// COOLING TIMES
							if(HOUR >= 22 || HOUR < 10) {					// base time
								if(reldose <= 1 && relexp <= .8)
									rivecOn = 0;
								//IF nonRivecVentSum < aeq THEN rivecOn = 1 (removed from revised algorithm)
							} else if(HOUR >= 10 && HOUR < 14) { 			// prepeak
								if(reldose <= 1 && relexp <= 1)
									rivecOn = 0;
								if(nonRivecVentSum > 1.25 * aeq)
									rivecOn = 0;
							} else if(HOUR >= 14 && HOUR < 18) {			// peak
									rivecOn = 0;							// always off
							} else {										// post peak
								if(reldose <= 1 && relexp <= 1)
									rivecOn = 0;
								if(nonRivecVentSum > 1.25 * aeq)
									rivecOn = 0;
							}
						
						}
					}
					
					//RIVEC decision after all other fans added up
					if(rivecOn == 0)										//RIVEC has turned off this fan
						fan[i].on = 0;
					else {													//limited set of fans controlled by RIVEC. Only supply, exhaust and HRV are controlled
						fan[i].on = 1;
						mvpow = mvpow + fan[i].power;						// vent fan power
						if(fan[i].q > 0) {									// supply fan - its heat needs to be added to the internal gains of the house
							fanheat = fan[i].power * .84;					// 16% efficiency for the particular fan used in this study.
							ventsumin = ventsumin + abs(fan[i].q) * 3600 / volume;
						} else												// ex fan
							ventsumout = ventsumout + abs(fan[i].q) * 3600 / volume;
					}
				}
				//---------------------FAN 31 end----------------------



				// FF: FAN 32, 33, 34, 35 commented due to dependency on tdv/td file
				/*
				//---------------------FAN 32---------------------- Uses an input file to determine the base/pre/post and peak time periods
				if(dcv == 1 && fan[i].oper == 32) {
					if(daymin == 1)
						//INPUT #8, tdv 									//variable tdv is used for tdv 4hr data
					if(hourmin == 1 || hourmin == 11 || hourmin == 21 || hourmin == 31 || hourmin == 41 || hourmin == 51) {
						rivecOn = 1; 										// this variable tells us if rivec is calling for the whole house fan to be turned off
						// start with controlled fan on as default
						if(HOUR >= (tdv - 4) && HOUR < tdv) {     			// prepeak
							if(reldose <= 1 && relexp <= 1)
								rivecOn = 0;
							if(nonRivecVentSum > 1.25 * aeq)
								rivecOn = 0;
						} else if(HOUR >= tdv && HOUR < (tdv + 4)) { 		//peak
							rivecOn = 0;   									//always off
						} else if(HOUR >= (tdv + 4) && HOUR < (tdv + 8)) { 	// post peak
							if(reldose <= 1 && relexp <= 1)
								rivecOn = 0;
							if(nonRivecVentSum > 1.25 * aeq)
								rivecOn = 0;
						} else {
							if(reldose <= 1 && relexp <= .8)
								rivecOn = 0;
						}
					}
					//RIVEC decision after all other fans added up
					if(rivecOn == 0) {  										//rivec has turned off this fan
						fan[i].on = 0;
					} else {  												//limited set of fans controllede by RIVEC only supply, exhaust and HRV are controlled
						fan[i].on = 1;
						mvpow = mvpow + fan[i].power;						// vent fan power
						if(fan[i].q > 0) { 									// supply fan - its heat needs to be added to the internal gains of the house
							fanheat = fan[i].power * .84;					// 16% efficiency for the particular fan used in this study.
							ventsumin = ventsumin + abs(fan[i].q) * 3600 / volume;
						} else {											// ex fan
							ventsumout = ventsumout + abs(fan[i].q) * 3600 / volume;
						}
					}
				}
				//---------------------FAN 32 end-----------------------------
				
				//---------------------FAN 33 begin / TD ---------------------
				if(dcv == 1 && fan[i].oper == 33) {
					
					if(hourmin == 1)
					//	INPUT #9, td					
					if(hourmin == 1 || hourmin == 11 || hourmin == 21 || hourmin == 31 || hourmin == 41 || hourmin == 51) {
						rivecOn = 1;								//start with controlled fan on as default
						if(td == 1) {								//peak
							rivecOn = 0;							//always off
						} else {
							if(reldose <= 1 && relexp <= 1)
								rivecOn = 0;
						}
					}
				
					if(rivecOn == 0)								//rivec has turned off this fan
						fan[i].on = 0;
					else {
						fan[i].on = 1;
						mvpow = mvpow + fan[i].power;				//vent fan power
						if(fan[i].q > 0) {							//supply fan - its heat needs to be added to the internal gains of the house
							fanheat = fan[i].power * .84;			// 16% efficiency for the particular fan used in this study.
							ventsumin = ventsumin + abs(fan[i].q) * 3600 / volume;
						} else {									//ex fan
							ventsumout = ventsumout + abs(fan[i].q) * 3600 / volume;
						}
					}
				}
				//---------------------FAN 33 end-----------------------------
				
				//---------------------FAN 34 begin / TD ---------------------
				if(dcv == 1 && fan[i].oper == 34) {
					if(hourmin == 1)
						// INPUT #9, td
					if(hourmin == 1 || hourmin == 11 || hourmin == 21 || hourmin == 31 || hourmin == 41 || hourmin == 51) {
						rivecOn = 1; //start with controlled fan on as default
						if(td == 1) { //peak
							rivecOn = 0;   //always off
							if(hcflag == 1 && (Tpredhouse - Tout) < 5) {
								rivecOn = 1;
							} else if(hcflag == 2 && (Tpredhouse - Tout) > 0) {
								rivecOn = 1;
							}
						} else {
							if(reldose <= 1 && relexp <= 1)
								rivecOn = 0;
						}
					}
				
					if(rivecOn == 0)		//rivec has turned off this fan
						fan[i].on = 0;
					else {
						fan[i].on = 1;
						mvpow = mvpow + fan[i].power;				//vent fan power
						if(fan[i].q > 0) {							//supply fan - its heat needs to be added to the internal gains of the house
							fanheat = fan[i].power * .84;			// 16% efficiency for the particular fan used in this study.
							ventsumin = ventsumin + abs(fan[i].q) * 3600 / volume;
						} else {									//ex fan
							ventsumout = ventsumout + abs(fan[i].q) * 3600 / volume;
						}
					}
				}
				//---------------------FAN 34 end------------------------------
				
				//---------------------FAN 35 begin / TDe ---------------------
				if(dcv == 1 && fan[i].oper == 35) {
					if(hourmin == 1)
						// INPUT #9, td
					if(hourmin == 1 || hourmin == 11 || hourmin == 21 || hourmin == 31 || hourmin == 41 || hourmin == 51) {
						rivecOn = 1;           							//start with controlled fan on as default
						if(td == 1) {          							//peak
							rivecOn = 0;    							//always off
							if(hcflag == 1 && (Tpredhouse - Tout) < 5) {
								rivecOn = 1;
							} else if(hcflag == 2 && (Tpredhouse - Tout) > 5) {
								rivecOn = 1;
							}
						} else {
							if(reldose <= 1 && relexp <= 1) {
								rivecOn = 0;
							}
						}
					}

					if(rivecOn == 0)									//rivec has turned off this fan
						fan[i].on = 0;
					else {
						fan[i].on = 1;
						mvpow = mvpow + fan[i].power;    				// vent fan power
						if(fan[i].q > 0) {            					// supply fan - its heat needs to be added to the internal gains of the house
							fanheat = fan[i].power * .84;				// 16% efficiency for the particular fan used in this study.
							ventsumin = ventsumin + abs(fan[i].q) * 3600 / volume;
						} else { 										//ex fan
							ventsumout = ventsumout + abs(fan[i].q) * 3600 / volume;
						}
					}
				}
				//---------------------FAN 35 end----------------------
				*/

				//---------------------FAN 36 rivec(4h - no shoulder periods)----------------------
				if(dcv == 1 && fan[i].oper == 36) {
					if(hourmin == 1 || hourmin == 11 || hourmin == 21 || hourmin == 31 || hourmin == 41 || hourmin == 51) {
						rivecOn = 1; 									// this variable tells us if rivec is calling for the whole house fan to be turned off
						// start with controlled fan on as default
					
						if(hcflag == 1) {              					// heating times
							if(HOUR >= 2 && HOUR < 6)					// peak
								rivecOn = 0;            				// always off
							else {                    			// post peak
								if(reldose <= 1 && relexp <= 1)
									rivecOn = 0;
							}
						} else {										// cooling times
							if(HOUR >= 14 && HOUR < 18)					// peak
								rivecOn = 0;                            // always off
							else {                            			// post peak
								if(reldose <= 1 && relexp <= 1)
									rivecOn = 0;
							}
						}
					}
					//RIVEC decision after all other fans added up
					//RIVEC MUST be the third fan
					//i = 3
					if(rivecOn == 0)	  								//rivec has turned off this fan
						fan[i].on = 0;
					else {  											//limited set of fans controllede by RIVEC only supply, exhaust and HRV are controlled
						//IF fan(i).oper = 1 THEN   // fan always on
						fan[i].on = 1;
						mvpow = mvpow + fan[i].power;					// vent fan power
						if(fan[i].q > 0) { 								// supply fan - its heat needs to be added to the internal gains of the house
							fanheat = fan[i].power * .84;				// 16% efficiency for the particular fan used in this study.
							ventsumin = ventsumin + abs(fan[i].q) * 3600 / volume;
						} else {										// ex fan
							ventsumout = ventsumout + abs(fan[i].q) * 3600 / volume;
						}
					//END IF
					}
				}
				//---------------------FAN 36 end----------------------

				//---------------------FAN 37 - RIVEC new heating hours----------------------
				if(dcv == 1 && fan[i].oper == 37) {
					if(hourmin == 1 || hourmin == 11 || hourmin == 21 || hourmin == 31 || hourmin == 41 || hourmin == 51) {
						rivecOn = 1;									// this variable tells us if rivec is calling for the whole house fan to be turned off
						// start with controlled fan on as default

						if(hcflag == 1) {								//heating times
							if(HOUR >= 1 && HOUR < 5) {					// prepeak
								if(reldose <= 1 && relexp <= 1)
									rivecOn = 0;
								if(nonRivecVentSum > 1.25 * aeq)
									rivecOn = 0;
							} else if( HOUR >= 5 && HOUR < 9) {			//peak
								rivecOn = 0;							//always off
							} else if( HOUR >= 9 && HOUR < 13) {		// post peak
								if(reldose <= 1 && relexp <= 1)
									rivecOn = 0;
								if(nonRivecVentSum > 1.25 * aeq)
									rivecOn = 0;
							} else {
								if(reldose <= 1 && relexp <= .8) {
									rivecOn = 0;
								}
								//IF nonRivecVentSum < aeq THEN rivecOn = 1
							}

						} else {										//cooling times
							if(HOUR >= 22 || HOUR < 10) {				// base time
								if(reldose <= 1 && relexp <= .8)
									rivecOn = 0;
								//IF nonRivecVentSum < aeq THEN rivecOn = 1
							} else if( HOUR >= 10 && HOUR < 14) {		// prepeak
								if(reldose <= 1 && relexp <= 1)
									rivecOn = 0;
								if(nonRivecVentSum > 1.25 * aeq)
									rivecOn = 0;
							} else if( HOUR >= 14 && HOUR < 18) {		//peak
								rivecOn = 0;							//always off
							} else {									// post peak
								if(reldose <= 1 && relexp <= 1)
									rivecOn = 0;
								if(nonRivecVentSum > 1.25 * aeq)
									rivecOn = 0;
							}
						}
					}
					//RIVEC decision after all other fans added up
					if(rivecOn == 0)									//rivec has turned off this fan
						fan[i].on = 0;
					else {												//limited set of fans controlled by RIVEC only supply, exhaust and HRV are controlled
						//IF fan(i).oper = 1 THEN   // fan always on
						fan[i].on = 1;
						mvpow = mvpow + fan[i].power;					// vent fan power
						if(fan[i].q > 0) {								// supply fan - its heat needs to be added to the internal gains of the house
							fanheat = fan[i].power * .84;				// 16% efficiency for the particular fan used in this study.
							ventsumin = ventsumin + abs(fan[i].q) * 3600 / volume;
						} else { // ex fan
							ventsumout = ventsumout + abs(fan[i].q) * 3600 / volume;
						}
						//END IF
					}
				}
				//---------------------FAN 37 end----------------------

				//---------------------FAN 40 - dcv by dose----------------------
				if(dcv == 1 && fan[i].oper == 40) {
					if(hourmin == 1 || hourmin == 11 || hourmin == 21 || hourmin == 31 || hourmin == 41 || hourmin == 51) {
						rivecOn = 1; 								//start with controlled fan on as default
						if(reldose <= 1)
							rivecOn = 0;
					}

					if(rivecOn == 0)  								//rivec has turned off this fan
						fan[i].on = 0;
					else {  										//limited set of fans controlled by RIVEC only supply, exhaust and HRV are controlled
						fan[i].on = 1;
						mvpow = mvpow + fan[i].power;				// vent fan power
						if(fan[i].q > 0) { 							// supply fan - its heat needs to be added to the internal gains of the house
							fanheat = fan[i].power * .84;			// 16% efficiency for the particular fan used in this study.
							ventsumin = ventsumin + abs(fan[i].q) * 3600 / volume;
						} else { 									// ex fan
							ventsumout = ventsumout + abs(fan[i].q) * 3600 / volume;
						}
					}
				}
				//---------------------FAN 40 end----------------------

				//---------------------FAN 41 - dcv by exposure----------------------
				if(dcv == 1 && fan[i].oper == 41) {
					if(hourmin == 1 || hourmin == 11 || hourmin == 21 || hourmin == 31 || hourmin == 41 || hourmin == 51) {
						rivecOn = 1; 								//start with controlled fan on as default
						if(relexp <= 1)
							rivecOn = 0;
					}

					if(rivecOn == 0)  								//rivec has turned off this fan
						fan[i].on = 0;
					else {  										//limited set of fans controlled by RIVEC only supply, exhaust and HRV are controlled
						fan[i].on = 1;
						mvpow = mvpow + fan[i].power;   	 		// vent fan power
						if(fan[i].q > 0) {                    		// supply fan - its heat needs to be added to the internal gains of the house
							fanheat = fan[i].power * .84;			// 16% efficiency for the particular fan used in this study.
							ventsumin = ventsumin + abs(fan[i].q) * 3600 / volume;
						} else { 									// ex fan
							ventsumout = ventsumout + abs(fan[i].q) * 3600 / volume;
						}
					}
				}
				//---------------------FAN 41 end----------------------

				//---------------------FAN 50---------- RIVEC OPERATION BASED ON CONTROL ALGORITHM v3 WJNT May 2011
				if(dcv == 1 && fan[i].oper == 50) {
					if(hourmin == 1 || hourmin == 11 || hourmin == 21 || hourmin == 31 || hourmin == 41 || hourmin == 51) {
						rivecOn = 1; 												// this variable tells us if rivec is calling for the whole house fan to be turned on
						// start with controlled fan on (1) as default

						if(hcflag == 1) {      										// HEATING TIMES
							if((HOUR >= 0 && HOUR < 4) || (HOUR >= 12)) {        	// BASE Time Period
								if(occupied[HOUR] == 1) {                        	// Occupied
									if(reldose < 1 && relexp < .8)
										rivecOn = 0;
									if(nonRivecVentSum < aeq)
										rivecOn = 1;
								} else if(occupied[HOUR] == 0) {                	// Unoccupied
									if(reldose < 1 && relexp < 2)
										rivecOn = 0;
									if(reldose > 1.2 || relexp > 2)
										rivecOn = 1;
									if(nonRivecVentSum > aeq)
										rivecOn = 0;
								}
							} else if(HOUR >= 4 && HOUR < 8) {						// PEAK Time Period
								rivecOn = 0;										// Always off
							} else if(HOUR >= 8 && HOUR < 12) {         			// RECOVERY Time Period
								if(reldose < 1 && relexp < 1)
									rivecOn = 0;
								if(nonRivecVentSum > 1.25 * aeq)
									rivecOn = 0;
							}

						} else {    												// COOLING TIMES
							if(HOUR >= 22 || (HOUR >= 0 && HOUR < 14)) {           	// BASE Time Period
								if(occupied[HOUR]) {                            	// Occupied
									if(reldose < 1 && relexp < .8)
										rivecOn = 0;
									if(nonRivecVentSum < aeq)
										rivecOn = 1;
								} else {                                           	// Unoccupied
									if(reldose < 1 && relexp < 2)
										rivecOn = 0;
									if(reldose > 1.2 || relexp > 2)
										rivecOn = 1;
									if(nonRivecVentSum > aeq)
										rivecOn = 0;
								}
							} else if(HOUR >= 14 && HOUR < 18) {                 	// PEAK Time Period
								rivecOn = 0;                                		// Always off
							} else if(HOUR >= 18 && HOUR < 22) {         			// RECOVERY Time Period
								if(reldose < 1 && relexp < 1)
									rivecOn = 0;
								if(nonRivecVentSum > 1.25 * aeq)
									rivecOn = 0;
							}
						}
					}
					// RIVEC decision after all other fans added up
					if(rivecOn == 0)	  											// RIVEC has turned off this fan
						fan[i].on = 0;
					else {  														// limited set of fans controlled by RIVEC. Only supply, exhaust and HRV are controlled
						fan[i].on = 1;
						mvpow = mvpow + fan[i].power;								// vent fan power
						if(fan[i].q > 0) { 											// supply fan - its heat needs to be added to the internal gains of the house
							fanheat = fan[i].power * .84;							// 16% efficiency for the particular fan used in this study.
							ventsumin = ventsumin + abs(fan[i].q) * 3600 / volume;
						} else { 													// ex fan
							ventsumout = ventsumout + abs(fan[i].q) * 3600 / volume;
						}
					}
				}
				// ---------------------FAN 50 end----------------------
			} // this is end of FOR loop for RIVEC fans statements

			for(int i=0; i < Nfan; i++) {
				// ************************** HRV code for HRV synched to air handler
				if(fan[i].oper == 16) {												// HRV+Air Handler
					if(rivecOn == 1) {												// rivec controller allows HRV operation
						ventsumin = ventsumin + abs(fan[i].q) * 3600 / volume;
						ventsumout = ventsumout + abs(fan[i].q) * 3600 / volume;
						fan[i].on = 1;
						mvpow = mvpow + fan[i].power;
						if(ahflag != 1 && ahflag != 2 && ahflag != 102) {			// ah is off for heat/cool so turn it on for venting
							ahflag = 100;
							ahminutes = ahminutes + 1;
							qah = cqah;
							SUPVELah = qah / (pow(supdiam,2) * pi / 4);
							retvelah = qah / (pow(retdiam,2) * pi / 4);
							QsupREG = -qah * supLF + qah;
							QretREG = qah * retLF - qah;
							qretleak = -qah * retLF;
							qsupleak = qah * supLF;
							MSUPLEAK = qsupleak * rhoin;
							mretleak = qretleak * rhoin;

							//*** HRV air is from outside so not simply added to qretleak or mretleak
							//*** but needs to be subtracted from mretreg1
							MAH = qah * rhoin;
							mhrv = fan[i].q * rhoin;								// mhrv needs to be a negative number - so does fan(i).q to get right mretreg and then fan(i).q is also the exhaust from house flow
							MSUPREG = QsupREG * rhoin;
							mretreg = QretREG * rhoin - mhrv;
							mvpow = mvpow + cfanpower;								// note ah operates at cooling speed
							fanpower = cfanpower * .85;								// for heat
						} else {													// open the outside air vent
							mhrv = fan[i].q * rhoin;
							mretreg = QretREG * rhoin - mhrv;
						}
					}
				}
			}

			// ***********************************************************************
			// ********************* equipment capacity calculations
			// NEW cooling EQUIPMENT MODEL 2005 - with 2006/7 additions for SHR
			evapcap = 0;
			latcap = 0;
			capacity = 0;
			capacityh = 0;
			powercon = 0;
			powercomp = 0;

			if(ahflag == 0) {														// AH OFF
				capacityc = 0;
				capacityh = 0;
				powercon = 0;
				powercomp = 0;
			} else if(ahflag == 100) {												// Fan on no heat/cool
				capacityh = fanpower;												// Include fan heat
				capacityc = 0;
			} else if(ahflag == 102) {												// End if heat cycle
				capacityh = .25 * hcapacity + fanpower;								// NOTE 0.25 refers to cycling for 15 minutes out of 60?
				capacityc = 0;
			} else if(ahflag == 1) {												// Heating
				capacityh = hcapacity + fanpower;									// include fan heat in furnace capacity (converted to W in main program)
				capacityc = 0;
			} else {																// we have cooling
				Toutf = (Tout - 273.15) * 1.8 + 32;
				tretf = (told[11] - 273.15 + fanpower / 1005 / MAH) * 1.8 + 32;		// return in F used for capacity includes fan heat
				dhret = fanpower / MAH / 2326;										// added to hret in capacity caluations for wet coil - converted to Btu/lb

				// wret=(hret - 1006*(told(12)-273.15))/(2501000+1860*(told(12)-273.15))  //enthalpy at return
				// the following corerctions are for TXV only
				// test for wet/dry coil using SHR calculation
				
				SHR = 1 - 50 * (HR[1] - .005);										// using humidity ratio - note only because we have small indoor dry bulb range
				
				if(SHR < .25)
					SHR = .25;  //setting low limit
				if(SHR > 1)
					SHR = 1;
				if(SHR == 1) { // dry coil
					// charge correction for capacity (dry coil)
					if(charge < .725) {
						chargecapd = 1.2 + (charge - 1);
					} else if(charge >= .725) {
						chargecapd = .925;
					}
					// charge correction for EER (dry coil)
					if(charge < 1) {
						chargeeerd = 1.04 + (charge - 1) * .65;
					} else if(charge >= 1) {
						chargeeerd = 1.04 - (charge - 1) * .35;
					}
					capacity = (capacityari * 1000) * .91 * qahcorrec * chargecapd * ((-.00007) * pow((Toutf - 82),2) - .0067 * (Toutf - 82) + 1);
					// EER for dry coil
					EER = EERari * qahcorrec * chargeeerd * ((-.00007) * pow((Toutf - 82),2) - .0085 * (Toutf - 82) + 1);
				} else {															// wet coil
					// charge correction for capacity (wet coil)
					if(charge < .85) {
						chargecapw = 1 + (charge - .85);
					} else if(charge >= .85) {
						chargecapw = 1;
					}
					// charge correction for EER(wet coil)
					if(charge < .85) {
						chargeeerw = 1 + (charge - .85) * .9;
					} else if(charge >= .85 && charge <= 1) {
						chargeeerw = 1;
					} else if(charge > 1) {
						chargeeerw = 1 - (charge - 1) * .35;
					}
					capacity = capacityari * 1000 * (1 + ((hret + dhret) - 30) * .025) * qahcorrec * chargecapw * ((-.00007) * pow((Toutf - 95),2) - .0067 * (Toutf - 95) + 1);
					// EER for wet coil
					EER = EERari * qahcorrec * chargeeerw * ((-.00007) * pow((Toutf - 95),2) - .0085 * (Toutf - 95) + 1);
				}

				// split sensible and latent capacities and convert capacity Btu/h to W
				// calculation of power used and total cacapity for air conditioning
				powercomp = capacity / EER; //compressor power

				// SHR ramps up for first three minutes of operation
				// comptime tracks the number of minutes the compressor is on for
				if(comptime == 1) {
					SHR = SHR + 2 / 3 * (1 - SHR);
				} else if(comptime == 2) {
					SHR = SHR + 1 / 3 * (1 - SHR);
				}

				capacityc = SHR * capacity / 3.413;									// sensible (equals total if SHR=1)
				// correct the sensible cooling for fan power heating
				capacityc = capacityc - fanpower;
				// tracking mass on coil (Mcoil) including condensation and evaporation - do this in main program
				evapcap = 0;
				latcap = 0;
			}

			// **************  tracking coil moisture
			if(ahflag == 2) {
				if(SHR < 1) {
					latcap = (1 - SHR) * capacity / 3.413;							// latent capacity J/s
					Mcoil = Mcoilprevious + latcap / 2501000 * dtau;				// condenstion in timestep kg/s * time
				} else {
					if(Mcoil > 0) {													// if moisture on coil but no latcap, then we evaporate coil moisture until Mcoil = zero
						latcap = -.3 * capacityraw / 1800 * 2501000;				// evaporation capacity J/s - this is negative latent capacity
						evapcap = latcap;
						Mcoil = Mcoilprevious - .3 * capacityraw / 1800 * dtau;		// evaporation in timestep kg/s * time
					} else {														// no latcap and dry coil
						latcap = 0;
					}
				}
			} else if(ahflag == 100) {												// then we have no cooling, but the fan is on
				if(Mcoil > 0) {														// if moisture on coil but no latcap, then we evaporate coil moisture until Mcoil = zero
					latcap = -.3 * capacityraw / 1800 * 2501000;					// evaporation capacity J/s - this is negative latent capacity
					evapcap = latcap;
					Mcoil = Mcoilprevious - .3 * capacityraw / 1800 * dtau;			// evaporation in timestep kg/s * time
				} else {															// no latcap and dry coil
					latcap = 0;
				}
			}

			if(Mcoil < 0)
				Mcoil = 0;
			if(Mcoil > .3 * capacityraw)
				Mcoil = .3 * capacityraw;											// maximum mass on coil is 0.3 kg per ton of cooling
			Mcoilprevious = Mcoil;													// maybe put this at top of the hour
			// ********************** end of equipment model


			// ********************************** moisture balance
			for(int i=0; i < 5; i++) {
				hrold[i] = HR[i];
			}

			moisture(HR, hrold, M1, M12, M15, M16, Mw5, dtau, matticenvout, Mceil, msupahoff, mretahoff,
				matticenvin, HROUT, MSUPLEAK, MAH, mretreg, mretleak, MSUPREG, latcap, mhousein, mhouseout,
				latload, mhrv, mfancyc, merv, MWha);

			if(HR[1] == 0)
				HR[1] = HROUT;
			if(HR[3] > .02)
				HR[3] = .02;

			// hret is in btu/lb and is used in capacity calculations
			hret = .24 * ((b[11] - 273.15) * 9 / 5 + 32) + HR[1] * (1061 + .444 * ((b[11] - 273.15) * 9 / 5 + 32));

			// *************** begin heat and mass transport calculations
			Mceilold = -1000;														// inital guess
			itercount = 0;
			limit = C / 10;
			if(limit <= .00001)
				limit = .00001;

			while(1) {												// this loop does the ventilation and heat transfer calculations
				itercount = itercount + 1;							// counting # of temperature/ventilation iterations
				flag = 0;

				while(1) {
					// house air flow calculations:
					houseleak(ahflag, flag, speed, direction, Tpredhouse, tpredattic, Tout, C, n, h, R, X, Nrflue,
						Chim, Wallfraction, FloorFraction, Sw, Swrfl, Nwindoor, Windoor, Nfan, fan, Npipe,
						Pipe, min, mout, Pint, Mflue, Mceil, Mfloor, Cattic, dPflue, dPceil, dPfloor, Crawl,
						Hfloor, row, Soffitfraction, Patticint, Cpwall, rhoref, MSUPREG, MAH, mretleak, MSUPLEAK,
						mretreg, mhousein, mhouseout, supC, supn, retC, retn, msupahoff, mretahoff, aeq);

					flag = flag + 1;

					if(abs(Mceilold - Mceil) < limit || flag > 5)
						break;
					else
						Mceilold = Mceil;
										
					// attic air flow caculations:
					Atticleak(flag, speed, direction, Tpredhouse, Tout, tpredattic, Cattic, Nattic, h, Hpeak,
						Swrfl, Sw, Navent, Avent, Soffit, Matticin, Matticout, Patticint, Mceil, row,
						Soffitfraction, pitch, RR, Nafans, AFan, rhoref, MSUPREG, mretleak, MSUPLEAK, matticenvin,
						matticenvout, dtau, msupahoff, mretahoff);
				}

				// HEAT TRANSFER PART: MCEIL is calculated in houseleak and used in Atticleak
				RHOATTIC = rhoref * 293 / tpredattic;

				// adding fan heat for supply fans, intgain1 is from input file, fanheat reset to zero each minute, intgain is common
				intgain = intgain1 + fanheat;

				bsize = sizeof(b)/sizeof(b[0]);

				heat(Tout, rhoref, Mceil, AL4, speed, ssolrad, nsolrad, told, ATTICVOL, housevol, sc, b, ERRCODE, TSKY,
					FLOORAREA, pitch, location, MSUPREG, mretreg, mretleak, MSUPLEAK, MAH, supRval, retRval, supdiam,
					retdiam, suparea, retarea, supthick, retthick, supvol, retvol, supcp, retcp, SUPVEL, retvel, suprho,
					retrho, pref, HROUT, diffuse, UA, matticenvin, matticenvout, mhousein, mhouseout, planarea, msupahoff,
					mretahoff, solgain, windowS, windowN, windowWE, shadcoef, mfancyc, Hpeak, h, Whouse, retlength, suplength,
					rooftype, M1, M12, M15, M16, rroof, rceil, ahflag, dtau, merv, mhrv, SBETA, CBETA, L, dec, Csol, idirect,
					capacityc, capacityh, evapcap, intgain, bsize);

				if(abs(b[0] - tpredattic) < .2) {										// Testing for convergence
					Tpredceil = b[6];
					TPREDFLOOR = b[7];
					tpredattic = b[0];
					TPREDSI = b[3];
					TPREDSO = b[4];
					TPREDNI = b[1];
					TPREDNO = b[2];
					TPREDENDI = b[8];
					TPREDENDO = b[9];
					TPREDWOOD = b[5];
					TPREDRETSURF = b[10];
					TPREDSUPSURF = b[13];
					tpredret = b[11];
					tpredsup = b[14];
					Tpredhouse = b[15];
					tpredhousemass = b[12];
					break;
				}

				if(itercount > 10) {													// Assume convergence
					Tpredceil = b[6];
					TPREDFLOOR = b[7];
					tpredattic = b[0];
					TPREDSI = b[3];
					TPREDSO = b[4];
					TPREDNI = b[1];
					TPREDNO = b[2];
					TPREDENDI = b[8];
					TPREDENDO = b[9];
					TPREDWOOD = b[5];
					TPREDRETSURF = b[10];
					TPREDSUPSURF = b[13];
					tpredret = b[11];
					tpredsup = b[14];
					Tpredhouse = b[15];
					tpredhousemass = b[12];
					break;
				}
				Tpredhouse = b[15];
				tpredattic = b[0];
			}
			// PRINT Mflue

			// setting "old" temps for next timestep to be current temps:
			told[0] = tpredattic;
			told[1] = TPREDNI;
			told[2] = TPREDNO;
			told[3] = TPREDSI;
			told[4] = TPREDSO;
			told[5] = TPREDWOOD;
			told[6] = Tpredceil;
			told[7] = TPREDFLOOR;
			told[8] = TPREDENDI;
			told[9] = TPREDENDO;
			told[10] = TPREDRETSURF;
			told[11] = tpredret;
			told[12] = tpredhousemass;
			told[13] = TPREDSUPSURF;
			told[14] = tpredsup;
			told[15] = Tpredhouse;


		
			// ************** house ventilation rate  - what would be measured with a tracer gas i.e., not just envelope and vent fan flows
			// min has msupreg added in mass balance calculations and mretleak contributes to house ventilation rate
			mhouse = min - mretleak * (1 - supLF) - MSUPREG;
			mhouse = mhouse - mfancyc - merv - mhrv;
			qhouse = abs(mhouse / rhoin); // m^3/s
			achhouse = qhouse / housevol * 3600; //ach
		
			if(Mflue < 0)
				achflue = (Mflue / rhoin) / housevol * 3600; //ach
			else
				achflue = (Mflue / rhoout) / housevol * 3600; //ach
		
			// ****** IAQ calculations **************************************************************************************************************************************************
			if(ventsumin > ventsumout) //ventsum based on largest of inflow or outflow
				ventsum = ventsumin;
			else
				ventsum = ventsumout;
			
			// -------dkm: To add flue in redose/relexp calc the following two IF sentences have been added------
			if(achflue <= 0)
				ventsumoutD = ventsumout + abs(achflue);
			else
				ventsuminD = ventsumin + achflue;
			
			if(ventsuminD > ventsumoutD)	//ventsum based on largest of inflow or outflow
				ventsumD = ventsuminD;
			else
				ventsumD = ventsumoutD;
		
			if(ventsumD <= 0)
				ventsumD = .000001;
			
			turnover = (1 - exp(-ventsumD * rivecdt)) / ventsumD + turnoverold * exp(-ventsumD * rivecdt);
			relexp = aeq * turnover;
			reldose = aeq * turnover * (1 - exp(-rivecdt / 24)) + reldoseold * exp(-rivecdt / 24);
		
			// ****** end IAQ calculations **************************************************************************************************************************************************
			
			// Writing results to a csv file
			double setpoint;
		
			if(hcflag == 1)
				setpoint = theat[ttime-1];
			else
				setpoint = tcool[ttime-1];

			outputfile << HOUR << ", " << hrtot << ", " << Tout << ", " << tpredattic << ", " << Tpredhouse << ", ";
			outputfile << tpredsup << ", " << tpredret << ", " << ahflag << ", " << fanpow << ", " << hcap << ", ";
			outputfile << powercomp << ", " << mvpow << ", " << HR[3] * 1000 << ", " << qhouse << ", " << achhouse << ", ";
			outputfile << capacityc << ", " << SHR << ", " << Mcoil << ", " << setpoint << ", " << Pint << ", ";
			outputfile << turnover << ", " << relexp << ", " << reldose << ", " << fan[0].on << ", " << fan[1].on << ", ";
			outputfile << fan[2].on << ", " << fan[3].on << ", " << fan[4].on << ", " << fan[5].on << ", " << ventsum << ", ";
			outputfile << rivecOn << ", " << nonRivecVentSum << ", " << achflue << ", " << occupied[HOUR] << endl;
		
			//****** WRITING MOISTURE DATA FILE ******
			moisturefile << HR[0] << ", " << HR[1] << ", " << HR[2] << ", " << HR[3] << ", " << HR[4] << endl;
			
			// Calculating sums for electrical and gas energy use
			gastherm = gastherm + hcap * 60 / 1000000 / 105.5;
			faneleckWh = faneleckWh + fanpow / 60000;
			compeleckWh = compeleckWh + powercomp / 60000;
			mveleckWh = mveleckWh + mvpow / 60000;
			numt = numt + 1;
			toutavg = toutavg + Tout;
			tatticavg = tatticavg + tpredattic;
			thouseavg = thouseavg + Tpredhouse;
			achavg = achavg + achhouse;
			daymin++;									// Minute count (of day)
			
			// To loop for a specific day count:
			//- set the ammount of days in variable "totaldays" at the start of the code
			//- comment the "} while (ws4file);" that is locates after this batch of comments
			//- uncomment the following 3 lines:
			//
			//		if(hrtot >= (totaldays) * 1440 + 1)
			//			break;
			//	} while (day >= totaldays);
			//
			// To loop until end of weather file:
			//- comment the previously mentioned lines
			//- uncomment the following line:
			//	} while (ws4file);

			if(hrtot >= (totaldays) * 1440 + 1)
				break;

		} while (day <= totaldays);


		//----------dkm: TIME for simulation--------
		time(&endtime);
		string te = ctime(&endtime);

		cout << "\nStart of simulation = " << ts;
		cout << "End of simulation = " << te << endl;
		//cout << "Time for one year = " << (ts - te) * 365 / day << endl;		// doesnt work in C++
		//-----------------------------------------

		outputfile.close();
		ws4file.close();

		if(dynamicScheduleFlag == 1)								// Fan Schedule
			fanschedulefile.close();

		//CLOSE #9 //tdv 4hr-file
		toutavg = toutavg / numt;
		tatticavg = tatticavg / numt;
		thouseavg = thouseavg / numt;
		achavg = achavg / numt;

		_chdir("..\\output\\");
		ofstream ou2file(Outfile + ".ou2"); 
		if(!ou2file) { 
			cout << "Cannot open: " << Outfile + ".ou2" << endl;
			system("pause");
			return 1; 
		}
		
		ou2file << "* Input file = " << f << " Weather file = " << weatherfile << endl;
		ou2file << "Tout, Tattpred, Thousepred, AHkWh, therms, compkWh, ventkWh, avgach" << endl;
		ou2file << toutavg << ", " << tatticavg << ", " << thouseavg << ", " << faneleckWh << ", ";
		ou2file << gastherm << ", " << compeleckWh << ", " << mveleckWh << ", " << achavg << endl;
		ou2file << "Number of days = " << day - 1 << endl;
		ou2file << "Start time = " << ts;
		ou2file << "End time = " << te << endl;

		ou2file.close();

		// ****** CLOSING MOISTURE DATA FILE ******
		moisturefile.close();	
	}
	system("pause");

	return 0; 
}
