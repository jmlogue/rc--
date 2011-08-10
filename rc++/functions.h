#pragma once
#ifndef functions_h
#define functions_h

#include <iostream>
#include <string>

using namespace std;

// additional data types
struct AventType {
	int wall;
	double h;
	double A;
	double n;
	double m;
	double dP;
};

struct SoffitType {
	double h;
	double m;
	double dP;
};

struct windoortype {
	int wall;
	double Bottom;
	double Top;
	double High;
	double Wide;
	double m;
	double min;
	double mout;
	double dPtop;
	double dPbottom;
};

struct fantype {
	double power;
	double q;
	double m;
	double on; // Set on=1 in main program based on .oper
	double oper; // 1= const operation, 2= fixed schedule bath, 3= fixed sched kitchen, 4 = fixed sched ex, 5 = hrv, 6 = fan cycler exhaust
};

struct pipetype {
	int wall;
	double h;
	double A;
	double n;
	double m;
	double dP;
	double Swf;
	double Swoff;
};

struct Chimtype {
	double Cflue;
	double Hflue;
	double Fluetemp;
};

// Additional functions

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
);

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
);

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
);

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
);

#endif
