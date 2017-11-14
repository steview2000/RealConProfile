/**********************************************************************
sf6_critical.c copied from Jae 12-15-01 fixed bug in CP 12-2-02 GA

This program calculates thermodynamic and kinetic properties of SF6 gas (SI unit.) near the critical region.

Input should be T in deg C and P in psi and all returned values (to the main program - noise.c or RBC.c) are in cgs.
But sf6_critical.c plays in SI unit.

The thermodynamic properties are from the equation of state (the crossover equation) in:
A. K. Wyczalkowska and J. V. Sengers, Journal of Chemical Physics, V.111, p.1551 (1999)

The thermodynamic properties in this program are valid in the following region:
	(1./_xi) <2.0,
    and 336kg/m^3 < rho < 1170kg/m^3,
    and 318.717K < T < 390K at rho=RHOc
    On the critical isochore, this corresponds to 544.55 < P < 1456.7 psi
    
If the fluid properties in the cell are out of this region, the program stops.

The viscosity is calculated by fitting experimental data to a formula which is a function of temperature and density.
The detailed are shown in Note.1 (see the last lines of the program).

The thermal conductivity is summation of the regular part and the singular critical part.
The regular (background) thermal conductivity is calculated by the concept of the excess thermal conductivity.
The singular (critical) thermal conductivity near critical point is basically obtained from the critical isochoric thermal conductivity data since no off-density thermal conductivity data are available.
The details are shown in Note.2 in the last lines of this program.

by Jaechul Oh, 2/10/2000

*  	rho     =	density			g/cm3                           *
*	alfa 	=	thermal expansion coeff	1/K                             *
*	comp	=	isothermal compressibility cm^2/dyne       		*
*	lambda	= 	thermal conductivity 	erg/s cm K                      * 
*	kappa 	=	thermal diffusivity	cm2/s                           *
*	nu 	=	viscosity 		g /s cm                         * 
*	cp 	=   	spec heat const press	erg/g K				*

***********************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define ZERO_C 		273.15
#define R 		8314.33    	/*  J/kmole*K  */

#define Tc	318.717 /* critical temperature of SF6, Kelvin = 45.567 deg C */
#define Pc	3754500. /* critical perssure of SF6 in Pa = 544.55 psi, (1bar=1.0e5Pa, 1Pa=1N/M^2, 1atm=101325Pa) */
#define RHOc	742. /* critical density of SF6, kg/m^3) */
#define SF6_MW	146.054 //9-25-2000  //146.066 /* molecular weight of SF6, kg/kmol */

/* Critical exponents */
#define alpha_	0.110 //+-0.003
#define beta_	0.326 //+-0.002
#define gamma_	1.239 //+-0.002
#define nu_	0.630 //+-0.001
#define eta_	0.03333 //0.033 //+-0.004
#define DelS_	0.50 //9-25-2000  //0.51  //0.54  //+-0.04
#define DelA_	1.32  //+-0.2

/* Crossover parameters */
#define ubar_	0.50043
#define uaster_	0.472    //a universal fixed-point coupling constant
//#define u_	0.236203 //=ubar*uaster
#define GAMMA_	0.80621

/* Scaling-field parameters */
#define ct_	1.73790
#define crho_	2.40061
#define c_	(-0.06904)
#define d1_	(-0.79271)

/* Classical parameters */
#define a05_	0.16838
#define a06_	0.73251
#define a14_	0.55023
#define a22_	1.22624

/* Equation of state background parameters */
#define _A0_	(-1.00000)
#define _A1_	(-6.06433)
#define _A2_	7.71501
#define _A3_	(-3.08501)
#define _A4_	6.70090

/* Caloric background parameters */
#define _mu2_	(-34.21210)
#define _mu3_	0.37334
#define _mu4_	(-59.99810)

/* ************************************************************************ */
double getd_P__d_T(double _T, double _rho, double dDel_A__dDel_T, double ddDel_A__dDel_rhodDel_T);
double getd_P__d_rho(double _rho, double Del_mu, double dDel_A__dDel_rho, double ddDel_A__ddDel_rho);
double get_cv(double _T, double _rho, double ddDel_A__ddDel_T);
double get_cp(double _T, double _P, double _rho, double _xi, double _cv, double d_P__d_T);
double getalfa(double T, double P, double rho, double d_P__d_T, double d_P__d_rho);
double get_W(double _rho, double _xi, double _cv, double _cp);
/* ************************************************************************ */
double gett(double Del_T, double dDel_Ax__dM);
double getM(double Del_T, double Del_rho, double dDel_Ax__dt);
double getDel_A(double Del_Ax, double dDel_Ax__dt, double dDel_Ax__dM);
double getDel_Ax(double t, double M, double Y);
double getdDel_Ax__dt(double t, double M, double Y, double dY__dt);
double getddDel_Ax__ddt(double t, double M, double Y, double dY__dt, double ddY__ddt);
double getdDel_Ax__dM(double t, double M, double Y, double dY__dM);
double getddDel_Ax__ddM(double t, double M, double Y, double dY__dM, double ddY__ddM);
double getddDel_Ax__dtdM(double t, double M, double Y, double dY__dt, double dY__dM, double ddY__dtdM);
double getkapp_2(double t, double M, double Y);
void bisection_Y(double *solution, double lY, double rY, double t, double M);
//double bisection_Y(double lY, double rY, double t, double M);
double getY(double t, double M);
double getdDel_A__dDel_rho(double dDel_Ax__dM);
double getdDel_A__dDel_T(double dDel_Ax__dt, double dDel_Ax__dM);
double getG(double ddDel_Ax__ddt, double ddDel_Ax__ddM, double ddDel_Ax__dtdM);
double getddDel_A__ddDel_rho(double ddDel_Ax__ddt, double ddDel_Ax__ddM, double ddDel_Ax__dtdM);
double getddDel_A__dDel_rhodDel_T(double ddDel_Ax__ddt, double ddDel_Ax__ddM, double ddDel_Ax__dtdM);
double getddDel_A__ddDel_T(double ddDel_Ax__ddt, double ddDel_Ax__ddM, double ddDel_Ax__dtdM);
double getf1(double Y, double kapp_2);
double getf3(double Y, double kapp_2);
double getf5(double t, double M, double Y);
double getf6(double t, double M, double Y);
double getg5(double t, double M, double Y);
double getdY__dt(double f1, double f3, double f5, double f6);
double getdY__dM(double f1, double f3, double g5, double f6);
double getdkapp_2__dt(double f1, double f3, double f5, double f6);
double getdkapp_2__dM(double f1, double f3, double g5, double f6);
double getdf1__dt(double Y, double kapp_2, double dY__dt, double dkapp_2__dt);
double getdf3__dt(double Y, double kapp_2, double dY__dt, double dkapp_2__dt);
double getdf5__dt(double t, double M, double Y, double dY__dt);
double getdf6__dt(double t, double M, double Y, double dY__dt);
double getdg5__dt(double t, double M, double Y, double dY__dt);
double getdf1__dM(double Y, double kapp_2, double dY__dM, double dkapp_2__dM);
double getdf3__dM(double Y, double kapp_2, double dY__dM, double dkapp_2__dM);
double getdg5__dM(double t, double M, double Y, double dY__dM);
double getdf6__dM(double t, double M, double Y, double dY__dM);
double getddY__ddt(double f1, double f3, double f5, double f6, 
			double df1__dt, double df3__dt, double df5__dt, double df6__dt);
double getddY__ddM(double f1, double f3, double g5, double f6,
			double df1__dM, double df3__dM, double dg5__dM, double df6__dM);
double getddY__dtdM(double f1, double f3, double g5, double f6,
			double df1__dt, double df3__dt, double dg5__dt, double df6__dt);

double getlambda_rho(double T, double rho);
double getcp_isochoric(double T);

/* *** for guessing rho *** */
double sf6_pressure(double temp, double rho);
double sf6_rho(double rho1, double rho2, double temp, double press);
double set_endrho_gasregion(double temp, double press, double rho_acc);

/* *** check if temp and press are valid * ***/
int istpOK(double temp, double press);

/* ******************************************************************** */

void sf6_crit(tt, pp, prho, palfa, pcomp, plambda, pkappa, pnu, pcp)
double tt, pp, *prho, *palfa, *pcomp, *plambda, *pkappa, *pnu, *pcp;
{
	double lambda_t, lambda_rho, lambda_trho;	// thermal conductivity
	double lambda_regular, lambda_singular, lambda_isochoric, cp_isochoric;
	double mu_tp;           			// viscosity
	double mu0, mu1, mu2, mu3, mu4, mu5, mu6, mu7;	// viscosity coefficients
	double coexistdelrho_t;		// rho_liquid - rho_gas at the coexistence curve at given t;
	double LHW;
	double xi, comp;			// suseptibility and compressibility

	double tk, T, _T, Del_T;           	// temperature (in Kelvin), reduced, reduced difference
	double rho, _rho, Del_rho;		// density, reduced, reduced difference
	double _rho1;
	double rrd;				// density in amagat
	double ppa, P, _P;			// pressure in Pa, reduced
	double _P1;
	double t;				// temperaturelike variable
	double t1;
	double M;				// densitylike order parameter
	double M1;
	double Del_A, Del_Ax, _A0;
	double dDel_A__dDel_rho, dDel_A__dDel_T;
	double ddDel_A__ddDel_rho, ddDel_A__ddDel_T;
	double ddDel_A__dDel_rhodDel_T;
	double dDel_Ax__dM, dDel_Ax__dt;
	double ddDel_Ax__ddM, ddDel_Ax__ddt;
	double ddDel_Ax__dtdM;
	double Del_mu;
	double Y;				// crossover function
	double kapp_2;				// crossover function
	double dY__dt, dY__dM;
	double dkapp_2__dt, dkapp_2__dM;
	double ddY__ddt, ddY__ddM, ddY__dtdM;
	double f1, f3, f6, f5, g5;
	double df1__dt, df3__dt, df6__dt, df5__dt, dg5__dt;
	double df1__dM, df3__dM, df6__dM, dg5__dM;
	double cv, cp, _cv, _cp;		// heat capacity and the reduced
	double alfa;				// thermal expansion
	double _xi;				// reduced xi
	double W, _W;				// speed of sound and the reduced
	double d_P__d_T, d_P__d_rho;
	int count, cnt;

	double b_rho, rho_acc;

	double _rhoarr[3], _rhoflag, rcorr;
	int cnt2;

if (!istpOK(tt,pp)) {
	printf("temp=%.3lfC and press=%.3lfpsi are not in the valid region.\n",tt,pp);
	exit(-1);
}

//if (tt<=Tc-ZERO_C) { printf("Out of valid temperature range. //temp=%.3lf<=%.3lf.\n",tt,Tc-ZERO_C); exit(-1); }
//if (tt>=390.-ZERO_C) { printf("Out of valid temperature range. //temp=%.3lf>=%.3lf\n",tt,390.-ZERO_C); exit(-1); }

	pp = pp/14.504;		/* convert pressure from psi to bar */ 
	tk = tt + ZERO_C;	/*  convert C into Kelvin */
	ppa = pp*1.0e+5;	/*  convert bar into Pa=N/M^2  */

	T=tk; _T=-Tc/T; Del_T=_T+1.;
	P=ppa; _P=Tc/Pc*P/T;
	_A0=_A0_+_A1_*Del_T+_A2_*Del_T*Del_T+_A3_*pow(Del_T,3.)+_A4_*pow(Del_T,4.);

/* *** guess rho using the equation of state from T. E. Morsy, J. Chem. and Engin. Data, v.15, 256 (1970) */

	rho_acc=40.;
	b_rho=set_endrho_gasregion(T,P,rho_acc);
	rho=sf6_rho(b_rho-rho_acc, b_rho, T, P);
	if (rho<336.) rho=336.;
	if (rho>1170.) rho=1170.;
	_rho=rho/RHOc;
cnt=0;
rcorr=20.;
for (;;) { // iteration for rho
	Del_rho=_rho-1.;

	t=gett(Del_T,0.);//ct_*Del_T;		//guess t
	M=getM(Del_T,Del_rho,0.);//crho_*(Del_rho-d1_*Del_T);	//guess M
    count=0;
    for (;;) { // iteration for t and M
	count++;
	Y=getY(t,M);
	kapp_2=getkapp_2(t,M,Y);
	f1=getf1(Y,kapp_2);
	f3=getf3(Y,kapp_2);
	f5=getf5(t,M,Y);
	f6=getf6(t,M,Y);
	g5=getg5(t,M,Y);
	dY__dt=getdY__dt(f1,f3,f5,f6);
	dY__dM=getdY__dM(f1,f3,g5,f6);
	dDel_Ax__dM=getdDel_Ax__dM(t,M,Y,dY__dM);
	t1=gett(Del_T,dDel_Ax__dM);
	dDel_Ax__dt=getdDel_Ax__dt(t,M,Y,dY__dt);
	M1=getM(Del_T,Del_rho,dDel_Ax__dt);
	if (fabs(t1-t)<1.0e-6*fabs(t) && fabs(M1-M)<1.0e-6*fabs(M)) {
		break;
	}
	else if (count>200) {
		printf("count>200\n");
		break;
	}
	else {
		t+=(t1-t)/10.; M+=(M1-M)/10.;
	}
    }

	Y=getY(t,M);
	kapp_2=getkapp_2(t,M,Y);
	f1=getf1(Y,kapp_2);
	f3=getf3(Y,kapp_2);
	f5=getf5(t,M,Y);
	f6=getf6(t,M,Y);
	g5=getg5(t,M,Y);
	dY__dt=getdY__dt(f1,f3,f5,f6);
	dY__dM=getdY__dM(f1,f3,g5,f6);
	dDel_Ax__dM=getdDel_Ax__dM(t,M,Y,dY__dM);
	t1=gett(Del_T,dDel_Ax__dM);
	dDel_Ax__dt=getdDel_Ax__dt(t,M,Y,dY__dt);
	M1=getM(Del_T,Del_rho,dDel_Ax__dt);
	if (fabs(t1-t)>1.0e-5*fabs(t) || fabs(M1-M)>1.0e-5*fabs(M)) {
		printf("More iterations are needed.\n");
		printf("t1/t-1=%lf M1/M-1=%lf\n",fabs(t1/t-1.),fabs(M1/M-1.));
		exit(-1);
	}
	Del_Ax=getDel_Ax(t,M,Y);
	dDel_Ax__dt=getdDel_Ax__dt(t,M,Y,dY__dt);
	dDel_Ax__dM=getdDel_Ax__dM(t,M,Y,dY__dM);
	Del_A=getDel_A(Del_Ax,dDel_Ax__dt,dDel_Ax__dM);
	Del_mu=getdDel_A__dDel_rho(dDel_Ax__dM);
	_P1=_rho*Del_mu-Del_A-_A0;
	_rho1=(_P+Del_A+_A0)/Del_mu;
cnt++;
	if (fabs(_rho1-_rho)<1.0e-4*fabs(_rho)) { break; }
	else {
	    if (cnt>200) { 
		printf("\nRho-seeking-loop stops at 200th iteration!\n");
		break;
	    }

	/* for fast approach in the iteration (turned out to be very effective) */
	    cnt2=(cnt-1)%3;
	    _rhoarr[cnt2]=_rho;
	    if (cnt>3) {
		if (cnt2==2) _rhoflag=(_rhoarr[2]-_rhoarr[1])*(_rhoarr[1]-_rhoarr[0]);
		else if (cnt2==1) _rhoflag=(_rhoarr[1]-_rhoarr[0])*(_rhoarr[0]-_rhoarr[2]);
		else if (cnt2==0) _rhoflag=(_rhoarr[0]-_rhoarr[2])*(_rhoarr[2]-_rhoarr[1]);
		else { printf("Error in cnt2!!!\n"); exit(0); }
	    }
	    if (_rhoflag<0.) { rcorr*=2.; }
	    if (_P1>_P) { _rho-=fabs((_rho1-_rho)/rcorr); }
	    else {	_rho+=fabs((_rho1-_rho)/rcorr); }
	    if (_rho<336./RHOc) { _rho=336./RHOc; }
	    if (_rho>1170./RHOc) { _rho=1170./RHOc; }
	}
//	printf("iteration %d:",cnt);
//	printf("_P=%lf _P1=%lf _rho1=%lf _rho=%lf\n",_P,_P1,_rho1,_rho);
}
	rho=_rho*RHOc;
	_P1=_rho*Del_mu-Del_A-_A0;

	dkapp_2__dt=getdkapp_2__dt(f1,f3,f5,f6);
	dkapp_2__dM=getdkapp_2__dM(f1,f3,g5,f6);
	df1__dt=getdf1__dt(Y,kapp_2,dY__dt,dkapp_2__dt);
	df3__dt=getdf3__dt(Y,kapp_2,dY__dt,dkapp_2__dt);
	df5__dt=getdf5__dt(t,M,Y,dY__dt);
	df6__dt=getdf6__dt(t,M,Y,dY__dt);
	dg5__dt=getdg5__dt(t,M,Y,dY__dt);
	df1__dM=getdf1__dM(Y,kapp_2,dY__dM,dkapp_2__dM);
	df3__dM=getdf3__dM(Y,kapp_2,dY__dM,dkapp_2__dM);
	dg5__dM=getdg5__dM(t,M,Y,dY__dM);
	df6__dM=getdf6__dM(t,M,Y,dY__dM);
	ddY__ddt=getddY__ddt(f1,f3,f5,f6,df1__dt,df3__dt,df5__dt,df6__dt);
	ddY__ddM=getddY__ddM(f1,f3,g5,f6,df1__dM,df3__dM,dg5__dM,df6__dM);
	ddY__dtdM=getddY__dtdM(f1,f3,g5,f6,df1__dt,df3__dt,dg5__dt,df6__dt);
	Del_Ax=getDel_Ax(t,M,Y);
	dDel_Ax__dt=getdDel_Ax__dt(t,M,Y,dY__dt);
	ddDel_Ax__ddt=getddDel_Ax__ddt(t,M,Y,dY__dt,ddY__ddt);
	dDel_Ax__dM=getdDel_Ax__dM(t,M,Y,dY__dM);
	ddDel_Ax__ddM=getddDel_Ax__ddM(t,M,Y,dY__dM,ddY__ddM);
	ddDel_Ax__dtdM=getddDel_Ax__dtdM(t,M,Y,dY__dt,dY__dM,ddY__dtdM);
	Del_A=getDel_A(Del_Ax,dDel_Ax__dt,dDel_Ax__dM);
	Del_mu=getdDel_A__dDel_rho(dDel_Ax__dM);

	dDel_A__dDel_T=getdDel_A__dDel_T(dDel_Ax__dt,dDel_Ax__dM);
	ddDel_A__ddDel_T=getddDel_A__ddDel_T(ddDel_Ax__ddt,ddDel_Ax__ddM,ddDel_Ax__dtdM);
	ddDel_A__dDel_rhodDel_T=getddDel_A__dDel_rhodDel_T(ddDel_Ax__ddt,ddDel_Ax__ddM,ddDel_Ax__dtdM);
	d_P__d_T=getd_P__d_T(_T,_rho,dDel_A__dDel_T,ddDel_A__dDel_rhodDel_T);
	dDel_A__dDel_rho=getdDel_A__dDel_rho(dDel_Ax__dM);
	ddDel_A__ddDel_rho=getddDel_A__ddDel_rho(ddDel_Ax__ddt,ddDel_Ax__ddM,ddDel_Ax__dtdM);
	d_P__d_rho=getd_P__d_rho(_rho,Del_mu,dDel_A__dDel_rho,ddDel_A__ddDel_rho);
	_xi=1./ddDel_A__ddDel_rho;
	_cv=get_cv(_T,_rho,ddDel_A__ddDel_T);
	_cp=get_cp(_T,_P,_rho,_xi,_cv,d_P__d_T);
	alfa=-getalfa(T,P,rho,d_P__d_T,d_P__d_rho); //put '-' for positive number
	_W=get_W(_rho,_xi,_cv,_cp);

	cv=_cv*Pc/Tc/rho/1000.; /* the isochoric heat capacity !! divide by 1000 added 12-1-02 GA */
	cp=_cp*Pc/Tc/rho/1000.; /* the isobaric heat capacity !! divide by 1000 added 12-1-02 GA*/
	W=_W/sqrt(Tc*RHOc/Pc/T); /* the speed of sound, m/sec */
	xi=_xi*RHOc*RHOc*Tc/Pc/T; /* xi is inverse of the susceptibility */
	comp=xi/rho/rho; /* compressibility */

//printf("T=%.2lf K, rho=%.2lf kg/m^3, P=%.3lf MPa\n",T,_rho*RHOc,_P1*Pc/Tc*T/1.0e6);
//printf("cv=%.3lf kJ/kgK, cp=%.3lf kJ/kgK, W=%.2lf m/s\n",cv, cp, W);
//printf("_xi,alfa,dP__dT=%.3lf, %.3lf, %.3lf\n",_xi,alfa,P/T+Pc/T/T*d_P__d_T);
//printf("_rho,_rho1=%lf,%lf: t,t1=%lf,%lf M,M1=%lf,%lf dDel_Ax__dM=%lf\n",_rho,_rho1,t,t1,M,M1,dDel_Ax__dM);
//printf("rho=%lf, T=%lf, P=%lf, Pcal=%lf\n",_rho*RHOc,T,P,(_rho*Del_mu-Del_A-_A0)*Pc/Tc*T);

if (1./_xi>=2.0 || rho<=336. || rho>1170.) {
	printf("Out of range where crossover functions are valid.\n");
	printf("rho=%lf kg/m^3, (1/_xi)=%.5lf\n",rho,1./_xi);
	printf("*** Valid region: 1./_xi <2.0\n");
	printf("              and 336kg/m^3 < rho < 1170kg/m^3\n");
	printf("              and 318.717K < T < 390K at rho=RHOc ***\n");
	printf(" see A.K.Wyczalkowska & J.V.Sengers\n");
	exit(-1);
}

/* From J. H. B. Hoogland, H. R. Van den Berg, and N. J. Trappeniers, Physica A, V.134, 169 (1985) and
	T. Strehlow and E. Vogel, Physica A, V.161 (1989)
  See Note.1 in the end of this file for details.
*/

	mu0=15.262*exp(0.548599*log(T/298.15)-0.651683/T*298.15+0.158578/T/T*298.15*298.15+0.491296); // \mPa sec
	mu1=6.6159*exp(30.9709*log(T)+2.794e4/T-3.318e6/T/T-240.06);	// \mPa sec/amagat
//	Ref. [182]	=> 
	mu2=2.85099e-3;		// \mPa sec/amagat^2
	mu3=-4.16848e-5;	// \mPa sec/amagat^3
	mu4=6.446345e-7;	// \mPa sec/amagat^4
	mu5=-4.997041e-9;	// \mPa sec/amagat^4
	mu6=1.970048e-11;	// \mPa sec/amagat^6
	mu7=-2.961807e-14;	// \mPa sec/amagat^7

	rrd = rho/6.6159; // density in unit of amagat (1 amagat = 6.6159 kg/m^3)
	mu_tp=mu0+mu1*rrd+mu2*pow(rrd,2.)+mu3*pow(rrd,3.)+mu4*pow(rrd,4.)+mu5*pow(rrd,5.)+mu6*pow(rrd,6.)+mu7*pow(rrd,7.); // mmPa s
	mu_tp = mu_tp*1.0e-6;	    /*  times 1.0e-6 to get the right unit */
			
/* *** for thermal conductivity *** */
/* ***	lambda_rho = lambda(rho,T) - lambda(rho=0,T) :
		eq.(5) in H. L. Swinney and D. L. Henry, Phys. Rev. A, Vol.8, p.2586 (1973)

	lambda(rho=0,tt) = 0.01303*(1.0 + 0.00549*(tt - 27.5)) : property at low density
		 From J. Kestin and N. Imaishi, Int. J. Thermophys. V.6, 107 (1985)
		 Valid 0<tt<100 C, tt is temperature in Celcius.

	The following data in double getlambda_rho(T,rho) are from 
		J. Lis and P. O. Kellard, Brit. J. Appl. Phys., Vol.16, p.1099 (1965)
		V. V. Burinskii and E. E. Totskii, High Temp., Vol.19, p.366 (1981)
*** */

/* lambda_rho = lambda(rho,T) - lambda(rho=0,T) :
		eq.(5) in H. L. Swinney and D. L. Henry, Phys. Rev. A, Vol.8, p.2586 (1973)
*/
	lambda_rho=getlambda_rho(T,rho); // in W/meter/K

/* lambda_t=lambda(rho=0,T) : thermal conductivity at low density
	from J. Kestin and N. Imaishi, Int. J. Thermophys. V.6, 107 (1985)
	valid 0<tt<100 C
*/
	lambda_t = 0.01303*(1.0 + 0.00549*(tt - 27.5));

/* the regular part of thermal conductivity, eq.(6)  in H. L. Swinney and D. L. Henry, Phys. Rev. A, Vol.8, p.2586 (1973) */
	lambda_regular=lambda_rho+lambda_t;

/* kappa_isochoric = thermal diffusivity at the critical isochore with background correction
	kappa_isochoric=2.4*(T/Tc-1.)^0.62*1.0e-4 cm^2/sec :
		from T. K. Lim, PhD thesis, Johns Hopkins University (1973)
*/
	cp_isochoric=getcp_isochoric(T);

/* singular part of lambda at the critical isochore */
	lambda_isochoric=2.4*pow(T/Tc-1.,0.62)*1.0e-8*cp_isochoric*RHOc;

/* (inverted) coexistence curve, E. S. Wu and W. W. Webb, J. Phys. (Paris), Vol.33, C1~149 (1972) */
	coexistdelrho_t=3.62*pow(T/Tc-1.,0.333)*RHOc;

/* assume that the critical isochoric (singular) thermal conductivity follow Lorentzian distribution w.r.t. rho at fixed temperature. */
	LHW=coexistdelrho_t/2.; // set the half width at half maximum of Lorentzian profile on the (inverted) coexistence curve

		lambda_singular=lambda_isochoric/(1.+(rho-RHOc)*(rho-RHOc)/LHW/LHW); // 10-2-2000

/* total thermal conductivity by summing the regular part and the singular part */
	lambda_trho=lambda_regular+lambda_singular;

//printf("lambda_regular,singular,LHW=%lf,%lf,%lf\n",lambda_regular,lambda_singular,LHW);

//printf("lambda_: rho,t,regular,isochoric,singular,trho=%lf,%lf,%lf,%lf,%lf,%lf\n",lambda_rho,lambda_t,lambda_regular,lambda_isochoric,lambda_singular,lambda_trho);
//printf("cp_isochoric=%lf\n\n",cp_isochoric);

/* convert to cgs and assign values */

	*plambda = 1.e5*lambda_trho;
	*pnu = 10.*mu_tp;	
	*palfa = alfa;
	*prho = 0.001*rho;//rr;
	*pcp = 1.e7*cp; // SF6_MW;	changed 12-2-02 GA
	*pcomp = comp/10.; //cm^2/dyne
	*pkappa = *plambda/((*prho)*(*pcp));

/* Don't forget we have also calculated the isochoric heat capacity (cv) and the speed of sound (W) which are not returned to noise.c here. */

}

/*** check if temp and press are in the valid region ***/
int istpOK(double temp, double press) {
	double y1, y2, y3, y4, y5;
	double x;

	x=press;
	y1=-68.439+0.2333*x;
	y2=22.2424+0.0339*x;
	y3=45.568;
	y4=-162.19+0.7834*x-0.001011*x*x+(6.613e-7)*x*x*x-(1.697e-10)*x*x*x*x;
	y5=-1093.3+2.5449*x-0.0018*x*x+(3.368e-7)*x*x*x+(5.589e-11)*x*x*x*x;
	if (temp<y1 && temp>y2 && temp>=y3 && temp<y4 && temp>y5 && press<1539.) {
		return(1);
	}
	else {
		return(0);
	}
}

/* *** for guessing rho *** */
double sf6_pressure(double temp, double rho)
/* equation of state from T. E. Morsy, J. Chem. and Engin. Data, v.15, 256 (1970) */
{
	double tr, rhor, pr;
	double alphac, oneoverZc, a1, a2 ,a3, a4, a5, a6, a7, a8, a9, a10, beta;

	alphac=6.651;
	oneoverZc=3.4998179;
	a1=2.7175825;
	a2=-1.1209876;
	a3=-7.9533437;
	a4=2.0935047;
	a5=0.238431307;
	a6=0.350845138;
	a7=-1.1128842;
	a8=1.4530336;
	a9=0.896589421;
	a10=0.174569983;
	beta=0.5348;

	tr=temp/Tc;
	rhor=rho/RHOc;
	pr=tr*rhor*oneoverZc+(a1+a2*tr+a3/tr/tr+a4/tr/tr/tr/tr)*rhor*rhor+(a5+a6*tr+a7/tr/tr)*rhor*rhor*rhor+(a8+a9/tr/tr)*exp(-beta*rhor*rhor)/tr/tr*(1+beta*rhor*rhor)*rhor*rhor*rhor+a10*rhor*rhor*rhor*rhor*rhor*rhor;
	return(Pc*pr);
}

double sf6_rho(double rho1, double rho2, double temp, double press)
{
	double p1, p2, c;

	c=(rho1+rho2)/2.;
	if (rho2-c<1.0e-10) {return(c);}
	else {
		p1=sf6_pressure(temp,c);
		p2=sf6_pressure(temp,rho2);
//printf("%.10lf\n",rho2-c);
		if ( (p1-press)*(p2-press)<0. ) {rho1=c;}
		else {rho2=c;}
		sf6_rho(rho1,rho2,temp,press);
	}
}

double set_endrho_gasregion(double temp, double press, double rho_acc)
{
	double trho;

	for (trho=0.; trho<15390.; trho+=rho_acc) {
		if ( (0.-press)*(sf6_pressure(temp,trho)-press)<0. ) return(trho);
	}
	printf("Error in the searching rho range!\n");
	exit(-1);
}

/* ************ for thermal conductivity ************ */

double getlambda_rho(double T, double rho) {

	int i,j;
	double rhos[18], lambda_rhos[18], lambda_rho;

/* ***	lambda_rho = lambda(rho,T) - lambda(rho=0,T) :
		eq.(5) in H. L. Swinney and D. L. Henry, Phys. Rev. A, Vol.8, p.2586 (1973)

	lambda(rho=0,tt) = 0.01303*(1.0 + 0.00549*(tt - 27.5)) : property at low density
		 From J. Kestin and N. Imaishi, Int. J. Thermophys. V.6, 107 (1985)
		 Valid 0<tt<100 C, tt is temperature in Celcius.

	The following data are from 
		J. Lis and P. O. Kellard, Brit. J. Appl. Phys., Vol.16, p.1099 (1965)
		V. V. Burinskii and E. E. Totskii, High Temp., Vol.19, p.366 (1981)
*** */

	rhos[0]=1259.558914; lambda_rhos[0]=37.505630; // at T=335.1 K, rho in kg/m^3, lambda_rho in mW/m/K
	rhos[1]=1103.507220; lambda_rhos[1]=32.926088; // at T=324.05 K
	rhos[2]=1030.553735; lambda_rhos[2]=32.276088; // at T=324.05 K
	rhos[3]=923.565772; lambda_rhos[3]=29.029726; // at T=330.15 K
	rhos[4]=837.711879; lambda_rhos[4]=27.929726; // at T=330.15 K
	rhos[5]=800.294483; lambda_rhos[5]=27.309726; // at T=330.15 K
	rhos[6]=770.910691; lambda_rhos[6]=26.759726; // at T=330.15 K
	rhos[7]=714.968208; lambda_rhos[7]=26.649726; // at T=330.15 K
	rhos[8]=685.788605; lambda_rhos[8]=26.069726; // at T=330.15 K
	rhos[9]=647.035399; lambda_rhos[9]=26.719726; // at T=330.15 K
	rhos[10]=615.768265; lambda_rhos[10]=26.119726; // at T=330.15 K
	rhos[11]=611.500579; lambda_rhos[11]=26.209726; // at T=330.15 K
	rhos[12]=590.959511; lambda_rhos[12]=22.049726; // at T=330.15 K
	rhos[13]=568.071442; lambda_rhos[13]=20.809726; // at T=330.15 K
	rhos[14]=482.397293; lambda_rhos[14]=14.969726; // at T=330.15 K
	rhos[15]=465.323008; lambda_rhos[15]=13.409726; // at T=330.15 K
	rhos[16]=393.532482; lambda_rhos[16]=10.009726; // at T=330.15 K
	rhos[17]=344.985050; lambda_rhos[17]=8.039726; // at T=330.15 K

	j=-1;
	for (i=0; i<17; i++) {
		if (rho<=rhos[i] && rho>=rhos[i+1]) {
			j=i;
			break;
		}
	}
	if (j==-1) {
		printf("rho is not in the valid range.\n");
		exit(-1);
	}
	lambda_rho=(lambda_rhos[j]-lambda_rhos[j+1])/(rhos[j]-rhos[j+1])*(rho-rhos[j+1])+lambda_rhos[j+1];
	return(lambda_rho/1000.); // in W/m/K
}

/* ************************************************************************ */

double getd_P__d_T(double _T, double _rho, double dDel_A__dDel_T, double ddDel_A__dDel_rhodDel_T) {
	double d_P__d_T, d_A0__d_T;

	d_A0__d_T=_A1_+2.*_A2_*(_T+1.)+3.*_A3_*(_T+1.)*(_T+1.)+4.*_A4_*(_T+1.)*(_T+1.)*(_T+1.);
	d_P__d_T=-d_A0__d_T+_rho*ddDel_A__dDel_rhodDel_T-dDel_A__dDel_T;
	return(d_P__d_T);
}

double getd_P__d_rho(double _rho, double Del_mu, double dDel_A__dDel_rho, double ddDel_A__ddDel_rho) {
	double d_P__d_rho;

	d_P__d_rho=Del_mu+_rho*ddDel_A__ddDel_rho-dDel_A__dDel_rho;
	return(d_P__d_rho);
}
	
double get_cv(double _T, double _rho, double ddDel_A__ddDel_T) {
	double dd_A0__dd_T, dd_mu0__dd_T, _cv;

	dd_A0__dd_T=2.*_A2_+6.*_A3_*(_T+1.)+12.*_A4_*(_T+1.)*(_T+1.);
	dd_mu0__dd_T=2.*_mu2_+6.*_mu3_*(_T+1.)+12.*_mu4_*(_T+1.)*(_T+1.);
	_cv=(-dd_A0__dd_T-_rho*dd_mu0__dd_T-ddDel_A__ddDel_T)*_T*_T;
	return(_cv);
}

double get_cp(double _T, double _P, double _rho, double _xi, double _cv, double d_P__d_T) {
	double _cp;

	_cp=_cv+_xi/_rho/_rho*pow(_P-_T*d_P__d_T,2.);
	return(_cp);
}

double getalfa(double T, double P, double rho, double d_P__d_T, double d_P__d_rho) {
	double alfa;

//	alfa=-1./rho*(P+Pc/T*d_P__d_T)/(Pc/Tc*T*T/RHOc*d_P__d_rho);
        alfa=-1./rho*(P+Pc*d_P__d_T)/(Pc/Tc*T*T/RHOc*d_P__d_rho);	/* changed 12-02-02 GA as instructed by JO */
        return(alfa);
}


double get_W(double _rho, double _xi, double _cv, double _cp) {
	return(sqrt(_rho/_xi*_cp/_cv));
}

/* ************************************************************************ */

double gett(double Del_T, double dDel_Ax__dM) {
	double t;

	t=ct_*Del_T+c_*dDel_Ax__dM;
	return(t);
}

double getM(double Del_T, double Del_rho, double dDel_Ax__dt) {
	double M;

	M=crho_*(Del_rho-d1_*Del_T)+c_*dDel_Ax__dt;
	return(M);
}

double getDel_A(double Del_Ax, double dDel_Ax__dt, double dDel_Ax__dM) {
	double Del_A;

	Del_A=Del_Ax-c_*dDel_Ax__dM*dDel_Ax__dt;
	return(Del_A);
}

double getDel_Ax(double t, double M, double Y) {
	double Del_Ax;
	double c1, c2, c3, c4, c5, c6, c7;
	double e1, e2, e3, e4, e5, e6, e7;

	c1=0.5;
	e1=(2.*nu_-1.-eta_*nu_)/DelS_;
	c2=uaster_*ubar_*GAMMA_/24.;
	e2=(nu_-2.*eta_*nu_)/DelS_;
	c3=a05_/120.;
	e3=(DelA_+nu_/2.-2.5*eta_*nu_)/DelS_;
	c4=a06_/720.;
	e4=(1.5*nu_-3.*eta_*nu_)/DelS_;
	c5=a14_/24.;
	e5=(1.5*nu_-1.-2.*eta_*nu_)/DelS_;
	c6=a22_/4.;
	e6=(3.5*nu_-2.-eta_*nu_)/DelS_;
	c7=nu_/alpha_/ubar_/GAMMA_/2.;
	e7=-alpha_/DelS_;

//	Del_Ax=0.5*t*M*M*pow(Y,(2.*nu_-1.-eta_*nu_)/DelS_);
//	Del_Ax+=uaster_*ubar_*GAMMA_/24.*M*M*M*M*pow(Y,(-2.*eta_*nu_+nu_)/DelS_);
//	Del_Ax+=a05_/120.*M*M*M*M*M*pow(Y,(-2.5*eta_*nu_+DelA_+nu_/2.)/DelS_);
//	Del_Ax+=a06_/720.*M*M*M*M*M*M*pow(Y,(-3.*eta_*nu_+1.5*nu_)/DelS_);
//	Del_Ax+=a14_/24.*t*M*M*M*M*pow(Y,(1.5*nu_-1.-2.*eta_*nu_)/DelS_);
//	Del_Ax+=a22_/4.*t*t*M*M*pow(Y,(3.5*nu_-2.-eta_*nu_)/DelS_);
//	Del_Ax-=t*t/2.*nu_/alpha_/ubar_/GAMMA_*(pow(Y,-alpha_/DelS_)-1.);
	Del_Ax=c1*t*M*M*pow(Y,e1)+c2*M*M*M*M*pow(Y,e2)+c3*M*M*M*M*M*pow(Y,e3)+c4*M*M*M*M*M*M*pow(Y,e4)+c5*t*M*M*M*M*pow(Y,e5)+c6*t*t*M*M*pow(Y,e6)-c7*t*t*(pow(Y,e7)-1.);
	return(Del_Ax);
}

double getdDel_Ax__dt(double t, double M, double Y, double dY__dt) {
	double dDel_Ax__dt;
	double c1, c2, c3, c4, c5, c6, c7;
	double e1, e2, e3, e4, e5, e6, e7;

	c1=0.5;
	e1=(2.*nu_-1.-eta_*nu_)/DelS_;
	c2=uaster_*ubar_*GAMMA_/24.;
	e2=(nu_-2.*eta_*nu_)/DelS_;
	c3=a05_/120.;
	e3=(DelA_+nu_/2.-2.5*eta_*nu_)/DelS_;
	c4=a06_/720.;
	e4=(1.5*nu_-3.*eta_*nu_)/DelS_;
	c5=a14_/24.;
	e5=(1.5*nu_-1.-2.*eta_*nu_)/DelS_;
	c6=a22_/4.;
	e6=(3.5*nu_-2.-eta_*nu_)/DelS_;
	c7=nu_/alpha_/ubar_/GAMMA_/2.;
	e7=-alpha_/DelS_;

	dDel_Ax__dt=c1*M*M*(pow(Y,e1)+t*e1*pow(Y,e1-1.)*dY__dt);
	dDel_Ax__dt+=c2*M*M*M*M*e2*pow(Y,e2-1.)*dY__dt;
	dDel_Ax__dt+=c3*M*M*M*M*M*e3*pow(Y,e3-1.)*dY__dt;
	dDel_Ax__dt+=c4*M*M*M*M*M*M*e4*pow(Y,e4-1.)*dY__dt;
	dDel_Ax__dt+=c5*M*M*M*M*(pow(Y,e5)+t*e5*pow(Y,e5-1.)*dY__dt);
	dDel_Ax__dt+=c6*M*M*(2.*t*pow(Y,e6)+t*t*e6*pow(Y,e6-1.)*dY__dt);
	dDel_Ax__dt-=c7*(2.*t*(pow(Y,e7)-1.)+t*t*e7*pow(Y,e7-1.)*dY__dt);
	return(dDel_Ax__dt);
}

double getddDel_Ax__ddt(double t, double M, double Y, double dY__dt, double ddY__ddt) {
	double ddDel_Ax__ddt;
	double c1, c2, c3, c4, c5, c6, c7;
	double e1, e2, e3, e4, e5, e6, e7;

	c1=0.5;
	e1=(2.*nu_-1.-eta_*nu_)/DelS_;
	c2=uaster_*ubar_*GAMMA_/24.;
	e2=(nu_-2.*eta_*nu_)/DelS_;
	c3=a05_/120.;
	e3=(DelA_+nu_/2.-2.5*eta_*nu_)/DelS_;
	c4=a06_/720.;
	e4=(1.5*nu_-3.*eta_*nu_)/DelS_;
	c5=a14_/24.;
	e5=(1.5*nu_-1.-2.*eta_*nu_)/DelS_;
	c6=a22_/4.;
	e6=(3.5*nu_-2.-eta_*nu_)/DelS_;
	c7=nu_/alpha_/ubar_/GAMMA_/2.;
	e7=-alpha_/DelS_;

	ddDel_Ax__ddt=c1*M*M*(2.*e1*pow(Y,e1-1.)*dY__dt+t*e1*(e1-1.)*pow(Y,e1-2.)*dY__dt*dY__dt+t*e1*pow(Y,e1-1.)*ddY__ddt);
	ddDel_Ax__ddt+=c2*M*M*M*M*e2*((e2-1.)*pow(Y,e2-2.)*dY__dt*dY__dt+pow(Y,e2-1.)*ddY__ddt);
	ddDel_Ax__ddt+=c3*M*M*M*M*M*e3*((e3-1.)*pow(Y,e3-2.)*dY__dt*dY__dt+pow(Y,e3-1.)*ddY__ddt);
	ddDel_Ax__ddt+=c4*M*M*M*M*M*M*e4*((e4-1.)*pow(Y,e4-2.)*dY__dt*dY__dt+pow(Y,e4-1.)*ddY__ddt);
	ddDel_Ax__ddt+=c5*M*M*M*M*e5*(2.*pow(Y,e5-1.)*dY__dt+t*(e5-1.)*pow(Y,e5-2.)*dY__dt*dY__dt+t*pow(Y,e5-1.)*ddY__ddt);
	ddDel_Ax__ddt+=c6*M*M*(2.*pow(Y,e6)+4.*t*e6*pow(Y,e6-1.)*dY__dt+t*t*e6*(e6-1.)*pow(Y,e6-2.)*dY__dt*dY__dt+t*t*e6*pow(Y,e6-1.)*ddY__ddt);
	ddDel_Ax__ddt-=c7*(2.*(pow(Y,e7)-1.)+4.*t*e7*pow(Y,e7-1.)*dY__dt+t*t*e7*(e7-1.)*pow(Y,e7-2.)*dY__dt*dY__dt+t*t*e7*pow(Y,e7-1.)*ddY__ddt);
	return(ddDel_Ax__ddt);
}

double getdDel_Ax__dM(double t, double M, double Y, double dY__dM) {
	double dDel_Ax__dM;
	double c1, c2, c3, c4, c5, c6, c7;
	double e1, e2, e3, e4, e5, e6, e7;

	c1=0.5;
	e1=(2.*nu_-1.-eta_*nu_)/DelS_;
	c2=uaster_*ubar_*GAMMA_/24.;
	e2=(nu_-2.*eta_*nu_)/DelS_;
	c3=a05_/120.;
	e3=(DelA_+nu_/2.-2.5*eta_*nu_)/DelS_;
	c4=a06_/720.;
	e4=(1.5*nu_-3.*eta_*nu_)/DelS_;
	c5=a14_/24.;
	e5=(1.5*nu_-1.-2.*eta_*nu_)/DelS_;
	c6=a22_/4.;
	e6=(3.5*nu_-2.-eta_*nu_)/DelS_;
	c7=nu_/alpha_/ubar_/GAMMA_/2.;
	e7=-alpha_/DelS_;

	dDel_Ax__dM=c1*t*M*(2.*pow(Y,e1)+M*e1*pow(Y,e1-1.)*dY__dM);
	dDel_Ax__dM+=c2*M*M*M*(4.*pow(Y,e2)+M*e2*pow(Y,e2-1.)*dY__dM);
	dDel_Ax__dM+=c3*M*M*M*M*(5.*pow(Y,e3)+M*e3*pow(Y,e3-1.)*dY__dM);
	dDel_Ax__dM+=c4*M*M*M*M*M*(6.*pow(Y,e4)+M*e4*pow(Y,e4-1.)*dY__dM);
	dDel_Ax__dM+=c5*t*M*M*M*(4.*pow(Y,e5)+M*e5*pow(Y,e5-1.)*dY__dM);
	dDel_Ax__dM+=c6*t*t*M*(2.*pow(Y,e6)+M*e6*pow(Y,e6-1.)*dY__dM);
	dDel_Ax__dM-=c7*t*t*e7*pow(Y,e7-1.)*dY__dM;
	return(dDel_Ax__dM);
}

double getddDel_Ax__ddM(double t, double M, double Y, double dY__dM, double ddY__ddM) {
	double ddDel_Ax__ddM;
	double c1, c2, c3, c4, c5, c6, c7;
	double e1, e2, e3, e4, e5, e6, e7;

	c1=0.5;
	e1=(2.*nu_-1.-eta_*nu_)/DelS_;
	c2=uaster_*ubar_*GAMMA_/24.;
	e2=(nu_-2.*eta_*nu_)/DelS_;
	c3=a05_/120.;
	e3=(DelA_+nu_/2.-2.5*eta_*nu_)/DelS_;
	c4=a06_/720.;
	e4=(1.5*nu_-3.*eta_*nu_)/DelS_;
	c5=a14_/24.;
	e5=(1.5*nu_-1.-2.*eta_*nu_)/DelS_;
	c6=a22_/4.;
	e6=(3.5*nu_-2.-eta_*nu_)/DelS_;
	c7=nu_/alpha_/ubar_/GAMMA_/2.;
	e7=-alpha_/DelS_;

	ddDel_Ax__ddM=c1*t*(2.*pow(Y,e1)+4.*M*e1*pow(Y,e1-1.)*dY__dM+M*M*e1*(e1-1.)*pow(Y,e1-2.)*dY__dM*dY__dM+M*M*e1*pow(Y,e1-1.)*ddY__ddM);
	ddDel_Ax__ddM+=c2*M*M*(12.*pow(Y,e2)+8.*M*e2*pow(Y,e2-1.)*dY__dM+M*M*e2*(e2-1.)*pow(Y,e2-2.)*dY__dM*dY__dM+M*M*e2*pow(Y,e2-1.)*ddY__ddM);
	ddDel_Ax__ddM+=c3*M*M*M*(20.*pow(Y,e3)+10.*M*e3*pow(Y,e3-1.)*dY__dM+M*M*e3*(e3-1.)*pow(Y,e3-2.)*dY__dM*dY__dM+M*M*e3*pow(Y,e3-1.)*ddY__ddM);
	ddDel_Ax__ddM+=c4*M*M*M*M*(30.*pow(Y,e4)+12.*M*e4*pow(Y,e4-1.)*dY__dM+M*M*e4*(e4-1.)*pow(Y,e4-2.)*dY__dM*dY__dM+M*M*e4*pow(Y,e4-1.)*ddY__ddM);
	ddDel_Ax__ddM+=c5*t*M*M*(12.*pow(Y,e5)+8.*M*e5*pow(Y,e5-1.)*dY__dM+M*M*e5*(e5-1.)*pow(Y,e5-2.)*dY__dM*dY__dM+M*M*e5*pow(Y,e5-1.)*ddY__ddM);
	ddDel_Ax__ddM+=c6*t*t*(2.*pow(Y,e6)+4.*M*e6*pow(Y,e6-1.)*dY__dM+M*M*e6*(e6-1.)*pow(Y,e6-2.)*dY__dM*dY__dM+M*M*e6*pow(Y,e6-1.)*ddY__ddM);
	ddDel_Ax__ddM-=c7*t*t*(e7*(e7-1.)*pow(Y,e7-2.)*dY__dM*dY__dM+e7*pow(Y,e7-1.)*ddY__ddM);
	return(ddDel_Ax__ddM);
}

double getddDel_Ax__dtdM(double t, double M, double Y, double dY__dt, double dY__dM, double ddY__dtdM) {
	double ddDel_Ax__dtdM;
	double c1, c2, c3, c4, c5, c6, c7;
	double e1, e2, e3, e4, e5, e6, e7;

	c1=0.5;
	e1=(2.*nu_-1.-eta_*nu_)/DelS_;
	c2=uaster_*ubar_*GAMMA_/24.;
	e2=(nu_-2.*eta_*nu_)/DelS_;
	c3=a05_/120.;
	e3=(DelA_+nu_/2.-2.5*eta_*nu_)/DelS_;
	c4=a06_/720.;
	e4=(1.5*nu_-3.*eta_*nu_)/DelS_;
	c5=a14_/24.;
	e5=(1.5*nu_-1.-2.*eta_*nu_)/DelS_;
	c6=a22_/4.;
	e6=(3.5*nu_-2.-eta_*nu_)/DelS_;
	c7=nu_/alpha_/ubar_/GAMMA_/2.;
	e7=-alpha_/DelS_;

	ddDel_Ax__dtdM=c1*M*(2.*pow(Y,e1)+M*e1*pow(Y,e1-1.)*dY__dM);
	ddDel_Ax__dtdM+=c1*t*M*(2.*e1*pow(Y,e1-1.)*dY__dt+M*e1*(e1-1.)*pow(Y,e1-2.)*dY__dt*dY__dM+M*e1*pow(Y,e1-1.)*ddY__dtdM);
	ddDel_Ax__dtdM+=c2*M*M*M*(4.*e2*pow(Y,e2-1.)*dY__dt+M*e2*(e2-1.)*pow(Y,e2-2.)*dY__dt*dY__dM+M*e2*pow(Y,e2-1.)*ddY__dtdM);
	ddDel_Ax__dtdM+=c3*M*M*M*M*(5.*e3*pow(Y,e3-1.)*dY__dt+M*e3*(e3-1.)*pow(Y,e3-2.)*dY__dt*dY__dM+M*e3*pow(Y,e3-1.)*ddY__dtdM);
	ddDel_Ax__dtdM+=c4*M*M*M*M*M*(6.*e4*pow(Y,e4-1.)*dY__dt+M*e4*(e4-1.)*pow(Y,e4-2.)*dY__dt*dY__dM+M*e4*pow(Y,e4-1.)*ddY__dtdM);
	ddDel_Ax__dtdM+=c5*M*M*M*(4.*pow(Y,e5)+M*e5*pow(Y,e5-1.)*dY__dM);
	ddDel_Ax__dtdM+=c5*t*M*M*M*(4.*e5*pow(Y,e5-1.)*dY__dt+M*e5*(e5-1.)*pow(Y,e5-2.)*dY__dt*dY__dM+M*e5*pow(Y,e5-1.)*ddY__dtdM);
	ddDel_Ax__dtdM+=c6*2.*t*M*(2.*pow(Y,e6)+M*e6*pow(Y,e6-1.)*dY__dM);
	ddDel_Ax__dtdM+=c6*t*t*M*(2.*e6*pow(Y,e6-1.)*dY__dt+M*e6*(e6-1.)*pow(Y,e6-2.)*dY__dt*dY__dM+M*e6*pow(Y,e6-1.)*ddY__dtdM);
	ddDel_Ax__dtdM-=c7*2.*t*e7*pow(Y,e7-1.)*dY__dM;
	ddDel_Ax__dtdM-=c7*t*t*(e7*(e7-1.)*pow(Y,e7-2.)*dY__dt*dY__dM+e7*pow(Y,e7-1.)*ddY__dtdM);
	return(ddDel_Ax__dtdM);
}

double getkapp_2(double t, double M, double Y) {
	double kapp_2;

	kapp_2=t*pow(Y,(2.*nu_-1.)/DelS_);
	kapp_2+=uaster_*ubar_*GAMMA_/2.*M*M*pow(Y,(-eta_*nu_+nu_)/DelS_);
	kapp_2+=a05_/6.*M*M*M*pow(Y,(-1.5*eta_*nu_+DelA_+nu_/2.)/DelS_);
	kapp_2+=a06_/24.*M*M*M*M*pow(Y,(-2.*eta_*nu_+1.5*nu_)/DelS_);
	kapp_2+=a14_/2.*t*M*M*pow(Y,(2.5*nu_-1.-eta_*nu_)/DelS_);
	kapp_2+=a22_/2.*t*t*pow(Y,(3.5*nu_-2.)/DelS_);
	return(kapp_2);
}

void bisection_Y(double *solution, double lY, double rY, double t, double M) {
	double fY1, fY2, c;

	c=(lY+rY)/2.;
	if (rY-c<1.0e-10) {*solution=c;}
	else {
		fY1=1.-(1.-ubar_)*c-ubar_*pow(c,nu_/DelS_)*sqrt(1.+GAMMA_*GAMMA_/getkapp_2(t,M,c));
		fY2=1.-(1.-ubar_)*rY-ubar_*pow(rY,nu_/DelS_)*sqrt(1.+GAMMA_*GAMMA_/getkapp_2(t,M,rY));
		if (fY1*fY2<=0.) {lY=c;}
		else {rY=c;}
		bisection_Y(solution,lY,rY,t,M);
	}
}
	
double getY(double t, double M) {
	int found;
	double Y, Y0, fY0, fY, lY, rY, stepY;
	double kapp_2;

	/* determine the Y region for the bisection rootfinding */
	found=0;
	Y0=0.00000001;
    for (;;) {
	kapp_2=getkapp_2(t,M,Y0);
//printf("t=%lf,M=%lf,Y=%lf,k_2=%lf; ",t,M,Y,kapp_2);
	if (Y0>100.) {
		printf("No solution! Adjust the Y0 limit (100, currently) in getY(double t, double M).\n");
		exit(-1);
	}
//	else if (1.+GAMMA_*GAMMA_/kapp_2<=0.) {
	else if (kapp_2<=0.) {
		Y0+=0.001;
	}
	else {
		break;
	}
    }
//printf("Y0=%lf\n",Y0);
	fY0=1.-(1.-ubar_)*Y0-ubar_*pow(Y0,nu_/DelS_)*sqrt(1.+GAMMA_*GAMMA_/getkapp_2(t,M,Y0));
	stepY=1.;
	for (Y=Y0+stepY; Y<=1000.; Y+=stepY) {
		kapp_2=getkapp_2(t,M,Y);
if (1.+GAMMA_*GAMMA_/kapp_2<0.) {
	printf("Check program!\n");
	exit(-1);
//printf("k2=%lf\n",kapp_2);
}
		fY=1.-(1.-ubar_)*Y-ubar_*pow(Y,nu_/DelS_)*sqrt(1.+GAMMA_*GAMMA_/kapp_2);
//printf("fY: %lf %lf %lf\n",fY0,fY,1.+GAMMA_*GAMMA_/kapp_2);
		if (fY0*fY<=0.) {
			lY=Y-stepY;
			rY=Y;
			found=1;
			break;
		}
	}
	if (!found) {
		printf("Error in searching Y in the given Y range (Y<=1000)!\n");
		exit(-1);
	}
	bisection_Y(&Y,lY,rY,t,M);

if (Y<=0.) {printf("Y=%lf, kapp_2=%lf\n",Y,kapp_2);
	    printf("Negative Y?\n");
	    exit(-1);}

	return(Y);
}

double getdDel_A__dDel_rho(double dDel_Ax__dM) {
	double dDel_A__dDel_rho;

	dDel_A__dDel_rho=crho_*dDel_Ax__dM;
	return(dDel_A__dDel_rho);
}

double getdDel_A__dDel_T(double dDel_Ax__dt, double dDel_Ax__dM) {
	double dDel_A__dDel_T;

	dDel_A__dDel_T=ct_*dDel_Ax__dt-crho_*d1_*dDel_Ax__dM;
	return(dDel_A__dDel_T);
}

double getG(double ddDel_Ax__ddt, double ddDel_Ax__ddM, double ddDel_Ax__dtdM) {
	double G;

	G=(1.-c_*ddDel_Ax__dtdM)*(1.-c_*ddDel_Ax__dtdM);
	G-=c_*c_*ddDel_Ax__ddt*ddDel_Ax__ddM;
	return(G);
}

double getddDel_A__ddDel_rho(double ddDel_Ax__ddt, double ddDel_Ax__ddM, double ddDel_Ax__dtdM) {
	double ddDel_A__ddDel_rho;
	double G;

	G=getG(ddDel_Ax__ddt, ddDel_Ax__ddM, ddDel_Ax__dtdM);
	ddDel_A__ddDel_rho=crho_*crho_*ddDel_Ax__ddM/G;
	return(ddDel_A__ddDel_rho);
}

double getddDel_A__dDel_rhodDel_T(double ddDel_Ax__ddt, double ddDel_Ax__ddM, double ddDel_Ax__dtdM) {
	double ddDel_A__dDel_rhodDel_T;
	double G;

	G=getG(ddDel_Ax__ddt, ddDel_Ax__ddM, ddDel_Ax__dtdM);
	ddDel_A__dDel_rhodDel_T=crho_*ct_*(ddDel_Ax__dtdM-c_*(ddDel_Ax__dtdM*ddDel_Ax__dtdM-ddDel_Ax__ddt*ddDel_Ax__ddM))/G;
	ddDel_A__dDel_rhodDel_T-=crho_*crho_*d1_*ddDel_Ax__ddM/G;
	return(ddDel_A__dDel_rhodDel_T);
}

double getddDel_A__ddDel_T(double ddDel_Ax__ddt, double ddDel_Ax__ddM, double ddDel_Ax__dtdM) {
	double ddDel_A__ddDel_T;
	double G;

	G=getG(ddDel_Ax__ddt, ddDel_Ax__ddM, ddDel_Ax__dtdM);
	ddDel_A__ddDel_T=(ct_*ct_*ddDel_Ax__ddt+crho_*crho_*d1_*d1_*ddDel_Ax__ddM)/G;
	ddDel_A__ddDel_T-=2.*crho_*ct_*d1_*(ddDel_Ax__dtdM-c_*(ddDel_Ax__dtdM*ddDel_Ax__dtdM-ddDel_Ax__ddt*ddDel_Ax__ddM))/G;
	return(ddDel_A__ddDel_T);
}

double getf1(double Y, double kapp_2) {
	double f1;

	f1=1.-ubar_+ubar_*nu_/DelS_*pow(Y,nu_/DelS_-1.)*pow(1.+GAMMA_*GAMMA_/kapp_2,0.5);
	return(f1);
}

double getf3(double Y, double kapp_2) {
	double f3;

	f3=ubar_/2.*pow(Y,nu_/DelS_)*GAMMA_*GAMMA_/kapp_2/kapp_2*pow(1.+GAMMA_*GAMMA_/kapp_2,-0.5);
	return(f3);
}

double getf5(double t, double M, double Y) {
	double f5;

	f5=pow(Y,(2.*nu_-1.)/DelS_)+a14_/2.*M*M*pow(Y,(2.5*nu_-1.-eta_*nu_)/DelS_)+a22_*t*pow(Y,(3.5*nu_-2.)/DelS_);
	return(f5);
}

double getf6(double t, double M, double Y) {
	double f6;

	f6=(2.*nu_-1.)/DelS_*t*pow(Y,(2.*nu_-1.)/DelS_-1.);
	f6+=uaster_*ubar_*GAMMA_/2.*(1.-eta_)*nu_/DelS_*M*M*pow(Y,(-eta_*nu_+nu_)/DelS_-1.);
	f6+=a05_/6.*(DelA_+nu_/2.-3.*eta_*nu_/2.)/DelS_*M*M*M*pow(Y,(-1.5*eta_*nu_+DelA_+nu_/2.)/DelS_-1.);
	f6+=a06_/24.*(1.5*nu_-2.*eta_*nu_)/DelS_*M*M*M*M*pow(Y,(-2.*eta_*nu_+1.5*nu_)/DelS_-1.);
	f6+=a14_/2.*(2.5*nu_-1.-eta_*nu_)/DelS_*t*M*M*pow(Y,(2.5*nu_-1.-eta_*nu_)/DelS_-1.);
	f6+=a22_/2.*(3.5*nu_-2.)/DelS_*t*t*pow(Y,(3.5*nu_-2.)/DelS_-1.);
	return(f6);
}

double getg5(double t, double M, double Y) {
	double g5;

	g5=uaster_*ubar_*GAMMA_*M*pow(Y,(-eta_*nu_+nu_)/DelS_);
	g5+=a05_/2.*M*M*pow(Y,(-1.5*eta_*nu_+DelA_+nu_/2.)/DelS_);
	g5+=a06_/6.*M*M*M*pow(Y,(-2.*eta_*nu_+1.5*nu_)/DelS_);
	g5+=a14_*t*M*pow(Y,(2.5*nu_-1.-eta_*nu_)/DelS_);
	return(g5);
}

double getdY__dt(double f1, double f3, double f5, double f6) {
	return(f3/(f1-f3*f6)*f5);
}

double getdY__dM(double f1, double f3, double g5, double f6) {
	return(f3/(f1-f3*f6)*g5);
}

double getdkapp_2__dt(double f1, double f3, double f5, double f6) {
	return(f1/(f1-f3*f6)*f5);
}

double getdkapp_2__dM(double f1, double f3, double g5, double f6) {
	return(f1/(f1-f3*f6)*g5);
}

double getdf1__dt(double Y, double kapp_2, double dY__dt, double dkapp_2__dt) {
	double df1__dt;

	df1__dt=ubar_*nu_/DelS_*( (nu_/DelS_-1.)*pow(Y,nu_/DelS_-2.)*dY__dt*sqrt(1.+GAMMA_*GAMMA_/kapp_2)-pow(Y,nu_/DelS_-1.)/2./sqrt(1.+GAMMA_*GAMMA_/kapp_2)*GAMMA_*GAMMA_/kapp_2/kapp_2*dkapp_2__dt );
	return(df1__dt);
}

double getdf3__dt(double Y, double kapp_2, double dY__dt, double dkapp_2__dt) {
	double df3__dt;

	df3__dt=ubar_/2.*( nu_/DelS_*pow(Y,nu_/DelS_-1.)*dY__dt*GAMMA_*GAMMA_/kapp_2/kapp_2/sqrt(1.+GAMMA_*GAMMA_/kapp_2)+pow(Y,nu_/DelS_)*GAMMA_*GAMMA_*(-2./pow(kapp_2,3.)*dkapp_2__dt/sqrt(1.+GAMMA_*GAMMA_/kapp_2)+0.5/kapp_2/kapp_2*pow(1.+GAMMA_*GAMMA_/kapp_2,-1.5)*GAMMA_*GAMMA_/kapp_2/kapp_2*dkapp_2__dt) );
	return(df3__dt);
}

double getdf5__dt(double t, double M, double Y, double dY__dt) {
	double df5__dt;

	df5__dt=(2.*nu_-1.)/DelS_*pow(Y,(2.*nu_-1.)/DelS_-1.)*dY__dt;
	df5__dt+=a14_/2.*M*M*(2.5*nu_-1.-eta_*nu_)/DelS_*pow(Y,(2.5*nu_-1.-eta_*nu_)/DelS_-1.)*dY__dt;
	df5__dt+=a22_*pow(Y,(3.5*nu_-2.)/DelS_);
	df5__dt+=a22_*t*(3.5*nu_-2.)/DelS_*pow(Y,(3.5*nu_-2.)/DelS_-1.)*dY__dt;
	return(df5__dt);
}

double getdf6__dt(double t, double M, double Y, double dY__dt) {
	double df6__dt;

	df6__dt=(2.*nu_-1.)/DelS_*(pow(Y,(2.*nu_-1.)/DelS_-1.)+t*((2.*nu_-1.)/DelS_-1.)*pow(Y,(2.*nu_-1.)/DelS_-2.)*dY__dt);
	df6__dt+=uaster_*ubar_*GAMMA_/2.*(1.-eta_)*nu_/DelS_*M*M*((nu_-eta_*nu_)/DelS_-1.)*pow(Y,(nu_-eta_*nu_)/DelS_-2.)*dY__dt;
	df6__dt+=a05_/6.*(DelA_+nu_/2.-3.*eta_*nu_/2.)/DelS_*M*M*M*((DelA_+nu_/2.-1.5*eta_*nu_)/DelS_-1.)*pow(Y,(DelA_+nu_/2.-1.5*eta_*nu_)/DelS_-2.)*dY__dt;
	df6__dt+=a06_/24.*(1.5*nu_-2.*eta_*nu_)/DelS_*M*M*M*M*((1.5*nu_-2.*eta_*nu_)/DelS_-1.)*pow(Y,(1.5*nu_-2.*eta_*nu_)/DelS_-2.)*dY__dt;
	df6__dt+=a14_/2.*(2.5*nu_-1.-eta_*nu_)/DelS_*M*M*(pow(Y,(2.5*nu_-1.-eta_*nu_)/DelS_-1.)+t*((2.5*nu_-1.-eta_*nu_)/DelS_-1.)*pow(Y,(2.5*nu_-1.-eta_*nu_)/DelS_-2.)*dY__dt);
	df6__dt+=a22_/2.*(3.5*nu_-2.)/DelS_*(2.*t*pow(Y,(3.5*nu_-2.)/DelS_-1.)+t*t*((3.5*nu_-2.)/DelS_-1.)*pow(Y,(3.5*nu_-2.)/DelS_-2.)*dY__dt);
	return(df6__dt);
}

double getdg5__dt(double t, double M, double Y, double dY__dt) {
	double dg5__dt;

	dg5__dt=uaster_*ubar_*GAMMA_*M*(nu_-eta_*nu_)/DelS_*pow(Y,(nu_-eta_*nu_)/DelS_-1.)*dY__dt;
	dg5__dt+=a05_/2.*M*M*(DelA_+nu_/2.-1.5*eta_*nu_)/DelS_*pow(Y,(DelA_+nu_/2.-1.5*eta_*nu_)/DelS_-1.)*dY__dt;
	dg5__dt+=a06_/6.*M*M*M*(1.5*nu_-2.*eta_*nu_)/DelS_*pow(Y,(1.5*nu_-2.*eta_*nu_)/DelS_-1.)*dY__dt;
	dg5__dt+=a14_*M*(pow(Y,(2.5*nu_-1.-eta_*nu_)/DelS_)+t*(2.5*nu_-1.-eta_*nu_)/DelS_*pow(Y,(2.5*nu_-1.-eta_*nu_)/DelS_-1.)*dY__dt);
	return(dg5__dt);
}

double getdf1__dM(double Y, double kapp_2, double dY__dM, double dkapp_2__dM) {
	double df1__dM;

	df1__dM=ubar_*nu_/DelS_*( (nu_/DelS_-1.)*pow(Y,nu_/DelS_-2.)*dY__dM*sqrt(1.+GAMMA_*GAMMA_/kapp_2)-pow(Y,nu_/DelS_-1.)/2./sqrt(1.+GAMMA_*GAMMA_/kapp_2)*GAMMA_*GAMMA_/kapp_2/kapp_2*dkapp_2__dM );
	return(df1__dM);
}

double getdf3__dM(double Y, double kapp_2, double dY__dM, double dkapp_2__dM) {
	double df3__dM;

	df3__dM=ubar_/2.*( nu_/DelS_*pow(Y,nu_/DelS_-1.)*dY__dM*GAMMA_*GAMMA_/kapp_2/kapp_2/sqrt(1.+GAMMA_*GAMMA_/kapp_2)+pow(Y,nu_/DelS_)*GAMMA_*GAMMA_*(-2./pow(kapp_2,3.)*dkapp_2__dM/sqrt(1.+GAMMA_*GAMMA_/kapp_2)+0.5/kapp_2/kapp_2*pow(1.+GAMMA_*GAMMA_/kapp_2,-1.5)*GAMMA_*GAMMA_/kapp_2/kapp_2*dkapp_2__dM) );
	return(df3__dM);
}

double getdg5__dM(double t, double M, double Y, double dY__dM) {
	double dg5__dM;

	dg5__dM=uaster_*ubar_*GAMMA_*(pow(Y,(nu_-eta_*nu_)/DelS_)+M*(nu_-eta_*nu_)/DelS_*pow(Y,(nu_-eta_*nu_)/DelS_-1.)*dY__dM);
	dg5__dM+=a05_/2.*(2.*M*pow(Y,(DelA_+nu_/2.-1.5*eta_*nu_)/DelS_)+M*M*(DelA_+nu_/2.-1.5*eta_*nu_)/DelS_*pow(Y,(DelA_+nu_/2.-1.5*eta_*nu_)/DelS_-1.)*dY__dM);
	dg5__dM+=a06_/6.*(3.*M*M*pow(Y,(1.5*nu_-2.*eta_*nu_)/DelS_)+M*M*M*(1.5*nu_-2.*eta_*nu_)/DelS_*pow(Y,(1.5*nu_-2.*eta_*nu_)/DelS_-1.)*dY__dM);
	dg5__dM+=a14_*t*(pow(Y,(2.5*nu_-1.-eta_*nu_)/DelS_)+(2.5*nu_-1.-eta_*nu_)/DelS_*pow(Y,(2.5*nu_-1.-eta_*nu_)/DelS_-1.)*dY__dM);
	return(dg5__dM);
}

double getdf6__dM(double t, double M, double Y, double dY__dM) {
	double df6__dM;

	df6__dM=(2.*nu_-1.)/DelS_*t*((2.*nu_-1.)/DelS_-1.)*pow(Y,(2.*nu_-1.)/DelS_-2.)*dY__dM;
	df6__dM+=uaster_*ubar_*GAMMA_/2.*(1.-eta_)*nu_/DelS_*(2.*M*pow(Y,(nu_-eta_*nu_)/DelS_-1.)+M*M*((nu_-eta_*nu_)/DelS_-1.)*pow(Y,(nu_-eta_*nu_)/DelS_-2.)*dY__dM);
	df6__dM+=a05_/6.*(DelA_+nu_/2.-3.*eta_*nu_/2.)/DelS_*(3.*M*M*pow(Y,(DelA_+nu_/2.-1.5*eta_*nu_)/DelS_-1.)+M*M*M*((DelA_+nu_/2.-1.5*eta_*nu_)/DelS_-1.)*pow(Y,(DelA_+nu_/2.-1.5*eta_*nu_)/DelS_-2.)*dY__dM);
	df6__dM+=a06_/24.*(1.5*nu_-2.*eta_*nu_)/DelS_*(4.*M*M*M*pow(Y,(1.5*nu_-2.*eta_*nu_)/DelS_-1.)+M*M*M*M*((1.5*nu_-2.*eta_*nu_)/DelS_-1.)*pow(Y,(1.5*nu_-2.*eta_*nu_)/DelS_-2.)*dY__dM);
	df6__dM+=a14_/2.*(2.5*nu_-1.-eta_*nu_)/DelS_*t*(2.*M*pow(Y,(2.5*nu_-1.-eta_*nu_)/DelS_-1.)+M*M*((2.5*nu_-1.-eta_*nu_)/DelS_-1.)*pow(Y,(2.5*nu_-1.-eta_*nu_)/DelS_-2.)*dY__dM);
	df6__dM+=a22_/2.*(3.5*nu_-2.)/DelS_*t*t*((3.5*nu_-2.)/DelS_-1.)*pow(Y,(3.5*nu_-2.)/DelS_-2.)*dY__dM;
	return(df6__dM);
}

double getddY__ddt(double f1, double f3, double f5, double f6, 
			double df1__dt, double df3__dt, double df5__dt, double df6__dt) {
	double F3;
	double dF3__dt;

	F3=f3/(f1-f3*f6);
	dF3__dt=df3__dt/(f1-f3*f6)-f3/(f1-f3*f6)/(f1-f3*f6)*(df1__dt-df3__dt*f6-f3*df6__dt);
	return(dF3__dt*f5+F3*df5__dt);
}

double getddY__ddM(double f1, double f3, double g5, double f6,
			double df1__dM, double df3__dM, double dg5__dM, double df6__dM) {
	double F3;
	double dF3__dM;

	F3=f3/(f1-f3*f6);
	dF3__dM=df3__dM/(f1-f3*f6)-f3/(f1-f3*f6)/(f1-f3*f6)*(df1__dM-df3__dM*f6-f3*df6__dM);
	return(dF3__dM*g5+F3*dg5__dM);
}

double getddY__dtdM(double f1, double f3, double g5, double f6,
			double df1__dt, double df3__dt, double dg5__dt, double df6__dt) {
	double F3;
	double dF3__dt;

	F3=f3/(f1-f3*f6);
	dF3__dt=df3__dt/(f1-f3*f6)-f3/(f1-f3*f6)/(f1-f3*f6)*(df1__dt-df3__dt*f6-f3*df6__dt);
	return(dF3__dt*g5+F3*dg5__dt);
}

double getcp_isochoric(double T) {

	double _T, Del_T;           	// temperature (in Kelvin), reduced, reduced difference
	double rho, _rho, Del_rho;		// density, reduced, reduced difference
	double _rho1;
	double _P;			// pressure in Pa, reduced
	double t;				// temperaturelike variable
	double t1;
	double M;				// densitylike order parameter
	double M1;
	double Del_A, Del_Ax, _A0;
	double dDel_A__dDel_rho, dDel_A__dDel_T;
	double ddDel_A__ddDel_rho, ddDel_A__ddDel_T;
	double ddDel_A__dDel_rhodDel_T;
	double dDel_Ax__dM, dDel_Ax__dt;
	double ddDel_Ax__ddM, ddDel_Ax__ddt;
	double ddDel_Ax__dtdM;
	double Del_mu;
	double Y;				// crossover function
	double kapp_2;				// crossover function
	double dY__dt, dY__dM;
	double dkapp_2__dt, dkapp_2__dM;
	double ddY__ddt, ddY__ddM, ddY__dtdM;
	double f1, f3, f6, f5, g5;
	double df1__dt, df3__dt, df6__dt, df5__dt, dg5__dt;
	double df1__dM, df3__dM, df6__dM, dg5__dM;
	double cv, cp, _cv, _cp;		// heat capacity and the reduced
	double _xi;				// reduced xi
	double d_P__d_T, d_P__d_rho;
	int count;

	_T=-Tc/T; Del_T=_T+1.;
	_A0=_A0_+_A1_*Del_T+_A2_*Del_T*Del_T+_A3_*pow(Del_T,3.)+_A4_*pow(Del_T,4.);

	rho=RHOc;
	_rho=rho/RHOc;

	Del_rho=_rho-1.;

	t=gett(Del_T,0.);//ct_*Del_T;		//guess t
	M=getM(Del_T,Del_rho,0.);//crho_*(Del_rho-d1_*Del_T);	//guess M
    count=0;
    for (;;) {
	count++;
	Y=getY(t,M);
	kapp_2=getkapp_2(t,M,Y);
	f1=getf1(Y,kapp_2);
	f3=getf3(Y,kapp_2);
	f5=getf5(t,M,Y);
	f6=getf6(t,M,Y);
	g5=getg5(t,M,Y);
	dY__dt=getdY__dt(f1,f3,f5,f6);
	dY__dM=getdY__dM(f1,f3,g5,f6);
	dDel_Ax__dM=getdDel_Ax__dM(t,M,Y,dY__dM);
	t1=gett(Del_T,dDel_Ax__dM);
	dDel_Ax__dt=getdDel_Ax__dt(t,M,Y,dY__dt);
	M1=getM(Del_T,Del_rho,dDel_Ax__dt);
	if (fabs(t1-t)<1.0e-6*fabs(t) && fabs(M1-M)<1.0e-6*fabs(M)) {
		break;
	}
	else if (count>200) {
		printf("count>200\n");
		break;
	}
	else {
		t+=(t1-t)/10.; M+=(M1-M)/10.;
	}
    }

	Y=getY(t,M);
	kapp_2=getkapp_2(t,M,Y);
	f1=getf1(Y,kapp_2);
	f3=getf3(Y,kapp_2);
	f5=getf5(t,M,Y);
	f6=getf6(t,M,Y);
	g5=getg5(t,M,Y);
	dY__dt=getdY__dt(f1,f3,f5,f6);
	dY__dM=getdY__dM(f1,f3,g5,f6);
	dDel_Ax__dM=getdDel_Ax__dM(t,M,Y,dY__dM);
	t1=gett(Del_T,dDel_Ax__dM);
	dDel_Ax__dt=getdDel_Ax__dt(t,M,Y,dY__dt);
	M1=getM(Del_T,Del_rho,dDel_Ax__dt);
	if (fabs(t1-t)>1.0e-5*fabs(t) || fabs(M1-M)>1.0e-5*fabs(M)) {
		printf("More iterations are needed.\n");
		printf("t1/t-1=%lf M1/M-1=%lf\n",fabs(t1/t-1.),fabs(M1/M-1.));
		exit(-1);
	}
	Del_Ax=getDel_Ax(t,M,Y);
	dDel_Ax__dt=getdDel_Ax__dt(t,M,Y,dY__dt);
	dDel_Ax__dM=getdDel_Ax__dM(t,M,Y,dY__dM);
	Del_A=getDel_A(Del_Ax,dDel_Ax__dt,dDel_Ax__dM);
	Del_mu=getdDel_A__dDel_rho(dDel_Ax__dM);
	_P=_rho*Del_mu-Del_A-_A0;
	_rho1=(_P+Del_A+_A0)/Del_mu;

	rho=_rho*RHOc;

	dkapp_2__dt=getdkapp_2__dt(f1,f3,f5,f6);
	dkapp_2__dM=getdkapp_2__dM(f1,f3,g5,f6);
	df1__dt=getdf1__dt(Y,kapp_2,dY__dt,dkapp_2__dt);
	df3__dt=getdf3__dt(Y,kapp_2,dY__dt,dkapp_2__dt);
	df5__dt=getdf5__dt(t,M,Y,dY__dt);
	df6__dt=getdf6__dt(t,M,Y,dY__dt);
	dg5__dt=getdg5__dt(t,M,Y,dY__dt);
	df1__dM=getdf1__dM(Y,kapp_2,dY__dM,dkapp_2__dM);
	df3__dM=getdf3__dM(Y,kapp_2,dY__dM,dkapp_2__dM);
	dg5__dM=getdg5__dM(t,M,Y,dY__dM);
	df6__dM=getdf6__dM(t,M,Y,dY__dM);
	ddY__ddt=getddY__ddt(f1,f3,f5,f6,df1__dt,df3__dt,df5__dt,df6__dt);
	ddY__ddM=getddY__ddM(f1,f3,g5,f6,df1__dM,df3__dM,dg5__dM,df6__dM);
	ddY__dtdM=getddY__dtdM(f1,f3,g5,f6,df1__dt,df3__dt,dg5__dt,df6__dt);
	Del_Ax=getDel_Ax(t,M,Y);
	dDel_Ax__dt=getdDel_Ax__dt(t,M,Y,dY__dt);
	ddDel_Ax__ddt=getddDel_Ax__ddt(t,M,Y,dY__dt,ddY__ddt);
	dDel_Ax__dM=getdDel_Ax__dM(t,M,Y,dY__dM);
	ddDel_Ax__ddM=getddDel_Ax__ddM(t,M,Y,dY__dM,ddY__ddM);
	ddDel_Ax__dtdM=getddDel_Ax__dtdM(t,M,Y,dY__dt,dY__dM,ddY__dtdM);
	Del_A=getDel_A(Del_Ax,dDel_Ax__dt,dDel_Ax__dM);
	Del_mu=getdDel_A__dDel_rho(dDel_Ax__dM);

	dDel_A__dDel_T=getdDel_A__dDel_T(dDel_Ax__dt,dDel_Ax__dM);
	ddDel_A__ddDel_T=getddDel_A__ddDel_T(ddDel_Ax__ddt,ddDel_Ax__ddM,ddDel_Ax__dtdM);
	ddDel_A__dDel_rhodDel_T=getddDel_A__dDel_rhodDel_T(ddDel_Ax__ddt,ddDel_Ax__ddM,ddDel_Ax__dtdM);
	d_P__d_T=getd_P__d_T(_T,_rho,dDel_A__dDel_T,ddDel_A__dDel_rhodDel_T);
	dDel_A__dDel_rho=getdDel_A__dDel_rho(dDel_Ax__dM);
	ddDel_A__ddDel_rho=getddDel_A__ddDel_rho(ddDel_Ax__ddt,ddDel_Ax__ddM,ddDel_Ax__dtdM);
	d_P__d_rho=getd_P__d_rho(_rho,Del_mu,dDel_A__dDel_rho,ddDel_A__ddDel_rho);
	_xi=1./ddDel_A__ddDel_rho;
	_cv=get_cv(_T,_rho,ddDel_A__ddDel_T);
	_cp=get_cp(_T,_P,_rho,_xi,_cv,d_P__d_T);

	cv=_cv*Pc/Tc/rho;
	cp=_cp*Pc/Tc/rho;

	return(cp);
}

/* ******************

Note.1. Viscosity Formula
	(The critical part of viscosity has been neglected in this approach since it has very weak 
	scailing near critical point. The critical contribution rises to 22% of the background value 
	at T-Tc=0.001 K as mentioned in D. S. Cannell, Phys. Rev. A, Vol.12, p.225 (1975).)

	Ref.[182]: J. H. B. Hoogland, H. R. Van den Berg, and N. J. Trappeniers, Physica A, V.134, 169 (1985)
	Ref.[185]: T. Strehlow and E. Vogel, Physica A, V.161 (1989)

** Fitting initial density viscosities (=m0) at various temperatures in [182] and [185] to the formula 
** From Table III in Ref. [185], but from Eq.(32) in Ref. [182] for T=333.18K 
** 1 \mPa sec l/mol = 1/MW \mPa sec m^3/kg = 6.6159/MW \mPa sec /amagat
** ( 1 amagat = 6.6159 kg/m^3 )

	T (K)		mu1(T) (/mPa m^3/kg)		mu1{T) (/mPa l/mol)	mu1(T) (/mPa /amagat)
	298.15		-0.0001095			-0.016
	323.15		0.0018690			0.273
	373.15		0.0039297			0.574
	423.15		0.0050251			0.734
	473.15		0.0058741			0.858
	523.15		0.0071885			1.050
	573.15		0.0089891			1.313
	623.15		0.0115290			1.684
	673.15		0.0146920			2.146
	333.18		0.0020957						0.013865

 ===> mu1 (\mPa sec m^3/kg) = exp(30.9709*ln(T)+2.794e4/T-3.318e6/T/T-240.06)  !!!!!
 ===> mu1 (\mPa sec /amagat) = 6.6159*exp(30.9709*ln(T)+2.794e4/T-3.318e6/T/T-240.06) !!!!!

** Summary for viscosity calculation

	Ref. [185]	=> mu0=15.262*exp(A*ln(T_R)+B/T_R+C/T_R/T_R+D) \mPa sec
				A=0.548599;
				B=-0.651683;
				C=0.158578;
				D=0.491296;
				T_R=T/298.15;
	Fitting data in Ref. [185] & [182]
			=> mu1=6.6159*exp(30.9709*ln(T)+2.794e4/T-3.318e6/T/T-240.06) \mPa sec/amagat
	Ref. [182]	=> mu2=(2.85099e-3)	// \mPa sec/amagat^2
			   mu3=(-4.16848e-5)	// \mPa sec/amagat^3
			   mu4=(6.446345e-7)	// \mPa sec/amagat^4
			   mu5=(-4.997041e-9)	// \mPa sec/amagat^4
			   mu6=(1.970048e-11)	// \mPa sec/amagat^6
			   mu7=(-2.961807e-14)	// \mPa sec/amagat^7
	==> mu_tp =
mu0+mu1*rho_amagat+mu2*pow(rho_amagat,2.)+mu3*pow(rho_amagat,3.)+mu4*pow(rho_amagat,4.)+mu5*pow(rho_amagat,5.)+mu6*pow(rho_amagat,6.)+mu7*pow(rho_amagat,7.);

	where 1 amagat = 6.6159 kg/m^3


Note.2. Thermal Conductivity

The regular part of thermal conductivity is, see eq.(6)  in H. L. Swinney and D. L. Henry, Phys. Rev. A, Vol.8, p.2586 (1973),
	lambda_regular=lambda_rho+lambda_t

	where
	    lambda_rho = lambda(rho,T) - lambda(rho=0,T) :
		eq.(5) in H. L. Swinney and D. L. Henry, Phys. Rev. A, Vol.8, p.2586 (1973)

	    lambda_t=lambda(rho=0,T) = 0.01303*(1.0 + 0.00549*(tt - 27.5)) :
		thermal conductivity at low density
		from J. Kestin and N. Imaishi, Int. J. Thermophys. V.6, 107 (1985)
		valid 0<tt<100 C

	    and the following data are used, see getlambda_rho(double T, double rho) in this program,

		(rho in kg/m^3, lambda_rho in mW/m/K)
		rhos[0]=1259.558914; lambda_rhos[0]=37.505630; // at T=335.1 K
		rhos[1]=1103.507220; lambda_rhos[1]=32.926088; // at T=324.05 K
		rhos[2]=1030.553735; lambda_rhos[2]=32.276088; // at T=324.05 K
		rhos[3]=923.565772; lambda_rhos[3]=29.029726; // at T=330.15 K
		rhos[4]=837.711879; lambda_rhos[4]=27.929726; // at T=330.15 K
		rhos[5]=800.294483; lambda_rhos[5]=27.309726; // at T=330.15 K
		rhos[6]=770.910691; lambda_rhos[6]=26.759726; // at T=330.15 K
		rhos[7]=714.968208; lambda_rhos[7]=26.649726; // at T=330.15 K
		rhos[8]=685.788605; lambda_rhos[8]=26.069726; // at T=330.15 K
		rhos[9]=647.035399; lambda_rhos[9]=26.719726; // at T=330.15 K
		rhos[10]=615.768265; lambda_rhos[10]=26.119726; // at T=330.15 K
		rhos[11]=611.500579; lambda_rhos[11]=26.209726; // at T=330.15 K
		rhos[12]=590.959511; lambda_rhos[12]=22.049726; // at T=330.15 K
		rhos[13]=568.071442; lambda_rhos[13]=20.809726; // at T=330.15 K
		rhos[14]=482.397293; lambda_rhos[14]=14.969726; // at T=330.15 K
		rhos[15]=465.323008; lambda_rhos[15]=13.409726; // at T=330.15 K
		rhos[16]=393.532482; lambda_rhos[16]=10.009726; // at T=330.15 K
		rhos[17]=344.985050; lambda_rhos[17]=8.039726; // at T=330.15 K
	    from
		J. Lis and P. O. Kellard, Brit. J. Appl. Phys., Vol.16, p.1099 (1965)
		V. V. Burinskii and E. E. Totskii, High Temp., Vol.19, p.366 (1981)

Thermal diffusivity at the critical isochore with background correction, kappa_isochoric, is
	kappa_isochoric=2.4*(T/Tc-1.)^0.62*1.0e-4 cm^2/sec :
		from T. K. Lim, PhD thesis, Johns Hopkins University (1973)

Isobaric heat capacity:
	cp_isochoric=getcp_isochoric(T)

Singular part of lambda at the critical isochore is obtained from the thermal diffusivity:
	lambda_isochoric=2.4*pow(T/Tc-1.,0.62)*1.0e-8*cp_isochoric*RHOc;

(Inverted) coexistence curve, E. S. Wu and W. W. Webb, J. Phys. (Paris), Vol.33, C1~149 (1972), :
	coexistdelrho_t=3.62*pow(T/Tc-1.,0.333)*RHOc;

Assume that the critical isochoric (singular) thermal conductivity follow Lorentzian distribution w.r.t. rho at fixed temperature:
	LHW=coexistdelrho_t/2. : set the half width at half maximum of Lorentzian profile to the coexistence curve
	if (fabs(rho-RHOc)<LHW) {
		lambda_singular=lambda_isochoric/(1.+(rho-RHOc)*(rho-RHOc)/LHW/LHW); //-lambda_isochoric/2.;
	}
//else{//lambda_singular=0.;//}

Total thermal conductivity by summing the regular part and the singular part:

	lambda_trho=lambda_regular+lambda_singular;

******************* */
