/**********************************************************************
sf6.c

SF_6 gas thermal properties (SI unit.)  The program is good (<2%) under the following conditions (Cp may be within 4%):

p=1 bar t=0 -- 100 C
p=5 bar t=0 -- 100 C
p=10 bar t=0 -- 100 C
p=15 bar t=10 -- 100 C
p=20 bar t=20 (25 ?) -- 100 C

When P > 20 bar, the results have larger errors.

The experimental data and formula are compiled from over 20 papers.  No single
paper gives a complete documentation.  These papers are documented by Jun Liu.

11/21/94 by Jun Liu.

Use a new method to calculate density and thermal expansion coefficient.  The
method is based on a full state equation and Newtonian method.  The results 
are more precise than those obtained by using the second virial coefficient.

4/6/95, Jun Liu.
10/10/95, Last modified.

Modified by G.A. Feb. 4, 97. Input should now be T in deg C and P in psi.
All returned values are in cgs.

***********************************************************************/

#include "../include/header.h"
#define NP 		100
#define JMAX 	26
#define DDT 	0.001


/* ******************************************************************** */

sf6(tt, pp, prho, palfa, pcomp, plambda, pkappa, pnu, pcp, ppsi, plewis)
double tt, pp, *prho, *palfa, *pcomp, *plambda, *pkappa, *pnu, *pcp, *ppsi, *plewis;
{
	double lambda_t, lambda_tp;		/*  thermal conductivity  */
	double lomu_t, mu_t, mu_tp;		/*  viscosity  */
	double rr, rrd, B_t, C_t, D_t;		/*  virial coeff. */
	double tk;             				/*  temperature  */
	double cp_t, cp_tp;				/*  heat capacity */
	double ppa, alp;

	pp = pp/14.504;		/* convert pressure from psi to bar */ 

	tk = tt + ZERO_C;		/*  convert C into Kelvin */
	ppa = pp*1.0e+5;		/*  convert bar into Pa=N/M^2  */

/*  properties at low density */

/* From J. Kestin and N. Imaishi, Int. J. Thermophys. V.6, 107 (1985).  */

	lambda_t = 0.01303*(1.0 + 0.00549*(tt - 27.5));

/* From J. H. B. Hoogland, H. R. Van den Berg, and N. J. Trappeniers, Physica A, V.134, 169 (1985).  */

	lomu_t = 0.51460*log(tk) - 190.1/tk + 10500/(tk*tk) + 0.3114;
	mu_t = exp(lomu_t);

/* Fit data cp in F. J. Uribe, E. A. Mason, and J. Kestin, {\it J. Phys. Chem. Ref. Data}, {\bf 19}, 1123 (1990).  */

	cp_t = 91.0719 + 0.2581*tt - 4.504e-4*tt*tt;  

/*  Virial coefficients derived from Bulletin JSME, Vol. 26, 1590 (1983).  */

	B_t = MW_SF6*MW_SF6*(0.054820 - 35.874/tk - 1.1443e+6/pow(tk,3))/R;
	C_t = MW_SF6*MW_SF6*MW_SF6*(-1.64315e-4 + 0.100787/tk - (3.03945e+3)/pow(tk,3))/R;
	D_t = MW_SF6*MW_SF6*MW_SF6*MW_SF6*(6.674565e-7 - 3.714493e-4/tk +15.73237/pow(tk,3))/R;
		           
        expansion(tk, ppa, &rr, &alp); 
	
/* ---  properties as a function of both temperature and pressure.  ----  */

/* From J. Kestin and N. Imaishi, Int. J. Thermophys. V.6, 107 (1985).  */

	lambda_tp = lambda_t + 1.1e-5*rr - 7.615e-9*rr*rr + 2.887e-10*rr*rr*rr;

/* Fit to data in J. H. B. Hoogland, H. R. Van den Berg, and N. J. Trappeniers, Physica A, V.134, 169 (1985).  */

	rrd = rr/6.6159;
	mu_tp = mu_t + 1.3865e-2*rrd + 2.85099e-3*rrd*rrd - 4.16848e-5*pow(rrd, 3) + 6.446345e-7*pow(rrd, 4);
	mu_tp = mu_tp - 4.99704e-9*pow(rrd, 5) + 1.970048e-11*pow(rrd, 6) - 2.961807e-14*pow(rrd, 7);
	mu_tp = mu_tp*1.0e-6;	    /*  times 1.0e-6 to get the right unit */

/* pressure dependence of cp is derived from W. Braker and A. L. Mossman, {\sl The Matheson Unabridged Gas Data Book : a Compilation of Physical and Thermodynamic Properties of Gases},  (Matheson Gas Products, East Rutherford, 1974).  However, we do not have enough data to show reliable pressure dependence.  */

	cp_tp = cp_t*(1.0 + (0.0023 + 1.346e-5*tt)*(pp-1.0) + (63.56 - 0.7144*tt)*1.0e-5*(pp-1.0)*(pp-1.0));
			
/* convert to cgs and assign values */

	*plambda = 1.e5*lambda_tp;
	*pnu = 10.*mu_tp;	
	*palfa = alp;
	*prho = 0.001*rr;
	*pcp = 1.e7*cp_tp/MW_SF6;
	*pcomp = 0.;
	*pkappa = *plambda/((*prho)*(*pcp));

	return;
}

/* ************************************************************************ */

/*  The following subroutine uses the equation of state in Bulletin of the
    JSME, V.26, 1590 (1983) by Oda et al. to calculate rho and expansion
    coefficient.  The Newton-Raphson method is used to solve the equation.  */

expansion(tk, ppa, rho, expan)
double tk, ppa, *rho, *expan;
  {
double dtk;             		           /*  temperature  */
double rrho, drrho;
double all;
double sf6_super_vir();

	rrho = sf6_super_vir(tk, ppa);
	dtk = tk + DDT;
        drrho = sf6_super_vir(dtk, ppa);
	all = (rrho - drrho)/(DDT*rrho);

	*rho = rrho;
	*expan = all;
	return;
  }

double sf6_super_vir(double temp, double press)
 {
double rho;
double rho0;
double xacc,x1,x2;	
double sf6_rtnewt();

	xacc = 1.e-5;
	rho0 = press*MW_SF6/(R*temp);
	x1 = 0.1*rho0;
	x2 = 1.9*rho0;
	rho = sf6_rtnewt(x1,x2,xacc,temp,press);
	return rho;
  }

double sf6_rtnewt(x1,x2,xacc,temp,press)
double x1,x2,xacc,temp,press;
{
int j;
double df,dx,f,rtn;

	rtn=0.5*(x1+x2);
	for (j=1;j<=JMAX;j++) {
		sf6_funcd(rtn,temp,press,&f,&df);
		dx=f/df;
		rtn -= dx;
		if ((x1-rtn)*(rtn-x2) < 0.0){
			printf("Jumped out of brackets in RTNEWT\n");
			exit(0);
		}
		if (fabs(dx) < xacc) return rtn;
	}
	printf("Maximum number of iterations exceeded in RTNEWT\n");
}

sf6_funcd(rho,temp,press,func,dfunc)
double rho, temp, press, *func,*dfunc;
  {
double a2, a3, a4, a5, a6, a7, a8, a60, a01, a02, b, tt2;
double fun, dfun;
	
	tt2 = 1.0/(temp*temp);

	a2 = 5.482e-2*temp - 35.874 - (1.1443e+6)*tt2;
	a3 = -1.64315e-4*temp + 0.100787 - (3.03945e+3)*tt2;
	a4 = 6.674565e-7*temp - 3.714493e-4 +15.73237*tt2;
	a5 = -9.112486e-10*temp + 4.8984006e-7 - 1.8778235e-2*tt2;
	a6 = 7.3006586e-13*temp - 3.6544374e-10 + 1.0840107e-5*tt2;
	a7 = -2.707784e-16*temp + 1.3897509e-13 - 3.388251e-9*tt2;
	a8 = 3.884286e-20*temp - 2.045054e-17 + 4.666427e-13*tt2;
	a60 = -3.4941e-17;
	a01 = -1.126e+2;
	a02 = -5.176e-4;
	b = 7.0e-6;

	fun = R*temp*rho/MW_SF6 + a2*rho*rho + a3*rho*rho*rho - press;
	fun = fun + a4*pow(rho,4.) + a5*pow(rho,5.) + a6*pow(rho,6.);
	fun = fun + a7*pow(rho,7.) + a8*pow(rho,8.) + a60*temp*temp*pow(rho,6.);
	fun = fun + tt2*(a01*rho*rho*rho + a02*pow(rho,5.))*(1.0 + b*rho*rho)*exp(-b*rho*rho);

	dfun = R*temp/MW_SF6 + 2.*a2*rho + 3.*a3*rho*rho;
	dfun = dfun + 4.*a4*pow(rho,3.) + 5.*a5*pow(rho,4.) + 6.*a6*pow(rho,5.);
	dfun = dfun + 7.*a7*pow(rho,6.) + 8.*a8*pow(rho,7.) + 6.*a60*temp*temp*pow(rho,5.);
	dfun = dfun + tt2*(3.*a01*rho*rho + 5.*a02*pow(rho,4.))*(1.0 + b*rho*rho)*exp(-b*rho*rho);
	dfun = dfun - 2.0*tt2*(a01*rho*rho*rho + a02*pow(rho,5.))*b*b*rho*rho*rho*exp(-b*rho*rho);

	*func = fun;
	*dfunc = dfun;

	return;
}
