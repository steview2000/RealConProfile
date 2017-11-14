
/*******************************************************************************\
	RBC.c calls one of
		cdo()
		fn_5cb()
		water()
		acetone()
		nitrogen()
		helium()
		sf6()
                sf6_crit()
		xenon()
		ethanol()
		methanol()
		isopropanol()
		toluene2()
		glycerol()
		h2_xe()
		he_sf6()  (does not work yet !!! )
	to calculate properties of the respective fluids and of a cylindrical
	convection geometry. The syntax is, for instance,

	cdo(temp, press, &rho, &alfa, &comp, &lambda, &kappa, &nu, &cp); 

	The pressure press is in psi, and the temperature temp in deg. C. 
	The returned parameters are
				
 		rho     =	density ( g/cm3 )                   
		alfa 	=	thermal expansion coeff		(1/K )                       
		comp	=	compressibility	 		( cm^2/dyne )
						(this is sometimes set to zero)               
		lambda	= 	thermal conductivity 	( erg/s cm K )                    
		kappa 	=	thermal diffusivity 		( cm2/s )                         
		nu 	=	shear viscosity 			( g /s cm )                      
		cp 	=   	spec. heat at const. press. 	(erg/g K )

	For the mixtures, the syntax is slightly different, e.g.:
 
	h2_xe(temp, press, x, &rho, &alfa, &comp, &lambda, &kappa, &nu, &cp, &psi, &lewis);

	with
		x	=	molar concentration of heavier component 
		psi	=	separation ratio
		lewis	=	Lewis number

\*******************************************************************************/

#define DIAMETER  19  		/* diameter of cylindrical cell in cm */
#define HEIGHT 10.2
#define CHOICE 2			/* 0 = fix dtc (onset), 1 fix height (onset), 2 fix height (turbulence) */
#include "header.h"
#include "header.def"

main()
{
double temp, rho, alfa, comp, lambda, kappa, nu, cp, tscale;
double temp1, rho1, alfa1, comp1, lambda1, kappa1, nu1, cp1;
double temp2, rho2, alfa2, comp2, lambda2, kappa2, nu2, cp2;
double gamma[5], corr[5][5];
double p0, p1, p2, p3, p4;
double q, rh20, rr20, l2;
double epsa, epsb, epsr;
double X, Y, Z, D0, a, b, epsT;
double omeg, taylor, tv;
double a2;
double dtc, qcrit, g, sigma, sigma1, sigma2, dsigma, a_r, n_index, n_1, n_2, beta, shadow_coeff;
double height, diameter, aspect, press;
double g3, g3_tilde, eps_crit, qc, xi0, tau0, F_th;
int fix, gas, i, j;
double x, psi, psi1, psi2, lewis, lewis1, lewis2, rayleigh;
double dtemp, volume, nusselt, lambda_eff, current, t_b, x_WL, alf_dt, drc, rc_gamma, a_g;
double tempTop, tempBot;
	g = G;  /*  cm/s^2   */

/* decide which output option is desired */

	if(CHOICE == 0)
		printf("\nDelta T_c is user specified. \nOutput for pattern formation near onset.\n");
	if(CHOICE == 1)
		printf("\nCell thickness is fixed by #define HEIGHT ... . \nOutput for pattern formation near onset.\n");
	if(CHOICE == 2)
		printf("\nCell thickness is fixed by #define HEIGHT ... . \nFor turbulent flows at large Rayleigh numbers.\n");
	printf("\nIf you would like a different option, change #define CHOICE \nin RBC.c and re-compile.\n\n");
	fix = CHOICE;	/* 0 = fix dtc, all else fix height */

/* input */

	printf("Which fluid ? \n\t 1 = H2 \t 2 = He \t 3 = N2 \t 4 = CO2 \t 5 = Xe  \n\t 6 = SF6  \t 7 = SF6_crit ");
		printf("\n\t10 = 5CB \t11 = H2O \t12 = ACETONE \t13 = METHANOL ");
		printf("\t14 = ETHANOL \n\t15 = 2-PROPANOL\t16 = TOLUENE\t17 = GLYCEROL ");
		printf("\n\t20 = H2-Xe \t21 = He-SF6\n");
	scanf("%d", &gas);
//	if(gas == 21){
//		printf("\nSorry, He-SF6 is under construction.\n\n");
//		exit(0);
//	}

        if((gas < 10) || (gas > 19)){
                printf("mean temperature (deg C) and pressure (psi) ?   ");
                scanf("%lf %lf", &temp, &press);
        }
        else{
                printf("mean temperature (deg C) ?   ");
                scanf("%lf", &temp);
                press = 1.;
        }

	x = 0.; psi = 0.; lewis = 0.; 	/* for mixtures only : */
	if(gas >= 20){
		printf("Molar concentration of the heavier component ?  ");
		scanf("%lf", &x);
	}

	if(fix == 2){
		printf("Temperature difference ?  ");
		scanf("%lf", &dtemp);
	}

/* print out a general header */

	fn_header(gas, temp, press, x);

/* get properties at the mean temperature. NOTE: nu is SHEAR viscosity. */
	get_prop(gas, temp, press, x, &rho, &alfa, &comp, &lambda, &kappa, &nu, &cp, &psi, &lewis);
	sigma = nu/(rho*kappa);				/* Prandtl number */
/* adiabatic coefficient */

        a_g = alfa*(temp + T0)/(rho*cp);
        a_g = a_g*rho*G;

/* determine height and Delta T_c */

	diameter = DIAMETER;
	if(fix == 0){
		printf("Critical temperature difference (deg C) ?   ");
		scanf("%lf", &dtc);
		height = RC*kappa*nu/(rho*alfa*G*dtc);	
		height = pow(height, 1./3.);
	}
	else{
		height = HEIGHT;
		dtc = RC*kappa*nu/(rho*alfa*G*height*height*height);
	}
	aspect = diameter/height;
	tscale = dtc/RC;

/* compute properties of turbulent state */

	if(fix == 2){
		volume = 1.e-3*height*diameter*diameter*PI/4.;   /* cell volume in liter */
		rayleigh = dtemp*rho*alfa*G*height*height*height/(kappa*nu);
/*
		nusselt = 0.3265*pow(rayleigh, 0.25)*pow(sigma, -1./12.);		nusselt += 0.002352*pow(rayleigh, 3./7.)*pow(sigma, -1./7.);
		nusselt = nusselt/(1. + 0.181*log10(aspect));
*/
		nusselt = NofR_Oregon(sigma, rayleigh);
		lambda_eff = lambda*nusselt;
	}

/* parameters needed near onset of RBC */

	qc = 3.117;
	xi0 = sqrt(0.148);
	tau0 = 1./(19.65*sigma/(sigma+0.5117));
		
/* thermal noise parameters, see Hohenberg and Swift, Phys. Rev. A 46, 4773 (1992). */

	F_th = (KB*(temp + T0)*rho/(height*nu*nu))*2.*sigma*qc/(xi0*tau0*RC);  /* Eq. 2.17 */
	g3 = 1.4;
	g3_tilde = (2./3.)*g3;
	eps_crit = 3.*g3_tilde*F_th/4.;
	eps_crit = pow(eps_crit, 2./3.);

/* get properties at hot and cold end */

	if((fix == 0) || (fix == 1)){
		temp1 = temp + dtc/2.;
		temp2 = temp - dtc/2.;
	}
	else if(fix == 2){
		temp1 = temp + dtemp/2.;
		temp2 = temp - dtemp/2.;
	}
	else{
		printf("\nno such option.\n");
		exit(0);
	}

	get_prop(gas, temp1, press, x, &rho1, &alfa1, &comp1, &lambda1, &kappa1, &nu1, &cp1, &psi1, &lewis1);
	get_prop(gas, temp2, press, x, &rho2, &alfa2, &comp2, &lambda2, &kappa2, &nu2, &cp2, &psi2, &lewis2);
	sigma1 = nu1*cp1/lambda1;
	sigma2 = nu2*cp2/lambda2;

/* some refractive index stuff. Needed for shadowgraph sensitivity. */

    /* for refractive index calculation. OK only for CO2 and SF6.  */
        if(gas == 4)	
            a_r = 0.152;		/* CO2 */
        else if(gas == 6)
            a_r = 0.0777;		/* SF6 */
        else
            a_r = 0.;			/* all else */
	n_index = sqrt((1.+2.*a_r*rho)/(1.-a_r*rho));    /* refractive index */
	n_1 = sqrt((1.+2.*a_r*rho1)/(1.-a_r*rho1));
	n_2 = sqrt((1.+2.*a_r*rho2)/(1.-a_r*rho2));
	beta = (n_1 - n_2)/dtc;
	shadow_coeff = 2.*PI*PI/height*beta*dtc/1708.*385.28*sqrt(F_th)/2.;


/* non-Boussinesq things */
	
	if(fix == 0 || fix == 1)
		dsigma = 2./(sigma1+sigma2)*(sigma2-sigma1)/dtc;
	else if(fix == 2)
		dsigma = 2./(sigma1+sigma2)*(sigma2-sigma1)/dtemp;
	
	gamma[0] = -(rho1 - rho2)/rho;
	gamma[1] = (rho1*alfa1-rho2*alfa2)/(2.*rho*alfa);
	gamma[2] = (nu1/rho1 - nu2/rho2)/(nu/rho);
	gamma[3] = (lambda1 - lambda2)/lambda;
	gamma[4] = (cp1 - cp2)/cp;

	/* Busse coefficients. These are not used. We use the new Pesch coefficients.*/
/*	 
	printf("NOB coeffs from Busse 1967\n");
	p0 = 2.676 - 0.1258/sigma;
	p1 = -6.603 - 0.5023/sigma;
	p2 = 2.755;
	p3 = 2.917 - 0.5023/sigma;
	p4 = -6.229 + 0.2512/sigma;
*/
	/* Pesch coefficients */
	printf("NOB coeffs from BPA00, Annu. Rev. Fluid Mech. 32, 709 (2000), Sect. 6.5\n");
	p0 = 2.676 - 0.361/sigma;
	p1 = -6.631 - 0.772/sigma;
	p2 = 2.765;
	p3 = 9.540;
	p4 = -6.225 + 0.3857/sigma;

	l2 = 0.29127 + 0.08147/sigma + 0.08933/(sigma*sigma);
	rh20 = 0.89360 + 0.04959/sigma + 0.06787/(sigma*sigma);
	rr20 = 0.69942 - 0.00472/sigma + 0.00832/(sigma*sigma);
	q = gamma[0]*p0 + gamma[1]*p1 + gamma[2]*p2 + gamma[3]*p3 + gamma[4]*p4;

/* these formulas are from Busse 1967*/
	epsa = -q*q/(4.*rh20*RC); 
	epsb = q*q*(9.*rh20 - 3.*l2)/(l2*l2*RC);		
	epsr = q*q*3.*rr20/(l2*l2*RC);

/* the following is from BdAC91 [PRL 67, 3078 (1991)] */
	X = 1.39892 - 0.00944/sigma + 0.01665/(sigma*sigma);
	Y = 1.98148 + 0.15349/sigma + 0.19529/(sigma*sigma);
	a = q*sqrt(6./(RC*X));
	a2 = a*a;
	b = Y/X;
	D0 = 1. + b - 2.*b*b;
	Z = a*b/D0 - (a/fabs(a))*sqrt(a*a*b*b/(D0*D0) + a*a/(2.*D0));
		/* Z is corrected ala footnote on p. 744 of BPA00 */ 
	epsT = a*Z+(2.*b+1)*Z*Z; 	

/* the following formulas are from BdAC91 [PRL 67, 3078 (1991)], and give the same results */
/*  as the Busse formulas used above for epsa, etc. 
	epsa2 = -a2/(4.+8.*b);
	epsr2 = a2/((1.-b)*(1.-b));
	epsb2 = a2*(2.+b)/((1.-b)*(1.-b));
*/

/* second order non-OB shift of R_c */
/* coefficients corr[i][j] from W. Pesch, private commun. */

        bzero(corr,25*sizeof(double));
        corr[0][0] = 10.122;
        corr[1][0] = -12.504;
        corr[1][1] = -12.984;
        corr[2][0] = 0.;
        corr[2][1] = 76.833;
        corr[2][2] = -88.059;
        corr[3][0] = -0.00157;
        corr[3][1] = 153.738;
        corr[3][2] = -76.844;
        corr[3][3] = -76.923;
        corr[4][0] = 6.252;
        corr[4][1] = -51.861;
        corr[4][2] = 38.417;
        corr[4][3] = 32.436;
        corr[4][4] = -3.246;
        
        drc = 0.;
        for(i = 0; i <= 4; i++){
            for(j = 0; j <= i; j++){
                drc += gamma[i]*corr[i][j]*gamma[j];
            }
        }
        rc_gamma = 1707.8 + drc;

/* correct dT_c for NOB effects */
/* Note that the gammas's are based on R_c = 1707.8 when dT_c is calculated. */

	dtc = dtc*rc_gamma/1707.8;

/* Boundary-layer calculations, Eq. 8 of Wu and L., Phys. Rev. A 43, 2833 (1991) */

	if(fix == 2){
	/* I will use an estimate of the average temperature of each BL to get the properties : */
		temp1 = temp + dtemp/4.;
		temp2 = temp - dtemp/4.;
		get_prop(gas, temp1, press, x, &rho1, &alfa1, &comp1, &lambda1, &kappa1, &nu1, &cp1, &psi1, &lewis1);
		get_prop(gas, temp2, press, x, &rho2, &alfa2, &comp2, &lambda2, &kappa2, &nu2, &cp2, &psi2, &lewis2);
		x_WL = alfa1*(nu2/rho2)*kappa2/(alfa2*(nu1/rho1)*kappa1);
		x_WL = pow(x_WL, 1./3.)*lambda1/lambda2;
		alf_dt = alfa*dtemp;		/* NOB criterion used by the Oregon group */
	}

/* End of non-Boussinesq */

// Addition by S. Weiss (January 2007)
/* Lets calculate the real temperature difference between the top of top plate and bottom of bottom plate */

	double term1,term2,nenner,realdtc;
	double d,d1,d2;
	double k,k1,k2,faktor;
	k = lambda*1e-5;// thermal conductivity of the gas in J/(s m K)
        k1= 145.7;  	// thermal conductivity of SI in J/(s m K) (from: www.rfcafe.com/references/general/thermal_conductivity.htm )
	k2= 35	;	// thermal conductivity of Sapphire in J/(s m K) (from: www.rfcafe.com/references/general/thermal_conductivity.htm )
	d = height*1e-2;	// gap between both plates;
	d1= 9.525e-3;	// thickness of bottom plate in m
	d2= 9.525e-3;	// thickness of top plate in m
	
	// top of bottom plate and bottom of top plate:
	faktor=1+k*d1/(d*k1)+k*d2/(d*k2);
	// calculate the real temperature (Regarding S. Weiss Notebook #1,pg.:101)
	realdtc=dtc*faktor;
		
/* Output */

/* common to all options: cell geometry and fluid props at mean temp. */

	printf("\nHeight = %.4f cm\t\tAspect Ratio = %.3f\n", height, aspect);
	printf("rho = %.6f g/cm3\n",rho);	
	printf("alpha = %.6f 1/K		comp = %.4e cm^2/dyne\n",alfa, comp);
	printf("lambda = %.1f erg/(s cm K)	kappa = %.4e cm^2/s\n",lambda, kappa);
	printf("eta = %.7f g/s cm		nu = %.4e g/(s cm)\n",nu, nu/rho);
	printf("cp = %.3e erg/g K		sigma = %.3f 	(1/sigma)(dsigma/dT) = %.5f\n", cp, sigma, dsigma);
	if(gas >= 20)
		printf("x = %.3f\t\tPsi = %.3f\t\t\tLewis = %.3f\n\n", x, psi, lewis);
	printf("\n");

/* Output relevant to RBC onset */

	if(fix == 0 || fix == 1){
		qcrit = dtc*lambda*diameter*diameter/4.*PI/height*1.e-7; 
		tv = height*height*rho/nu;
		omeg = 2.*3.14159*tv*0.25;
		taylor = 4.*omeg*omeg;
                printf("gamma_i :	");
                for(i = 0; i < 5; i++)
                    printf("%.4f   ", gamma[i]);
                printf("\n");
		printf("F_th = %.3e\t\teps_crit = %.3e\n", F_th, eps_crit);
		printf("Q_Busse = %.3f\t\tEpsa  = %.3e\t\tEpsr  = %.3e\n", q, epsa, epsr);
		printf("EpsT' = %.3e\t\tEpsb  = %.3e\n", epsT, epsb);
                printf("R_c^gamma = %.1f	dR_c = %.3f	a_g = %.3e K/cm\n", rc_gamma, drc, a_g);
		printf("delta tc = %.4e C		Qcrit = %.4e W\n",dtc, qcrit);

		printf("visc. tv = %.3f s		therm. tv = %.3f s\n\n",tv, tv*sigma);	  
		printf("1/4 Hz yields Omega = %.2f and\tTaylor = %.1f\n", omeg, taylor);
		if((gas == 4) || (gas == 6)){
//			sens = -9.44e4*beta*height*dtc;
			printf("Refractive index n = %.5f   dn/dT = %.3e\n\n", n_index, beta);
		}
		printf("\nAssuming the bottom plate is %.3lf mm thick Si and the top plate is %.3lf mm thick Sapphire.\n",d1*1e3,d2*1e3);
		printf("The temperature difference you will measure is: %lf K\n\n",realdtc);
		printf("The conversion factor: %lf K\n\n",faktor);

	}

/* Output relevant to turbulent RBC */

	else if(fix == 2){
		printf("Q_Busse = %.3f\t\tx_WL = %.4f\t\talpha*DeltaT = %.4f\n", q, x_WL, alf_dt);
		qcrit = dtc*lambda*diameter*diameter/4.*PI/height*1.e-7; 
		printf("delta tc = %.4e C		Qcrit = %.4e W\n",dtc, qcrit);
		current = dtemp*lambda*nusselt*height*aspect*height*aspect/4.*PI/height*1.e-7; 
		printf("Rayleigh = %.2e\t\tNusselt = %.2f\t\tPrandtl = %.3f\n", rayleigh, nusselt, sigma);
		printf("Q = %.3f W\t\t\tVolume = %.3f liter\n", current, volume);

		tv = height*height*rho/nu;
		t_b = height*height*rho*cp/lambda_eff;

		printf("visc. tv = %.3f s\ttherm. tv = %.3f s\tlambda_eff = %.3e ergs/s cm K\n\n",tv, t_b, lambda_eff);	  }
	else{
		printf("\nno such option.\n");
		exit(0);
	}



}



void fn_header(int gas, double temp, double press, double x)
{

printf("gas = %d\n", gas);

	if(gas == 1)
		printf("\n\nH2\t\tT_bar = %6.3f deg C\t\tP = %6.2f psi\n\n",temp, press);
	else if(gas == 2)
		printf("\n\nHe\t\tT_bar = %6.3f deg C\t\tP = %6.2f psi\n\n",temp, press);
	else if(gas == 3)
		printf("\n\nN2\t\tT_bar = %6.3f deg C\t\tP = %6.2f psi\n\n",temp, press);
	else if(gas == 4)
		printf("\n\nCO2	\t\tT_bar = %6.3f deg C\t\tP = %6.2f psi\n\n",temp, press);
	else if(gas == 5)
		printf("\n\nXe\t\tT_bar = %6.3f deg C\t\tP = %6.2f psi\n\n",temp, press);
	else if(gas == 6)
		printf("\n\nSF6\t\tT_bar = %6.3f deg C\t\tP = %6.2f psi\n\n",temp, press);
	else if(gas == 7)
		printf("\n\nSF6 near the critical point\t\tT_bar = %6.3f deg C\t\tP = %6.2f psi\n\n",temp, press);
	else if(gas == 10)
		printf("\n\n5CB\t\tT_bar = %6.3f deg C\n\n",temp);
	else if(gas == 11)
		printf("\n\nH2O	\t\tT_bar = %6.3f deg C\n\n",temp);
	else if(gas == 12)
		printf("\n\nACETONE\t\tT_bar = %6.3f deg C\n\n",temp);
	else if(gas == 13)
		printf("\n\nMETHANOL\t\tT_bar = %6.3f deg C\n\n",temp);
	else if(gas == 14)
		printf("\n\nETHANOL\t\tT_bar = %6.3f deg C\n\n",temp);
	else if(gas == 15)
		printf("\n\n2-PROPANOL\t\tT_bar = %6.3f deg C\n\n",temp);
	else if(gas == 16)
		printf("\n\nTOLUENE\t\tT_bar = %6.3f deg C\n\n",temp);
	else if(gas == 17)
		printf("\n\nGLYCEROL\t\tT_bar = %6.3f deg C\n\n",temp);
	else if(gas == 20)
		printf("\n\nH2-Xe\t\tT_bar = %6.3f deg C\t\tP = %6.2f psi\t\tX_Xe = %6.3f\n\n",temp, press, x);
	else if(gas == 21)
		printf("\n\nHe-SF6\t\tT_bar = %6.3f deg C\t\tP = %6.2f psi\t\tX_SF6 = %6.3f\n\n",temp, press, x);
	else{
		printf("no such fluid\n");
		exit(0);
	}
}



void get_prop(gas, temp, press, x, prho, palfa, pcomp, plambda, pkappa, pnu, pcp, ppsi, plewis)
int gas;
double temp, press, x, *prho, *palfa, *pcomp, *plambda, *pkappa, *pnu, *pcp, *ppsi, *plewis;
{
	if(gas == 1){
		h2_xe(temp, press, x, prho, palfa, pcomp, plambda, pkappa, pnu, pcp, ppsi, plewis);
	}
	else if(gas == 2){
		helium(temp, press, prho, palfa, pcomp, plambda, pkappa, pnu, pcp);
	}
	else if(gas == 3){
		nitrogen(temp, press, prho, palfa, pcomp, plambda, pkappa, pnu, pcp);
	}
	else if(gas == 4){
		cdo(temp, press, prho, palfa, pcomp, plambda, pkappa, pnu, pcp);
	}
	else if(gas == 5){
		xenon(temp, press, prho, palfa, pcomp, plambda, pkappa, pnu, pcp);
	}
	else if(gas == 6){
		sf6(temp, press, prho, palfa, pcomp, plambda, pkappa, pnu, pcp);
	}
	else if(gas == 7){
		sf6_crit(temp, press, prho, palfa, pcomp, plambda, pkappa, pnu, pcp);
	}
	else if(gas == 10){
		fn_5cb(temp, press, prho, palfa, pcomp, plambda, pkappa, pnu, pcp);
	}
	else if(gas == 11){
		water2(temp, press, prho, palfa, pcomp, plambda, pkappa, pnu, pcp);
	}
	else if(gas == 12){
		acetone(temp, press, prho, palfa, pcomp, plambda, pkappa, pnu, pcp);
	}
	else if(gas == 13){
		methanol(temp, press, prho, palfa, pcomp, plambda, pkappa, pnu, pcp);
	}
	else if(gas == 14){
		ethanol(temp, press, prho, palfa, pcomp, plambda, pkappa, pnu, pcp);
	}
	else if(gas == 15){
		isopropanol(temp, press, prho, palfa, pcomp, plambda, pkappa, pnu, pcp);
	}
	else if(gas == 16){
		toluene2(temp, press, prho, palfa, pcomp, plambda, pkappa, pnu, pcp);
	}
	else if(gas == 17){
		glycerol(temp, press, prho, palfa, pcomp, plambda, pkappa, pnu, pcp);
	}
	else if(gas == 20){
		h2_xe(temp, press, x, prho, palfa, pcomp, plambda, pkappa, pnu, pcp, ppsi, plewis);
	}
	else if(gas == 21){
		he_sf6(temp, press, x, prho, palfa, pcomp, plambda, pkappa, pnu, pcp, ppsi, plewis);
	}
	else{
		printf("no such fluid available\n");
		exit(0);
	}
}



double NofR_Oregon(double sig, double r)
{

/* A fit to the Oregon data, after correction for wall conductance using Model 2. */

double n, ln, lr, a, b;
	a = -1.140182;
	b = 0.324721;
	lr = log10(r);
	ln = a + b*lr;
	n = pow(10., ln);
	return(n);
}
