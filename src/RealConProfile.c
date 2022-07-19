// This program calculates the real conducted heat in RBC, by integrating the heat equation (shooting method) taken into
// consideration different pressures (hydrostatic) and the fluid properties at any vertical position

#include <ctype.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "libheatcond.h"
#include "libFluidPropC.h"

// If FLAG is set to HEOS, fluid properties are directly calculated via CoolProp
#define FLAG "HEOS"
// If FLAG is set to REFPROP, fluid properties are calculated using REFPROP (NIST)  via the CoolProp
// interface
//#define FLAG "REFPROP"

double e,D,r,dT,Tt=999,Tb=999,height=999,P=999,P0;

extern double z[NX];
extern double q;
extern double h;
char fluid[100]="fluidnotyetdefined";

int main(int argc, char **argv){
	double dt,h,A,AOld,d;
	double tstart,tend,t,ratio,TDeltaOld,Tmax,Tc;
	double T0,lambda_0,rho_0,alpha_0,comp,kappa,nu,cp,psi,lewis;
	double z[NX],p[NX],pOld[NX]; // the vertical coordinate
	double rho[NX],alpha[NX],lambda[NX],lambdaOld[NX];
	double T[NX],T_bOld;
	double diffx,k1,k2,k3,k4;
	double q;
	double intensity,intensityOld;
	double Tm,dT;
	int count,i,j,low,high,k,N,entries;
	char opt;
	FILE *fp;
	
	// handingling arguments
	
	while ((opt=getopt(argc,argv,"f:b:t:P:H:h:")) != -1 ){
		switch (opt){
			case 'f':
				snprintf(fluid,20,"%s",optarg);
				printf("fluid: %s\n",fluid);
				break;
			case 'b':
				Tb = atof(optarg);
				printf("Tb: %lf degC\n",Tb);
				break;
			case 't':
				 Tt= atof(optarg);
				printf("Tt: %lf degC\n",Tt);
				break;
			case 'P':
				P = atof(optarg);
				printf("P: %lf bar\n",P);
				break;
			case 'H':
				height = atof(optarg);
				printf("H: %lf m\n",height);
				break;
			case 'h':
				printUsage(argv);
				break;
			default:
				printUsage(argv);
				exit(-1);
		}	
	}
	//printf("Fluid: %s\n",fluid);
	//printf("Compare: %d\n",strncmp(fluid,"fluidnotyetdefines",16));	
	if (strncmp(fluid,"fluidnotyetdefined",16)==0) {
		printf("Possbible fluids:\n");
		printf(\
	"Fluids:\n"
    "	Air     \n"
    "	Hydrogen\n"
    "	Helium  \n"
    "	Nitrogen\n"
    "	CO2     \n"
    "	Xenon   \n"
    "	SF6     \n"
    "	Ethane  \n"
    "	Water   \n"
    "	Acetone \n"
    "	Methanol\n"
    "	Ethanol \n"
		);
	};

	printf("Height [m]: ");
	scanf("%lf",&height);
	printf("\nP [bar]: ");	
	scanf("%lf",&P);
	printf("Tb [degC]: ");	
	scanf("%lf",&Tb);
	printf("\nTt [degC]: ");
	scanf("%lf",&Tt);
	printf("\n");	
	
	h = height/(NX-1); // there are NX points but only NX-1 intervalls	
	printf("#Tt: %.3f\tTb: %.3f\n",Tt,Tb);
	
	// Calculating starting parameters:
	P0 = 1.e5*P;
	Tb += 273.15;		
	Tt += 273.15;		
	T0 = (Tb+Tt)/2.;
	// Pobulating starting values:
	getCoolProp(fluid, T0, P0, &rho_0, &alpha_0, &comp, &lambda_0, &kappa, &nu, &cp, &psi, &lewis,FLAG);
	//printf("rho_0: %lf\t,lambda: %lf\n",rho_0,lambda_0);
	for (i=0;i<NX;i++){
		rho[i]    = rho_0;
		lambda[i] = lambda_0;
		alpha[i]  = alpha_0;
		z[i] = i*h;
		p[i] = P0;
	}	
	q = 0.01;
	
	for (k=0;k<7;k++){
		q= get_temperature(z,T,lambda,q);
		//get_temperature2(z,T,lambda,q);
    	get_pressure(rho,p,z);
		get_prop(T,p,lambda,rho,alpha);
	}	
	printf("q[W/m^2]: %.6g \n",q);	
	
	
	return 0;
}

double get_temperature(double *z,double *T,double *lambda,double q){
	double err,TDelta,T_tOld,Tmax,TDeltaOld,lowerT;
	double k1,h;
	int count,low,high,i;
	//printf("Tt: %.3f\tTb: %.3f\n",Tt,Tb);
	//printf("inside get_temperature \n");
	h = 1.*height/(NX-1);
	T[0] = Tb;
	err  = 1;
	TDelta=0.001;
	count = 0;
	
	// shotting method
	while (err>0.0001){
		count += 1;
		low    = 0;
		high   = 0;
		T_tOld = T[NX-1];
	    T[0] = Tb; 
		
		// Euler integration
		for (i=0;i<NX;i++) { 
			k1 = -q * h/lambda[i];
			T[i+1] = (T[i]+k1);//6+k2/3+k3/3+k4/6);
			// start value
		};

		//printf("Tb: %lf\tTb_fit: %lf\n",Tb,T[NX-1]);
		Tmax = 0;
	
		// now check how close T[NX-1] is to T_t
		err = (T[NX-1]-Tt)*(T[NX-1]-Tt);
		
		TDeltaOld = TDelta;

		if (T[NX-1]==T_tOld) TDelta = -0.1*(T[NX-1]-Tt);
		else TDelta= -0.1*(T[NX-1]-Tt)*(TDeltaOld)/(T[NX-1]-T_tOld);

		q = q + lambda[0]*TDelta;
	
		T_tOld = T[NX-1];

	}; // end while
	return q;
}

int get_temperature2(double *z, double *T, double *lambda, double q){
	double int_inf_lambda[NX];
	double q_array[NX];
	double h;
	int i;

	h = height/(NX-1);
	
	dT = Tb-Tt;
	int_inf_lambda[0] = lambda[0];
    printf("dT: %lf\tTb: %lf\n",dT,Tb);
	for (i=0;i<NX-1;i++){
	    int_inf_lambda[i+1] = int_inf_lambda[i]+h/lambda[i];
	}	
	
	for (i=0;i<NX;i++){
	    T[i] = Tt + dT*int_inf_lambda[i]/int_inf_lambda[NX-1];
	    if (i>0) q_array[i] = lambda[i]*(T[i]-T[i-1])/h;

	}	
	
	return 0;
}

int get_pressure(double *rho,double *p,double *z){
    double h,p0; //is pressure at half height
	int i;
    h = height/(NX-1);

    p[0]=P0;
    for (i=0;i<NX-1;i++) p[i+1] = p[i]-G*h*rho[i]	;
   	p0 = P0-p[NX/2]; 
    for (i=0;i<NX;i++){
    	p[i] = p[i]+ p0;//P0-p[(NX/2)]	;
   	} 
    
    return 0;
}

int get_prop(double *T, double *p, double *lambda,double *rho, double *alpha){
	int i;
	double temp,press,alpha_0,rho_0,comp,lambda_0,kappa,nu,cp,kappa_0,psi,lewis;

	rho_0   =1;
	lambda_0=1;
	
	for (i=0;i<NX;i++){ 
		temp  = T[i];
		press = p[i]; 
        //printf("%d: temp: %.4lf\t press: %.4lf\n",i,temp,press);
		getCoolProp(fluid, temp, press, &rho_0, &alpha_0, &comp, &lambda_0, &kappa, &nu, &cp, &psi, &lewis,FLAG);
		rho[i]    = rho_0;
		if (isnan(lambda_0)){
			lambda[i] = lambda[i-1];
			printf("Lambda isnan!!");
			printf("i: %d\tT: %.5f\tp: %.2f\n",i,T[i],p[i]);
		} else	lambda[i] = lambda_0;
		
		alpha[i]  = alpha_0;
	}
	return 0;
}

int printUsage(char **argv){
		printf("Usage: %s [option] [value] \n",argv[0]);
		printf("Calculates the conductive heat flux through a fluid of varying heat conductivity\n");
		printf(\
	"Options are optional and don't have to be provided.\n"
	"\t-h\t prints this help message\n"
	"\t-f\t <fluids>\t defines the fluid"
	"\t\n"
	"\t\t Possible fluids:\n"
    "\t\t	Air     \n"
    "\t\t	Hydrogen\n"
    "\t\t	Helium  \n"
    "\t\t	Nitrogen\n"
    "\t\t	CO2     \n"
    "\t\t	Xenon   \n"
    "\t\t	SF6     \n"
    "\t\t	Ethane  \n"
    "\t\t	Water   \n"
    "\t\t	Acetone \n"
    "\t\t	Methanol\n"
    "\t\t	Ethanol \n\n"
	"\t-b <value> bottom plate temperature in degree Celsius\n"
	"\t-t <value> top plate temperature in degree Celsius\n"
	"\t-P <value> pressure in bar \n"
	"\t-H <value> height of the cell in meter \n\n"
	);
	return 0;
}

