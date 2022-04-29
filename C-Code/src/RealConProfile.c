// This program calculates the real conducted heat in RBC, by integrating the heat equation (shooting method) taken into
// consideration different pressures (hydrostatic) and the fluid properties at any vertical position

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "libheatcond.h"
#include "libFluidPropC.h"

//#define FLUID "SF6"
//#define FLAG "HEOS"
#define FLAG "SLOW"

double e,D,r,dT,Tt,Tb,height,P,P0;

extern double z[NX];
extern double q;
extern double h;
char fluid[100];

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
	double Tm,dT,Ra_m,Ra_c,Nu_m,Nu_c;
	int count,i,j,low,high,k,N,entries;
	FILE *fp;
	
	if (argc != 2){
		printf("Usage: %s [fluid]\n",argv[0]);
		printf("\n");
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
		exit(-1);
	}
	snprintf(fluid,100,argv[1]);

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
	printf("q[W/m^2]: %.4lf \n",q);	
	//	Temperature profile (shooting method):
	//getchar();
	
	//fp = fopen("testout-c.csv","w");
	//fprintf(fp,"# z[m], T,p,rho,lambda\n");
	//for (i=0;i<NX;i++){
	//	fprintf(fp,"%.4lf,%.4lf,%.5g,%.5g,%.5g\n",z[i],T[i],p[i],rho[i],lambda[i]);
	//}
	//fclose(fp);
	
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
			//printf("lambda: %.4lf\n",lambda[i]);
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

		//printf("lambda: %lf\t",lambda[0]);	
		q = q + lambda[0]*TDelta;
	
		T_tOld = T[NX-1];
		//printf("Err: %lf\t q: %lf\n",err,q);	

	}; // end while
	//printf("Done get_temp\n");
	return q;
}

int get_temperature2(double *z, double *T, double *lambda, double q){
	double int_inf_lambda[NX];
	double q_array[NX];
	double h;
	int i;

	//qOut = open("q.txt","w");
	h = height/(NX-1);
	
	dT = Tb-Tt;
	int_inf_lambda[0] = lambda[0];
    printf("dT: %lf\tTb: %lf\n",dT,Tb);
	for (i=0;i<NX-1;i++){
	    //print("lambda: %lf\t%lf\n",lambda[i],h);
	    int_inf_lambda[i+1] = int_inf_lambda[i]+h/lambda[i];
        //printf("int: %lf\n",int_inf_lambda[i+1]);
	}	
	
	for (i=0;i<NX;i++){
        //printf("lambda: %lf\t%lf\n",int_inf_lambda[i],int_inf_lambda[NX-1]);
	    T[i] = Tt + dT*int_inf_lambda[i]/int_inf_lambda[NX-1];
	    if (i>0) q_array[i] = lambda[i]*(T[i]-T[i-1])/h;

	   // fprintf(qOut,"%lf\t%lf\n",T[i],q_array[i]);	
       //	printf("T: %lf\tq: %lf\tlambda: %lf\t\n",T[i],q_array[i],lambda[i]);	
	}	
	//fclose(qOut);
	
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



