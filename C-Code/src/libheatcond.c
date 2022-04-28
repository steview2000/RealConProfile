// This program calculates the real conducted heat in RBC, by integrating the heat equation (shooting method) taken into
// consideration different pressures (hydrostatic) and the fluid properties at any vertical position

// 12.2.2016: I will make the following changes

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "libheatcond.h"

#define NX 10000

// Function that actually calculates the real conduction profile
double get_RealCond(double Tt, double Tb, double P,double height,double *z,double *T,double *p,double *rho,double *lambda){
	double P0,Tm,lambda_0,rho_0,alfa,comp,kappa,nu,cp,psi,lewis;
	double q,h;
	int i;
	
	h = height/(NX-1);
	
	// Calculating starting parameters:
	P0 = 14.504e-5*P; // pa -> psi
	Tm = (Tb+Tt)/2.;
	// Pobulating starting values:
	sf6(Tm,P0, &rho_0, &alfa, &comp, &lambda_0, &kappa, &nu, &cp);
	printf("Height: %.3lf(m)\nTm: %.3lf (degC)\n P: %.3lf (bar)\n",height,Tm,1e-5*P)	;
	
	// Convert into SI units:
	rho_0    = rho_0*1e3;
	lambda_0 = lambda_0*1e-5;
	
	for(i=0;i<NX;i++){
		p[i]      = P0;
		lambda[i] = lambda_0;
		rho[i]    = rho_0;
		z[i]      = i*h;
		q         = 1;
	}
	// Iteration to find correct T(z), p(z), rho(z) and lambda(z)
	for(i=0;i<7;i++){
		q = get_temperature(z,T,lambda,q,Tt,Tb,P);
		get_pressure(rho,p,z,P);
		get_prop(T,p,lambda,rho);
	}	
	//	Temperature profile (shooting method):
	printf("cond heat transp: %.6lf [W/m^2]\n",i,q);	
	return q;
}

double get_temperature(double *z,double *T,double *lambda,double q,double Tt, double Tb, double P){
	double err,TDelta,T_bOld,Tmax,TDeltaOld,lowerT;
	double k1,h;
	int count,low,high,i;
	
	h      = z[1]-z[0];
	T[0]   = Tt;
	err    = 1;
	TDelta = 0.001;
	count  = 0;
	while(err>0.0001){
		count++;
		low    = 0;
		high   = 0;
		T_bOld = T[NX-1];
	    T[0]   = Tt;
		
		for(i=0;i<(NX-1);i++){
			k1 = q * h/lambda[i];
			T[i+1] = (T[i]+k1);
			// start value
		};
	
		Tmax=0;
	
		// now check how close T[NX-1] is to T_t
		err=sqrt((T[NX-1]-Tb)*(T[NX-1]-Tb));
		
		TDeltaOld = TDelta;

		if(T[NX-1]==T_bOld) TDelta = -0.1*(T[NX-1]-Tb);
		else{	
			TDelta= -0.1*(T[NX-1]-Tb)*(TDeltaOld)/(T[NX-1]-T_bOld);
		};

		q = q + lambda[0]*TDelta;
	
		T_bOld = T[NX-1];

	}; //end while

	return q;
}

//double get_temperature2(double *z,double *T,double *lambda,double q){
//	double int_inf_lambda[NX];
//	double q_array[NX];
//	double h,dT;
//	int i;
//	FILE *qOut;
//
//	qOut = fopen("q.txt","w");
//
//	h = HEIGHT/(NX-1);
//	
//	dT = Tb-Tt;
//	int_inf_lambda[0] = lambda[0];
////	printf("dT: %lf\tTb: %lf\n",dT,Tb);
//	for(i=0;i<(NX-1);i++){
//	//	printf("lambda: %lf\t%lf\n",lambda[i],h);
//		int_inf_lambda[i+1] = int_inf_lambda[i]+h/lambda[i];
////		printf("int: %lf\n",int_inf_lambda[i+1]);
//	};
//	
//	for(i=0;i<NX;i++){
////		printf("lambda: %lf\t%lf\n",int_inf_lambda[i],int_inf_lambda[NX-1]);
//		T[i] = Tt + dT*int_inf_lambda[i]/int_inf_lambda[NX-1];
//		if(i>0) q_array[i] = lambda[i]*(T[i]-T[i-1])/h;
//		fprintf(qOut,"%lf\t%lf\n",T[i],q_array[i]);	
////		printf("T: %lf\tq: %lf\tlambda: %lf\t\n",T[i],q_array[i],lambda[i]);	
//	};
//	fclose(qOut);
//	
//	return q_array[NX-1];
//}

int get_pressure(double *rho,double *p,double *z,double P){
	double h,p0;// p0 is pressure at half height
	int i;
	
	h = z[1]-z[0];

	p[0]=0;
	for(i=0;i<(NX-1);i++){
		p[i+1] = p[i]+G*h*rho[i]	;
	}
	p0 = P-p[NX/2];
	for(i=0;i<NX;i++){
		p[i] = p[i]+ p0	;
	}
	return 1;
}

int get_prop(double *T,double *p,double *lambda,double *rho){
	int i;
	double temp,press,alfa,comp,kappa,nu,cp,rho_0,lambda_0,alpha_0;
	double h;
	FILE *fp;
	rho_0=1;
	lambda_0=1;

	for(i=0;i<NX;i++){
		temp = T[i];
		press = 14.504e-5*p[i]; // pa -> psi
		sf6(temp,press, &rho_0, &alpha_0, &comp, &lambda_0, &kappa, &nu, &cp);
		rho[i] = rho_0*1e3;
		lambda[i] = lambda_0*1e-5;
	}
	return 1;
}

