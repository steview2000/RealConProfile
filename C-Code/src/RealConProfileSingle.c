// This program calculates the real conducted heat in RBC, by integrating the heat equation (shooting method) taken into
// consideration different pressures (hydrostatic) and the fluid properties at any vertical position

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gnuplot_i.h"
//#include "header.h"

#define NX 10000
#define HEIGHT 1.12 // in m
#define Tb (21.431+11.858/2)	// in degC
#define Tt (21.431-11.858/2)	// in degC
#define P 8.032e5 // in pa (at half-height)
#define PI 3.14159

double g=9.81; // gravity (in m/s^2)
double e,D,r;
double get_temperature(double *,double *,double *,double );
int get_pressure(double *,double *,double *);
int get_prop(double *,double *,double *,double *, double *);

int main(){
	double dt,h,A,AOld,d;
	double tstart,tend,t,ratio,TDeltaOld,Tmax;
	double P0,T0,lambda_0,rho_0,alfa,comp,kappa,nu,cp,psi,lewis;
	double z[NX]; // the vertical coordinate
	double T[NX],T_bOld;
	double p[NX],pOld[NX],rho[NX],alpha[NX], rhoOld[NX],lambda[NX],lambdaOld[NX];
	double diffx,k1,k2,k3,k4;
	double q;
	double intensity,intensityOld;
	int count,i,j,low,high,k;

		

	// Initialising graphical outpu
	gnuplot_ctrl *h1;
	h1 = gnuplot_init();
	
	h = HEIGHT/(NX-1); // there are NX points but only NX-1 intervalls	
	
	// Calculating starting parameters:
	P0 = 14.504e-5*P;
	T0 = (Tb+Tt)/2.;

	// Pobulating starting values:
	sf6(T0,P0, &rho_0, &alfa, &comp, &lambda_0, &kappa, &nu, &cp);
	printf("rho_0: %lf\t,lambda: %lf\n",rho_0,lambda_0);
	
	// Convert into SI units:
	rho_0 = rho_0*1e3;
	lambda_0 = lambda_0*1e-5;

	for(i=0;i<NX;i++){
		p[i] = P0;
		lambda[i] = lambda_0;
		rho[i] = rho_0;
		z[i] = i*h;
		q = 1;
	}

	for(k=0;k<7;k++){
		q = get_temperature(z,T,lambda,q);
    	get_pressure(rho,p,z);
		get_prop(T,p,lambda,rho,alpha);

	
		//	Temperature profile (shooting method):
		gnuplot_resetplot(h1);
	 	//gnuplot_cmd(h1,"load  \"plot.gnu\" ");
		printf("%i\t q: %.6lf [W/m^2]\n",k,q);	
		getchar();
	}

	return 0;

}//end main

double get_temperature(double *z,double *T,double *lambda,double q){
	double err,TDelta,T_bOld,Tmax,TDeltaOld,lowerT;
	double k1,h;
	int count,low,high,i;
	FILE *fp;

	fp = fopen("temp.csv","w");
	fprintf(fp,"# z[m]\t| T [degC]\n");
	h = 1.*HEIGHT/(NX-1);
	T[0] = Tt;
	err = 1;
	TDelta=0.001;
	count=0;
	while(err>0.0001){
		count++;
		//printf("err: %lf\n",err);
		low=0;
		high=0;
		T_bOld=T[NX-1];
	    T[0] = Tt; 
		
		for(i=0;i<NX;i++){
			k1 = q * h/lambda[i];
//			printf("lambda: %.4lf\n",lambda[i]);
			T[i+1] = (T[i]+k1);///6+k2/3+k3/3+k4/6);
			// start value
		};

	
		//if(count%1000==0) 
//		printf("Tb: %lf\tTb_fit: %lf\n",Tb,T[NX-1]);
		Tmax=0;
	
		// now check how close T[NX-1] is to T_t
		err=sqrt((T[NX-1]-Tb)*(T[NX-1]-Tb));
		
		TDeltaOld = TDelta;

		if(T[NX-1]==T_bOld) TDelta = -0.1*(T[NX-1]-Tb);
		else{	
			TDelta= -0.1*(T[NX-1]-Tb)*(TDeltaOld)/(T[NX-1]-T_bOld);
		};

	//	printf("lambda: %lf\t",lambda[0]);	
		q = q + lambda[0]*TDelta;
	
		T_bOld = T[NX-1];
	//	printf("Err: %lf\t q: %lf\n",err,q);	

	}; //end while
	for(i=0;i<NX;i++)	fprintf(fp,"%.4lf\t|%.4lf\n",i*h,T[i]);
	
	fclose(fp);

	return q;
}

int get_pressure(double *rho,double *p,double *z){
	double h,p0;// p0 is pressure at half height
	int i;
	FILE *fp;

	fp = fopen("pressure.csv","w");
	fprintf(fp,"#z [m]\t|p [pa]\n");
	h = HEIGHT/(NX-1);

	p[0]=0;
	for(i=0;i<(NX-1);i++){
		p[i+1] = p[i]+g*h*rho[i]	;
	}
	p0 = P-p[NX/2];
	for(i=0;i<NX;i++){
		p[i] = p[i]+ p0	;
		fprintf(fp,"%.4lf\t|%.4lf\n",(i*h),p[i]);
	}

	fclose(fp);
	return 1;
}

int get_prop(double *T,double *p,double *lambda,double *rho,double *alpha){
	int i;
	double temp,press,alfa,comp,kappa,nu,cp,rho_0,lambda_0,alpha_0;
	double h;
	FILE *fp;
	rho_0=1;
	lambda_0=1;
	h = 1.*HEIGHT/(NX-1);
	fp = fopen("prop.csv","w");

	fprintf(fp,"#z [m]\t|lamba [J/s m K]\t|rho [kg/m^3]\t|alpha [1/K]\n");
	for(i=0;i<NX;i++){
		temp = T[i];
		press = 14.504e-5*p[i]; // pa -> psi
//		printf("temp: %.4lf\t press: %.4lf\n",temp,press);
		sf6(temp,press, &rho_0, &alpha_0, &comp, &lambda_0, &kappa, &nu, &cp);
		rho[i] = rho_0*1e3;
		lambda[i] = lambda_0*1e-5;
		alpha[i] = alpha_0;
		fprintf(fp,"%lf \t| %lf \t|%lf \t|%lf \n",i*h,lambda[i],rho[i],alpha[i]);
	}
	fclose(fp);
	return 1;
}

