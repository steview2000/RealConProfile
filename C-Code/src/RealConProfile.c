// This program calculates the real conducted heat in RBC, by integrating the heat equation (shooting method) taken into
// consideration different pressures (hydrostatic) and the fluid properties at any vertical position

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gnuplot_i.h"
//#include "header.h"

#define NX 1000
#define HEIGHT 1.12 // in m
//#define Tb 40.	// in degC
//#define Tt 13.	// in degC
//#define P 19e5 // in pa (at half-height)
#define PI 3.14159
#define COLN 28


double g=9.81; // gravity (in m/s^2)
double e,D,r;
double get_RealCond(double,  double , double ,double, double *,double *,double *,double *,double *);
double get_temperature(double *,double *,double *, double, double, double, double );
double get_temperature2(double *,double *,double *, double, double, double, double );
int get_pressure(double *,double *,double *,double);
int get_prop(double *,double *,double *,double *);
void sf6(double, double, double *, double *,double *, double *, double *, double *, double *);


int main(){
	double Tb,Tt,P;
	double dt,h,A,AOld,d,Tm,dT,Ra_m,Ra_c,Nu_m,Nu_c;
	double tstart,tend,t,ratio,TDeltaOld,Tmax,Tc;
	double P0,lambda_0,rho_0,alfa,comp,kappa,nu,cp,psi,lewis,lambda_c;
	double z[NX]; // the vertical coordinate
	double T[NX],T_bOld;
	double p[NX],pOld[NX],rho[NX],rhoOld[NX],lambda[NX],lambdaOld[NX];
	double matrix[COLN];
	double diffx,k1,k2,k3,k4;
	double q;
	double intensity,intensityOld;
	char DataSet[6],junk[100],line[1000];
	int count,i,j,low,high,k,N,entries;

	FILE *fin,*fout;
	
	fin = fopen("nus_table_Tm_Tc_HPCF4b.csv","r");
	fout = fopen("q.csv","w");
	fprintf(fout,"#DataSet: \t|Ra_m:\t|Ra_c\t|Nu_m:\t|Nu_c:\t| q_av_m:\t|q_real:\t|q_av_c: \n");	
		
	// read the header - check how many columns there are
	fgets(line,1000,fin);
	printf("%s\n",line);

	printf("%i: ",i);
	entries = 1000;
	while(entries>10){
		// Read each line
		for(k=0;k<COLN;k++){
			if (EOF == fscanf(fin,"%s\t",junk)){
				entries=0;
				break;
			}
			else{
			//	printf("%lf\t",P0); 
			 	if(k==0) sprintf(DataSet,"%s",junk);
				if(k==2) P = 1e5*atof(junk);
				if(k==3) Tm=atof(junk);
				if(k==5) dT=atof(junk);
				if(k==6) Ra_m=atof(junk);
				if(k==8) Nu_m=atof(junk);
				if(k==14) Tc=atof(junk);
				if(k==15) Ra_c=atof(junk);
				if(k==17) Nu_c=atof(junk);
			}

		}
		printf("P: %.3lf (bar)",P*1e-5);
		fscanf(fin,"\n");
		printf("\n");
		
		Tb = Tm+dT/2.;
		Tt = Tm-dT/2.;
		q=get_RealCond(Tt,Tb,P,HEIGHT,z,T,p,rho,lambda);

		sf6(Tm,14.504e-5*P, &rho_0, &alfa, &comp, &lambda_0, &kappa, &nu, &cp);
		sf6(Tc,14.504e-5*P, &rho_0, &alfa, &comp, &lambda_c, &kappa, &nu, &cp);
		lambda_0 = lambda_0*1e-5;
		lambda_c = lambda_c*1e-5;
		printf("lambda: %lf (W/m s K)\t rho: %lf (g/cm^3)\n",lambda_0,rho_0);
		fprintf(fout,"%s\t|%lf\t|%lf\t|%lf\t|%lf\t|%lf\t|%lf\t|%lf\n",DataSet,Ra_m,Ra_c,Nu_m,Nu_c,dT*lambda_0/HEIGHT,q,dT*lambda_c/HEIGHT);	

	} // end while loop

	fclose(fout);
	fclose(fin);
	return 0;

}//end main


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
	printf("Height: %.3lf(m)\tTm: %lf(degC)\t P: %lf(bar)\n",height,Tm,1e-5*P)	;
	// Convert into SI units:
	rho_0 = rho_0*1e3;
	lambda_0 = lambda_0*1e-5;
	printf("rho_0: %lf\t,lambda: %lf\n",rho_0,lambda_0);
	
	for(i=0;i<NX;i++){
		p[i] = P0;
		lambda[i] = lambda_0;
		rho[i] = rho_0;
		z[i] = i*h;
		q = 1;
	}
	// Iteration to find correct T(z), p(z), rho(z) and lambda(z)
	for(i=0;i<7;i++){
		q = get_temperature(z,T,lambda,q,Tt,Tb,P);
		get_pressure(rho,p,z,P);
		get_prop(T,p,lambda,rho);
	
	}	
	//	Temperature profile (shooting method):
	printf("Iteration: %i\t q: %.6lf [W/m^2]\n",i,q);	
	return q;
}

double get_temperature(double *z,double *T,double *lambda,double q,double Tt, double Tb, double P){
	double err,TDelta,T_bOld,Tmax,TDeltaOld,lowerT;
	double k1,h;
	int count,low,high,i;
	FILE *fp;

	fp = fopen("temp.csv","w");
	fprintf(fp,"# z[m]\t| T [degC]\n");
	h = 1.*HEIGHT/NX;
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
		
		for(i=0;i<(NX-1);i++){
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

double get_temperature2(double *z,double *T,double *lambda,double q,double Tt, double Tb,double P){
	double int_inf_lambda[NX];
	double q_array[NX];
	double h,dT;
	
	int i;
	FILE *qOut;

	qOut = fopen("q.txt","w");

	h = HEIGHT/(NX-1);
	
	dT = Tb-Tt;
	int_inf_lambda[0] = lambda[0];
//	printf("dT: %lf\tTb: %lf\n",dT,Tb);
	for(i=0;i<(NX-1);i++){
	//	printf("lambda: %lf\t%lf\n",lambda[i],h);
		int_inf_lambda[i+1] = int_inf_lambda[i]+h/lambda[i];
//		printf("int: %lf\n",int_inf_lambda[i+1]);
	};
	
	for(i=0;i<NX;i++){
//		printf("lambda: %lf\t%lf\n",int_inf_lambda[i],int_inf_lambda[NX-1]);
		T[i] = Tt + dT*int_inf_lambda[i]/int_inf_lambda[NX-1];
		if(i>0) q_array[i] = lambda[i]*(T[i]-T[i-1])/h;
		fprintf(qOut,"%lf\t%lf\n",T[i],q_array[i]);	
//		printf("T: %lf\tq: %lf\tlambda: %lf\t\n",T[i],q_array[i],lambda[i]);	
	};
	fclose(qOut);
	
	return q_array[NX-1];
}

int get_pressure(double *rho,double *p,double *z,double P){
	double h,p0;// p0 is pressure at half height
	int i;
	FILE *fp;

	fp = fopen("pressure.csv","w");
	fprintf(fp,"#z [m]\t|p [pa]\n");
	h = HEIGHT/NX;

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

int get_prop(double *T,double *p,double *lambda,double *rho){
	int i;
	double temp,press,alfa,comp,kappa,nu,cp,rho_0,lambda_0;
	double h;
	FILE *fp;
	rho_0=1;
	lambda_0=1;
	h = 1.*HEIGHT/NX;
	fp = fopen("prop.csv","w");

	fprintf(fp,"#z [m]\t|lamba [J/s m K]\t|rho [kg/m^3]\n");
	for(i=0;i<NX;i++){
		temp = T[i];
		press = 14.504e-5*p[i]; // pa -> psi
//		printf("temp: %.4lf\t press: %.4lf\n",temp,press);
		sf6(temp,press, &rho_0, &alfa, &comp, &lambda_0, &kappa, &nu, &cp);
		rho[i] = rho_0*1e3;
		lambda[i] = lambda_0*1e-5;
		fprintf(fp,"%.4lf \t| %.4lf \t|%.4lf \n",i*h,lambda[i],rho[i]);
	}
	fclose(fp);
	return 1;
}

