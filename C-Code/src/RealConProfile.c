// This program calculates the real conducted heat in RBC, by integrating the heat equation (shooting method) taken into
// consideration different pressures (hydrostatic) and the fluid properties at any vertical position

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "libheatcond.h"

//#define HEIGHT 1.12 // in m
#define COLN 28

double e,D,r;

int main(){
	double Tb,Tt,P,height;
	double dt,h,A,AOld,d,Tm,dT,Ra_m,Ra_c,Nu_m,Nu_c;
	double tstart,tend,t,ratio,TDeltaOld,Tmax,Tc;
	double P0,lambda_0,rho_0,alfa,comp,kappa,nu,cp,psi,lewis,lambda_c;
	double matrix[COLN];
	double diffx,k1,k2,k3,k4;
	double q;
	int count,i,j,low,high,k,N,entries;

	printf("Height [m]: ");
	scanf("%lf",&height);
	printf("\nP [bar]: ");	
	scanf("%lf",&P);
	P = 1e5*P;
	printf("\nTb [degC]: ");	
	scanf("%lf",&Tb);
	printf("\nTt [degC]: ");
	scanf("%lf",&Tt);
	printf("\n");	
	q=get_RealCond(Tt,Tb,P,height,z,T,p,rho,lambda);

	return 0;

}

