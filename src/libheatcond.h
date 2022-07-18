#ifndef REALCONPROFILE_H
#define REALCONPROFILE_H

#define NX 1000
#define PI 3.14159
#define G 9.81

//double get_RealCond(double Tt, double Tb, double P,double height,double *z,double *T,double *p,double *rho,double *lambda);
double get_temperature(double *z,double *T,double *lambda,double q);
int get_temperature2(double *,double *,double *,double );
int get_pressure(double *rho,double *p,double *z);
int get_prop(double *T,double *p,double *lambda,double *rho,double *alpha);
void sf6(double, double, double *, double *,double *, double *, double *, double *, double *);

#endif
