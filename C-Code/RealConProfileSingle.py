# This program calculates the real conducted heat in RBC, by integrating the heat equation (shooting method) taken into
# consideration different pressures (hydrostatic) and the fluid properties at any vertical position
#
# 12.2.2016: I will make the following changes

import sys
import numpy as np
import matplotlib.pyplot as py
#sys.path.append('/home/sweiss/Convection/FluidsProp') # for the FluidProp module
import FluidProp as fluid

NX= 1000
HEIGHT= 1.12 # in m
Tb = 25.0 # in degC
Tt = 15.0 # in degC
P = 8.0e5 # in pa (at half-height)
PI= 3.14159

g=9.81; # gravity (in m/s^2)

#double get_temperature(double *,double *,double *,double );
#double get_temperature2(double *,double *,double *,double );
#int get_pressure(double *,double *,double *);
#int get_prop(double *,double *,double *,double *, double *);

def get_temperature(z,T,lambda,q):
#	double err,TDelta,T_bOld,Tmax,TDeltaOld,lowerT;
#	double k1,h;
#	int count,low,high,i;
#	FILE *fp;
	print("inside get_tem\n");
	fp = open("temp.csv","w");
	fp.write("# z[m]\t| T [degC]\n");
	h = 1.*HEIGHT/(NX-1);
	T[0] = Tt;
	err = 1;
	TDelta=0.001;
	count=0;
	while(err>0.0001):
		count+=1;
		low=0;
		high=0;
		T_bOld=T[NX-1];
	        T[0] = Tt; 
		
		for i in range(NX):
			k1 = q * h/lambda[i];
			#printf("lambda: %.4lf\n",lambda[i]);
			T[i+1] = (T[i]+k1);#6+k2/3+k3/3+k4/6);
			# start value
		};

	
		#printf("Tb: %lf\tTb_fit: %lf\n",Tb,T[NX-1]);
		Tmax=0;
	
		# now check how close T[NX-1] is to T_t
		err=sqrt((T[NX-1]-Tb)*(T[NX-1]-Tb));
		
		TDeltaOld = TDelta;

		if(T[NX-1]==T_bOld):
                    TDelta = -0.1*(T[NX-1]-Tb);
		else:
                    TDelta= -0.1*(T[NX-1]-Tb)*(TDeltaOld)/(T[NX-1]-T_bOld);
		

		#printf("lambda: %lf\t",lambda[0]);	
		q = q + lambda[0]*TDelta;
	
		T_bOld = T[NX-1];
		print("Err: %lf\t q: %lf\n",err,q);	

	}; # end while
	for i in range(NX):
            fp.write("%.4lf\t|%.4lf\n"%(i*h,T[i]));
	
	fp.close();

	return z,T,lambda,q;


def get_temperature2(z,T,lambda, q):
	int_inf_lambda = np.zeros(NX)
        q_array = np.zeros(NX);

	qOut = open("q.txt","w");
	h = HEIGHT/(NX-1);
	
	dT = Tb-Tt;
	int_inf_lambda[0] = lambda[0];
        #printf("dT: %lf\tTb: %lf\n",dT,Tb);
	for i in range(NX-1):
	    #print("lambda: %lf\t%lf\n",lambda[i],h);
	    int_inf_lambda[i+1] = int_inf_lambda[i]+h/lambda[i];
            #printf("int: %lf\n",int_inf_lambda[i+1]);
	
	
	for i in range(NX):
            #printf("lambda: %lf\t%lf\n",int_inf_lambda[i],int_inf_lambda[NX-1]);
	    T[i] = Tt + dT*int_inf_lambda[i]/int_inf_lambda[NX-1];
	    if (i>0):
                q_array[i] = lambda[i]*(T[i]-T[i-1])/h;
	    qOut.writw("%lf\t%lf\n"%(T[i],q_array[i]));	
            #	printf("T: %lf\tq: %lf\tlambda: %lf\t\n",T[i],q_array[i],lambda[i]);	
	
	qOut.close();
	
	return z,T,lambda,q_array[NX-1];


def get_pressure(rho,p,z):
    #	double h,p0;// p0 is pressure at half height
    fp = open("pressure.csv","w");
    fp.write(fp,"#z [m]\t|p [pa]\n");
    h = HEIGHT/(NX-1);

    p[0]=0;
    for i in range(NX-1):
        p[i+1] = p[i]+g*h*rho[i]	;
	p0 = P-p[NX/2];
    
    for i in range(NX):
    	p[i] = p[i]+ p0	;
    	fp.write("%.4lf\t|%.4lf\n"%((i*h),p[i]);
    
    
    fp.close();
    return rho,p,z;


def get_prop(T,p,lambda,rho,alpha):
	#int i;
	#double temp,press,alfa,comp,kappa,nu,cp,rho_0,lambda_0,alpha_0;
	#double h;
	rho_0=1;
	lambda_0=1;
	h = 1.*HEIGHT/(NX-1);
	fp = open("prop.csv","w");

	fp.write("#z [m]\t|lamba [J/s m K]\t|rho [kg/m^3]\t|alpha [1/K]\n");
	for i in range(NX):
		temp = T[i];
		press = 14.504e-5*p[i]; // pa -> psi
                #printf("temp: %.4lf\t press: %.4lf\n",temp,press);
		rho,lambda,nu,kappa,alpha = fluid.SF6_C(temp,press);
		rho[i] = rho_0*1e3;
		lambda[i] = lambda_0*1e-5;
		alpha[i] = alpha_0;
		fp.write(fp,"%lf \t| %lf \t|%lf \t|%lf \n"%(i*h,lambda[i],rho[i],alpha[i]));
	}
	fclose(fp);
	return T,p,lambda,rho,alpha;
}

if __name__== "__main__":
    z = np.zeros(NX) 
	T = np.zeros(NX) 
    p = np.zeros(NX)
    pOld=np.zeros(NX)
    rho = np.zeros(NX)
    alpha=np.zeros(NX)
    lambda=np.zeros(NX)
    lambdaOld=np.zeros(N)

	h = HEIGHT/(NX-1); # there are NX points but only NX-1 intervalls	
	
	# Calculating starting parameters:
	P0 = 14.504e-5*P;
	T0 = (Tb+Tt)/2.;

	# Pobulating starting values:
	rho,lambda,nu,kappa,alpha =fluid.SF6_C(T0,P0);
	print("rho_0: %lf\t,lambda: %lf\n"%(rho_0,lambda_0));
	
	# Convert into SI units:
	rho_0 = rho_0*1e3;
	lambda_0 = lambda_0*1e-5;

	for i in range(N):
		p[i] = P0;
		lambda[i] = lambda_0;
		rho[i] = rho_0;
		z[i] = i*h;

	q = 1;
	for k in range(7):
		q,T,lambda,q = get_temperature2(z,T,lambda,q);

    	rho,p,z = get_pressure(rho,p,z);
	T,p,lambda,rho,alpha=get_prop(T,p,lambda,rho,alpha);
	
	#	Temperature profile (shooting method):
	printf("%i\t q: %.6lf [W/m^2]\n",k,q);	
	getchar();

	return 0;


