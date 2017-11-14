
# This program calculates the real conducted heat in RBC, by integrating the heat equation (shooting method) taken into
# consideration different pressures (hydrostatic) and the fluid properties at any vertical position
#
# 12.2.2016: I will make the following changes
#
# see our paper: Shishkina, Weiss, Bodenschatz, PRF, 1, 062301 (R)

import sys
import numpy as np
sys.path.append('/home/sweiss/Convection/FluidsProp') # for the FluidProp module
import FluidProp as fluid
reload(fluid)

#NX= 1000
#HEIGHT= 1.12 # in m for the G=1 cell HPCF
#HEIGHT=0.3 # in m for the He cell of Urban et al.
#Tb= 30.#(21.431+11.858/2) # in degC
#Tt= 18.#(21.431-11.858/2)	# in degC
#P= 8.032e5 # in pa (at half-height)
g=9.81; # gravity (in m/s^2)

def get_temperature(z,T,Tt,Tb,lambd,q):
#	print("inside get_tem\n");
        NX = np.size(z)
        HEIGHT = z[-1]-z[0]
	#fp = open("temp.csv","w");
	#fp.write("# z[m]\t| T [degC]\n");
	h = 1.*HEIGHT/(NX-1);
	T[0] = Tb;
	err = 1;
	TDelta=0.001;
	count=0;
	while(err>0.0001):
		count+=1;
		low=0;
		high=0;
		T_tOld=T[NX-1];
	        T[0] = Tb; 
		
		for i in range(NX-1):
			k1 = -q * h/lambd[i];
			#printf("lambd: %.4lf\n",lambd[i]);
			T[i+1] = (T[i]+k1);#6+k2/3+k3/3+k4/6);
			# start value
	
		#printf("Tb: %lf\tTb_fit: %lf\n",Tb,T[NX-1]);
		Tmax=0;
	
		# now check how close T[NX-1] is to T_t
		err=np.sqrt((T[NX-1]-Tt)*(T[NX-1]-Tt));
		
		TDeltaOld = TDelta;

		if(T[NX-1]==T_tOld):
                    TDelta = -0.1*(T[NX-1]-Tt);
		else:
                    TDelta= -0.1*(T[NX-1]-Tt)*(TDeltaOld)/(T[NX-1]-T_tOld);
		

		#printf("lambd: %lf\t",lambd[0]);	
		q = q + lambd[0]*TDelta;
	
		T_tOld = T[NX-1];
		#print("Err: %lf\t q: %lf\n",err,q);	

	for i in range(NX):
            fp.write("%lf\t|%lf\n"%(i*h,T[i]));
	
	fp.close();

	return T,q;


def get_temperature2(z,T,Tt,Tb,lambd, q):
        NX = np.size(z)
        HEIGHT = z[-1]-z[0]
	int_inf_lambd = np.zeros(NX)
        q_array = np.zeros(NX);

	qOut = open("q.txt","w");
	h = HEIGHT/(NX-1);
	
	dT = Tb-Tt;
	int_inf_lambd[0] = lambd[0];
        #printf("dT: %lf\tTb: %lf\n",dT,Tb);
	for i in range(NX-1):
	    #print("lambd: %lf\t%lf\n",lambd[i],h);
	    int_inf_lambd[i+1] = int_inf_lambd[i]+h/lambd[i];
            #printf("int: %lf\n",int_inf_lambd[i+1]);
	
	for i in range(NX):
            #printf("lambd: %lf\t%lf\n",int_inf_lambd[i],int_inf_lambd[NX-1]);
	    T[i] = Tb - dT*int_inf_lambd[i]/int_inf_lambd[NX-1];
	    if (i>0):
                q_array[i] = -lambd[i]*(T[i]-T[i-1])/h;
            
	    qOut.write("%lf\t%lf\t%lf\n"%(T[i],q_array[i],(Tb-Tt)/int_inf_lambd[NX-1]));	
            #	printf("T: %lf\tq: %lf\tlambd: %lf\t\n",T[i],q_array[i],lambd[i]);	
	
	qOut.close();
	return T,(Tb-Tt)/int_inf_lambd[NX-1];


def get_pressure(rho,p,Pin,z):
    #	double h,p0;// p0 is pressure at half height
    fp = open("pressure.csv","w");
    fp.write("#z [m]\t|p [pa]\n");
    NX = np.size(z)
    h = (z[-1]-z[0])/(NX-1);

    p[0]=0;
    for i in range(NX-1):
        p[i+1] = p[i]-g*h*rho[i]	;
   
    p0 = Pin-p[NX/2];
    
    for i in range(NX):
    	p[i] = p[i]+ p0	
        fp.write("%lf\t|%lf\n"%(i*h,p[i]));
    fp.close();
    return p;


def get_prop_C(T,p):
        NX = np.size(T)
        rho = np.zeros(NX)
        lambd=np.zeros(NX)
        alpha=np.zeros(NX)
	rho_0=1;
	lambd_0=1;
	fp = open("prop.csv","w");
        #print("size alpha: %i"%np.size(alpha)) 

	#fp.write("#z [m]\t|lamba [J/s m K]\t|rho [kg/m^3]\t|alpha [1/K]\n");
	for i in range(NX):
		temp = T[i];
		press = p[i]; # pa -> bar
		rho_0,lambd_0,nu,kappa,alpha_0 = fluid.SF6_C(temp,press*1e-5);
		rho[i] = rho_0;
		lambd[i] = lambd_0;
		alpha[i] = alpha_0;
	#	fp.write("%lf \t| %lf \t|%lf \t|%lf \n"%(i*h,lambd[i],rho[i],alpha[i]));
#                print "Got something"	
	#fp.close();
	return lambd,rho,alpha;

def get_prop_NIST(fluidType,T,p):
        NX = np.size(T)
        rho = np.zeros(NX)
        lambd=np.zeros(NX)
        alpha=np.zeros(NX)
	rho_0=1;
	lambd_0=1;
	#fp = open("prop.csv","w");
        #print("size alpha: %i"%np.size(alpha)) 

        if fluidType=="SF6":
	    rho,lambd,nu,kappa,alpha = fluid.SF6_NIST(T,p*1e-5);
        elif fluidType=="He":
	    rho,lambd,nu,kappa,alpha = fluid.He_NIST(T,p*1e-5);
        else:
            print("EROROR!\nFluid: "+fluidType+"unkown!\n")
            return

	return lambd,rho,alpha;

def get_RealCond(Tt,Tb,P,height,fluidType,database="NIST"):
    NX = 1000 
    print("Using: %s\n\n"%database)
    print("Height: %lf"%height)
    z = np.zeros(NX) 
    T = np.zeros(NX) 
    p = np.zeros(NX)
    pOld=np.zeros(NX)
    rho = np.zeros(NX)
    alpha = np.zeros(NX)
    lambd=np.zeros(NX)
    lambdOld=np.zeros(NX)
    h = height/(NX-1)
    
    q=1
    #Calculating starting parameters:
    P0 = P; # 
    Tm = (Tb+Tt)/2.;

    # Pobulating starting values:
    if fluidType=="SF6":
        if database=="C":
       	    rho_0,lambd_0,nu,kappa,alpha_0 = fluid.SF6_C(Tm,P0*1e-5);
        else:
            rho_0,lambd_0,nu,kappa,alpha_0 = fluid.SF6_NIST(Tm,P0*1e-5);
            print("Height: %.3lf (m)\tT: %lf %lf (degC)\t P: %lf (bar)\n"%(height,Tt,Tb,1e-5*P))  ;
    elif fluidType=="He":
        rho_0,lambd_0,nu,kappa,alpha_0 = fluid.He_NIST(Tm,P0*1e-5);
    else:
        print("ERROR!\nFluid "+fluidType+"is not in the database!")
    
    # Convert into SI units:
    for i in range(NX):
    	p[i] = P0;
    	lambd[i] = lambd_0;
    	rho[i] = rho_0;
    	z[i] = i*h;

    for k in range(6):
        T,q = get_temperature2(z,T,Tt,Tb,lambd,q);
        p = get_pressure(rho,p,P0,z);

    	if database=="C":
            lambd,rho,alpha=get_prop_C(T,p);
    	else:
	    lambd,rho,alpha=get_prop_NIST(fluidType,T,p);

        #Temperature profile (shooting method):
#        print("%i\t q: %.6lf [W/m^2]\n"%(k,q));	
 
    return q,z,T,p,rho,lambd



if __name__== "__main__":
    fluid.PrepNISTData("He",0,7,0.,6) 
    #qN,zN,TN,pN,rhoN,lambdN = get_RealCond(np.array(Tt),np.array(Tb),np.array(P),HEIGHT,'NIST')
    qN,zN,TN,pN,rhoN,lambdN = get_RealCond(Tt,Tb,P,HEIGHT,"SF6",'NIST')
