#!/usr/bin/python3
# This program calculates the real conducted heat in RBC, by integrating the heat equation (shooting method) taken into
# consideration different pressures (hydrostatic) and the fluid properties at any vertical position
#
# 12.2.2016: I will make the following changes
#
# see our paper: Shishkina, Weiss, Bodenschatz, PRF, 1, 062301 (R)

import sys
import numpy as np
import PyFluidProp as fluid

#NX= 1000

g=9.81; # gravity (in m/s^2)


def get_temperature(z,Tt,Tb,lambd):
	"""
	Calculates the vertical temperature profile T(z) at a given, top plate temperature (Tt), bottom
	plate temperature (Tb), profile of the thermal conductivity (lambda)

	Usage:
	------
	T,q = get_temperature2(z,Tt,Tb,lambd)
	
	Arguments:
	---------
		z     ... vertical co-ordinates
		Tt    ... top plate temperature
		Tb    ... bottom plate temperature
		lamba ... vertical profile of the heat conductivity

	Returns:
	--------
		T ... vertical temperature profile
		q ... conductive heat flux
	"""

	NX            = np.size(z)
	HEIGHT        = z[-1]-z[0]
	int_inf_lambd = np.zeros(NX)
	q_array       = np.zeros(NX);
	T			  = np.zeros(NX)

	h = HEIGHT/(NX-1);
	
	dT = Tb-Tt;
	int_inf_lambd[0] = lambd[0];
	for i in range(NX-1):
		int_inf_lambd[i+1] = int_inf_lambd[i]+h/lambd[i];
	
	for i in range(NX):
		T[i] = Tb - dT*int_inf_lambd[i]/int_inf_lambd[NX-1];
	#	if (i>0):
	#		q_array[i] = -lambd[i]*(T[i]-T[i-1])/h;
	
	return T,(Tb-Tt)/int_inf_lambd[NX-1];

def get_temperature2(z,T,Tt,Tb,lambd,q):
	"""	
	Calculates the temperature profile based on a z-array using an iterative approach method
	"""

	NX     = np.size(z)
	HEIGHT = z[-1]-z[0]
	h      = 1.*HEIGHT/(NX-1);
	T[0]   = Tb;
	err    = 1;
	TDelta = 0.001;
	count  = 0;

	while(err>0.0001):
		count += 1;
		low    = 0;
		high   = 0;
		T_tOld = T[NX-1];
		T[0]   = Tb;
		
		for i in range(NX-1):
			k1     = -q * h/lambd[i];
			T[i+1] = (T[i]+k1);
			# start value
	
		Tmax=0;
	
		# now check how close T[NX-1] is to T_t
		err=np.sqrt((T[NX-1]-Tt)*(T[NX-1]-Tt));
		
		TDeltaOld = TDelta;

		if(T[NX-1]==T_tOld):
			TDelta = -0.1*(T[NX-1]-Tt);
		else:
			TDelta= -0.1*(T[NX-1]-Tt)*(TDeltaOld)/(T[NX-1]-T_tOld);
		
		q = q + lambd[0]*TDelta;
	
		T_tOld = T[NX-1];

	return T,q;

def get_pressure(rho,p,Pin,z):
	"""
	This function calculates the hydrostatic pressure (in pa)
	"""
	NX = np.size(z)
	h = (z[-1]-z[0])/(NX-1);

	p[0]=0;
	for i in range(NX-1):
		p[i+1] = p[i]-g*h*rho[i]	;
   
	p0 = Pin-p[int(NX/2)];
	
	for i in range(NX):
		p[i] = p[i]+ p0	
	
	return p;


def get_prop_C(T,p):
	"""
	This function calculates fluid properties based on the NIST data at every vertical position
	"""

	NX = np.size(T)
	rho   = np.zeros(NX)
	lambd = np.zeros(NX)
	alpha = np.zeros(NX)
	rho_0   = 1;
	lambd_0 = 1;

	for i in range(NX):
		temp = T[i];
		press = p[i]; # pa -> bar
		rho_0,lambd_0,nu,kappa,alpha_0 = fluid.SF6_CoolProp(temp,press*1e-5);
		rho[i] = rho_0;
		lambd[i] = lambd_0;
		alpha[i] = alpha_0;
	
	return lambd,rho,alpha;

def get_prop_NIST(fluidType,T,p):
	"""
	This function calculates fluid properties based on the NIST data at every vertical position
	"""

	NX      = np.size(T)
	rho     = np.zeros(NX)
	lambd   = np.zeros(NX)
	alpha   = np.zeros(NX)
	rho_0   = 1;
	lambd_0 = 1;
	
	if fluidType=="SF6":
		rho,lambd,nu,kappa,alpha = fluid.SF6_NIST(T,p*1e-5);
	elif fluidType=="He":
		rho,lambd,nu,kappa,alpha = fluid.He_NIST(T,p*1e-5);
	elif fluidType=="N2":
		rho,lambd,nu,kappa,alpha = fluid.N2_NIST(T,p*1e-5);
	else:
	    print("EROROR!\nFluid: "+fluidType+" unkown!\n")
	    return

	return lambd,rho,alpha;

def get_RealCond(Tt,Tb,P,height,fluidType,database="NIST"):
	"""
	Calculates the conductive heat transport (important to calculate the Nusselt number.

	q,z,T,p,rho,lambd = get_RealCond(Tt,Tb,P,height,fluidType,database=\"NIST\")
	
	Arguments:
	----------
		Tt			... top plate temperature (in degC)
		Tb 			... bottom plate temperature (in degC)
		P  			... pressure (in pa)
		height  	... cell height (in m)
		fluidType	... fluid (currently only He and SF6 available)
		database	... can be either \"NIST\" for the NIST database or \"C\" for Guenters C-code
	
	Returns:
	--------
		q		... conductive heat flux (in W)
		z 		... array with the vertical coordinates
		T 		... vertical temperature profile
		p 		... vertical pressure profile (in pa)
		rho		... vertical density profile
		lambd	... vertical profile of the heat conductivity
	
	"""
	
	NX = 1000  # spatial resolution in z-direction (number of sampling points)
	print("Using: %s\n\n"%database)
	print("Height: %lf"%height)
	print("Fluid: "+fluidType)

	z		= np.zeros(NX) # z-coordinates
	T		= np.zeros(NX) # array with vertical temperature values
	p		= np.zeros(NX) # array with vertical pressure values 
	pOld    = np.zeros(NX)
	rho	    = np.zeros(NX) # array with vertical density distribution
	alpha	= np.zeros(NX) # array with veritcal expansivity
	lambd	= np.zeros(NX) # array with therm conductivity
	lambdOld = np.zeros(NX)
	h	  	= height/(NX-1) # vertical steps
	
	q=1

	#Calculating starting parameters:
	P0 = P; # 
	Tm = (Tb+Tt)/2.;

	# Pobulating starting values:
	if fluidType=="SF6":
		if database=="C":
	   		rho_0,lambd_0,nu,kappa,alpha_0 = fluid.SF6_CoolProp(Tm,P0*1e-5);
		else:
			rho_0,lambd_0,nu,kappa,alpha_0 = fluid.SF6_NIST(Tm,P0*1e-5);
			print("Height: %.3lf (m)\tT: %lf %lf (degC)\t P: %lf (bar)\n"%(height,Tt,Tb,1e-5*P))  ;
	elif fluidType=="He":
		rho_0,lambd_0,nu,kappa,alpha_0 = fluid.He_NIST(Tm,P0*1e-5);
		print("Height: %.3lf (m)\tT: %lf %lf (degC)\t P: %lf (bar)\n"%(height,Tt,Tb,1e-5*P))  ;
	elif fluidType=="N2":
		rho_0,lambd_0,nu,kappa,alpha_0 = fluid.N2_NIST(Tm,P0*1e-5);
		print("Height: %.3lf (m)\tT: %lf %lf (degC)\t P: %lf (bar)\n"%(height,Tt,Tb,1e-5*P))  ;
	else:
		print("ERROR!\nFluid "+fluidType+"is not in the database!")
	
	# Convert into SI units:
	for i in range(NX):
		p[i] = P0;
		lambd[i] = lambd_0;
		rho[i] = rho_0;
		z[i] = i*h;

	# now iterate (6 times is usually enough)
	for k in range(6):
		T,q = get_temperature(z,Tt,Tb,lambd);
		p = get_pressure(rho,p,P0,z);

		if database=="C":
			lambd,rho,alpha = get_prop_C(T,p);
		else:
			lambd,rho,alpha = get_prop_NIST(fluidType,T,p);

	return q,z,T,p,rho,lambd



if __name__== "__main__":
    #fluid.PrepNISTData("He",0,7,0.,6) 
    #qN,zN,TN,pN,rhoN,lambdN = get_RealCond(np.array(Tt),np.array(Tb),np.array(P),HEIGHT,'NIST')
	#HEIGHT= 2.24 # in m for the G=1 cell HPCF
	#Tb = 30.     #(21.431+11.858/2) # in degC
	#Tt = 18.     #(21.431-11.858/2)	# in degC
	#P  = 8.032e5 # in pa (at half-height)
	HEIGHT = float(input("Height: "))
	P = float(input("Press: "))*1e5
	Tb =float(input("Tb: "))
	Tt =float(input("Tt: "))
	if sys.argv[1]=='SF6':
		qN,zN,TN,pN,rhoN,lambdN = get_RealCond(Tt,Tb,P,HEIGHT,"SF6",'NIST')
	elif sys.argv[1]=='Helium':
		qN,zN,TN,pN,rhoN,lambdN = get_RealCond(Tt,Tb,P,HEIGHT,"He",'NIST')
	else:
		print(sys.argv[1]+" not defined!")
		sys.exit(1)

	print("q[W/m^2]: ",qN)
	#header_string = "z,Tm,p,rho,lambda"
	#np.savetxt('testout-py.csv',np.array([zN,TN,pN,rhoN,lambdN]).transpose(),header=header_string,fmt="%.4g",delimiter=",")

	
