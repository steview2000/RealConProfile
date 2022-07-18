/* ****************************************************************************
    gasindex.c

    This program estimates the refractive index of three gases (He, CO2, SF6)
and two gas mixtures (He-CO2, He-SF6).  The refractivity virial coefficients
of pure gases are obtained from the following papers:

   Proc. R. Soc. Lond. A. 336, 275 (1974)
   J. Chem. Phys. 71, 4951 (1979)
   J. Chem. Phys. 94, 5669 (1991)

The temperature is assumed to be 20 - 50 C.

   The definition of refractivity virial coefficients for pure gases is:

      Rm = (n^2-1)/((n^2+2)*rho) = Ar + Br*rho + ...

where rho is the gas density in _moles_ per unit volume (unit: mol/m^3).
It is easy to get n from this equation for pure gases.

For gas mixtures, we assume that (n - 1) is additive, therefore, we calculate
(n - 1) for two components respectively and then add them together to get
(n_mix - 1).


by Jun Liu, 3/17/95.

modified to include H2-Xe mixture
by Kapil Bajaj 3/14/00
For Xenon: H. J. Achtermann et. al J. Chem Phys 98, 1 Feb 1993, pp2307
For Hydrogen: Achtermann et. al., J. Chem Phys. 94, 15 April 1991, pp 5674

*************************************************************************** */

#include <stdio.h>
#include <math.h>
#define AR_HE	5.2e-7      /* unit:  m^3/mol   */
#define AR_CO2	6.65e-6
#define AR_SF6	1.134e-5
#define BR_HE	-6.0e-14    /* unit:  m^6/mol^2   */
#define BR_CO2	3.0e-12
#define BR_SF6	2.8e-11
#define AR_XE   1.036e-5 /* at 298 K */	
#define BR_XE	2.58e-11 /* at 298 K deg C */
#define AR_H2	2.0713e-6      /* unit:  m^3/mol   */
#define BR_H2	1.3e-13      /* unit:  m^6/mol^2   */
#define MW_HE 	4.0026
#define MW_CO2 	44.010
#define MW_SF6	146.06
#define MW_H2   2.0156	
#define MW_XE   131.29	

char *prgname;
void main(argc, argv)
int argc;
char *argv[];
{
	int mixflag, gasflag=0, biflag=0;
        double ind, a, b, middle, middle2, ind2;
        double rho, cmhe, c, rho1, rho2;
	double sqrt();

	fprintf(stderr, "Pure gas (enter 0) or mixture (enter 1)? \n");
	scanf("%d", &mixflag);

	if(mixflag == 0)  
	  {
	   fprintf(stderr, "Which gas? (He=1, CO2=2, SF6=3)\n");
		scanf("%d", &gasflag);
	  }
	else
          {
	   fprintf(stderr, "Which mixture? (He+CO2=1, He+SF6=2, H2-Xe=3)\n");
		scanf("%d", &biflag);
          }

	fprintf(stderr, "Gas density (unit: kg/m^3)? \n");
	scanf("%lf", &rho);
        
	rho = rho*1000.0;   	/* convert to g/m^3  */

/*   Pure gas   */

	if (gasflag !=0)
         {
          if (gasflag == 1)
            {  
	    	a = AR_HE;
		b = BR_HE;
		rho = rho/MW_HE;	/* molar density */
	    }
          else if (gasflag == 2)
            {
  		a = AR_CO2;
		b = BR_CO2;
		rho = rho/MW_CO2;
	    }
          else if (gasflag == 3)
            {
  		a = AR_SF6;
		b = BR_SF6;
		rho = rho/MW_SF6;
	    }
	  
          middle = a*rho + b*rho*rho;
          ind = sqrt( 3/(1-middle) - 2 );
         }

/*  Gas mixture  */

	if (biflag != 0) 
	 {
	  if (biflag == 1)
       	    {  
	  fprintf(stderr, "Molar concentration of He? \n");
	  scanf("%lf", &cmhe);
	     	a = AR_CO2;
		b = BR_CO2;
                c = cmhe*MW_HE/(cmhe*MW_HE + (1 - cmhe)*MW_CO2);  /*  mass 
						concentration of He */
		rho2 = c*rho/MW_HE;		/*  molar density of He */
		rho1 = (1-c)*rho/MW_CO2;        /*  molar density of CO2 */ 
            }
          else
	  if (biflag == 2)
            {
	  fprintf(stderr, "Molar concentration of He? \n");
	  scanf("%lf", &cmhe);
 		a = AR_SF6;
		b = BR_SF6;
                c = cmhe*MW_HE/(cmhe*MW_HE + (1 - cmhe)*MW_SF6);  /*  mass 
						concentration of He */
		rho2 = c*rho/MW_HE;		/*  molar density of He */
		rho1 = (1-c)*rho/MW_SF6;        /*  molar density of CO2 */ 
            }
          else
            {
	  fprintf(stderr, "Molar concentration of H2? \n");
	  scanf("%lf", &cmhe);
 		a = AR_XE;
		b = BR_XE;
                c = cmhe*MW_H2/(cmhe*MW_H2 + (1 - cmhe)*MW_XE);  /*  mass 
						concentration of H2 */
		rho2 = c*rho/MW_H2;		/*  molar density of H2 */
		rho1 = (1-c)*rho/MW_XE;        /*  molar density of Xe */ 
	}
           middle = a*rho1 + b*rho1*rho1;
	   ind = sqrt( 3/(1-middle) - 2 );

           middle2 = AR_H2*rho2 + BR_H2*rho2*rho2;
	   ind2 = sqrt( 3/(1-middle2) - 2 );

           ind = ind + ind2 - 1;
         }

 	fprintf(stdout, "The refractive index you ask for is --> %8g\n", ind);
}
          
