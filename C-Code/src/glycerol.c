/* On 04/24/00, modify from guenter/projects2/gasprop/acetone.c by Xiaochao Xu */
/* Ref. Physical and Thermodynamic properties of pure chemicals, Data Compilation. 
UCSB Lib# TP200 .D39 1989 */

#include "header.h"

glycerol(temp, press, prho, palfa, pcomp, plambda, pkappa, peta, pcp, ppsi, plewis)
double temp, press, *prho, *palfa, *pcomp, *plambda, *pkappa, *peta, *pcp, *ppsi, *plewis;
{
double tempa,temp_2,temp_3,conversion,MolWt,term,power;

	tempa = 273.15 +temp;
	temp_2 = tempa * tempa;
	temp_3 = temp_2 * tempa;
	MolWt = MW_GLYCEROL;              /*Molecular weight of glycerol (gm/mole) */
	conversion = 1.0e-03;     /* conversion factor for Kg/m^3 to gm/cm^3  */
/*
	Density Calculation (in gm/cm^3) 
*/
	term = (1-tempa/723.00);
	power=1+pow(term,0.15410);
	*prho = 0.94390/pow(0.24902,power);
	*prho *= MolWt * conversion;
/*
	Thermal expansion coefficient (1/K)
*/
	*palfa = pow((1-tempa/723.00),(0.15410-1.0))*(-0.15410*log(0.24902)/723.00);
/*	
	Heat Capacity (J/g K). Multiply by 10^7 to convert to erg/g K
*/
	*pcp = 6.8230e4 + 5.0520e2 * tempa;
	*pcp /= MolWt;
	*pcp *= conversion;
	*pcp *= 1.0e7;
/*
	Shear viscosity in gm/cm s
*/
	power = -237.03 + 1.6739e4/tempa +31.734*log(tempa);
	*peta = exp(power);
	*peta *= 10.0;
/*
	Thermal conductivity ( erg/(s cm K) )
*/
		*plambda = 2.5800e-01 + 1.1340e-04*tempa;
		*plambda *= 1.0e-02;
		*plambda *= 1.0e7;
/*
	Thermal diffusivity
*/
	*pkappa = ( *plambda / (*prho) / (*pcp) );

	return;

}
