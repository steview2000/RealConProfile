/* On 12/30/99, modified from guenter/projects2/gasprop/acetone.c by Xiaochao Xu */
/* Ref.: Physical and Thermodynamic properties of pure chemicals, 
Data Compilation. UCSB Lib# TP200 .D39 1989 */

#include "header.h"

ethanol(temp, press, prho, palfa, pcomp, plambda, pkappa, peta, pcp, ppsi, plewis)
double temp, press, *prho, *palfa, *pcomp, *plambda, *pkappa, *peta, *pcp, *ppsi, *plewis;
{
double tempa,temp_2,temp_3,conversion,MolWt,term,power;

	tempa = 273.15 +temp;
	temp_2 = tempa * tempa;
	temp_3 = temp_2 * tempa;
	MolWt = MW_ETHANOL;
	conversion = 1.0e-03;     /* Kg/m^3 to gm/cm^3  */
/*
	Density Calculation (in gm/cm^3) 
*/
	term = (1-tempa/516.25);
	power=1+pow(term,0.23670);
	*prho = 1.5223/pow(0.26395,power);
	*prho *= MolWt * conversion;
/*
	Thermal expansion coefficient (1/K)
*/	
	*palfa = 6.1072e-4/pow((1.0-0.00193705*tempa),0.76330);
/*	
	Heat Capacity (J/g K). Multiply by 10^7 to convert to erg/g K
*/
	*pcp = 9.4560e4 - 5.6200e1 * tempa -3.2900e-01 * temp_2 + 2.3980e-3 * temp_3;
	*pcp /= MolWt;
	*pcp *= conversion;
	*pcp *= 1.0e7;
/*
	Shear viscosity in gm/cm s
*/
	power = 8.0490 + 7.7600e2/tempa -3.0680*log(tempa);
	*peta = exp(power);
	*peta *= 10.0;
/*
	Thermal conductivity ( erg/(s cm K) )
*/
	*plambda = 2.5300e-01 - 2.8100e-04*tempa;
	*plambda *= 1.0e-02;
	*plambda *= 1.0e7;
/*
	Thermal diffusivity
*/
	*pkappa = ( *plambda / (*prho) / (*pcp) );

	return;

}
