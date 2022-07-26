# RealConProfile

This code calculates the conductive heat transport between two plates within a fluid of varying heat conductivity (non-Oberbeck-Boussinesq
 as descriped in  O. Shishkina, S. Weiss, and E. Bodenschatz, , Phys. Rev. Fluids, 1(6), 062301(R)

## Installation

1. Make sure libFluidPropC is installed (which itself depends on CoolProp)
	

3. Clone this repository: 

```
git clone https://github.com/steview2000/RealConProfile.git
```

3. Enter directory and compile:

```
cd RealConProfile 
make
```	
4. Install in ${HOME}/bin
```
make install
```

## Dependencies:
	* libFluidPropC (https://github.com/steview2000/libFluidProp.git)

## Usage: 
Usage: ./RealConProfile \[option\] \[value\] 
Calculates the conductive heat flux through a fluid of varying heat conductivity
Options are optional and don't have to be provided.
        -h       prints this help message
        -f       <fluids>        defines the fluid
                 Possible fluids:
                        Air     
                        Hydrogen
                        Helium  
                        Nitrogen
                        CO2     
                        Xenon   
                        SF6     
                        Ethane  
                        Water   
                        Acetone 
                        Methanol
                        Ethanol 

        -b <value> bottom plate temperature in degree Celsius
        -t <value> top plate temperature in degree Celsius
        -P <value> pressure in bar 
        -H <value> height of the cell in meter 

## Example
The following code calculates the heat flux for SF6 at P=10bar with bottom plate temperate 15 C and
top plate temperature 35 C in a 2m tall cell:
```
./RealConProfile -f SF6 -b 15 -t 35 -P 10 -H 2
```

