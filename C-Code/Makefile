all: RealConProfile 

S=./src
O=./obj

RealConProfile: $O/RealConProfile.o #$O/sf6.o $O/libheatcond.o
	cc $O/RealConProfile.o -lFluidPropC -ldl -lCoolProp -o RealConProfile -lm -lstdc++ 
	#cc -o RealConProfile ${O}/RealConProfile.o ${O}/libheatcond.o ${O}/sf6.o -lm

#${O}/libheatcond.o: $(S)/libheatcond.c
#	cc -c $S/libheatcond.c
#	mv libheatcond.o $O/


$(O)/RealConProfile.o: $(S)/RealConProfile.c
	cc -c $S/RealConProfile.c -fPIE
	mv RealConProfile.o ${O}/

#$(O)/sf6.o: $(S)/sf6.c
#	cc -c $S/sf6.c -fPIE 
#	mv sf6.o ${O}/
