all: RealConProfile 

S=./src
O=./obj

RealConProfile: $O/RealConProfile.o #$O/sf6.o $O/libheatcond.o
	cc $O/RealConProfile.o -lFluidPropC -ldl -lCoolProp -o RealConProfile -lm -lstdc++ 

$(O)/RealConProfile.o: $(S)/RealConProfile.c
	mkdir -p $(O)
	cc -c $S/RealConProfile.c -fPIE
	mv RealConProfile.o ${O}/

install:
	cp RealConProfile ${HOME}/bin
