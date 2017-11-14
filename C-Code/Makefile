all: RealConProfile RealConProfileSingle RealConProfileSingle2

S=./src
O=./obj


RealConProfile: $O/RealConProfile.o
	cc -o RealConProfile RealConProfile.o gnuplot_i.o sf6.o -lm

RealConProfileSingle: $O/RealConProfileSingle.o
	cc -o RealConProfileSingle RealConProfileSingle.o gnuplot_i.o sf6.o -lm

RealConProfileSingle2: $(O)/RealConProfileSingle2.o
	cc -o RealConProfileSingle2 $(O)/RealConProfileSingle2.o gnuplot_i.o sf6.o -lm

RealConProfile.o: $(S)/RealConProfile.c
	cc -c $S/RealConProfile.c -o $O/RealConProfile.o

RealConProfileSingle.o: $(S)/RealConProfileSingle.c
	cc -c $S/RealConProfileSingle.c

$(O)/RealConProfileSingle2.o: $(S)/RealConProfileSingle2.c
	cc -c $S/RealConProfileSingle2.c -o ${O}/RealConProfileSingle2.o

$(O)/RealConProfile.o: $(S)/RealConProfile.c
	cc -c $S/RealConProfile.c -o ${O}/RealConProfile.o
