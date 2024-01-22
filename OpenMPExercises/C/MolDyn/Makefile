# Makefile for the molecular dynamics code

#
# C compiler and options for Intel
#
#CC=     icc -O3 -qopenmp -std=c99
#LIB=    -lm

#
# C compiler and options for GNU 
#
CC=     gcc -O3 -fopenmp 
LIB=	-lm
CPPFLAGS = -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp

#
# Object files
#
OBJ=    main.o \
	dfill.o \
	domove.o \
        dscal.o \
	fcc.o \
	forces.o \
	mkekin.o \
	mxwell.o \
	prnout.o \
	velavg.o

#
# Compile
#
md:	$(OBJ)
	$(CC)  -o $@ $(OBJ)  $(LIB) $(CPPFLAGS)

.c.o:
	$(CC) -c $< $(CPPFLAGS)

#
# Clean out object files and the executable.
#
clean:
	rm *.o md
