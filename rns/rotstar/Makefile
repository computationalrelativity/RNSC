#######################################################################
#
#                     MAKEFILE FOR ROTSTAR.C
#
#######################################################################

#######################################################################
#                     COMPILING FLAGS
#######################################################################

# Specify C compiler and ANSI option:
CC=gcc

# DEBUGGING OPTION
#CFLAGS =-g3 -Wall

# OPTIMIZATION
#CFLAGS =-O1 
CFLAGS =-O2 -fPIE
#CFLAGS =-O3 -fPIE

#######################################################################
#                     SPECIFY GRID SIZE
#######################################################################

#STANDARD
#SIZE=-DMDIV=65 -DSDIV=129

#HIGH
#SIZE=-DMDIV=101 -DSDIV=201

#VERY HIGH
#SIZE=-DMDIV=151 -DSDIV=301

#VERY VERY HIGH
#SIZE=-DMDIV=201 -DSDIV=401

#LOW
#SIZE=-DMDIV=51 -DSDIV=101

#VERY LOW
#SIZE=-DMDIV=41 -DSDIV=71

# USED IN WHISKY:
SIZE=-DMDIV=301 -DSDIV=601

#######################################################################
#                       OBJECT AND TARGETS 
#######################################################################

OBJ=equil.o equil_util.o nrutil.o rotstar.o 
EXE=rotstar.x

rotstar: $(OBJ)
	$(CC) $(CFLAGS) $(SIZE) -o $(EXE) $(OBJ) -lm 

rotstar.o: equil.h  nrutil.h consts.h rotstar.c
	$(CC) -c $(CFLAGS) $(SIZE) rotstar.c 

kepler: kepler.o equil.o equil_util.o nrutil.o 
	$(CC) $(CFLAGS) -lm  $(SIZE) -o kepler.x kepler.o equil.o equil_util.o nrutil.o -lm 

kepler.o: equil.h  nrutil.h consts.h kepler.c
	$(CC) -c $(CFLAGS) $(SIZE) kepler.c 

equil.o: equil.h equil_util.h nrutil.h consts.h equil.c
	$(CC) -c $(CFLAGS) $(SIZE) equil.c 

equil_util.o:equil_util.h nrutil.h consts.h equil_util.c
	$(CC) -c $(CFLAGS) $(SIZE) equil_util.c 

nrutil.o:nrutil.h nrutil.c
	$(CC) -c $(CFLAGS) $(SIZE) nrutil.c

clean:
	rm -rf $(EXE) kepler.x *.o *~

# 
