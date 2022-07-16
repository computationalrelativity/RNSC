# Makefile for RNS

BASE=$(shell /bin/pwd)
SRCD=$(BASE)/src
INCD=$(BASE)/include
OBJD=$(BASE)/obj
LIBD=$(BASE)/lib


NAME=RNSRun
LIBNAME=RNS
EXE=$(NAME).x
SRC=$(wildcard $(SRCD)/*.c)
OBJ=$(patsubst $(SRCD)/%.c,$(OBJD)/%.o,$(SRC))
INC=$(wildcard $(INCD)/*.c)
INC_PARAMS=$(foreach d, $(INCD), -I$d)

OBJ2=$(wildcard $(OBJD)/*.o)
LIB=$(LIBD)/lib$(LIBNAME).so


# mandatory flags

CC = gcc
LD = ld
AR = ar

CFLAGS = -std=c99 -fPIC -pedantic
LFLAGS=
LDLFLAGS= -lm # -lgsl -lgslcblas

CFLAGS += -O3

all: $(EXE) $(LIB)
	@echo "All done"

$(EXE): $(OBJ)
	@echo "Building $@ ..."
	@echo

	$(CC) $(CFLAGS) $(LFLAGS) $(OBJ) -o $@ $(LDLFLAGS)

$(OBJD)/%.o: $(SRCD)/%.c
	@echo "Building objects ..."
	@echo

	@mkdir -p $(dir $@)
#	@echo "Compiling $< ... -> $@"
	$(CC) $(CFLAGS) $(INC_PARAMS) -c $< -o $@

$(LIB): $(OBJ)
	@echo "Making libraries... "
	@echo

	@mkdir -p $(LIBD)
	$(CC) $(CFLAGS) -shared $(OBJ) -o $(LIBD)/lib$(LIBNAME).so

	$(AR) rcs $(LIBD)/libRNS.a $(OBJ)
#	$(AR) rcs $(LIBD)/libRNS.a $@ $^

clean:
	@echo "Cleaning ..."
	@rm -rf $(OBJD)
	@rm -rf $(LIBD)
	@rm -rf $(EXE)
	@echo "... done"

.PHONY: all clean
