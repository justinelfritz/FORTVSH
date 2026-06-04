NAME = vsh_test.exe

BASE = $(HOME)/Desktop/magfric/fortran_VSH
EXEC = $(BASE)/$(NAME)
OBJD = $(BASE)/obj
SRCD = $(BASE)/src

F77 = gfortran -O -I $(SRCD) $(SHTOOLSMODPATH) -L $(SHTOOLSLIBPATH) -lSHTOOLS -lfftw3 -lm -llapack -lblas -O3 -o M

FFLAGS = 

VPATH = $(OBJD):$(SRCD)

OBJS = main.o globals.o kinds.o vsh.o tests.o

$(NAME): $(OBJS)
	@echo Building executable
	cd $(OBJD); $(F77) -o $(EXEC) $(OBJS)

.f.o:
	@echo Building $*.o 
	cd $(OBJD); $(F77) -c $<

clean:
	rm -f terminal* makefile~ *.dat~ obj/*.o obj/*.mod $(EXEC)
