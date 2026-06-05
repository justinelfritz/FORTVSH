NAME = vsh_test.exe

BASE = $(HOME)/Desktop/magfric/fortran_VSH
EXEC = $(BASE)/$(NAME)
OBJD = $(BASE)/obj
SRCD = $(BASE)/src

FC = gfortran -O3 -I $(OBJD)

VPATH = $(OBJD):$(SRCD)

OBJS = kinds.o globals.o vsh.o tests.o main.o

$(NAME): $(OBJS)
	@echo Building executable
	cd $(OBJD); $(FC) -o $(EXEC) $(OBJS)

.f.o:
	@echo Building $*.o
	cd $(OBJD); $(FC) -J $(OBJD) -c $<

clean:
	rm -f obj/*.o obj/*.mod src/*.mod validation/*.dat $(EXEC)
