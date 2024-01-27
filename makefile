
FC = gfortran
RM = rm -f

#FFLAGS = -O3
FFLAGS = -g


OBJS = main_box.o  problem.o

all:  $(OBJS) 
	$(FC) -o dfbox $(OBJS)

.SUFFIXES : .f90 .o

.f90.o: $* parameter.f ; $(FC) $(FFLAGS) -c $*.f90

clean: 
	$(RM) *.o
	$(RM) fort.*
	$(RM) solution.dat
	$(RM) dfbox

