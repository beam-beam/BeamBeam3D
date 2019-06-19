#*****************************************************
#  General Makefile
#
#*****************************************************

#**************************************************************************
# Macros defining the Fortran, C/C++ compiler and linker.

CC = ftn
LINK = ftn
FFLAGS = -O3

#**************************************************************************
# List of .o files that EXENAME depends on.  Edit as appropriate for MP.

OBJS = \
	Timer.o Pgrid.o FFT.o Commun.o CompDom.o Utility.o Input.o Output.o Distribution.o \
	Dampfluc.o Transfer.o Transpose.o Orbit.o Linearmap.o Extmap2nd.o \
	Extmap3rd.o Extmap4th.o DepoScat.o FieldSolver.o Beambeam.o \
	Feedback.o AccSimulator.o main.o

#**************************************************************************
# Change this line if you don't like 'a.out'.

EXENAME = xmain

#************************************************************************
# disable predefined suffixes and define your own set of allowed suffixes
 .SUFFIXES:
 .SUFFIXES: .o .f .F .c .f90 .F90

#*************************************************************************
# inference rules (how to compile object files that have no explicit rules)
#  $* = name part of target
#  $@ = full target name
#  $< = dependent name

.f90.o:
	$(CC) -c $(FFLAGS) $<

#**************************************************************************
# Rules for building EXENAME from OBJS and OBJS from your source.

$(EXENAME):  $(OBJS) 
	$(LINK) -o $(EXENAME) $(OBJS) 

#************************************************************************
# if you wish to compile a certain object with different flags
# or in some special way, then specify the target & dependency explicitly
# the following line say Timer.o is depended on Timer.f90
#Timer.o: Timer.f90
#	$(CC) -c -O3 Timer.f90

#***********************************************************************
clean:
	@if [ -f $(EXENAME) ]; then rm $(EXENAME); fi
	@for file in `ls *.o *.mod`; do rm "$$file"; done

