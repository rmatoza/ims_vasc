
# Makefile for IMS_VASC
# set fortran compiler and flags
FC = gfortran
FFLAG = -O3

# bin directory for executables
BINDIR=../bin

# subroutine object files to link
OBJS  =  sphere_subs.o \
datetime.o \


all: $(OBJS) ims_vasc ims_gridproc readwrite_imsvascbin_gridproc readwrite_imsvascbin_basic

# make clean
clean:
	\rm -f ${BINDIR}/* *.o *.mod

# make .o subroutine files (%.o: %.f)
%.o: %.f
	$(FC) $(FFLAG) -c $< -o $@

ims_vasc: ims_vasc.f90
	$(FC) $(FFLAG) ims_vasc.f90 $(OBJS) -o $(BINDIR)/ims_vasc

ims_gridproc: ims_gridproc.f90
	$(FC) $(FFLAG) ims_gridproc.f90 $(OBJS) -o $(BINDIR)/ims_gridproc

readwrite_imsvascbin_basic: readwrite_imsvascbin_basic.f90
	$(FC) $(FFLAG) readwrite_imsvascbin_basic.f90 -o $(BINDIR)/readwrite_imsvascbin_basic
	
readwrite_imsvascbin_gridproc: readwrite_imsvascbin_gridproc.f90
	$(FC) $(FFLAG) readwrite_imsvascbin_gridproc.f90 -o $(BINDIR)/readwrite_imsvascbin_gridproc


