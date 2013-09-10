#Makefile for parallel compiling
#=====================================================
# Note: use -C=all to catch array out of bounds errors
#============================================================================ 
# Compiler options
#============================================================================
FC = ifort -O3 -fpp -warn all -g -pg

LIBPATH=./cpp
#LIBPATH=./ispc
#LIBPATH=./cuda

ifeq ($(AVX), 1)
	FC += -xavx
else 
	FC += -xsse
endif

ifeq ($(SINGLE),1)
	FC += -DSINGLE
else
	FC += -r8
endif

ifeq ($(OMP),1)
	FC += -openmp
endif

ifeq ($(FORTRAN),1)
	FC += -D_USEFORTRAN_
endif

#=======================================================================
# Library
#======================================================================
#LINKS = /usr/lib64/liblapack.so.3 /usr/lib64/libblas.so.3 \
	-L$(FFTWLIB) -lfftw3 -lfftw3f

ifeq ($(SINGLE),1)
LINKS = -lfftw3f -lblas -llapack
else
LINKS = -lfftw3  -lblas -llapack
endif
LINKS += -lstdc++
# LINKS += -lpapi -L$(HOME)/usr/lib
ifeq ($(OMP),1)
	LINKS += -openmp
endif
LINKS += -L$(LIBPATH) -lrbflow
#============================================================================ 
# make PROGRAM   
#============================================================================


PROGRAM = boutnp

SRC = fortran/cfl.f \
			 fortran/coetar.f \
			 fortran/continua.f \
			 fortran/cordin.f \
			 fortran/densbo.f \
			 fortran/densmc.f \
			 fortran/divg.f \
			 fortran/divgck.f \
			 fortran/gcurv.f \
			 fortran/hdnl2.f \
			 fortran/hdnl3.f \
			 fortran/hdnlro.f \
			 fortran/indic.f \
			 fortran/inirea.f \
			 fortran/initia.f \
			 fortran/inqpr.f \
			 fortran/invtr2.f \
			 fortran/invtr3.f \
			 fortran/invtrro.f \
			 fortran/meshes.f \
			 fortran/ntrvb_line.f \
			 fortran/openfi.f \
			 fortran/outh.f \
			 fortran/outpf.f \
			 fortran/papero.f \
			 fortran/param.f \
			 fortran/phcalcTP.f \
			 fortran/prcalc.f \
			 fortran/radialgrid.f \
			 fortran/solq2j.f \
			 fortran/solq2k.f \
			 fortran/solq3j.f \
			 fortran/solq3k.f \
			 fortran/solroj.f \
			 fortran/solrok.f \
			 fortran/tripvmy_line.f \
			 fortran/tsch.f \
			 fortran/updvp.f \
			 fortran/vmaxv.f \
			 fortran/wson.f \
			 fortran/phiniTP_intel.f  


OBJ = $(SRC:%.f=%.o) 

MODULES = fortran/param.o

all: $(MODULES) $(PROGRAM)
#============================================================================ 
# Linking    
#============================================================================

$(PROGRAM) : $(OBJ)
#	make -C $(LIBPATH) 
	$(FC) $(OBJ) $(LINKS) -o $@ 

#============================================================================
#  Dependencies
#============================================================================

param.o: fortran/param.f
	$(FC) -c $< -o $@ 

%.o:   fortran/%.f 
	$(FC) -c $< -o $@
        

#============================================================================
#  Clean up
#============================================================================

clean:
	/bin/rm -rf  $(PROGRAM) fortran/*.o  cpp/*.o *.mod *__genmod.f90 *__genmod.mod

clean_all: clean
	/bin/rm -rf *~ fortran/*~ cpp/*~ fortran/.??* cpp/.??* *.out
	ls -d .??*|grep -v '.git'|xargs /bin/rm -rf 
	make -C $(LIBPATH) clean_all
	


$(OBJ): $(MODULES) ./dims.h
