#------------------------------------------------------------------------------
ifeq ($(TARGET_MIC),)
CXX=nvcc
LD=g++
else
# Intel Xeon Phi/MIC requires using Intel C++ Compiler (ICC)
CXX=icpc
LD=icpc
CXXFLAGS=-mmic -x c++
endif

CXXFLAGS += -O3
DEFINEFLAGS = -DDUMMY=dummy

UNAME=$(shell uname)
ifeq ($(UNAME), Darwin)
CXXFLAGS+=-m64
endif

ifneq ($(CUDAPRINT),)
DEFINEFLAGS += -DCUDAPRINT=yes
endif

ifneq ($(PRINTCALLS),)
DEFINEFLAGS += -DPRINTCALLS=yes
endif

ifneq ($(PROFILE),)
DEFINEFLAGS += -DPROFILING=yes
endif

ifeq ($(TARGET_OMP),)
# nvcc (CUDA)
CXXFLAGS += -arch=sm_35
else
# OpenMP common flags
DEFINEFLAGS += -fno-inline -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_BACKEND_OMP
ifeq ($(TARGET_MIC),)
# GCC/Clang
DEFINEFLAGS += -fopenmp
else
# Intel C++ Compiler (ICC)
DEFINEFLAGS += -openmp
endif
endif

ifeq ($(CUDALOCATION), )
CUDALOCATION = /usr/local/cuda-7.5/bin
#CUDALOCATION = /cmshome/adrianodif/CUDAs/cuda-6.0/bin
endif
CUDAHEADERS = $(CUDALOCATION)/include/

PWD = $(shell /bin/pwd)
SRCDIR = $(PWD)/PDFs

INCLUDES += -I$(SRCDIR) -I$(PWD) -I$(CUDAHEADERS) -I$(ROOTSYS)/include -I$(PWD)/rootstuff 

# GooPdf must be first in CUDAglob, as it defines global variables.
FUNCTORLIST    = $(SRCDIR)/GooPdf.cu
FUNCTORLIST    += $(filter-out $(SRCDIR)/GooPdf.cu, $(wildcard $(SRCDIR)/*Pdf.cu))
FUNCTORLIST   += $(wildcard $(SRCDIR)/*Aux.cu)
HEADERLIST     = $(patsubst %.cu,%.hh,$(SRCFILES))
WRKFUNCTORLIST = $(patsubst $(SRCDIR)/%.cu,wrkdir/%.cu,$(FUNCTORLIST))

#NB, the above are used in the SRCDIR Makefile.

THRUSTO		= wrkdir/Variable.o wrkdir/FitManager.o wrkdir/GooPdfCUDA.o wrkdir/Faddeeva.o wrkdir/FitControl.o wrkdir/PdfBase.o wrkdir/DataSet.o wrkdir/BinnedDataSet.o wrkdir/UnbinnedDataSet.o wrkdir/FunctorWriter.o
ROOTRIPDIR	= $(PWD)/rootstuff
ROOTRIPOBJS	= $(ROOTRIPDIR)/TMinuit.o $(ROOTRIPDIR)/TRandom.o $(ROOTRIPDIR)/TRandom3.o #$(ROOTRIPDIR)/TVirtualFitter.o
ROOTUTILLIB	= $(ROOTRIPDIR)/libRootUtils.so

.SUFFIXES:
.PHONY:		goofit clean

goofit:		$(THRUSTO) $(ROOTUTILLIB)
		@echo "Built GooFit objects"

# One rule for GooFit objects.
wrkdir/%.o:	%.cc %.hh
		@mkdir -p wrkdir
		$(CXX) $(INCLUDES) $(CXXFLAGS) $(DEFINEFLAGS) -c -o $@ $<

# A different rule for user-level objects. Notice ROOT_INCLUDES.
%.o:	%.cu
	$(CXX) $(INCLUDES) $(ROOT_INCLUDES) $(DEFINEFLAGS) $(CXXFLAGS) -c -o $@ $<

# Still a third rule for the ROOT objects - these have their own Makefile.
$(ROOTRIPDIR)/%.o:	$(ROOTRIPDIR)/%.cc
			rm -f $@
			@echo "Postponing $@ for separate Makefile"

$(ROOTUTILLIB):	$(ROOTRIPOBJS)
		@cd rootstuff; $(MAKE)

include $(SRCDIR)/Makefile

FitManager.o:		FitManager.cc FitManager.hh wrkdir/ThrustFitManagerCUDA.o Variable.o
			$(CXX) $(DEFINEFLAGS) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

##FitManager1.o:	FitManager1.cc FitManager1.hh wrkdir/ThrustFitManagerCUDA.o Variable.o
##	$(CXX) $(DEFINEFLAGS) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

wrkdir/GooPdfCUDA.o:	wrkdir/CUDAglob.cu PdfBase.cu
			$(CXX) $(CXXFLAGS) $(INCLUDES) -I. $(DEFINEFLAGS) -c $< -o $@
			@echo "$@ done"

##wrkdir/GooPdfCUDA1.o:	wrkdir/CUDAglob1.cu PdfBase1.cu
##			$(CXX) $(CXXFLAGS) $(INCLUDES) -I. $(DEFINEFLAGS) -c $< -o $@
##			@echo "$@ done"



clean:
		@rm -f *.o wrkdir/*
		cd rootstuff; $(MAKE) clean
