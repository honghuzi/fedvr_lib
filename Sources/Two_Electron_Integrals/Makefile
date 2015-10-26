# @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f .f90 .mod $(SUFFIXES)
TOP_DIR = $(SOURCES)/fedvr_lib
include $(TOP_DIR)/Makefile.inc
SRCS = \
	Two_Electron_Integrals.f90 

all: $(TWO_ELECTRON_DRIVER) 

$(TWO_ELECTRON_DRIVER) : $(SRCS:.f90=.o) 
	$(LOAD) $(TWO_ELECTRON_DRIVER)

files: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod *~

clean_rep:
	rm -fR RCS* svn* .svn*
