#@(#)Makefile	1.4   8/2/91
#
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
TOP_DIR = $(SOURCES)/fedvr_lib
include $(TOP_DIR)/Makefile.inc
SRCS = Solve_Diffeq.f90
all: $(DIFFEQ_DRIVER) 

$(DIFFEQ_DRIVER): $(SRCS:.f90=.o)
$(DIFFEQ_DRIVER): 
	$(LOAD) $(DIFFEQ_DRIVER)


files: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod *~

clean_rep:
	rm -fR RCS* svn* .svn*
