#@(#)Makefile	1.4   8/2/91
#
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
include $(FEDVR_LIB)/Makefile.inc
SRCS = FEDVR_Matrix_Elements.f90
all: $(MATRIX_ELEMENT_DRIVER) 

$(MATRIX_ELEMENT_DRIVER): $(SRCS:.f90=.o)
	$(LOAD) $(MATRIX_ELEMENT_DRIVER)

files: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod *~

clean_rep:
	rm -fR RCS* svn* .svn*
