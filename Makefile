#  @(#)Makefile	5.1 11/6/94
#
# This assumes the following directory system as well as makefiles
#
# You should set your environment variables for the fortran compiler
# and the various flags and libraries needed.
#
TOP_DIR = $(FEDVR_LIB)
include $(TOP_DIR)/Makefile.inc
#
all: $(FEDVR)
.RECURSIVE: $(FEDVR)
$(FEDVR): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

make: FORCE
	cd $(FEDVR_MODULES) ; $(MAKE) $(MFLAGS)
	cd $(MATRIX_ELEMENTS) ; $(MAKE) $(MFLAGS)
	cd $(DIFFEQ) ; $(MAKE) $(MFLAGS)
	cd $(TWO_ELECTRON_INTEGRALS) ; $(MAKE) $(MFLAGS)

clean: FORCE
	cd $(FEDVR_MODULES) ; $(MAKE) $(MFLAGS) clean
	cd $(MATRIX_ELEMENTS) ; $(MAKE) $(MFLAGS) clean
	cd $(DIFFEQ) ; $(MAKE) $(MFLAGS) clean
	cd $(TWO_ELECTRON_INTEGRALS) ; $(MAKE) $(MFLAGS) clean
	rm -f *.a *~

clean_rep: FORCE
	rm RCS* -fR svn* .svn*
	cd $(FEDVR_MODULES) ; rm RCS* -fR svn* .svn* ; $(MAKE) $(MFLAGS) clean_rep
	cd $(MATRIX_ELEMENTS) ; rm RCS* -fR svn* .svn* ; $(MAKE) $(MFLAGS) clean_rep
	cd $(DIFFEQ) ; rm RCS* -fR svn* .svn* ; $(MAKE) $(MFLAGS) clean_rep
	cd $(TWO_ELECTRON_INTEGRALS) ; rm RCS* -fR svn* .svn* ; $(MAKE) $(MFLAGS) clean_rep


FORCE:
