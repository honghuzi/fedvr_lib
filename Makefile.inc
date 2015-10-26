#
# Including Makefiles must define TOP_DIR appropriately.
#
# $Id: Makefile.inc,v 1.1 2008/11/22 19:43:25 bschneid Exp bschneid $
#
TOP_DIR = $(FEDVR_LIB)
SOURCE_DIR = $(TOP_DIR)/Sources
MODULE_DIR = $(TOP_DIR)/.mod_dir
OBJECT_DIR = $(TOP_DIR)/.o_dir
LIB_DIR = $(TOP_DIR)/Library
RUN_DIR = $(TOP_DIR)/Run
#
FEDVR_MODULE_LIB = $(LIB_DIR)/fedvr_mod_lib.a
LIBS = \
     $(FEDVR_MODULE_LIB) \
     $(MOD_LIB) \
     $(POTENTIAL_LIB) \
     $(UTILITY_LIB) \
     $(MD_LIB)
#
FEDVR_MODULES = $(SOURCE_DIR)/Modules
MATRIX_ELEMENTS = $(SOURCE_DIR)/Matrix_Elements
TWO_ELECTRON_INTEGRALS = $(SOURCE_DIR)/Two_Electron_Integrals
DIFFEQ = $(SOURCE_DIR)/Diffeq
#
FEDVR = \
      $(FEDVR_MODULES) \
      $(MATRIX_ELEMENTS) \
      $(TWO_ELECTRON_INTEGRALS) \
      $(DIFFEQ)
#
# Driver
#
MATRIX_ELEMENT_DRIVER = $(RUN_DIR)/fedvr_me.x
TWO_ELECTRON_DRIVER = $(RUN_DIR)/fedvr_two_el.x
DIFFEQ_DRIVER = $(RUN_DIR)/diffeq.x
# System-dependent definitions
#
#
# Suffix rules
#
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
ALL_SRCS= $(SRCS.F90:.F90=.o) $(SRCS.f90:.f90=.o) $(SRCS:.f=.o) $(SRCS.f:.f=.o) $(SRCS.F:.F=.o)
ALL_INCLUDE = $(INCLUDE_LIST) -I$(MODULE_DIR)
NOLINK = -c
COMPILE = $(F90) $(ALL_INCLUDE) $(NOLINK) $(F90FLAGS)
.F90.o:
	$(COMPILE) $(CPPFLAGS) $<
.f90.o:
	$(COMPILE) $<
.f.o:
	$(COMPILE) $<
.F.o:
	$(COMPILE) $(CPPFLAGS) $<

LOAD = $(F90) $(ALL_INCLUDE) $(ALL_SRCS) $(LIBS) -o
#
# Miscellaneous
#
MAKE = make
MFLAGS =
RCS = RCS
CO = co -u
RANLIB = ranlib
