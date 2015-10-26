#   @(#)Makefile	1.4   8/2/91
#
include $(FEDVR_LIB)/Makefile.inc

SRCS = \
	Utility_Module.f90 \
	Read_DVR_Module.f90 \
	Renormalization_Module.f90 \
	DVR_Polynomials_Module.f90 \
	DVR_Kinetic_Energy_Module.f90 \
	DVR_H_0_Module.f90 \
	DVR_Module.f90 \
	Matrix_Scale_and_Assemble.f90 \
	Matrix_Diagonalization.f90 \
	Poisson_Module.f90 \
	P_LM.f90 \
	Two_Electron_FEDVR_Module.f90 \
	Diffeq_Module.f90 
all: $(FEDVR_MODULE_LIB) 
$(FEDVR_MODULE_LIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(FEDVR_MODULE_LIB) $(?:.f90=.o)
	$(RANLIB) $(FEDVR_MODULE_LIB)
	cp *.mod $(MODULE_DIR)/.
	cp *.o $(OBJECT_DIR)/.
	cd $(MATRIX_ELEMENTS) ; $(MAKE) $(MFLAGS) clean
	cd $(DIFFEQ) ; $(MAKE) $(MFLAGS) clean
	cd $(TWO_ELECTRON_INTEGRALS) ; $(MAKE) $(MFLAGS) clean
clean:
	rm -f *.o 
	rm -f *.mod 
	rm -f $(OBJECT_DIR)/*.o 
	rm -f $(MODULE_DIR)/*.mod
	rm -f $(LIB_DIR)/*.a
	rm -f *~
movein:
	mv *.mod $(MODULE_DIR)/.
	mv *.o $(OBJECT_DIR)/.

moveout:
	mv $(MODULE_DIR)/*.mod .
	mv $(OBJECT_DIR)/*.o .

clean_rep:
	rm -fR RCS* svn* .svn*
