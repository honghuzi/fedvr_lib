!***********************************************************************
! DVR_Module
!**begin prologue     DVR_Module
!**date written       090219   (yymmdd)
!**revision date               (yymmdd)
!**keywords           DVR, FEDVR
!**
!**author             schneider, b. i.(nsf)
!**source             DVR Library
!**purpose            This is the driver module to calculate the FEDVR functions 
!***                  and matrix elements.  The actual work is done in the 
!***                  subroutines that are called.
!***                  in a FEDVR basis
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***                  
!***                  
!***end prologue      DVR_Module
!***********************************************************************
!***********************************************************************
                           MODULE DVR_Module
                           USE Data_Module
                           USE FEDVR_Shared
                           USE FEDVR_Derived_Types
                           USE DVR_Polynomials_Module
                           USE DVR_Kinetic_Energy_Module
                           USE DVR_H_0_Module
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck KE_DVR_Matrices
!***begin prologue     KE_DVR_Matrices
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Driving routine to compute the raw sector matrix elements
!***                   of the kinetic energy operator.  These are done for each sector
!***                   with no regard for boundary conditions or normalization.
!***                   For certain coordinates even and odd m must be treated differently.
!***                   The real work is done in the DVR_Kinetic_Energy_Module where the
!***                   quadrature is performed by Kinetic_Energy.
!***  
!***references
!***routines called    iosys, util and mdutil
!***end prologue       
!
  SUBROUTINE KE_DVR_Matrices(grid)
  IMPLICIT NONE
  TYPE (coordinates)             :: grid
!
  len=lenth(grid%label)
!
  IF (keyword == 'cartesian') THEN
!
      IF(typwt == 'one' .or. typwt == 'legendre') THEN      
         ALLOCATE(grid%reg_mat(1:nreg))
         Call Kinetic_Energy(grid, grid%reg_mat, grid%reg_poly)
      ELSE IF (typwt == 'hermite') THEN
         ALLOCATE(grid%reg_mat_hermite(1:nreg))
         Call Kinetic_Energy(grid, grid%reg_mat_hermite, grid%reg_poly)
      END IF
!
  ELSE IF(keyword =='spherical') THEN
!
      IF(grid%label(1:len) == 'r') THEN
         IF(typwt == 'one' .or. typwt == 'legendre') THEN
            ALLOCATE(grid%reg_mat(1:nreg))
            Call Kinetic_Energy(grid, grid%reg_mat, grid%reg_poly)
         ELSE IF (typwt == 'laguerre') THEN
            ALLOCATE(grid%reg_mat_laguerre(1))
            Call Kinetic_Energy(grid, grid%reg_mat_laguerre, grid%reg_poly)
         ELSE IF (typwt == 'spherical_hermite') THEN        
            ALLOCATE(grid%reg_mat_hermite(1:nreg))
            Call Kinetic_Energy(grid, grid%reg_mat_hermite, grid%reg_poly)
         END IF
!
      ELSE IF( grid%label(1:len) == 'theta') THEN
!
!        Compute the Even KE.
!
         pre_factor = - one
         ALLOCATE(grid%reg_mat(1:nreg))
         Call Kinetic_Energy(grid, grid%reg_mat, grid%reg_poly)
!
!            Compute the Odd KE.
!
         IF (m_max > 0 ) THEN
             ALLOCATE(grid%reg_mat_odd(1:nreg))
             Call Kinetic_Energy(grid, grid%reg_mat_odd, grid%reg_poly)
         END IF
      END IF
!
  ELSE IF(keyword =='cylindrical') THEN
!
      IF(grid%label(1:len) == 'rho') THEN
         ALLOCATE(grid%reg_mat(1:nreg))
         Call Kinetic_Energy(grid, grid%reg_mat, grid%reg_poly)
      END IF
!
  ELSE IF (keyword == 'spheroidal') THEN
!
      ALLOCATE(grid%reg_mat(1:nreg))
      IF (grid%label(1:len) == 'eta') THEN
          pre_factor = - one
      ELSE IF(grid%label(1:len) == 'xi') THEN
          pre_factor = one
      END IF
!
!     Compute the even kinetic energy matrix elements.
!
      Call Kinetic_Energy(grid, grid%reg_mat, grid%reg_poly)
      IF ( m_max > 0 ) THEN
           ALLOCATE(grid%reg_mat_odd(1:nreg))
!
!      Compute the odd kinetic energy matrix elements.
!
           Call Kinetic_Energy(grid, grid%reg_mat_odd, grid%reg_poly)
      END IF
  ELSE IF (keyword == 'fourier') THEN
!
      ALLOCATE( grid%reg_mat_fourier(1))
      Call Kinetic_Energy(grid, grid%reg_mat_fourier)
!
!
  END IF
!
!
END SUBROUTINE KE_DVR_Matrices
!***********************************************************************
!***********************************************************************
!deck Final_Functions
!***begin prologue     Final_Functions
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Take the basic interpolation polynomials in each sector
!***                   and renormalize them.  The renormalization factors are computed
!***                   so that the polyomials at beginning and end of the sector, the
!***                   bridge functions, are properly normalized. See the normaliztion
!***                   routine for details. 
!***references
!***routines called    iosys, util and mdutil
!***end prologue       
!
  SUBROUTINE Final_Functions(grid)
  IMPLICIT NONE
  TYPE(REAL_MATRIX)              :: type_real_matrix
  TYPE (coordinates)             :: grid
  INTEGER                        :: ir
  INTEGER                        :: i
!
!
!                  Normalize the functions.
!
  DO ir = 1, nreg 
     write(iout,1) ir
     DO i = 1, npt(ir)     
        grid%reg_poly(ir)%pr(:,i)                                                    &
                                  =                                                  &
        grid%reg_poly(ir)%pr(:,i) * grid%reg_poly(ir)%normalization(i)
        grid%reg_poly(ir)%dpr(:,i)                                                   &
                                  =                                                  &
        grid%reg_poly(ir)%dpr(:,i) * grid%reg_poly(ir)%normalization(i)
        grid%reg_poly(ir)%ddpr(:,i)                                                  &
                                  =                                                  &
        grid%reg_poly(ir)%ddpr(:,i) * grid%reg_poly(ir)%normalization(i)
     END DO
     IF (prn(3) == .true. ) THEN

         Call Print_Matrix(type_real_matrix,grid%reg_poly(ir)%pr,npt(ir),npt(ir),    &
                           title = 'normalized polynomials')
         Call Print_Matrix(type_real_matrix,grid%reg_poly(ir)%dpr,npt(ir),npt(ir),   &
                           title = 'first derivative of normalized polynomials')

         Call Print_Matrix(type_real_matrix,grid%reg_poly(ir)%ddpr,npt(ir),npt(ir),  &   
                           title = 'second derivative of normalized polynomials')
     END IF
  END DO
!
1 Format(/,10x,'Region = ',i4)
END SUBROUTINE Final_Functions
!***********************************************************************
!***********************************************************************
!deck Final_KE_DVR_Matrices
!***begin prologue     Final_KE_DVR_Matrices
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Compute the renormalized sector matrix elements
!***                   of the kinetic energy operator.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Final_KE_DVR_Matrices
!
  SUBROUTINE Final_KE_DVR_Matrices(grid)
  IMPLICIT NONE
  TYPE (coordinates)             :: grid
!
!
  len=lenth(grid%label)
!
  IF (keyword == 'cartesian') THEN
!
      IF(typwt == 'one' .or. typwt == 'legendre') THEN      
         Call Renormalization(grid,grid%reg_mat)
      ELSE IF (typwt == 'hermite') THEN
         Call Renormalization(grid,grid%reg_mat_hermite)
      END IF
!
  ELSE IF(keyword =='spherical') THEN
!
      IF(grid%label(1:len) == 'r') THEN
!
         IF(typwt == 'one' .or. typwt == 'legendre') THEN
            Call Renormalization(grid,grid%reg_mat)
         ELSE IF (typwt == 'laguerre') THEN
            Call Renormalization(grid,grid%reg_mat_laguerre)
         ELSE IF (typwt == 'hermite') THEN
            Call Renormalization(grid,grid%reg_mat_hermite)
         ELSE IF (typwt == 'spherical_hermite') THEN
            Call Renormalization(grid,grid%reg_mat_hermite)
         END IF
!
      ELSE IF( grid%label(1:len) == 'theta') THEN
!
!        Compute the Even KE.
!
            Call Renormalization(grid,grid%reg_mat)
!
!        Compute the Odd polynomials.
!
         IF (m_max > 0 ) THEN
             Call Renormalization(grid,grid%reg_mat_odd)
         END IF
      END IF
!
  ELSE IF(keyword =='cylindrical') THEN
!
      IF(grid%label(1:len) == 'rho') THEN
         Call Renormalization(grid,grid%reg_mat)
      END IF
!
  ELSE IF (keyword == 'spheroidal') THEN
!
!
!     Compute the even kinetic energy matrix elements.
!
         Call Renormalization(grid,grid%reg_mat)
      IF ( m_max > 0 ) THEN
!
!      Compute the odd kinetic energy matrix elements.
!
           Call Renormalization(grid, grid%reg_mat_odd)
      END IF
  ELSE IF (keyword == 'fourier') THEN
           Call Renormalization(grid, grid%reg_mat_fourier)
  END IF
!
END SUBROUTINE Final_KE_DVR_Matrices
!***********************************************************************
!***********************************************************************
!deck H_0.f
!***begin prologue     H_0
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            The renormalized sector matrix elements of the kinetic energy 
!***                   operator with proper boundary conditions are combined with
!***                   any one body potential.  Then the angular momentum parts are
!***                   added and the matrices "trimmed" if there are boundary conditions
!***                   at the end points which need to be satisfied.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       
!
  SUBROUTINE H_0(grid)
  IMPLICIT NONE
  TYPE (coordinates)                   :: grid
  TYPE (even)                          :: type_even
  TYPE (odd)                           :: type_odd
  TYPE (nabla)                         :: type_nabla
  TYPE (ham)                           :: type_hamiltonian
  INTEGER                              :: i
  INTEGER                              :: j
  INTEGER                              :: lm
  CHARACTER (LEN=3)                    :: itoc
!****************************************************************************************************
!
!          The first step is to take the sector matrix elements which have been computed
!          without regard to angular momenta or potential energy terms and to store them
!          as specific KE matrices.  The subroutines K_0_Matrix are respeonsible for that.  
!          These matrices are then trimmed by removing any leading or trailing basis functions
!          to reflect the specific boundary conditions. The trimmed basic matrices are then 
!          processed by adding appropriate terms for the potential and angular momenta and stored 
!          in the final arrays. 
!
!          Store the original size of the first and last element in these two variables.  They will be
!          reset in Fix_BC if that routine is called and then reset later.
!
!****************************************************************************************************
!                   Begin with seeting up some variables used in printing and labelling files
!****************************************************************************************************
  ALLOCATE(sub_title(1:10),title_len(1:10))
!
!                    Set up some variable to use for printing and writing files
!
  sub_title(1) = 'Final '//grid%label(1:len)//' Hamiltonian Region = '
  title_len(1) = lenth(sub_title(1)) 
  sub_title(2) = 'Final '//grid%label(1:len)//' Nabla Region = '
  title_len(2) = lenth(sub_title(2)) 
  sub_title(3) = '"basic '//grid%label(1:len)//' even hamiltonian region = '
  title_len(3) = lenth(sub_title(3)) 
  Call IOsys('write character file_title_1 to '//FEDVR_File,title_len(3),sub_title(3),0,sub_title(3))
  sub_title(4) = '"basic '//grid%label(1:len)//' odd hamiltonian region = '
  title_len(4) = lenth(sub_title(4)) 
  Call IOsys('write character file_title_2 to '//FEDVR_File,title_len(4),sub_title(4),0,sub_title(4))
  IF (form_nabla) THEN
     sub_title(5) = '"basic '//grid%label(1:len)//' even nabla region = '
     title_len(5) = lenth(sub_title(5)) 
     Call IOsys('write character file_title_3 to '//FEDVR_File,title_len(5),sub_title(5),0,sub_title(5))
     sub_title(6) = '"basic '//grid%label(1:len)//' odd nabla region = '
     title_len(6) = lenth(sub_title(6)) 
     Call IOsys('write character file_title_4 to '//FEDVR_File,title_len(6),sub_title(6),0,sub_title(6))
  END IF
!******************************************************************************************************n
!
! This section of code computes two basic matrices which are then used to put together all the rest
! The basic matrices are the pieces of the kinetic energy that depend on what coordinate system is
! being used but have no added angular momentum nor any added potential matrix elements.
! That will come later.  Before we construct those matrices we need to impose boundary conditions which
! mean trimming the regional first and last matrices.  This is done after this section is concluded.
! K_0_Matrix is the generic routine(s) used to compute the basic matrices.  The routines are located
! in the module DVR_H_0_module.
!******************************************************************************************************
  lm_max=max(l_max,m_max)
  ALLOCATE(grid%reg_type_op(1:nreg,0:lm_max))
!
!       Here we allocate the basic matrices we need.  The actual Hamiltonian and Nabla matrices for
!       each angular momentum and with needed potentials are allocated and computed later.
!
  key_to_allocate = int_zero  ! This variable will be set to one if there is a need to allocate and
                              ! compute the odd matrices.  This will only happen when computing the 
                              ! angular variable theta or the spheroidal variable xi and eta and the only
                              ! if m_max is greater than zero
!
  write(iout,1)
  DO i = 1, nreg
    ALLOCATE(grid%reg_type_op(i,0)%tr(1:npt(i),1:npt(i))) ! Allocate the even tr array it is always needed
  END DO
  IF (keyword == 'cartesian') THEN                        ! Compute grid%reg_type_op(i,0)%tr 
      Write(iout,*) 'Computing Cartesian Coordinate K_0 Matrix = '//grid%label(1:len)
      IF(typwt == 'one' .or. typwt == 'legendre') THEN      
         Call K_0_Matrix(grid, grid%reg_mat)  
      ELSE IF(typwt == 'hermite') THEN
         Call K_0_Matrix(grid, grid%reg_mat_hermite)
      END IF
      Write(iout,*) 'Finished Cartesian Coordinate K_0 Matrix = '//grid%label(1:len)
  ELSE IF(keyword =='spherical') THEN
      IF(grid%label(1:len) == 'r') THEN                   ! Compute grid%reg_type_op(i,0)%tr
         Write(iout,*) 'Computing Spherical Coordinate K_0 Matrix = '//grid%label(1:len)
         IF(typwt == 'one' .or. typwt == 'legendre') THEN
            Call K_0_Matrix(grid, grid%reg_mat)
         ELSE IF (typwt == 'laguerre') THEN
            Call K_0_Matrix(grid, grid%reg_mat_laguerre)
         ELSE IF (typwt == 'hermite') THEN
            Call K_0_Matrix(grid, grid%reg_mat_hermite)
         ELSE IF (typwt == 'spherical_hermite') THEN
            Call K_0_Matrix(grid, grid%reg_mat_hermite)
         END IF
         Write(iout,*) 'Finished Spherical Coordinate K_0 Matrix = '//grid%label(1:len)
      ELSE IF( grid%label(1:len) == 'theta') THEN 
         Write(iout,*) 'Computing Spherical Coordinate K_0 Matrix = '//grid%label(1:len)//' even'
         Call K_0_Matrix(grid, grid%reg_mat)              ! Compute grid%reg_type_op(i,0)%tr 
         Write(iout,*) 'Finished Spherical Coordinate K_0 Matrix = '//grid%label(1:len)//' even'
         IF (m_max > int_zero ) THEN  
             write(iout,2)
             Write(iout,*) 'Computing Spherical Coordinate K_0 Matrix = '//grid%label(1:len)//' odd'
             key_to_allocate = int_one                    ! Set the variable to one so other subroutines know that
             DO i = 1, nreg
                ALLOCATE(grid%reg_type_op(i,1)%tr(1:npt(i),1:npt(i)))  ! Allocate the odd tr array it is needed.
             END DO
             Call K_0_Matrix(grid, grid%reg_mat_odd)      ! Compute grid%reg_type_op(i,1)%tr  
             Write(iout,*) 'Finished Spherical Coordinate K_0 Matrix = '//grid%label(1:len)//' odd'
         END IF
      END IF
  ELSE IF(keyword =='cylindrical') THEN
      IF(grid%label(1:len) == 'rho') THEN
         Write(iout,*) 'Computing Cylindrical Coordinate K_0 Matrix = '//grid%label(1:len)
         Call K_0_Matrix(grid, grid%reg_mat)              ! Compute grid%reg_type_op(i,0)%tr 
         Write(iout,*) 'Finished Cylindrical Coordinate K_0 Matrix = '//grid%label(1:len)
      END IF
  ELSE IF (keyword == 'spheroidal') THEN
      IF (grid%label(1:len) == 'eta') THEN
          pre_factor = - one
      ELSE IF(grid%label(1:len) == 'xi') THEN
          pre_factor = one
      END IF
      Write(iout,*) 'Computing Spheroidal Coordinate K_0 Matrix = '//grid%label(1:len)//' even'
      Call K_0_Matrix(grid, grid%reg_mat)                 ! Compute grid%reg_type_op(i,0)%tr 
      write(iout,*) 'Finished spheroidal K_0 Matrix = '//grid%label(1:len)//' kind = '//' even'
      IF (m_max > int_zero )  THEN   
          key_to_allocate = int_one                       ! Set the variable to one so other subroutine know that
          write(iout,2)
          Write(iout,*) 'Computing Spheroidal Coordinate K_0 Matrix '//grid%label(1:len)//' odd'
          DO i = 1, nreg
             ALLOCATE(grid%reg_type_op(i,1)%tr(1:npt(i),1:npt(i)))  ! Allocate the odd tr array it is needed.
          END DO
          Call K_0_Matrix(grid, grid%reg_mat_odd)        ! Compute grid%reg_type_op(i,1)%tr 
          Write(iout,*) 'Finished Spheroidal Coordinate K_0 Matrix '//grid%label(1:len)//' odd'
      END IF
  ELSE IF (keyword == 'fourier') THEN                     ! Compute grid%reg_type_op(i,0)%tr 
      Call K_0_Matrix(grid, grid%reg_mat_fourier)
  END IF
!******************************************************************************************************
!
!      Now that we have computed the basic tr matrices we can trim them if required
!
  IF (drop(1) == .true. .or. drop(2) == .true. ) THEN  
!
      Call Fix_BC(grid)     ! Trim basic matrices if necesssary and re-store them with their proper dimensions
                            ! Note that the size of the arrays in each region are now set at their proper
                            ! trimmed values.
  END IF
!
!******************************************************************************************************
!
!
! Add the potential and angular momenta to form the Hamiltonian and Nabla matrices.
!
!
  Call Form_Final_Sector_Matrices(grid, grid%label(1:len))
!
!
! Before we construct the Hamiltonian and Nabla matrices, we need to write some information to the disk.
!
  Call IOsys('write integer "number of regions" to '//FEDVR_File,1,nreg,0,' ')
  Call IOsys('write integer "number of points in each region" to '//FEDVR_File,nreg,              &
              npt,0,' ')
  DO i  = 1, nreg
     Call IOsys('write real "points region = '//itoc(i)//' " to '//FEDVR_File,                    &
                 npt(i),grid%reg_pt_wt(i)%qr,0,' ')
     Call IOsys('write real "weights region = '//itoc(i)//' " to '//FEDVR_File,                   &
                 npt(i),grid%reg_pt_wt(i)%wtr,0,' ')
     Call IOsys('write real "qr_fac region = '//itoc(i)//' " to '//FEDVR_File,                    &
                 npt(i),grid%reg_pt_wt(i)%qr_fac,0,' ')
     Call IOsys('write real "inv_qr_fac region = '//itoc(i)//' " to '//FEDVR_File,                &
                 npt(i),grid%reg_pt_wt(i)%inv_qr_fac,0,' ')
     Call IOsys('write real "inv_sqrt_qr_fac region = '//itoc(i)//' " to '//FEDVR_File,           &
                 npt(i),grid%reg_pt_wt(i)%inv_sqrt_qr_fac,0,' ')
  END DO
!
  DEALLOCATE(sub_title,title_len)
!
  DO i = 1, nreg
     DEALLOCATE(reg_pot(i)%vec)
  END DO
  DEALLOCATE(reg_pot)
!
1 Format(/,30x,'Allocating the even tr sector matrices.  These are always required',/)
2 Format(/,30x,'Allocating the odd tr sector matrices.   These are needed for this coordinate',/)
END SUBROUTINE H_0
!***********************************************************************
!***********************************************************************
           END MODULE DVR_Module
!***********************************************************************
!***********************************************************************
