!***********************************************************************
! Diffeq_Module
!**begin prologue     Diffeq_Module
!**date written       090219   (yymmdd)
!**revision date               (yymmdd)
!**keywords           DVR, FEDVR
!**
!**author             schneider, b. i.(nsf)
!**source             DVR Library
!**purpose            Calculates two electron integrals in DVR basis
!***                  
!***description       Compute  < rho_ik | 1/r_12 | rho_jl > where
!***                        rho_ij(r) = psi_i(r) * psi_j(r) and 
!***                        psi_i(r_j) = delta_ij/sqrt(w_i)
!***                        The points and weights are those defined
!***                        in the FEDVR quadrature basis. 
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***                  
!***                  
!***end prologue      Diffeq_Module
!***********************************************************************
!***********************************************************************
                           MODULE Diffeq_Module
                            USE Matrix_Scale_and_Assemble
                            USE Matrix_Diagonalization
                            USE Poisson_Module
                            USE FEDVR_Shared
                            USE FEDVR_Derived_Types
                            USE Data_Module
                            USE Legendre_Data
!
                              IMPLICIT NONE
  INTEGER                                          :: maximum_orbital_l
  INTEGER                                          :: maximum_orbital_m
  INTEGER                                          :: minimum_orbital_m
  INTEGER                                          :: maximum_total_L
  INTEGER                                          :: maximum_total_M
  INTEGER                                          :: minimum_total_M
!***********************************************************************
!                           Explicit Interfaces
!***********************************************************************
!
!
!
!
!
  INTEGER                                          :: nen
  INTEGER                                          :: ndim
  INTEGER                                          :: lm
  INTEGER                                          :: begin
  INTEGER                                          :: end
  INTEGER, DIMENSION(10)                           :: lchr
  CHARACTER(LEN=16)                                :: type_eq
  CHARACTER(LEN=16)                                :: type_operator
  CHARACTER(LEN=3)                                 :: new_old
  REAL(idp), DIMENSION(:), ALLOCATABLE             :: energy
  REAL(idp)                                        :: R_factor  
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck Read_Diffeq
!***begin prologue     Read_Diffeq
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            data input for spherical or spheroidal radial coordinate
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       Read_Diffeq
  SUBROUTINE Read_Diffeq
  IMPLICIT NONE
  REAL(idp)                                        :: fpkey
  LOGICAL                                          :: dollar
  LOGICAL                                          :: logkey
  INTEGER                                          :: intkey  
  INTEGER                                          :: i

  CHARACTER(LEN=80)                                :: chrkey
  CHARACTER(LEN=8)                                 :: itoc
  TYPE(spherical)                                  :: name_spherical
  TYPE(spheroidal)                                 :: name_spheroidal
!
  IF ( dollar('$pde',card,cpass,inp) ) THEN
       keyword = chrkey(card,'coordinate_system','spherical',' ')
       type_eq = chrkey(card,'type_differential_equation','eigenvalue',' ')
       form_hamiltonian=logkey(card,'form_hamiltonian',.false.,' ')
       IF ( form_hamiltonian == .true. ) THEN
            read_hamiltonian=.true.
            diagonalize_hamiltonian=logkey(card,'diagonalize_hamiltonian',.false.,' ')
            print_hamiltonian=logkey(card,'print_hamiltonian',.false.,' ')
       END IF
       form_nabla=logkey(card,'form_nabla',.false.,' ')
       IF ( form_nabla == .true. ) THEN
            read_nabla=.true.
            diagonalize_nabla=logkey(card,'diagonalize_nabla',.false.,' ')
            print_nabla=logkey(card,'print_nabla',.false.,' ')
       END IF
       nen = intkey(card,'number_of_energies',0,' ')
       ndim = intkey(card,'number_of_dimensions',1,' ')
       ALLOCATE(reg_grid(1:ndim))
       lchr(1)=lenth(keyword)
       write(iout,1) keyword(1:lchr(1)), type_eq, ndim, nen, form_hamiltonian, diagonalize_hamiltonian 
       IF ( nen > int_zero) THEN
            ALLOCATE (energy(1:nen))
            Call fparr( card,'energies',energy,nen,' ')
       END IF
       maximum_orbital_l = intkey(card,'maximum_orbital_l',0,' ')
       maximum_orbital_m = intkey(card,'maximum_orbital_m',maximum_orbital_l,' ')
       minimum_orbital_m = - maximum_orbital_l 
       maximum_total_L = intkey(card,'maximum_total_L',maximum_orbital_l + maximum_orbital_l,' ') 
       maximum_total_M = intkey(card,'maximum_total_M',maximum_total_L,' ')  
       minimum_total_M = - maximum_total_M
       write(iout,2) maximum_orbital_l, maximum_orbital_m, maximum_total_L, maximum_total_M
       DO i = 1, 12
          prn(i) = logkey(card,prnkey(i),.false.,' ')
       END DO
       DO i=1,ndim
          reg_grid(i)%label = chrkey(card,'coordinate_label_'//itoc(i),'r',' ')
       END DO
  END IF
  IF (keyword(1:lchr(1)) /= 'spherical'.and.keyword(1:lchr(1))/='spheroidal') THEN
      Call lnkerr(' Quit.  Coordinate System not spherical or spheroidal ')
  END IF
  new_old='old'  
1 FORMAT(/,25x,'Coordinate System  = ',a16,5x,'Type Equation           = ',a16,      &
         /25x,'number of dimension = ',i2, 5x,'number of energies      = ',i3,       &
         /25x,'form_Hamiltonion    = ',l1, 5x,'diagonalize_Hamiltonian = ',l1)
2 Format(/,15x,'Maximum Orbital l Read In = ',i4,/,15x,                              &
               'Maximum Orbital m Read In = ',i4,/,15x,                              &
               'Maximum Total L Read In   = ',i4,/,15x,                              &
               'Maximum Total M Read In   = ',i4 )
  END SUBROUTINE Read_Diffeq
!***********************************************************************
!***********************************************************************
!deck Read_FEDVR_Information
!***begin prologue     Read_FEDVR_Information
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            data input for spherical or spheroidal radial coordinate
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       Spherical_Data
!  SUBROUTINE Read_FEDVR_Information(grid)
  SUBROUTINE Read_FEDVR_Information(grid)
  IMPLICIT NONE
  TYPE(coordinates)                             :: grid
  TYPE(Real_Vector)                             :: type_real_vector
  CHARACTER(LEN=3)                              :: itoc
  INTEGER                                       :: i
  INTEGER                                       :: lm
!
  Write(iout,1) FEDVR_File
  Call IOsys('read integer "number of physical grid points" from '//FEDVR_File,1,physical_points,0,0)
  Call IOsys('read integer "drop left" from '//FEDVR_File,1,drop(1),0,0)
  Call IOsys('read integer "drop right" from '//FEDVR_File,1,drop(2),0,0)
  ALLOCATE(grid%grid_points(1:physical_points), grid%grid_weights(1:physical_points))
  Call IOsys('read real "physical grid points" from '//FEDVR_File,physical_points,grid%grid_points,0,' ')
  Call IOsys('read real "physical grid weights" from '//FEDVR_File,physical_points,grid%grid_weights,0,' ')
  Call Print_Matrix(type_real_vector,grid%grid_points,title=grid%label(1:end)//' Points')
  Call Print_Matrix(type_real_vector,grid%grid_weights,title=grid%label(1:end)//' Weights')
  Call IOsys('read real "last grid point" from '//FEDVR_File,1,R_max,0,' ')
  IF (drop(2) == .true. ) THEN
      Write(iout,2) R_max
  ELSE
      Write(iout,3) R_max
  END IF
!
!           This will read in the global Hamiltonian or Nabla matrices.
!
  Call IOsys('read integer "L Max" from '//FEDVR_File,1,l_max,0,' ')
  Call IOsys('read integer "M Max" from '//FEDVR_File,1,m_max,0,' ')
  lm_max=max(l_max,m_max)
  ALLOCATE(dvr_mat(0:lm_max))
  Write(iout,4)
  IF (read_nabla == .true.) THEN
      DO lm = 0, lm_max
         ALLOCATE(dvr_mat(lm)%tr(1:physical_points,1:physical_points))
      END DO
  END IF
  IF (read_hamiltonian == .true.) THEN
      DO lm = 0, lm_max
         ALLOCATE(dvr_mat(lm)%ham(1:physical_points,1:physical_points))
      END DO
  END IF
  Call Read_Global_Matrices(grid)
  Write(iout,6) FEDVR_File
  Call IOsys('close '//FEDVR_File,0,0,0,' ')
1 Format(/15x,'Opening ',a24,' and reading in all the grid information')
2 Format(/15x,'The last grid point = 'f12.6,' has been dropped due to boundary conditions')
3 Format(/15x,'The last grid point = 'f12.6,' has not been dropped due to boundary conditions')
4 Format(/15x,'Reading in all the matrices')
5 Format(/15x,'Operator Selected = ',a24)
6 Format(/15x,'Closing ',a24)
  END SUBROUTINE Read_FEDVR_Information
!
!***********************************************************************
!***********************************************************************
! Read_Global_Matrices
!***begin prologue     Read Global_Matrices
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Read_Global_Matrices
!***references
!***routines called    iosys, util and mdutil
!***end prologue       
!
  SUBROUTINE Read_Global_Matrices(grid)
  IMPLICIT NONE
  TYPE (coordinates)                      :: grid
  TYPE(Real_Matrix)                       :: type_real_matrix
  INTEGER                                 :: lm
  CHARACTER(LEN=3)                        :: itoc
!  
  write(iout,*)
  IF (read_hamiltonian) THEN
      DO lm = 0, lm_max
         title='global H-mat l/m = '//itoc(lm)
         len=lenth(title)
         write(iout,*)
         write(iout,*) '               Reading '//title(1:len)//' from disk'
         Call IOsys('read real "'//title//' " from '//FEDVR_File,                         &
                     physical_points*physical_points,dvr_mat(lm)%ham,0,' ')
         IF (print_hamiltonian) THEN
             Call Print_Matrix(type_real_matrix,dvr_mat(lm)%ham,physical_points,          &
                               physical_points,title=title)
         END IF
      END DO
  END IF
  IF (read_nabla) THEN
      DO lm = 0, lm_max
         title = 'global T-mat l/m = '//itoc(lm)
         len=lenth(title)
         write(iout,*)
         write(iout,*) '               Read '//title(1:len)//' from disk'
         Call IOsys('read real "'//title//' " from '//FEDVR_File,                         &
                     physical_points*physical_points,dvr_mat(lm)%tr,0,' ')
         IF (print_nabla) THEN
             Call Print_Matrix(type_real_matrix,dvr_mat(lm)%tr,physical_points,           &
                               physical_points,title=title)
         END IF
      END DO
  END IF
END SUBROUTINE Read_Global_Matrices
!***********************************************************************
!***********************************************************************
!deck PDE
!***begin prologue     PDE
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Solve the appropriate differential equation in FEDVR basis
!***description        
!***                   
!***                   
!***references
!***routines called
!***end prologue       PDE
  SUBROUTINE PDE (grid)
  IMPLICIT NONE
  TYPE(coordinates)                        :: grid
  TYPE(spherical)                          :: name_spherical
!
!            It should be noted here again the operators that come in are SYMMETRIZED.  That means that in
!            their construction there were factors that were required to turn what would be a generalized
!            eigenvalue problem into a simple eigenvalue problem.
!              
!
  IF (type_eq == 'eigenvalue') THEN
      Call Diagonalize_Global_Matrix(grid)
      Write(iout,1)
  ELSE IF (type_eq == 'linear_system') THEN
      Call Poisson_Equation(grid)
  Write(iout,2)
  ELSE
      Call lnkerr('error in type equation.  Quit')
  END IF
!
1 Format(/,25x,'Finished Diagonalizing the Global Matrix')
2 Format(/,25x,'Finished Solving the Linear System for Global Matrix')
END SUBROUTINE PDE
!***********************************************************************
!***********************************************************************
           END MODULE Diffeq_Module
!***********************************************************************
!***********************************************************************
