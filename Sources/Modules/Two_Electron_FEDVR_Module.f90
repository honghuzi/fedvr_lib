!***********************************************************************
! Two_Electron_FEDVR_Module
!**begin prologue     Two_Electron_FEDVR_Module
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
!***end prologue      Two_Electron_FEDVR_Module
!***********************************************************************
!***********************************************************************
                           MODULE Two_Electron_FEDVR_Module
                            USE Data_Module
                            USE FEDVR_Shared
                            USE FEDVR_Derived_Types
                            USE Poisson_Module
  TYPE(REAL_MATRIX)                              :: type_real_matrix
  TYPE(REAL_VECTOR)                              :: type_real_vector
  INTEGER, DIMENSION(10)                         :: strlen
  INTEGER                                        :: number_of_radial_points
  INTEGER                                        :: number_of_angular_points
  INTEGER                                        :: number_of_eta_points
  INTEGER                                        :: number_of_xi_points
  INTEGER                                        :: n_tri
  INTEGER                                        :: maximum_orbital_l
  INTEGER                                        :: maximum_orbital_m
  INTEGER                                        :: minimum_orbital_m
  INTEGER                                        :: maximum_total_L
  INTEGER                                        :: maximum_total_M
  INTEGER                                        :: minimum_total_M
  CHARACTER(LEN=8)                               :: Key
!***********************************************************************
!                           Explicit Interfaces
!***********************************************************************
!
!
!
!
                           INTERFACE V_ijkl                        
                       MODULE PROCEDURE V_ijkl_Spherical,          &
                                        V_ijkl_Spheroidal  
                            END INTERFACE V_ijkl
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck Setup_2e
!***begin prologue     Setup_2e
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keyword           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            setup for two electron integral calculation
!***description        
!***references
!***routines called    
!***end prologue       Setup_2e
  SUBROUTINE Setup_2e  
  IMPLICIT NONE
  TYPE(spherical)                                  :: name_spherical
  TYPE(spheroidal)                                 :: name_spheroidal
  TYPE(spherical_harmonics)                        :: harmonics
  REAL(idp)                                        :: fpkey
  LOGICAL                                          :: dollar
  LOGICAL                                          :: logkey
  LOGICAL                                          :: get_angular_files
  INTEGER                                          :: intkey  
  INTEGER                                          :: i
  CHARACTER(LEN=80)                                :: chrkey
!
!
  IF ( dollar('$two_electron_integrals',card,cpass,inp) ) THEN
       keyword = chrkey(card,'coordinate_system','spherical',' ')
       representation = chrkey(card,'representation','spherical_harmonics',' ')
       strlen(1)=lenth(keyword)
       write(iout,1) keyword(1:strlen(1))
       IF (keyword(1:strlen(1)) /= 'spherical'.and.keyword(1:strlen(1))/='spheroidal') THEN
           Call lnkerr(' Quit.  Coordinate System not spherical or spheroidal ')
       END IF
       DO i = 1, 12
          prn(i) = logkey(card,prnkey(i),.false.,' ')
       END DO
  END IF
  maximum_orbital_l = intkey(card,'maximum_orbital_l',0,' ')
  maximum_orbital_m = intkey(card,'maximum_orbital_m',maximum_orbital_l,' ')
  minimum_orbital_m = - maximum_orbital_l 
  maximum_total_L = intkey(card,'maximum_total_L',maximum_orbital_l + maximum_orbital_l,' ') 
  maximum_total_M = intkey(card,'maximum_total_M',maximum_total_L,' ')  
  minimum_total_M = - maximum_total_M
  ALLOCATE(reg_grid(3))
  write(iout,2) maximum_orbital_l, maximum_orbital_m, maximum_total_L, maximum_total_M
!
!   This reads in and manipulates either the spherical or spheroidal radial variable
!   information providing it from the solution to the appropriate differential equation in a
!   FEDVR basis.  From this and the asymptotic form of the irregular solution, it is
!   possible to solve the radial Poisson equation in the same basis.
!
  IF(keyword(1:strlen(1)) == 'spherical') THEN
     reg_grid(1)%label='r'
     FEDVR_File = 'spherical'//'_'//reg_grid(1)%label(1:1)
  ELSE IF (keyword(1:strlen(1)) == 'spheroidal') THEN
      reg_grid(1)%label='xi'
      FEDVR_File = 'spheroidal'//'_'//reg_grid(1)%label(1:2)
  ELSE
      Call lnkerr('quit. Not a spherical or spheroidal system')
  END IF
  strlen(1)=lenth(FEDVR_File)
  write(iout,*) 'Opening Data File = '//FEDVR_File(1:strlen(1))
  file_loc=File_Directory(4)(1:len_dir(4))//'/'//FEDVR_File(1:strlen(1))
  strlen(2)=lenth(file_loc)
  write(iout,*) 'File Location = ',file_loc(1:strlen(2))
  Call IOsys('open '//FEDVR_File//' as old',0,0,0,file_loc(1:strlen(2)))
  Call IOsys('read integer "L Max" from '//FEDVR_File,1,l_max,0,' ')
  Call IOsys('read integer "M Max" from '//FEDVR_File,1,m_max,0,' ')
  Call Read_Input_Data(reg_grid(1),keyword(1:strlen(1)),reg_grid(1)%label)
  Call IOsys('rewind all on '//FEDVR_File//' read-and-write',0,0,0,' ')
  Call IOsys('close '//FEDVR_File,0,0,0,' ')
  Call Poisson_Equation(reg_grid(1))
!
!  Now we do the angular variable either using the spherical harmonic representation
!  which translates into Clebsch-Gordan coefficients or from the FEDVR solution to the
!  angular differential equation.  In the latter approach we directly use the expansion
!  of the angular functions in the FEDVR basis and it then basically simple to get what
!  is needed.
!
  IF(keyword(1:strlen(1)) == 'spherical') THEN
     Call V_ijkl (reg_grid(1), name_spherical)
     write(iout,*) 'Opening Data File = Spherical_2_Electron_Integrals as new'
  ELSE IF (keyword(1:strlen(1)) == 'spheroidal') THEN
     Call V_ijkl (reg_grid(1), name_spheroidal)
     write(iout,*) 'Opening Data File = Spheroidal_2_Electron_Integrals as new'
  END IF
  get_angular_files = .true.
  IF(keyword(1:strlen(1)) == 'spherical') THEN
     reg_grid(2)%label='theta'
     FEDVR_File = 'spherical'//'_'//reg_grid(2)%label(1:5)
     IF (representation == 'spherical_harmonics' ) THEN
         get_angular_files = .false.
     END IF
  ELSE IF (keyword(1:strlen(1)) == 'spheroidal') THEN
      reg_grid(2)%label='eta'
      FEDVR_File = 'spheroidal'//'_'//reg_grid(2)%label(1:3)
  ELSE
      Call lnkerr('quit. Not a spherical or spheroidal system')
  END IF
  IF (get_angular_files == .true.) THEN
      strlen(1)=lenth(FEDVR_File)
      write(iout,*) 'Opening Data File = '//FEDVR_File(1:strlen(1))
      file_loc=File_Directory(7)(1:len_dir(7))//'/'//FEDVR_File(1:strlen(1))
      strlen(2)=lenth(file_loc)
      write(iout,*) 'File Location = ',file_loc(1:strlen(2))
      Call IOsys('open '//FEDVR_File//' as old',0,0,0,file_loc(1:strlen(2)))
      Call IOsys('read integer "L Max" from '//FEDVR_File,1,l_max,0,' ')
      Call IOsys('read integer "M Max" from '//FEDVR_File,1,m_max,0,' ')
      Call Read_Input_Data(reg_grid(2),keyword(1:strlen(1)),reg_grid(2)%label)
      Call IOsys('rewind all on '//FEDVR_File//' read-and-write',0,0,0,' ')
      Call IOsys('close '//FEDVR_File,0,0,0,' ')
  ELSE
      Call D_LM
  END IF
  Call IOsys('open FEDVR_Two_Electron_Integral_File as new',0,0,0,                    &
             'Spherical_2_Electron_Integrals')
1 FORMAT(/,25x,'Coordinate System = ',a32)
2 Format(/,15x,'Maximum Orbital l Read In = ',i4,/,15x,                               &
               'Maximum Orbital m Read In = ',i4,/,15x,                               &
               'Maximum Total L Read In   = ',i4,/,15x,                               &
               'Maximum Total M Read In   = ',i4 )
  END SUBROUTINE Setup_2e
!***********************************************************************
!***********************************************************************
!deck Read_Input_Data
!***begin prologue     Read_Input_Data
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
  SUBROUTINE Read_Input_Data(grid,keyword,label)
  IMPLICIT NONE
  TYPE(coordinates)                             :: grid
  TYPE(spherical_fedvr)                         :: fedvr
  CHARACTER(*)                                  :: keyword
  CHARACTER(*)                                  :: label
  CHARACTER(LEN=3)                              :: itoc
  INTEGER                                       :: i
  INTEGER                                       :: lm
  INTEGER                                       :: l_1
  INTEGER                                       :: l_
  INTEGER                                       :: begin
  INTEGER                                       :: end
!
  Call IOsys('read integer "number of regions" from '//FEDVR_File,1,grid%num_reg,0,' ')
  ALLOCATE(grid%num_pts_reg(1:grid%num_reg),grid%reg_pt_wt(1:grid%num_reg))
  Call IOsys('read integer "number of points in each region" from '//FEDVR_File,grid%num_reg,         &
              grid%num_pts_reg,0,' ')
  nreg=grid%num_reg
  npt(1:nreg) =  grid%num_pts_reg(1:nreg)
  write(iout,*) nreg, npt(1:nreg)
  DO i = 1, grid%num_reg
     ALLOCATE(grid%reg_pt_wt(i)%qr(1:grid%num_pts_reg(i)),                                            &
              grid%reg_pt_wt(i)%wtr(1:grid%num_pts_reg(i)),                                           &
              grid%reg_pt_wt(i)%qr_fac(1:grid%num_pts_reg(i)),                                        &
              grid%reg_pt_wt(i)%inv_qr_fac(1:grid%num_pts_reg(i)),                                    &
              grid%reg_pt_wt(i)%inv_sqrt_qr_fac(1:grid%num_pts_reg(i)))
     Call IOsys('read real "points region = '//itoc(i)//' " from '//FEDVR_File,                       &
                 grid%num_pts_reg(i),grid%reg_pt_wt(i)%qr,0,' ')
     Call IOsys('read real "weights region = '//itoc(i)//' " from '//FEDVR_File,                      &
                 grid%num_pts_reg(i),grid%reg_pt_wt(i)%wtr,0,' ')
     Call IOsys('read real "qr_fac region = '//itoc(i)//' " from '//FEDVR_File,                      &
                 grid%num_pts_reg(i),grid%reg_pt_wt(i)%qr_fac,0,' ')
     Call IOsys('read real "inv_qr_fac region = '//itoc(i)//' " from '//FEDVR_File,                   &
                 grid%num_pts_reg(i),grid%reg_pt_wt(i)%inv_qr_fac,0,' ')
     Call IOsys('read real "inv_sqrt_qr_fac region = '//itoc(i)//' " '//'from '//FEDVR_File,         &
                 grid%num_pts_reg(i),grid%reg_pt_wt(i)%inv_sqrt_qr_fac,0,' ')
  END DO  
  lm_max = max(l_max,m_max)
!
! Get the total number of points, make the total final point and weight arrays
!
  physical_points = grid%num_pts_reg(1)
  DO i = 2, grid%num_reg 
     physical_points = physical_points + grid%num_pts_reg(i) - 1
  END DO
  ALLOCATE(grid%grid_points(1:physical_points), grid%grid_weights(1:physical_points) )
  grid%grid_points(1:grid%num_pts_reg(1)) = grid%reg_pt_wt(1)%qr(1:grid%num_pts_reg(1))
  grid%grid_weights(1:grid%num_pts_reg(1)) = grid%reg_pt_wt(1)%wtr(1:grid%num_pts_reg(1))
  begin = grid%num_pts_reg(1)
  DO i = 2, grid%num_reg 
     begin = begin + 1
     end = begin + grid%num_pts_reg(i) - 2
     grid%grid_points(begin:end) = grid%reg_pt_wt(i)%qr(2:grid%num_pts_reg(i))
     grid%grid_weights(begin:end) = grid%reg_pt_wt(i)%wtr(2:grid%num_pts_reg(i))
     begin = end
  END DO
  n_final = end
  n_tri = n_final * (n_final + int_one ) / int_two
  Call IOsys('read real "last grid point" from '//FEDVR_File,1,R_max,0,' ')
  write(iout,1) l_max, m_max, physical_points, R_max
  Call Print_Matrix(type_real_vector,grid%grid_points,title = label//' Points')
  title = label//' Weights'
  Call Print_Matrix(type_real_vector,grid%grid_weights,title = label//' Weights')
!
!           This will read in the regional Nabla matrices.
!
  ALLOCATE(grid%reg_type_op(1:grid%num_reg,0:lm_max))
  DO lm = 0, lm_max
     DO i = 1,grid%num_reg
        ALLOCATE ( grid%reg_type_op(i,lm)%tr(1:grid%num_pts_reg(i),1:grid%num_pts_reg(i)) )
        Call IOsys('read real "scaled nabla lm = '//itoc(lm)//' region = '//itoc(i)//' " from '//FEDVR_File,       &
                   grid%num_pts_reg(i)*grid%num_pts_reg(i),grid%reg_type_op(i,lm)%tr,0,' ')
     END DO
  END DO
  IF (prn(4) == .true. ) THEN
      DO lm = 0, lm_max
         DO i = 1, grid%num_reg
            Call Print_Matrix(type_real_matrix,grid%reg_type_op(i,lm)%tr,grid%num_pts_reg(i),grid%num_pts_reg(i),  &
                              title = 'Scaled '//label//' Nabla matrix LM = '//itoc(lm)//' Region '//itoc(i) )
         END DO
      END DO
  END IF
!
!           Form the gloabl Nabla matrix
!
  ALLOCATE(dvr_mat(0:l_max))
  Call Form_Global_Nabla(grid,grid%label)
1 Format(/,15x,'maximum orbital l from disk       = ',i4,/,15x,                       &
               'maximum orbital m from disk       = ',i4,/,15x,                       &
               'number of points from disk = ',i5,/,15x,                              &
               'box size from disk                = ',f15.8)
  END SUBROUTINE Read_Input_Data
!***********************************************************************
!***********************************************************************
!deck Form_Global_Nabla.f
!***begin prologue     Form_Global_Nabla
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Form the global matrices.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       
!
  SUBROUTINE Form_Global_Nabla(grid,label)
  IMPLICIT NONE
  TYPE (coordinates)                      :: grid
  CHARACTER(LEN=*)                        :: label
  INTEGER                                 :: i
  INTEGER                                 :: lm
  INTEGER                                 :: first
  INTEGER                                 :: last
  CHARACTER(LEN=3)                        :: itoc
!  
DO lm = 0, lm_max
   ALLOCATE(dvr_mat(lm)%tr(1:physical_points,1:physical_points))
   dvr_mat(lm)%tr(:,:) = zero
   first = 1
   DO i = 1, grid%num_reg
      last = first + grid%num_pts_reg(i) - 1
      dvr_mat(lm)%tr(first:last,first:last) =                       &
                  grid%reg_type_op(i,lm)%tr(1:grid%num_pts_reg(i),1:grid%num_pts_reg(i))         
      first = last
      DEALLOCATE(grid%reg_type_op(i,lm)%tr)
   END DO    
   IF (prn(9) == .true. ) THEN
       Call Print_Matrix(type_real_matrix,dvr_mat(lm)%tr,physical_points,physical_points,           &
                         title = 'Global Nabla L or M = '//itoc(lm))
   END IF
END DO
END SUBROUTINE Form_Global_Nabla
!***********************************************************************
!***********************************************************************
!deck V_ijkl_Spherical
!***begin prologue     V_ijkl_Spherical
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            V_ijkl_Spherical
!***description        The radial parts of the two-electron integrals
!***                   in spherical coordinates are computed for all
!***                   required L values.  They do not explicitly depend on M.
!***references
!***routines called
!***end prologue       V_ijkl_Spherical
  SUBROUTINE V_ijkl_Spherical (grid, name_spherical)
  IMPLICIT NONE
  TYPE(coordinates)                        :: grid  
  TYPE(spherical)                          :: name_spherical
  REAL(idp)                                :: R_factor  
  REAL(idp), DIMENSION(:), ALLOCATABLE     :: factr
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: l_factor  
  INTEGER                                  :: l_two
  INTEGER                                  :: info  
  INTEGER                                  :: lm  
  INTEGER                                  :: count
  INTEGER                                  :: n_tri
  INTEGER                                  :: maximum_total_L
  CHARACTER(LEN=3)                         :: itoc
  REAL(idp)                                :: value
!
!       V_ijkl = delta_ik delta_jl ( r_i )^(L+2)*( r_j )^(L+2) /( r_n)^(2L+1) 
!                                          - 
!                         (2L+1) * ( r_i )^2 * ( T_ij )^-1 /sqrt(w_i*w_j}
!
!
!      Compute the inverse from the LU decomposition.  The last point
!      is excluded since the solution of the Poisson equation we are
!      computing is zero at the last point.  We will add a solution of the
!      homogeneous equation to take care of the boundary conditions later. 
!      
!         
  ALLOCATE( dvr_mat(0)%ipvt( n_final ), dvr_mat(0)%lower( n_tri ),                  &
            dvr_mat(0)%Inverse(n_final,n_final),                                    &
            dvr_mat(0)%tr(physical_points,physical_points), factr( n_total ) )
  factr(1:n_total) = grid%grid_points(1:n_total) / sqrt ( grid%grid_weights(1:n_total) ) 
  keyword = 'spherical'
  strlen(1)=lenth(keyword)
  grid%label = 'r'
  strlen(2)=lenth(grid%label)
!
!    Set up the matrix to be used to solve the linear equations
!
  DO lm = 0, maximum_total_L
     title=keyword(1:strlen(1))//'_'//grid%label(1:strlen(2))//'_'//itoc(lm)
     strlen(3)=lenth(title)
     write(iout,*) 'Reading Global Nabla = ',title(1:strlen(3))//' from disk'
     Call IOsys('read real "'//title(1:strlen(3))//' nabla" from '//FEDVR_File,           &
                 physical_points*physical_points,dvr_mat(0)%tr,0,' ')
     count = 0
     DO i = 1, n_final
        DO j = 1, i
           count = count + 1
           dvr_mat(0)%lower(count) = dvr_mat(lm)%tr(i,j)
        END DO
     END DO
!
!    Factor the matrix
!
     Call DSPTRF ('u',n_final,dvr_mat(0)%lower,dvr_mat(0)%ipvt,info)

!
!    Compute the inverse
!
     dvr_mat(0)%Inverse(:,:) = zero
     DO i = 1, n_final
        dvr_mat(0)%Inverse(i,i) = one
     END DO
!
     Call DSPTRS('u',n_final,n_final,dvr_mat(0)%lower,dvr_mat(0)%ipvt,              &
                                     dvr_mat(0)%Inverse,n_total,info)
!
!         
     DO i = 1, n_final
        dvr_mat(0)%Inverse(1:n_final,i) = factr(1:n_final) * dvr_mat(0)%Inverse(1:n_final,i) * factr(i) 
     END DO
     l_factor = - ( lm + lm + int_one )
     dvr_mat(0)%Inverse(:,:) = l_factor * dvr_mat(0)%Inverse(:,:)
     l_two = lm + int_two     
     factr(1:n_total) = grid%grid_points(1:n_total) ** l_two
!
!    Add the solution to the homogeneous equation
!                                                                           
     ALLOCATE( dvr_mat(lm)%Q(n_total,n_total) )
     write(iout,*) 'Begin Calculation of Basic Radial Two-Electron Integrals L = '//itoc(lm)
     dvr_mat(lm)%Q(:,:) = zero
     dvr_mat(lm)%Q(1:n_final,1:n_final) = dvr_mat(0)%Inverse(1:n_final,1:n_final)
     R_factor = one / R_max**l_factor
     DO i = 1, n_total
        dvr_mat(lm)%Q(1:n_total,i) = dvr_mat(lm)%Q(1:n_total,i) + R_factor * factr(1:n_total) * factr(i)
     END DO
     write(iout,*) 'End Calculation of Basic Radial Two-Electron Integrals L = '//itoc(lm)
     Call Print_Matrix(type_real_matrix,dvr_mat(lm)%Q,n_total, n_total,                    &
                       title='Q Function for Two Electron Integrals_'//itoc(lm))
  END DO
  DEALLOCATE( dvr_mat(0)%ipvt, dvr_mat(0)%lower, dvr_mat(0)%Inverse, dvr_mat(0)%tr )
END SUBROUTINE V_ijkl_Spherical
!***********************************************************************
!***********************************************************************
!deck V_ijkl_Spheroidal
!***begin prologue     V_ijkl_Spheroidal
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Q function for computation of two-electron integrals.
!***
!***references
!***routines called
!***end prologue       V_ijkl_Spheroidal
  SUBROUTINE V_ijkl_Spheroidal(grid, name_spheroidal)
  IMPLICIT NONE
  TYPE(coordinates)                        :: grid  
  TYPE(spheroidal)                         :: name_spheroidal
  REAL(idp)                                :: scale  
  INTEGER                                  :: info
  INTEGER                                  :: count
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: val
  CHARACTER(LEN=3)                         :: itoc
!
END SUBROUTINE V_ijkl_Spheroidal
!***********************************************************************
!***********************************************************************
!deck D_LM
!***begin prologue     D_LM
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Calculate the D_LM Coefficients
!***description        The D_LM coefficients are products of Wigner 3J coefficients
!***                   (-1)^m * C_3J(l,l',L|0,0,0) * C_3J(l,l',L|-m,m',M)
!***                   and arise from the integation of three spherical harmonics one being
!***                   complex conjugated which leads to the (-1)^m prefactor.
!***references
!***routines called
!***end prologue       D_LM
  SUBROUTINE D_LM
  IMPLICIT NONE
  INTEGER                                  :: l_1
  INTEGER                                  :: l_2
  INTEGER                                  :: m_1
  INTEGER                                  :: m_2  
  INTEGER                                  :: m_fac  
  INTEGER                                  :: l_tot  
  INTEGER                                  :: m_tot
  INTEGER                                  :: l_upper  
  INTEGER                                  :: l_lower
  INTEGER                                  :: count
  INTEGER                                  :: possible
  INTEGER                                  :: maximum_total_l
  REAL(idp)                                :: pre_factor
  REAL(idp)                                :: multiplier
  REAL(idp)                                :: F_3J
  REAL(idp)                                :: coef
  REAL(idp)                                :: value
  DO l_1 = 0, l_max
     DO l_2 = 0, l_1
        IF (  ( l_1 + l_2 ) > maximum_total_l ) THEN
           exit
        ELSE
           l_upper = l_1 + l_2
           l_lower = abs(l_1-l_2)
           pre_factor = sqrt ( dfloat ( (l_1 + l_1 + int_one) * ( l_2 + l_2 + int_one ) ) )
           possible = ( l_1 + l_1 + int_one ) * ( l_2 + l_2 + int_one )
           ALLOCATE( ang_mat(l_tot:l_tot,-l_tot:l_tot) )
           DO l_tot = l_lower, l_upper
              coef = pre_factor * F_3J(l_1,int_zero,l_2,int_zero,l_tot,int_zero,.false.) 
              ALLOCATE( ang_mat(l_tot:l_tot,-l_tot:l_tot) )
              DO m_tot = -l_tot, l_tot
                 ALLOCATE( ang_mat(l_tot,m_tot)%D_LM_Coef(1:possible) )
                 m_fac = (-int_one)**(-l_1)
                 multiplier = m_fac * pre_factor
                 count = int_zero
                 DO m_1 = -l_1, l_1
                    DO m_2 = -l_2, l_2
                       IF ( m_tot /= ( m_2-m_1) ) THEN
                            exit
                       ELSE
                          value = multiplier * F_3J(l_1,-m_1,l_2,m_2,l_tot,m_tot,.false.) 
                          IF (value /= zero) THEN
                              count = count + int_one
                              ang_mat(l_tot,m_tot)%D_LM_Coef(count) = value
                              ang_mat(l_tot,m_tot)%D_LM_Index(count,1) = l_1
                              ang_mat(l_tot,m_tot)%D_LM_Index(count,2) = m_1
                              ang_mat(l_tot,m_tot)%D_LM_Index(count,3) = l_2
                              ang_mat(l_tot,m_tot)%D_LM_Index(count,4) = m_2
                          END IF
                       END IF 
                    END DO
                    m_fac = - m_fac
                 END DO
              END DO
           END DO
        END IF
     END DO
  END DO
END SUBROUTINE D_LM
!***********************************************************************
!***********************************************************************
           END MODULE Two_Electron_FEDVR_Module
!***********************************************************************
!***********************************************************************
