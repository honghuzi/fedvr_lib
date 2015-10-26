!***********************************************************************
! DVR_H_0_Module
!**begin prologue     DVR_H_0_Module
!**date written       090119   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Calculate the FEDVR one body matrix elements
!***description       The raw matrix KE elements are used to assemble
!***                  full one-body Kinetic energy and/or Hamiltonian
!***                  sector matrices.  Thus angular momentum terms
!***                  and any one-body potentials are added in where required.
!***                  
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***                  
!***                  
!***end prologue      DVR_H_0_Module
!***********************************************************************
!***********************************************************************
                           MODULE DVR_H_0_Module
                           USE Data_Module
                           USE FEDVR_Shared
                           USE FEDVR_Derived_Types
                           USE Legendre_Data
                           USE Matrix_Print
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!
                            INTERFACE K_0_Matrix                       
                       MODULE PROCEDURE K_0_DVR,                                &
                                        K_0_DVR_Odd,                            &
                                        K_0_DVR_Fourier,                        &  
                                        K_0_DVR_Hermite,                        &  
                                        K_0_DVR_Laguerre  
                            END INTERFACE K_0_Matrix
!

                            INTERFACE Add_m_Angular_Momentum
                       MODULE PROCEDURE Add_m_Even_Angular_Momentum,            &
                                        Add_m_Odd_Angular_Momentum
                            END INTERFACE Add_m_Angular_Momentum
!
           REAL(idp), Private, DIMENSION(:,:), ALLOCATABLE       :: matrix
           TYPE(REAL_MATRIX)                                     :: type_real_matrix
           TYPE(REAL_VECTOR)                                     :: type_real_vector
!
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck K_0_DVR.f
!***begin prologue     K_0_DVR
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Form the basic nabla matrix elements for even FE regional matrices.
!***                   This corresponds to the l=0 angular momentum matrix.
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE K_0_DVR(grid,reg_mat)
  IMPLICIT NONE
  TYPE (coordinates)                               :: grid
  TYPE (matrices), DIMENSION(:)                    :: reg_mat
  INTEGER                                          :: i
  CHARACTER (LEN=3)                                :: itoc
!
  DO i = 1, nreg
     grid%reg_type_op(i,0)%tr(:,:) = grid%reg_mat(i)%ham(:,:) 
     DEALLOCATE(grid%reg_mat(i)%ham)         
  END DO
  IF (prn(4) == .true. ) THEN
      DO i = 1, nreg
        Call Print_Matrix(type_real_matrix,grid%reg_type_op(i,0)%tr,   &
                          npt(i),npt(i),title='Even matrix before trim Region = '//itoc(i))
      END DO
  END IF
END SUBROUTINE K_0_DVR
!***********************************************************************
!***********************************************************************
!deck K_0_DVR_Odd.f
!***begin prologue     K_0_DVR_Odd
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Form basic nabla matrix elements for odd FE regional matrices.
!***                   This corresponds to the l=1 angular momentum matrix.
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE K_0_DVR_Odd(grid,reg_mat_odd)
  IMPLICIT NONE
  TYPE (coordinates)                               :: grid
  TYPE (odd_matrices), DIMENSION(:)                :: reg_mat_odd
  INTEGER                                          :: i
  CHARACTER (LEN=3)                                :: itoc
!
  DO i = 1, nreg
     grid%reg_type_op(i,1)%tr(:,:) = grid%reg_mat_odd(i)%ham(:,:) 
     DEALLOCATE(grid%reg_mat_odd(i)%ham)         
  END DO
  IF (prn(4) == .true. ) THEN
      DO i = 1, nreg
         Call Print_Matrix(type_real_matrix,grid%reg_type_op(i,1)%tr,  &
                           npt(i),npt(i),                               &
                           title = 'Odd matrix before trim Region = '//itoc(i) )
      END DO
  END IF
END SUBROUTINE K_0_DVR_Odd
!***********************************************************************
!***********************************************************************
!deck K_0_DVR_Fourier.f
!***begin prologue     K_0_DVR_Fourier
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Form basic nabla matrix elements for Fourier DVR.
!***                   
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE K_0_DVR_Fourier(grid,reg_mat_fourier)
  IMPLICIT NONE
  TYPE (coordinates)                               :: grid
  TYPE (fourier_matrices), DIMENSION(:)            :: reg_mat_fourier
  INTEGER                                          :: i
  CHARACTER (LEN=3)                                :: itoc
!
  DO i = 1, nreg
     grid%reg_type_op(i,0)%tr(:,:) = grid%reg_mat_fourier(i)%ham(:,:) 
     DEALLOCATE(grid%reg_mat_fourier(i)%ham) 
  END DO
  IF (prn(4) == .true. ) THEN
      DO i = 1, nreg
         Call Print_Matrix(type_real_matrix,grid%reg_type_op(i,0)%tr,  &
                           npt(i),npt(i),                               &
                           title = 'Fourier matrix before trim Region = '//itoc(i))
      END DO
  END IF 
END SUBROUTINE K_0_DVR_Fourier
!***********************************************************************
!***********************************************************************
!deck K_0_DVR_Hermite.f
!***begin prologue     K_0_DVR_Hermite
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Form basic nabla matrix elements for Hermiten FE regional matrices.
!***                   
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE K_0_DVR_Hermite(grid,reg_mat_hermite)
  IMPLICIT NONE
  TYPE (coordinates)                               :: grid
  TYPE (hermite_matrices), DIMENSION(:)            :: reg_mat_hermite
  INTEGER                                          :: i
  CHARACTER (LEN=3)                                :: itoc
!
 DO i = 1, nreg
     grid%reg_type_op(i,0)%tr(:,:) = grid%reg_mat_hermite(i)%ham(:,:) 
     DEALLOCATE(grid%reg_mat_hermite(i)%ham) 
  END DO
  IF (prn(4) == .true. ) THEN
      DO i = 1, nreg
         Call Print_Matrix(type_real_matrix,grid%reg_type_op(i,0)%tr,  &
                           npt(i),npt(i),                               &
                           title = 'Hermite matrix before trim Region = '//itoc(i) )
      END DO
  END IF
END SUBROUTINE K_0_DVR_Hermite
!***********************************************************************
!***********************************************************************
!deck K_0_DVR_Laguerre.f
!***begin prologue     K_0_DVR_Laguerre
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Form basic nabla matrix elements for Laguerre regional matrices.
!***                   
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE K_0_DVR_Laguerre(grid,reg_mat_laguerre)
  IMPLICIT NONE
  TYPE (coordinates)                               :: grid
  TYPE (laguerre_matrices), DIMENSION(:)           :: reg_mat_laguerre
  INTEGER                                          :: i
  CHARACTER (LEN=3)                                :: itoc
!
  DO i = 1, nreg
     grid%reg_type_op(i,0)%tr(:,:) =  grid%reg_mat_laguerre(i)%ham(:,:) 
     DEALLOCATE(grid%reg_mat_laguerre(1)%ham) 
  END DO 
  IF (prn(4) == .true. ) THEN
      DO i = 1, nreg

         Call Print_MAtrix(type_real_matrix,grid%reg_type_op(i,0)%tr,  &
                          npt(i),npt(i),                                &
                          title = 'Laguerre matrix before trim Region = '//itoc(i) )
      END DO
  END IF
END SUBROUTINE K_0_DVR_Laguerre
!***********************************************************************
!***********************************************************************
!deck Form_Final_Sector_Matrices
!***begin prologue     Form_Final_Sector_Matrices
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Form the sector matrices from the raw sector matrix 
!***                   containing the second derivative term appropriate to the designated
!***                   coordinate plus the potential.  At the end of the routine the full
!***                   sector matrix including thr angular momentum is constructed.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Form_Final_Sector_Matrices

  SUBROUTINE Form_Final_Sector_Matrices(grid, label)
  IMPLICIT NONE
  TYPE (coordinates)                               :: grid
  TYPE (ham)                                       :: type_hamiltonian
  TYPE (nabla)                                     :: type_nabla
  TYPE (even)                                      :: type_even
  TYPE (odd)                                       :: type_odd
  INTEGER                                          :: i
  INTEGER                                          :: j
  INTEGER                                          :: lm
  INTEGER                                          :: l_num
  CHARACTER(LEN=*)                                 :: label
  CHARACTER(LEN=3)                                 :: itoc
!
!     Allocate and fill them
!
  write(iout,*)
  IF ( key_to_allocate == int_zero ) THEN
    write(iout,*) '          If we are dealing with the radial coordinate in spherical coordinates'
    write(iout,*) '          or the m even angular coordinate theta or eta, we need to allocate and fill'
    write(iout,*) '          all of the  matrices grid%reg_type_op(i,lm)%tr with '
    write(iout,*) '          grid%reg_type_op(i,0)%tr'
    DO lm = int_one, lm_max
       DO i = 1, nreg
          ALLOCATE(grid%reg_type_op(i,lm)%tr(1:npt(i),1:npt(i)))
          grid%reg_type_op(i,lm)%tr(:,:) = grid%reg_type_op(i,0)%tr(:,:) ! All can be filled with the even matrix        
       END DO
    END DO
  ELSE IF ( key_to_allocate == int_one ) THEN
    write(iout,*) '          If we are dealing with the m odd angular coordinate in theta or eta,' 
    write(iout,*) '          we need to allocate and fill the odd Tr matrices'
    write(iout,*) '          grid%reg_type_op(i,lm)%tr with grid%reg_type_op(i,1)%tr'
    DO lm = int_two, lm_max, int_two
       DO i = 1, nreg
          ALLOCATE(grid%reg_type_op(i,lm)%tr(1:npt(i),1:npt(i)))
          grid%reg_type_op(i,lm)%tr(:,:) = grid%reg_type_op(i,0)%tr(:,:) ! Fill with the even matrix        
       END DO
    END DO
    DO lm = int_three, lm_max, int_two
       DO i = 1, nreg
          ALLOCATE(grid%reg_type_op(i,lm)%tr(1:npt(i),1:npt(i)))
          grid%reg_type_op(i,lm)%tr(:,:) = grid%reg_type_op(i,1)%tr(:,:) ! Fill with the odd matrix        
       END DO
    END DO
  END IF
!
!    Now add in the angular momentum dealing with all the special cases.
!
  IF(keyword =='spherical') THEN
     IF(label == 'r') THEN     
!
!
        write(iout,*) '          Add the orbital angular momentum to the radial operators'
        Call Add_l_Angular_Momentum ( grid ) 
!
        write(iout,*)'          Finished Construction of Regional Tr Matrices for Coordinate = ',label
        write(iout,*)'          for  All Angular Momenta'      
!
    ELSE IF( label == 'theta') THEN
!
        IF ( m_max >= int_two) THEN
             Call Add_m_Angular_Momentum(grid, type_even, int_two, int_two)
        END IF
        IF ( m_max >= int_one) THEN
             Call Add_m_Angular_Momentum(grid, type_odd, int_one, int_two)
        END IF
        write(iout,*)'          Finished Construction of Regional Tr Matrices for Coordinate = ',label
        write(iout,*)'          for  All Angular Momenta'
     END IF
  ELSE IF(keyword =='cylindrical') THEN
!
        IF(label == 'rho') THEN
           IF (m_max > int_zero ) THEN
               Call Add_m_Angular_Momentum(grid, type_even, int_one, int_one)
           END IF
        END IF
        write(iout,*)'          Finished Construction of Regional Tr Matrices for Coordinate = ',label
        write(iout,*)'          for  All Angular Momenta'
  ELSE IF (keyword == 'spheroidal') THEN
        IF ( m_max >= int_two ) THEN
            Call Add_m_Angular_Momentum(grid, type_even, int_two, int_two)
        END IF
        IF ( m_max >= int_one) THEN
             Call Add_m_Angular_Momentum(grid, type_odd,  int_one, int_two)
        END IF
        write(iout,*)'          Finished Construction of Regional Tr Matrices for Coordinate = ',label
        write(iout,*)'          for  All Angular Momenta'      
  END IF
  IF (prn(4) == .true. ) THEN
      DO lm = int_zero, lm_max
         write(iout,*)
         write (iout,*) '           Results for symmetry = '//itoc(lm)
         write(iout,*)
         DO i = 1, nreg
            Call Print_MAtrix(type_real_matrix,grid%reg_type_op(i,lm)%tr,npt(i),npt(i),              &
                              title = sub_title(2)(1:title_len(2))//' '//itoc(i)//' with proper BC')
         END DO
      END DO
  END IF
!
!       Now form the Hamiltonian by adding in the potential.
!
  dscale = - half * hbar * hbar / mass
  IF(units == 'atomic_units') then
     dscale = - half
  END IF     
!
  IF(reg_pot(i)%type == 'spheroidal') then
     dscale = - R_ab * quarter
  END IF
!
  Call IOsys('write real dscale to '//FEDVR_File,1,dscale,0,' ')
  IF ( form_hamiltonian == .true. ) THEN
       write(iout,*)
       write(iout,*) '          Allocating the Hamiltonian matrices and filling them '
       write(iout,*) '                     with the nabla matrices'
       DO lm = int_zero, lm_max
          DO i = 1, nreg
             ALLOCATE(grid%reg_type_op(i,lm)%ham(1:npt(i),1:npt(i)))  ! Allocate the Hamiltonian matrices
             grid%reg_type_op(i,lm)%ham(:,:) = dscale * grid%reg_type_op(i,lm)%tr(:,:)  ! Fill with the Nabla matrices
                                                                                        ! and scale Nabla 
                                                                                        ! to get Hamiltonian
          END DO
       END DO
       write(iout,*)
       write(iout,*) '          Adding in the potential'
       Call Add_Potential( grid )  ! Add the the potential
       IF (prn(4) == .true. ) THEN
           DO lm = int_zero, lm_max
              write(iout,*)
              write (iout,*) '           Results for symmetry = '//itoc(lm)
              write(iout,*)
              DO i = 1, nreg
                 Call Print_Matrix(type_real_Matrix,grid%reg_type_op(i,lm)%ham,npt(i),npt(i),       &
                                   title = sub_title(1)(1:title_len(1))//' '//itoc(i)//' with proper BC' )
              END DO
           END DO
       END IF
  END IF
  IF ( form_nabla == .false. ) THEN
       DO lm = int_zero, lm_max
          DO i = 1, nreg
             DEALLOCATE(grid%reg_type_op(i,lm)%tr)  ! Deallocate the tr matrices
          END DO
       END DO
  END IF
!
END SUBROUTINE Form_Final_Sector_Matrices
!***********************************************************************
!***********************************************************************
!deck Add_Potential
!***begin prologue     Add_Potential
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Add in the potential
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Add_Potential
  SUBROUTINE Add_Potential( grid )
  IMPLICIT NONE
  TYPE(coordinates)              :: grid
  INTEGER                        :: key
  INTEGER                        :: i
  INTEGER                        :: j
!
  DO key = int_zero, lm_max
     DO i = 1, nreg
        DO j = 1, npt(i)
           grid%reg_type_op(i,key)%ham(j,j) = grid%reg_type_op(i,key)%ham(j,j)       &
                                                      +                              &
                                                 reg_pot(i)%vec(j)                   &
                                                      *                              &
                                              grid%reg_pt_wt(i)%qr_fac(j)
        END DO
     END DO     
  END DO
!
END SUBROUTINE Add_Potential
!***********************************************************************
!***********************************************************************
!deck Add_l_Angular_Momentum
!***begin prologue     Add_l_Angular_Momentum
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Add in the angular momentum
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Add_l_Angular_Momentum
  SUBROUTINE Add_l_Angular_Momentum( grid )
  IMPLICIT NONE
  TYPE(coordinates)              :: grid
  INTEGER                        :: lm
  INTEGER                        :: l_num
  INTEGER                        :: i
  INTEGER                        :: j
  DO lm = int_one, l_max
     l_num = lm * ( lm + int_one )
     DO i = 1, nreg
        DO j = 1, npt(i)
           grid%reg_type_op(i,lm)%tr(j,j) = grid%reg_type_op(i,lm)%tr(j,j)         &
                                          -                                        &
           l_num / ( grid%reg_pt_wt(i)%qr(j) * grid%reg_pt_wt(i)%qr(j) )           &
                                          *                                        &
                         grid%reg_pt_wt(i)%qr_fac(j)
        END DO
     END DO
  END DO
!
END SUBROUTINE Add_l_Angular_Momentum
!***********************************************************************
!***********************************************************************
!deck Add_m_Even_Angular_Momentum
!***begin prologue     Add_m_Even_Angular_Momentum
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Add in the m angular momentum to the nabla.
!***                   The matrix with m=zero is already formed.
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Add_m_Even_Angular_Momentum
  SUBROUTINE Add_m_Even_Angular_Momentum(grid, type_even, m_start, m_skip)
  IMPLICIT NONE
  TYPE (coordinates)                               :: grid
  TYPE (even)                                      :: type_even
  INTEGER                                          :: m
  INTEGER                                          :: m_num
  INTEGER                                          :: m_start
  INTEGER                                          :: m_skip
  INTEGER                                          :: ir
  INTEGER                                          :: i
!
  type_even%kind = 'even m'
  DO m = m_start, m_max, m_skip
     m_num = m * m
     DO ir = 1, nreg
        DO i = 1, npt(ir)
           grid%reg_type_op(ir,m)%tr(i,i) = grid%reg_type_op(ir,m)%tr(i,i)                   &
                                                      -                                      &
                                            m_num * grid%reg_pt_wt(ir)%inv_qr_fac(i)
        END DO 
     END DO
  END DO
END SUBROUTINE Add_m_Even_Angular_Momentum
!***********************************************************************
!***********************************************************************
!deck Add_m_Odd_Angular_Momentum
!***begin prologue     Add_m_Odd_Angular_Momentum
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Add in the m angular momentum to the nabla.
!***                   The matrix with m=one is already formed.
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Add_m_Odd_Angular_Momentum
  SUBROUTINE Add_m_Odd_Angular_Momentum(grid, type_odd, m_start, m_skip)
  IMPLICIT NONE
  TYPE (coordinates)                               :: grid
  TYPE (odd)                                       :: type_odd
  INTEGER                                          :: ir
  INTEGER                                          :: i
  INTEGER                                          :: m
  INTEGER                                          :: m_num
  INTEGER                                          :: m_start
  INTEGER                                          :: m_skip
!
  type_odd%kind = 'odd m'
  DO m = m_start, m_max, m_skip
     m_num = m * m
     DO ir = 1, nreg
        DO i = 1, npt(ir)
           grid%reg_type_op(ir,m)%tr(i,i) = grid%reg_type_op(ir,m)%tr(i,i)                   &
                                                      -                                      &
                                            m_num * grid%reg_pt_wt(ir)%inv_qr_fac(i)
        END DO 
     END DO
  END DO
END SUBROUTINE Add_m_Odd_Angular_Momentum
!***********************************************************************
!***********************************************************************
!deck Fix_BC
!***begin prologue     Fix_BC
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Fixing boundary conditions by removing certain basis functions
!***                   from either the first or last finite element.
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Fix_BC

  SUBROUTINE Fix_BC(grid)
  IMPLICIT NONE
  TYPE (coordinates)                               :: grid
  INTEGER                                          :: lm
  INTEGER                                          :: n_max
  INTEGER                                          :: upper
  INTEGER, DIMENSION(2)                            :: from
  INTEGER                                          :: id
  INTEGER, DIMENSION(2)                            :: reg
  CHARACTER (LEN=3)                                :: itoc
!
  write(iout,*)
  write(iout,*)
  write(iout,*) '          **** Enforcing the Boundary Conditions on the Sector Matrices ****'
  n_max = max(int_zero,npt(1),npt(nreg))
  upper = int_zero
  IF (m_max > int_zero ) THEN
      upper = int_one
  END IF
  ALLOCATE(matrix(1:n_max,1:n_max), n_tmp(1:nreg))
  n_tmp(1:nreg) = npt(1:nreg)  ! Use n_tmp to store npt temporarily
  reg(1) = int_one
  reg(2) = nreg
  from(1) = int_two
  from(2) = int_one
  DO id = int_one, int_two
     IF ( drop(id) == .true.) THEN
          write(iout,1) reg(id)
          n_tmp(reg(id)) = n_tmp(reg(id)) - int_one
          DO lm = int_zero, upper
             write(iout,2) from(id),n_tmp(reg(id))
             matrix(1:n_tmp(reg(id)),1:n_tmp(reg(id)))                                     &
                                             =                                             &
             grid%reg_type_op(reg(id),lm)%tr(from(id):npt(reg(id)),from(id):npt(reg(id)))
             DEALLOCATE(grid%reg_type_op(reg(id),lm)%tr)
             ALLOCATE(grid%reg_type_op(reg(id),lm)%tr(1:n_tmp(reg(id)),1:n_tmp(reg(id))))
             grid%reg_type_op(reg(id),lm)%tr(1:n_tmp(reg(id)),1:n_tmp(reg(id)))            &
                                             =                                             &
             matrix(1:n_tmp(reg(id)),1:n_tmp(reg(id)))
          END DO
          IF (prn(4) == .true. ) THEN
              DO lm = 0, upper
                 Call Print_Matrix(type_real_matrix,grid%reg_type_op(reg(id),lm)%tr,n_tmp(reg(id)),n_tmp(reg(id)),  &
                                   title='region = '//itoc(reg(id))//' basic LM = '//itoc(lm)//' with proper boundary conditions')
              END DO
          END IF
          write(iout,3)
          grid%reg_pt_wt(reg(id))%qr(1:n_tmp(reg(id)))              =                                            &       
                                                           grid%reg_pt_wt(reg(id))%qr(from(id):npt(reg(id)))
          grid%reg_pt_wt(reg(id))%wtr(1:n_tmp(reg(id)))              =                                           &       
                                                           grid%reg_pt_wt(reg(id))%wtr(from(id):npt(reg(id)))
          grid%reg_pt_wt(reg(id))%qr_fac(1:n_tmp(reg(id)))          =                                            &
                                                           grid%reg_pt_wt(reg(id))%qr_fac(from(id):npt(reg(id)))
          grid%reg_pt_wt(reg(id))%inv_qr_fac(1:n_tmp(reg(id)))      =                                            &
                                                           grid%reg_pt_wt(reg(id))%inv_qr_fac(from(id):npt(reg(id)))
          grid%reg_pt_wt(reg(id))%inv_sqrt_qr_fac(1:n_tmp(reg(id))) =                                            &
                                                           grid%reg_pt_wt(reg(id))%inv_sqrt_qr_fac(from(id):npt(reg(id)))
          reg_pot(reg(id))%vec(1:n_tmp(reg(id))) = reg_pot(reg(id))%vec(from(id):npt(reg(id)))
     END IF 
  END DO
  npt(1:nreg) = n_tmp(1:nreg)
  DEALLOCATE(matrix,n_tmp)
1 Format(/,20x,'Trimming Region = ', i3)
2 Format(/,20x,'New Matrix Extends from ',i3,' to ',i3 )
3 Format(/,20x,'Trimming some additional vectors needed to form the final and composite matrices')
END SUBROUTINE Fix_BC
!***********************************************************************
!***********************************************************************
!deck PE_DVR_Matrix.f
!***begin prologue     PE_DVR_Matrix
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Calculate the local potential matrix for a number of
!***                   fairly standard potentials.
!***                   
!***                   
!***                   
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       PE_DVR_Matrix
  SUBROUTINE PE_DVR_Matrix(grid)
  IMPLICIT NONE
  TYPE (coordinates)                           :: grid
  INTEGER                                      :: i
  INTEGER                                      :: intkey
  CHARACTER (LEN=3)                            :: itoc
  REAL(idp)                                    :: fpkey
!
!
  DO i = 1, nreg
     ALLOCATE( reg_pot(i)%vec(npt(i)) )
     reg_pot(i)%vec(:) = zero     
!
     IF(reg_pot(i)%type == 'none') then
!
        call none(reg_pot(i)%vec,                                   &
                  grid%reg_pt_wt(i)%qr,                             &
                  npt(i),                                           &
                  prnt)
!
     ELSE IF(reg_pot(i)%type == 'well') then
!
        reg_pot(i)%well%depth=fpkey(card,'well_depth',zero,' ')
        call vwell(reg_pot(i)%vec,                                  &
                   reg_pot(i)%well%depth,                           &
                   npt(i),                                          &
                   prnt)
!
     ELSE IF(reg_pot(i)%type == 'exponential') then
!
        reg_pot(i)%exponential%amp(1)=fpkey(card,'amplitude',-1.d0,' ')
        reg_pot(i)%exponential%expnt(1)=fpkey(card,'exponent',reg_pot(i)%exponential%expnt,' ')
        call vexp(reg_pot(i)%vec,                                   &
                  grid%reg_pt_wt(i)%qr,                             &
                  reg_pot(i)%exponential%amp,                       &
                  reg_pot(i)%exponential%expnt,                     &
                  npt(i),                                           &
                  prnt)
!
     ELSE IF(reg_pot(i)%type == 'yukawa') then
!
        reg_pot(i)%yukawa%amp(1)=fpkey(card,'amplitude',-1.d0,' ')
        reg_pot(i)%yukawa%expnt(1)=fpkey(card,'exponent',reg_pot(i)%exponential%expnt,' ')
        call vyukawa(reg_pot(i)%vec,                                &
                     grid%reg_pt_wt(i)%qr,                          &
                     reg_pot(i)%yukawa%amp,                         &
                     reg_pot(i)%yukawa%expnt,                       &
                     npt(i),                                        &
                     prnt)
!
     ELSE IF(reg_pot(i)%type == 'power_exponential') then
!
        reg_pot(i)%power_exp%amp(1)=fpkey(card,'amplitude',1.d0,' ')
        reg_pot(i)%power_exp%expnt(1)=fpkey(card,'exponent',reg_pot(i)%power_exp%expnt,' ')
        reg_pot(i)%power_exp%n_p=intkey(card,'power',0,' ')
        call v_pow_exp(reg_pot(i)%vec,                              &
                       grid%reg_pt_wt(i)%qr,                        &
                       reg_pot(i)%power_exp%amp,                    &
                       reg_pot(i)%power_exp%expnt,                  &
                       reg_pot(i)%power_exp%n_p,                    &
                       npt(i),                                      &
                       prnt)
!
     ELSE IF(reg_pot(i)%type == 'sum_exponential') then
!
        call fparr(card,'amplitudes',reg_pot(i)%sum_exp%amp,2,' ')
        call fparr(card,'exponents',reg_pot(i)%sum_exp%expnt,2,' ')
        call vexp_sum(reg_pot(i)%vec,                               &
                      grid%reg_pt_wt(i)%qr,                         &
                      reg_pot(i)%sum_exp%amp,                       &
                      reg_pot(i)%sum_exp%expnt,                     &
                      npt(i),                                       &
                      prnt)
!
     ELSE IF(reg_pot(i)%type == 'coulomb') then
!
        reg_pot(i)%coulomb%charge=fpkey(card,'charge',-one,' ')
        call vcoul(reg_pot(i)%vec,                                  &
                   grid%reg_pt_wt(i)%qr,                            &
                   reg_pot(i)%coulomb%charge,                       &
                   npt(i),                                          &
                   prnt)
!
     ELSE IF(reg_pot(i)%type == 'eberlonium') then
!
        reg_pot(i)%eberlonium%charge=fpkey(card,'charge',-one,' ')
        reg_pot(i)%eberlonium%n_p=intkey(card,'power',0,' ')
        reg_pot(i)%eberlonium%amp(1)=fpkey(card,'a',1.d0,' ')
        reg_pot(i)%eberlonium%amp(2)=fpkey(card,'b',1.d0,' ')
        call v_eberlonium(reg_pot(i)%vec,                           &
                          grid%reg_pt_wt(i)%qr,                     &
                          reg_pot(i)%eberlonium%charge,             &
                          reg_pot(i)%eberlonium%amp(1),             &
                          reg_pot(i)%eberlonium%amp(2),             &
                          reg_pot(i)%eberlonium%n_p,                &
                          npt(i),                                   &
                          prnt)
!
     ELSE IF(reg_pot(i)%type == 'inverse_r4') then
!
        call vir4(reg_pot(i)%vec,                                   &
                  grid%reg_pt_wt(i)%qr,                             &
                  npt(i),                                           &
                  prnt)
!
     ELSE IF(reg_pot(i)%type == 'rounded_well') then
!
        reg_pot(i)%well%nwell=intkey(card,'n_well',ten,' ')
        reg_pot(i)%well%awell=fpkey(card,'a_well',14.d0,' ')
        call vrwell(reg_pot(i)%vec,                                 &
                    grid%reg_pt_wt(i)%qr,                           &
                    reg_pot(i)%well%awell,                          &
                    reg_pot(i)%well%nwell,                          &
                    npt(i),                                         &
                    prnt)
!
     ELSE IF(reg_pot(i)%type == 'harmonic_oscillator') then
!
!        enter the mass and frequency in atomic units.
!         
        reg_pot(i)%harmonic_oscillator%mass=fpkey(card,'mass',1.d0,' ')
        reg_pot(i)%harmonic_oscillator%omega=fpkey(card,'omega',1.d0,' ')
        factor= reg_pot(i)%harmonic_oscillator%mass * reg_pot(i)%harmonic_oscillator%omega &
                                                    * reg_pot(i)%harmonic_oscillator%omega &
                                                    * half
        hbar=one
        write(iout,1) mass, reg_pot(i)%harmonic_oscillator%omega
        call vhmo(reg_pot(i)%vec,                                   &
                  grid%reg_pt_wt(i)%qr,                             &
                  factor,                                           &
                  npt(i),                                           &
                  prnt)
!
     ELSE IF(reg_pot(i)%type == 'anharmonic_oscillator') then
!
        call vanhmo(reg_pot(i)%vec,                                 &
                    grid%reg_pt_wt(i)%qr,                           &
                    npt(i),                                         &
                    prnt)
!
     ELSE IF(reg_pot(i)%type == 'expres') then
!
        call fparr(card,'amplitude',reg_pot(i)%expres%amp,2,' ')
        call fparr(card,'exponent',reg_pot(i)%expres%expnt,2,' ')
        reg_pot(i)%expres%shift=fpkey(card,'exponent_shift',zero,' ')
        call vres(reg_pot(i)%vec,                                   &
                  grid%reg_pt_wt(i)%qr,                             &
                  reg_pot(i)%expres%amp,                            &
                  reg_pot(i)%expres%expnt,                          &
                  reg_pot(i)%expres%shift,                          &
                  npt(i),                                           &
                  prnt)
!
     ELSE IF(reg_pot(i)%type == 'periodic') then
!
        reg_pot(i)%periodic%n_scale=intkey(card,'n_i',10,' ')         
        reg_pot(i)%periodic%e_c=fpkey(card,'e_c',.001d0,' ')         
        reg_pot(i)%periodic%awell=reg_pot(i)%periodic%n_scale/reg_pot(i)%periodic%e_c
        call vperiod(reg_pot(i)%vec,                                &
                     grid%reg_pt_wt(i)%qr,                          &
                     reg_pot(i)%periodic%awell,                     &
                     npt(i),                                        &
                     prnt)
!
     ELSE IF(reg_pot(i)%type == 'spheroidal_coulomb') then
!
!    In the spheroidal case, this is not actually the potential but a scaled version which
!    takes into account the volume factor in the integral which is non-separable in the two
!    coordinates.
!
        IF ( grid%label == 'eta') THEN
             reg_pot(i)%vec(:) = R_ab * ( Z_b - Z_a ) * grid%reg_pt_wt(i)%qr(:) * grid%reg_pt_wt(i)%inv_qr_fac(:)
        ELSE IF ( grid%label == 'xi' ) THEN
             reg_pot(i)%vec(:) = R_AB * ( Z_a + Z_b ) * grid%reg_pt_wt(i)%qr(:) * grid%reg_pt_wt(i)%inv_qr_fac(:)
        END IF
     ELSE IF(reg_pot(i)%type == 'general_spheroidal') then
             reg_pot(i)%vec(:) = (                                                                                &
                                   reg_pot(i)%c(0)                                                                &
                                                  +                                                               &
                                   reg_pot(i)%c(1) * grid%reg_pt_wt(i)%qr(:)                                      &
                                                  +                                                               &
                                   reg_pot(i)%c(2) * grid%reg_pt_wt(i)%qr(:) * grid%reg_pt_wt(i)%qr(:)            &
                                                                                                        )         & 
                                                  *                                                               &
                                   grid%reg_pt_wt(i)%inv_qr_fac(:)
     ELSE
        Call lnkerr('screwed up potential. quit')
     END IF
!
  END DO
  IF(prn(3)) then
     DO i = 1, nreg
        write(iout,2) reg_pot(i)%type

        call Print_Matrix(type_real_vector,reg_pot(i)%vec,    &
                          title='potential matrix elements for region = '//itoc(i) )
     END DO
  END IF
!
!
 1 Format(/,1x,'oscillator mass      = ',e15.8, &
          /,1x,'oscillator-frequency = ',e15.8)
 2 FORMAT(/,20x,'Potential Type = ',a32)
!
END SUBROUTINE PE_DVR_Matrix
!***********************************************************************
!***********************************************************************
           END MODULE DVR_H_0_Module
!***********************************************************************
!***********************************************************************
