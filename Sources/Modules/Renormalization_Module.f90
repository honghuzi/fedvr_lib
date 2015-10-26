!***********************************************************************
!Renormalization_Module
!**begin prologue     Renormalization_Module
!**date written       090119   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Compute renormalization factors required for the
!***                  renormalization of the sector polynomials and kinetic
!***                  matrix elements in a FEDVR basis.
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***end prologue      Renormalization_Module
!***********************************************************************
!***********************************************************************
                           MODULE Renormalization_Module
                           USE Data_Module
                           USE FEDVR_Shared
                           USE FEDVR_Derived_Types
                           USE Matrix_Print
  TYPE(REAL_MATRIX)                              :: type_real_matrix
  TYPE(REAL_VECTOR)                              :: type_real_vector
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!
                            INTERFACE Renormalization
                       MODULE PROCEDURE Cartesian_Renormalization,       &
                                        Odd_Renormalization,             &
                                        Laguerre_Renormalization,        &
                                        Hermite_Renormalization,         &
                                        Fourier_Renormalization
                            END INTERFACE Renormalization                                     
!
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck Norm.f
!***begin prologue     Norm
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Modify the weights at the ends of the interval to
!***                   reflect that there are bridge functions present. 
!***                   
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Norm

  SUBROUTINE Norm(wt,normalization,n,region,wt_left_end,wt_right_end)
  IMPLICIT NONE
  REAL(idp), DIMENSION (:)             :: wt
  REAL(idp), DIMENSION (:)             :: normalization
  INTEGER                              :: n
  INTEGER                              :: region
  REAL(idp), OPTIONAL                  :: wt_right_end
  REAL(idp), OPTIONAL                  :: wt_left_end
!
!
  IF ( .not.present(wt_left_end).and..not.present(wt_right_end) ) THEN
!
       normalization(:) = Sqrt ( 1.d0 / wt(:) )
!
  ELSE IF ( .not.present(wt_left_end).and.present(wt_right_end) ) THEN
!
!      Modify the last weight and then get the inverse square roots.
!
       normalization(1:n-1) = 1.d0 / sqrt ( wt(:) )     
       normalization(n) = 1.d0 / sqrt ( wt(n) + wt_right_end )     
!
  ELSE IF ( present(wt_left_end).and..not.present(wt_right_end) ) THEN
!
!      
       normalization(1) = 1.d0 / sqrt ( wt_left_end + wt(1) )
       normalization(2:n) = 1.d0 / sqrt ( wt(2:n) )
!
!
  ELSE IF ( present(wt_left_end).and.present(wt_right_end) ) THEN
!      
!      Modify the last weight and then get the inverse square roots.
!
       normalization(1) = 1.d0 / Sqrt ( wt_left_end + wt(1) )
       normalization(n)  = 1.d0 / Sqrt ( wt_right_end + wt(n) )
       normalization(2:n-1)  =  1.d0 / sqrt ( wt(2:n-1) )
  END IF
END SUBROUTINE Norm
!***********************************************************************
!***********************************************************************
!deck Re_Poly.f
!***begin prologue     Re_Poly
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            normalize the lobatto dvr functions
!***                   again, this does not account for any
!***                   difference due to (l,m)
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Re_Poly
  SUBROUTINE Re_Poly(p,dp,ddp,inv_sqrt_wt,n)
  IMPLICIT NONE
  REAL(idp), DIMENSION(:,:)       :: p
  REAL(idp), DIMENSION(:,:)       :: dp
  REAL(idp), DIMENSION(:,:)       :: ddp
  REAL(idp), DIMENSION (:)        :: inv_sqrt_wt
  INTEGER                         :: n
  INTEGER                         :: i
!
!
  DO i = 1, n
       p(:, i ) =    p(:, i )  * inv_sqrt_wt(i)
      dp(:, i ) =   dp(:, i )  * inv_sqrt_wt(i)
     ddp(:, i ) =  ddp(:, i )  * inv_sqrt_wt(i)
  END DO
END SUBROUTINE Re_Poly
!***********************************************************************
!***********************************************************************
!deck Cartesian_Renormalization.f
!***begin prologue     Cartesian_Renormalization
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Renormalize the kinetic energy matrix elements so that they
!***                   are correct across the joining regions.

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Cartesian_Renormalization
  SUBROUTINE Cartesian_Renormalization(grid,reg_mat)
  IMPLICIT NONE
  TYPE (coordinates)                  :: grid
  TYPE(matrices),  DIMENSION(:)       :: reg_mat
  INTEGER                             :: i
!
!
  IF ( nreg == 1 ) THEN
!
!      One Region only.  No change needed except normalization.
!
       i = 1
       ALLOCATE(grid%reg_mat(1)%ham(npt(1),npt(1)))
!
       Call Re_KE (grid%reg_mat(i)%ham,                                       &
                   grid%reg_mat(i)%tr,                                        &
                   grid%reg_poly(i)%normalization,                            &
                   npt(i),                                                    &
                   i )
       IF (prn(4) == .true. ) THEN
           write(iout,1) i
           Call Print_Matrix(type_real_matrix,grid%reg_mat(i)%ham,npt(i),npt(i),              &
                             title='renormalized KE even FEDVR matrix')
       END IF
  ELSE
       i = 1
!
!      Region one is a special case.  Only the last right hand element needs to be changed
!      and the elements renormalized.
!
       ALLOCATE(grid%reg_mat(i)%ham(npt(i),npt(i)))
       Call Re_KE (grid%reg_mat(i)%ham,                                       &
                   grid%reg_mat(i)%tr,                                        &
                   grid%reg_poly(i)%normalization,                            &
                   npt(i),                                                    &
                   i,                                                         &
                   right_end_mat_el=grid%reg_mat(i+1)%tr(1,1) )
       IF (prn(4) == .true. ) THEN
           write(iout,1) i
           Call Print_Matrix(type_real_matrix,grid%reg_mat(i)%ham,npt(i),npt(i),   &
                             title='renormalized KE even FEDVR matrix')
       END IF
!
!      Now do the general case.  Both the left and right elemenents are changed.
! 
       DO i = 2 , nreg - 1
!
          ALLOCATE(grid%reg_mat(i)%ham(npt(i),npt(i)))
          Call Re_KE (grid%reg_mat(i)%ham,                                     &
                      grid%reg_mat(i)%tr,                                      &
                      grid%reg_poly(i)%normalization,                          &
                      npt(i),                                                  &
                      i,                                                       &
                      left_end_mat_el=grid%reg_mat(i-1)%tr(npt(i-1),npt(i-1)), &
                      right_end_mat_el=grid%reg_mat(i+1)%tr(1,1) )
          IF (prn(4) == .true. ) THEN
              write(iout,1) i

             Call Print_Matrix(type_real_matrix,grid%reg_mat(i)%ham,npt(i),npt(i),           &
                               title='renormalized KE even FEDVR matrix')
          END IF
       END DO
       i = nreg
!
!      Special case of last region.  Only the left hand elements needs to be changed.
!
!
          ALLOCATE(grid%reg_mat(i)%ham(npt(i),npt(i)))
          Call Re_KE (grid%reg_mat(i)%ham,                                     &
                      grid%reg_mat(i)%tr,                                      &
                      grid%reg_poly(i)%normalization,                          &
                      npt(i),                                                  &
                      i,                                                       &
                      left_end_mat_el=grid%reg_mat(i-1)%tr(npt(i-1),npt(i-1)) )
          IF (prn(4) == .true. ) THEN
              write(iout,1) i
             Call Print_Matrix(type_real_matrix,grid%reg_mat(i)%ham,npt(i),npt(i),           &
                               title='renormalized KE FEDVR matrix')
          END IF
  END IF
  DO i = 1, nreg
     DEALLOCATE(grid%reg_mat(i)%tr)
  END DO
1 Format(/,10x,'Region = ',i4)
END SUBROUTINE Cartesian_Renormalization
!***********************************************************************
!***********************************************************************
!deck Odd_Renormalization.f
!***begin prologue     Odd_Renormalization
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Renormalization
!***                   
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Odd_Renormalization

  SUBROUTINE Odd_Renormalization(grid,reg_mat_odd)
  IMPLICIT NONE
  TYPE(coordinates)                    :: grid  
  TYPE(odd_matrices),  DIMENSION(:)    :: reg_mat_odd
  INTEGER                              :: i
!
!
!
!

  IF ( nreg == 1 ) THEN
!
!      One Region only.  No change needed except normalization.
!
       i = 1
       ALLOCATE(grid%reg_mat_odd(i)%ham(npt(i),npt(i)))
!
       Call Re_KE (grid%reg_mat_odd(i)%ham,                                         &
                   grid%reg_mat_odd(i)%tr,                                          &
                   grid%reg_poly(i)%normalization,                                  &
                   npt(i),                                                          &
                   i )
       write(iout,1) i
       IF (prn(4) == .true. ) THEN
           title='renormalized KE odd FEDVR matrix'
           Call Print_Matrix(type_real_matrix,grid%reg_mat_odd(i)%ham,npt(i),npt(i), &
                             title='renormalized KE odd FEDVR matrix')
       END IF
  ELSE
       i = 1
!
!      Region one is a special case.  Only the last right hand element needs to be changed
!      and the elements renormalized.
!
       ALLOCATE(grid%reg_mat_odd(i)%ham(npt(i),npt(i)))
       Call Re_KE (grid%reg_mat_odd(i)%ham,                                         &
                   grid%reg_mat_odd(i)%tr,                                          &
                   grid%reg_poly(i)%normalization,                                  &
                   npt(i),                                                          &
                   i,                                                               &
                   right_end_mat_el=grid%reg_mat_odd(i+1)%tr(1,1) )
       write(iout,1) i
       IF (prn(4) == .true. ) THEN

           Call Print_Matrix(type_real_matrix,grid%reg_mat_odd(i)%ham,npt(i),npt(i), &
                             title='renormalized KE odd FEDVR matrix')     
       END IF
!
!      Now do the general case.  Both the left and right elemenents are changed.
! 
       DO i = 2 , nreg - 1
!
          ALLOCATE(grid%reg_mat_odd(i)%ham(npt(i),npt(i)))
          Call Re_KE (grid%reg_mat_odd(i)%ham,                                      &
                      grid%reg_mat_odd(i)%tr,                                       &
                      grid%reg_poly(i)%normalization,                               &
                      npt(i),                                                       &
                      i,                                                            &
                      left_end_mat_el=grid%reg_mat_odd(i-1)%tr(npt(i-1),npt(i-1)),  &
                      right_end_mat_el=grid%reg_mat_odd(i+1)%tr(1,1) )
          write(iout,1) i
          IF (prn(4) == .true. ) THEN
              Call Print_Matrix(type_real_matrix,grid%reg_mat_odd(i)%ham,npt(i),npt(i), &
                                title='renormalized KE FEDVR matrix')
          END IF
       END DO
       i = nreg
!
!      Special case of last region.  Only the left hand elements needs to be changed.
!
!
          ALLOCATE(grid%reg_mat_odd(i)%ham(npt(i),npt(i)))
          Call Re_KE (grid%reg_mat_odd(i)%ham,                                      &
                      grid%reg_mat_odd(i)%tr,                                       &
                      grid%reg_poly(i)%normalization,                               &
                      npt(i),                                                       &
                      i,                                                            &
                      left_end_mat_el=grid%reg_mat_odd(i-1)%tr(npt(i-1),npt(i-1)) )
          IF (prn(4) == .true. ) THEN
              Call Print_Matrix(type_real_matrix,grid%reg_mat_odd(i)%ham,npt(i),npt(i), &
                                title='renormalized KE odd FEDVR matrix')
          END IF
  END IF
  DO i = 1, nreg
     DEALLOCATE(grid%reg_mat_odd(i)%tr)
  END DO
1 Format(/,10x,'Region = ',i4)
END SUBROUTINE Odd_Renormalization
!***********************************************************************
!***********************************************************************
!deck Fourier_Renormalization.f
!***begin prologue     Fourier_Renormalization
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Renormalization.
!***

!***references
!***routines called    iosys, util and mdutil
!***end prologue       Fourier_Renormalization
  SUBROUTINE Fourier_Renormalization(grid,reg_mat_fourier)
  IMPLICIT NONE  
  TYPE (coordinates)                        :: grid
  TYPE(fourier_matrices),  DIMENSION(:)     :: reg_mat_fourier
  INTEGER                                   :: i
!
!
  DO i = 1, nreg
!
       ALLOCATE(grid%reg_mat_fourier(i)%ham(1:npt(i),1:npt(i)))
       grid%reg_mat_fourier(i)%ham(:,:) = grid%reg_mat_fourier(i)%tr(:,:)
       DEALLOCATE(grid%reg_mat_fourier(i)%tr)
!
  END DO
  IF (prn(4) == .true. ) THEN
      DO i = 1, nreg
         write(iout,1) i
         Call Print_Matrix(type_real_matrix,grid%reg_mat_fourier(i)%ham,npt(i),npt(i), &
                           title='renormalized fourier kinetic energy mat')
      END DO
  END IF
1 Format(/,10x,'Region = ',i4)
END SUBROUTINE Fourier_Renormalization
!***********************************************************************
!***********************************************************************
!deck Laguerre_Renormalization.f
!***begin prologue     Laguerre_Renormalization
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Renormalization.
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Laguerre_Renormalization
  SUBROUTINE Laguerre_Renormalization(grid,reg_mat_laguerre)
  IMPLICIT NONE  
  TYPE (coordinates)                        :: grid
  TYPE(laguerre_matrices),  DIMENSION(:)    :: reg_mat_laguerre
  INTEGER                                   :: i
!
!

  IF ( nreg == 1 ) THEN
!
!      One Region only.  No change needed except normalization.
!
       i = 1
       ALLOCATE(grid%reg_mat_laguerre(i)%ham(npt(i),npt(i)))
!
       Call Re_KE (grid%reg_mat_laguerre(i)%ham,                                       &
                   grid%reg_mat_laguerre(i)%tr,                                        &
                   grid%reg_poly(i)%normalization,                                     &
                   npt(i),                                                             &
                   i )
       write(iout,1) i
       IF (prn(4) == .true. ) THEN
           Call Print_Matrix(type_real_matrix,grid%reg_mat_laguerre(i)%ham,npt(i),npt(i), &
                             title='renormalized laguerre kinetic energy mat')
       END IF
  ELSE
       i = 1
!
!      Region one is a special case.  Only the last right hand element needs to be changed
!      and the elements renormalized.
!
       ALLOCATE(grid%reg_mat_laguerre(i)%ham(npt(i),npt(i)))
       Call Re_KE (grid%reg_mat_laguerre(i)%ham,                                       &
                   grid%reg_mat_laguerre(i)%tr,                                        &
                   grid%reg_poly(i)%normalization,                                     &
                   npt(i),                                                             &
                   i,                                                                  &
                   right_end_mat_el=grid%reg_mat_laguerre(i+1)%tr(1,1) )
       write(iout,1) i
       IF (prn(4) == .true. ) THEN
           Call Print_Matrix(type_real_matrix,grid%reg_mat_laguerre(i)%ham,npt(i),npt(i), &
                             title='renormalized laguerre kinetic energy mat')
       END IF
!
!      Now do the general case.  Both the left and right elemenents are changed.
! 
       DO i = 2 , nreg - 1
!
          ALLOCATE(grid%reg_mat_laguerre(i)%ham(npt(i),npt(i)))
          Call Re_KE (grid%reg_mat_laguerre(i)%ham,                                     &
                      grid%reg_mat_laguerre(i)%tr,                                      &
                      grid%reg_poly(i)%normalization,                                   &
                      npt(i),                                                           &
                      i,                                                                &
                      left_end_mat_el=grid%reg_mat(i-1)%tr(npt(i-1),npt(i-1)),          &
                      right_end_mat_el=grid%reg_mat(i+1)%tr(1,1) )
          write(iout,1) i
          IF (prn(4) == .true. ) THEN
           Call Print_Matrix(type_real_matrix,grid%reg_mat_laguerre(i)%ham,npt(i),npt(i), &
                             title='renormalized laguerre kinetic energy mat')
          END IF
       END DO
       i = nreg
!
!      Special case of last region.  Only the left hand elements needs to be changed.
!
!
          ALLOCATE(grid%reg_mat(i)%ham(npt(i),npt(i)))
          Call Re_KE (grid%reg_mat_laguerre(i)%ham,                                     &
                      grid%reg_mat_laguerre(i)%tr,                                      &
                      grid%reg_poly(i)%normalization,                                   &
                      npt(i),                                                           &
                      i,                                                                &
                      left_end_mat_el=grid%reg_mat_laguerre(i-1)%tr(npt(i-1),npt(i-1)) )
          IF (prn(4) == .true. ) THEN
              write(iout,1) i
           Call Print_Matrix(type_real_matrix,grid%reg_mat_laguerre(i)%ham,npt(i),npt(i), &
                             title='renormalized laguerre kinetic energy mat')
          END IF
  END IF
  DO i = 1, nreg
     DEALLOCATE(grid%reg_mat_laguerre(i)%tr)
  END DO
1 Format(/,10x,'Region = ',i4)
END SUBROUTINE Laguerre_Renormalization
!***********************************************************************
!***********************************************************************
!deck Hermite_Renormalization.f
!***begin prologue     Hermite_Renormalization
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Renormalization.
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Hermite_Renormalization
  SUBROUTINE Hermite_Renormalization(grid,reg_mat_hermite)
  IMPLICIT NONE  
  TYPE (coordinates)                        :: grid
  TYPE(hermite_matrices),  DIMENSION(:)     :: reg_mat_hermite
  INTEGER                                   :: i
!
!

  IF ( nreg == 1 ) THEN
!
!      One Region only.  No change needed except normalization.
!
       i = 1
       ALLOCATE(grid%reg_mat_hermite(i)%ham(npt(i),npt(i)))
!
       Call Re_KE (grid%reg_mat_hermite(i)%ham,                                       &
                   grid%reg_mat_hermite(i)%tr,                                        &
                   grid%reg_poly(i)%normalization,                                     &
                   npt(i),                                                             &
                   i )
       IF (prn(4) == .true. ) THEN
           write(iout,1) i
           Call Print_Matrix(type_real_matrix,grid%reg_mat_hermite(i)%ham,npt(i),npt(i), &
                            title='renormalized hermite kinetic energy mat')
       END IF
  ELSE
       i = 1
!
!      Region one is a special case.  Only the last right hand element needs to be changed
!      and the elements renormalized.
!
       ALLOCATE(grid%reg_mat_hermite(i)%ham(npt(i),npt(i)))
       Call Re_KE (grid%reg_mat_hermite(i)%ham,                                       &
                   grid%reg_mat_hermite(i)%tr,                                        &
                   grid%reg_poly(i)%normalization,                                    &
                   npt(i),                                                            &
                   i,                                                                 &
                   right_end_mat_el=grid%reg_mat_hermite(i+1)%tr(1,1) )
       IF (prn(4) == .true. ) THEN
           write(iout,1) i
           Call Print_Matrix(type_real_matrix,grid%reg_mat_hermite(i)%ham,npt(i),npt(i), &
                            title='renormalized hermite kinetic energy mat')
       END IF
!
!      Now do the general case.  Both the left and right elemenents are changed.
! 
       DO i = 2 , nreg - 1
!
          Call Re_KE (grid%reg_mat_hermite(i)%ham,                                      &
                      grid%reg_mat_hermite(i)%tr,                                       &
                      grid%reg_poly(i)%normalization,                                   &
                      npt(i),                                                           &
                      i,                                                                &
                      left_end_mat_el=grid%reg_mat(i-1)%tr(npt(i-1),npt(i-1)),          &
                      right_end_mat_el=grid%reg_mat(i+1)%tr(1,1) )
          write(iout,1) i
          IF (prn(4) == .true. ) THEN
           Call Print_Matrix(type_real_matrix,grid%reg_mat_hermite(i)%ham,npt(i),npt(i), &
                            title='renormalized hermite kinetic energy mat')
          END IF
       END DO
       i = nreg
!
!      Special case of last region.  Only the left hand elements needs to be changed.
!
!
          ALLOCATE(grid%reg_mat_hermite(i)%ham(npt(i),npt(i)))
          Call Re_KE (grid%reg_mat_hermite(i)%ham,                                     &
                      grid%reg_mat_hermite(i)%tr,                                      &
                      grid%reg_poly(i)%normalization,                                   &
                      npt(i),                                                           &
                      i,                                                                &
                      left_end_mat_el=grid%reg_mat_hermite(i-1)%tr(npt(i-1),npt(i-1)) )
          IF (prn(4) == .true. ) THEN
              write(iout,1) i
              Call Print_Matrix(type_real_matrix,grid%reg_mat_hermite(i)%ham,npt(i),npt(i), &
                                title='renormalized hermite kinetic energy mat')
          END IF
  END IF
  DO i = 1, nreg
     DEALLOCATE(grid%reg_mat_hermite(i)%tr)
  END DO
1 Format(/,10x,'Region = ',i4)
END SUBROUTINE Hermite_Renormalization
!***********************************************************************
!***********************************************************************
!deck Re_KE.f
!***begin prologue     Re_KE
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            1. calculate piecewise lobatto dvr functions and
!***                      their one-body matrices
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Re_KE

  SUBROUTINE Re_KE( ham, tr, norm, n , region, left_end_mat_el, right_end_mat_el )
  IMPLICIT NONE
  REAL(idp),  DIMENSION(:,:)      :: ham
  REAL(idp),  DIMENSION(:,:)      :: tr
  REAL(idp), DIMENSION (:)        :: norm
  REAL(idp), OPTIONAL             :: left_end_mat_el
  REAL(idp), OPTIONAL             :: right_end_mat_el
  INTEGER                         :: n
  INTEGER                         :: region
  INTEGER                         :: i
  INTEGER                         :: j
!
!
!
  IF ( .not.present(left_end_mat_el).and..not.present(right_end_mat_el) ) THEN
!
!      No bridge functions, so its straightforward.
!
       DO i = 1, n
          ham(i,1:i) = norm(i) * tr(i,1:i) * norm(1:i)
          ham(1:i,i) = ham(i,1:i)
       END DO
!
!       return
!  END IF
!
  ELSE IF ( .not.present(left_end_mat_el).and.present(right_end_mat_el) ) THEN
!
!      There is a bridge function at the right hand boundary that requires the
!      the matrix element involving the first function in the next sector.  This
!      is in the right_end_mat_el.
!      
       DO i = 1, n - 1
          ham(i,1:i) = norm(i) * tr(i,1:i) * norm(1:i)
          ham(1:i,i) = ham(i,1:i)
       END DO
       ham(n,1:n-1) = norm(n) * tr(n,1:n-1) * norm(1:n-1)
       ham(1:n-1,n) = ham(n,1:n-1)
       ham(n,n) = norm(n) * ( tr(n,n) + right_end_mat_el ) * norm(n)
!         
  ELSE IF ( present(left_end_mat_el).and..not.present(right_end_mat_el) ) THEN
!
!      This is the last sector.  Here we require information involving the bridge function
!      at the left hand boundary.
!
       ham(1,1) = norm(1) * ( tr(1,1) + left_end_mat_el ) * norm(1)       
       DO i = 2, n
          ham(i,1:i) = norm(i) * tr(i,1:i) * norm(1:i)
          ham(1:i,i) = ham(i,1:i)
       END DO
!
  ELSE IF ( present(left_end_mat_el).and.present(right_end_mat_el) ) THEN
!
!      This is the general case when there are contributions at both ends of the sector.
!      at the left hand boundary.
!       
       ham(1,1) = norm(1) * ( tr(1,1) + left_end_mat_el ) * norm(1)       
       DO i = 2, n - 1
          ham(i,1:i) = norm(i) * tr(i,1:i) * norm(1:i)
          ham(1:i,i) = ham(i,1:i)
       END DO
       ham(n,1:n-1) = norm(n) * tr(n,1:n-1) * norm(1:n-1)
       ham(1:n-1,n) = ham(n,1:n-1)
       ham(n,n) = norm(n) * ( tr(n,n) + right_end_mat_el ) * norm(n)
  END IF
END SUBROUTINE Re_KE
!***********************************************************************
!***********************************************************************
           END MODULE Renormalization_Module
!***********************************************************************
!***********************************************************************
