!***********************************************************************
! Matrix_Diagonalization
!**begin prologue     Matrix_Diagonalization
!**date written       090219   (yymmdd)
!**revision date               (yymmdd)
!**keywords           DVR, FEDVR
!**
!**author             schneider, b. i.(nsf)
!**source             DVR Library
!**purpose            Diagonalize the final DVR matrix
!***                  
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***                  
!***                  
!***end prologue      Matrix_Diagonalization
!***********************************************************************
!***********************************************************************
                           MODULE Matrix_Diagonalization
                           USE DVR_H_0_Module

!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck Diagonalize_Global_Matrix
!***begin prologue     Diagonalize_Global_Matrix
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Diagonalize the global matrix
!***                   
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Diagonalize_Global_Matrix
!
  SUBROUTINE Diagonalize_Global_Matrix(grid)
  IMPLICIT NONE
  TYPE(REAL_MATRIX)                              :: type_real_matrix
  TYPE(REAL_VECTOR)                              :: type_real_vector
  TYPE (coordinates)                             :: grid
  REAL(idp), ALLOCATABLE, DIMENSION(:,:)         :: test_matrix
  INTEGER                                        :: count
  INTEGER                                        :: i
  INTEGER                                        :: j
  INTEGER                                        :: lm
  INTEGER                                        :: info
  CHARACTER(LEN=3)                               :: itoc
  LOGICAL                                        :: logkey
!  
  write(iout,*)
  write(iout,*) '                              Diagonalizing the Full Matrix'
  write(iout,*)
  ALLOCATE(dvr_mat(0)%eigenvectors(1:physical_points,1:physical_points),                &
           dvr_mat(0)%eigenvalues(1:physical_points),                                   &
           dvr_mat(0)%work(3*physical_points),                                          &
           dvr_mat(0)%lower( physical_points*(physical_points+1)/2 ) )
  IF (grid%label(1:5) == 'theta' .or. grid%label(1:3) == 'eta') THEN
      write(iout,*) ' Will test for orthonormality of eigenvectors'
      ALLOCATE (test_matrix(1:physical_points,1:physical_points))
  END IF
  IF (diagonalize_nabla) THEN
      DO lm = 0, the_size
         count = 0
         DO i = 1, physical_points
            DO j = 1, i
               count = count + 1
               dvr_mat(0)%lower(count) = dvr_mat(lm)%tr(i,j)
            END DO
         END DO
         IF (.not.poisson) THEN
             DEALLOCATE(dvr_mat(lm)%tr)
         END IF
!
!        Used packed form of diagonalizer
!
         Call dspev('v','u',physical_points,dvr_mat(0)%lower,                                  &
                    dvr_mat(0)%eigenvalues,dvr_mat(0)%eigenvectors,                            &
                    physical_points,dvr_mat(0)%work,info)
         CALL Print_Matrix(type_real_vector,dvr_mat(0)%eigenvalues,                            &
                           title='eigenvalues of Nabla matrix angular quantum number = '//itoc(lm) )
         IF(prn(10)) THEN
            CALL Print_Matrix(type_real_matrix,dvr_mat(0)%eigenvectors,physical_points,        &
                              physical_points,title='eigenvectors of Nabla matrix')
         END IF
         IF (grid%label(1:5) == 'theta' .or. grid%label(1:3) == 'eta') THEN
             Call ebtc(test_matrix,dvr_mat(0)%eigenvectors,dvr_mat(0)%eigenvectors,            &
                       physical_points,physical_points,physical_points)
             CALL Print_Matrix(type_real_matrix,test_matrix,physical_points,physical_points,   & 
                               title='orthonormality matrix')
         END IF
         IF (grid%label(1:5) == 'theta' .or. grid%label(1:3) == 'eta') THEN
             Call ebtc(test_matrix,dvr_mat(0)%eigenvectors,dvr_mat(0)%eigenvectors,            &
                       physical_points,physical_points,physical_points)
             CALL Print_Matrix(type_real_matrix,test_matrix,physical_points,physical_points,   &
                               title='orthonormality matrix')
         END IF
      END DO
  END IF
  IF (diagonalize_hamiltonian) THEN
      DO lm = 0, the_size
         count = 0
         DO i = 1, physical_points
            DO j = 1, i
               count = count + 1
               dvr_mat(0)%lower(count) = dvr_mat(lm)%ham(i,j)
            END DO
         END DO
         DEALLOCATE(dvr_mat(lm)%ham)
!
!        Used packed form of diagonalizer
!
         Call dspev('v','u',physical_points,dvr_mat(0)%lower,                                  &
                    dvr_mat(0)%eigenvalues,dvr_mat(0)%eigenvectors,                            &
                    physical_points,dvr_mat(0)%work,info)
         CALL Print_Matrix(type_real_vector,dvr_mat(0)%eigenvalues,                            &
                           title='eigenvalues of Hamiltonian matrix l = '//itoc(lm) )
         IF(prn(10)) THEN
            CALL Print_Matrix(type_real_matrix,dvr_mat(0)%eigenvectors,physical_points,        &
                              physical_points,title='eigenvectors of Hamiltonian matrix l = '//itoc(lm)    )
         END IF
         IF (grid%label(1:5) == 'theta' .or. grid%label(1:3) == 'eta') THEN
            Call ebtc(test_matrix,dvr_mat(0)%eigenvectors,dvr_mat(0)%eigenvectors,            &
                      physical_points,physical_points,physical_points)
            CALL Print_Matrix(type_real_matrix,test_matrix,physical_points,physical_points,   &
                              title='orthonormality matrix')
         END IF
      END DO
  ELSE
      Call lnkerr('bad call in matrix diagonalization routine.')
  END IF
  IF (grid%label(1:5) == 'theta' .or. grid%label(1:3) == 'eta') THEN
      DEALLOCATE(test_matrix)
  END IF
  DEALLOCATE(dvr_mat(0)%eigenvectors,                                                     &
             dvr_mat(0)%eigenvalues,                                                      &
             dvr_mat(0)%work,                                                             &
             dvr_mat(0)%lower )
END SUBROUTINE Diagonalize_Global_Matrix
!***********************************************************************
!***********************************************************************
           END MODULE Matrix_Diagonalization
!***********************************************************************
!***********************************************************************
