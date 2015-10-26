!***********************************************************************
! Matrix_Scale_and_Assemble
!**begin prologue     Matrix_Scale_and_Assemble
!**date written       090219   (yymmdd)
!**revision date               (yymmdd)
!**keywords           DVR, FEDVR
!**
!**author             schneider, b. i.(nsf)
!**source             DVR Library
!**purpose            Assemble the final DVR matrix
!***                  in a FEDVR basis
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***                  
!***                  
!***end prologue      Matrix_Scale_and_Assemble
!***********************************************************************
!***********************************************************************
                           MODULE Matrix_Scale_and_Assemble
                           USE DVR_H_0_Module
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck Form_Matrix
!***begin prologue     Form_Matrix
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Assemble the full Hamiltonian matrix elements from the sector matrices
!***description        The routine assembles the full matrix from the sector matrices after transforming 
!***                   them to standard form.  The final result is the global matrix ready for further use.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Form_Matrix
!
  SUBROUTINE Form_Matrix(grid)
  IMPLICIT NONE
  TYPE (coordinates)                             :: grid
!
!         Some coordinate systems and coordinates have a metric which needs to scale the Hamiltonian
!
  write(iout,*)
  write(iout,*) '                              Forming the Final Full Matrix'
  write(iout,*)
  len=lenth(grid%label)
!
!       For the pure radial spherical case and the pure rho cylindrical case, we scale the
!       appropriate Laplacian so that the metric is accounted for.
!
  Call Write_Unscaled_Regional_Matrices(grid)
  IF ( keyword == 'spherical' ) THEN
       IF ( grid%label(1:len) == 'r') THEN
            Call Transform_Matrix_to_Standard_Form(grid)          
       END IF
  END IF
  IF ( keyword == 'cylindrical' ) THEN
       IF ( grid%label(1:len) == 'rho') THEN
            Call Transform_Matrix_to_Standard_Form(grid)          
       END IF
  END IF
  IF (prn(4) == .true. ) THEN
      Call Print_Scaled_Matrices(grid)
  END IF
  Call Write_Scaled_Regional_Matrices(grid)
!  Call Read_Scaled_Matrices(grid)
  Call Form_Global_Matrix(grid)
  Call Write_Global_Matrices(grid)
  IF (prn(4) == .true. ) THEN
      Call Print_Global_Matrices(grid)
  END IF
!
!
END SUBROUTINE Form_Matrix
!***********************************************************************
!***********************************************************************
!deck Transform_Matrix_to_Standard_Form
!***begin prologue     Transform_Matrix_to_Standard_Form
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            In certain coordinate systems there is a simple metric
!                      due to the structure of the kinetic energy matrix
!                      which turn a standard eigenvalue problem to what looks
!                      like a general eigenvalue problem.  This is trivially
!                      removed by a simple diagonal pre and post multiplication.
!                      This unitary transformation also makes the definition of
!                      the vectors different and this needs to be accounted for later.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       
!
  SUBROUTINE Transform_Matrix_to_Standard_Form(grid)
  IMPLICIT NONE
  TYPE (coordinates)             :: grid
  INTEGER                        :: lm
  INTEGER                        :: i
  INTEGER                        :: j
!
  IF (form_hamiltonian) THEN
      write(iout,*)
      write(iout,*) 'Forming Scaled Hamiltonian'
      DO lm = 0, lm_max
         DO i = 1, nreg
            DO j = 1, npt(i)
!
!              Pre and post multiplication
!
               grid%reg_type_op(i,lm)%ham(j,1:j) =                                       &
                                                   grid%reg_pt_wt(i)%inv_sqrt_qr_fac(j)  &
                                                           *                             &
                                                   grid%reg_type_op(i,lm)%ham(j,1:j)     &
                                                           *                             &
                                                   grid%reg_pt_wt(i)%inv_sqrt_qr_fac(1:j)
               grid%reg_type_op(i,lm)%ham(1:j,j) = grid%reg_type_op(i,lm)%ham(j,1:j) 
            END DO
         END DO
      END DO
  END IF
!
  IF (form_nabla) THEN
      write(iout,*)
      write(iout,*) 'Forming Scaled Nabla'
      DO lm = 0, lm_max
         DO i = 1, nreg
            DO j = 1, npt(i)
!
!              Pre and post multiplication
!
               grid%reg_type_op(i,lm)%tr(j,1:j) =                                        &
                                                   grid%reg_pt_wt(i)%inv_sqrt_qr_fac(j)  &
                                                           *                             &
                                                   grid%reg_type_op(i,lm)%tr(j,1:j)      &
                                                           *                             &
                                                   grid%reg_pt_wt(i)%inv_sqrt_qr_fac(1:j)
               grid%reg_type_op(i,lm)%tr(1:j,j) = grid%reg_type_op(i,lm)%tr(j,1:j) 
            END DO
         END DO
      END DO
  END IF
!
END SUBROUTINE Transform_Matrix_to_Standard_Form
!***********************************************************************
!***********************************************************************
!deck Print_Scaled_Matrices
!***begin prologue     Print_Scaled_Matrices
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Print final scaled regional matrices to disk
!***references
!***routines called    iosys, util and mdutil
!***end prologue       
!
  SUBROUTINE Print_Scaled_Matrices(grid)
  IMPLICIT NONE
  TYPE (coordinates)             :: grid
  INTEGER                        :: lm
  INTEGER                        :: i
  CHARACTER (LEN=3)              :: itoc
!
  IF (form_hamiltonian) THEN
      DO lm = 0, lm_max
         DO i = 1, nreg
            Call Print_Matrix(type_real_matrix,grid%reg_type_op(i,lm)%ham,npt(i),npt(i),   &
                              title = 'Scaled Hamiltonian matrix LM = '//itoc(lm)//' Region '//itoc(i) )
         END DO
      END DO
  END IF
  IF (form_nabla) THEN
      DO lm = 0, lm_max
         DO i = 1, nreg
            Call Print_Matrix(type_real_matrix,grid%reg_type_op(i,lm)%tr,npt(i),npt(i),    &
                              title = 'Scaled Nabla matrix LM = '//itoc(lm)//' Region '//itoc(i) )
         END DO
      END DO
  END IF
!
END SUBROUTINE Print_Scaled_Matrices
!***********************************************************************
!***********************************************************************
!deck Write_Unscaled_Regional_Matrices                                                                                        
!***begin prologue     Write_Unscaled_Regional_Matrices                                                                       
!***date written       960718   (yymmdd)                                                                                      
!***revision date               (yymmdd)                                                                                      
!***keywords                                                                                                                  
!***author             schneider, b. i.(nsf)                                                                                  
!***source                                                                                                                    
!***purpose            Write unscaled regional matrices to disk                                                               
!***references                                                                                                                
!***routines called    iosys, util and mdutil                                                                                 
!***end prologue                                                                                                              
!                                                                                                                             
  SUBROUTINE Write_Unscaled_Regional_Matrices(grid)
  IMPLICIT NONE
  TYPE (coordinates)             :: grid
  INTEGER                        :: lm
  INTEGER                        :: i
  INTEGER                        :: len
  INTEGER                        :: lenth
  CHARACTER (LEN=3)              :: itoc
!                                                                                                                             
  IF (form_hamiltonian) THEN
      write(iout,*) '         Writing Unscaled Regional Hamiltonian Matrices to disk'
      DO lm = 0, lm_max
         title = 'Unscaled H-mat l/m = '//itoc(lm)
         len=lenth(title)
         write(iout,*) title(1:len)
         DO i = 1, nreg
            title = 'Unscaled H-mat l/m = '//itoc(lm)//' reg = '//itoc(i)
            len=lenth(title)
            Call iosys ('write real "'//title(1:len)//' " to '                     &
                        //FEDVR_File,npt(i)*npt(i),grid%reg_type_op(i,lm)%ham,0,' ')
         END DO
      END DO
  END IF
  IF (form_nabla) THEN
      write(iout,*) '         Writing Unscaled Nabla Matrices to disk'
      DO lm = 0, lm_max
         title = 'Unscaled T-mat l/m = '//itoc(lm)
         len=lenth(title)
         write(iout,*) title(1:len)
         DO i = 1, nreg
            title = 'Unscaled T-mat l/m = '//itoc(lm)//' reg = '//itoc(i)
            len=lenth(title)
            Call iosys ('write real "'//title(1:len)//' " to '                      &
                        //FEDVR_File,npt(i)*npt(i),grid%reg_type_op(i,lm)%tr,0,' ')
         END DO
      END DO
  END IF
! 
  END SUBROUTINE Write_Unscaled_Regional_Matrices
!***********************************************************************
!***********************************************************************
!deck Write_Scaled_Regional_Matrices                                                                                        
!***begin prologue     Write_Scaled_Regional_Matrices                                                                       
!***date written       960718   (yymmdd)                                                                                      
!***revision date               (yymmdd)                                                                                      
!***keywords                                                                                                                  
!***author             schneider, b. i.(nsf)                                                                                  
!***source                                                                                                                    
!***purpose            Write scaled regional matrices to disk                                                               
!***references                                                                                                                
!***routines called    iosys, util and mdutil                                                                                 
!***end prologue                                                                                                              
!                                                
  SUBROUTINE Write_Scaled_Regional_Matrices(grid)
  IMPLICIT NONE
  TYPE (coordinates)             :: grid
  INTEGER                        :: lm
  INTEGER                        :: i
  INTEGER                        :: len
  INTEGER                        :: lenth
  CHARACTER (LEN=3)              :: itoc
!                                                                                                                             
  IF (form_hamiltonian) THEN
      write(iout,*) '         Writing Scaled Regional Hamiltonian Matrices to disk'
      DO lm = 0, lm_max
         title = 'Scaled H-mat l/m = '//itoc(lm)
         len=lenth(title)
         write(iout,*) title(1:len)
         DO i = 1, nreg
            title = 'Scaled H-mat l/m = '//itoc(lm)//' reg = '//itoc(i)
            len=lenth(title)
            Call iosys ('write real "'//title(1:len)//' " to '                     &
                        //FEDVR_File,npt(i)*npt(i),grid%reg_type_op(i,lm)%ham,0,' ')
         END DO
      END DO
  END IF
  IF (form_nabla) THEN
      write(iout,*) '         Writing Scaled Nabla Matrices to disk'
      DO lm = 0, lm_max
         title = 'Scaled T-mat l/m = '//itoc(lm)
         len=lenth(title)
         write(iout,*) title(1:len)
         DO i = 1, nreg
            title = 'Scaled T-mat l/m = '//itoc(lm)//' reg = '//itoc(i)
            len=lenth(title)
            Call iosys ('write real "'//title(1:len)//' " to '                      &
                        //FEDVR_File,npt(i)*npt(i),grid%reg_type_op(i,lm)%tr,0,' ')
         END DO
      END DO
  END IF
! 
  END SUBROUTINE Write_Scaled_Regional_Matrices
!***********************************************************************
!***********************************************************************
!deck Read_Scaled_Matrices
!***begin prologue     Read_Scaled_Matrices
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Read final scaled regional matrices to disk
!***references
!***routines called    iosys, util and mdutil
!***end prologue       
!
  SUBROUTINE Read_Scaled_Matrices(grid)
  IMPLICIT NONE
  TYPE (coordinates)             :: grid
  INTEGER                        :: lm
  INTEGER                        :: i
  INTEGER                        :: len
  INTEGER                        :: lenth
  CHARACTER (LEN=3)              :: itoc
!
  write(iout,*)
  IF (read_hamiltonian) THEN
      DO lm = 0, lm_max
         title = 'scaled H-mat l/m = '//itoc(lm)
         len=lenth(title)
         write(iout,*) '               Reading '//title(1:len)//' to disk'
         DO i = 1, nreg
            title = 'scaled H-mat l/m = '//itoc(lm)//' reg = '//itoc(i)
            len=lenth(title)
            Call iosys ('read real "'//title(1:len)//' " from '                     &
                        //FEDVR_File,npt(i)*npt(i),grid%reg_type_op(i,lm)%ham,0,' ')
         END DO
      END DO
  END IF
  IF (read_nabla) THEN
      DO lm = 0, lm_max
         title = 'scaled T-mat l/m = '//itoc(lm)
         len=lenth(title)
         write(iout,*) '               Reading '//title(1:len)//' to disk'
         DO i = 1, nreg
            title = 'scaled T-mat l/m = '//itoc(lm)//' reg = '//itoc(i)
            len=lenth(title)
            Call iosys ('read real "'//title(1:len)//' " from '                      &
                        //FEDVR_File,npt(i)*npt(i),grid%reg_type_op(i,lm)%tr,0,' ')
         END DO
      END DO
  END IF
!
END SUBROUTINE Read_Scaled_Matrices
!***********************************************************************
!***********************************************************************
!deck Form_Global_Matrix.f
!***begin prologue     Form_Global_Matrix
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
  SUBROUTINE Form_Global_Matrix(grid)
  IMPLICIT NONE
  TYPE (coordinates)                      :: grid
  INTEGER                                 :: i
  INTEGER                                 :: lm
  INTEGER                                 :: first
  INTEGER                                 :: last
!  
  IF (form_hamiltonian) THEN
      DO lm = 0, lm_max
         ALLOCATE(dvr_mat(lm)%ham(1:physical_points,1:physical_points))
         dvr_mat(lm)%ham(:,:) = zero
         first = 1
         DO i = 1, nreg
            last = first + npt(i) - 1
            dvr_mat(lm)%ham(first:last,first:last) = grid%reg_type_op(i,lm)%ham(1:npt(i),1:npt(i))         
            first = last
            DEALLOCATE(grid%reg_type_op(i,lm)%ham)
         END DO    
      END DO
  END IF
  IF (form_nabla) THEN
      DO lm = 0, lm_max
         ALLOCATE(dvr_mat(lm)%tr(1:physical_points,1:physical_points))
         dvr_mat(lm)%tr(:,:) = zero
         first = 1
         DO i = 1, nreg
            last = first + npt(i) - 1
            dvr_mat(lm)%tr(first:last,first:last) = grid%reg_type_op(i,lm)%tr(1:npt(i),1:npt(i))         
            first = last
            DEALLOCATE(grid%reg_type_op(i,lm)%tr)
         END DO    
      END DO
  END IF
  DEALLOCATE(grid%reg_type_op)
END SUBROUTINE Form_Global_Matrix
!***********************************************************************
!***********************************************************************
!Write_Global_Matrices
!***begin prologue     Write_Global_Matrices
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Write_Global_Matrices
!***references
!***routines called    iosys, util and mdutil
!***end prologue       
!
  SUBROUTINE Write_Global_Matrices(grid)
  IMPLICIT NONE
  TYPE (coordinates)                      :: grid
  INTEGER                                 :: lm
  INTEGER                                 :: len
  INTEGER                                 :: lenth
  CHARACTER(LEN=3)                        :: itoc
!  
  IF (write_hamiltonian) THEN
      DO lm = 0, lm_max
         title = 'global H-mat l/m = '//itoc(lm)
         len=lenth(title)
         write(iout,*) '               Writing '//title(1:len)//' to disk'
         Call IOsys('write real "'//title//' " to '//FEDVR_File,       &
                     physical_points*physical_points,dvr_mat(lm)%ham,0,' ')
      END DO
  END IF
  IF (write_nabla) THEN
      DO lm = 0, lm_max
         title = 'global T-mat l/m = '//itoc(lm)
         len=lenth(title)
         write(iout,*) '               Writing '//title(1:len)//' to disk'
         Call IOsys('write real "'//title//' " to '//FEDVR_File,        &
                     physical_points*physical_points,dvr_mat(lm)%tr,0,' ')

      END DO
  END IF
END SUBROUTINE Write_Global_Matrices
!***********************************************************************
!***********************************************************************
! Print_Global_Matrices
!***begin prologue      Print_Global_Matrices
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Matrix_Output
!***references
!***routines called    iosys, util and mdutil
!***end prologue       
!
  SUBROUTINE  Print_Global_Matrices(grid)
  IMPLICIT NONE
  TYPE (coordinates)                      :: grid
  INTEGER                                 :: i
  INTEGER                                 :: lm
  CHARACTER(LEN=3)                        :: itoc
!  
  IF (form_hamiltonian) THEN
      DO lm = 0, lm_max
         Call Print_Matrix(type_real_matrix,dvr_mat(lm)%ham,physical_points,physical_points,           &
                           title = 'global H-mat l/m = '//itoc(lm) )
      END DO
  END IF
  IF (form_nabla) THEN
      DO lm = 0, lm_max
         Call Print_Matrix(type_real_matrix,dvr_mat(lm)%tr,physical_points,physical_points,           &
                           title = 'global T-mat l/m = '//itoc(lm) )
      END DO
  END IF
END SUBROUTINE  Print_Global_Matrices
!***********************************************************************
!***********************************************************************
           END MODULE Matrix_Scale_and_Assemble
!***********************************************************************
!***********************************************************************
