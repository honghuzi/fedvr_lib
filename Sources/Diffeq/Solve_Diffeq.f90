!***********************************************************************
! Solve_Diffeq
!**begin prologue     Solve_Diffeq
!**date written       090219   (yymmdd)
!**revision date               (yymmdd)
!**keywords           DVR, FEDVR
!**
!**author             schneider, b. i.(nsf)
!**source             DVR Library
!**purpose            Driver to solve second order differentail eautions using FEDVR
!***                  
!***description       
!***                  
!***                  
!***                  
!***                  
!***references
!***modules needed    
!***comments          
!***                  
!***                  
!***                  
!***                  
!***end prologue      Solve_Diffeq
  PROGRAM Solve_Diffeq
  USE Read_DVR_Module
  USE Matrix_Scale_and_Assemble
  USE Matrix_Diagonalization
  USE Poisson_Module
  USE Diffeq_Module
!  USE Two_Electron_FEDVR_Module
  IMPLICIT NONE
  CHARACTER(LEN=80)                                :: home_directory
  REAL(idp)                                        :: fpkey
  LOGICAL                                          :: dollar
  LOGICAL                                          :: logkey
  INTEGER                                          :: intkey  
  INTEGER                                          :: i
  CHARACTER(LEN=80)                                :: chrkey
  CHARACTER(LEN=8)                                 :: itoc
!
  Call GetEnv('SOURCES',home_directory)
  len = lenth(home_directory)
  home_directory = home_directory(1:len)//'/fedvr_lib'
  nunits=10
  n_dir=6
  ALLOCATE(names(1:nunits),len_name(1:nunits),File_Directory(1:n_dir),len_dir(1:n_dir))
  Call Command_Line(home_directory=home_directory)
!
  write(iout,1)
  write(iout,2)
  write(iout,1)
!
!
  Call Read_Diffeq
  DO i = 1, ndim  
     write(iout,3) reg_grid(i)%label
     Call FEDVR_File_Open(reg_grid(i),new_old)
     Call Read_FEDVR_Information(reg_grid(i))
     IF ( dollar('$poisson_'//reg_grid(i)%label,card,cpass,inp) ) THEN
          number_of_right_hand_sides=intkey(card,'number_of_right_hand_sides',            &
                                            physical_points,' ')
          type_inhomo=chrkey(card,'inhomogeneity','linear',' ') 
          last_value = fpkey(card,'last_value',one,' ')
     END IF
     Call PDE(reg_grid(i))
  END DO
1 Format('           **************************************************************' &
         '****************')
2 Format(15x,'Solve Second Order Differential Equations using an FEDVR Points')
3 Format(/,15x,'Processing Coordinate = ',a8)
 stop
END PROGRAM Solve_Diffeq
!***********************************************************************
!***********************************************************************
