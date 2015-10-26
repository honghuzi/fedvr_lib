!deck FEDVR_Matrix_Elements
!**begin prologue     FEDVR_Matrix_Elements
!**date written       060711   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           
!**author             schneider, barry (nsf)
!**source
!**purpose            test FEDVR_Matrix_Elements code
!**description        calls routines to input data, construct dvr matrices
!**                   
!**references
!**routines called        
!                       ----                    --------
!
!**end prologue       FEDVR_Matrix_Elements
  PROGRAM FEDVR_Matrix_Elements
  USE DVR_Module
  USE Read_DVR_Module
  USE Matrix_Scale_and_Assemble
  USE Matrix_Diagonalization
  USE Poisson_Module
  IMPLICIT NONE
!  INTEGER                          :: input
!  INTEGER                          :: output
  INTEGER                          :: i
  INTEGER                          :: sets
  INTEGER                          :: intkey
  CHARACTER(LEN=80)                :: chrkey
  LOGICAL                          :: logkey
  REAL(idp)                        :: fpkey
  INTEGER                          :: no_sets
  LOGICAL                          :: test_key
  LOGICAL                          :: dollar
  CHARACTER(LEN=80)                :: home_directory
  CHARACTER(LEN=3)                 :: itoc

!  COMMON /io/ input, output
!
!  Get the input and output file numbers which appear in the 
!  Module input_output.f90 and put them into input and output which appears
!  in the ONLY common block in the code.  This is needed to pass into library
!  routines and for no other purpose.
!
  Call GetEnv('SOURCES',home_directory)
  len = lenth(home_directory)
  home_directory = home_directory(1:len)//'/fedvr_lib'
  nunits=10
  n_dir=6
  ALLOCATE(names(1:nunits),len_name(1:nunits),File_Directory(1:n_dir),len_dir(1:n_dir))
  Call Command_Line(home_directory=home_directory)
!  input = inp
!  output = iout
!
!  Open the input and output files
!
!  len=lenth(File_Directory(3))
!  OPEN(input,file=File_Directory(3)(1:len)//'/FEDVR_Matrix_Elements.inp',status='old')
!  len=lenth(File_Directory(4))
!  OPEN(output,file=File_Directory(4)(1:len)//'/FEDVR_Matrix_Elements.out',status='unknown')
!
  write(iout,1)
  write(iout,2)
  write(iout,1)
  IF ( dollar('$begin_calculation',card,cpass,inp) )THEN
     no_sets = intkey(card,'Number_of_Data_Sets',1,' ') ! number of data sets
     Write(iout,3) no_sets
  ELSE
     Call lnkerr('No initialization of data')
  END IF
Loop_over_Data_Sets :  &
  DO sets = 1, no_sets
     write(iout,*)
     Write(iout,4) sets
     write(iout,*)
     Call Set_General_Keywords(sets)  ! here we just read in some general keyword data
                                      ! see module Read_DVR_Module
!
     write(iout,*)
     Write(iout,5)
     write(iout,*)
     Loop_Over_Coordinates : &
     DO i = 1, spdim
        len=lenth(reg_grid(i)%label)
        write(iout,*)
        write(iout,6) i, reg_grid(i)%label(1:len)
        write(iout,*)
!------------------------------------------------------------------------------------------------------------------------
        Call Read_Data(reg_grid(i)) !  Read input data for 
                                    !  each type of coordinate
                                    !  or DVR. Contained in module Read_DVR_Module
!------------------------------------------------------------------------------------------------------------------------
        write(iout,*)
        write(iout,7)
        write(iout,*)
        Call Read_Potential_Data (reg_grid(i))    ! read data for one-body 
                                    ! potential. Contained in module Read_DVR_Module
!------------------------------------------------------------------------------------------------------------------------
        write(iout,*)
        write(iout,8)
        write(iout,*)
        Call Lobatto_Functions (reg_grid(i))  ! see module DVR_Polynomials_Module, calls
                                              ! either a lobatto or fourier routine
!------------------------------------------------------------------------------------------------------------------------
        write(iout,*)
        write(iout,9)
        write(iout,*)
        Call Coordinate_Factors(reg_grid(i))  ! also in DVR_Module.  routine computes
                                              ! the grid and a number of functions
                                              ! of the grid that are used in the
                                              ! construction of the matrix elements
!------------------------------------------------------------------------------------------------------------------------
        write(iout,*)
        write(iout,10)
        write(iout,*)
        Call KE_DVR_Matrices(reg_grid(i))    ! Here is where we calculate the
                                             ! the raw sector functions and the kinetic energy matrices.
                                             ! For some types of DVR basis sets there
                                             ! is only one sector.  The result is 
                                             ! dependent on both the type of DVR and
                                             ! the coordinate system.  The reg_mat arrays contain
                                             ! the results and these arrays have diffrent names.
!------------------------------------------------------------------------------------------------------------------------
        write(iout,*)
        write(iout,11)
        write(iout,*)
        Call PE_DVR_Matrix(reg_grid(i))       ! Calculates the sector potential energy matrix elements.
!------------------------------------------------------------------------------------------------------------------------
        write(iout,*)
        write(iout,12)
        write(iout,*)
        write(iout,13)
        write(iout,*)
        Call Final_KE_DVR_Matrices(reg_grid(i))
        write(iout,*)
        write(iout,14) 
        write(iout,*)
        Call H_0(reg_grid(i))
        write(iout,*)
        write(iout,15)
        write(iout,*)
        IF(diagonalize_hamiltonian .or. diagonalize_nabla .or. poisson) THEN
           the_size=max(0,l_max,m_max) 
           ALLOCATE(dvr_mat(0:the_size))            
        END IF
        Call Form_Matrix(reg_grid(i))
        IF(diagonalize_hamiltonian ==.true. .or. diagonalize_nabla == .true.) THEN
           write(iout,*)
           write(iout,16)
           write(iout,*)
           Call Diagonalize_Global_Matrix(reg_grid(i))
        END IF
        IF (poisson) THEN
            write(iout,*)
            write(iout,17)
            write(iout,*)
            IF ( dollar('$poisson_'//reg_grid(i)%label(1:len),card,cpass,inp) ) THEN
                 number_of_right_hand_sides=intkey(card,'number_of_right_hand_sides',            &
                                                   physical_points,' ')
                 type_inhomo=chrkey(card,'inhomogeneity','linear',' ') 
                 last_value = fpkey(card,'last_value',one,' ')
            END IF
            Call Poisson_Equation(reg_grid(i))
        END IF
        write(iout,*)
        write(iout,18)
        write(iout,*)
        IF(diagonalize_hamiltonian .or. diagonalize_nabla .or. poisson) THEN
            DEALLOCATE(dvr_mat)           
        END IF
        Call IOsys('rewind all on '//FEDVR_File//' read-and-write', &
                    0,0,0,' ')
        Call IOsys('close '//FEDVR_File,0,0,0,' ')
        write(iout,*)
        write(iout,19) sets, i
     END DO Loop_Over_Coordinates
        DEALLOCATE(reg_grid)
  END DO Loop_over_Data_Sets
  CLOSE(input)
  CLOSE(output)
!
1 Format('           **************************************************************' &
         '****************')
2 Format(15x,'Calculation of FEDVR Points, Weights, Polynomials and Matrix Elements')
3 Format(30x,'** Number of Data Sets ** = ',i3)
4 Format(30x,'***    Processing Data Set *** = ',i3)
5 Format(30x,'***    Beginning Loop over Coordinates    ***')
6 Format(30x,'***    Reading data for Coordinate = ',i3,' Coordinate label = ',a8,' ***')
7 Format(30x,'***    Reading Potential Data    ***')
8 Format(30x,'***    Compute the required FEDVR Basis set information    ***')
9 Format(30x,'***    Compute some additional grid information and a number of functions',      &
        /30x '       that depend on the grid and are used later    ***')
10 Format(30x,'***    Compute the primitive raw DVR sector matrices for the kinetic energy    ***')
11 Format(30x,'***    Compute the diagonal potential matrix elements    ***')
12 Format(30x,'***    Compute the final FEDVR sector matrix elements: This is done in few steps')
13 Format(30x,'***    First step, Renormalize the basic FEDVR sector matrices    ***')
14 Format(30x,'***    Second step, Form the final FEDVR sector matrices with angular momenta',  &
         /30x,'       and potential added    ***')
15 Format(30x,'***    Third step, Take the final FEDVR sector matrices,',                      &
         /30x,'       transform to standard form and assemble the global matrix ',              &
         /30x,'       for further processing    ***')
16 Format(30x,'***    Fourth step. Diagonalize the global matrix if requested    ***')
17 Format(30x,'***    Fifth step. Solve the Poisson equation if requested    ***')
18 Format(30x,'***    Clean up.  We are done with this coordinate    ***')
19 Format(30x,'***    Computation finished for set = 'i3,' and dimension = ',i2)
  stop
END PROGRAM FEDVR_Matrix_Elements
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
