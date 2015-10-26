!deck Two_Electron_Integrals
!**begin prologue     Two_Electron_Integrals
!**date written       060711   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           
!**author             schneider, barry (nsf)
!**source
!**purpose            Two_Electron_Integrals
!**description        Computes two electron integrals using the Poisson equation
!**                   
!**references
!**routines called        
!                       ----                    --------
!
!**end prologue       
  PROGRAM Two_Electron_Integrals
  USE Two_Electron_FEDVR_Module
  USE input_output
  IMPLICIT NONE
!  INTEGER                          :: input
!  INTEGER                          :: output
  INTEGER                          :: intkey
  LOGICAL                          :: dollar
  CHARACTER (LEN=240)              :: chrkey
  CHARACTER (LEN=80)               :: home_directory
!  COMMON /io/ input, output
!
!  Get the input and output file numbers which appear in the 
!  Module input_output.f90 and put them into input and output which appears
!  in the ONLY common block in the code.  This is needed to pass into library
!  routines and for no other purpose.
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
!  OPEN(input,file=File_Directory(3)(1:len)//'/FEDVR_Two_Electron_Integrals.inp',status='old')
!  len=lenth(File_Directory(4))
!  OPEN(output,file=File_Directory(4)(1:len)//'/FEDVR_Two_Electron_Integrals.out',status='unknown')
!
  write(iout,1)
  write(iout,2)
  write(iout,1)
  Call Setup_2e 
  Call IOsys('rewind all on FEDVR_Two_Electron_Integral_File read-and-write',0,0,0,' ')
  Call IOsys('close Two_Electron_Integral_File',0,0,0,' ')
  CLOSE(input)
  CLOSE(output)
!
1 Format('           **************************************************************' &
         '****************')
2 Format(25x,'Calculation of FEDVR Two Electron Integrals')
END PROGRAM Two_Electron_Integrals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
