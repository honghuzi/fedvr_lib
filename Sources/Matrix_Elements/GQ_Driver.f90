!deck GQ_Driver
!**begin prologue     GQ_Driver
!**date written       060711   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           
!**author             schneider, barry (nsf)
!**source
!**purpose            test GQ_Driver code
!**description        calls routines to input data, construct dvr matrices
!**                   
!**references
!**routines called        
!                       ----                    --------
!
!**end prologue       GQ_Driver
  PROGRAM GQ_Driver
  USE input_output
  USE  Generalized_Gaussian_Quadrature_Module
  USE DVR_Module
  IMPLICIT NONE
  INTEGER                        :: input
  INTEGER                        :: output
  COMMON /io/ input, output
!
  input = inp
  output = iout
!
!  Open the input and output files
!
  OPEN(input,file='GQ.inp',status='old')
  OPEN(output,file='GQ.out',status='unknown')
  write(iout,1)
  write(iout,2)
  Call Drive_Gauss 
  write(iout,1)
  CLOSE(input)
  CLOSE(output)
!
1 Format('           **************************************************************' &
         '****************')
2 Format(15x,'Calculation of FEDVR Points, Weights, Polynomials')
  stop
END PROGRAM GQ_Driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
