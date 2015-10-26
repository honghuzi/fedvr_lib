!***********************************************************************
! FEDVR_Grid_Functions
!**begin prologue     FEDVR_Grid_Functions
!**date written       090119   (yymmdd)
!**revision date               (yymmdd)
!**keywords           FEDVR, grid, polynomials
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Calculate the FEDVR points, weights, polynomials, first and second derivatives
!***                  of the polynomials and the sector renormaliztion factors.
!***                  
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***end prologue      FEDVR_Grid_Functions
!***********************************************************************
!***********************************************************************
                           MODULE FEDVR_Grid_Functions
                           USE Data_Module
                           USE FEDVR_Shared
                           USE FEDVR_Derived_Types
                           USE Renormalization_Module
!***********************************************************************
!***********************************************************************
!***********************************************************************                          
!                            Explicit Interfaces

                            INTERFACE Coordinate_Functions
                       MODULE PROCEDURE Cartesian_Grid,                    &
                                        Spherical_Grid,                    &
                                        Cylindrical_Grid,                  &
                                        Spheroidal_Grid
                            END INTERFACE Coordinate_Functions
!**********************************************************************
!**********************************************************************   

                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck Grid_Functions
!***begin prologue     Grid_Functions
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Compute sector grid, weights and coordinate factors
!***                   needed to construct the KE matrix elements.

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE Grid_Functions(grid,reg)
  IMPLICIT NONE         
  TYPE(coordinates)                          :: grid
  TYPE(regional)                             :: reg
!
!
  IF (grid%type_coordinates == 'cartesian') THEN
      Call Coordinate_Finctions(grid%xyz,reg)
  ELSE IF(grid%type_coordinates == 'spherical') THEN
      Call Coordinate_Finctions(grid%r_theta,reg)
  ELSE IF(grid%type_coordinates == 'cylindrical') THEN
      Call Coordinate_Finctions(grid%rho_z,reg)
  ELSE IF(grid%type_coordinates == 'spheroidal') THEN
      Call Coordinate_Finctions(grid%xi_eta,reg)
  END IF
  Call Global_Points_Weights(grid,reg)
  Call Normalization(reg)
END SUBROUTINE Grid_Function
!***********************************************************************
!***********************************************************************
!deck Cartesian_Grid
!***begin prologue     Cartesian_Grid
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Compute sector grid, weights and coordinate factors
!***                   needed to construct the KE matrix elements.

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE Cartesian_Grid(xyz,reg)
  IMPLICIT NONE
  TYPE(cartesian)               :: xyz
  TYPE(regional)                :: reg
!
!
  DO i = 1, nreg
     ALLOCATE( reg(i)%q(1:reg(i)%n), reg(i)%wt(1:reg(i)%n),                                        &
               reg(i)%p(1:reg(i)%n,1:reg(i)%n), reg(i)%dp(1:reg(i)%n,1:reg(i)%n),                  &
               reg(i)%ddp(1:reg(i)%n,1:reg(i)%n), reg(i)%q_fac(1:reg(i)%n),                        &
               reg(i)%inv_q_fac(1:reg(i)%n), reg(i)%inv_sqrt_q_fac(1:reg(i)%n) )
!
     Call gauss(reg(i)%q, reg(i)%wt, reg%(i)%edge, reg(i)%p, reg(i)%dp, reg(i)%ddp,                &
                reg(i)%type_quadrature, reg(i)%fixed_point, reg(i)%n)
!
!         For cartesian coordinates they are all one.
!
     reg(i)%q_fac(1:reg(i)%n) = one
     reg(i)%inv_reg(i)%q_fac(1:reg(i)%n) = one 
     reg(i)%inv_sqrt_reg(i)%q_fac(1:reg(i)%n) = one 
!
!
  END DO
  IF (prn == .true.) THEN
      DO i = 1, nreg
         call prntfmn('q-'//itoc(i),reg(i)%q,reg(i)%n,1,reg(i)%n,1,iout,'e')
         call prntfmn('wt-'//itoc(i),reg(i)%wt,reg(i)%n,1,reg(i)%n,1,iout,'e')
         call prntfmn('qr_factor-'//itoc(i),reg(i)%q_fac,reg(i)%n,1,reg(i)%n,1,iout,'e')
         call prntfmn('inverse_qr_factor-'//itoc(i),reg(i)%inv_q_fac,reg(i)%n,1,reg(i)%n,1,iout,'e')
         call prntfmn('inverse_sqrt_qr_factor-'//itoc(i),reg(i)%inv_sqrt_q_fac,                   &
                                                         reg(i)%n,1,reg(i)%n,1,iout,'e'))
         Call prntfmn('unnormalized polynomials-'//itoc(i),reg(i)%pr,reg(i)%n,reg(i)%n,           &
                                                           reg(i)%n,reg(i)%n,iout,'e'
         Call prntfmn('first derivative of unnormalized polynomials-'//itoc(i),reg(i)%dpr,        &
                                                                               reg(i)%n,          &
                                                                               reg(i)%n,          &
                                                                               reg(i)%n,          &
                                                                               reg(i)%n,iout,'e')
         Call prntfmn('second derivative of unnormalized polynomials-'//itoc(i),reg(i)%ddpr,      &
                                                                               reg(i)%n,          &
                                                                               reg(i)%n,          &
                                                                               reg(i)%n,          &
                                                                               reg(i)%n,iout,'e')
     END DO
  END IF  
!
END SUBROUTINE Cartesian_Grid

!***********************************************************************
!***********************************************************************
!Deck Spherical_Grid
!***begin prologue     Spherical_Grid
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Compute sector grid, weights and coordinate factors
!***                   needed to construct the KE matrix elements.

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE Spherical_Grid(r_theta,reg)
  IMPLICIT NONE
  TYPE(spherical)              :: r_theta
  TYPE(regional)               :: reg
  CHARACTER(LEN=3)             :: itoc
!
!
  IF (r_theta%axis == 'r') THEN
      IF ( grid%drop(1) == .true.) THEN
           DO i = 1, nreg
              ALLOCATE( reg(i)%q(1:reg(i)%n), reg(i)%wt(1:reg(i)%n),                                        &
                        reg(i)%p(1:reg(i)%n,1:reg(i)%n), reg(i)%dp(1:reg(i)%n,1:reg(i)%n),                  &
                        reg(i)%ddp(1:reg(i)%n,1:reg(i)%n), reg(i)%q_fac(1:reg(i)%n),                        &
                        reg(i)%inv_q_fac(1:reg(i)%n), reg(i)%inv_sqrt_q_fac(1:reg(i)%n) )
!
              Call gauss(reg(i)%q, reg(i)%wt, reg%(i)%edge, reg(i)%p, reg(i)%dp, reg(i)%ddp,                &
                         reg(i)%type_quadrature, reg(i)%fixed_point, reg(i)%n)
!
                         reg(i)%q_fac(1:reg(i)%n) = one
                         reg(i)%inv_reg(i)%q_fac(1:reg(i)%n) = one 
                         reg(i)%inv_sqrt_reg(i)%q_fac(1:reg(i)%n) = one 
!
           END DO
      ELSE
!
           DO i = 1, nreg
!
              ALLOCATE( reg(i)%q(1:reg(i)%n), reg(i)%wt(1:reg(i)%n),                                        &
                        reg(i)%p(1:reg(i)%n,1:reg(i)%n), reg(i)%dp(1:reg(i)%n,1:reg(i)%n),                  &
                        reg(i)%ddp(1:reg(i)%n,1:reg(i)%n), reg(i)%q_fac(1:reg(i)%n),                        &
                        reg(i)%inv_q_fac(1:reg(i)%n), reg(i)%inv_sqrt_q_fac(1:reg(i)%n) )
!          
              Call gauss(reg(i)%q, reg(i)%wt, reg%(i)%edge, reg(i)%p, reg(i)%dp, reg(i)%ddp,                &
                         reg(i)%type_quadrature, reg(i)%fixed_point, reg(i)%n)
                         reg(i)%q_fac(1:reg(i)%n) = reg(i)%q(1:reg(i)%n) * reg(i)%q(1:reg(i)%n) 
                         reg(i)%inv_q_fac(1:reg(i)%n)= one  / reg(i)%q_fac(1:reg(i)%n)
!
                         reg(i)%q_fac(1:reg(i)%n) = reg(i)%q(1:reg(i)%n) * reg(i)%q(1:reg(i)%n) 
                         reg(i)%inv_q_fac(1:reg(i)%n)= one  / reg(i)%q_fac(1:reg(i)%n)
                         reg(i)%inv_sqrt_q_fac(1:reg(i)%n) = one / reg(i)%q(1:reg(i)%n) 
!
           END DO
      END IF
  ELSE IF (r_theta%axis == 'theta') THEN
       DO i = 1, nreg
!
          ALLOCATE( reg(i)%q(1:reg(i)%n), reg(i)%wt(1:reg(i)%n),                                        &
                    reg(i)%p(1:reg(i)%n,1:reg(i)%n), reg(i)%dp(1:reg(i)%n,1:reg(i)%n),                  &
                    reg(i)%ddp(1:reg(i)%n,1:reg(i)%n), reg(i)%q_fac(1:reg(i)%n),                        &
                    reg(i)%inv_q_fac(1:reg(i)%n), reg(i)%inv_sqrt_q_fac(1:reg(i)%n) )
!
          Call gauss(reg(i)%q, reg(i)%wt, reg%(i)%edge, reg(i)%p, reg(i)%dp, reg(i)%ddp,                &
                     reg(i)%type_quadrature, reg(i)%fixed_point, reg(i)%n)
!
                    reg(i)%q_fac(1:reg(i)%n) = one - reg(i)%q(1:reg(i)%n) * reg(i)%q(1:reg(i)%n)
                    reg(i)%inv_q_fac(1:reg(i)%n) = one / reg(i)%q_fac(1:reg(i)%n)
                    reg(i)%inv_sqrt_q_fac(1:reg(i)%n) = Sqrt ( reg(i)%inv_q_fac(1:reg(i)%n) )
!
       END DO
  ELSE
       Call lnkerr('error in axis variable')
!
  END IF
!
  IF (prn == .true.) THEN
      DO i = 1, nreg
         call prntfmn('q-'//itoc(i),reg(i)%q,reg(i)%n,1,reg(i)%n,1,iout,'e')
         call prntfmn('wt-'//itoc(i),reg(i)%wt,reg(i)%n,1,reg(i)%n,1,iout,'e')
         call prntfmn('qr_factor-'//itoc(i),reg(i)%q_fac,reg(i)%n,1,reg(i)%n,1,iout,'e')
         call prntfmn('inverse_qr_factor-'//itoc(i),reg(i)%inv_q_fac,reg(i)%n,1,reg(i)%n,1,iout,'e')
         call prntfmn('inverse_sqrt_qr_factor-'//itoc(i),reg(i)%inv_sqrt_q_fac,                   &
                                                         reg(i)%n,1,reg(i)%n,1,iout,'e'))
         Call prntfmn('unnormalized polynomials-'//itoc(i),reg(i)%pr,reg(i)%n,reg(i)%n,           &
                                                           reg(i)%n,reg(i)%n,iout,'e'
         Call prntfmn('first derivative of unnormalized polynomials-'//itoc(i),reg(i)%dpr,        &
                                                                               reg(i)%n,          &
                                                                               reg(i)%n,          &
                                                                               reg(i)%n,          &
                                                                               reg(i)%n,iout,'e')
         Call prntfmn('second derivative of unnormalized polynomials-'//itoc(i),reg(i)%ddpr,      &
                                                                               reg(i)%n,          &
                                                                               reg(i)%n,          &
                                                                               reg(i)%n,          &
                                                                               reg(i)%n,iout,'e')
  END IF  
!
END SUBROUTINE Spherical_Grid
!***********************************************************************
!***********************************************************************
!deck Cylindrical_Grid
!***begin prologue     Cylindrical_Grid
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Compute sector grid, weights and coordinate factors
!***                   needed to construct the KE matrix elements.

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE Cylindrical_Grid(rho_z,reg)
  IMPLICIT NONE
  TYPE(cylindical)              :: rho_z
  TYPE(regional)                :: reg
!
  IF (rho_z%axis == 'z') THEDN
!
      DO i = 1, nreg
         ALLOCATE( reg(i)%q(1:reg(i)%n), reg(i)%wt(1:reg(i)%n),                                        &
                   reg(i)%p(1:reg(i)%n,1:reg(i)%n), reg(i)%dp(1:reg(i)%n,1:reg(i)%n),                  &
                   reg(i)%ddp(1:reg(i)%n,1:reg(i)%n), reg(i)%q_fac(1:reg(i)%n),                        &
                   reg(i)%inv_q_fac(1:reg(i)%n), reg(i)%inv_sqrt_q_fac(1:reg(i)%n) )
!
!
         Call gauss(reg(i)%q, reg(i)%wt, reg%(i)%edge, reg(i)%p, reg(i)%dp, reg(i)%ddp,                &
                    reg(i)%type_quadrature, reg(i)%fixed_point, reg(i)%n)
!
                   reg(i)%q_fac(1:reg(i)%n) = one
                   reg(i)%inv_q_fac(1:reg(i)%n) = one 
                   reg(i)%inv_sqrt_q_fac(1:reg(i)%n) = one
     END DO
  ELSE IF (rho_z%axis == 'rho') THEN
!
     DO i = 1, nreg
        ALLOCATE( reg(i)%q(1:reg(i)%n), reg(i)%wt(1:reg(i)%n),                                        &
                  reg(i)%p(1:reg(i)%n,1:reg(i)%n), reg(i)%dp(1:reg(i)%n,1:reg(i)%n),                  &
                  reg(i)%ddp(1:reg(i)%n,1:reg(i)%n), reg(i)%q_fac(1:reg(i)%n),                        &
                  reg(i)%inv_q_fac(1:reg(i)%n), reg(i)%inv_sqrt_q_fac(1:reg(i)%n) )
!
       Call gauss(reg(i)%q, reg(i)%wt, reg%(i)%edge, reg(i)%p, reg(i)%dp, reg(i)%ddp,                &
                  reg(i)%type_quadrature, reg(i)%fixed_point, reg(i)%n)
!
       reg(i)%q_fac(1:reg(i)%n) = reg(i)%q(1:reg(i)%n)
       reg(i)%inv_q_fac(1:reg(i)%n) = one / reg(i)%q(1:reg(i)%n)
       reg(i)%inv_sqrt_q_fac(1:reg(i)%n) = Sqrt ( reg(i)%inv_q_fac(1:reg(i)%n) )
!
     END DO
!
  ELSE
     Call lnkerr('error in axis variable')
  END IF
!
  IF (prn == .true.) THEN
      DO i = 1, nreg
         call prntfmn('q-'//itoc(i),reg(i)%q,reg(i)%n,1,reg(i)%n,1,iout,'e')
         call prntfmn('wt-'//itoc(i),reg(i)%wt,reg(i)%n,1,reg(i)%n,1,iout,'e')
         call prntfmn('qr_factor-'//itoc(i),reg(i)%q_fac,reg(i)%n,1,reg(i)%n,1,iout,'e')
         call prntfmn('inverse_qr_factor-'//itoc(i),reg(i)%inv_q_fac,reg(i)%n,1,reg(i)%n,1,iout,'e')
         call prntfmn('inverse_sqrt_qr_factor-'//itoc(i),reg(i)%inv_sqrt_q_fac,                   &
                                                                reg(i)%n,1,reg(i)%n,1,iout,'e'))
         Call prntfmn('unnormalized polynomials-'//itoc(i),reg(i)%pr,reg(i)%n,reg(i)%n,           &
                                                           reg(i)%n,reg(i)%n,iout,'e'
         Call prntfmn('first derivative of unnormalized polynomials-'//itoc(i),reg(i)%dpr,        &
                                                                               reg(i)%n,          &
                                                                               reg(i)%n,          &
                                                                               reg(i)%n,          &
                                                                               reg(i)%n,iout,'e')
         Call prntfmn('second derivative of unnormalized polynomials-'//itoc(i),reg(i)%ddpr,      &
                                                                               reg(i)%n,          &
                                                                               reg(i)%n,          &
                                                                               reg(i)%n,          &
                                                                               reg(i)%n,iout,'e')
      END DO
  END IF  
END SUBROUTINE Cylindrical_Grid
!***********************************************************************
!***********************************************************************
!deck Spheroidal_Grid
!***begin prologue     Spheroidal_Grid
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Compute sector grid, weights and coordinate factors
!***                   needed to construct the KE matrix elements.

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE Spheroidal_Grid(xi_eta,reg)
  IMPLICIT NONE
  TYPE(spheroidal)              :: xi_eta
  TYPE(regional)                :: reg
  CHARACTER(LEN=3)                           :: itoc
!
!
  IF (xi_eta%axis == 'eta') THEN
!
      DO i = 1, nreg
         ALLOCATE( reg(i)%q(1:reg(i)%n), reg(i)%wt(1:reg(i)%n),                                        &
                   reg(i)%p(1:reg(i)%n,1:reg(i)%n), reg(i)%dp(1:reg(i)%n,1:reg(i)%n),                  &
                   reg(i)%ddp(1:reg(i)%n,1:reg(i)%n), reg(i)%q_fac(1:reg(i)%n),                        &
                   reg(i)%inv_q_fac(1:reg(i)%n), reg(i)%inv_sqrt_q_fac(1:reg(i)%n) )
!
         Call gauss(reg(i)%q, reg(i)%wt, reg%(i)%edge, reg(i)%p, reg(i)%dp, reg(i)%ddp,                &
                   reg(i)%type_quadrature, reg(i)%fixed_point, reg(i)%n)
!
                   reg(i)%q_fac(1:reg(i)%n) = one - reg(i)%q(1:reg(i)%n) * reg(i)%q(1:reg(i)%n)
                   reg(i)%inv_q_fac(1:reg(i)%n) = one / reg(i)%q_fac(1:reg(i)%n)
                   reg(i)%inv_sqrt_q_fac(1:reg(i)%n) = Sqrt ( reg(i)%inv_q_fac(1:reg(i)%n) )
      END DO
!
  ELSE IF (xi_eta%axis == 'xi' ) THEN
      DO i = 1, nreg
         ALLOCATE( reg(i)%q(1:reg(i)%n), reg(i)%wt(1:reg(i)%n),                                        &
                   reg(i)%p(1:reg(i)%n,1:reg(i)%n), reg(i)%dp(1:reg(i)%n,1:reg(i)%n),                  &
                   reg(i)%ddp(1:reg(i)%n,1:reg(i)%n), reg(i)%q_fac(1:reg(i)%n),                        &
                   reg(i)%inv_q_fac(1:reg(i)%n), reg(i)%inv_sqrt_q_fac(1:reg(i)%n) )
!
         Call gauss(reg(i)%q, reg(i)%wt, reg%(i)%edge, reg(i)%p, reg(i)%dp, reg(i)%ddp,                &
                    reg(i)%type_quadrature, reg(i)%fixed_point, reg(i)%n)
!
         reg(i)%q_fac(1:reg(i)%n) = reg(i)%q(1:reg(i)%n) * reg(i)%q(1:reg(i)%n) - one
         reg(i)%inv_q_fac(1:reg(i)%n) = one / reg(i)%q_fac(1:reg(i)%n)
         reg(i)%inv_sqrt_q_fac(1:reg(i)%n) = Sqrt ( reg(i)%inv_q_fac(1:reg(i)%n) )
      END DO
!
  ELSE
      Call lnkerr('error in axis variable')
  END IF
  IF (prn == .true.) THEN
      DO i = 1, nreg
         call prntfmn('q-'//itoc(i),reg(i)%q,reg(i)%n,1,reg(i)%n,1,iout,'e')
         call prntfmn('wt-'//itoc(i),reg(i)%wt,reg(i)%n,1,reg(i)%n,1,iout,'e')
         call prntfmn('qr_factor-'//itoc(i),reg(i)%q_fac,reg(i)%n,1,reg(i)%n,1,iout,'e')
         call prntfmn('inverse_qr_factor-'//itoc(i),reg(i)%inv_q_fac,reg(i)%n,1,reg(i)%n,1,iout,'e')
         call prntfmn('inverse_sqrt_qr_factor-'//itoc(i),reg(i)%inv_sqrt_q_fac,                   &
                                                         reg(i)%n,1,reg(i)%n,1,iout,'e'))
         Call prntfmn('unnormalized polynomials-'//itoc(i),reg(i)%pr,reg(i)%n,reg(i)%n,           &
                                                           reg(i)%n,reg(i)%n,iout,'e'
         Call prntfmn('first derivative of unnormalized polynomials-'//itoc(i),reg(i)%dpr,        &
                                                                               reg(i)%n,          &
                                                                               reg(i)%n,          &
                                                                               reg(i)%n,          &
                                                                               reg(i)%n,iout,'e')
         Call prntfmn('second derivative of unnormalized polynomials-'//itoc(i),reg(i)%ddpr,      &
                                                                               reg(i)%n,          &
                                                                               reg(i)%n,          &
                                                                               reg(i)%n,          &
                                                                               reg(i)%n,iout,'e')
     END DO
  END IF  
END SUBROUTINE Spheroidal_Grid
!***********************************************************************
!***********************************************************************
!deck Global_Points_Weights
!***begin prologue     Global_Points_Weights
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Compute sector grid, weights and coordinate factors
!***                   needed to construct the KE matrix elements.

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE Global_Points_Weights(grid,reg)
  IMPLICIT NONE
  TYPE(coordinates)        :: grid
  TYPE(regional)           :: reg
!
  nphy=0
  DO i = 1, nreg
     nphy = nphy + reg(i)%n - 2
  END DO
!                                                                                                                         
! add one bridge function between each interval and the extreme left and right points.                                                                      
!                                                                                                                         
  nphy = nphy + nreg + 1
!                                                                                                                      
  nglobal = nphy
  IF(reg(1)%drop(1) == .true.) THEN
        nphy=nphy-1
  END IF
  IF(reg(nreg)%drop(2) == .true.) THEN
        nphy=nphy-1
  END IF
  ALLOCATE(grid%points(1:nphy),grid%weights(1:nphy))
  grid%weights(1:nphy) = zero
  i = 1
  begin = 1
  end = reg(i)%n
  first = 1
  last = reg(i)%n
  IF ( reg(i)%drop(1) == .true. ) THEN
       first = 2
  END IF
  end = last - first + 1
  IF ( nreg == 1) THEN
       IF ( reg(i)%drop(2) == .true.) THEN
            last = reg(i)%n - 1
            end = last - first + 1
       END IF
       grid%points(begin:end)  = reg(nreg)%q(first:last) 
       grid%weights(begin:end) = grid%weights(begin:end) + reg(i)%wt(first:last)  
  ELSE
       grid%points(begin:end)  = reg(i)%q(first:last) 
       grid%weights(begin:end) = grid%weights(begin:end) + reg(i)%wt(first:last)  
       DO i = 2, nreg - 1
          begin = end
          end =  end + reg(i)%n - 1
          grid%points(begin:end) = reg(i)%q(1:reg(i)%n) 
          grid%weights(begin:end) = grid%weights(begin:end) + reg(i)%wt(1:reg(i)%n)  
       END DO
       i = nreg
       begin = end
       end = end + reg(i)%n - 1
       last = reg(i)%n
       IF (reg(i)%drop(2) == .true. ) THEN
           last = reg(i)%n - 1  
           end = end - 1
       END IF
       grid%points(begin:end) = reg(i)%q(1:reg(i)%n) 
       grid%weights(begin:end) = grid%weights(begin:end) + reg(i)%wt(1:reg(i)%n)  
  END IF
  call prntfmn('final global grid points',grid%points,nphy,1,nphy,1,iout,'e') 
  call prntfmn('final global grid weights',grid%weights,nphy,1,nphy,1,iout,'e') 
  R_max = reg(nreg)%q(reg(nreg)%n)
  write(iout,2) R_max
1 Format(/,10x,'Region = ',i4)
2 Format(/,10x,'Last Grid Point = ',e15.8)
END SUBROUTINE Global_Points_Weights
!***********************************************************************
!***********************************************************************
!deck Normalization
!***begin prologue     Normalization
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Compute sector grid, weights and coordinate factors
!***                   needed to construct the KE matrix elements.

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE Normalization(reg)
  IMPLICIT NONE
  TYPE(regional)                :: reg
!
!
  IF ( nreg == 1) THEN
!
       i = 1
       ALLOCATE(reg(i)%normalization(1:reg(i)%n) )
!
!      Only one region.  No endpoint corrections required.
!
       Call Norm ( reg(i)%wt, reg(i)%normalization, reg(i)%n, i)
  ELSE
!

!                  First region of multi-region domain.  Correction  at  
!                  right endpoint needed from first function in second region.
       i = 1
       ALLOCATE( reg(i)%normalization(1:reg(i)%n) )
!
       Call Norm ( reg(i)%wt, reg(i)%normalization, reg(i)%n, i, wt_right_end=reg(i+1)%wt(1))
!       
       DO i = 2, nreg - 1
          ALLOCATE( reg(i)%normalization(1:reg(i)%n) )
!
!                  General case.  Correction at both the left and right
!                  endpoints needed.  The correction at the left enpoint
!                  requires the last weight from the previous region while
!                  the right endpoint correction requires the first weight
!                  from the next region.
!
          Call Norm ( reg(i)%wt, reg(i)%normalization, reg(i)n, i,          &      
                      wt_left_end=reg(i-1)%wt(reg(i-1)%n)),                 &
                      wt_right_end=reg(i+1)%wt(1))                     
!
       END DO
!
!                  Last region.  Correct the left end point using the last
!                  weight from the previous region.
!
       i = nreg
!
       Call Norm ( reg(i)%wt, reg(i)%normalization, reg(i)%n, i             &
                   wt_left_end=reg(i-1)%wt(reg(i-1)%n) )
  END IF
  IF (prn == .true. ) THEN
      DO i = 1, nreg
          Call prntfmn('normalization-'//itoc(i)title,reg(i)%normalization,reg(i)%n,1,        &
                                                                           reg(i)%n,1,iout,'e')      
      END DO
  END IF
END SUBROUTINE Normalization
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
!deck gauss.f
!***begin prologue     gauss
!***date written       000702   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr
!***
!***author             schneider, b. i.(nsf)
!***source             gauss
!***purpose            Points, weights and coordinate functions for generalized
!***                   Gauss quadratures.
!***description        Lanczos recursion using a reference weight function
!***                   is used to generate the points and weights of Gauss quadratures
!***                   for generalized weight functions.  The eigenvector matrix
!***                   of the tridiagonal matrix is used to compute the
!***                   coordinate functions, their first and second derivatives.

!***references         see papers and notes appended.

!***routines called    iosys, util and mdutil
!***end prologue       gauss

!     This is the main library routine to compute the orthogonal and
!     coordinate DVR functions for general weight functions.  The approach
!     is to use a reference quadrature to compute the $\alpha$ and $\beta$
!     recursion coefficients satisfied by the orthogonal polynonials.  The
!     three term recursion relationship is then diagonalized to obtain the
!     generalized points and weights.  The eigenvectors are used to transform
!     from the orthogonal polynomial to coordinate representation.  First
!     and second derivatives are also computed by taking the derivative of the
!     recursion relationship.

  SUBROUTINE gauss(q,wt,edge,p,dp,ddp,type_quadrature,fixed_point,n)
  IMPLICIT NONE
  REAL(idp), DIMENSION(:)                :: q
  REAL(idp), DIMENSION(:)                :: wt
  REAL(idp), DIMENSION(:)                :: edge
  REAL(idp), OPTIONAL, DIMENSION(:,:)    :: p
  REAL(idp), OPTIONAL, DIMENSION(:,:)    :: dp
  REAL(idp), OPTIONAL, DIMENSION(:,:)    :: ddp
  CHARACTER(LEN=*)                       :: type_quadrature
  INTEGER                                :: fixed_point
  INTEGER                                :: n
  REAL(idp), DIMENSION(2)                :: endpts
  REAL(idp), DIMENSION(2)                :: ptfix
  DATA ptfix / -1.d0, 1.d0 /
!
  write(iout,1)
!
! If the weight function is a one of the classical weight functions the
! points and weights are known analytically and after computing them we
! go directly to getting the coordinate functions.
!
  ALLOCATE( b(1:n) )
  endpts(:)=edge(:)
  IF (type_quadrature == "gauss") THEN
      CALL gaussq('one',n,0.d0,0.d0,0,ptfix,b,q,wt)
  ELSE IF ( type_quadrature == "radau") THEN
      IF (fixed_point = 1)
          ptfix(1) = -1.d0
      ELSE IF ( fixed_point = 2)
          ptfix(1) = 1.d0
      END IF
      CALL gaussq('one',n,0.d0,0.d0,1,ptfix,b,q,wt)
  ELSE IF ( type_quadrature == "lobatto" ) THEN
      ptfix(1) = -1.d0
      ptfix(2) = 1.d0
      CALL gaussq('one',n,0.d0,0.d0,2,ptfix,b,q,wt)
  END IF
  CALL cnvtpt(q,wt,edge,n)
  IF(prn) THEN
     CALL prntfm('final nodes',q,n,1,n,1,iout)
     CALL prntfm('final weights',wt,n,1,n,1,iout)
  END IF
  DEALLOCATE(b)
!  
! Generate the needed functions at all required points.
!  
  CALL cpoly(p,dp,ddp,q,n-1,n,prn)
!  
! The DVR library assumes that the polynomials are $\delta$
! functions at the quadrature points.  Convert to this normalization
!
!
1    FORMAT(/,20X,'Computing DVR Points and Weights')
END SUBROUTINE gauss
!***********************************************************************
!***********************************************************************
!*deck lgngr
   SUBROUTINE LGNGR(p,dp,ddp,x,y,nx,ny,type,drctv,prnt) 
!***begin prologue     lgngr
!***date written       940504   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G% 
!***purpose            lagrange polynomials at arbitrary points.
!***description
!***            
!               
!               
!***references
!
!***routines called
!
!***end prologue       lgngr
!
  USE input_output
  IMPLICIT NONE
  REAL*8, DIMENSION (ny,nx)           :: p
  REAL*8, DIMENSION (ny,nx)           :: dp
  REAL*8, DIMENSION (ny,nx)           :: ddp
  REAL*8, DIMENSION(nx)               :: x
  REAL*8, DIMENSION(ny)               :: y
  REAL*8, DIMENSION(:), ALLOCATABLE   :: xt
  REAL*8, DIMENSION(:), ALLOCATABLE   :: yt
  REAL*8                              :: sn
  REAL*8                              :: ssn
  REAL*8                              :: fac
  LOGICAL, OPTIONAL                   :: prnt
  CHARACTER (LEN = 80)                :: title 
  CHARACTER (LEN = *), OPTIONAL       :: drctv
  CHARACTER (LEN = *), OPTIONAL       :: type
  INTEGER                             :: nx
  INTEGER                             :: ny
  INTEGER                             :: i
  INTEGER                             :: j
  INTEGER                             :: k
  INTEGER                             :: first
  INTEGER                             :: second
  INTEGER                             :: zerfac
  INTEGER                             :: inp
  INTEGER                             :: iout
!
!     generate polynomials and derivatives with respect to x
!
  p(:,:) = 1.d0
  IF (present(type) ) THEN 
      ALLOCATE(xt(nx),yt(ny))
      xt(:) = x(:)
      yt(:) = y(:)
      x(:) = x(:) * x(:)
      y(:) = y(:) * y(:)
  END IF
  DO i=1,ny
     zerfac = 0
     DO j=1,nx
        fac =  y(i) - x(j) 
        IF(abs(fac) <= 1.d-10) THEN
           zerfac = j
        ENDIF  
     END DO
     DO j=1,nx
        DO k = 1, j-1
           p(i,j) = p(i,j) * ( y(i) - x(k) )   &
                           / ( x(j) - x(k) )
        END DO
        DO k=j+1,nx
           p(i,j) = p(i,j) * ( y(i) - x(k) )   &
                           / ( x(j) - x(k) )
        END DO
        IF(present(drctv) ) THEN
            IF ( abs(p(i,j)) > 1.d-10) THEN
                 sn = 0.d0
                 ssn = 0.d0
                 DO k=1,j-1
                    fac = 1.d0/( y(i) - x(k) )
                    sn = sn + fac
                    ssn = ssn + fac*fac
                 END DO
                 DO k=j+1,nx
                    fac = 1.d0/( y(i) - x(k) )
                    sn = sn + fac
                    ssn = ssn + fac*fac
                 END DO                                 
                 dp(i,j) = sn*p(i,j)               
                 ddp(i,j) = sn*dp(i,j) - ssn*p(i,j)
            ELSE
                 first=j
                 second=zerfac
                 IF(first > second) THEN
                    first=zerfac
                    second=j
                 END IF
                 sn = 1.d0
                 ssn = 0.d0
                 DO k=1,first-1
                    fac = 1.d0/( x(j) - x(k) )
                    sn = sn*fac*( y(i) - x(k) )
                    ssn = ssn + 1.d0/(y(i) - x(k))
                 END DO
                 DO k=first+1,second-1
                    fac = 1.d0/( x(j) - x(k) )
                    sn = sn*fac*( y(i) - x(k) )
                    ssn = ssn + 1.d0/( y(i) - x(k) )             
                 END DO
                 DO k=second+1,nx
                    fac = 1.d0/( x(j) - x(k) )
                    sn = sn*fac*( y(i) - x(k) )
                    ssn = ssn + 1.d0/( y(i) - x(k) )             
                 END DO
                 dp(i,j) = sn/( x(j) - x(zerfac) )
                 ddp(i,j) = 2.d0*ssn*dp(i,j)
            END IF                    
        END IF
!
     END DO
  END DO
!
  IF (present(type)) THEN 
      DO i=1,ny
         ddp(i,:) = 2.d0*dp(i,:) + 4.d0 * yt(i) * yt(i) * ddp(i,:) 
         dp(i,:) = 2.d0 * yt(i) * dp(i,:)
      END DO
      x(:) = xt(:)
      y(:) = yt(:)
      DEALLOCATE(xt,yt)
!
  END IF
  IF(present(prnt)) THEN
     title='polynomials'
     call prntfm(title,p,ny,nx,ny,nx,iout)
     IF(present(drctv)) then
        title='derivative of polynomials'
        call prntfm(title,dp,ny,nx,ny,nx,iout)
        title='second derivative of polynomials'
        call prntfm(title,ddp,ny,nx,ny,nx,iout)
     END IF
  END IF
  END SUBROUTINE Lgngr
!***********************************************************************
!***********************************************************************
!deck cnvtpt.f
  SUBROUTINE cnvtpt(pt,wt,endpts,n)
  USE input_output
  IMPLICIT NONE
  INTEGER                                :: n
  REAL(idp), DIMENSION(n)                :: pt
  REAL(idp), DIMENSION(n)                :: wt
  REAL(idp), DIMENSION(2)                :: endpts
  REAL(idp)                              :: f1
  REAL(idp)                              :: f2
  f1 = ( endpts(2)-endpts(1) )*.5D0
  f2 = ( endpts(1) + endpts(2) )*.5D0
  pt =  f1*pt + f2
  wt = wt*f1
END SUBROUTINE cnvtpt
!**********************************************************************
!**********************************************************************
!deck cpoly.f
!***begin prologue     cpoly
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           coordinate eigenfunctions
!***author             schneider, barry (nsf)
!***source
!***purpose            generate coordinate functions.
!***
!***description
!***references
!***routines called
!***end prologue       cpoly

  SUBROUTINE cpoly(cp,dcp,ddcp,pt,n,npt,prn)
  USE input_output
  IMPLICIT NONE
  INTEGER                             :: n
  INTEGER                             :: npt
  REAL(idp), DIMENSION(npt,0:n)       :: cp
  REAL(idp), DIMENSION(npt,0:n)       :: dcp
  REAL(idp), DIMENSION(npt,0:n)       :: ddcp
  REAL(idp), DIMENSION(npt)           :: pt
  LOGICAL                             :: prn
  CHARACTER (LEN=80)                  :: title
  CALL lgngr(cp,dcp,ddcp,pt,pt,npt,npt,drctv='on')
  IF(prn) THEN
     title='coordinate function'
     CALL prntfm(title,cp(1,0),npt,n+1,npt,n+1,iout)
     title='first derivative of coordinate function'
     CALL prntfm(title,dcp(1,0),npt,n+1,npt,n+1,iout)
     title='second derivative of coordinate function'
     CALL prntfm(title,ddcp(1,0),npt,n+1,npt,n+1,iout)
   END IF
END SUBROUTINE cpoly
!***********************************************************************
!***********************************************************************
          END MODULE FEDVR_Grid_Functions
!***********************************************************************
!***********************************************************************
