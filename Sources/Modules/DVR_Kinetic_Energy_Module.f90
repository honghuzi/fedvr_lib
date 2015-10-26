!***********************************************************************
! DVR_Kinetic_Energy_Module
!**begin prologue     DVR_Kinetic_Energy_Module
!**date written       090119   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Calculate the regional kinetic energy matrix elements
!***                  in a FEDVR basis.  It is the basic routine that is first
!***                  called before anything else is done.
!***description       The routines in this module compute the regional kinetic
!***                  energy matrix elements for many DVR basis sets.  For
!***                  the FEDVR basis sets based on Gauss quadatures defined
!***                  on (-1,1) there are cases where one needs to explicitly
!***                  extract a square root singularity before using the FEDVR.
!***                  To do that we define a KE and and Odd_KE type.  The other
!***                  FEDVRs actually are only defined in a single region.
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***                  
!***                  
!***end prologue      DVR_Kinetic_Energy_Moduls
!***********************************************************************
!***********************************************************************
                           MODULE DVR_Kinetic_Energy_Module
                           USE Data_Module
                           USE FEDVR_Shared
                           USE FEDVR_Derived_Types
                           USE Matrix_Print
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!
                            INTERFACE Kinetic_Energy                       
                       MODULE PROCEDURE KE,                             &
                                        Odd_KE,                         &
                                        Fourier_KE,                     &  
                                        Hermite_KE,                     &  
                                        Laguerre_KE  
                            END INTERFACE Kinetic_Energy
!
                      TYPE(REAL_MATRIX)              :: type_real_matrix
                      TYPE(REAL_VECTOR)              :: type_real_vector
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck KE.f
!***begin prologue     KE
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Calculate the parts of the regional matrices
!***                   not depending on the (l,m) quantum numbers.
!***                   Also, the formulas programmed assume that 1) the kinetic energy
!***                   is left in its standard form and no first derivative has been transformed
!***                   away and 2) the volume element is included in the definition of the the matrix
!***                   element.  The programmed expression is what is obtained after 
!***                   integrating by parts.  if there is a singular point at the boundaries, the
!***                   integrand remains well behaved and it is not necessary to make that a quadrature
!***                   point.  This allows us to use Gauss-Radau rules for those elements instead of
!***                   Gauss-Lobatto and you do not have to remove the first basis function.  This was
!***                   first pointed out to me by Brett Esry.  this routine will handle cases where
!***                   the FEDVR basis functions are polynomials.
!***references

!***routines called    iosys, util and mdutil
!***end prologue       Ke
  SUBROUTINE KE(grid,reg_mat,reg_poly)
  IMPLICIT NONE
  TYPE(coordinates)                    :: grid  
  TYPE(matrices), DIMENSION(:)         :: reg_mat  
  TYPE(functions), DIMENSION(:)        :: reg_poly  
  INTEGER                              :: ir
  INTEGER                              :: i
  INTEGER                              :: k
!
!
!    Loop over regions
!
  DO ir = 1, nreg
     ALLOCATE( grid%reg_mat(ir)%tr( npt(ir), npt(ir) ) )
     grid%reg_mat(ir)%tr(:,:) = zero
!
!    Loop over (i,j) 
!
     DO i = 1, npt(ir)
!
!       Sum over quadrature points
!
        DO k = 1, npt(ir)
           grid%reg_mat(ir)%tr(i,1:i) = grid%reg_mat(ir)%tr(i,1:i)   &
                                        -                            &
           grid%reg_pt_wt(ir)%qr_fac(k) * grid%reg_pt_wt(ir)%wtr(k)  &
                                        *                            &
           grid%reg_poly(ir)%dpr(k,i)   * grid%reg_poly(ir)%dpr(k,1:i) 
        END DO
        grid%reg_mat(ir)%tr(1:i,i) = grid%reg_mat(ir)%tr(i,1:i)
     END DO
  END DO
  IF ( drop(2) == .false. ) THEN
       ir = nreg
       i  = npt(ir)
       grid%reg_mat(ir)%tr(i,i) = grid%reg_mat(ir)%tr(i,i)            &
                                +                                     &
       bc * grid%reg_pt_wt(ir)%qr_fac(i) * grid%reg_poly(ir)%pr(i,i)  &
                                         * grid%reg_poly(ir)%pr(i,i) 
  END IF
  IF (prn(4) == .true. ) THEN
      DO i = 1, nreg
         write(iout,1) i
         Call Print_Matrix(type_real_matrix,grid%reg_mat(i)%tr,npt(i),npt(i),        &
                           title = 'raw kinetic energy matrix for even FEDVR')
      END DO
  END IF
1 FORMAT(/,10x,'Region = ',i4)
!
END SUBROUTINE KE
!***********************************************************************
!***********************************************************************
!deck Odd_KE.f
!***begin prologue     Odd_KE
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            this routine will handle FEDVR basis functions when there
!***                   is a square root singularity that needs to be removed explicitly.
!***                   once the square root factor is removed the remaining part of
!***                   each basis function may be taken to be of polynomial form.
!***                   this is needed for certain variables such as the odd m legendre
!***                   or spheroidal functions.  again, the programmed formulas are
!***                   obtained after integration by parts.
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Odd_KE

  SUBROUTINE Odd_KE(grid,reg_mat_odd,reg_poly)
  IMPLICIT NONE
  TYPE(coordinates)                    :: grid  
  TYPE(functions), DIMENSION(:)        :: reg_poly  
  TYPE(odd_matrices), DIMENSION(:)     :: reg_mat_odd  
  INTEGER                              :: ir
  INTEGER                              :: i
  INTEGER                              :: j
  INTEGER                              :: k
!
!
! Loop over the regions
!
  DO ir = 1, nreg
     ALLOCATE( grid%reg_mat_odd(ir)%tr( npt(ir), npt(ir) ) )
     grid%reg_mat_odd(ir)%tr(:,:) = zero
!
!    Loop over the (i,j) with j an implicit loop
!
     DO  i = 1, npt(ir)
!
!        Sum over the quadrature points
!
         DO k = 1, npt(ir)
            grid%reg_mat_odd(ir)%tr(i,1:i) = grid%reg_mat_odd(ir)%tr(i,1:i)    &
                                           -                                   &
            grid%reg_pt_wt(ir)%qr_fac(k)   * grid%reg_pt_wt(ir)%qr_fac(k)      &
                                           *                                   &
            grid%reg_pt_wt(ir)%wtr(k)      * grid%reg_poly(ir)%dpr(k,i)        &
                                           * grid%reg_poly(ir)%dpr(k,1:i)
         END DO
         grid%reg_mat_odd(ir)%tr(i,1:i) = grid%reg_mat_odd(ir)%tr(i,1:i)       &
                                         -                                     &
                                      pre_factor                               &
                                         *                                     &
                         ( grid%reg_pt_wt(ir)%qr(i)                            &
                                         *                                     &
                           grid%reg_pt_wt(ir)%qr_fac(i)                        &
                                         *                                     &
                           grid%reg_poly(ir)%pr(i,i)                           &
                                         *                                     &
                           grid%reg_poly(ir)%dpr(i,1:i) )                      &
                                         *                                     &
                           grid%reg_pt_wt(ir)%wtr(i)
         DO j = 1 , i
            grid%reg_mat_odd(ir)%tr(i,j) = grid%reg_mat_odd(ir)%tr(i,j)        &
                                         -                                     &
                                      pre_factor                               &
                                         *                                     &
                         ( grid%reg_pt_wt(ir)%qr(j)                            &
                                      *                                        &
                           grid%reg_pt_wt(ir)%qr_fac(j)                        &
                                      *                                        &
                           grid%reg_poly(ir)%pr(j,j)                           &
                                      *                                        &
                           grid%reg_poly(ir)%dpr(j,i) )                        &
                                      *                                        &
                           grid%reg_pt_wt(ir)%wtr(j)         
         END DO
         grid%reg_mat_odd(ir)%tr(i,i)   = grid%reg_mat_odd(ir)%tr(i,i)         &
                                        -                                      &
         grid%reg_pt_wt(ir)%qr(i)       * grid%reg_pt_wt(ir)%qr(i)             &
                                        *                                      &
         grid%reg_poly(ir)%pr(i,i)      * grid%reg_poly(ir)%pr(i,i)            &
                                        *                                      &
                             grid%reg_pt_wt(ir)%wtr(i)         
     END DO
!
!    Symmetrize
!
     DO  i = 1, npt(ir)
         grid%reg_mat_odd(ir)%tr(i,1:i) =                                      &
                          grid%reg_pt_wt(ir)%inv_sqrt_qr_fac(i)                &
                                        *                                      &
                          grid%reg_mat_odd(ir)%tr(i,1:i)                       &
                                        *                                      &
                          grid%reg_pt_wt(ir)%inv_sqrt_qr_fac(1:i)
         grid%reg_mat_odd(ir)%tr(1:i,i) = grid%reg_mat_odd(ir)%tr(i,1:i)
!         write(iout,*) (grid%reg_mat_odd(ir)%tr(i,j), j=1,i)
     END DO
  END DO
!
  IF ( drop(2) == .false. ) THEN
       ir = nreg
       i = npt(ir)
       grid%reg_mat_odd(ir)%tr(i,i)    = grid%reg_mat_odd(ir)%tr(i,i)             &
                                       +                                          &
       bc * grid%reg_pt_wt(ir)%qr(i)   * grid%reg_pt_wt(ir)%qr(i)                 &
                                       *                                          &
             grid%reg_poly(ir)%pr(i,i) * grid%reg_poly(ir)%pr(i,i)               
  END IF

  IF (prn(4) == .true. ) THEN
      DO i = 1, nreg
         Call Print_Matrix(type_real_matrix,grid%reg_mat_odd(i)%tr,npt(i),npt(i), &
                           title = 'raw kinetic energy matrix for odd FEDVR')
      END DO
  END IF
!
END SUBROUTINE Odd_KE
!***********************************************************************
!***********************************************************************
!***deck Fourier_KE.f90
!***begin prologue     Fourier_KE
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)                                                               
!***keywords           kinetic energy, hamiltonian, fourier basis                                                    
!***author             schneider, barry (nsf)                                                                        
!***source                                                                                                           
!***purpose            generate matrix elements                                                                      
!***                   for a fourier expansion. only a single element
!***                   is considered. there is no integration by parts and
!***                   what is produced is the discretized second derivative
!***                   matrix. 
!***references                                                                                                       
  SUBROUTINE Fourier_KE(grid,reg_mat_fourier)
  IMPLICIT NONE
  TYPE(coordinates)                      :: grid  
  TYPE(fourier_matrices), DIMENSION(:)   :: reg_mat_fourier  
  INTEGER                                :: i
  INTEGER                                :: j
  INTEGER                                :: k
  INTEGER                                :: l
  INTEGER                                :: m
  INTEGER                                :: ir
  INTEGER                                :: delim
  INTEGER                                :: pre
  INTEGER                                :: ifac
  REAL(idp)                              :: kii
  REAL(idp)                              :: kij
  REAL(idp)                              :: cosfac
  REAL(idp)                              :: sinfac_2
  DO ir = 1, nreg
     ALLOCATE( grid%reg_mat_fourier(ir)%tr( npt(ir), npt(ir) ) ) 
     j=( npt(ir) - 1 )/2
     kii=-pi*pi*( npt(ir)*npt(ir) - one)/(three*box*box)
     kij=-two*pi*pi/(box*box)
     DO  i=1,npt(ir)
         grid%reg_mat_fourier(ir)%tr(i,i) = kii
     END DO
     k=-j
     DO i=1,npt(ir)
        l=-j
        DO m=1,i-1
           delim=k-l
           cosfac=cos(pi*delim/npt(ir))
           sinfac_2=sin(pi*delim/npt(ir))
           sinfac_2=sinfac_2 * sinfac_2
           pre=1
           ifac=delim - 2*(delim/2)
           if(ifac == 1) THEN
              pre=-1
           END IF
           grid%reg_mat_fourier(ir)%tr(i,m) = kij * pre * cosfac / sinfac_2
           grid%reg_mat_fourier(ir)%tr(m,i) = grid%reg_mat_fourier(ir)%tr(i,m)
           l = l + 1
        END DO
        k = k + 1
     END DO
     IF (prn(4) == .true. ) THEN
         Call Print_Matrix(type_real_matrix,grid%reg_mat_fourier(ir)%tr,npt(ir),npt(ir), &
                           title = 'Normalized fourier kinetic energy matrix')
     END IF
  END DO
1 Format(/,5x,'fourier weight = ',e15.8)
END SUBROUTINE Fourier_KE
!***********************************************************************
!***********************************************************************
!deck Hermite_KE.f
!***begin prologue     Hermite_KE
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           kinetic energy with Hermite weight function
!***author             schneider, barry (nsf)
!***source
!***purpose            generate kinetic energy matrix elements
!***                   for a hermite weight function.  only one region and
!***                   no integration by parts has been used.
!***
!***description
!***references
!***routines called
!***end prologue       Hermite_KE
  SUBROUTINE Hermite_KE(grid,reg_mat_hermite,reg_poly)
  IMPLICIT NONE
  TYPE(coordinates)                      :: grid  
  TYPE(hermite_matrices), DIMENSION(:)   :: reg_mat_hermite  
  TYPE(functions), DIMENSION(:)          :: reg_poly  
  INTEGER                                :: i 
  INTEGER                                :: j
  INTEGER                                :: ir
  DO ir = 1, nreg
     ALLOCATE( grid%reg_mat_hermite(ir)%tr( npt(ir), npt(ir) ) )
     grid%reg_mat_hermite(ir)%tr(:,:) = zero
     DO  i=1,npt(ir)
         DO  j=1,npt(ir)
             grid%reg_mat_hermite(ir)%tr(i,:)                                     &
                                    =                                             &
             grid%reg_mat_hermite(ir)%tr(i,:)                                     &
                                    -                                             &
                          grid%reg_pt_wt(ir)%qr_fac(j)                            &
                                        *                                         &
                          grid%reg_pt_wt(ir)%wtr(j)                               &
                                        *                                         &
             grid%reg_poly(ir)%dpr(j,i) * grid%reg_poly(ir)%dpr(j,1:i)
         END DO
         grid%reg_mat_hermite(ir)%tr(i,1:i)  = grid%reg_mat_hermite(ir)%tr(i,1:i) &
                                          +                                       &
                         grid%reg_poly(ir)%pr(i,i)                                &
                                          *                                       &
                         grid%reg_poly(ir)%dpr(i,1:i)                             &
                                          *                                       &
                         grid%reg_pt_wt(ir)%qr_fac(i)                             &
                                          *                                       &
                         grid%reg_pt_wt(ir)%qr(i)                                 &
                                          *                                       &
                         grid%reg_pt_wt(ir)%wtr(i)                                 
         DO j = 1, i
            grid%reg_mat_hermite(ir)%tr(i,j) =  grid%reg_mat_hermite(ir)%tr(i,j)  &
                                         +                                        &
                            grid%reg_poly(ir)%pr(j,j)                             &
                                         *                                        &
                            grid%reg_poly(ir)%dpr(j,i)                            &
                                         *                                        &
                            grid%reg_pt_wt(ir)%qr_fac(j)                          &
                                         *                                        &
                            grid%reg_pt_wt(ir)%qr(j)                              &
                                         *                                        &
                           grid%reg_pt_wt(ir)%wtr(j)                              
         END DO       
         grid%reg_mat_hermite(ir)%tr(1:i,i)  = grid%reg_mat_hermite(ir)%tr(i,1:i)
         grid%reg_mat_hermite(ir)%tr(i,i)    = grid%reg_mat_hermite(ir)%tr(i,i)   &
                                           -                                      &
                            grid%reg_pt_wt(ir)%qr(i)                              &
                                        *                                         &
                            grid%reg_pt_wt(ir)%qr(i)                              &
                                        *                                         &
                            grid%reg_pt_wt(ir)%qr_fac(i)                          &
                                        *                                         &
                            grid%reg_pt_wt(ir)%wtr(i)                             &
                                        *                                         &
                            grid%reg_poly(ir)%pr(i,i)                             &
                                        *                                         &
                            grid%reg_poly(ir)%pr(i,i)                              
     END DO
  END DO
  DO ir = 1, nreg
     write(iout,1) ir
     IF (prn(4) == .true. ) THEN
         Call Print_Matrix(type_real_matrix,grid%reg_mat_hermite(ir)%tr,npt(ir),npt(ir), &
                           title = 'raw hermite kinetic energy matrix')
     END IF
  END DO
1 FORMAT(/,10x,'Region = ',i4)
END SUBROUTINE Hermite_KE
!***********************************************************************
!**********************************************************************
!deck LaGuerre_KE.f
!***begin prologue     LaGuerre_KE
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           kinetic energy with Laguerre weight function
!***author             schneider, barry (nsf)
!***source
!***purpose            generate kinetic energy matrix elements
!***                   for a laguerre weight function. only one region
!***                   and no integration by parts has been used.
!***description
!***references
!***routines called
!***end prologue       LaGuerre_KE
  SUBROUTINE LaGuerre_KE(grid,reg_mat_laguerre,reg_poly)
  IMPLICIT NONE
  TYPE(coordinates)                       :: grid  
  TYPE(laguerre_matrices), DIMENSION(:)   :: reg_mat_laguerre  
  TYPE(functions), DIMENSION(:)           :: reg_poly  
  INTEGER                                 :: i
  INTEGER                                 :: j
  ALLOCATE( grid%reg_mat_laguerre(1)%tr( npt(1), npt(1) ) )
  grid%reg_mat_laguerre(1)%tr(:,:) = zero
  DO  i=1,npt(1)
      DO j = 1, npt(1)
         grid%reg_mat_laguerre(1)%tr(i,1:i)    = grid%reg_mat_laguerre(1)%tr(i,1:i)     &
                                      -                                                 &
                        grid%reg_poly(1)%dpr(j,i)                                       &
                                      *                                                 &
                        grid%reg_poly(1)%dpr(j,1:i)                                     &
                                      *                                                 &
                        grid%reg_pt_wt(1)%qr_fac(j)                                     &
                                      *                                                 &
                        grid%reg_pt_wt(1)%wtr(j)                                        
      END DO
      grid%reg_mat_laguerre(1)%tr(i,1:i)  = grid%reg_mat_laguerre(1)%tr(i,1:i)          &
                                      +                                                 &
                     grid%reg_poly(1)%pr(i,i)                                           &
                                      *                                                 &
                     grid%reg_poly(1)%dpr(i,1:i)                                        &
                                      *                                                 &
                     grid%reg_pt_wt(1)%qr_fac(i)                                        &
                                      *                                                 &
                     grid%reg_pt_wt(1)%wtr(i)                                           &
                                      *                                                 &
                                     half      
      DO j = 1, i
         grid%reg_mat_laguerre(1)%tr(i,j)     =  grid%reg_mat_laguerre(1)%tr(i,j)       &
                                     +                                                  &
                        grid%reg_pt_wt(1)%wtr(j)                                        &
                                     *                                                  &
                        grid%reg_pt_wt(1)%qr_fac(j)                                     &
                                      *                                                 &
                        grid%reg_poly(1)%pr(j,j)                                        &
                                      *                                                 &
                        grid%reg_poly(1)%dpr(j,i)                                       &
                                      *                                                 &
                                     half      
      END DO
      grid%reg_mat_laguerre(1)%tr(1:i,i)    = grid%reg_mat_laguerre(1)%tr(i,1:i)
      grid%reg_mat_laguerre(1)%tr(i,i)      = grid%reg_mat_laguerre(1)%tr(i,i)          &
                                   -                                                    &
                        grid%reg_pt_wt(1)%qr_fac(i)                                     &
                                   *                                                    &
                        grid%reg_pt_wt(1)%wtr(i)                                        &
                                   *                                                    &
                        grid%reg_poly(1)%pr(i,i)                                        &
                                   *                                                    &
                        grid%reg_poly(1)%pr(i,i)                                        &
                                   *                                                    &
                                quarter
  END DO
  IF (prn(4) == .true. ) THEN
      Call Print_Matrix(type_real_matrix,grid%reg_mat_laguerre(1)%tr,npt(1),npt(1),     &
                        title = 'raw laguerre kinetic energy matrix')
  END IF
END SUBROUTINE LaGuerre_KE
!***********************************************************************
!***********************************************************************
           END MODULE DVR_Kinetic_Energy_Module
!***********************************************************************
!***********************************************************************
