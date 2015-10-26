!***********************************************************************
! DVR_Polynomials_Module
!**begin prologue     DVR_Polynomials_Module
!**date written       090119   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Calculate the piecewise FEDVR functions.
!***                  
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***end prologue      DVR_Polynomials_Module
!***********************************************************************
!***********************************************************************
                           MODULE DVR_Polynomials_Module
                           USE Data_Module
                           USE FEDVR_Shared
                           USE FEDVR_Derived_Types
                           USE Renormalization_Module
                           USE input_output
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
                            INTERFACE LaGrange_Polynomials                             
                       MODULE PROCEDURE Polynomials,                    &
                                        Fourier_Polynomials             
                            END INTERFACE LaGrange_Polynomials
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck Lobatto_Functions.f
!***begin prologue     Lobatto_Functions
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Calculate piecewise dvr functions
!***description        DVR polynomials are computed for a number                   
!***                   of coordinate systems.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Lobatto_Functions

  SUBROUTINE Lobatto_Functions(grid)
  IMPLICIT NONE
  TYPE (coordinates)            :: grid
!
  IF (keyword == 'fourier') THEN
      Call LaGrange_Polynomials(grid,grid%reg_poly_fourier)
!
  ELSE
      Call LaGrange_Polynomials(grid,grid%reg_poly)
!
  END IF
1 FORMAT(/,10x,'Calculate the Regional Basis Functions for Coordinate = ',a4)
2 FORMAT(/,10x,'Region = ',i3,2x,'Number of Polynomials = ', i3)
END SUBROUTINE Lobatto_Functions
!***********************************************************************
!***********************************************************************
!deck Coordinate_Factors.f
!***begin prologue     Coordinate_Factors
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Compute needed prefactors for various
!***                   coordinate systems.  These appear in the KE
!***                   matrix elements.

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE Coordinate_Factors(grid)
  IMPLICIT NONE
  TYPE(coordinates)                    :: grid
  INTEGER                              :: i
  INTEGER                              :: first
  INTEGER                              :: last
  INTEGER                              :: begin
  INTEGER                              :: end
!
  len=lenth(grid%label)
  ALLOCATE(grid%grid_points(1:physical_points),grid%grid_weights(1:physical_points))
  grid%grid_weights(:) = zero
  i = 1
  begin = 1
  end = npt(1)
  first = 1
  last = npt(1)
  IF ( drop(1) == .true. ) THEN
       first = 2
  END IF
  end = last - first + 1
  IF ( nreg == 1) THEN
       IF ( drop(2) == .true.) THEN
            last = npt(1) - 1
            end = last - first + 1
       END IF
       grid%grid_points(begin:end)  = grid%reg_pt_wt(i)%qr(first:last) 
       grid%grid_weights(begin:end) = grid%grid_weights(begin:end) + grid%reg_pt_wt(i)%wtr(first:last)  
  ELSE
       grid%grid_points(begin:end)  = grid%reg_pt_wt(i)%qr(first:last) 
       grid%grid_weights(begin:end) = grid%grid_weights(begin:end) + grid%reg_pt_wt(i)%wtr(first:last)  
       DO i = 2, nreg - 1
          begin = end
          end =  end + npt (i) - 1
          grid%grid_points(begin:end) = grid%reg_pt_wt(i)%qr(1:npt(i)) 
          grid%grid_weights(begin:end) = grid%grid_weights(begin:end) + grid%reg_pt_wt(i)%wtr(1:npt(i))  
       END DO
       i = nreg
       begin = end
       end = end + npt(i) - 1
       last = npt(i)
       IF (drop(2) == .true. ) THEN
           last = npt(i) - 1  
           end = end - 1
       END IF
       grid%grid_points(begin:end) = grid%reg_pt_wt(i)%qr(1:npt(i)) 
       grid%grid_weights(begin:end) = grid%grid_weights(begin:end) + grid%reg_pt_wt(i)%wtr(1:npt(i))  
  END IF
  call Print_Matrix(type_real_vector,grid%grid_points, title = 'final grid points'  )
  Call IOsys('write real "physical grid points" to '//FEDVR_File,physical_points,grid%grid_points,0,' ')
  Call IOsys('write real "physical grid weights" to '//FEDVR_File,physical_points,grid%grid_weights,0,' ')
  R_max = grid%reg_pt_wt(nreg)%qr(npt(nreg))
  Call IOsys('write real "last grid point" to '//FEDVR_File,1,R_max,0,' ')
  write(iout,2) R_max
!
!         Here we compute once and for all some variable that are used over and over again
!         This is for convenience and nothing else and does not take a lot of memory
!         The qr_fac arrays are related to the facors appearing in front of the derivatives
!         in the one-dimensional Schreodinger equation.  For example, r*r, 1/(r*r) and 1/r
!

!
  IF (grid%label(1:len) /= 'fourier' ) THEN  
!
      DO i = 1, nreg
         ALLOCATE ( grid%reg_pt_wt(i)%qr_fac( npt(i) ),             &
                    grid%reg_pt_wt(i)%inv_qr_fac( npt(i) ),         &
                    grid%reg_pt_wt(i)%inv_sqrt_qr_fac( npt(i) ) )
!         For cartesian coordinates they are all one.
!
         IF (grid%label(1:len) == 'x' .or. grid%label(1:len) == 'y' .or. grid%label(1:len) == 'z' ) THEN  
             grid%reg_pt_wt(i)%qr_fac(:) = one
             grid%reg_pt_wt(i)%inv_qr_fac(:) = one 
             grid%reg_pt_wt(i)%inv_sqrt_qr_fac(:) = one 
!
!     ELSE IF (grid%label(1:len) == 'fourier' ) THEN  
!         grid%reg_pt_wt(i)%qr_fac(:) = one
!         grid%reg_pt_wt(i)%inv_qr_fac(:) = one 
!         grid%reg_pt_wt(i)%inv_sqrt_qr_fac(:) = one 
!
         ELSE IF (grid%label(1:len) == 'r' ) THEN  
!          For the radial spherical variable the volume element
!          is r*r.
!
              IF ( grid%drop_pt(1) == .true.) THEN
                   grid%reg_pt_wt(i)%qr_fac(:) = one
                   grid%reg_pt_wt(i)%inv_qr_fac(:) = one 
                   grid%reg_pt_wt(i)%inv_sqrt_qr_fac(:) = one 
              ELSE
                   grid%reg_pt_wt(i)%qr_fac(:) = grid%reg_pt_wt(i)%qr(:) * grid%reg_pt_wt(i)%qr(:) 
                   grid%reg_pt_wt(i)%inv_qr_fac(:) = one  / grid%reg_pt_wt(i)%qr_fac(:)
                   grid%reg_pt_wt(i)%inv_sqrt_qr_fac(:) = one / grid%reg_pt_wt(i)%qr(:) 
              END IF
!
         ELSE IF (grid%label(1:len) == 'theta'.or.grid%label(1:len) == 'eta') THEN
!
!          For both the spherical and the spheroidal angular variable
!          we need ( 1 - eta*eta) as the factor where we are using the
!          cosine of the angle as the relevant variable.
!
              grid%reg_pt_wt(i)%qr_fac(:) = one - grid%reg_pt_wt(i)%qr(:) * grid%reg_pt_wt(i)%qr(:)
              grid%reg_pt_wt(i)%inv_qr_fac(:) = one / grid%reg_pt_wt(i)%qr_fac(:)
              grid%reg_pt_wt(i)%inv_sqrt_qr_fac(:) = Sqrt ( grid%reg_pt_wt(i)%inv_qr_fac(:) )
!
         ELSE IF (grid%label(1:len) == 'xi' ) THEN
!
!          For the spheroidal radial variable we need ( xi*xi - 1) 
!          as the factor.
!
              grid%reg_pt_wt(i)%qr_fac(:) = grid%reg_pt_wt(i)%qr(:) * grid%reg_pt_wt(i)%qr(:) - one
              grid%reg_pt_wt(i)%inv_qr_fac(:) = one / grid%reg_pt_wt(i)%qr_fac(:)
              grid%reg_pt_wt(i)%inv_sqrt_qr_fac(:) = Sqrt ( grid%reg_pt_wt(i)%inv_qr_fac(:) )
!
         ELSE IF (grid%label(1:len) == 'rho') THEN
!
!          For cylindrical coordinates we need rho as the factor.
!
              grid%reg_pt_wt(i)%qr_fac(:) = grid%reg_pt_wt(i)%qr(:)
              grid%reg_pt_wt(i)%inv_qr_fac(:) = one / grid%reg_pt_wt(i)%qr(:)
              grid%reg_pt_wt(i)%inv_sqrt_qr_fac(:) = Sqrt ( grid%reg_pt_wt(i)%inv_qr_fac(:) )
!
         END IF
         IF (prn(2) == .true.) THEN
             write(iout,1) i
             call Print_Matrix(type_real_vector,grid%reg_pt_wt(i)%qr_fac,title='qr_factor')
             call Print_Matrix(type_real_vector,grid%reg_pt_wt(i)%inv_qr_fac,title='inverse_qr_factor')
             call Print_Matrix(type_real_vector,grid%reg_pt_wt(i)%inv_sqrt_qr_fac,title='inverse_sqrt_qr_factor')
         END IF  
      END DO
  ELSE
      DO i = 1, nreg
         ALLOCATE ( grid%reg_pt_wt(i)%qr_fac( npt(i) ),             &
                    grid%reg_pt_wt(i)%inv_qr_fac( npt(i) ),         &
                    grid%reg_pt_wt(i)%inv_sqrt_qr_fac( npt(i) ) )
         grid%reg_pt_wt(i)%qr_fac(:) = one
         grid%reg_pt_wt(i)%inv_qr_fac(:) = one 
         grid%reg_pt_wt(i)%inv_sqrt_qr_fac(:) = one 
      END DO
  END IF
1 Format(/,10x,'Region = ',i4)
2 Format(/,10x,'Last Grid Point = ',e15.8)
END SUBROUTINE Coordinate_Factors
!***********************************************************************
!***********************************************************************
!deck Polynomials.f
!***begin prologue     Polynomials
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate piecewise lobatto dvr functions
!***                   these are just the simple interpolation
!***                   polynomials.  any factors that are needed
!**                    to distinguish particular (l,m) aspects
!**                    of the matrix elements are done there.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Polynomials

  SUBROUTINE Polynomials(grid,reg_poly)
  IMPLICIT NONE
  TYPE(coordinates)                   :: grid
  TYPE(functions),  DIMENSION(:)      :: reg_poly
  REAL(idp)                           :: alpha
  REAL(idp)                           :: beta
  INTEGER                             :: i
  CHARACTER (LEN=80 )                 :: title
  LOGICAL                             :: dollar
!
!
!           
 ALLOCATE ( grid%reg_poly(1:nreg) )
  DO  i=1,nreg
      Write(iout,1) i
!
!        calculate the sector functions and their derivatives.
!
      ALLOCATE( grid%reg_pt_wt(i)%qr(npt(i)),                               &
                grid%reg_pt_wt(i)%wtr(npt(i)),                              &
                grid%reg_poly(i)%pr(npt(i),npt(i)),                         &
                grid%reg_poly(i)%dpr(npt(i),npt(i)),                        &
                grid%reg_poly(i)%ddpr(npt(i),npt(i)),                       &
                grid%reg_poly(i)%normalization(npt(i)) )              
!
!         compute the reference points and weights
!
     IF (refwt == 'none') THEN
         Call gauss(grid%reg_pt_wt(i)%qr,                                   &
                    grid%reg_pt_wt(i)%wtr,                                  &
                    edge(i),                                                &
                    typwt,                                                  &
                    grid%reg_poly(i)%pr,                                    &
                    grid%reg_poly(i)%dpr,                                   &
                    grid%reg_poly(i)%ddpr,                                  &
                    npt(i),                                                 &
                    i)
     ELSE
         Call gauss(grid%reg_pt_wt(i)%qr,                                   &
                    grid%reg_pt_wt(i)%wtr,                                  &
                    edge(i),                                                &
                    typwt,                                                  &
                    grid%reg_poly(i)%pr,                                    &
                    grid%reg_poly(i)%dpr,                                   &
                    grid%reg_poly(i)%ddpr,                                  &
                    npt(i),                                                 &
                    i,                                                      &
                    refwt,                                                  &
                    alpha,                                                  &
                    beta)
      END IF
      IF (prn(1) == .true. ) THEN
          Call Print_Matrix(type_real_vector,grid%reg_pt_wt(i)%qr,title = 'points')
          Call Print_Matrix(type_real_vector,grid%reg_pt_wt(i)%wtr,title = 'weights')
      END IF
      IF (prn(3) == .true. ) THEN
          Call Print_Matrix(type_real_matrix,grid%reg_poly(i)%pr,npt(i),npt(i),   &
                            title = 'unnormalized polynomials')
          Call Print_Matrix(type_real_matrix,grid%reg_poly(i)%dpr,npt(i),npt(i),  &
                            title = 'first derivative of unnormalized polynomials')

          Call Print_Matrix(type_real_matrix,grid%reg_poly(i)%ddpr,npt(i),npt(i), &
                            title = 'second derivative of unnormalized polynomials')
      END IF
  END DO
!
!                  To compute most of what is required, it is not necessary
!                  to construct anything else. Since the bridge functions 
!                  span two elements, one can define them at the grid points
!                  but their derivatives are discontinuous across the
!                  sector boundaries. The matrix elements can be constructed 
!                  entirely from re-normalized sector quantities. It is
!                  convenient to have the normalization factors which involve
!                  more than one sector available.  They are computed now.
!
  IF ( nreg == 1) THEN
!
!                  Only one region.  No endpoint corrections required.
!
       i = 1
       Call Norm ( grid%reg_pt_wt(i)%wtr,                                   &
                   grid%reg_poly(i)%normalization,                          &
                   npt(i),                                                  &
                   i)

  ELSE
!
!                  First region of multi-region domain.  Correction  at  
!                  right endpoint needed from first function in second region.
       i = 1
       Call Norm ( grid%reg_pt_wt(i)%wtr,                                   &
                   grid%reg_poly(i)%normalization,                          &
                   npt(i),                                                  &
                   i,                                                       &      
                   wt_right_end=grid%reg_pt_wt(i+1)%wtr(1))
!       
       DO i = 2, nreg - 1
!
!                  General case.  Correction at both the left and right
!                  endpoints needed.  The correction at the left enpoint
!                  requires the last weight from the previous region while
!                  the right endpoint correction requires the first weight
!                  from the next region.
!
          Call Norm ( grid%reg_pt_wt(i)%wtr,                                &
                      grid%reg_poly(i)%normalization,                       &
                      npt(i),                                               &
                      i,                                                    &      
                      wt_left_end=grid%reg_pt_wt(i-1)%wtr(npt(i-1)),        &
                      wt_right_end=grid%reg_pt_wt(i+1)%wtr(1))                     
!
       END DO
!
!                  Last region.  Correct the left end point using the last
!                  weight from the previous region.
!
       i = nreg
!
          Call Norm ( grid%reg_pt_wt(i)%wtr,                                &
                      grid%reg_poly(i)%normalization,                       &
                      npt(i),                                               &
                      i,                                                    &
                      wt_left_end=grid%reg_pt_wt(i-1)%wtr(npt(i-1)))
  END IF
  IF (prn(3) == .true. ) THEN
      DO i = 1, nreg
         write(iout,1) i 
         Call Print_Matrix(type_real_vector,grid%reg_poly(i)%normalization, & 
                           title = 'normalization')      
      END DO
  END IF
1 FORMAT(/,10x,'Region = ',i4)
END SUBROUTINE Polynomials
!***********************************************************************
!***********************************************************************
!***deck Fourier_Polynomials.f                                                                          
!***begin prologue     Fourier_Polynomials.f
!***date written       000702   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr
!***
!***author             schneider, b. i.(nsf)
!***source             Fourier_Points_Weights
!***purpose            Points, and weights for fourier dvr
!***description
!***
!***
!***
!*** 
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Fourier_Polynomials
  SUBROUTINE Fourier_Polynomials(grid,reg_poly_fourier)
  IMPLICIT NONE
  TYPE(coordinates)                      :: grid
  TYPE(fourier_functions),  DIMENSION(:) :: reg_poly_fourier
  INTEGER                                :: i
  INTEGER                                :: j
  INTEGER                                :: k
  INTEGER                                :: ir
  REAL(idp)                              :: fac
  CHARACTER (LEN=80)                     :: title
  ALLOCATE(grid%reg_poly_fourier(1:nreg))
  DO ir = 1, nreg
     ALLOCATE(grid%reg_pt_wt(ir)%qr(1:npt(ir)), grid%reg_pt_wt(ir)%wtr(1:npt(ir)),            &
              grid%reg_poly_fourier(ir)%pr(1:npt(ir),1:npt(ir)),                              &
              grid%reg_poly_fourier(ir)%normalization(1:npt(ir)) )
     box=edge(ir+1)-edge(ir)
     deltax=box/npt(ir)
     grid%reg_pt_wt(ir)%wtr(:) = deltax
     j = ( npt(ir) - 1)/2
     fac = edge(ir) + .5d0*box
     k=-j
     DO i=1,npt(ir)
        grid%reg_pt_wt(ir)%qr(i) = deltax * k + fac
        k=k+1
     END DO
     IF (prn(1)  == .true.) THEN
         CALL Print_Matrix(type_real_vector,grid%reg_pt_wt(ir)%qr,title='fourier points')
         CALL Print_Matrix(type_real_vector,grid%reg_pt_wt(ir)%wtr,title='fourier weights')
     END IF
     Call Norm ( grid%reg_pt_wt(ir)%wtr,                                                 &
                 grid%reg_poly_fourier(ir)%normalization,                                &
                 npt(ir),                                                                &
                 ir)
     grid%reg_poly_fourier(ir)%pr(:,:) = zero
     DO i=1,npt(ir)
        grid%reg_poly_fourier(ir)%pr(i,i)  = one
     END DO
     IF (prn(3) == .true. ) THEN
         Call Print_Matrix(type_real_matrix,grid%reg_poly_fourier(ir)%pr,npt(ir),npt(ir), &
                           title = 'unnormalized polynomials')
     END IF
  END DO
  IF (prn(3) == .true. ) THEN
      DO i = 1, nreg
         write(iout,1) i 
         Call Print_Matrix(type_real_vector,grid%reg_poly_fourier(i)%normalization,title = 'normalization')      
      END DO
  END IF
1 FORMAT(/,10x,'Region = ',i4)
  END SUBROUTINE Fourier_Polynomials
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

  SUBROUTINE gauss(q,wt,edge,typwt,p,dp,ddp,nord,reg_number,         &
                   refwt,alpha,beta)
  IMPLICIT NONE
  INTEGER                                :: nord
  REAL(idp), DIMENSION(:)                :: q
  REAL(idp), DIMENSION(:)                :: wt
  REAL(idp), OPTIONAL, DIMENSION(:,:)    :: p
  REAL(idp), OPTIONAL, DIMENSION(:,:)    :: dp
  REAL(idp), OPTIONAL, DIMENSION(:,:)    :: ddp
  REAL(idp), OPTIONAL                    :: alpha
  REAL(idp), OPTIONAL                    :: beta
  REAL(idp), DIMENSION(nreg+1)           :: edge
  CHARACTER (LEN=*)                      :: typwt
  CHARACTER (LEN=*), OPTIONAL            :: refwt
  INTEGER                                :: n
  INTEGER                                :: i
  INTEGER                                :: j
  INTEGER                                :: nq
  INTEGER                                :: reg_number
  INTEGER                                :: n_fixed
  INTEGER                                :: max_num = 10
  CHARACTER (LEN=80)                     :: title
  REAL(idp), DIMENSION(2)                :: endpts
  REAL(idp), DIMENSION(2)                :: ptse
  REAL(idp), DIMENSION(2)                :: ptfix
  REAL(idp)                              :: nrm
  REAL(idp), DIMENSION(:), ALLOCATABLE   :: a
  REAL(idp), DIMENSION(:), ALLOCATABLE   :: b
  REAL(idp), DIMENSION(:), ALLOCATABLE   :: r
  REAL(idp), DIMENSION(:), ALLOCATABLE   :: rwt
  REAL(idp), DIMENSION(:), ALLOCATABLE   :: wtfn
  REAL(idp), DIMENSION(:,:), ALLOCATABLE :: ply
  REAL(idp), DIMENSION(:), ALLOCATABLE   :: scr
  REAL(idp), DIMENSION(:), ALLOCATABLE   :: arg
  REAL(idp), DIMENSION(:), ALLOCATABLE   :: scrat
  REAL(idp), DIMENSION(:), ALLOCATABLE   :: q_diff
  REAL(idp), DIMENSION(:,:), ALLOCATABLE :: eigv
  REAL(idp)                              :: max_diff
  REAL(idp)                              :: ddot
  REAL(idp)                              :: mu
  REAL(idp)                              :: tol=1.d-10
  DATA ptfix / -1.d0, 1.d0 /
!
  write(iout,1)
!
! If the weight function is a one of the classical weight functions the
! points and weights are known analytically and after computing them we
!  go directly to getting the coordinate functions.
!
  endpts=edge
  ptse = endpts
  IF (nreg == 1) THEN
      n_fixed=0
      IF(nfix == 1) THEN
         n_fixed=1       
         ptfix(1) = -1.d0
         ptse(1) = endpts(1)
         IF(fix(2)) THEN
            ptfix(1) = 1.d0
            ptse(1) = endpts(2)
         END IF
      ELSE IF(nfix == 2) THEN
         n_fixed=2      
         ptfix(1) = -1.d0
         ptfix(2) = 1.d0
         ptse = endpts
      END IF
  ELSE
      IF(reg_number == 1) THEN 
         n_fixed=1       
         ptfix(1)=1.d0
         ptse(1) = endpts(2)
         IF(fix(1)) THEN
            n_fixed=2
            ptfix(1)=-1.d0
            ptfix(2)=1.d0
            ptse = endpts
         END IF
      ELSE IF(reg_number == nreg) THEN 
            n_fixed=1   
            ptfix(1)=-1.d0  
            ptse(1) = endpts(1)  
            IF(fix(2)) THEN
               n_fixed=2
               ptfix(1)=-1.d0
               ptfix(2)=1.d0
               ptse = endpts
            END IF
      ELSE
         n_fixed=2
         ptfix(1)=-1.d0
         ptfix(2)=1.d0
         ptse = endpts
      END IF
  END IF
! For finite regions
!
  IF ( .not.present(refwt)) THEN  
       write(iout,2) typwt
       ALLOCATE(b(nord))
       IF(typwt == 'one'.or.         &
          typwt == 'chebyshev-1'.or. &
          typwt == 'chebyshev-2'.or. &
          typwt == 'jacobi'.or.      &
          typwt == 'spherical'.or.   &
          typwt == 'cylindrical'.or. &
          typwt == 'legendre') THEN
          CALL gaussq(typwt,nord,0.d0,0.d0,n_fixed,ptfix,b,q,wt)
          CALL cnvtpt(q,wt,edge,nord)
          IF(prn(1)) THEN
             CALL Print_Matrix(type_real_vector,q,title='final nodes')
             CALL Print_Matrix(type_real_vector,wt,title='final weights')
          END IF
!
!       For Infinite Regions
!
       ELSE IF( typwt == 'hermite') THEN
          n_fixed = 0
          CALL gaussq(typwt,nord,0.d0,0.d0,n_fixed,ptfix,b,q,wt)
       ELSE IF(typwt == 'laguerre' .or.                          &
               typwt=='spherical_hermite') THEN
          n_fixed=0
          IF(nfix == 1) THEN
             n_fixed=1       
             ptfix(1) = 1.0d-20
             ptse(1) = ptfix(1)    
          END IF
          CALL gaussq('laguerre',nord,0.d0,0.d0,n_fixed,ptfix,b,q,wt)
       END IF
       DEALLOCATE(b)
  ELSE
!  
!     Compute the non-classical recursion coefficients using the Lanczos method
!     based on computing the needed integrals via a reference quadrature.
     ALLOCATE(a(nord),b(nord),q_diff(nord))
     write(iout,3) typwt, refwt
     nq = nord + nord
     q_diff = 0.d0
     DO i = 1,max_num
        ALLOCATE(r(nq),rwt(nq),wtfn(nq),scr(nq),arg(nq),ply(nq,0:nord))
!  
!        Here we compute the reference nodes and weights.
!  
        CALL gaussq(refwt,nq,0.d0,0.d0,0,ptfix,scr,r,rwt)
!  
!       Convert the Gauss points and weights on [-1.,1.] to [edge(1),edge(2)]
!
        IF (refwt == 'legendre') THEN
            CALL cnvtpt(r,rwt,edge,nq)
        END IF
        IF(prn(2)) THEN
           CALL Print_Matrix(type_real_vector,r,title='reference nodes')
           CALL Print_Matrix(type_real_vector,rwt,title='reference weights')
        END IF
!  
!       Compute the effective weight function.
!
        CALL genrwt(wtfn,r,typwt,nq,refwt,alpha,beta)
        IF(prn(3)) THEN
           CALL Print_Matrix(type_real_vector,wtfn,title='ratio weight factor')
        END IF
!  
!        Generate the recursion coefficients numerically.
!        Initialize the first function and normalize.
!  
!        ply(:,0) = 1.d0
!        arg=ply(:,0)*ply(:,0)*rwt(:)*wtfn(:)
        mu=ddot(nq,rwt,1,wtfn,1)
        WRITE(iout,4) mu
!        arg(:) = r(:)
!  
!        Note that lancz will return the recursion coefficients on
!        the actual interval.  This is consistent with what is done for
!        the known cases.
!  
        CALL lancz(r,a,b,rwt,wtfn,nq,nord)
        CALL modab(a,b,n_fixed,ptse,nord)
        IF(prn(4)) THEN
           CALL Print_Matrix(type_real_vector,a,title='lanczos a coefficients')
           CALL Print_Matrix(type_real_vector,b,title='lanczos b coefficients')
        END IF
!  
!        Get the points and weights and then compute the coordinate
!        functions
!  
        ALLOCATE(eigv(nord,nord),scrat(nord))
        q = a
        scrat=b
!  
!        Generate the non-classical points and weights and the
!        transformation matrix from the orthogonal polynomials
!        to the co-ordinate functions.
!
        CALL genq(q,scrat,wt,eigv,fix,endpts,mu,nord)
        DEALLOCATE(r,rwt,wtfn,arg,scr,ply,eigv,scrat)
        IF(prn(1)) THEN
          CALL Print_Matrix(type_real_vector,q,title='final nodes')
          CALL Print_Matrix(type_real_vector,wt,title='final weights')
        END IF
        q_diff = abs ( q - q_diff)
        max_diff=0.d0
        DO j=1,nord
           max_diff = max(max_diff,q_diff(j))
        END DO
        write(iout,5) i, max_diff
        IF (max_diff <= tol) THEN
            EXIT
        END IF
        q_diff = q
        nq = nq + nq
     END DO
     DEALLOCATE(a,b,q_diff)
  END IF
!  
! Generate the needed functions at all required points.
!  
  CALL cpoly(p,dp,ddp,q,nord-1,nord,prn(5))
!  
! The DVR library assumes that the polynomials are $\delta$
! functions at the quadrature points.  Convert to this normalization
!
!
1    FORMAT(/,20X,'Computing DVR Points and Weights')
2    FORMAT(/,1X,'Standard Weight Function = ',a20)
3    FORMAT(/,1X,'Non-Standard Weight Function = ',a20,1x,'Reference Weight Function= ',a20)
4    FORMAT(/,1X,'weight integral = ',e15.8)
5    FORMAT(/,1X,'Iteration = ',i2,1x,'Maximum_Difference = ',e15.8)
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
     call Print_Matrix(type_real_matrix,p,ny,nx,title='polynomials')
     IF(present(drctv)) then
        call Print_Matrix(type_real_matrix,dp,ny,nx,title='derivative of polynomials')
        call Print_Matrix(type_real_matrix,ddp,ny,nx,title='second derivative of polynomials')
     END IF
  END IF
  END SUBROUTINE Lgngr
!***********************************************************************
!***********************************************************************
!deck cnvtpt.f
  SUBROUTINE cnvtpt(pt,wt,endpts,n)
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
!***********************************************************************
!***********************************************************************
!deck genq.f
  SUBROUTINE genq( a, b, w, z, fix, edge, muzero, n)
  IMPLICIT NONE
!        a, b     on input the recursion coefficients; will be destroyed
!     output parameters (both arrays of length n)
!        a        the desired points
!        w        the desired weights
  INTEGER                                :: n
  REAL(idp), DIMENSION(n)                :: a
  REAL(idp), DIMENSION(n)                :: b
  REAL(idp), DIMENSION(n)                :: w
  REAL(idp), DIMENSION(n,n)              :: z
  LOGICAL, DIMENSION(2)                  :: fix
  REAL(idp), DIMENSION(2)                :: edge
  REAL(idp)                              :: muzero
  REAL(idp)                              :: sumwt
  INTEGER                                :: i
  INTEGER                                :: ierr
!     the method used is a ql-type method with origin shifting
  z(:,1)=b
  b(2:n)=z(1:n-1,1)
  w=0.d0
  z=0.d0
  DO  i=1,n
      z(i,i)=1.d0
  END DO
  CALL imtql2 (n, n, a, b, z, ierr)
  w=muzero*z(1,:)*z(1,:)
  sumwt=sum(w,1)
  IF(fix(1)) THEN
     a(1)=edge(1)
  END IF
  IF(fix(2)) THEN
     a(n)=edge(2)
  END IF
  WRITE(iout,1) sumwt
1    FORMAT(/,1X,'sum of the weights = ',e15.8)
END SUBROUTINE genq
!***********************************************************************
!***********************************************************************
!deck genrwt.f
!***begin prologue     genrwt
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           reference weight
!***author             schneider, barry (nsf)
!***source
!***purpose            calculate weight functions and derivatives.
!***
!***description
!***
!***
!***references
!***routines called
!***end prologue       genrwt
  SUBROUTINE genrwt(rwt,pt,wt,n,refwt,alpha,beta)
  IMPLICIT NONE
  INTEGER                                :: n
  REAL(idp), DIMENSION(:)                :: rwt
  REAL(idp), DIMENSION(:)                :: pt
  CHARACTER (LEN=*)                      :: wt
  CHARACTER (LEN=*)                      :: refwt
  REAL(idp), OPTIONAL                    :: alpha
  REAL(idp), OPTIONAL                    :: beta
  IF (refwt == 'legendre' .or. refwt == 'one') Then
      rwt(:) =  1.d0
  ELSE IF(refwt == 'laguerre') THEN
      rwt(:) = exp( - pt(:) )
  ELSE IF(refwt == 'hermite') THEN
      rwt(:) = exp( - pt(:) * pt(:) )
  ELSE
      call lnkerr('error in reference weight function')
  END IF
!     Weight functions and their derivatives are computed for a variety of
!     cases
  IF(wt == 'r') THEN
     rwt=pt/rwt
  ELSE IF(wt == 'rr') THEN
     rwt=pt*pt/rwt 
  ELSE IF(wt == 'hermite') THEN
    rwt=EXP(-pt*pt)/rwt
  ELSE IF(wt == 'spherical_hermite') THEN
    rwt=EXP(-pt*pt)/rwt
  ELSE IF(wt == 'chebyshev-1') THEN
    rwt= SQRT ( 1.d0/( 1.d0 - pt*pt ) )/rwt
  ELSE IF(wt == 'chebyshev-2') THEN
    rwt= SQRT (  1.d0 -pt*pt )/rwt
  ELSE IF(wt == 'laguerre') THEN
    rwt= pt**alpha*EXP(-pt)/rwt
  ELSE IF(wt == 'jacobi') THEN
    rwt= (1.d0-pt)**alpha * (1.d0+pt)**beta/rwt
  ELSE IF(wt == 'rys') THEN
    rwt=EXP(-alpha*pt*pt)/rwt
  END IF
END SUBROUTINE genrwt
!***********************************************************************
!***********************************************************************
!deck lancz.f
!***begin prologue     lancz
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           eigenvalues, eigenvectors
!***author             schneider, barry (nsf)
!***source
!***purpose            generate orthogonal polynomials using recursion.
!***                   using the generated alpha and beta coefficients find
!***                   the points and weights of the generalized gauss
!***                   quadrature.
!***
!***description
!***references
!***routines called
!***end prologue       lancz
!\begin{eqnarray}
!  \beta_{j} P_{j}(x) &=& \big ( x - \alpha_{j} \big ) P_{j-1}(x) - \beta_{j-1} P_{j-2}(x)
!                             \\ \nonumber
!  \alpha_{j} &=& \langle P_{j-1} \mid x \mid P_{j-1} \rangle \\ \nonumber
!  \beta_{j} &=& \langle P_{j-1} \mid x \mid P_{j} \rangle \\ \nonumber
!  \beta_{j} P_{j}^{\prime}(x) &=& \big ( x - \alpha_{j} \big ) P_{j-1}^{\prime}(x)
!                             - \beta_{j-1} P_{j-2}^{\prime}(x) + P_{j-1}(x)
!                                 \\ \nonumber
!  \beta_{j} P_{j}^{\prime \prime}(x) &=& \big ( x - \alpha_{j} \big )
!                                      P_{j-1}^{\prime \prime}(x)
!                            - \beta_{j-1} P_{j-2}^{\prime \prime}(x)
!                            + 2 P_{j-1}^{\prime}(x)
!\end{eqnarray}
  SUBROUTINE lancz(r,a,b,wt,refwt,n,iter)
  IMPLICIT NONE
  INTEGER                                :: n
  INTEGER                                :: iter
  REAL(idp), DIMENSION(:,:), ALLOCATABLE :: v
  REAL(idp), DIMENSION(:)                :: r
  REAL(idp), DIMENSION(:)                :: wt
  REAL(idp), DIMENSION(:)                :: refwt
  REAL(idp), DIMENSION(:)                :: a
  REAL(idp), DIMENSION(:)                :: b
  REAL(idp)                              :: anorm
  INTEGER                                :: i
  ALLOCATE(v(n,0:iter))
  v(:,0) = 1.d0
  anorm=SQRT(1.d0/lanczos_dot(n,v(:,0),v(:,0),wt,refwt=refwt))
  v(:,0)=anorm*v(:,0)
!     we now have the first function
  IF (iter > 1) THEN
!         form argument times the first function
!         calculate a(1)
      a(1)=lanczos_dot(n,v(:,0),v(:,0),wt,refwt,r)
!         form mat times the first function - a(1) times the first function
!         and store it in the next polynomial
      v(:,1) = ( r(:) - a(1) ) * v(:,0)
!         calculate b(1)
      b(1)=SQRT( lanczos_dot(n,v(:,1),v(:,1),wt,refwt=refwt) )
!         normalize the second polynomial
      v(:,1)=(1.d0/b(1))*v(:,1)
      DO  i=2,iter
!            multiply the last calculated polynomial by mat
!            orthogonalize it to the two previous polynomials
!            calculating a(i) as we go
          a(i)=lanczos_dot(n,v(:,i-1),v(:,i-1),wt,refwt,r)
          v(:,i)= ( r(:) - a(i) ) * v(:,i-1) - b(i-1)*v(:,i-2)
!            calculate b(i)
          b(i)=SQRT( lanczos_dot(n,v(:,i),v(:,i),wt,refwt=refwt) )
!            normalize the polynomial and we are done
          v(:,i)=(1.d0/b(i))*v(:,i)
      END DO
  END IF
  DEALLOCATE(v)
END SUBROUTINE lancz
!***********************************************************************
!**********************************************************************
!deck modab.f
  SUBROUTINE modab(a,b,kpts,endpts,n)
  IMPLICIT NONE
  INTEGER                                :: n
  REAL(idp), DIMENSION(n)                :: a
  REAL(idp), DIMENSION(n)                :: b
  INTEGER                                :: kpts
  REAL(idp), DIMENSION(2)                :: endpts
  REAL(idp)                              :: gbslve
  REAL(idp)                              :: gam
  REAL(idp)                              :: t1
  IF (kpts == 0) THEN
      RETURN
  ELSE IF (kpts == 1) THEN
!         only a(n) must be changed
      a(n) =gbslve(endpts(1), n, a, b) * b(n-1)**2 + endpts(1)
      RETURN
  ELSE IF(kpts == 2) THEN
!         a(n) and b(n-1) must be recomputed
      gam =gbslve(endpts(1), n, a, b)
      t1 = ((endpts(1) - endpts(2))/(gbslve(endpts(2), n, a, b) - gam))
      b(n-1) =  SQRT(t1)
      a(n) = endpts(1) + gam*t1
      RETURN
  END IF
END SUBROUTINE modab
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
!     CALL Print_Matrix(type_real_matrix,cp(1,0),npt,n+1,npt,n+1,iout)
     CALL Print_Matrix(type_real_matrix,cp,npt,n+1,title='coordinate function')
     CALL Print_Matrix(type_real_matrix,dcp,npt,n+1,title='first derivative of coordinate function')
     CALL Print_Matrix(type_real_matrix,ddcp,npt,n+1,title='second derivative of coordinate function')
   END IF
END SUBROUTINE cpoly
!***********************************************************************
!***********************************************************************
!deck lanczos_dot.f
!***begin prologue     lanczos_dot
!***date written       090117   (yymmdd)
!***revision date      yymmdd   (yymmdd)                                       
!***keywords           
!***author             schneider, barry (nsf)
!***source
!***purpose            special dot product needed for lanczos routine.
!***
!***
!***
!***
!***description
!***references
!***routines called
!***end prologue       lanczos_dot
  REAL(idp) FUNCTION lanczos_dot(n,v_a,v_b,wt,pt,refwt)
  IMPLICIT NONE
  INTEGER                                :: n
  REAL(idp), DIMENSION(:)                :: v_a
  REAL(idp), DIMENSION(:)                :: v_b
  REAL(idp), DIMENSION(:), OPTIONAL      :: pt
  REAL(idp), DIMENSION(:)                :: wt
  REAL(idp), DIMENSION(:)                :: refwt
  INTEGER                                :: i
  lanczos_dot = 0.d0
  IF (present(pt)) THEN
      DO i=1,n
         lanczos_dot = lanczos_dot + v_a(i) * wt(i) * pt(i) * refwt(i) * v_b(i)
      END DO
  ELSE
      DO i=1,n
         lanczos_dot = lanczos_dot + v_a(i) * wt(i) * refwt(i) * v_b(i)
      END DO
  END IF
END FUNCTION lanczos_dot
!***********************************************************************
!***********************************************************************
           END MODULE DVR_Polynomials_Module
!***********************************************************************
!***********************************************************************
