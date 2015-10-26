!***********************************************************************
! Utility_Module
!**begin prologue     Utility_Module
!**date written       090219   (yymmdd)
!**revision date               (yymmdd)
!**keywords           Utility
!**
!**author             schneider, b. i.(nsf)
!**source             
!**purpose            Utility routines for FEDVR Library
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***                  
!***                  
!***end prologue      Utility_Module
!***********************************************************************
!***********************************************************************
                           MODULE Utility_Module
              USE FEDVR_Shared
              USE Matrix_Print
              TYPE(REAL_MATRIX)              :: type_real_matrix
              TYPE(REAL_VECTOR)              :: type_real_vector
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!*deck lgngr
      SUBROUTINE LGNGR(p,dp,ddp,x,y,nx,ny,prnt,drctv,type) 
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
  REAL*8, DIMENSION (:,:)    :: p
  REAL*8, DIMENSION (:,:)    :: dp
  REAL*8, DIMENSION (:,:)    :: ddp
  REAL*8, DIMENSION(:)       :: x
  REAL*8, DIMENSION(:)       :: y
  REAL*8                     :: sn
  REAL*8                     :: ssn
  REAL*8                     :: fac
  LOGICAL                    :: prnt
  CHARACTER (LEN = 80)       :: title 
  CHARACTER (LEN = *)        :: drctv
  CHARACTER (LEN = *)        :: type
  INTEGER                    :: nx
  INTEGER                    :: ny
  INTEGER                    :: i
  INTEGER                    :: j
  INTEGER                    :: k
  INTEGER                    :: first
  INTEGER                    :: second
  INTEGER                    :: zerfac
  INTEGER                    :: inp
  INTEGER                    :: iout
  common /io/ inp, iout

!
!     generate polynomials and derivatives with respect to x
!
  p(:,:) = 1.d0
  IF (type /= 'even') THEN 
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
            IF( drctv /= 'functions-only') THEN
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
         END DO
!
      END DO
  ELSE 
      DO i=1,ny
         zerfac = 0
         DO j=1,nx
            fac =  y(i)*y(i) - x(j)*x(j) 
            IF(abs(fac) <= 1.d-10) THEN
               zerfac = j
            ENDIF  
         END DO
         DO j=1,nx
            DO k=1,j-1
               p(i,j) = p(i,j) * ( y(i)*y(i) - x(k)*x(k) )   &
                            / ( x(j)*x(j) - x(k)*x(k) )
            END DO
            DO k=j+1,nx
               p(i,j) = p(i,j) * ( y(i)*y(i) - x(k)*x(k) )  &
                            / ( x(j)*x(j) - x(k)*x(k) )
            END DO
            IF( drctv /= 'functions-only') THEN
                IF ( abs(p(i,j)) > 1.d-10) THEN
                     sn = 0.d0
                     ssn = 0.d0
                     DO k=1,j-1
                        fac = 1.d0/( y(i)*y(i) - x(k)*x(k) )
                        sn = sn + fac
                        ssn = ssn + fac*fac
                     END DO
                     DO k=j+1,nx
                        fac = 1.d0/( y(i)*y(i) - x(k)*x(k) )
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
                        fac = 1.d0/( x(j)*x(j) - x(k)*x(k) )
                        sn = sn*fac*( y(i)*y(i) - x(k)*x(k) )
                        ssn = ssn + 1.d0/(y(i)*y(i) - x(k)*x(k))
                     END DO
                     DO k=first+1,second-1
                        fac = 1.d0/( x(j)*x(j) - x(k)*x(k) )
                        sn = sn*fac*( y(i)*y(i) - x(k)*x(k) )
                        ssn = ssn + 1.d0/( y(i)*y(i) - x(k)*x(k) )             
                     END DO
                     DO k=second+1,nx
                        fac = 1.d0/( x(j)*x(j) - x(k)*x(k) )
                        sn = sn*fac*( y(i)*y(i) - x(k)*x(k) )
                        ssn = ssn + 1.d0/( y(i)*y(i) - x(k)*x(k) )             
                     END DO
                     dp(i,j) = sn/( x(j)*x(j) - x(zerfac)*x(zerfac) )
                     ddp(i,j) = 2.d0*ssn*dp(i,j)
                     ddp(i,j) = 2.d0*dp(i,j) + 4.d0 * y(i) * y(i) * ddp(i,j) 
                     dp(i,j) = 2.d0 * y(i) * dp(i,j)
                END IF                    
            END IF
         END DO
      END DO
!
  END IF
  IF(prnt) THEN
     Call Print_Matrix(type_real_matrix,p,ny,nx,title='polynomials')
     IF(drctv.ne.'functions-only') then
        title='derivative of polynomials'
        Call Print_Matrix(type_real_matrix,dp,ny,nx,title='first derivative of polynomials')
        Call Print_Matrix(type_real_matrix,ddp,ny,nx,title='second derivative of polynomials')
    END IF
  END IF
  END SUBROUTINE Lgngr
!***********************************************************************
!***********************************************************************
!deck wddtot.f
!***begin prologue     wddtot
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           eigenvalues, eigenvectors
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***                   
!***                   
!***                   
!***
!***description
!***references
!***routines called
!***end prologue       wddtot
  FUNCTION WDDOT(n,v_a,v_b,wt,pt)
  USE accuracy
  IMPLICIT NONE
  INTEGER                                :: n
  REAL(idp), DIMENSION(:)                :: v_a
  REAL(idp), DIMENSION(:)                :: v_b
  REAL(idp), DIMENSION(:), OPTIONAL      :: pt
  REAL(idp), DIMENSION(:)                :: wt
  REAL(idp)                              :: WDDOT
  INTEGER                                :: i
  WDDOT = 0.D0
  IF (present(pt)) THEN
      DO i=1,n
         WDDOT = WDDOT + v_a(i) * wt(i) * pt(i) * v_b(i)
      END DO 
  ELSE
      DO i=1,n
         WDDOT = WDDOT + v_a(i) * wt(i) * v_b(i)
      END DO 
  END IF
END FUNCTION WDDOT
!***********************************************************************
!***********************************************************************
           END MODULE Utility_Module
!***********************************************************************
!***********************************************************************
