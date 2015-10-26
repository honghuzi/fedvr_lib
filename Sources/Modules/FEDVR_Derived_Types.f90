!
MODULE FEDVR_Derived_Types
!***begin prologue     FEDVR_Derived_Types
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           time, finite element dvr, orthogonal polynomial
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            global shared variables for dvr library
!***description        this routine defines the global variables
!***                   and data needed for the dvrlib
!
!***references

!***routines called    
!***end prologue       FEDVR_Derived_Types
  USE accuracy
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!---------------------------------------------------------------------
!                    Derived types and Allocated variables for 
!                             the dvr quantities
!---------------------------------------------------------------------
!

  TYPE spherical_harmonics
    CHARACTER(LEN=3)                             :: kind
  END TYPE spherical_harmonics

  TYPE spherical_fedvr
    CHARACTER(LEN=3)                             :: kind
  END TYPE spherical_fedvr

  TYPE cartesian
    CHARACTER(LEN=3)                             :: kind
  END TYPE cartesian

  TYPE spherical
    TYPE(spherical_harmonics)                    :: spherical_harmonics
    TYPE(spherical_fedvr)                        :: spherical_fedvr
    CHARACTER(LEN=3)                             :: kind
  END TYPE spherical

  TYPE spheroidal
    CHARACTER(LEN=3)                             :: kind
  END TYPE spheroidal

  TYPE Even
    Character(LEN=32)                            :: kind
  END TYPE Even  

  TYPE Odd
    Character(LEN=32)                            :: kind
  END TYPE Odd 

  TYPE Ham
    LOGICAL                                      :: true_false
  END TYPE Ham  

  TYPE Nabla
    LOGICAL                                      :: true_false
  END TYPE Nabla 

  TYPE dvr_matrices
    REAL(idp), DIMENSION(:,:),     ALLOCATABLE   :: ham
    REAL(idp), DIMENSION(:,:),     ALLOCATABLE   :: tr
    REAL(idp), DIMENSION(:,:),     ALLOCATABLE   :: Q
    REAL(idp), DIMENSION(:,:),     ALLOCATABLE   :: eigenvectors
    REAL(idp), DIMENSION(:),       ALLOCATABLE   :: eigenvalues
    REAL(idp), DIMENSION(:,:),     ALLOCATABLE   :: RHS
    REAL(idp), DIMENSION(:,:),     ALLOCATABLE   :: Inverse
    REAL(idp), DIMENSION(:,:),     ALLOCATABLE   :: Exact_Solution
    REAL(idp), DIMENSION(:),       ALLOCATABLE   :: work
    REAL(idp), DIMENSION(:),       ALLOCATABLE   :: points
    REAL(idp), DIMENSION(:),       ALLOCATABLE   :: lower
    INTEGER,   DIMENSION(:),       ALLOCATABLE   :: ipvt
  END TYPE dvr_matrices

  TYPE angular_matrices
    REAL(idp), DIMENSION(:),       ALLOCATABLE   :: D_LM_Coef
    INTEGER,   DIMENSION(:,:),     ALLOCATABLE   :: D_LM_Index
  END TYPE angular_matrices

  TYPE pt_wt
       REAL(idp), DIMENSION(:),    ALLOCATABLE   :: qr
       REAL(idp), DIMENSION(:),    ALLOCATABLE   :: wtr
       REAL(idp), DIMENSION(:),    ALLOCATABLE   :: qr_fac
       REAL(idp), DIMENSION(:),    ALLOCATABLE   :: inv_sqrt_qr_fac
       REAL(idp), DIMENSION(:),    ALLOCATABLE   :: inv_qr_fac
  END TYPE pt_wt

  TYPE functions
       REAL(idp), DIMENSION(:,:),  ALLOCATABLE   :: pr
       REAL(idp), DIMENSION(:,:),  ALLOCATABLE   :: dpr
       REAL(idp), DIMENSION(:,:),  ALLOCATABLE   :: ddpr
       REAL(idp), DIMENSION(:),    ALLOCATABLE   :: normalization
  END TYPE functions

  TYPE fourier_functions
       REAL(idp), DIMENSION(:,:),  ALLOCATABLE    :: pr
       REAL(idp), DIMENSION(:,:),  ALLOCATABLE    :: dpr
       REAL(idp), DIMENSION(:,:),  ALLOCATABLE    :: ddpr
       REAL(idp), DIMENSION(:),    ALLOCATABLE    :: normalization
  END TYPE fourier_functions

  TYPE odd_matrices
       REAL(idp), DIMENSION(:,:),  ALLOCATABLE    :: tr
       REAL(idp), DIMENSION(:,:),  ALLOCATABLE    :: ham
       CHARACTER (LEN=4)                          :: type
  END TYPE odd_matrices

  TYPE fourier_matrices
       REAL(idp), DIMENSION(:,:),  ALLOCATABLE    :: tr
       REAL(idp), DIMENSION(:,:),  ALLOCATABLE    :: ham
  END TYPE fourier_matrices

  TYPE laguerre_matrices
       REAL(idp), DIMENSION(:,:),  ALLOCATABLE    :: tr
       REAL(idp), DIMENSION(:,:),  ALLOCATABLE    :: ham
  END TYPE laguerre_matrices

  TYPE hermite_matrices
       REAL(idp), DIMENSION(:,:),  ALLOCATABLE    :: tr
       REAL(idp), DIMENSION(:,:),  ALLOCATABLE    :: ham
  END TYPE hermite_matrices

  TYPE matrices
       TYPE(Even)                                 :: type_even
       TYPE(Odd)                                  :: type_odd
       TYPE(Ham)                                  :: type_hamiltonian
       TYPE(Nabla)                                :: type_nabla
       REAL(idp), DIMENSION(:,:),  ALLOCATABLE    :: tr
       REAL(idp), DIMENSION(:,:),  ALLOCATABLE    :: ham
       CHARACTER (LEN=4)                          :: type
  END TYPE matrices

  TYPE vectors
       REAL(idp), DIMENSION(:),    ALLOCATABLE    :: vr
       REAL(idp), DIMENSION(:),    ALLOCATABLE    :: s
  END TYPE vectors

  TYPE well
       REAL(idp)                                  :: depth
       INTEGER                                    :: nwell
       REAL(idp)                                  :: awell
  END TYPE well

  TYPE exponential
       REAL(idp), DIMENSION(2)                    :: expnt
       REAL(idp), DIMENSION(2)                    :: amp
  END TYPE exponential

  TYPE sum_exponential
       REAL(idp), DIMENSION(2)                    :: expnt
       REAL(idp), DIMENSION(2)                    :: amp
  END TYPE sum_exponential

  TYPE power_exponential
       REAL(idp), DIMENSION(2)                    :: expnt
       REAL(idp), DIMENSION(2)                    :: amp
       INTEGER                                    :: n_p
  END TYPE power_exponential

  TYPE yukawa
       REAL(idp), DIMENSION(2)                    :: expnt
       REAL(idp), DIMENSION(2)                    :: amp
  END TYPE yukawa

  TYPE coulomb
       REAL(idp)                                  :: charge
  END TYPE coulomb

  TYPE eberlonium
       REAL(idp)                                  :: charge
       REAL(idp), DIMENSION(2)                    :: amp
       INTEGER                                    :: n_p
  END TYPE eberlonium

  TYPE harmonic_oscillator
       REAL(idp)                                  :: mass
       REAL(idp)                                  :: omega
  END TYPE harmonic_oscillator

  TYPE expres
       REAL(idp), DIMENSION(2)                    :: expnt
       REAL(idp), DIMENSION(2)                    :: amp
       REAL(idp)                                  :: shift
  END TYPE expres

  TYPE periodic
       INTEGER                                    :: n_scale
       REAL(idp)                                  :: e_c
       REAL(idp)                                  :: awell
  END TYPE periodic

  TYPE inverse_r4
       REAL(idp)                                  :: dummy
  END TYPE inverse_r4 

  TYPE potential
       TYPE (well)                                :: well
       TYPE (yukawa)                              :: yukawa
       TYPE (exponential)                         :: exponential
       TYPE (coulomb)                             :: coulomb
       TYPE (expres)                              :: expres
       TYPE (harmonic_oscillator)                 :: harmonic_oscillator
       TYPE (eberlonium)                          :: eberlonium
       TYPE (power_exponential)                   :: power_exp
       TYPE (sum_exponential)                     :: sum_exp
       TYPE (inverse_r4)                          :: inv_r4
       TYPE (periodic)                            :: periodic
       CHARACTER(LEN=80)                          :: type
       REAL(idp), DIMENSION(:),    ALLOCATABLE    :: vec
       REAL(idp), DIMENSION(0:2)                  :: c
  END TYPE potential

  TYPE coordinates
       CHARACTER(LEN=24)                          :: label
       INTEGER                                    :: num_reg
       INTEGER, DIMENSION(:),      ALLOCATABLE    :: num_pts_reg
       INTEGER                                    :: num_points
       LOGICAL, DIMENSION(2)                      :: drop_pt
       LOGICAL, DIMENSION(2)                      :: fix_pt
       INTEGER                                    :: num_fixed
       REAL(idp), DIMENSION(:),    ALLOCATABLE    :: grid_points
       REAL(idp), DIMENSION(:),    ALLOCATABLE    :: grid_weights
       TYPE(cartesian)                            :: name_cartesian
       TYPE(spherical)                            :: name_spherical
       TYPE(spheroidal)                           :: name_spheroidal
       TYPE(pt_wt),                DIMENSION(:),                                &
                                   ALLOCATABLE    :: reg_pt_wt
       TYPE(functions),            DIMENSION(:),                                &
                                   ALLOCATABLE    :: reg_poly
       TYPE(fourier_functions),    DIMENSION(:),                                &
                                   ALLOCATABLE    :: reg_poly_fourier
       TYPE(matrices),             DIMENSION(:),                                &
                                   ALLOCATABLE    :: reg_mat
       TYPE(matrices),             DIMENSION(:,:),                              &
                                   ALLOCATABLE    :: reg_type_op
       TYPE(odd_matrices),         DIMENSION(:),                                &
                                   ALLOCATABLE    :: reg_mat_odd
       TYPE(fourier_matrices),     DIMENSION(:),                                &
                                   ALLOCATABLE    :: reg_mat_fourier
       TYPE(hermite_matrices),     DIMENSION(:),                                &
                                   ALLOCATABLE    :: reg_mat_hermite
       TYPE(laguerre_matrices),    DIMENSION(:),                                &
                                   ALLOCATABLE    :: reg_mat_laguerre
       TYPE(vectors),              DIMENSION(:),                                &
                                   ALLOCATABLE    :: reg_vec
  END TYPE coordinates

  TYPE (coordinates),              DIMENSION(:),                                &
                                   ALLOCATABLE    :: reg_grid

  TYPE (potential),                DIMENSION(:),                                &
                                   ALLOCATABLE    :: reg_pot

  TYPE(dvr_matrices),              DIMENSION(:),                                &
                                   ALLOCATABLE    :: dvr_mat

  TYPE(angular_matrices),          DIMENSION(:,:),                              &
                                   ALLOCATABLE    :: ang_mat
  TYPE buffer
    REAL(idp),                     DIMENSION(:),                                &
                                   ALLOCATABLE    :: d, hbuf

    REAL(idp),                     DIMENSION(:,:),                              &
                                   ALLOCATABLE    :: hibuf
  END TYPE buffer

  TYPE(buffer),                    DIMENSION(:),                                &
                                   ALLOCATABLE    :: buf    
  REAL(idp),                       DIMENSION(:),                                &
                                   ALLOCATABLE    :: v_pert
!
!
END MODULE FEDVR_Derived_Types
