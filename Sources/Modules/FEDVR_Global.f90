!
MODULE FEDVR_Global
!***begin prologue     FEDVR_Global
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
!***end prologue       FEDVR_Global
  USE accuracy
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!---------------------------------------------------------------------
!                    shared variables used in lots of routines
!---------------------------------------------------------------------
  CHARACTER (LEN=80), DIMENSION(2)        ::          vtyp
  CHARACTER (LEN=80), DIMENSION(30)       ::          coord
  CHARACTER (LEN=128)                     ::          filbec
  CHARACTER (LEN=8)                       ::          typke
  CHARACTER (LEN=24)                      ::          kinetic_energy_type
  CHARACTER (LEN=24)                      ::          type_calculation
  CHARACTER (LEN=80)                      ::          title
  CHARACTER (LEN=128)                     ::          matrix_title
  CHARACTER (LEN=80)                      ::          keyword
  CHARACTER (LEN=80)                      ::          units
  CHARACTER (LEN=80)                      ::          typpot
  CHARACTER (LEN=80)                      ::          typwt
  CHARACTER (LEN=80)                      ::          typarg
  CHARACTER (LEN=80)                      ::          dentyp
  CHARACTER (LEN=80)                      ::          typint
  CHARACTER (LEN=80)                      ::          refwt
  CHARACTER (LEN=80)                      ::          cpass
  CHARACTER (LEN=80)                      ::          drctv
  CHARACTER (LEN=2400)                    ::          card
  CHARACTER (LEN=4)                       ::          parity
  CHARACTER (LEN=16)                      ::          range
  INTEGER                                 ::          physical_points
  INTEGER                                 ::          global_points
  INTEGER, DIMENSION(30)                  ::          nphy
  INTEGER, DIMENSION(30)                  ::          nglobal
  INTEGER, DIMENSION(30)                  ::          row
  INTEGER, DIMENSION(30)                  ::          nonz
  INTEGER, PARAMETER                      ::          maxreg=5000
  INTEGER, DIMENSION(maxreg)              ::          n
  INTEGER, DIMENSION(maxreg)              ::          npt
  INTEGER, DIMENSION(maxreg)              ::          tri_reg
  INTEGER, DIMENSION(maxreg)              ::          nrq
  INTEGER                                 ::          spdim
  INTEGER                                 ::          nfix
  INTEGER                                 ::          bcl
  INTEGER                                 ::          bcr
  INTEGER                                 ::          nreg
  INTEGER                                 ::          l_orb_max 
  INTEGER                                 ::          l_max 
  INTEGER                                 ::          m_max 
  INTEGER                                 ::          l_val 
  INTEGER                                 ::          m_val 
  INTEGER                                 ::          nphy_tri 
  INTEGER                                 ::          angmom
  INTEGER                                 ::          size
  INTEGER, DIMENSION(:,:), ALLOCATABLE    ::          nfun_reg
  INTEGER, DIMENSION(:), ALLOCATABLE      ::          num_reg
  LOGICAL                                 ::          diag
  LOGICAL                                 ::          poisson
  LOGICAL                                 ::          two_electron
  LOGICAL                                 ::          nlse
  LOGICAL                                 ::          genpts
  LOGICAL                                 ::          reuse
  LOGICAL                                 ::          nodiag
  LOGICAL                                 ::          reuse_sector_information
  LOGICAL                                 ::          atomic
  LOGICAL                                 ::          proj
  LOGICAL, DIMENSION(2)                   ::          drop
  LOGICAL, DIMENSION(2)                   ::          fix
  LOGICAL                                 ::          skip
  REAL*4,  DIMENSION(20)                  ::          del
  REAL*4,  DIMENSION(20)                  ::          tim
  REAL(idp),  DIMENSION(2)                ::          endpts
  REAL(idp)                               ::          mass
  REAL(idp),  DIMENSION(maxreg+1)         ::          edge
  REAL(idp)                               ::          box
  REAL(idp)                               ::          deltax
  REAL(idp)                               ::          dscale
  REAL(idp)                               ::          Z_a
  REAL(idp)                               ::          Z_b
  REAL(idp)                               ::          R_ab
  REAL(idp)                               ::          pre_factor
  REAL(idp)                               ::          factor
  REAL(idp)                               ::          last_value
  LOGICAL                                 ::          toau
  LOGICAL                                 ::          useau
  LOGICAL                                 ::          prnt
!
!
!---------------------------------------------------------------------
!                    Keywords for Print Options
!---------------------------------------------------------------------
!
  CHARACTER (LEN=80), DIMENSION(12) :: prnkey
  CHARACTER (LEN=80)                :: secprn
  DATA prnkey /'sector_points', 'sector_factors',                         &
               'sector_polynomials', 'sector_matrices',                   &
               'global_points', 'global_polynomials',                     &
               'potential','global_matrices',                             &
               'hamiltonian', 'eigenvalues','eigenvectors', 'all' /
  DATA secprn / 'none' /
!
!
 CHARACTER (LEN=80),DIMENSION(12)         :: prloc
 LOGICAL, DIMENSION(12)                   :: prn
!
!
!---------------------------------------------------------------------
!                    Derived types and Allocated variables for 
!                             the dvr quantities
!---------------------------------------------------------------------
  TYPE dvr_matrices
    REAL(idp), DIMENSION(:,:), ALLOCATABLE       :: ham
    REAL(idp), DIMENSION(:,:), ALLOCATABLE       :: eigenvectors
    REAL(idp), DIMENSION(:),   ALLOCATABLE       :: eigenvalues
    REAL(idp), DIMENSION(:),   ALLOCATABLE       :: work
    REAL(idp), DIMENSION(:,:), ALLOCATABLE       :: rhs
    REAL(idp), DIMENSION(:),   ALLOCATABLE       :: points
    REAL(idp), DIMENSION(:),   ALLOCATABLE       :: upper
    INTEGER,   DIMENSION(:),   ALLOCATABLE       :: ipvt
    INTEGER                                      :: number_of_right_hand_sides
    CHARACTER(LEN=80)                            :: type_inhomo
  END TYPE dvr_matrices

  TYPE pt_wt
       REAL(idp), DIMENSION(:), ALLOCATABLE      :: qr
       REAL(idp), DIMENSION(:), ALLOCATABLE      :: wtr
       REAL(idp), DIMENSION(:), ALLOCATABLE      :: qr_fac
       REAL(idp), DIMENSION(:), ALLOCATABLE      :: inv_sqrt_qr_fac
       REAL(idp), DIMENSION(:), ALLOCATABLE      :: inv_qr_fac
  END TYPE pt_wt

  TYPE functions
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: pr
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: dpr
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: ddpr
       REAL(idp), DIMENSION(:),   ALLOCATABLE    :: normalization
  END TYPE functions

  TYPE odd_functions
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: pr
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: dpr
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: ddpr
  END TYPE odd_functions

  TYPE fourier_functions
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: pr
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: dpr
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: ddpr
  END TYPE fourier_functions

  TYPE matrices
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: aa
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: tr
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: ham
  END TYPE matrices

  TYPE fourier_matrices
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: tr
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: ham
  END TYPE fourier_matrices

  TYPE laguerre_matrices
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: tr
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: ham
  END TYPE laguerre_matrices

  TYPE hermite_matrices
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: tr
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: ham
  END TYPE hermite_matrices

  TYPE odd_matrices
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: tr
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: ham
  END TYPE odd_matrices

  TYPE hamiltonian
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: tr
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: ham
  END TYPE hamiltonian

  TYPE vectors
       REAL(idp), DIMENSION(:),   ALLOCATABLE    :: vr
       REAL(idp), DIMENSION(:),   ALLOCATABLE    :: s
  END TYPE vectors

  TYPE well
       REAL(idp)                                 :: depth
       INTEGER                                   :: nwell
       REAL(idp)                                 :: awell
  END TYPE well

  TYPE exponential
       REAL(idp), DIMENSION(2)                   :: expnt
       REAL(idp), DIMENSION(2)                   :: amp
  END TYPE exponential

  TYPE sum_exponential
       REAL(idp), DIMENSION(2)                   :: expnt
       REAL(idp), DIMENSION(2)                   :: amp
  END TYPE sum_exponential

  TYPE power_exponential
       REAL(idp), DIMENSION(2)                   :: expnt
       REAL(idp), DIMENSION(2)                   :: amp
       INTEGER                                   :: n_p
  END TYPE power_exponential

  TYPE yukawa
       REAL(idp), DIMENSION(2)                   :: expnt
       REAL(idp), DIMENSION(2)                   :: amp
  END TYPE yukawa

  TYPE coulomb
       REAL(idp)                                 :: charge
  END TYPE coulomb

  TYPE eberlonium
       REAL(idp)                                 :: charge
       REAL(idp), DIMENSION(2)                   :: amp
       INTEGER                                   :: n_p
  END TYPE eberlonium

  TYPE harmonic_oscillator
       REAL(idp)                                 :: mass
       REAL(idp)                                 :: omega
  END TYPE harmonic_oscillator

  TYPE expres
       REAL(idp), DIMENSION(2)                   :: expnt
       REAL(idp), DIMENSION(2)                   :: amp
       REAL(idp)                                 :: shift
  END TYPE expres

  TYPE periodic
       INTEGER                                   :: n_scale
       REAL(idp)                                 :: e_c
       REAL(idp)                                 :: awell
  END TYPE periodic

  TYPE inverse_r4
       REAL(idp)                                 :: dummy
  END TYPE inverse_r4 

  TYPE potential
       TYPE (well)                               :: well
       TYPE (yukawa)                             :: yukawa
       TYPE (exponential)                        :: exponential
       TYPE (coulomb)                            :: coulomb
       TYPE (expres)                             :: expres
       TYPE (harmonic_oscillator)                :: harmonic_oscillator
       TYPE (eberlonium)                         :: eberlonium
       TYPE (power_exponential)                  :: power_exp
       TYPE (sum_exponential)                    :: sum_exp
       TYPE (inverse_r4)                         :: inv_r4
       TYPE (periodic)                           :: periodic
       CHARACTER(LEN=80)                         :: type
       REAL(idp), DIMENSION(:),   ALLOCATABLE    :: vec
  END TYPE potential

  TYPE coordinates
       CHARACTER(LEN=24)                         :: label
       INTEGER                                   :: num_reg
       INTEGER, DIMENSION(:),     ALLOCATABLE    :: num_pts_reg
       LOGICAL, DIMENSION(2)                     :: drop_pt
       LOGICAL, DIMENSION(2)                     :: fix_pt
       INTEGER                                   :: num_fixed
       REAL(idp), DIMENSION(:),   ALLOCATABLE    :: grid_points
       REAL(idp), DIMENSION(:),   ALLOCATABLE    :: grid_weights
       TYPE(pt_wt),               DIMENSION(:),                                 &
                                  ALLOCATABLE    :: reg_pt_wt
       TYPE(functions),           DIMENSION(:),                                 &
                                  ALLOCATABLE    :: reg_poly
       TYPE(odd_functions),       DIMENSION(:),                                 &
                                  ALLOCATABLE    :: reg_poly_odd
       TYPE(fourier_functions),   DIMENSION(:),                                 &
                                  ALLOCATABLE    :: reg_poly_fourier
       TYPE(matrices),            DIMENSION(:),                                 &
                                  ALLOCATABLE    :: reg_mat
       TYPE(odd_matrices),        DIMENSION(:),                                 &
                                  ALLOCATABLE    :: reg_mat_odd
       TYPE(fourier_matrices),    DIMENSION(:),                                 &
                                  ALLOCATABLE    :: reg_mat_fourier
       TYPE(hermite_matrices),    DIMENSION(:),                                 &
                                  ALLOCATABLE    :: reg_mat_hermite
       TYPE(laguerre_matrices),   DIMENSION(:),                                 &
                                  ALLOCATABLE    :: reg_mat_laguerre
       TYPE(vectors),             DIMENSION(:),                                 &
                                  ALLOCATABLE    :: reg_vec
       TYPE(hamiltonian),         DIMENSION(:,:),                               &
                                  ALLOCATABLE    :: reg_ham
  END TYPE coordinates

  TYPE (coordinates),           DIMENSION(:),                                   &
                                ALLOCATABLE      :: reg_grid

  TYPE (potential),             DIMENSION(:),                                   &
                                ALLOCATABLE      :: reg_pot

  TYPE(dvr_matrices),           DIMENSION(:),                                   &
                                ALLOCATABLE      :: dvr_mat

  TYPE buffer
    REAL(idp),    DIMENSION(:), ALLOCATABLE                                     &
                                                 :: d, hbuf
    REAL(idp),    DIMENSION(:,:), ALLOCATABLE                                   &
                                                 :: hibuf
  END TYPE buffer

  TYPE(buffer), DIMENSION(:),                                                   &
                ALLOCATABLE                      :: buf    
  REAL(idp), DIMENSION(:), ALLOCATABLE           :: v_pert
!
!
END MODULE FEDVR_Global
