!
MODULE Two_Electron_Shared
!***begin prologue     Two_Electron_Shared
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           time, finite element dvr, orthogonal polynomial
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            Global shared variables for two electron
!***description        integral code
!***                   
!
!***references

!***routines called    
!***end prologue       Two_Electron_Shared
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
  INTEGER                                 ::          number_of_radial_points
  INTEGER                                 ::          number_of_angular_points
  INTEGER                                 ::          number_of_eta_points
  INTEGER                                 ::          number_of_xi_points
  INTEGER                                 ::          n_tri
  INTEGER                                 ::          maximum_orbital_l
  INTEGER                                 ::          maximum_orbital_m
  INTEGER                                 ::          minimum_orbital_m
  INTEGER                                 ::          maximum_total_L
  INTEGER                                 ::          maximum_total_M
  INTEGER                                 ::          minimum_total_M
  CHARACTER(LEN=8)                        ::          Key

!
!
!---------------------------------------------------------------------
!                    Keywords for Print Options
!---------------------------------------------------------------------
!
!
!
END MODULE Two_Electron_Shared
