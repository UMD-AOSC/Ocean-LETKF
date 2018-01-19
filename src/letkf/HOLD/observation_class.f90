MODULE observation_class

USE common

TYPE observation

  REAL(r_size) :: elem
  REAL(r_size) :: lon
  REAL(r_size) :: lat
  REAL(r_size) :: lev
  REAL(r_size) :: time
  REAL(r_size) :: dat
  REAL(r_size) :: err
  REAL(r_size) :: hour
  REAL(r_size) :: qc

  REAL(r_size) :: 

CONTAINS

  PROCEDURE :: initialize
  PROCEDURE :: finalize

END TYPE observation

TYPE, EXTENDS (observation) :: observation_operator

  REAL(r_size) :: hx

END TYPE hx


TYPE, EXTENDS (observation_operator) :: innovation

  REAL(r_size) :: omf
  REAL(r_size) :: sigma_radius0
  REAL(r_size) :: sigma_radius
  REAL(r_size) :: sigma_radiusv

END TYPE hx

INTERFACE 
  SUBROUTINE initialize(yo)

    CLASS(observation), ALLOCATABLE, DIMENSION(:) :: yo


  END SUBROUTINE initialize 

END INTERFACE


END MODULE observation_class
