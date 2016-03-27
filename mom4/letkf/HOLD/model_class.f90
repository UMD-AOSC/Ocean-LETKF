MODULE observation_class

USE common

TYPE model

  TYPE(model_grid)       :: x
  TYPE(model_parameters) :: p

CONTAINS

  PROCEDURE :: initialize
  PROCEDURE :: finalize

END TYPE model

TYPE model_grid

  REAL(r_size), DIMENSION(:,:) :: lon
  REAL(r_size), DIMENSION(:,:) :: lat
  REAL(r_size), DIMENSION(:) :: lev
  REAL(r_size) :: time
  INTEGER, DIMENSION(:,:) :: mask
  TYPE(vertical_coordinate) :: v_coord

END TYPE model_grid

TYPE model_parameters

  REAL(r_size)                   :: p0d
  REAL(r_size), DIMENSION(:)     :: p1d
  REAL(r_size), DIMENSION(:,:)   :: p2d
  REAL(r_size), DIMENSION(:,:,:) :: p3d

END TYPE model_parameters

TYPE vertical_coordinate

  CHARACTER(slen) i              :: vc
  REAL(r_size), DIMENSION(:)     :: height

END TYPE vertical_coordinate

TYPE, EXTENDS (model) :: innovation

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
