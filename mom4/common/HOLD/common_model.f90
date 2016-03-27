MODULE common_model
!STEVE: The goal of this module is to generalize the common_<model>
!       modules to a single module that can switch between desired models.

INTEGER, PARAMETER :: hycom = 1
INTEGER, PARAMETER :: roms = 2
INTEGER, PARAMETER :: sis  = 3
INTEGER, PARAMETER :: mom4 = 4
INTEGER, PARAMETER :: mom5 = 5
INTEGER, PARAMETER :: mom6 = 6
INTEGER, PARAMETER :: model_type = mom4

CONTAINS

SUBROUTINE set_common_model()
USE common_mom4, ONLY: set_common_mom4

SELECT(model_type)
CASE(mom4)
  CALL set_common_mom4
DEFAULT
  WRITE(6,*) "model_type not available :: ", model_type
  STOP(75)
END SELECT

SUBROUTINE set_common_model

SUBROUTINE read_diag
END SUBROUTINE read_diag

SUBROUTINE read_restart
END SUBROUTINE read_restart

SUBROUTINE write_restart
END SUBROUTINE write_restart

END MODULE common_model
