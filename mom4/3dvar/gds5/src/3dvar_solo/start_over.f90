PROGRAM start_over

! STEVE: start from scratch and try to create a comparable 3DVar program to godas


! Read in model state
!CALL read_grd(infile,v3d,v2d)

! Read in observation innovations (need to compute obsop.f90 on mean state prior to running this)
!CALL read_obs2

! Read in the background error covariance inputs

! Compute analysis
!CALL setup_pcg
!CALL pcg
!CALL finish_pcg


! Write analyzed model state
!CALL write_grd(outfile,v3d,v2d)


END PROGRAM start_over
