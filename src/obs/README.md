# Observation Operator

  This directory is the initiation of a set of observation operators for each new
  data 'class'. I am assuming that it will be most convenient and efficient to 
  process all similar data type (e.g. profiles, including temperature and salinity
  from Argo, XBTs, CTDs, etc.) in its own observation operator program.

  Ultimately, this should allow distributed work effort on each new data type,
  possibily even allowing the satellite data providers to compute their own
  observation operators for rapid inclusion into the operational ocean data
  assimilation system.

  Future:
  I'd like to make this more object oriented (as much as fortran can handle), e.g.
  with inheritance for an observation class down to specific instrument types.


# Add to obsop_sst, obsop_adt:
  Subroutines to read in raw data, e.g. read_aviso_adt.f90, read_avhrr_sst.f90
  Looking for L2 AVHRR and Microwave

# sample run lines:
./obsop.TEST5.adt.x -obsin 20070102adt.nc -gues gues -obsout obsout.dat -day 20820 -debug .true. -superob .true.