DSET ^%y4%m2%d2.000000.ocean_sbc.res.nc
*DTYPE netcdf
TITLE Ocean Surface Data restart
OPTIONS template
*UNDEF 99999.0

XDEF xaxis_1
YDEF yaxis_1

TDEF Time 1000 linear 0z2jan1980 1dy

@ xaxis_1 String units degree_east 
@ xaxis_1 String long_name Nominal Longitude of T-cell center
@ yaxis_1 String units degree_north
@ yaxis_1 String long_name Nominal Latitude of T-cell center
@ Time String units days since 0001-01-01 00:00:00
@ Time String long_name Time

@ t_surf String long_name Surface Temperature
@ t_surf String units K
@ s_surf String long_name Surface Salinity
@ s_surf String units psu
@ u_surf String long_name Surface Zonal Velocity
@ u_surf String units m/s
@ v_surf String long_name Surface Meridional Velocity
@ v_surf String units m/s
@ sea_lev String long_name Sea Level
@ sea_lev String units m
@ frazil String long_name Frazil
@ frazil String units units

