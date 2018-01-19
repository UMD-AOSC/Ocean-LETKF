DSET ^%y4%m2%d2.000000.ocean_temp_salt.res.nc
*DTYPE netcdf
TITLE Ocean Temperature/Salinity restart
OPTIONS template
*UNDEF 99999.0

XDEF xaxis_1
YDEF yaxis_1
ZDEF zaxis_1

TDEF Time 1000 linear 0z2jan1980 1dy

@ xaxis_1 String units degree_east 
@ xaxis_1 String long_name Nominal Longitude of T-cell center
@ yaxis_1 String units degree_north
@ yaxis_1 String long_name Nominal Latitude of T-cell center
@ zaxis_1 String units meters
@ zaxis_1 String long_name tz
@ Time String units days since 0001-01-01 00:00:00
@ Time String long_name Time

@ temp String long_name initial potential temp
@ temp String units deg_c
@ salt String long_name initial salinity
@ salt String units psu

