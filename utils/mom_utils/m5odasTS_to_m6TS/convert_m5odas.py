#!/usr/bin/env python3

"""
Usage:
./convert_m5odas.py ts_out.nc --m5Ana ocean_temp_salt.res.nc --m5Grid grid_spec.nc.mom5_LL1440x721 --m5Date 2015080900
"""

import numpy as np
from netCDF4 import Dataset
import os, sys, argparse
import datetime as dt

MISSING_VALUE = 9.96921e36

## need to be consistent with UMD_rc/input.nml: &ocean_tempsalt_nml
T_MIN_LIMIT = -5.0
T_MAX_LIMIT = 32.0
S_MIN_LIMIT =  2.0
S_MAX_LIMIT = 42.0

def read_m5ana(fnin="ocean.temp_salt.res.latlon.nc"):
    """
    read in salt & temp from lat-lon m5 odas analysis
    """
    f = Dataset(fnin,"r")
    f.set_auto_mask(False)
    temp =  f.variables['temp'][:] 
    salt =  f.variables['salt'][:]  
    print("temp: range = ", temp.min(),temp.max())
    print("salt: range = ", salt.min(),salt.max())
    f.close()

    return temp, salt

def read_m5llgrid(fngrd="grid_spec.latlon.nc"):
    """
    read in grid info related to lat-lon m5 odas analysis
    """
    f = Dataset(fngrd,"r")
    f.set_auto_mask(False)  
    lat   = f.variables['grid_y_T'][:]
    lon   = f.variables['grid_x_T'][:]
    depth = f.variables['zt'][:]
    num_levels = f.variables['num_levels'][:]
    wet        = f.variables['wet'][:]
    olmask     = f.variables['olmask'][:]
    f.close()

    vnames = ["lat","lon","depth","wet","num_levels","olmask"]
    vs = [lat,lon,depth, wet, num_levels, olmask ]
    for i, v in enumerate(vs):
        print("var, range, shape=", vnames[i], v.min(), v.max(), "|", v.shape)

    return lat, lon, depth, olmask, num_levels, wet

"""
netcdf WOA05_pottemp_salt {
dimensions:
	TIME = UNLIMITED ; // (12 currently)
	DEPTH = 33 ;
	LAT = 180 ;
	LON = 360 ;
variables:
	double TIME(TIME) ;
		TIME:units = "days since 0001-01-01 00:00:00" ;
		TIME:calendar = "noleap" ;
		TIME:cartesian_axis = "T" ;
	double DEPTH(DEPTH) ;
		DEPTH:units = "m" ;
		DEPTH:direction = -1 ;
		DEPTH:cartesian_axis = "Z" ;
	double LAT(LAT) ;  -89.5 ->  89.5
		LAT:units = "degrees_north" ;
		LAT:cartesian_axis = "Y" ;
	double LON(LON) ; 0.5 -> 359.5
		LON:units = "degrees_east" ;
		LON:cartesian_axis = "X" ;
	float PTEMP(TIME, DEPTH, LAT, LON) ;
		PTEMP:_FillValue = -1.e+34f ;
	float SALT(TIME, DEPTH, LAT, LON) ;
		SALT:_FillValue = -1.e+34f ;
"""
def create_m6_ts(fnout, t, s, lat, lon, depth, olmask, num_levels, wet, anadt=dt.datetime(1,1,1,0,0,0)):
    """
    write out T/S grid file which can be used by the MOM6
    """
    if os.path.exists(fnout):
        raise RuntimeError(f"file {fnout} already exits.")
        sys.exit(1)
    #
    # Flag out grids below ocean basins.
    #
    GRID_USED   = 1
    GRID_UNUSED = 0
    levmask4d = GRID_UNUSED * np.ones(s.shape,dtype=int) 
    wetmask4d = np.ones(s.shape,dtype=int)

    print("levmask4d: ",levmask4d.min(),levmask4d.max(),levmask4d.shape)
    print("num_levels: ", num_levels.min(), num_levels.max(), num_levels.shape)
    for j in range(lat.size):
        #print(j, j*100.0/lat.size,"%")
        for i in range(lon.size):
            nlev = round(num_levels[j,i])

            # note nlev is for Fortran index: we let levmask4d(1:nlev,j,i) = 1
            # corresponding python: levmask4d[0:nlev,j,i] = 1, (nlev excluded in 0:nlev)

            levmask4d[0,:nlev,j,i] = GRID_USED
            wetmask4d[:,:,j,i] = wet[j,i]

    olmask4d = np.zeros(s.shape)
    olmask4d[0,:,:,:] = olmask

    # [FIXME] CDA: need to tune T/S range for MOM6
    missing_grids = (t < -10.0) | (s < 3.0) | (levmask4d==GRID_UNUSED) | (olmask4d < 0.5)
    #missing_grids = (levmask4d==GRID_UNUSED) | (olmask4d < 0.5) | (wetmask4d < 0.5)

    # olmask: 0-land, 1-sea
    #print("mismatch:", np.sum((levmask4d==GRID_UNUSED) != (olmask4d < 0.5)) )
    t[missing_grids] = MISSING_VALUE
    s[missing_grids] = MISSING_VALUE


    #
    # Create T/S grid file
    #
    fnout_tmp = fnout+"_tmp"
    f = Dataset(fnout_tmp,"w")
    # create dim
    d_time  = f.createDimension("time",  None)
    d_lev   = f.createDimension("depth", depth.size)
    d_lat   = f.createDimension("lat",   lat.size)
    d_lon   = f.createDimension("lon",   lon.size)

    # create vars
    v_time  = f.createVariable("time", "f8", ("time") )
    v_depth = f.createVariable("depth","f8", ("depth") )
    v_lat   = f.createVariable("lat",  "f8", ("lat") )
    v_lon   = f.createVariable("lon",  "f8", ("lon") )

    v_olmask     = f.createVariable("olmask",    "f8", ("depth", "lat", "lon") )
    v_num_levels = f.createVariable("num_levels","f8", ("lat", "lon") )
    v_wet        = f.createVariable("wet",       "f8", ("lat", "lon") )
  
    # CDA: MOM6 uses _FillValue as missing_value to flag out grids in the Z file
    v_s   = f.createVariable("salt", "f4", ("time","depth", "lat", "lon"),fill_value = MISSING_VALUE )
    v_t   = f.createVariable("temp", "f4", ("time","depth", "lat", "lon"),fill_value = MISSING_VALUE )

    # add attributes
    v_time.units = anadt.strftime("days since %Y-%m-%d %H:00:00")
    v_time.calendar = "julian"
    v_time.cartesian_axis = "T"

    v_depth.units = "m"
    v_depth.direction = np.int32(-1)
    v_depth.cartesian_axis = "Z"

    v_lat.units = "degrees_north"
    v_lat.cartesian_axis = "Y"

    v_lon.units = "degrees_east"
    v_lon.cartesian_axis = "X"

    v_olmask.long_name = "3d land/sea mask, 0=land"

    # fill in values
    v_time[:]   = 0.0
    v_depth[:]  = depth
    v_lat[:]    = lat
    v_lon[:]    = lon

    v_olmask[:]     = olmask
    v_num_levels[:] = num_levels
    v_wet[:]        = wet

    v_s[:]      = s.astype(np.float32)
    v_t[:]      = t.astype(np.float32)

    f.close()

    # check bounds
    f = Dataset(fnout_tmp,"r")
    t = f.variables['temp'][:]
    s = f.variables['salt'][:]
    f.close()

    if not pass_bounds_check(t.min(),t.max(),s.min(),s.max()):
        raise RuntimeError(f"failed to pass bounds check: Check intermedidate file {fnout_tmp}")
        sys.exit(2)
    else:
        os.rename(fnout_tmp, fnout)
        print(f"Finished genearting T/S file {fnout}")

def pass_bounds_check(t_min, t_max, s_min, s_max):
    passed = True
    print("Final Check: temp range=", t_min, t_max)
    print("Final Check: salt range=", s_min, s_max)
    if (t_min < T_MIN_LIMIT or t_max > T_MAX_LIMIT ):
        passed = False
        print("temp out of permitted range [{T_MIN_LIMIT}, {T_MAX_LIMIT}]")
    elif (s_min < S_MIN_LIMIT or s_max > S_MAX_LIMIT ):
        passed = False
        print("salt out of permitted range [{S_MIN_LIMIT}, {S_MAX_LIMIT}]")
    
    return passed

def convert_m5odas_ts(fnin  = "ocean_temp_salt.res.nc",
                      fngrd = "grid_spec_latlon.nc",
                      fnout = "m5_ocean_temp_salt.res.nc",
                      anadt = dt.datetime(1,1,1,0,0,0)):

    t, s = read_m5ana(fnin)
    lat,lon,depth, olmask, num_levels, wet = read_m5llgrid(fngrd)
    create_m6_ts(fnout, t, s, lat, lon, depth, olmask, num_levels, wet, anadt)

    return None

def parse_args():
    parser = argparse.ArgumentParser(description=("convert MOM5 ODAS T/S analysis to a file to initialize T/S for MOM6"))
    parser.add_argument("fnOut",default="MOM6_TS_Z_profile.nc", type=str, help=("output T/S file to initialize MOM6"))
    parser.add_argument("--m5Ana", default="ANA-ocean_temp_salt.res.nc" ,type=str, required=True, help=("MOM5 ODAS T/S analysis on lat-lon grid"))
    parser.add_argument("--m5Grid", default="grid_spec.nc.mom5_LL1440x721", type=str, required=True, help=("grid_spec.nc in latlon grid for MOM5 ODAS"))
    parser.add_argument("--m5Date", default="0001100100", metavar="YYYYMMDDHH", required=False, help=("Analysis date time"))
    parser.add_argument("--checkOnly", action="store_true", required=False, help=("run check only but no exeuction"))

    args = parser.parse_args()

    args.fnOut  = os.path.abspath(args.fnOut)
    args.m5Ana  = os.path.abspath(args.m5Ana)
    args.m5Grid = os.path.abspath(args.m5Grid)
    args.m5Date = dt.datetime.strptime(args.m5Date,"%Y%m%d%H")
    print(args)

    if not args.checkOnly:
        convert_m5odas_ts(args.m5Ana, args.m5Grid, args.fnOut, args.m5Date)


if __name__ == '__main__':
    print(" ".join(sys.argv[:]))
    parse_args()
