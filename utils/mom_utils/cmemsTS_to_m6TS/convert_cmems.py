#!/usr/bin/env python3

import sys, os, numpy as np, datetime as dt
import argparse
from netCDF4 import Dataset
import multiprocessing
from multiprocessing import shared_memory
from numba import jit
from regrid_tools import iterative_fill_POP_core, iterative_fill_sor



missing_pts_ic_value = 1.e10 

def load_cmems_oana( fnin  = "cmems_mod_glo_phy-all_my_0.25deg_P1d-m-20230306.nc", 
                     oname = "oras"):
    """
    75-level, 0.25deg, daily ocean analysis ensemble 
    ref: 
    - https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_ENS_001_031/description
    - http://marine.copernicus.eu/documents/QUID/CMEMS-GLO-QUID-001_031.pdf

    input file name example:
    cmems_mod_glo_phy-all_my_0.25deg_P1d-m-20230306.nc
    """
    
    if oname not in ["oras","glor","cglo"]:
        raise RuntimeError(f"oname (f{oname}) not found in the dataset (oras, glor, cglo)")
        sys.exit(1)

    f = Dataset(fnin)
    f.set_auto_mask(False)
    z1d   = f.variables["depth"][:]
    lat1d = f.variables["latitude"][:]
    lon1d = f.variables["longitude"][:]

    tname = f"thetao_{oname}"
    sname = f"so_{oname}"
    t = f.variables[tname][:].squeeze()
    tmin, tmax = f.variables[tname].valid_min, f.variables[tname].valid_max
    t_pts_invalid = (t<tmin) | (t>tmax)
    t[t_pts_invalid] = np.nan

    s = f.variables[sname][:].squeeze()
    smin, smax = f.variables[sname].valid_min, f.variables[sname].valid_max
    s_pts_invalid = (s<tmin) | (s>tmax)
    s[s_pts_invalid] = np.nan

	#int64 time(time) ;
	#	time:axis = "T" ;
	#	time:long_name = "Time" ;
	#	time:standard_name = "time" ;
	#	time:units = "seconds since 1950-01-01" ;
	#	time:calendar = "gregorian" ;
    delT = dt.timedelta(seconds=int(f.variables['time'][0]))
    T0   = dt.datetime(1950,1,1,0,0,0)
    odate = delT + T0

    f.close()
    print(f"t: shape, min, max = {t.shape}, {np.nanmin(t)}, {np.nanmax(t)}")
    print(f"s: shape, min, max = {s.shape}, {np.nanmin(s)}, {np.nanmax(s)}")
    print(f"lon1d: shape = {lon1d.shape}, val = {lon1d}")
    print(f"lat1d: shape = {lat1d.shape}, val = {lat1d}")
    print(f"z1d:   shape = {z1d.shape},   val = {z1d}")
    print(odate)

    return lon1d, lat1d, z1d, t, s, odate

def flood_level(ilev, shm_name, shape, use_fill_pop=False,tol=1.e-2,nitermax=100):
    shm = shared_memory.SharedMemory(name=shm_name)
    t_f = np.ndarray(shape, dtype=np.float64, buffer=shm.buf)
    t_f_lev = t_f[ilev,:,:]
    ## flood land and -90 -> -80
    pts_to_fill = np.isnan(t_f_lev)

    use_fill_pop = False
    if use_fill_pop:
        t_f_lev[pts_to_fill] = missing_pts_ic_value
        iterative_fill_POP_core(var           = t_f_lev,
                                fillmask      = pts_to_fill,
                                missing_value = missing_pts_ic_value,
                                tol           = tol, ltripole = False, nitermax = nitermax)
    else: 
        iterative_fill_sor(nlat = t_f_lev.shape[0],
                           nlon = t_f_lev.shape[1], 
                           var  = t_f_lev, 
                           fillmask = pts_to_fill, 
                           tol = tol, rc = 1.5, ltripole = False, max_iter = nitermax)

    t_f_lev[0,:]  = t_f_lev[1,:].copy() # set South-most row
    print(f"finish flooding at lev={ilev}: min, max=",t_f_lev.min(),t_f_lev.max())
    shm.close()
    

def flood_globe(lat1d, t, s, npes=4):
    
    # redefine the grid
    lat1d_f = np.arange(-90,90+0.25,0.25)
    print(t.shape)
    ix = np.argwhere( (lat1d_f-lat1d[0]) == 0.0)[0]
    #[nz, nlat, nlon]
    t_f = np.ones((t.shape[0], lat1d_f.size, t.shape[2])) * np.nan
    s_f = np.ones((t.shape[0], lat1d_f.size, t.shape[2])) * np.nan
    print(f"t_f: shape = {t_f.shape}")
 
    ib = np.argwhere( (lat1d_f-lat1d[0]) == 0.0)[0][0]
    ie = lat1d_f.size + 1
    print(ib,ie)

    # fill -80 -> 90
    t_f[:,ib:ie,:] = t[:,:,:].copy()
    s_f[:,ib:ie,:] = s[:,:,:].copy()

    flood_t = True
    flood_s = True
    use_fill_pop = False

    shm = shared_memory.SharedMemory(create=True, size = t_f.size * 8) #size in bytes
    var_f_shared = np.ndarray(t_f.shape,dtype=np.float64, buffer=shm.buf)

    if flood_t:
       print("="*80+"\nFlood T")
       var_f_shared[:] = t_f
       tol = 5.e-3; nitermax = 2000
       with multiprocessing.Pool(processes=npes) as pool:
            pool.starmap(flood_level, [ (ilev, shm.name, t_f.shape, use_fill_pop, tol, nitermax) for ilev in range(t.shape[0]) ]   )
       t_f[:] = var_f_shared
       t_f[-1,:,:] = missing_pts_ic_value   # deepest level has no valid value

    if flood_s:
       print("="*80+"\nFlood S")
       var_f_shared[:] = s_f
       tol = 5.e-3; nitermax = 2000
       with multiprocessing.Pool(processes=npes) as pool:
            pool.starmap(flood_level, [ (ilev, shm.name, t_f.shape, use_fill_pop, tol, nitermax) for ilev in range(s.shape[0]) ]   )
       s_f[:] = var_f_shared
       s_f[-1,:,:] = missing_pts_ic_value   # deepest level has no valid value
    
    shm.close()
    shm.unlink()

    return lat1d_f, t_f, s_f

def write_netcdf(fnout, lon1d, lat1d, z1d, t, s, odate):

    nz, ny, nx = t.shape
    t = t.reshape((1,nz,ny,nx))
    s = s.reshape((1,nz,ny,nx))

    f = Dataset(fnout,"w")
    # create dim
    d_time  = f.createDimension("time",  None)
    d_lev   = f.createDimension("depth", z1d.size)
    d_lat   = f.createDimension("lat",   lat1d.size)
    d_lon   = f.createDimension("lon",   lon1d.size)

    # create vars
    v_time  = f.createVariable("time", "f8", ("time") )
    v_depth = f.createVariable("depth","f8", ("depth") )
    v_lat   = f.createVariable("lat",  "f8", ("lat") )
    v_lon   = f.createVariable("lon",  "f8", ("lon") )

    # CDA: MOM6 uses _FillValue as missing_value to flag out grids in the Z file
    v_s   = f.createVariable("salt", "f4", ("time","depth", "lat", "lon"),fill_value = missing_pts_ic_value )
    v_t   = f.createVariable("temp", "f4", ("time","depth", "lat", "lon"),fill_value = missing_pts_ic_value )

    # add attributes
    v_time.units = odate.strftime("days since %Y-%m-%d %H:00:00")
    v_time.calendar = "gregorian" # from cmems time attribute
    v_time.cartesian_axis = "T"

    v_depth.units = "m"
    v_depth.direction = np.int32(-1)
    v_depth.cartesian_axis = "Z"

    v_lat.units = "degrees_north"
    v_lat.cartesian_axis = "Y"

    v_lon.units = "degrees_east"
    v_lon.cartesian_axis = "X"

    # fill in values
    v_time[:]   = 0.0
    v_depth[:]  = z1d
    v_lat[:]    = lat1d
    v_lon[:]    = lon1d
    v_s[:]      = s.astype(np.float32)
    v_t[:]      = t.astype(np.float32)
    f.close()

def parse_args():
    parser = argparse.ArgumentParser(description=("convert CMEMS ocean analysis to a file to initialize T/S for MOM6"))
    parser.add_argument("fnout", default="IC_TS_mom6.nc", type=str, help=("output T/S file to initialize MOM6"))
    parser.add_argument("--fnin", default="cmems_mod_glo_phy-all_my_0.25deg_P1d-m-20230306.nc" ,type=str, required=True, help=("input cmems file name"))
    parser.add_argument("--oname", default="oras", type=str, required=True, help=("ocean analysis to use (oras, glor, cglo)"))
    parser.add_argument("--npes", default=8, type=int, required=False, help=("num of procceses for parallel-processing"))

    args = parser.parse_args()

    args.fnout = os.path.abspath(args.fnout)
    args.fnin  = os.path.abspath(args.fnin)
    print(args)

    lon1d, lat1d, z1d, t, s, odate = load_cmems_oana(args.fnin, args.oname)
    lat1d_f, t_f, s_f = flood_globe(lat1d, t, s, npes=args.npes)
    write_netcdf(args.fnout, lon1d, lat1d_f, z1d, t_f, s_f, odate)

if __name__ == '__main__':
    args = parse_args()

