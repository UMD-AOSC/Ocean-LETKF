#!/usr/bin/env python3
import os, sys, argparse, yaml
import numpy as np
import datetime as dt
import xesmf as xe
from netCDF4 import Dataset, num2date
from regrid_tools import iterative_fill_POP_core
from mom6_tools import RestoreTemplate, parse_time_units

class ConfigFill:
    def __init__(self):
        # static config applied for each processing
        self.regridder_file_path = "xesmf_wts.nc"
        self.restore_template_path = "T_restore_file.yaml"
        self.ocean_static_path = "ocean_static.nc"
        self.ocean_hgrid_path = "ocean_hgrid.nc"
        self.mom_lat_var = "geolat"  
        self.mom_lon_var = "geolon"
        self.mom_wet_var = "wet"
        self.sst_lat_var = "lat"
        self.sst_lon_var = "lon"
        self.sst_sst_var = "analysed_sst"
        self.sst_var_renamed = "SST"
        self.new_time_units = "days since 1900-01-01 00:00:00"

    def __str__(self):
        info=f"""
        regridder_file_path = {self.regridder_file_path}
        restore_template_path = {self.restore_template_path}
        ocean_static_path = {self.ocean_static_path}
        ocean_hgrid_path = {self.ocean_hgrid_path}
        mom_lat_var = {self.mom_lat_var}
        mom_lon_var = {self.mom_lon_var}
        mom_wet_var = {self.mom_wet_var}
        sst_lat_var = {self.sst_lat_var}
        sst_lon_var = {self.sst_lon_var}
        sst_sst_var = {self.sst_sst_var}
        sst_var_renamed = {self.sst_var_renamed}
        new_time_units = {self.new_time_units}
        """
        return info

    def read_cfg(self, fnin, namelist="fill"):
        with open(fnin,"r") as f:
            self.cfg = yaml.safe_load(f)
        
        for att in self.cfg[namelist].keys(): 
            if hasattr(self, att):
               setattr(self, att, self.cfg[namelist][att])
            else:
               raise RuntimeError("ConfigFill does not have attribute ({})".format(att)) 
               exit(11)
        print(self)

def parseCommandLine():
    
    parser = argparse.ArgumentParser(description=(""))
    parser.add_argument("--config_path", required=True, help=(
                "yaml file for storing all configure files"))
    parser.add_argument("--sst_file_path",required=True,type=str,help=(
                "original L4 SST data"))
    parser.add_argument("--remapped_file_path", required=True, type=str, help=(
                "output SST remapped to MOM6 tripolar grids"))
    parser.add_argument("--debug_file_path", required=False, default=None, help=(
                "file with addtional debug info"))
    args =parser.parse_args()
    print(args)

    return args

def write_nc_file(fnout, lat1d_l4sst, lon1d_l4sst, lat2d_mom, lon2d_mom, \
                    sst_unfilled, sst_filled, sst_remapped_mom, wet_mom):
    #print("write out netcdf file: {}".format(fnout))
    if os.path.exists(fnout):
        raise Exception("fnout already exists at: {}".format(fnout))
        sys.exit(2)

    fnout = os.path.abspath(fnout)
    f = Dataset(fnout, mode='w', format='NETCDF4_CLASSIC')
    f.createDimension('lat_l4sst', lat1d_l4sst.size)
    f.createDimension('lon_l4sst', lon1d_l4sst.size)
    nlat_mom, nlon_mom = lat2d_mom.shape
    f.createDimension('lat_mom', nlat_mom)
    f.createDimension('lon_mom', nlon_mom)

    lat1d_l4sst_to_file = f.createVariable('lat1d_l4sst',np.float32, ('lat_l4sst',))
    lon1d_l4sst_to_file = f.createVariable('lon1d_l4sst',np.float32, ('lon_l4sst',))
    lat2d_mom_to_file = f.createVariable('lat2d_mom',np.float32, ('lat_mom','lon_mom'))
    lon2d_mom_to_file = f.createVariable('lon2d_mom',np.float32, ('lat_mom','lon_mom'))

    unfill_to_file    = f.createVariable('sst_unfilled',    np.float32, ('lat_l4sst','lon_l4sst'))
    filled_to_file = f.createVariable('sst_filled',np.float32, ('lat_l4sst','lon_l4sst'))
    remapped_to_file  = f.createVariable('sst_remapped_mom',np.float32, ('lat_mom','lon_mom'))
    wet_to_file       = f.createVariable('wet',             np.float32, ('lat_mom','lon_mom'))

    lat1d_l4sst_to_file[:] = lat1d_l4sst
    lon1d_l4sst_to_file[:] = lon1d_l4sst
    lat2d_mom_to_file[:]   = lat2d_mom
    lon2d_mom_to_file[:]   = lon2d_mom

    unfill_to_file[:]   = sst_unfilled
    filled_to_file[:]   = sst_filled
    remapped_to_file[:] = sst_remapped_mom
    wet_to_file[:]      = wet_mom

    f.close()


def main(args):

    #
    # read configurations 
    #
    cfg = ConfigFill()
    cfg.read_cfg(args.config_path,namelist="fill")

    #
    # create the regridder info
    #
    print("="*80+"\ncreate the regridder info\n"+"-"*80)

    # read in output tripolar MOM lon, lat, land/ocean mask (ocean: wet=1)
    f = Dataset(cfg.ocean_static_path)
    lat2d_grd_out = f.variables[cfg.mom_lat_var][:].astype(np.float32)
    lon2d_grd_out = f.variables[cfg.mom_lon_var][:].astype(np.float32)
    wet2d_grd_out = f.variables[cfg.mom_wet_var][:].astype(np.float32)
    f.close()

    print("lat2d_grd_out: shape, min, max=",lat2d_grd_out.shape, np.min(lat2d_grd_out), np.max(lat2d_grd_out))
    print("lon2d_grd_out: shape, min, max=",lon2d_grd_out.shape, np.min(lon2d_grd_out), np.max(lon2d_grd_out))
    print("wet2d_grd_out: shape, min, max=", wet2d_grd_out.shape, np.min(wet2d_grd_out), np.max(wet2d_grd_out))

    # read in the input lon, lat
    f = Dataset(args.sst_file_path)
    lat1d_grd_in = f.variables[cfg.sst_lat_var][:]
    lon1d_grd_in = f.variables[cfg.sst_lon_var][:]
    sst_unfilled = np.squeeze( f[cfg.sst_sst_var][:] )
    time = num2date(f.variables['time'],f.variables['time'].units, \
                     only_use_cftime_datetimes=False, \
                     only_use_python_datetimes=True)[0]
    f.close()
    print("sst time={}".format(time))

    dtime_units, base_time = parse_time_units(cfg.new_time_units)
    new_dtime = (time-base_time)/dtime_units

    print("lat1d_grd_in: shape, min, max=",lat1d_grd_in.shape, np.min(lat1d_grd_in), np.max(lat1d_grd_in))
    print("lon1d_grd_in: shape, min, max=",lon1d_grd_in.shape, np.min(lon1d_grd_in), np.max(lon1d_grd_in))


    # read in the regridder
    grd_in = {"lon": lon1d_grd_in, "lat": lat1d_grd_in}
    grd_out = {"lon": lon2d_grd_out, "lat": lat2d_grd_out}
    interp_method = "bilinear"
    periodic = True
    fn_regridder_in = os.path.abspath(cfg.regridder_file_path)
    regridder = xe.Regridder(grd_in, grd_out, interp_method, weights=fn_regridder_in, periodic=periodic)
    regridder.filename = fn_regridder_in  # a bug of xesmf
    print(regridder)


    #
    # flood the globle
    #   
    print("="*80+"\nflood the land\n"+"-"*80)

    wk2d = sst_unfilled.data
    missing_value = sst_unfilled.fill_value
    fillmask = sst_unfilled.mask
    iterative_fill_POP_core(var=wk2d, fillmask=fillmask, missing_value=missing_value, tol=1.e-4, ltripole=False, nitermax=1000)
    sst_filled  = wk2d.copy()

    #
    # interp from uniform latlon to mom tripolar
    #
    print("="*80+"\ninterpolate to MOM tripolar grids\n"+"-"*80)

    print("BEFORE remapped (global): min, max=", np.min(wk2d), np.max(wk2d))
    print("remap SST from latlon to mom tripolar")
    sst_remapped = regridder(wk2d)
    print("AFTER remapped (global): min, max=", np.min(sst_remapped), np.max(sst_remapped))
    print("AFTER remapped (ocean): min, max=", np.min(sst_remapped[wet2d_grd_out>0.5]), np.max(sst_remapped[wet2d_grd_out>0.5]))
    n_pts_missing = np.sum( np.isnan(sst_remapped) & (wet2d_grd_out>0.5))
    print("n_pts_missing=",n_pts_missing)
    if n_pts_missing != 0:
        raise Exception("Have {} NaN in remapped tripolar grids over ocean when interpolating to tripolar grids".format(n_pts_missing))
        sys.exit(3)

    sst_remapped -= 273.15
    print("AFTER K->degC (global): min, max=", np.min(sst_remapped), np.max(sst_remapped))
    print("AFTER K->degC (ocean): min, max=", np.min(sst_remapped[wet2d_grd_out>0.5]), np.max(sst_remapped[wet2d_grd_out>0.5]))

    #
    # write out the daily restore file
    #
    print("="*80+"\nwrite out\n"+"-"*80)

    tpl = RestoreTemplate()
    tpl.read(os.path.abspath(cfg.restore_template_path))
    input_dims, input_vars = tpl.get_restore_static(\
                               ocean_static_file = os.path.abspath(cfg.ocean_static_path), \
                               ocean_hgrid_file = os.path.abspath(cfg.ocean_hgrid_path))

    input_vars[cfg.sst_var_renamed] = np.expand_dims(sst_remapped,axis=0)
    input_vars['time'] = np.array([new_dtime])
    tpl.create_restore(args.remapped_file_path, input_dims, input_vars)

    #
    # write out the daily original/filled/remapped sst
    #
    if args.debug_file_path is not None:
        write_nc_file(args.debug_file_path, lat1d_grd_in, lon1d_grd_in, lat2d_grd_out, lon2d_grd_out, \
                  sst_unfilled, sst_filled, sst_remapped, wet2d_grd_out)


if __name__ == '__main__':
    args = parseCommandLine()
    main(args)
