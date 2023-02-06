#!/usr/bin/env python3

import argparse, os, sys
import yaml
import datetime as dt
import numpy as np
import xesmf as xe
from netCDF4 import Dataset
from mom6_tools import RestoreTemplate, parse_time_units


class ConfigRemap:
    def __init__(self):
        # static config applied for each processing
        self.regridder_file_path = "xesmf_wts.nc"
        self.restore_template_path = "S_restore_file.yaml"
        self.ocean_static_path = "ocean_static.nc"
        self.ocean_hgrid_path = "ocean_hgrid.nc"
        self.mom_lat_var = "geolat"
        self.mom_lon_var = "geolon"
        self.mom_wet_var = "wet"
        self.woa_grid_path = "woa_grid.nc"
        self.woa_lat_var = "lat"
        self.woa_lon_var = "lon"
        self.woa_var = 'SALT'
        self.woa_var_renamed = 'SALT'
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
        woa_grid_path = {self.woa_grid_path}
        woa_lat_var = {self.woa_lat_var}
        woa_lon_var = {self.woa_lon_var}
        woa_var     = {self.woa_var}
        woa_var_renamed = {self.woa_var_renamed}
        new_time_units = {self.new_time_units}
        """
        return info

    def read_cfg(self, fnin, namelist="remap"):
        with open(fnin,"r") as f:
            self.cfg = yaml.safe_load(f)

        for att in self.cfg[namelist].keys():
            if hasattr(self, att):
               setattr(self, att, self.cfg[namelist][att])
            else:
               raise RuntimeError("ConfigRemap does not have attribute ({})".format(att))
               exit(11)
        print(self)

def parseCommandLine():
    parser = argparse.ArgumentParser(description=("Generate regrid weight file by xesmf"))
    parser.add_argument("--remapped_file_path",required=True,help=("output file path"))
    parser.add_argument("--config_path", required=True, help=())
    parser.add_argument("--sss_filled_path", required=True, help=(
                "path of the flooded source grid file (e.g., WOA18)"))
    parser.add_argument("--year",required=True, type=int, help=())

    args = parser.parse_args()
    print(args)
    return args

def main(args):


    cfg = ConfigRemap()
    cfg.read_cfg(args.config_path,namelist="remap")
    #
    # create regridder info
    #

    # read in output tripolar lon, lat, land/ocean mask (ocean: wet==1)
    f = Dataset(cfg.ocean_static_path)
    lat2d_grd_out = f.variables[cfg.mom_lat_var][:].astype(np.float32)
    lon2d_grd_out = f.variables[cfg.mom_lon_var][:].astype(np.float32)
    wet2d_grd_out = f.variables[cfg.mom_wet_var][:].astype(np.float32)
    f.close()

    print("lat2d_grd_out: shape, min, max=",lat2d_grd_out.shape, np.min(lat2d_grd_out), np.max(lat2d_grd_out))
    print("lon2d_grd_out: shape, min, max=",lon2d_grd_out.shape, np.min(lon2d_grd_out), np.max(lon2d_grd_out))
    print("wet2d_grd_out: shape, min, max=", wet2d_grd_out.shape, np.min(wet2d_grd_out), np.max(wet2d_grd_out))

    # read in the input lon, lat
    f = Dataset(cfg.woa_grid_path)
    lat1d_grd_in = f.variables[cfg.woa_lat_var][:]
    lon1d_grd_in = f.variables[cfg.woa_lon_var][:]
    f.close()

    print("lat1d_grd_in: shape, min, max=",lat1d_grd_in.shape, np.min(lat1d_grd_in), np.max(lat1d_grd_in))
    print("lon1d_grd_in: shape, min, max=",lon1d_grd_in.shape, np.min(lon1d_grd_in), np.max(lon1d_grd_in))

    # generate the regridder 
    grd_in = {"lon": lon1d_grd_in, "lat": lat1d_grd_in}  
    grd_out = {"lon": lon2d_grd_out, "lat": lat2d_grd_out}
    interp_method = "bilinear"
    periodic = True
    fn_regridder_in = os.path.abspath(cfg.regridder_file_path)
    regridder = xe.Regridder(grd_in, grd_out, interp_method, weights=fn_regridder_in, periodic=periodic)
    regridder.filename = fn_regridder_in  # a bug of xesmf 
    print(regridder)

    #
    # load SSS_filled file
    #
    f = Dataset(args.sss_filled_path)
    sss_filled = f.variables[cfg.woa_var][:]
    print("shape(sss_filled)=",sss_filled.shape, np.min(sss_filled), np.max(sss_filled))
    f.close()

    sss_filled[sss_filled.mask] = np.nan
    sss_filled.mask = False
    sss_filled.fill_value = 0

    #
    # interpolate SSS from latlon to tripolar grids
    #
    remapped_sss_filled = np.zeros((sss_filled.shape[0],lat2d_grd_out.shape[0],lat2d_grd_out.shape[1]))
    for i in range(sss_filled.shape[0]):
        wk2d = regridder(sss_filled[i,:,:])
        print("remapped: i, min, max=", i, np.min(wk2d), np.max(wk2d))
        if np.sum( np.isnan(wk2d) & (wet2d_grd_out>0.5)) != 0:
            raise Exception("Have NaN in remapped tripolar grids over ocean when filling remapped_sss_filled[{},:,:]".format(i))
            sys.exit(3)
        remapped_sss_filled[i,:,:] = wk2d.copy()

    #
    # write out the restore file
    #
    dtime_units, base_time = parse_time_units(cfg.new_time_units)
    dtime = np.zeros((12))
    for i in range(12):
        month_time = dt.datetime(args.year,i+1,1,0,0,0)
        dtime[i] = (month_time-base_time)/dtime_units

    tpl = RestoreTemplate()
    tpl.read(os.path.abspath(cfg.restore_template_path))
    input_dims, input_vars = tpl.get_restore_static(\
                               ocean_static_file = os.path.abspath(cfg.ocean_static_path), \
                               ocean_hgrid_file = os.path.abspath(cfg.ocean_hgrid_path))

    input_vars[cfg.woa_var_renamed] = remapped_sss_filled
    input_vars['time'] = dtime
    tpl.create_restore(args.remapped_file_path, input_dims, input_vars)

if __name__ == '__main__':
    args = parseCommandLine()
    main(args)

