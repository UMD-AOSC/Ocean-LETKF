#!/usr/bin/env python3

import argparse, os
import numpy as np
import xesmf as xe
from netCDF4 import Dataset
from config_tools import Config

class ConfigRegrid(Config):
    def __init__(self):
        # static config applied for each processing
        self.woa_grid_path = "OSTIA.nc"
        self.woa_lat_var = "lat"
        self.woa_lon_var = "lon"
        self.l4_grid_path = "L4SSS.nc"
        self.l4_lat_var = "lat"
        self.l4_lon_var = "lon"
        self.regridder_file_path = "xesmf_wts.nc"


def parseCommandLine():
    parser = argparse.ArgumentParser(description=("Generate regrid weight file by xesmf"))
    parser.add_argument("config_path", help=("config file path"))
    args = parser.parse_args()
    print(args)
    return args



def main(args):

    cfg = ConfigRegrid()
    cfg.read_cfg(args.config_path, "regrid2")

    # WOA18 latlon grids as output grids
    fn_grd_out = os.path.abspath(cfg.woa_grid_path)
    f = Dataset(fn_grd_out)
    lat1d_grd_out = f.variables[cfg.woa_lat_var][:]
    lon1d_grd_out = f.variables[cfg.woa_lon_var][:]
    f.close()

    print("lat1d_grd_out: shape, min, max=",lat1d_grd_out.shape, np.min(lat1d_grd_out), np.max(lat1d_grd_out))
    print("lon1d_grd_out: shape, min, max=",lon1d_grd_out.shape, np.min(lon1d_grd_out), np.max(lon1d_grd_out))

    # L4 latlon grids as input grids
    fn_grd_in = os.path.abspath(cfg.l4_grid_path)
    f = Dataset(fn_grd_in)
    lat1d_grd_in = f.variables[cfg.l4_lat_var][:]
    lon1d_grd_in = f.variables[cfg.l4_lon_var][:]
    f.close()

    print("lat1d_grd_in: shape, min, max=",lat1d_grd_in.shape, np.min(lat1d_grd_in), np.max(lat1d_grd_in))
    print("lon1d_grd_in: shape, min, max=",lon1d_grd_in.shape, np.min(lon1d_grd_in), np.max(lon1d_grd_in))

    # we only generate wts to WOA18 grids within the L4 grid range; otherwise, xesmf fails to work
    lat_max_allowed = lat1d_grd_out[lat1d_grd_out<lat1d_grd_in.max()][-1]
    idx_lat_max_allowed = np.where(lat1d_grd_out==lat_max_allowed)[0][0]

    lat_min_allowed = lat1d_grd_out[lat1d_grd_out>lat1d_grd_in.min()][0]
    idx_lat_min_allowed = np.where(lat1d_grd_out==lat_min_allowed)[0][0]

    print("min_WOA18_lat can be mapped from L4 ESA: idx_woa18, val_woa18=", idx_lat_min_allowed, lat1d_grd_out[idx_lat_min_allowed])
    print("max_WOA18_lat can be mapped from L4 ESA: idx_woa18, val_woa18=", idx_lat_max_allowed, lat1d_grd_out[idx_lat_max_allowed])

    # generate the regridder 
    grd_in = {"lon": lon1d_grd_in, "lat": lat1d_grd_in}  
    grd_out = {"lon": lon1d_grd_out, "lat": lat1d_grd_out[idx_lat_min_allowed:idx_lat_max_allowed+1]}
    interp_method = "bilinear"
    periodic = True
    regridder = xe.Regridder(grd_in, grd_out, interp_method, periodic=periodic)
  
    fn_regridder_out = os.path.abspath(cfg.regridder_file_path)
    print("write regrid info into the file: ", fn_regridder_out)
    regridder.filename = fn_regridder_out
    regridder.to_netcdf()

    print(regridder)

    return None

if __name__ == '__main__':
    args = parseCommandLine()
    main(args)

