#!/usr/bin/env python3

import argparse
import numpy as np
from gradsio import gradsctl
from netCDF4 import Dataset
import os, sys

DEFAULT_MISSING_VALUE = 1.e20

def write_ncdiag(fn, ctl, v3d, v2d):
    fn = os.path.abspath(fn)

    if os.path.exists(fn):
        raise RuntimeError("{} already exists there. Exit...".format(fn))
        sys.exit(1)

    ds = Dataset(fn, mode='w' , format="NETCDF4_CLASSIC")

    # Note: actually U/V does not use xh, yh due to Arakawa-C stagged grid, but it should not
    # matter as long as our validation does not use u,v related obs
    dim_x_name = "xh"
    dim_y_name = "yh"
    dim_z_name = "zl"
    dim_t_name = "Time"
    ds.createDimension(dim_x_name, ctl.nx)
    ds.createDimension(dim_y_name, ctl.ny)
    ds.createDimension(dim_z_name, ctl.nz)
    ds.createDimension(dim_t_name, None)

    var_u_name = "u"
    var_v_name = "v"
    var_s_name = "salt"
    var_t_name = "temp"
    var_h_name = "h"
    v3d_shape = (dim_t_name, dim_z_name,dim_y_name,dim_x_name)
    var_u = ds.createVariable(var_u_name, "f4", v3d_shape)
    var_v = ds.createVariable(var_v_name, "f4", v3d_shape)
    var_s = ds.createVariable(var_s_name, "f4", v3d_shape)
    var_t = ds.createVariable(var_t_name, "f4", v3d_shape)
    var_h = ds.createVariable(var_h_name, "f4", v3d_shape)
        
    var_ssh_name = "ssh"
    var_sst_name = "sst"
    var_sss_name = "sss"
    var_eta_name = "eta"
    var_mld_name = "mld"
    var_ssu_name = "ssu" # additional variables
    var_ssv_name = "ssv"
    v2d_shape = (dim_t_name, dim_y_name,dim_x_name)
    var_ssh = ds.createVariable(var_ssh_name, "f4", v2d_shape)
    var_sst = ds.createVariable(var_sst_name, "f4", v2d_shape)
    var_sss = ds.createVariable(var_sss_name, "f4", v2d_shape)
    var_eta = ds.createVariable(var_eta_name, "f4", v2d_shape)
    var_mld = ds.createVariable(var_mld_name, "f4", v2d_shape)

    var_ssu = ds.createVariable(var_ssu_name, "f4", v2d_shape)
    var_ssv = ds.createVariable(var_ssv_name, "f4", v2d_shape)
   
    print("writing 3d vars")
    var_u[0,:,:,:] = v3d['u'].astype("f4")
    var_v[0,:,:,:] = v3d['v'].astype("f4")
    var_s[0,:,:,:] = v3d['s'].astype("f4")
    var_t[0,:,:,:] = v3d['t'].astype("f4")
    var_h[0,:,:,:] = v3d['h'].astype("f4")

    print("writing 2d vars")
    var_ssh[0,:,:] = v2d['ssh'].astype("f4")
    var_sst[0,:,:] = v2d['sst'].astype("f4")
    var_sss[0,:,:] = v2d['sss'].astype("f4")
    var_eta[0,:,:] = v2d['eta'].astype("f4")
    var_mld[0,:,:] = v2d['mld'].astype("f4")
    var_ssu[0,:,:] = v2d['ssu'].astype("f4")
    var_ssv[0,:,:] = v2d['ssv'].astype("f4")

    ds.close()

    # add missing values, since they will be read by read_diag for obsop
    ds = Dataset(fn, mode='r+' , format="NETCDF4_CLASSIC")
    varlist = ["u","v","salt","temp","h","ssh","sst","sss","eta","mld","ssu","ssv"]
    for var in varlist:
        ds.variables[var].missing_value = np.array(DEFAULT_MISSING_VALUE, dtype="f4")
    ds.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=("convert binary ens mean/sprd files into NetCDF files"))
    parser.add_argument("fnin",help=("input binary file"))
    parser.add_argument("fnout",help=("output NetCDF file"))
    args = parser.parse_args()
    print(args)

    fnin = os.path.abspath(args.fnin)
    fnout = os.path.abspath(args.fnout)

    # CTL description
    # https://github.com/UMD-AOSC/Ocean-LETKF/blob/83d53a28b102cce47c99e0d78a9ca43340f1f51e/src/model_specific/mom6/params_model.f90
    mom = gradsctl(nv3d = 5, nv2d = 5, nx = 1440, ny = 1080, nz = 75, \
                   v3dnames = ['u','v','t','s','h'], \
                   v2dnames = ['ssh', 'sst', 'sss', 'eta', 'mld'])
    print(mom)
    print(mom.nx, mom.ny, mom.nz)

    v3d, v2d = mom.readVars3d2d(fnin, vtype='<f4')
    for var in v3d:
        print( "Var, min, max=",v3d[var].shape, np.min(v3d[var]), np.max(v3d[var]) )
    for var in v2d:
        print( "Var, min, max=",v2d[var].shape, np.min(v2d[var]), np.max(v2d[var]) )

    v2d['sst'][:] = v3d['t'][0,:,:].copy()
    v2d['sss'][:] = v3d['s'][0,:,:].copy()

    # append additional 2d vars
    v2d['ssu'] = v3d['u'][0,:,:].copy()
    v2d['ssv'] = v3d['v'][0,:,:].copy()

    write_ncdiag(fnout, mom, v3d, v2d)
    print("FINISHED!")

