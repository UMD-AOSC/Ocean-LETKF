#!/usr/bin/env python3

import datetime as dt
import calendar, argparse, os, sys
import numpy as np
import xesmf as xe
from netCDF4 import Dataset, num2date
from regrid_tools import fill_nan_grds_aux, iterative_fill_POP_core
from mom6_tools import RestoreTemplate, parse_time_units
from config_tools import Config

class ConfigMerge(Config):
    def __init__(self):
        # static config applied for each processing
        self.regridder_l4_to_woa_path  = "wts_l4sss_to_woa.nc"
        self.regridder_woa_to_mom_path = "wts_woa_to_mom.nc"
        self.restore_template_path = "T_restore_file.yaml"
        self.ocean_static_path = "ocean_static.nc"
        self.ocean_hgrid_path = "ocean_hgrid.nc"
        self.mom_lat_var = "geolat"  
        self.mom_lon_var = "geolon"
        self.mom_wet_var = "wet"
        self.woa_grid_path = "woa_grid.nc"
        self.woa_lat_var = "lat"
        self.woa_lon_var = "lon"
        self.l4_lat_var = "lat"
        self.l4_lon_var = "lon"
        self.sss_var_renamed = "SALT"
        self.new_time_units = "days since 1900-01-01 00:00:00"


def parseCommandLine():
    parser = argparse.ArgumentParser(description=("Generate regrid weight file by xesmf"))
    parser.add_argument("--config_path", required=True, help=())
    parser.add_argument("--woa18_filled_path", required=True, help=(
                    "path of the source grid file (e.g., WOA18)"))
    parser.add_argument("--l4_sss_path", required=True, help=("Level-4 SSS file"))
    parser.add_argument("--radius", required=True, type=int, help=("man-made gap between L4 & WOA"))
    parser.add_argument("--output_file_path", required=True, type=str, help=("output file name"))
    parser.add_argument("--debug_file_path", default=None, type=str, help=("inmediate file for debug"))
    args = parser.parse_args()
    print(args)

    return args


def write_nc_file(fnout, lat1d_woa, lon1d_woa, lat2d_mom, lon2d_mom, \
                    sss_unfilled, sss_filled_step1, sss_filled_step2, sss_remapped_mom, wet_mom):
    print("write out netcdf file: {}".format(fnout))
    if os.path.exists(fnout):
        raise Exception("fnout already exists at: {}".format(fnout))
        sys.exit(2)

    fnout = os.path.abspath(fnout)
    f = Dataset(fnout, mode='w', format='NETCDF4_CLASSIC')
    f.createDimension('lat_woa', lat1d_woa.size)
    f.createDimension('lon_woa', lon1d_woa.size)
    f.createDimension('lat_mom', lat2d_mom.shape[0])
    f.createDimension('lon_mom', lat2d_mom.shape[1])

    lat1d_woa_to_file = f.createVariable('lat1d_woa',np.float32, ('lat_woa',))
    lon1d_woa_to_file = f.createVariable('lon1d_woa',np.float32, ('lon_woa',))
    lat2d_mom_to_file = f.createVariable('lat2d_mom',np.float32, ('lat_mom','lon_mom'))
    lon2d_mom_to_file = f.createVariable('lon2d_mom',np.float32, ('lat_mom','lon_mom'))

    unfill_to_file    = f.createVariable('sss_unfilled',    np.float32, ('lat_woa','lon_woa'))
    filled_s1_to_file = f.createVariable('sss_filled_step1',np.float32, ('lat_woa','lon_woa'))
    filled_s2_to_file = f.createVariable('sss_filled_step2',np.float32, ('lat_woa','lon_woa'))
    remapped_to_file  = f.createVariable('sss_remapped_mom',np.float32, ('lat_mom','lon_mom'))
    wet_to_file       = f.createVariable('wet',             np.float32, ('lat_mom','lon_mom'))

    lat1d_woa_to_file[:] = lat1d_woa
    lon1d_woa_to_file[:] = lon1d_woa
    lat2d_mom_to_file[:] = lat2d_mom
    lon2d_mom_to_file[:] = lon2d_mom

    unfill_to_file[:]    = sss_unfilled
    filled_s1_to_file[:] = sss_filled_step1
    filled_s2_to_file[:] = sss_filled_step2
    remapped_to_file[:]  = sss_remapped_mom
    wet_to_file[:]       = wet_mom
    f.close() 


def read_l4_sss_esacci(ds_path):
    f = Dataset(ds_path)
    sss = np.squeeze( f.variables['sss'][:] )
    sss_qc = np.squeeze( f.variables['sss_qc'][:] )
    lsc_qc = np.squeeze( f.variables['lsc_qc'][:] )
    isc_qc = np.squeeze( f.variables['isc_qc'][:] )
    #delta_t  = f.variables['time'][:][0]
    t = num2date(f.variables['time'],f.variables['time'].units, \
                     calendar=f.variables['time'].calendar, \
                     only_use_cftime_datetimes=False, \
                     only_use_python_datetimes=True)[0]
    f.close()

    return sss, t, sss_qc, lsc_qc, isc_qc


def read_monthly_woa18_filled(ds_path,varname='SALT'):
    f = Dataset(ds_path)
    sss_12mn = f.variables[varname][:]
    f.close()

    return sss_12mn

def woa_month_to_day(woa_12mn, date):
    month_prev = date.month 
    month_next = date.month + 1 if date.month < 12 else 1
    imn_prev = month_prev - 1
    imn_next = month_next - 1

    ndays_total = calendar.monthrange(date.year,date.month)[1]
    wts_prev = 1.0 - (date.day-1)/ndays_total
    wts_next = 1.0 - wts_prev
    print("odate, wts_prev, wts_next=", date, wts_prev, wts_next)

    woa_daily = wts_prev*woa_12mn[imn_prev,:,:] + wts_next*woa_12mn[imn_next,:,:]

    return woa_daily    


def main(args):

    debug = False

    cfg = ConfigMerge()
    cfg.read_cfg(args.config_path, namelist="merge")
    
    #
    # load grids info 
    #

    # read in WOA18 lon, lat
    f = Dataset(cfg.woa_grid_path)
    lat1d_woa = f.variables[cfg.woa_lat_var][:]
    lon1d_woa = f.variables[cfg.woa_lon_var][:]
    f.close()

    # read in L4SSS lon, lat
    f = Dataset(args.l4_sss_path)
    lat1d_l4 = f.variables[cfg.l4_lat_var][:]
    lon1d_l4 = f.variables[cfg.l4_lon_var][:]
    f.close()

    # read in tripolar mom lon, lat
    f = Dataset(cfg.ocean_static_path)
    lat2d_mom = f.variables[cfg.mom_lat_var][:].astype(np.float32)
    lon2d_mom = f.variables[cfg.mom_lon_var][:].astype(np.float32)
    wet2d_mom = f.variables[cfg.mom_wet_var][:].astype(np.float32)
    f.close()

    # we only generate wts to WOA18 grids within the L4 grid range; otherwise, xesmf fails to work
    lat_max_allowed = lat1d_woa[lat1d_woa<lat1d_l4.max()][-1]
    idx_lat_max_allowed = np.where(lat1d_woa==lat_max_allowed)[0][0]

    lat_min_allowed = lat1d_woa[lat1d_woa>lat1d_l4.min()][0]
    idx_lat_min_allowed = np.where(lat1d_woa==lat_min_allowed)[0][0]
    lat1d_woa_allowed = lat1d_woa[idx_lat_min_allowed:idx_lat_max_allowed+1]

    print("min_WOA18_lat can be mapped from L4 ESA: idx_woa18, val_woa18=", idx_lat_min_allowed, lat1d_woa[idx_lat_min_allowed])
    print("max_WOA18_lat can be mapped from L4 ESA: idx_woa18, val_woa18=", idx_lat_max_allowed, lat1d_woa[idx_lat_max_allowed])

    #
    # load regridder
    #

    interp_method = "bilinear"
    periodic = True

    # load the regridder: L4SSS => WOA
    grd_in = {"lon": lon1d_l4, "lat": lat1d_l4}
    grd_out = {"lon": lon1d_woa, "lat": lat1d_woa[idx_lat_min_allowed:idx_lat_max_allowed+1]}
    fn_regridder_in = os.path.abspath(cfg.regridder_l4_to_woa_path)
    regridder_l4_to_woa = xe.Regridder(grd_in, grd_out, interp_method, weights=fn_regridder_in, periodic=periodic)
    regridder_l4_to_woa.filename = fn_regridder_in  # a bug of xesmf 
    if debug: print(regridder_l4_to_woa)

    # load the regridder: WOA => MOM
    grd_in = {"lon": lon1d_woa, "lat": lat1d_woa}
    grd_out = {"lon": lon2d_mom, "lat": lat2d_mom}
    fn_regridder_in = os.path.abspath(cfg.regridder_woa_to_mom_path)
    regridder_woa_to_mom = xe.Regridder(grd_in, grd_out, interp_method, weights=fn_regridder_in, periodic=periodic)
    regridder_woa_to_mom.filename = fn_regridder_in  # a bug of xesmf 
    if debug: print(regridder_woa_to_mom)
 

    #
    # interpolate L4SSS to WOA grids
    #

    # read in monthly WOA files
    sss_12mn_woa = read_monthly_woa18_filled(args.woa18_filled_path)
    if debug:
        for i in range(0,12):
            print(i, sss_12mn_woa[i,:,:].min(), sss_12mn_woa[i,:,:].max())
 
    # read in daily L4 SSS
    sss_l4, time_l4, sss_qc, lsc_qc, isc_qc = read_l4_sss_esacci(args.l4_sss_path)
    if debug:
        print("sss_l4: shape, min, max", sss_l4.shape, sss_l4.min(), sss_l4.max())
        print("time_l4", time_l4)
    all_qc = sss_qc + lsc_qc + isc_qc

    # calculate new time based on new time units from user input
    dtime_units, base_time = parse_time_units(cfg.new_time_units)
    new_dtime = (time_l4 - base_time)/dtime_units


    # create daily WOA18 SSS using monthly WOA18 (maskarray)
    sss_daily_woa = woa_month_to_day(sss_12mn_woa, time_l4)

    # interp L4 SSS to WOA18 grids using regridder (sss_l4 destroyed)
    # output is interped L4 SSS with land unfilled
    QC_GOOD = 0
    sss_l4[all_qc!=QC_GOOD] = np.nan
    sss_l4[sss_l4.mask] = np.nan
    sss_l4.mask = False
    sss_l4.fill_value = 0

    wk2d = np.empty((lat1d_woa.size,lon1d_woa.size))
    wk2d[:,:] = np.nan
    wk2d[idx_lat_min_allowed:idx_lat_max_allowed+1,:] = regridder_l4_to_woa(sss_l4)

    sss_unfilled = wk2d.copy()

    # fill the grids whose nearby have no L4SSS 
    wk2d = fill_nan_grds_aux(wk2d, sss_daily_woa.data, args.radius)
    sss_filled_step1 = wk2d.copy()

    # merge 2 datasets
    missing_value = 10.0e23
    fillmask = np.isnan(wk2d)
    wk2d[fillmask] = missing_value
    iterative_fill_POP_core(var=wk2d, fillmask=fillmask, missing_value=missing_value, tol=1.e-4, ltripole=False)

    sss_filled_step2 = wk2d.copy()

    if debug: print("range of merged salinity:", wk2d.min(), wk2d.max() )


    # 
    # interpolate to MOM tripolar grids
    #
    sss_remapped_mom = regridder_woa_to_mom(wk2d)
    ocn_pts_invalid = np.isnan(sss_remapped_mom) & (wet2d_mom>0.5)
    ocn_pts_valid = (~np.isnan(sss_remapped_mom)) & (wet2d_mom>0.5)
    if np.sum(ocn_pts_invalid) > 0:
        raise Exception("invalid ocean pts found: #={}".format(np.sum(ocn_pts_invalid)))
        sys.exit(1)

    print("salinity range after remapping: ", np.min(sss_remapped_mom[ocn_pts_valid]), np.max(sss_remapped_mom[ocn_pts_valid]))
   
    #
    # write out the daily restore file
    #
    tpl = RestoreTemplate()
    tpl.read(os.path.abspath(cfg.restore_template_path))
    input_dims, input_vars = tpl.get_restore_static(\
                               ocean_static_file = os.path.abspath(cfg.ocean_static_path), \
                               ocean_hgrid_file = os.path.abspath(cfg.ocean_hgrid_path))

    input_vars[cfg.sss_var_renamed] = np.expand_dims(sss_remapped_mom,axis=0)
    input_vars['time'] = np.array([new_dtime])
    tpl.create_restore(args.output_file_path, input_dims, input_vars)

    if args.debug_file_path is not None:
        write_nc_file(args.debug_file_path, lat1d_woa, lon1d_woa, lat2d_mom, lon2d_mom, \
                    sss_unfilled, sss_filled_step1, sss_filled_step2, sss_remapped_mom, wet2d_mom) 

if __name__ == '__main__':
    args = parseCommandLine()
    main(args)

