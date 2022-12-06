#/usr/bin/env python3
import argparse
import os, sys
import numpy as np
from netCDF4 import Dataset
from regrid_tools import iterative_fill_POP_core

def parseCommandLine():
    
    parser = argparse.ArgumentParser(description=(""))
    parser.add_argument("woa18_dir",type=str,help=(
                "directory of WOA18 monthly files such as woa18_decav_s11_04.nc"))
    parser.add_argument("--start_month",required=True,default=1,type=int,help=(
                ""))
    parser.add_argument("--end_month",required=True,default=12,type=int,help=(
                ""))
    args =parser.parse_args()
    print(args)
    return args

def write_nc_file(fnout, lat1d, lon1d, v2d_list, month_list, v2d_name="unfilled"):
    print("write out netcdf file: {}".format(fnout))
    if os.path.exists(fnout):
        raise Exception("fnout already exists at: {}".format(fnout))
        sys.exit(1)

    fnout = os.path.abspath(fnout)
    f = Dataset(fnout, mode='w', format='NETCDF4_CLASSIC')
    f.createDimension('lat', lat1d.size)
    f.createDimension('lon', lon1d.size)
    f.createDimension('time', len(v2d_list))

    lat_to_file = f.createVariable('lat',np.float32, ('lat',))
    lon_to_file = f.createVariable('lon',np.float32, ('lon',))
    v2d_to_file = f.createVariable(v2d_name,np.float32, ('time','lat','lon'))
    month_to_file = f.createVariable('month',np.float32,('time'))

    lat_to_file[:] = lat1d
    lon_to_file[:] = lon1d

    v2d_as_array = np.array(v2d_list)
    month_as_array = np.array(month_list,dtype=np.float32)

    v2d_to_file[:,:,:] = v2d_as_array
    month_to_file[:] = month_as_array

    f.close()


def main(args):

    #
    # read in the lat & lon of the input WOA18 fields, & generate unfilled monthly file
    #

    sss_unfilled_list = []
    month_list = []
    for month in range(args.start_month, args.end_month+1):
        fin_name = 'woa18_decav_s{:02d}_04.nc'.format(month)
        fin_path = os.path.join(args.woa18_dir, fin_name)
        print("read in monthly file: {}".format(fin_path))
        if not os.path.exists(fin_path):
            raise Exception("WOA monthly file does not exist at: {}".format(fin_path))
            sys.exit(2)
        f = Dataset(fin_path)
        sss_grd_in = np.squeeze(f.variables['s_an'][:])[0,:,:] # use the 1st layer
        if month == args.start_month:
            lat1d_grd_in = f.variables['lat'][:]
            lon1d_grd_in = f.variables['lon'][:]
            depth_grd_in = f.variables['depth'][:]

        sss_unfilled_list.append(sss_grd_in.copy())
        month_list.append(month)
        f.close()

    print("month_list=",month_list)
    print("sss_unfilled_list=",len(sss_unfilled_list))
    fout_name_unfilled = "unfilled_sss_s{:02d}_e{:02d}.nc".format(args.start_month,args.end_month)
    fout_name_filled = "filled_sss_s{:02d}_e{:02d}.nc".format(args.start_month,args.end_month)
   
    # write out monthly-aggregated file out with missing values on land
    write_nc_file(fout_name_unfilled, lat1d_grd_in, lon1d_grd_in, sss_unfilled_list,month_list)

    #
    # fill the land
    #
    sss_filled_list = []
    for month, sss_unfilled in zip(month_list, sss_unfilled_list):
        print("fill month: {}".format(month))
        wk2d = sss_unfilled.data
        missing_value = sss_unfilled.fill_value
        #wk2d[sss_unfilled.mask] = missing_value

        fillmask = sss_unfilled.mask
        iterative_fill_POP_core(var=wk2d, fillmask=fillmask, missing_value=missing_value, tol=1.e-4, ltripole=True)
        sss_filled_list.append(wk2d.copy())

    # write out monthly-aggregated file out with extrapolation value on land
    write_nc_file(fout_name_filled, lat1d_grd_in, lon1d_grd_in, sss_filled_list, month_list)



if __name__ == '__main__':
    args = parseCommandLine()
    main(args)
