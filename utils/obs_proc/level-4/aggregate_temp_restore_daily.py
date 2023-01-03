#!/usr/bin/env python3

import datetime as dt
import argparse
import os, sys
import numpy as np
from netCDF4 import Dataset

def parseCommandLine():
    parser = argparse.ArgumentParser(description=("merge daily temperature restore file"))
    parser.add_argument("--template_file_path", required=True, help=(
                "MOM6 template file to which we write daily sst"))
    parser.add_argument("--start_date", required=True, help=("YYYYMMDDHH"))
    parser.add_argument("--end_date", required=True, help=("YYYYMMDDHH"))
    parser.add_argument("--daily_temp_dir", required=True, help=(
                "directory name to store daily sst file"))
    parser.add_argument("--daily_temp_file_name", required=True, help=(
                "something like temp_YYYYMMDD.nc where YYYYMMDD will be replaced with actual dates"))
    parser.add_argument("--output_file_name", required=False, type=str, help=(
                    "rename template file after finishing output"))
    args = parser.parse_args()

    args.start_date = dt.datetime.strptime(args.start_date,"%Y%m%d%H")
    args.end_date = dt.datetime.strptime(args.end_date,"%Y%m%d%H")


    args.template_file_path = os.path.abspath(args.template_file_path)
    if not os.path.exists(args.template_file_path):
        raise Exception("template_file_path ({}) does not exist".format(args.template_file_path))
        sys.exit(1)

    args.daily_temp_dir = os.path.abspath(args.daily_temp_dir)
    if not os.path.exists(args.daily_temp_dir):
        raise Exception("daily_temp_dir ({}) does not exist".format(args.daily_temp_dir))
        sys.exit(2)

    if args.output_file_name is not None:
        if os.path.exists(args.output_file_name):   
            raise Exception("output_file_name ({}) already exists".format(args.output_file_name))
            sys.exit(3)


    print(args)
    return args

def main(args):
    ndays = (args.end_date - args.start_date).days + 1
    cdate = args.start_date

    while cdate <= args.end_date:
        if cdate.day == 1: print("reading daily file: "+cdate.strftime("%Y%m%d"))

        temp_daily_file = args.daily_temp_file_name.replace("YYYYMMDD",cdate.strftime("%Y%m%d"))
        temp_daily_path = os.path.join(args.daily_temp_dir, temp_daily_file)
        if not os.path.exists(temp_daily_path):
            raise Exception("daily temp file ({}) does not exist".format(temp_daily_path))
            sys.exit(4)

        f = Dataset(temp_daily_path)
        wk2d = f.variables['sst_remapped_mom'][:]
        f.close()

        if cdate == args.start_date:
            nlat, nlon = wk2d.shape
            print("nlat,nlon=", nlat,nlon)
            ndays = (args.end_date - args.start_date).days + 1
            temp_aggr = np.zeros((ndays, nlat, nlon))
            day_aggr = np.zeros(ndays)
            base_date = dt.datetime(1900,1,1,0,0,0)
       
        idx_day = (cdate - args.start_date).days
        temp_aggr[idx_day,:,:] = wk2d
        day_aggr[idx_day] = (cdate - base_date).days

        cdate += dt.timedelta(days=1)

    f = Dataset(args.template_file_path,"r+")
    f.variables["SST"][:] = temp_aggr.astype(np.float32)
    f.variables["time"][:] = day_aggr.astype(np.float32)
    f.close()

    if args.output_file_name is not None:
        os.rename(args.template_file_path, args.output_file_name)


if __name__ == '__main__':
    args = parseCommandLine()
    main(args)
