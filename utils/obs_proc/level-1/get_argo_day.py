#!/usr/bin/env python3

"""
To download global daily files of ARGO data from several servers.
This script works on Argo data python library
Documentation: https://argopy.readthedocs.io/en/latest/#
Requires module with the argopy library. Ej:
    module use -a ~sakella/modulefiles/
    module load python/argopy-22Feb2023
"""

from argopy import DataFetcher as ArgoDataFetcher
import pandas as pd
import argparse
import sys
import os
from tqdm import tqdm

def parseCommandLine():

    parser = argparse.ArgumentParser(description=("Script to download daily global ARGO data"))
    parser.add_argument("-s", "--start_date", default="20190101", metavar="YYYY-MM-DD",required=True, help=("Start date"))
    parser.add_argument("-e", "--end_date", default=None, metavar="YYYY-MM-DD", required=False, help=("End date"))
    parser.add_argument("-d", "--domain", default=[-180, 180, -90, 90], metavar="[lon_l, lon_r, lat_b, lat_t]", required=False, help=("domain [lon_left, lon_right, lat_bot, lat_top]"))
    parser.add_argument("-p", "--pres", default=[0, 6000], metavar="[p0, ..]", required=False, help=("Pressure levels"))
    parser.add_argument("-o", "--outdir", default="./", required=False, help=("Output directory. Data will be stored under outdir/YYYY/MM/"))
    parser.add_argument("--show", action="store_true", required=False, help=("display all the input arguments w/o running the script"))

    args = parser.parse_args()

    args.start_date = pd.to_datetime(args.start_date)

    if args.end_date is None:
        args.end_date = args.start_date
    else:
        args.end_date = pd.to_datetime(args.end_date)

    if args.end_date < args.start_date:
        raise Exception("end_date ({}) should be after start_date ({}).".format(args.end_date, args.start_date))
        sys.exit(1)

    args.outdir = os.path.abspath(args.outdir)

    if args.show:
        print(args)
        sys.exit(0)

    return args

argo_loader_p=ArgoDataFetcher(parallel=True)

def main(args):

    cdate = args.start_date

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # n = max((args.end_date - args.start_date).days, 1)

    # pbar = tqdm(total=n)

    while cdate <= args.end_date:
        edate = cdate + pd.to_timedelta('1 d')
        cdir = args.outdir + "/" + cdate.strftime("%Y/%m/")
        
        if not os.path.exists(cdir):
            os.makedirs(cdir)
    
        fpath = cdir+'argo_global_{}.nc'.format(cdate.strftime('%Y-%m-%d'))

        if not os.path.isfile(fpath):
            roi = args.domain + args.pres + [cdate.strftime('%Y-%m-%d'), edate.strftime('%Y-%m-%d')]
            try:
                ds = argo_loader_p.region(roi).to_xarray()
                ds.to_netcdf(fpath)
            except Exception as e:
                print(f'caught {type(e)}: {e}')
        # pbar.update(1)
        cdate = edate



if __name__== "__main__":
    args = parseCommandLine()
    main(args)
