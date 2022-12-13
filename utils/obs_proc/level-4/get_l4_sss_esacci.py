#!/usr/bin/env python3

"""
ftp://anon-ftp.ceda.ac.uk/neodc/esacci/sea_surface_salinity/data/v03.21/7days//2019/ESACCI-SEASURFACESALINITY-L4-SSS-MERGED_OI_7DAY_RUNNINGMEAN_DAILY_25km-20190101-fv3.21.nc
"""

import common
import argparse
import datetime as dt
import subprocess as sp
import os, sys


def check_wget(logging):
    """
        check if wget is available in the current system
    """
    cmd_wget = 'wget --version'
    proc = sp.Popen(cmd_wget, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = proc.communicate()
    if proc.returncode != 0:
         errmsg = "wget not available on the current system"
         logging.critical("{}\n\nSTDOUT: {}\nERR: {}\n".format(errmsg, out.decode(), err.decode()))
         sys.exit(1)


def parseCommandLine():

    parser = argparse.ArgumentParser(description=("Script to download ESACCI 50km daily SSS"))
    parser.add_argument("--start_date", default="20190101", metavar="YYYYMMDD",required=True, help=("start date"))
    parser.add_argument("--end_date", default="20190101", metavar="YYYYMMDD", required=False, help=("end date"))
    parser.add_argument("--outdir", default="./", required=False, help=("output directory. downloaded data will be stored under outdir/topdir_name/YYYY/YYYYMM/YYYYMMDD"))
    parser.add_argument("--topdir_name", default="l4_sss_esacci", required=False, help=("top subdirectory name under outdir"))
    parser.add_argument("--show", action="store_true", required=False, help=("display all the input arguments w/o running the script"))
    parser.add_argument("--log", action="store_true", required=False, help=("generate a text log file"))
    parser.add_argument("--verbose", action="store_true", required=False, help=("print more diag info"))

    args = parser.parse_args()

    args.start_date = dt.datetime.strptime(args.start_date,"%Y%m%d")
    if args.end_date is None:
        args.end_date = args.start_date
    else:
        args.end_date = dt.datetime.strptime(args.end_date,"%Y%m%d")
    if args.end_date < args.start_date:
        raise Exception("end_date ({}) should be after start_date ({}).".format(args.end_date, args.start_date))
        sys.exit(1)

    args.outdir = os.path.abspath(args.outdir)

    if args.show:
        print(args)
        sys.exit(0)

    return args


def main(args):

    # initialize log file
    if args.log:
        logging = common.setupEasyLogging("get_l4_sss_esacci")
    else:
        logging = common.setupEasyLogging("get_l4_sss_esacci", fileLog=False)
    logging.info(" ".join(sys.argv))

    # check if wget exists
    check_wget(logging)

    # loop over days
    cdate = args.start_date
    while cdate <= args.end_date:
        logging.info("start downloading {}".format(cdate.strftime("%Y-%m-%d")))

        # create local directory
        cdir = args.outdir + "/" + args.topdir_name + cdate.strftime("/%Y/%Y%m/%Y%m%d/")
        if args.verbose:
            logging.debug(cdir)
        if not os.path.exists(cdir):
            if args.verbose: logging.debug("directory ({}) not found. create it.".format(cdir))
            os.makedirs(cdir)

        # start to download files
        # example: ftp://anon-ftp.ceda.ac.uk/neodc/esacci/sea_surface_salinity/data/v03.21/7days//2019/ESACCI-SEASURFACESALINITY-L4-SSS-MERGED_OI_7DAY_RUNNINGMEAN_DAILY_25km-20190101-fv3.21.nc
        cmd_get = 'wget -c ftp://anon-ftp.ceda.ac.uk/neodc/esacci/sea_surface_salinity/data/v03.21/7days/{}/ESACCI-SEASURFACESALINITY-L4-SSS-MERGED_OI_7DAY_RUNNINGMEAN_DAILY_25km-{}-fv3.21.nc'.format(cdate.year,cdate.strftime("%Y%m%d"))
        if args.verbose:
            logging.debug(cmd_get)
        proc = sp.Popen(cmd_get, shell=True, cwd=cdir, stdout=sp.PIPE, stderr=sp.PIPE)
        out, err = proc.communicate()
        if proc.returncode != 0:
            errmsg = "encounter errors when downloading data from remote server"
            logging.critical("{}\nSTDOUT: {}\nERR: {}\n".format(errmsg, out.decode(), err.decode()))
            sys.exit(3)
        
        logging.info("finish downloading {}".format(cdate.strftime("%Y-%m-%d")))
        cdate +=dt.timedelta(days=1)

    logging.info("FINISHED. YEAH!!!")


if __name__ == '__main__':
    args = parseCommandLine()
    main(args)

