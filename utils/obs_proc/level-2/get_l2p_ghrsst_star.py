#!/usr/bin/env python3

"""
Example:

PYTHONPATH=../../pycommon ./get_l2p_ghrsst_star.py --start_date 20220110 --end_date 20220111 --sis ABI_G17 --outdir ./ --topdir_name test_l2p_sst --skip_remote_check  --get_nav_file --verbose --log

"""

import common
import argparse
import datetime as dt
import subprocess as sp
import os, sys

def get_sis_info(sis):
    sis_dict={"ABI_G16":{"sat":"GOES16", "ver":"v2.70", "geodir":"docs", "geofile":"G16_075_0_W.nc"}, 
              "ABI_G17":{"sat":"GOES17", "ver":"v2.71", "geodir":"nav", "geofile":"G17_137_0_W.nc"},
              "H8":     {}
    }
    if sis in sis_dict.keys():
        info = sis_dict[sis]
        ret  = 0
    else:
        info = None
        ret  = 255

    return info, ret

parser = argparse.ArgumentParser(description=("Script to download STAR-L2P_GHRSST data"),
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--sis", default="ABI_G16", metavar="SENSOR_SATELLITE", required=True,help=("sensor/satellite name"))
parser.add_argument("--start_date",default="20220803", metavar="YYYYMMDDHH", required=True,help=("start date of data"))
parser.add_argument("--end_date",  default=None, metavar="YYYYMMDDHH", required=False,help=("end date of data"))
parser.add_argument("--outdir", default="./", required=False, help=("output directory. downloaded data will be stored under outdir/topdir_name/sis/YYYY/YYYYMM/YYYYMMDD"))


parser.add_argument("--topdir_name", default="l2p_ghrsst_star", required=False, help=("top directory name for saving all data"))
parser.add_argument("--get_nav_file", action="store_true", required=False, help=("get constant navigation .nc file (e.g., lat, long)"))
parser.add_argument("--show", action="store_true", required=False, help=("display all the input arguments w/o running the scripts"))
parser.add_argument("--skip_remote_check", action="store_true",required=False, help=("skip file checks in the remote server"))
parser.add_argument("--verbose", action="store_true", required=False, help=("show more debug info"))
parser.add_argument("--log", action="store_true", required=False, help=("generate a text log file"))

args = parser.parse_args()

args.start_date = dt.datetime.strptime(args.start_date,"%Y%m%d")
if args.end_date is None:
    args.end_date = args.start_date
else:
    args.end_date = dt.datetime.strptime(args.end_date,"%Y%m%d")

args.outdir = os.path.abspath(args.outdir)

if args.show:
    print(args)
    sys.exit(0)

if args.log:
    logging = common.setupEasyLogging("get_l2p_ghrsst_star")
else:
    logging = common.setupEasyLogging("get_l2p_ghrsst_star",fileLog=False)
logging.info(" ".join(sys.argv))

# check if the sensor mapping
sis_info, returncode = get_sis_info( args.sis )
if returncode != 0:
    logging.critical("Fail to find the corresponding satellite for sis={}".format(args.sis))
    sys.exit(4)
if args.verbose:
    logging.debug("sis_info = {}".format(sis_info))

# check if wget available in the current system
cmd_wget = 'wget --version'
proc = sp.Popen(cmd_wget, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
out, err = proc.communicate()
if proc.returncode != 0:
    errmsg = "wget not available on the current system"
    logging.critical("{}\n\nSTDOUT: {}\nERR: {}\n".format(errmsg, out.decode(), err.decode()))
    sys.exit(1)

if args.get_nav_file:
    navdir = args.outdir + "/" + args.topdir_name + "/" + args.sis + "/nav"
    if args.verbose:
        logging.debug(navdir)
    if not os.path.exists(navdir):
        logging.info("navigation directory ({}) not found. create it".format(navdir))
        os.makedirs(navdir)
 
    cmd_get_nav='wget -r -l1 -np -nH --cut-dirs 20 "https://opendap.jpl.nasa.gov/opendap/hyrax/allData/ghrsst/data/GDS2/L2P/{}/STAR/{}/" -A "*.nc" -P ./ -e robots=off'.format(sis_info["sat"], sis_info["geodir"] )
    if args.verbose:
        logging.debug(cmd_get_nav)
    proc = sp.Popen(cmd_get_nav, shell=True, cwd=navdir, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = proc.communicate()
    if proc.returncode != 0:
        errmsg = "fail to get navigation *.nc file for sis=".format(args.sis)
        logging.critical("{}\n\nSTDOUT: {}\nERR: {}\n".format(errmsg, out.decode(), err.decode()))
        sys.exit(5)
 
# loop over inquired days
cdate = args.start_date
while cdate <= args.end_date:
    logging.info("start downloading {}".format(cdate.strftime("%Y-%m-%d")))
    julian_day = cdate.timetuple().tm_yday


    if not args.skip_remote_check:
        # check if remote files exist. Raise error if remote files do not exist
        cmd_inq = 'wget -r -l1 -np -nH --spider --cut-dirs 20 "https://opendap.jpl.nasa.gov/opendap/hyrax/allData/ghrsst/data/GDS2/L2P/{}/STAR/{}/{:04d}/{:03d}/" -A "{:04d}*.nc" -P ./ -e robots=off'.format(sis_info["sat"], sis_info["ver"], cdate.year, julian_day, cdate.year)

        if args.verbose:
            logging.debug(cmd_inq)
        proc = sp.Popen(cmd_inq, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
        out, err = proc.communicate()
        if proc.returncode != 0:
            errmsg = "remote files on this day do not exist"
            logging.critical("{}\n\nSTDOUT: {}\nERR: {}\n".format(errmsg, out.decode(), err.decode()))
            sys.exit(2)
       
    # create local dirs to store daily downloads 
    cdir = args.outdir + "/" + args.topdir_name + "/" + args.sis + cdate.strftime("/%Y/%Y%m/%Y%m%d/")
    if args.verbose:
        logging.debug(cdir)
    if not os.path.exists(cdir):
        logging.info("directory ({}) not found. create it".format(cdir))
        os.makedirs(cdir)

    # start daily download
    cmd_get = 'wget -r -l1 -np -nH --cut-dirs 20 "https://opendap.jpl.nasa.gov/opendap/hyrax/allData/ghrsst/data/GDS2/L2P/{}/STAR/{}/{:04d}/{:03d}/" -A "{:04d}*.nc" -P ./ -e robots=off'.format(sis_info["sat"], sis_info["ver"], cdate.year, julian_day, cdate.year)
    if args.verbose:
        logging.debug(cmd_get)
    proc = sp.Popen(cmd_get, shell=True, cwd = cdir, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = proc.communicate()
    if proc.returncode != 0:
        errmsg = "encounter errors when downloading data from remote server"
        logging.critical("{}\nSTDOUT: {}\nERR: {}\n".format(errmsg, out.decode(), err.decode()))
        sys.exit(3)

    logging.info("finish downloading {}".format(cdate.strftime("%Y-%m-%d")))
    cdate += dt.timedelta(days=1)

logging.info("FINISHED. YEAH!!!")
