#!/usr/bin/env python3
import sys, os
import common
import subprocess as sp
import argparse

def parseCommandLine():

    parser = argparse.ArgumentParser(description=("script to download WOA18 datasets"))
    parser.add_argument("outdir",type=str,help=("where to download the WOA18 data"))

    parser.add_argument("--dvar", choices=["salinity","temperature"], required=True, help=("var to download"))
    parser.add_argument("--start_month", type=int, required=True, help=("start month to download"))
    parser.add_argument("--end_month", type=int, required=True, help=("end month to download"))
    parser.add_argument("--res", choices=["5","1","0.25"], required=True, help=("resolution in deg (options: 5, 1, 0.25)"))

    parser.add_argument("--annual", action="store_true",required=False, help=("get annual mean"))
    parser.add_argument("--seasonal", action="store_true",required=False, help=("get seasonal mean"))
    parser.add_argument("--verbose",action="store_true",required=False,help=("print out wget command"))

    args = parser.parse_args()
    print(args)
    return args


def check_wget(logging):
    # check if wget is available in the current system
    cmd_wget = 'wget --version'
    proc = sp.Popen(cmd_wget, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = proc.communicate()
    if proc.returncode != 0:
         errmsg = "wget not available on the current system"
         logging.critical("{}\n\nSTDOUT: {}\nERR: {}\n".format(errmsg, out.decode(), err.decode()))
         sys.exit(1)

def main(args):

    logging = common.setupEasyLogging("get_monthly_woa18",fileLog=False)

    # check if wget exists
    check_wget(logging)

    # check if the output directory exists. Create it if not
    if not os.path.exists(args.outdir):
        logging.info("directory ({}) not found. create it".format(args.outdir))
        os.makedirs(args.outdir)

    # start to download required data
    month_list = list(range(args.start_month, args.end_month+1))
    if args.seasonal: month_list += [13,14,15,16]
    if args.annual:   month_list += [0]
    for month in month_list:
        if args.dvar == 'salinity':
            dvar_s = 's'
        elif args.dvar == 'temperature':
            dvar_s = 't'

        if args.res == '5':
            dres  = '5deg'
            dres_s = '5d'
        elif args.res == '1':
            dres = '1.00'
            dres_s = '01'
        elif args.res == '0.25':
            dres = '0.25'
            dres_s = '04'
        dname = 'woa18_decav_{}{:02d}_{}.nc'.format(dvar_s, month, dres_s)
        cmd_get = "wget -c https://www.ncei.noaa.gov/thredds-ocean/fileServer/ncei/woa/{}/decav/{}/{}".format(args.dvar, dres, dname)
        if args.verbose: logging.debug(cmd_get)

        proc = sp.Popen(cmd_get, shell=True, cwd=args.outdir, stdout=sp.PIPE, stderr=sp.PIPE)
        out, err = proc.communicate()
        if proc.returncode != 0:
            errmsg = "encounter errors when downloading data from remote server"
            logging.critical("{}\nSTDOUT: {}\nERR: {}\n".format(errmsg, out.decode(), err.decode()))
            sys.exit(2)
        logging.info("finish downloading {}".format(dname))
   
    logging.info("FINISHED. YEAH!!!")


if __name__ == '__main__':
    args = parseCommandLine()
    main(args)
