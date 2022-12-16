#!/usr/bin/env python3

import os, sys
import subprocess as sp
import datetime as dt

sdate = dt.datetime(2019,6,2,0,0,0)
edate = dt.datetime(2019,12,31,0,0,0)
radius = 10

cdate = sdate
while cdate <= edate:
    print("-"*80+"\n"+cdate.strftime("%Y-%m-%d"))
    print("time now:", dt.datetime.now())

    output_dir = "remapped_l4_sss"
    output_file_name = "sss_remapped_{}_radius{:03d}.nc".format(cdate.strftime("%Y%m%d"), radius)
    output_path  = output_dir+"/"+output_file_name

    l4_sss_daily_dir = "l4_sss_esacci/{}".format(cdate.strftime("%Y/%Y%m/%Y%m%d"))
    l4_sss_daily_file_name = "ESACCI-SEASURFACESALINITY-L4-SSS-MERGED_OI_7DAY_RUNNINGMEAN_DAILY_25km-{}-fv3.21.nc".format(cdate.strftime("%Y%m%d"))
    l4_sss_path = l4_sss_daily_dir + "/" + l4_sss_daily_file_name

    cmd = \
    f""" time PYTHONPATH=../../pycommon ./merge_l4_woa18.py \
        --regridder_l4_to_woa_path wts_l4sss_to_woa18.nc \
        --regridder_woa_to_mom_path wts_woa18_to_mom_p25.nc \
        --woa_grid_path ../../test_data/woa18_decav_s01_04.nc \
        --l4_grid_path {l4_sss_path} \
        --mom_grid_path ../../test_data/flood_salt_restore_PHC2.1440x1080.v20180405.nc \
        --mom_lsmsk_path ../../test_data/ocean_topog.nc \
        --woa18_filled_path ../woa/sss_woa18_10m/filled_sss_s01_e12.nc \
        --l4_sss_path {l4_sss_path} \
        --radius {radius} \
        --output_file_name {output_path}
    """
    print(cmd)
    proc = sp.Popen(cmd, shell=True, cwd="./", stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = proc.communicate()
    print(out.decode())
    if proc.returncode != 0:
        raise Exception("Fail to generate daily restore file ERROR: {}\n".format(err.decode()))
        sys.exit(3)



    cdate += dt.timedelta(days=1)
print("finished")
