import datetime as dt        # Python standard library datetime  module
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import ntpath

# STEVE: use better way to input list of files:
#obsfiles=['../OBS/00/0001_nrt_5_2_20161226_20170101_profb_01_fdbk.nc']
#obsfiles=['../OBS/00/0001_nrt_5_2_20161226_20170101_profb_01_fdbk.nc','../OBS/01/0001_nrt_5_2_20161226_20170101_profb_01_fdbk.nc','../OBS/02/0001_nrt_5_2_20161226_20170101_profb_01_fdbk.nc','../OBS/03/0001_nrt_5_2_20161226_20170101_profb_01_fdbk.nc','../OBS/04/0001_nrt_5_2_20161226_20170101_profb_01_fdbk.nc']

# Input a list of files
obsfilelist = 'observation_infile_list.txt'
olist_fid = open(obsfilelist,'r')
obsfiles = olist_fid.read().splitlines()
print(obsfiles)

# Variables to keep, dimensioned a: (nobs), and b: (nobs,depth)
varlist_a = ['LONGITUDE','LATITUDE','JULD','POTM_QC','PSAL_QC','ORIGINAL_FILE_INDEX']
varlist_b = ['DEPTH','POTM_OBS','PSAL_OBS','POTM_Hx','PSAL_Hx']

# Loop through each observation file
idx=0
for nc_f in obsfiles:
  # Read in observation file (netcdf)
  nc_fid = Dataset(nc_f, 'r')  # Dataset is the class behavior to open the file
                               # and create an instance of the ncCDF4 class

  # Sort on the original file index
  oid = np.array(nc_fid.variables['ORIGINAL_FILE_INDEX'][:])
  sort_idx = np.argsort(oid)
  print(sort_idx[0:20])

  # Write out temporary replacement file that is properly sorted
  outfile = 'sorted'+'_'+str(idx).zfill(2)+'_'+ntpath.basename(nc_f)
  print(outfile)
  out_fid = Dataset(outfile, 'w', format='NETCDF4')
  
  for name, dimension in nc_fid.dimensions.iteritems():
    out_fid.createDimension(name, len(dimension) if not dimension.isunlimited() else None)

  for name, variable in nc_fid.variables.iteritems():

    if any(name in s for s in varlist_a):
      x = out_fid.createVariable(name, variable.datatype, variable.dimensions)
      print(name)
      x[:] = variable[sort_idx]
      print(x)

    if any(name in s for s in varlist_b):
      x = out_fid.createVariable(name, variable.datatype, variable.dimensions)
      print(name)
      x[:] = variable[sort_idx,:]
      print(x)

  out_fid.close()
  nc_fid.close()
  idx += 1
