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
varlist_b = ['DEPTH','POTM_OBS','PSAL_OBS','POTM_Hx0','PSAL_Hx0']

# Reference NEMOVAR qc'd feedback file:
qc_f  = '/home/rd/dasp/perm/DATA/OBS/00/0001_nrt_5_2_20161226_20170101_profbqc_01_fdbk.nc'
DO_QC = True

# Open qc'd file with index of observations to keep:
if (DO_QC):

  print('Opening file: ', qc_f)
  qc_fid = Dataset(qc_f, 'r')
  qc_oid = np.array(qc_fid.variables['ORIGINAL_FILE_INDEX'][:])
  qc_depth = np.array(qc_fid.variables['DEPTH'][:,:])

  qc_dim_a = len(qc_oid)
  print(qc_dim_a)

  # Sort new list
  sort_idx = np.argsort(qc_oid)
  print('original qc list (unsorted)')
  print(qc_oid[0:30])
  print('qc sort index:')
  print(sort_idx[0:30])
  print('original qc list (sorted)')
  print(qc_oid[sort_idx[0:30]])
  qc_fid.close()
  # Convert sort index to apply to ORIGINAL_FILE_INDEX:
  sort_idx0 = qc_oid[sort_idx]-1

else:
  print('Using raw observations (no qc from feedback files)')

# Loop through each observation file
idx=0
for nc_f in obsfiles:
  # Read in observation file (netcdf)
  nc_fid = Dataset(nc_f, 'r')  # Dataset is the class behavior to open the file
                               # and create an instance of the ncCDF4 class

  # Sort on the original file index
  oid = np.array(nc_fid.variables['ORIGINAL_FILE_INDEX'][:])
  sort_idx_i = np.argsort(oid)
  print('New file (unsorted):')
  print(oid[0:30])
  print('New file sort index:')
  print(sort_idx_i[0:30])

  if (DO_QC):
    print('Using qc sort instead of full dataset...')
    sort_idx = sort_idx_i[sort_idx0]
  else:
    sort_idx = sort_idx_i

  print('New file qc sort index:')
  print(sort_idx[0:30])

# quit()

  # Write out temporary replacement file that is properly sorted
  outfile = 'sorted'+'_'+str(idx).zfill(2)+'_'+ntpath.basename(nc_f)
  print(outfile)
  out_fid = Dataset(outfile, 'w', format='NETCDF4')

  print ('Cycle dimensions...')
  for name, dimension in nc_fid.dimensions.iteritems():
    print (dimension)
    if (name in ['N_OBS']):
#     print ('Setting N_OBS based on N_OBS from qc file...')
#     print (qc_dim_a)
      out_fid.createDimension(name, qc_dim_a if not dimension.isunlimited() else None)
    else:
      out_fid.createDimension(name, len(dimension) if not dimension.isunlimited() else None)

  print ('Cycle variables...')
  for name, variable in nc_fid.variables.iteritems():
    if any(name in s for s in varlist_a):
#     print ('------------------------------------------')
#     print(name)
#     print ('------------------------------------------')
      x = out_fid.createVariable(name, variable.datatype, variable.dimensions)
      x[:] = variable[sort_idx]
#     print(x)
#     print(variable[0:30])
#     print(x[0:30])
#     print(x.size)

    if any(name in s for s in varlist_b):
#     print ('------------------------------------------')
#     print(name)
#     print ('------------------------------------------')
      x = out_fid.createVariable(name, variable.datatype, variable.dimensions)
      x[:] = variable[sort_idx,:]
#     print(x)
#     print(x[0:10])
#     print(x.size)

  out_fid.close()
  nc_fid.close()
  idx += 1
