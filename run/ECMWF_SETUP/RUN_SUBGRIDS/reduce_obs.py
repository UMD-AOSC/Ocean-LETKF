import datetime as dt        # Python standard library datetime  module
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import ntpath
import sys

# Author - Professor Stephen G. Penny
#          University of Maryland, College Park
#          ECMWF Visiting Scientist 1 June - 28 Nov, 2017
#          contact: Steve.Penny@noaa.gov
#
# Purpose:
#          Sort the observations based on the ORIGINAL_FILE_INDEX.
#          Apply quality control by removing any observations discarded by the QC
#          operations that have been applied to the control member.
#          Replace observation locations of each member with the location of the
#          control member (in case perturbations have been applied to the obs locations).

def fill_idx (sort_idx_i,oid):
  # Check to make sure there are no missing indices, if so, fill them in with a -1 flag
  # So the ORIGINAL_FILE_INDEX and sort_idx_i match correctly.
  j = 0
  for i in xrange(0,max(oid)-1):
    if (j<oid[sort_idx_i[i]]-1):
      # If there is a missing index, put a -1
#     print('Observation ORIGINAL_FILE_INDEX number oid[sort_idx_i[i]]=',i+1,' is missing')
#     print('i = ',i)
#     print('j = ', j)
#     print('Pre-insert')
#     print('sort_idx_i[i-2] = ', sort_idx_i[i-2])
#     print('sort_idx_i[i-1] = ', sort_idx_i[i-1])
#     print('sort_idx_i[i]   = ', sort_idx_i[i])
#     print('sort_idx_i[i+1] = ', sort_idx_i[i+1])
#     print('sort_idx_i[i+2] = ', sort_idx_i[i+2])
#     print('oid[sort_idx_i[i-2]] = ', oid[sort_idx_i[i-2]])
#     print('oid[sort_idx_i[i-1]] = ', oid[sort_idx_i[i-1]])
#     print('oid[sort_idx_i[i]]   = ', oid[sort_idx_i[i]])
#     print('oid[sort_idx_i[i+1]] = ', oid[sort_idx_i[i+1]])
#     print('oid[sort_idx_i[i+2]] = ', oid[sort_idx_i[i+2]])
#     print('len(sort_idx_i) = ', len(sort_idx_i))

      # (1) Put an missing value in oid
#     print('Inserting ',i+1,' before oid[sort_idx_i[i]] = ', oid[sort_idx_i[i]])
      oid = np.insert(oid,sort_idx_i[i],i+1)

#     print('Post-insert (a)')
#     print('oid[sort_idx_i[i-2]] = ', oid[sort_idx_i[i-2]])
#     print('oid[sort_idx_i[i-1]] = ', oid[sort_idx_i[i-1]])
#     print('oid[sort_idx_i[i]]   = ', oid[sort_idx_i[i]])
#     print('oid[sort_idx_i[i+1]] = ', oid[sort_idx_i[i+1]])
#     print('oid[sort_idx_i[i+2]] = ', oid[sort_idx_i[i+2]])
#     print('len(sort_idx_i) = ', len(sort_idx_i))

      # (2) increment the sort index for all remaining points
#     print('Inserting ',sort_idx_i[i],' before sort_idx_i[i] = ', sort_idx_i[i])
      sort_idx_i = np.insert(sort_idx_i,sort_idx_i[i],sort_idx_i[i])
      sort_idx_i[i+1:] = sort_idx_i[i+1:]+1

#     print('Post-insert (b)')
#     print('sort_idx_i[i-2] = ', sort_idx_i[i-2])
#     print('sort_idx_i[i-1] = ', sort_idx_i[i-1])
#     print('sort_idx_i[i]   = ', sort_idx_i[i])
#     print('sort_idx_i[i+1] = ', sort_idx_i[i+1])
#     print('sort_idx_i[i+2] = ', sort_idx_i[i+2])
#     print('len(sort_idx_i) = ', len(sort_idx_i))
    j=j+1

  return sort_idx_i, oid

# Input a list of files
obsfilelist = 'observation_infile_list.txt'
olist_fid = open(obsfilelist,'r')
obsfiles = olist_fid.read().splitlines()
print(obsfiles)

# Variables to keep, dimensioned a: (nobs), and b: (nobs,depth)
varlist_a = ['LONGITUDE','LATITUDE','JULD','POTM_QC','PSAL_QC','ORIGINAL_FILE_INDEX']
varlist_b = ['DEPTH','POTM_OBS','PSAL_OBS','POTM_Hx','PSAL_Hx']

# Reference NEMO qc'd feedback file:
qc_f  = sys.argv[1] #'/home/rd/dasp/perm/DATA/OBS/00/0001_nrt_5_2_20161226_20170101_profbqc_01_fdbk.nc'
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
  print(sort_idx[-30:])
  print('original qc list (sorted)')
  print(qc_oid[sort_idx[0:30]])
  print(qc_oid[sort_idx[-30:]])
  qc_fid.close()
  # Convert sort index to apply to ORIGINAL_FILE_INDEX:
  sort_idx0 = qc_oid[sort_idx]-1
  print('sort_idx0 = ')
  print(sort_idx0[0:30])
  print(sort_idx0[-30:])

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

# sort_idx_i, oid = fill_idx(sort_idx_i,oid)

  print('New file (unsorted):')
  print(oid[0:30])
  print(oid[-30:])
  print('New file sort index:')
  print(sort_idx_i[0:30])
  print(sort_idx_i[-30:])

  if (DO_QC):
    # ORIGINAL_FILE_INDEX should not exceed the number of observations in the control file to use this program.
    # (ISSUE)
    print('Using qc sort instead of full dataset...')
    sort_idx = sort_idx_i[sort_idx0]
  else:
    sort_idx = sort_idx_i

  print('New file qc sort index:')
  print(sort_idx[0:30])
  print(sort_idx[-30:])

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
      print ('------------------------------------------')
      print(name)
      print ('------------------------------------------')
      x = out_fid.createVariable(name, variable.datatype, variable.dimensions)
      print('sort_idx = ')
      print(sort_idx[0:30])
      print(sort_idx[-30:])
      print('len(sort_idx) = ', len(sort_idx))
      print('variable[:] = ', variable[:])
      print('len(variable[:]) = ', len(variable[:]))

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
