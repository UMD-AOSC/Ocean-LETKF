#!/usr/bin/env python3

from netCDF4 import Dataset

f = Dataset("temp_restore_template_p25.nc","r+")
f.renameVariable("SALT","SST")
f['SST'].long_name = "sea surface temperature"
f['SST'].units = "degC"
f.close()
