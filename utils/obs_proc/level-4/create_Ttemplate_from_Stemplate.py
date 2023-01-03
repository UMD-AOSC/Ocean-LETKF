from netCDF4 import Dataset

f = Dataset("temp_restore_template.nc","r+")
f.renameVariable("SALT","TEMP")
f['TEMP'].long_name = "sea surface temperature"
f['TEMP'].units = "deg C"
f.close()
