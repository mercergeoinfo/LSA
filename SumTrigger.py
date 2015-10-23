#!/Users/andrew/anaconda/bin/python
import AM_Func as AM
import os
filename = '../Edit/2012mb/sg_mb_modbw_Bw/sg_mb_modbw_Bw_kriged_data.tif'
results = AM.geopixsum(filename)
gdalcom = 'gdalinfo -stats ' + filename
os.system(gdalcom)
