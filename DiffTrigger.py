#!/Users/andrew/anaconda/bin/python
import AM_Func as AM
import os
Afile = '/Users/andrew/Google Drive/Work/MeltModel/Edit/2012mb/sg_mb_obsbw_Bw/sg_mb_obsbw_Bw_kriged_data.tif'
Bfile = '/Users/andrew/Google Drive/Work/MeltModel/Edit/2012mb/sg_mb_modbw_Bw/sg_mb_modbw_Bw_kriged_data.tif'
outfolder = '/Users/andrew/Google Drive/Work/MeltModel/Edit/'

AM.rasterDiff(Afile, Bfile, outfolder)