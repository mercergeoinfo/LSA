#!/Users/andrew/anaconda/bin/python
import spatialfiles as AM
import os
Afile = '/Users/andrew/Google Drive/Work/MeltModel/Output/Surface/2012/Kriging/melt_e_Org_Bw/melt_e_Org_Bw_kriged_data.tif'
Bfile = '/Users/andrew/Google Drive/Work/MeltModel/Output/Surface/2012/Kriging/melt_e_Mod1_Bw_186/melt_e_Mod1_Bw_186_kriged_data.tif'
outfolder = '/Users/andrew/Google Drive/Work/MeltModel/Output/Surface/2012/Kriging/'

AM.rasterDiff(Afile, Bfile, outfolder)


#gdalwarp -s_srs EPSG:3006 -t_srs EPSG:3006 -te 647984.320 7535662.820 651339.320 7537607.820 -tr 5 5 -r near -dstnodata -9999 -of GTiff -cutline  /Users/andrew/Google\ Drive/Work/MeltModel/InData/DEM/Outline_2010.shp /Users/andrew/Documents/Work/TRSGeoData/StorglacData/massbalance/mbgrids_sweref99/2010bw.tif TRS_Bw_2010.tif

