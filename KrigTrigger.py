#!/Users/andrew/anaconda/bin/python
import AM_Func as AM
dataFile = '/Users/andrew/Google Drive/Work/MeltModel/Output/Surface/2012b2/melt_e.csv'
DEM = '../InData/DEM/SG_DEM_2010.tif'
AM.kriging(dataFile,DEM)