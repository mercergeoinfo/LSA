#!/Users/andrew/anaconda/bin/python
import spatialfiles
dataFile = '/Users/andrew/Google Drive/Work/MeltModel/Output/Surface/2010/melt_e.csv'
DEM = '../InData/DEM/SG_DEM_2010.tif'
output = spatialfiles.idw(dataFile, DEM, epsg='3006')