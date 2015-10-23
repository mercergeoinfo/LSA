#!/Users/andrew/anaconda/bin/python
import AM_Func as AM
dataFile = '../Edit/2009mb/model.csv'
DEM = '../InData/DEM/SG_DEM_2010.tif'
AM.kriging(dataFile,DEM)