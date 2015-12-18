#!/Users/andrew/anaconda/bin/python
# -*- coding: utf-8 -*-
#from __future__ import division
"""Shading model"""
__author__ = 'Andrew Mercer'
__contact__ = 'mercergeoinfo@gmail.com'
__maintainer__ = 'Andrew Mercer'
__version__ = '1.0'
__date__ = '25/04/2014'
#
# Created: 25/01/2014
# Version: 1.0
# Andrew Mercer
# Calculate Sun's position from Latitude, Longitude and Time/Date then calculate shade on Storglaciären
#
import ephem
from datetime import datetime
import sys
import os
import math
import numpy as np
from scipy.stats.stats import nanmean
import scipy.stats as stats
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
#
# FUNCTIONS
#
def rasterImport(file):
	'''Import a raster file and return grid plus meta data.
	To use: data, meta, metadata = rasterImport(file)'''
	## Called by:
	print '\nImporting ',file
	#print time_one.strftime('at day:%j %H:%M:%S')
	# register all of the GDAL drivers
	gdal.AllRegister()
	# Open file
	try:
		raster = gdal.Open(file)
	except:
		print "gdal import error"
	if raster is None:
		print 'Could not open ',file,'\n'
		sys.exit(1)
	# Get coordinate system parameters
	projec = raster.GetProjection()
	srs=osr.SpatialReference(wkt=projec)
	transf = raster.GetGeoTransform()
	ul_x = transf[0]
	ul_y = transf[3]
	xres = transf[1]
	yres = transf[5]
	# get image size
	rows = raster.RasterYSize
	cols = raster.RasterXSize
	dims = {'xres':xres,'yres':yres,'rows':rows,'cols':cols}
	# Calculate corners
	ll_x = ul_x
	ll_y = ul_y + (rows * yres)
	ur_x = ul_x + (cols * xres)
	ur_y = ul_y
	lr_x = ur_x
	lr_y = ll_y
	corners = {'ll':(ll_x,ll_y),'ul':(ul_x,ul_y),'ur':(ur_x,ur_y),'lr':(lr_x,lr_y)}
	#
	driveshrt = raster.GetDriver().ShortName
	driver = gdal.GetDriverByName(driveshrt)
	metadata = driver.GetMetadata()
	metakeys = metadata.keys()
	# Read the file band to a matrix called band_1
	band = raster.GetRasterBand(1)
	# Access data in rastern band as array
	data = band.ReadAsArray(0,0,cols,rows)
	#
	# gdal interpolation creates "upside down files, hence this
	Yres=yres
	if Yres/np.fabs(Yres)==-1:
		Yres=-1*Yres
		data = np.flipud(data)
	#
	# get nodata value
	nandat = band.GetNoDataValue()
	# get minimum and maximum value
	mindat = band.GetMinimum()
	maxdat = band.GetMaximum()
	if mindat is None or maxdat is None:
			(mindat,maxdat) = band.ComputeRasterMinMax(1)
	dats = [nandat,mindat,maxdat]
	sumvals = data[np.where(np.logical_not(data == nandat))]
	sumval = sumvals.sum()
	numvals = len(sumvals)
	avval = sumval/numvals
	volvals = sumval*(math.fabs((transf[1]*transf[5])))
	meta = {'transform':transf,'projection':projec,'corners':corners,'dimension':dims,'dats':dats,'sumval':sumval,'numvals':numvals,'avval':avval,'volvals':volvals}
#	if srs.IsProjected:
# 		print 'projcs 1: ',srs.GetAttrValue('projcs')
# 		print 'projcs 2: ',srs.GetAuthorityCode('projcs')
# 		print 'projcs 3: ',srs.GetAuthorityName('projcs')
# 	print 'geogcs 1:: ',srs.GetAttrValue('geogcs')
# 	print 'Sum of pixel values = ',sumval,' in range ',mindat,' to ',maxdat
# 	print 'From ',numvals,' of ',rows,' X ',cols,' = ',rows*cols,' pixels, or ',((1.0*numvals)/(rows*cols))*100,'%'
# 	print 'Average then is ',avval
# 	print 'Pixels size, x:',transf[1],' y:',transf[5]
# 	print '"Volume" is then the sum * 1 pixel area: ',volvals
	return data, meta, metadata
#
#
# MAIN
#
gdal.AllRegister()
# Give year (not needed after this)
year = 2012
# Set up vectors for sun position data
Clock = []
JDay = []
Azim = []
Alti = []
days = range(1,366)
#days = range(100,260)
#days = range(150,160)
hours = range(6,19)
#hours = range(6,9)
# Loop through days of year and daylight hours (06.00 to 18.00)
for jday in days:
	for hour in hours:
		dtin = str(year) + str(jday) + str(hour)
		dateTime = datetime.strptime(dtin, '%Y%j%H').strftime('%Y/%m/%d %H:%M:%S')
		Latitude = 67.8974541
		Longitude = 18.5728927
		gatech = ephem.Observer()
		gatech.lon, gatech.lat = Longitude, Latitude
		gatech.date = dateTime	 # 12:22:56 EDT
		sun = ephem.Sun()
		sun.compute(gatech)
		#print dateTime
		#print("%s %s" % (sun.alt, sun.az))
		Clock.append(dateTime)
		JDay.append(jday)
		Azim.append(math.degrees(sun.az))
		Alti.append(math.degrees(sun.alt))
#
# Get DEM of entire region
RegDEM = '../InData/DEM/1999_KebDEM.tiff'

# Get DEM of Storglaciären
SG = '../InData/DEM/SG_DEM_2010.tif'
SGdata, SGmeta, SGmetadata = rasterImport(SG)
ul = SGmeta['corners']['ul']
lr = SGmeta['corners']['lr']
SGmask = '../InData/DEM/SG_DEM_2010_mask.tif'
#
# Output folder
outDir = '../InData/Shades/'
# Loop through days
for day in days:
	print "Day: ", day
	DailyName = str(day)+'.tif'
	DailyFile = os.path.join(outDir,DailyName)
	bashCommand0 = 'gdal_calc.py -A %s -B %s --outfile=%s --calc="A-B" --NoDataValue=-9999' % (SGmask, SGmask, DailyFile)
	bashRun = os.system(bashCommand0)
	if bashRun != 0:
		print "bash 0: %s failed\n" % bashCommand0
		sys.exit()
	for i in np.where(np.array(JDay)==day)[0]:
		# Calculate Shade at time stamp
		az = Azim[i]
		alt = Alti[i]
		dateTime = datetime.strptime(Clock[i], '%Y/%m/%d %H:%M:%S').strftime('%j%H')
		FileName = dateTime+'.tif'
		OutFile = os.path.join(outDir,FileName)
		bashCommand1 = 'gdaldem hillshade -az %s -alt %s %s %s' % (str(az), str(alt), RegDEM, OutFile)
		bashRun = os.system(bashCommand1)
		if bashRun != 0:
			print "bash 1: %s failed\n" % bashCommand1
			sys.exit()
		#
		# Clip to Storglaciären
		FileName2 = dateTime+'_c.tif'
		ClippedFile = os.path.join(outDir,FileName2)
		bashCommand2 = 'gdal_translate -projwin %s %s %s %s %s %s' % (ul[0], ul[1], lr[0], lr[1], OutFile, ClippedFile)
		bashRun = os.system(bashCommand2)
		if bashRun != 0:
			print "bash 2: %s failed\n" % bashCommand2
			sys.exit()
		#
		# Delete regional shade map
		bashCommand3 = 'rm %s' % (OutFile)
		bashRun = os.system(bashCommand3)
		if bashRun != 0:
			print "bash 3: %s failed\n" % bashCommand3
			sys.exit()
		#
		# Add hourly shade map to daily sum
		bashCommand4 = 'gdal_calc.py -A %s -B %s --outfile=%s --calc="A+B" --NoDataValue=-9999' % (DailyFile, ClippedFile, DailyFile)
		bashRun = os.system(bashCommand4)
		if bashRun != 0:
			print "bash 4: %s failed\n" % bashCommand4
			sys.exit()		#
		# Erase hourly map
		bashCommand5 = 'rm %s' % (ClippedFile)
		bashRun = os.system(bashCommand5)
		if bashRun != 0:
			print "bash 5: %s failed\n" % bashCommand5
			sys.exit()
		#
		# Convert to %
		# Open file
		raster = gdal.Open( DailyFile)
		band = raster.GetRasterBand(1)
		mindat = band.GetMinimum()
		maxdat = band.GetMaximum()
		if mindat is None or maxdat is None:
			(mindat,maxdat) = band.ComputeRasterMinMax(1)
		DailyPName = '{:03d}'.format(day)+ '_perc.tif'
		DailyPFile = os.path.join(outDir,DailyPName)
		bashCommand6 = 'gdal_calc.py -A %s --outfile=%s --calc="1-A/%s" --NoDataValue=-9999' % (DailyFile, DailyPFile, str(maxdat))
		bashRun = os.system(bashCommand6)
		if bashRun != 0:
			print "bash 6: %s failed\n" % bashCommand6
			sys.exit()
		#
	# Erase Daily Sum
	bashCommand7 = 'rm %s' % (DailyFile)
	bashRun = os.system(bashCommand7)
	if bashRun != 0:
		print "bash 7: %s failed\n" % bashCommand7
		sys.exit()
# Merge files
bashCommand8 = 'gdal_merge.py -separate -a_nodata -9999 -o %sSG_shade.tif %s*.tif' % (outDir, outDir)
bashRun = os.system(bashCommand8)
if bashRun != 0:
	print "bash 8: %s failed\n" % bashCommand8
