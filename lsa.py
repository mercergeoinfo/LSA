#!/Users/andrew/anaconda/bin/python
#from __future__ import division
"""Degree Day Melt Model"""
__author__ = 'Andrew Mercer'
__contact__ = 'mercergeoinfo@gmail.com'
__maintainer__ = 'Andrew Mercer'
__version__ = '3.0'
__date__ = '15/12/2015'
#
# Created: 16/04/2014
# Edits:
#		16/01/2015
#		20/01/2015
#		21/01/2015
#		26/01/2015
#		03/02/2015
#		05/02/2015
#		16/02/2015
#		17/02/2015
#		17/03/2015
# This file runs a degree day melt model on test data from Storglaciaren
# It is a "working" version, running data from the probing grid
# Reworking of ddf_model.py
#
from datetime import datetime
import sys
import os
import shutil
import fnmatch
import csv
import struct
import copy
import itertools
## SCIENTIFIC
#import math
import numpy as np
#import numpy.ma as ma
#import scipy.stats as stats
from scipy.stats.stats import nanmean
from osgeo import gdal
#from osgeo import ogr
from osgeo import osr
from osgeo.gdalconst import *
## PLOT TOOLS
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
## OWN MODULES
import ParameterAnalysis
from spatialfiles import pt2fmt
from standard import namer, makeDir, makeDirTS, filelist, import2vector
#
# CLASSES
#
class DdfCell:
	# Set a counter for how many objects of this class we have created.
	cellCount = 0
	def __init__(self, east,north,elev,snowmwe,jday, shadevals, paramDict):
		'''Initialise the DdfCell class object,
		setting parameters required by class definition'''
		# Give the values passed during initialisation to internal variables.
		self.easting = east # Easting coordinate of point
		self.northing = north # Northing coordinate of point
		self.elevation = elev # Elevation of point
		self.initSnowMass = snowmwe # Snow mass in m w.e. at beginning of melt model period, Bw
		self.currSnowMass = snowmwe # The snow in m w.e. left at location for current calculation
		self.initjTime = jday # The julian day at beginning of melt model period
		self.jTime = jday # The julian day for current calculation
		self.jTimeSeries = [] # Stores the sequence of days
		self.shadevals = shadevals # The shading factor of point for each day
		self.ddfSnow = paramDict['ddfSnow'] # The degree day factor for snow
		self.ddfIce = paramDict['ddfIce'] # The degree day factor for ice
		self.ddfFirn = paramDict['ddfFirn'] # The degree day factor for firn
		self.ddfSi = paramDict['ddfSi'] # The degree day factor for superimposed ice
		self.ELA = paramDict['ELA'] # The estimated average ELA for determining firn or ice under snow
		self.lapse = paramDict['lapse'] # The adiabatic lapse rate
		self.elevLapse = paramDict['elevLapse'] # An elevation dependent lapse rate adjustment for transition from dry to saturated
		self.sfe = paramDict['sfe'] # An exponential to adjust the effect of shading
		self.refElev = paramDict['refElev'] # Elevation of the reference weather station
		self.meltSum = 0 # The sum of melt
		self.meltSeries = [] # Vectors for storing each melt value
		self.meltSumSeries = []
		self.Bn = 0 # Net Balance
		self.BnSeries = []
		DdfCell.cellCount += 1
	#
	def spatFactor(self):
		'''This offsets the temperature at the weather station by a spatially
		derived factor.'''
		# self.tempDiff is the amount by which the temperature registered at the AWS will be altered
		# Elevation difference to SMHI station.
		self.elevDiff = self.elevation - self.refElev # Fix station elevation
		self.vlr = 1 + self.elevDiff/(self.elevLapse) # Try to capture effect of moisture entering air mass.
		self.tempDiff = self.elevDiff * (self.lapse * self.vlr) # Elevation difference determines temperature difference
		# Shading factor: 1 = no shade, 0 = complete shade throughout the day
		self.shadeFactor = self.shadevals[self.workTime]
		self.tempDiff = self.tempDiff / self.shadeFactor**self.sfe
	#
	def meltModel(self,degC,ddf,jPeriod):
		'''This is the degree day melt model. Temperature passed as mean temperature
		of time interval. Time interval given in julian days.
		It is called by meltInst'''
		# First calculate the temperature at the location in question.
		localTemp = (degC - self.tempDiff)
		# Then run model if the temperature is above 0
		if localTemp > 0:
			melt = jPeriod * ddf * localTemp
		else:
			melt = 0
		return melt
	#
	def meltInst(self,degC,t):
		'''This performs melt calculation for the period in question and then subtracts
		appropriate snow cover and adds melt to the cumulative melt series.
		This method is triggered externally by the script.'''
		self.workTime = t
		# Check if tempDif created, if not then create it
		try:
			self.tempDiff
		except:
			self.spatFactor()
		# Choose DDF based on snow cover
		si_depth = 0.3 # Depth of superimposed ice and fudge value for missing snow
		if self.currSnowMass > si_depth :
			ddf = self.ddfSnow
		elif self.currSnowMass < si_depth and self.currSnowMass > 0.0:
			ddf = self.ddfSi
		elif self.currSnowMass < 0.0 and self.elevation > self.ELA:
			ddf = self.ddfFirn
		else:
			ddf = self.ddfIce
		# Calculate time interval in julian days from last calculation
		jPeriod = self.workTime - self.jTime
		# If time is odd or date prior to Bw survey return 0
		if jPeriod < 0: return 0
		# Send average temperature for elapsed time period and the degree day factor
		# for surface type to the melt model method
		melt = self.meltModel(degC,ddf,jPeriod)
		# Track snow cover for DDF model selection
		self.currSnowMass -= melt
		# Avoid negative snow cover
		if self.currSnowMass < 0:
			self.currSnowMass = 0
		# Sum the total melt
		self.meltSum += melt
		# Add sum to cumulative melt series
		self.meltSumSeries.append(self.meltSum)
		# Add melt to melt series
		self.meltSeries.append(melt)
		# Calculate net balance
		self.Bn = self.initSnowMass - self.meltSum
		# Add net balance to series (cumulative)
		self.BnSeries.append(self.Bn)
		# Add time marker to time series
		self.jTimeSeries.append(self.workTime)
		# Set last calculated time to present
		self.jTime = self.workTime
		return melt
	#
#
# FUNCTIONS
#
def getSettings(DataLoc, lastDate):
	'''Set weather station elevation, julian day of last winter probing and list of julian days for which melt sum to be exported and first day in shade model.
	to use: refElev, jdayBw, jdatelist, startday = getSettings(DataLoc)
	File (settings.txt) to be placed as so ../InData/yyyy/settings.txt, where yyyy is the relevant year.
	Example:
	Elevation=1150
	jdayBw=115
	ExportDate=236,256
	ShadeStart=100
	'''
	# Called by:
	#
	settingsFile = os.path.join(DataLoc, 'settings.txt')
	settings = {}
	InFile = open(settingsFile,'rb')
	# Check contents and set up dictionary
	for row in InFile:
		line = row.strip().split('=')
		print line
		settings[line[0].strip()]=line[1].strip()
	# Set elevation of weather station
	if 'Elevation' not in settings:
		settings['Elevation'] = 1150
		print "Elevation not found. Default value %s used" %(settings['Elevation'])
	refElev = int(settings['Elevation'])
	# Set julian day for probing data
	if 'jdayBw' not in settings:
		settings['jdayBw'] = 115
		print "jdayBw not found. Default value %s used" %(settings['jdayBw'])
	jdayBw = int(settings['jdayBw'])
	# Set julian days for which modelled values are to be exported
	jdatelist = []
	if 'ExportDates' not in settings:
		dates = range(jdayBw+1, lastDate+1)
		for date in dates:
			jdatelist.append(date)
		print "ExportDates not found.\nDefaults to exporting all dates between last probing %s and last temperature reading %s" %(jdayBw, lastDate)
	else:
		expDays = settings['ExportDates'].split(',')
		for j in expDays:
			jdatelist.append(int(j))
	# Give the Julian day at which the shade raster starts. This should not vary from year to year.
	if 'ShadeStart' not in settings:
		startday = 100
	else:
		startday = int(settings['ShadeStart'])
	return refElev, jdayBw, jdatelist, startday
#
def getShadeFile(file):
	'''Get shade map for creating shade value vectors
	To use: raster, transf, bandcount = getShadeFile(file)'''
	# Called by:
	# register all of the GDAL drivers
	gdal.AllRegister()
	# Open file
	raster = gdal.Open(file, GA_ReadOnly)
	if raster is None:
		print 'Could not open ',file,'\n'
		sys.exit(1)
	# Get coordinate system parameters
	projec = raster.GetProjection()
	srs=osr.SpatialReference(wkt=projec) # pyflakes - unused
	transf = raster.GetGeoTransform()
	bandcount = raster.RasterCount
	#
	return raster, transf, bandcount
#
def GetShadeVals(x,y, raster, transf, bandcount, vals, startday):
	'''Create vector of shade factors.
	To use: vals = GetShadeVals(x,y, raster, transf, bandcount, startday)'''
	# Called by:
	# get image size
	success, transfInv = gdal.InvGeoTransform(transf)
	if not success:
		print "Failed InvGeoTransform()"
		sys.exit(1)
	xpix, ypix = gdal.ApplyGeoTransform(transfInv, x, y)
	# Read the file band to a matrix called band_1
	for i in range(1,bandcount+1):
		band = raster.GetRasterBand(i)
		bandtype = gdal.GetDataTypeName(band.DataType)
		if band is None:
			continue
		structval = band.ReadRaster(int(xpix), int(ypix), 1,1, buf_type = band.DataType )
		try:
			fmt = pt2fmt(bandtype)
			print fmt
		except:
			print "fmt error: ", fmt
		try:
			intval = struct.unpack(fmt , structval)
			print intval
		except:
			print "intval error: ", intval
		try:
			vals[i] = intval[0]
		except:
			print "vals error: ", vals[i], intval[0]
	return vals
#
def meltDataOut(pointKeys, points, outDir):
	'''Write results of melt model to csv file
	To use: meltDataOut(pointKeys, points, outDir)'''
	# Called by:
	outFile = os.path.join(outDir,'melt.csv')
	# Create string vector from output dates for Bs and Bn
	outString = []
	for outd in jdatelist:
		outString.append('Mod_Bs_'+str(outd))
		outString.append('Mod_Bn_'+str(outd))
	with open(outFile, 'wb') as csvfile:
		writer = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
		writer.writerow(['Stake','Easting','Northing','Elevation','Bw']+outString)
		for j in pointKeys:
			# list to write as row. Start with point name and elevation
			output =[j,round(points[j]['easting'],2),round(points[j]['northing'],2),int(points[j]['elevation']),points[j]['Bw']]
			# then add modelled data
			for d in jdatelist:
				loc = points[j]['MeltModel'].jTimeSeries.index(d)
				output.append(round(points[j]['MeltModel'].meltSumSeries[loc],3))
				output.append(round(points[j]['MeltModel'].BnSeries[loc],3))
			writer.writerow(output)
	writer = None
#
def plotDifElev(outputname,outDir,x,y,colour):
	'''Plot modelled data.
	To use: plotDifElev(outputname, outDir, x, y, colour)'''
	# Called by:
	#
	matplotlib.rcParams['axes.grid'] = True
	matplotlib.rcParams['legend.fancybox'] = True
	# matplotlib.rcParams['figure.figsize'] = 18, 9 # Mine
	# matplotlib.rcParams['figure.figsize'] = 16.54, 11.69 # A3
	matplotlib.rcParams['figure.figsize'] = 11.69, 8.27 # A4
	matplotlib.rcParams['savefig.dpi'] = 300
	plotName = outputname + '.pdf'
	pp1 = PdfPages(os.path.join(outDir,plotName))
	fig1 = plt.figure(1)
	ax1 = fig1.add_subplot(111)
	# lnclr = ['k','r','g','b','y','c','m'] # pyflakes - unused
	# lnsty = ['-','--','-.',':'] # pyflakes - unused
	# mrsty = ['o','s','v','*','x','+','1','2','3','4'] # pyflakes - unused
	xmax = max(x)
	xmin = min(x)
	# ymax = 3 # pyflakes - unused
	# ymin = -3 # pyflakes - unused
	labelString = outputname.split("_")[0] + " " + outputname.split("_")[-1]
	ax1.plot(x,y, color=colour, marker='o', linestyle='None', label=labelString)
	matplotlib.pyplot.axes().set_position([0.04, 0.065, 0.8, 0.9])
	ax1.legend(bbox_to_anchor=(0.0, 0), loc=3, borderaxespad=0.1, ncol=3, title = "Component, Julian Day")
	ax1.plot([xmin,xmax],[0,0],'-k')
	#lastind = str(jdatelist[-1])
	#ax1.plot([0,5],np.multiply(dbvecs[lastind]['slope'],[0,5]) + dbvecs[lastind]['intercept'],'-b')
	#plt.axis([xmin, xmax, ymin, ymax*1.2])
	plt.axis([xmin, xmax, -2, 2])
	# for pnt in range(len(dbvecs[jd]['Stake'])):
	# ax1.annotate(dbvecs[jd]['Stake'][pnt],xy=(dbvecs[jd]['Elevation'][pnt],1.8), rotation=90)
	plt.xlabel('Elevation (m.a.s.l.')
	plt.ylabel('Measured - Modelled Melt (m w.e.)')
	plt.title('Measured - Modelled against Elevation')
	# plt.show()
	pp1.savefig(bbox_inches='tight')
	pp1.close()
	plt.close()
	return 0
#
#
################################### MAIN ###################################
#
# To run this you will need the following:
# A shading file where each pixel represents shading on the glacier at each julian day to be modelled
# A temperature data file for each julian day to be modelled
# A stake data file containing winter balance values and coordinates for each stake
# A settings file
# The directory structure I have used here is as follows ('yyyy' should be replaced by the year):
# /InData/yyyy
# /Indata/yyyy/weatheryyyy.csv
# /Indata/yyyy/StakeReadings.csv
# /Indata/settings.txt
# /Output/
# /Output/Shades/SG_shade.tif 		This file is here as it is the output of another script
# /Output/yyyy/		These are created by this script as needed
# /Scripts/
#
# Format examples:
	# settings.txt (note that if 'ExportDate' omitted then all dates between bW and final temperature reading exported):
		# Elevation=1150
		# jdayBw=115
		# ExportDate=236,256
		# ShadeStart=100
	# weatheryyyy.csv:
		# Date,Temp
		# 2010-01-25,-8.3
	# StakeDatayyyy.csv:
		# Stake,Easting,Northing,Elevation,Bw,Bs,Bn,Surface
		# 04C,651103.586397,7536381.86553,1219,0.334,2.53,-2.196,ice
#
time_zero = datetime.now()
time_one = time_zero
#
# Get shading data
# Location of file containing a multiband raster, each band represents the shading on one day. 1 = no shade, 0 = really quite dark
shadefile = '../InData/Shades/SG_shade.tif'
# Read the shade factor raster in to memory
raster, transf, bandcount = getShadeFile(shadefile)
#
# Set test to 0 to run a single set of parameters (2009 to 2013 best parameters) and don't compare results with measured data
# Set test to 1 to run all parameters and compare results with field data
# Set test to 2 to run 2005 to 2008 best parameters and compare results with field data
# Set test to 3 (or higher) to run 2005 to 2013 best parameters and compare results with field data
test = 1
#
# [2009,2010,2011,2012,2013]
# [2005,2006,2007,2008]
# [2005,2006,2007,2008,2009,2010,2011,2012,2013]
# Choose which data set is to be run.
for year in [2013]:
	strYear = str(year)
	dataLoc = '../InData/' + strYear
	# Temperature data. Following two lines are example of input format:
	# Date,Temp
	# 2010-01-25,-8.3
	weather = 'weather' + strYear + '.csv'
	TfileName = os.path.join(dataLoc, weather)
	# Read Temperature data from csv file and convert dates to julian days. Date format '%Y-%m-%d' is SMHI's
	TinFile = open(TfileName,'rb')
	dates = []
	times = []
	temps = []
	for line in csv.DictReader(TinFile, delimiter=','):
		dates.append(line['Date'].strip())
		date = datetime.strptime(line['Date'].strip(),'%Y-%m-%d')
		jdate = datetime.strftime(date,'%j')
		times.append(int(jdate))
		temps.append(float(line['Temp'].strip()))
	TinFile.close()
	# Stake data. Following two lines are example of input format:
	# Stake,Easting,Northing,Elevation,Bw,Bs,Bn,Surface
	# 04C,651103.586397,7536381.86553,1219,0.334,2.53,-2.196,ice
	stakeFileName = 'StakeData' + strYear + '.csv'
	SfileName =os.path.join(dataLoc, stakeFileName)
	# File not read until loop through parameters. This could be improved by creating point object first and then inserting parameters on each run.
	#
	# Get settings for model: AWS elevation, date of snow probing, dates for model export, first date in shading file (could start this at 1 by default but
	# shading file created for limited range of dates to reduce file size)
	refElev, jdayBw, jdatelist, startday = getSettings(dataLoc, times[-1])
	print "For year %s following settings used: " %(strYear)
	print "refElev set to %s" %(refElev)
	print "jdayBw set to %s" %(jdayBw)
	#
	# Set parameters for the melt model
	if test == 1: # Block 2, 2009 to 2013 best parameters
		ddfSnow = 0.0038
		ddfSi = 0.0051
		ddfFirn = 0.0050
		ddfIce = 0.0054
		lapse = 0.0055
		elevLapse = (2100 - 1150) # Elevation dependant lapse rate
		sfe = 1.5 # Shading factor exponent (adjusts the shading value at each point)
		ELA = 1500 # Equilibrium line, for firn or ice under snow
	elif test == 2: # Block 1, 2005 to 2008 best parameters
		ddfSnow=0.0040
		ddfSi=0.0048
		ddfFirn=0.0050
		ddfIce=0.0050
		lapse=0.0064
		elevLapse = (2100 - 1150) # Elevation dependant lapse rate
		sfe = 1.5 # Shading factor exponent (adjusts the shading value at each point)
		ELA = 1500 # Equilibrium line, for firn or ice under snow
	else: #All, 2005 to 2013 best parameters
		ddfSnow = 0.0039
		ddfSi = 0.0049
		ddfFirn = 0.0050
		ddfIce = 0.0052
		lapse = 0.0059
		elevLapse = (2100 - 1150) # Elevation dependant lapse rate
		sfe = 1.5 # Shading factor exponent (adjusts the shading value at each point)
		ELA = 1500 # Equilibrium line, for firn or ice under snow
	paramDict = {}
	paramDict['refElev'] = refElev
	paramDict['elevLapse'] = elevLapse
	paramDict['sfe'] = sfe
	paramDict['ELA'] = ELA
	paramDict['ddfSnow'] = ddfSnow
	paramDict['ddfSi'] = ddfSi
	paramDict['ddfFirn'] = ddfFirn
	paramDict['ddfIce'] = ddfIce
	paramDict['lapse'] = lapse
	#
	# Directory for output
	outputDir = os.path.join('../Output/', strYear)
	if not os.path.exists(outputDir):
		os.makedirs(outputDir)
	outputname = strYear +'_DDM'
	truthDir = os.path.join(dataLoc,"truthing")
	truthfiles = filelist(truthDir,'csv')
	print "Truthing files: "
	for tf in truthfiles: print tf
	#
	#
	trudb = {} # FOR VERIFICATION OF RESULTS ONLY
	for file in truthfiles: # FOR VERIFICATION OF RESULTS ONLY
		trudb[file.split('.')[0]] = import2vector(os.path.join(truthDir,file)) # FOR VERIFICATION OF RESULTS ONLY
	counter = 0
	# Read stake data
	SinFile = open(SfileName,'rb')
	stakeData = {}
	# Read point data file with position and winter balance as of last winter probing (jdayBw) and send to model
	for line in csv.DictReader(SinFile, delimiter=','):
		stakeName = line['Stake'].strip()
		stakeData[stakeName] = {}
		# Coordinates
		stakeData[stakeName]['easting'] = float(line['Easting'].strip())
		stakeData[stakeName]['northing'] = float(line['Northing'].strip())
		stakeData[stakeName]['elevation'] = float(line['Elevation'].strip())
		# Get shading factor for location
		# Create vector for shade values
		vals = []
		for d in range(366):
			vals.append(1)
		try:
			stakeData[stakeName]['Shadevals'] = GetShadeVals(stakeData[stakeName]['easting'], stakeData[stakeName]['northing'], raster, transf, bandcount, vals, startday)
		except:
			stakeData[stakeName]['Shadevals'] = vals
			print "No shade value obtained for ", stakeName
		# Get the measured winter balance
		try:
			stakeData[stakeName]['Bw'] = float(line['Bw'].strip())
		except:
			print "No winter balance data found (Bw column)"
			break
		# Get the measured summer balance
		try:
			stakeData[stakeName]['Bs'] = float(line['Bs'].strip())
		except:
			pass
		# Get the measured net balance
		try:
			stakeData[stakeName]['Bn'] = float(line['Bn'].strip())
		except:
			pass
	#
	#
	points = {}
	data = copy.deepcopy(stakeData)
	stakeNames = stakeData.keys()
	# Read stake data
	for stake in stakeNames:
		# Send input data to Degree Day Model object
		data[stake]['MeltModel'] = DdfCell(data[stake]['easting'], data[stake]['northing'], data[stake]['elevation'], data[stake]['Bw'], jdayBw, data[stake]['Shadevals'], paramDict)
		points[stake] = copy.deepcopy(data[stake])
	pointKeys = points.keys()
	pointKeys.sort()
	#
	# For each julian day in the "times" vector call the meltInst method for each point object, passing the temperature and the day number.
	# This is what runs the model at each time step in the temperature time series file
	for i in range(len(temps)):
		for j in pointKeys:
			points[j]['MeltModel'].meltInst(temps[i],times[i])
	#
	# Create database of vectors for assessment of model output
	dbcomp = {}
	# Get number of modelled objects
	pntsz = len(pointKeys)
	# Create empty arrays for data
	Stk_mod = np.zeros((pntsz,))
	Stk_mod.fill(np.nan)
	dbcomp['Stake'] = [None]*pntsz
	dbcomp['Easting'] = np.copy(Stk_mod)
	dbcomp['Northing'] = np.copy(Stk_mod)
	dbcomp['Elevation'] = np.copy(Stk_mod)
	dbcomp['Bw'] = np.copy(Stk_mod)
	for d in jdatelist:
		bsn = "Mod_Bs_" + str(d)
		bnn = "Mod_Bn_" + str(d)
		dbcomp[bsn] = np.copy(Stk_mod)
		dbcomp[bnn] = np.copy(Stk_mod)
	for pki in range(pntsz):
		Stk = pointKeys[pki]
		dbcomp['Stake'][pki] = pointKeys[pki]
		dbcomp['Easting'][pki] = points[Stk]['easting']
		dbcomp['Northing'][pki] = points[Stk]['northing']
		dbcomp['Elevation'][pki] = points[Stk]['elevation']
		dbcomp['Bw'][pki] = points[Stk]['Bw']
		for d in jdatelist:
			bsn = "Mod_Bs_" + str(d)
			bnn = "Mod_Bn_" + str(d)
			loc = points[j]['MeltModel'].jTimeSeries.index(d)
			dbcomp[bsn][pki] = round(points[Stk]['MeltModel'].meltSumSeries[loc],3)
			dbcomp[bnn][pki] = round(points[Stk]['MeltModel'].BnSeries[loc],3)
	# Create lists for R-squared values
	bsR2list = []
	bnR2list = []
	# Loop through model results (in vectors)
	for d in jdatelist:
		obsBs = []
		obsBn = []
		for stk in dbcomp['Stake']:
			try:
				ind = np.where(trudb[str(d)]['Stake'] == stk)
				bs = trudb[str(d)]['Bs'][ind][0]
				bn = trudb[str(d)]['Bn'][ind][0]
			except:
				bs = np.nan
				bn = np.nan
			obsBs.append(bs)
			obsBn.append(bn)
		bsn = "Mod_Bs_" + str(d)
		bnn = "Mod_Bn_" + str(d)
		# Calculate 'R-squared'
		# Bs
		modBs = dbcomp[bsn]
		obsBsmean = nanmean(obsBs)
		obsBsMinModBs = obsBs - modBs
		obsBsMinMean = obsBs - obsBsmean
		BsR2 = 1 - ((np.nansum(obsBsMinModBs**2)) / (np.nansum(obsBsMinMean**2)))
		paramDict[(bsn+'_R2')] = BsR2
		# Bn
		modBn = dbcomp[bnn]
		obsBnmean = nanmean(obsBn)
		obsBnMinModBn = obsBn - modBn
		obsBnMinMean = obsBn - obsBnmean
		BnR2 = 1 - ((np.nansum(obsBnMinModBn**2)) / (np.nansum(obsBnMinMean**2)))
		paramDict[(bnn+'_R2')] = BnR2
		# Output model data to file
		flnm = str(counter)
		outDir = makeDir(outputDir, flnm)
		meltDataOut(pointKeys, points, outDir)
		# Write parameters used to text file
		paramFile = os.path.join(outDir,'Parameters.txt')
		with open (paramFile, 'w') as fp:
			# for p in paramDict.items():
			for p in sorted(paramDict.keys()):
				fp.write("%s:%2.6f\n" % (p, paramDict[p]))
		# Plot model results
		x = dbcomp['Elevation']
		for d in jdatelist:
			obsBs = []
			obsBn = []
			for stk in dbcomp['Stake']:
				try:
					ind = np.where(trudb[str(d)]['Stake'] == stk)
					bs = trudb[str(d)]['Bs'][ind][0]
					bn = trudb[str(d)]['Bn'][ind][0]
				except:
					bs = np.nan
					bn = np.nan
				obsBs.append(bs)
				obsBn.append(bn)
			bsn = "Mod_Bs_" + str(d)
			bnn = "Mod_Bn_" + str(d)
			modBs = dbcomp[bsn]
			modBn = dbcomp[bnn]
			bsDiff = obsBs - modBs
			bnDiff = obsBn - modBn
			pltnmBs = outputname + str(d) + '_Bs_diff_'
			pltnmBn = outputname + str(d) + '_Bn_diff_'
			plotDifElev(pltnmBs,outDir,x,bsDiff,'r')
			plotDifElev(pltnmBn,outDir,x,bnDiff,'k')
	counter = counter +1

