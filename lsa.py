#!/Users/andrew/anaconda/bin/python
#from __future__ import division
"""Late Season Accumulation Model"""
__author__ = 'Andrew Mercer'
__contact__ = 'mercergeoinfo@gmail.com'
__maintainer__ = 'Andrew Mercer'
__version__ = '1.0'
__date__ = '16/12/2015'
#
# Created: 16/12/2015
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
		startday = 1
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
		except:
			print "GetShadeVals fmt error: ", fmt
		try:
			intval = struct.unpack(fmt , structval)
		except:
			print "GetShadeVals intval error: ", intval
		try:
			vals[i] = intval[0]
		except:
			print "GetShadeVals vals error: ", vals[i], intval[0]
	return vals
#
def meltDataWrite(points, outDir):
	'''Write results of melt model to csv file
	To use: meltDatatWrite(points, outDir)'''
	outFile = os.path.join(outDir,'melt.csv')
	stakes = points.keys()
	stakes.remove('DataSets')
	stakes.sort()
	headers = points[stakes[0]]['Headers']
	headers.remove('MeltModel')
	headers.remove('Shadevals')
	headers.insert(0,'Stake')
	# Write to file
	with open(outFile, 'wb') as csvfile:
		writer = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
		# Write header to top of file
		writer.writerow(headers)
		# Write data row by row (stake by stake)
		headers.remove('Stake')
		for stake in stakes:
			outputRow = [stake]
			for column in headers:
				outputRow.append(points[stake][column])
			writer.writerow(outputRow)
	writer = None
#
def reportWrite(report, outDir):
	reportFile = os.path.join(outDir,'Report.txt')
	with open (reportFile, 'w') as fp:
		for p in sorted(report.keys()):
			fp.write("%s:%2.4f\n" % (p, report[p]))
	return 0
#
def plotDifElev(outputname,outDir, title, x, y, colour):
	'''Plot modelled data.
	To use: plotDifElev(outputname,outDir, title, x, y, colour)'''
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
	# labelString = outputname.split("_")[0] + " " + outputname.split("_")[-1]
	labelString = title
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
# Set test to 1 to run 2005 to 2008 best parameters
# Set test to 2 to run2009 to 2008 best parameters
# Set test to 3 (or higher) to run 2005 to 2013 best parameters
test = 3
year = 2013
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
# Get settings for model: AWS elevation, date of snow probing, dates for model export, first date in shading file (could start this at 1 by default but
# shading file may be created for limited range of dates to reduce file size)
refElev, jdayBw, jdatelist, startday = getSettings(dataLoc, times[-1])
print "For year %s following settings used: " %(strYear)
print "refElev set to %s" %(refElev)
print "jdayBw set to %s" %(jdayBw)
#
# Set parameters for the melt model
if test == 1: # Block 1, 2005 to 2008 best parameters
	ddfSnow=[0.0040]
	ddfSi=[0.0048]
	ddfFirn=[0.0050]
	ddfIce=[0.0050]
	lapse=[0.0064]
	rangeZ = (2100 - 1150)
	elevLapse =[rangeZ]# Elevation dependant lapse rate
	sfe = [1.5] # Shading factor exponent (adjusts the shading value at each point)
	ELA = [1500] # Equilibrium line, for firn or ice under snow
elif test == 2: # Block 2, 2009 to 2013 best parameters
	ddfSnow = [0.0038]
	ddfSi = [0.0051]
	ddfFirn = [0.0050]
	ddfIce = [0.0054]
	lapse = [0.0055]
	rangeZ = (2100 - 1150)
	elevLapse = [rangeZ] # Elevation dependant lapse rate
	sfe = [1.5] # Shading factor exponent (adjusts the shading value at each point)
	ELA = [1500] # Equilibrium line, for firn or ice under snow
elif test == 3:
	ddfSnow = [0.0030,0.0032,0.0034,0.0036,0.0038,0.0040,0.0042,0.0044,0.0046,0.0048]
	ddfSi = [0.0045,0.0047,0.0049,0.0051]
	ddfFirn = [0.0043,0.0045,0.0048,0.0050,0.0053,0.0055,0.0058]
	ddfIce = [0.0043,0.0045,0.0047,0.0049,0.0051,0.0053,0.0055,0.0057]
	lapse = [0.0042,0.0044,0.0046,0.0048,0.0050,0.0055,0.0057,0.0059,0.0061,0.0063,0.0065,0.0068,0.0070,0.0072,0.0074,0.0076]
	rangeZ = (2100 - 1150)
	elevLapse = [rangeZ] # Elevation dependant lapse rate
	sfe = [1.5] # Shading factor exponent (adjusts the shading value at each point)
	ELA = [1500] # Equilibrium line, for firn or ice under snow
else: #All, 2005 to 2013 best parameters
	ddfSnow = [0.0039]
	ddfSi = [0.0049]
	ddfFirn = [0.0050]
	ddfIce = [0.0052]
	lapse = [0.0059]
	rangeZ = (2100 - 1150)
	elevLapse = [rangeZ] # Elevation dependant lapse rate
	sfe = [1.5] # Shading factor exponent (adjusts the shading value at each point)
	ELA = [1500] # Equilibrium line, for firn or ice under snow
counter = 0
bestBsR2 = -9999
writeTest = 0
for it1, it2, it3, it4, it5, it6, it7, it8 in itertools.product(ddfSnow, ddfSi, ddfFirn, ddfIce, lapse, elevLapse, sfe, ELA):
	paramDict = {}
	paramDict['ddfSnow'] = it1
	paramDict['ddfSi'] =it2
	paramDict['ddfFirn'] = it3
	paramDict['ddfIce'] = it4
	paramDict['lapse'] = it5
	paramDict['elevLapse'] = it6
	paramDict['sfe'] = it7
	paramDict['ELA'] = it8
	paramDict['refElev'] = refElev
	#
	# Directory for output
	outputDir = os.path.join('../Output/', strYear)
	if not os.path.exists(outputDir):
		os.makedirs(outputDir)
	outputname = strYear +'_DDM_'
	#
	# Truthing of data against field survey data. Each survey stored in sperate csv file.
	# The trudb is used to store both data and assessment of results
	truthDir = os.path.join(dataLoc,"truthing")
	try:
		truthfiles = filelist(truthDir,'csv')
		trudb = {}
		print "Truthing files: "
		for file in truthfiles:
			if int(file.split('.')[0]) not in jdatelist:
				print "%s does not match date given in settings file: %s" % (file, jdatelist)
			else:
				print file
			trudb[file.split('.')[0]] = import2vector(os.path.join(truthDir,file))
		trudbKeys = trudb.keys()
		trudbKeys.sort()
	except:
		print "No truthing data found."
	#
	# Read stake data
	SinFile = open(SfileName,'rb')
	stakeData = {}
	# Read point data file with position and winter balance as of last winter probing (jdayBw) and send to model
	for line in csv.DictReader(SinFile, delimiter=','):
		stakeName = line['Stake'].strip()
		stakeData[stakeName] = {}
		# Coordinates
		stakeData[stakeName]['Easting'] = float(line['Easting'].strip())
		stakeData[stakeName]['Northing'] = float(line['Northing'].strip())
		stakeData[stakeName]['Elevation'] = float(line['Elevation'].strip())
		# Get shading factor for location
		# Create vector for shade values
		vals = []
		for d in range(366):
			vals.append(1)
		try:
			stakeData[stakeName]['Shadevals'] = GetShadeVals(stakeData[stakeName]['Easting'], stakeData[stakeName]['Northing'], raster, transf, bandcount, vals, startday)
		except:
			stakeData[stakeName]['Shadevals'] = vals
			print "No shade value obtained for ", stakeName
		# Get the measured winter balance
		try:
			stakeData[stakeName]['Org_Bw'] = float(line['Bw'].strip())
		except:
			print "No winter balance data found (Bw column)"
			break
		# Get the measured summer balance
		try:
			stakeData[stakeName]['Org_Bs'] = float(line['Bs'].strip())
		except:
			pass
		# Get the measured net balance
		try:
			stakeData[stakeName]['Org_Bn'] = float(line['Bn'].strip())
		except:
			pass
	#
	# 'points' stores each survey point and the model results, all inherited from 'data'
	# 'data' is a copy of the original 'stakeData'
	points = {}
	data = copy.deepcopy(stakeData)
	stakeNames = stakeData.keys()
	stakeNames.sort()
	data['DataSets'] = {}
	for stake in stakeNames:
		# For ordered headers/keys
		data[stake]['Headers'] = ['MeltModel', 'Shadevals', 'Easting', 'Northing', 'Elevation', 'Org_Bw']
		if 'Org_Bn' in data[stake].keys():
			data[stake]['Headers'].append('Org_Bn')
		# Send input data to Degree Day Model object
		data[stake]['MeltModel'] = DdfCell(data[stake]['Easting'], data[stake]['Northing'], data[stake]['Elevation'], data[stake]['Org_Bw'], jdayBw, data[stake]['Shadevals'], paramDict)
		# For each julian day in the "times" vector call the meltInst method for each point object, passing the temperature and the day number.
		# This is what runs the model at each time step in the temperature time series file
		for i in range(len(temps)):
			data[stake]['MeltModel'].meltInst(temps[i],times[i])
		for day in jdatelist:
			# Fetch modelled melt and net balance for each julian day specific in settings and create new entry for each
			loc = data[stake]['MeltModel'].jTimeSeries.index(day)
			data[stake]['Mod_Bs_' + str(day)] =  round(data[stake]['MeltModel'].meltSumSeries[loc],3)
			data[stake]['Mod_Bn_' + str(day)] =  round(data[stake]['MeltModel'].BnSeries[loc],3)
			data[stake]['Headers'].append('Mod_Bs_' + str(day))
			data[stake]['Headers'].append('Mod_Bn_' + str(day))
			# Fetch any truthing data available
			if 'trudbKeys' in locals():
				if str(day) in trudbKeys:
					loc = np.where(trudb[str(day)]['Stake']==stake)[0][0]
					data[stake]['Org_Bs_' + str(day)] = round(trudb[str(day)]['Bs'][loc],3)
					data[stake]['Org_Bn_' + str(day)] = round(trudb[str(day)]['Bn'][loc],3)
					data[stake]['Mod_Bw_' + str(day)] = round((data[stake]['Org_Bn_' + str(day)] +data[stake]['Mod_Bs_' + str(day)]), 3)
					data[stake]['Headers'].insert(-2, 'Org_Bs_' + str(day))
					data[stake]['Headers'].insert(-2, 'Org_Bn_' + str(day))
					data[stake]['Headers'].insert(-2, 'Mod_Bw_' + str(day))
					# Add values to lists for calculating R2 later
					if 'Mod_Bs_' + str(day) not in data['DataSets'].keys():
						data['DataSets']['Mod_Bs_' + str(day)] = []
					data['DataSets']['Mod_Bs_' + str(day)].append(data[stake]['Mod_Bs_' + str(day)])
					if 'Org_Bs_' + str(day) not in data['DataSets'].keys():
						data['DataSets']['Org_Bs_' + str(day)] = []
					data['DataSets']['Org_Bs_' + str(day)].append(data[stake]['Org_Bs_' + str(day)])
		dataKeys = data.keys()
		dataKeys.sort()
	#
	if len(data['DataSets']) > 0:
		report = copy.deepcopy(paramDict)
		setKeys = data['DataSets'].keys()
		# Order all Mod first, then all Org
		setKeys.sort()
		start = 0
		end = len(setKeys)
		middle = end/2
		i = start
		while i < end/2:
			# Calculate R2
			modBs = np.array(data['DataSets'][setKeys[i]])
			obsBs = np.array(data['DataSets'][setKeys[middle]])
			modBsmean = nanmean(modBs)
			obsBsmean = nanmean(obsBs)
			obsBsMinModBs = obsBs - modBs
			obsBsMinMean = obsBs - obsBsmean
			BsR2 = 1 - ((np.nansum(obsBsMinModBs**2)) / (np.nansum(obsBsMinMean**2)))
			report[(setKeys[i]+'_R2')] = BsR2
			if i == 1:
				if BsR2 > bestBsR2:
					bestBsR2 = BsR2
					writeTest = 1
			i = i+1
			middle = middle+1

	if writeTest == 1:
		# Output model data to file
		flnm = str(counter)
		outDir = makeDir(outputDir, flnm)
		meltDataWrite(data, outDir)
		# 	Write report to text file
		reportWrite(report, outDir)
		# 	Plot model results
		x = []
		for stake in stakeNames:
			x.append(data[stake]['Elevation'])
		if len(data['DataSets']) > 0:
			setKeys = data['DataSets'].keys()
			# Order all Mod first, then all Org
			setKeys.sort()
			start = 0
			end = len(setKeys)
			middle = end/2
			i = start
			while i < end/2:
				modBs = np.array(data['DataSets'][setKeys[i]])
				obsBs = np.array(data['DataSets'][setKeys[middle]])
				bsDiff = obsBs - modBs
				pltnmBs = outputname + setKeys[i] + '_measured'
				plotDifElev(pltnmBs, outDir, setKeys[i], x, bsDiff, 'r')
				i = i+1
				middle = middle+1
	writeTest = 0
	counter = counter+1
