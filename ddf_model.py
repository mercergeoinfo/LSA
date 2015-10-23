#!/Users/andrew/anaconda/bin/python
#from __future__ import division
"""Degree Day Melt Model"""
__author__ = 'Andrew Mercer'
__contact__ = 'mercergeoinfo@gmail.com'
__maintainer__ = 'Andrew Mercer'
__version__ = '2.0'
__date__ = '17/03/2015'
#
# Created: 16/04/2014
# Edits: 16/01/2015
#		20/01/2015
#		21/01/2015
#		26/01/2015
#		03/02/2015
#		05/02/2015
#		16/02/2015
#		17/02/2015
# Version: 2.0
# Andrew Mercer
# This file runs a degree day melt model on test data from Storglaciaren
# It is a "working" version, running data from the probing grid
#
from datetime import datetime
import sys
import os
import shutil
import fnmatch
import csv
import struct
## SCIENTIFIC
#import math
import numpy as np
#import numpy.ma as ma
#import scipy.stats as stats
from scipy.stats.stats import nanmean
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
from osgeo.gdalconst import *
## PLOT TOOLS
import matplotlib
import matplotlib.pyplot as plt
#import matplotlib.dates as dates
#from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
import ParameterAnalysis
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
def namer(pathstring):
	'''Get name and extension of a file
	To use: name,ext,path,namefull = namer(pathstring) '''
	## Called by:
	namefull = os.path.basename(pathstring)
	path = os.path.dirname(pathstring)
	ext = namefull.split('.')[1]
	name = namefull.split('.')[0]
	return name,ext,path,namefull
#
def makeDir(path,name):
	'''Make a directory at given location
	To use: outDir = makeDir(path, name)'''
	## Called by:
	targetDir = os.path.join(path,name)
	if not os.path.exists(targetDir):
		os.makedirs(targetDir)
	return targetDir
#
def makeDirTS(path,name):
	'''Make a directory at given location with current time stamp added to name
	To use: outDir, containDir = makeDirTS(path, name)'''
	## Called by:
	tempDir = os.path.join(path,name)
	namstr = datetime.now().strftime('%j_%H%M%S')
	targetDir = os.path.join(tempDir,namstr)
	if not os.path.exists(targetDir):
		os.makedirs(targetDir)
	return targetDir, tempDir
#
def filelist(folder,ending):
	"""Return list of specific file type in folder other than current folder.
	To use: filelist = filelist(folder,ending)"""
	# Called by:
	matchstring = '*.' + ending
	filelist = fnmatch.filter(os.listdir(folder),matchstring)
	print 'From ', folder, ' the following matched ', matchstring, '\n', filelist
	return filelist
#
def pt2fmt(pt):
	'''Format types for data extraction from shade file'''
	fmttypes = {
		GDT_Byte: 'B',
		GDT_Int16: 'h',
		GDT_UInt16: 'H',
		GDT_Int32: 'i',
		GDT_UInt32: 'I',
		GDT_Float32: 'f',
		GDT_Float64: 'f'
		}
	return fmttypes.get(pt, 'x')
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
	srs=osr.SpatialReference(wkt=projec)
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
	#print x, y
	success, transfInv = gdal.InvGeoTransform(transf)
	if not success:
		print "Failed InvGeoTransform()"
		sys.exit(1)
	rows = raster.RasterYSize
	cols = raster.RasterXSize
	xpix, ypix = gdal.ApplyGeoTransform(transfInv, x, y)
	# Read the file band to a matrix called band_1
	for i in range(1,bandcount+1):
		band = raster.GetRasterBand(i)
		bandtype = gdal.GetDataTypeName(band.DataType)
		if band is None:
			continue
		structval = band.ReadRaster(int(xpix), int(ypix), 1,1, buf_type = band.DataType )
		fmt = pt2fmt(band.DataType)
		intval = struct.unpack(fmt , structval)
		vals[i] = intval[0]
	return vals
#
def import2vector(fileName):
	"""Imports the data as vectors in a dictionary."""
	InFile = open(fileName,'rb')
	Name = os.path.basename(fileName).split('.')[0]
	# Check if file has iButton header
	line = InFile.next()
	Headers = line.strip().split(',')
	data = {}
	data['Source'] = Name
	i=0
	#
	nanlist = ['NAN','NaN','nan','NULL','null','-9999',-9999,'']
	for i in Headers:
		data[i] = []
	for row in InFile:
		if row != "\n":
			# Split read line of data into list
			dataIn = row.strip().split(',')
			for i in range(len(dataIn)):
				if dataIn[i] in nanlist:
				#if dataIn[i] == 'nan' or dataIn[i] == 'NaN' or dataIn[i] == 'NULL' :
					 dataIn[i] = np.nan
					 #data[Headers[i]].append(dataIn[i])
				elif dataIn[i] == "":
					dataIn[i] = np.nan
					#data[Headers[i]].append(dataIn[i])
				else:
					try:
						dataIn[i] = datetime.strptime(dataIn[i],'%d/%m/%y %H:%M:%S')
					except:
						try:
							dataIn[i] = float(dataIn[i])
						except:
							dataIn[i] = dataIn[i]
				data[Headers[i]].append(dataIn[i])
	for i in Headers:
		data[i] = np.array(data[i])
	return data
#
def meltDataOut(pointKeys, points, outDir):
	'''Write results of melt model to csv file
	To use: meltDataOut(pointKeys, points, outDir)'''
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
	'''Plot modelled data '''
	# Plot the modelled data
	matplotlib.rcParams['axes.grid'] = True
	matplotlib.rcParams['legend.fancybox'] = True
	#matplotlib.rcParams['figure.figsize'] = 18, 9 #Mine
	#matplotlib.rcParams['figure.figsize'] = 16.54, 11.69 #A3
	matplotlib.rcParams['figure.figsize'] = 11.69, 8.27 #A4
	matplotlib.rcParams['savefig.dpi'] = 300
	plotName = outputname + '.pdf'
	pp1 = PdfPages(os.path.join(outDir,plotName))
	fig1 = plt.figure(1)
	ax1 = fig1.add_subplot(111)
	lnclr = ['k','r','g','b','y','c','m']
	lnsty = ['-','--','-.',':']
	mrsty = ['o','s','v','*','x','+','1','2','3','4']
	xmax = max(x)
	xmin = min(x)
	ymax = 3
	ymin = -3
	ax1.plot(x,y, color=colour, marker='o', linestyle='o', label=outputname)
	matplotlib.pyplot.axes().set_position([0.04, 0.065, 0.8, 0.9])
	ax1.legend(bbox_to_anchor=(0.0, 0), loc=3, borderaxespad=0.1, ncol=3, title = "Zone, Day")
	ax1.plot([xmin,xmax],[0,0],'-k')
	#lastind = str(jdatelist[-1])
	#ax1.plot([0,5],np.multiply(dbvecs[lastind]['slope'],[0,5]) + dbvecs[lastind]['intercept'],'-b')
	#plt.axis([xmin, xmax, ymin, ymax*1.2])
	plt.axis([xmin, xmax, -2, 2])
	# for pnt in range(len(dbvecs[jd]['Stake'])):
# 		ax1.annotate(dbvecs[jd]['Stake'][pnt],xy=(dbvecs[jd]['Elevation'][pnt],1.8), rotation=90)
	plt.xlabel('Elevation (m.a.s.l.')
	plt.ylabel('Measured - Modelled Melt (m w.e.)')
	plt.title('Measured - Modelled against Elevation')
	#plt.show()
	pp1.savefig(bbox_inches='tight')
	pp1.close()
	plt.close()
	return 0
#
#
# MAIN
#
time_zero = datetime.now()
time_one = time_zero
# Get shading data
# Location of file containing a multiband raster, each band represents the shading on one day. 1 = no shade, 0 = really quite dark
shadefile = '../Output/Shades/SG_shade.tif'
# Read the shade factor raster in to memory
raster, transf, bandcount = getShadeFile(shadefile)
# Give the Julian day at which the shade raster starts
startday = 100
# [2009,2010,2011,2012,2013]
# Choose which data set is to be run. This section could be simplified and shortened by introducing a naming convention for the raw data file
for year in [2012] :
	#
	if year == 2009:
		# Temperature data: Following two lines are example of input format
		# Date,Temp
		# 2010-01-25,-8.3
		TfileName = '../InData/SMHI/SMHI2009.csv' # 2009
		refElev = 1150
		# Stake data: Following two lines are example of input format
		# Stake,Easting,Northing,Elevation,Bw,Bs,Bn,Surface
		# 04C,651103.586397,7536381.86553,1219,0.334,2.53,-2.196,ice
		SfileName = '../InData/2009/256.csv' # 2010
		# Give day of last winter probing (julian day for Bw)
		jdayBw = 115 # 2009
		# Output melt sum at julian days to csv
		jdatelist = [256] # 2009
		# Directory for output
		outputDir = "../Output/2009/" # 2009
		outputname = '2009_DDM' # 2009
		truthDir = '../InData/2009/' # 2009
		truthfiles = filelist(truthDir,'csv')
	elif year == 2010:
		# Temperature data: Following two lines are example of input format
		# Date,Temp
		# 2010-01-25,-8.3
		TfileName = '../InData/SMHI/SMHI2010.csv' # 2010
		refElev = 1150
		# Stake data: Following two lines are example of input format
		# Stake,Easting,Northing,Elevation,Bw,Bs,Bn,Surface
		# 04C,651103.586397,7536381.86553,1219,0.334,2.53,-2.196,ice
		SfileName = '../InData/2010/256.csv' # 2010
		# Give day of last winter probing (julian day for Bw)
		jdayBw = 105 # 2010
		# Output melt sum at julian days to csv
		jdatelist = [256] # 2010
		# Directory for output
		outputDir = "../Output/2010/" # 2010
		outputname = '2010_DDM'
		truthDir = '../InData/2010/' # 2010
		truthfiles = filelist(truthDir,'csv')
	elif year == 2011:
		# Temperature data: Following two lines are example of input format
		# Date,Temp
		# 2010-01-25,-8.3
		TfileName = '../InData/SMHI/TRS2011.csv' # 2011
		refElev = 1150
		# Stake data: Following two lines are example of input format
		# Stake,Easting,Northing,Elevation,Bw,Bs,Bn,Surface
		# 04C,651103.586397,7536381.86553,1219,0.334,2.53,-2.196,ice
		SfileName = '../InData/2011/257.csv' # 2011
		# Give day of last winter probing (julian day for Bw)
		jdayBw = 115 # 2010
		# Output melt sum at julian days to csv
		jdatelist = [257] # 2010
		# Directory for output
		outputDir = "../Output/2011/" # 2011
		outputname = '2011_DDM'
		truthDir = '../InData/2011/' # 2011
		truthfiles = filelist(truthDir,'csv')
	elif year == 2012:
		# Temperature data:  Following two lines are example of input format
		# Date,Temp
		# 2010-01-25,-8.3
		TfileName = '../InData/SMHI/SMHI2012.csv' # 2012
		refElev = 1150
		# Stake data: Following two lines are example of input format
		# Stake,Easting,Northing,Elevation,Bw,Bs,Bn,Surface
		# 04C,651103.586397,7536381.86553,1219,0.334,2.53,-2.196,ice
		#SfileName = '../InData/2012/257.csv' # 2012
		SfileName = '../InData/2012/2012melt_r4.csv' # 2012
		# Give day of last winter probing (julian day for Bw)
		jdayBw = 115 # 2012
		# Output melt sum at julian days to csv
		#jdatelist = [247] # 2012
		jdatelist = [186,204,226,247,257] # 2012
		# Directory for output
		outputDir = "../Output/2012/" # 2012
		outputname = '2012_DDM'
		truthDir = '../InData/2012/' # 2012
		truthfiles = filelist(truthDir,'csv')
	elif year == 2013:
		# Temperature data:  Following two lines are example of input format
		# Date,Temp
		# 2010-01-25,-8.3
		TfileName = '../InData/SMHI/SMHI2013.csv' # 2013
		refElev = 1150
		# Stake data: Following two lines are example of input format
		# Stake,Easting,Northing,Elevation,Bw,Bs,Bn,Surface
		# 04C,651103.586397,7536381.86553,1219,0.334,2.53,-2.196,ice
		SfileName = '../InData/2013/256.csv' # 2013
		# Give day of last winter probing (julian day for Bw)
		jdayBw = 114 # 2013
		# Output melt sum at julian days to csv
		#jdatelist = [256] # 2013
		jdatelist = [194,198,211,224,256] # 2013
		# Directory for output
		outputDir = "../Output/2013/" # 2013
		outputname = '2013_DDM'
		truthDir = '../InData/2013/' # 2013
		truthfiles = filelist(truthDir,'csv')
	#
	if not os.path.exists(outputDir):
		os.makedirs(outputDir)
	#
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
	#
	#
	# Set parameters for the melt model
	# Set test to zero to run a single set of parameters and don't compare results with measured data
	test = 0
	if test == 0:
		ddfSnow = [0.0038]
		ddfSi = [0.0051]
		ddfFirn = [0.0050]
		ddfIce = [0.0054]
		lapse = [0.0055]
	elif test == 1:
		ddfSnow = [0.0030,0.0032,0.0034,0.0036,0.0038,0.0040,0.0042,0.0044,0.0046,0.0048]
		ddfSi = [0.0045,0.0047,0.0049,0.0051]
		ddfFirn = [0.0043,0.0045,0.0048,0.0050,0.0053,0.0055,0.0058]
		ddfIce = [0.0043,0.0045,0.0047,0.0049,0.0051,0.0053,0.0055,0.0057]
		lapse = [0.0042,0.0044,0.0046,0.0048,0.0050,0.0055,0.0057,0.0059,0.0061,0.0063,0.0065,0.0068,0.0070,0.0072,0.0074,0.0076]
	else:
		ddfSnow = [0.0038]
		ddfSi = [0.0051]
		ddfFirn = [0.0050]
		ddfIce = [0.0054]
		lapse = [0.0055]
	# The parameters below are fixed for this study
	elevLapse = (2100 - 1150) # Elevation dependant lapse rate
	sfe = 1.5 # Shading factor exponent (adjusts the shading value at each point)
	ELA = 1500 # Equilibrium line, for firn or ice under snow
	# Calculate how many iterations are to be performed (feedback for user)
	Snow = range(len(ddfSnow))
	Super = range(len(ddfSi))
	Firn = range(len(ddfFirn))
	Ice = range(len(ddfIce))
	Lapse = range(len(lapse))
	pltot = len(ddfSnow)*len(ddfSi)*len(ddfFirn)*len(ddfIce)*len(lapse)
	print "Year: %d Iterations: %d" % (year,pltot)
	if test != 0:
		trudb = {} # FOR VERIFICATION OF RESULTS ONLY
		for file in truthfiles: # FOR VERIFICATION OF RESULTS ONLY
			trudb[file.split('.')[0]] = import2vector(os.path.join(truthDir,file)) # FOR VERIFICATION OF RESULTS ONLY
	bestBsR2 = np.nan
	bestBnR2 = np.nan
	rsq1 = []
	rsq2 = []
	counterv = []
	bestDir = []
	parUsage = {}
	parUsage['ddfSnow'] = []
	parUsage['ddfSi'] = []
	parUsage['ddfFirn'] = []
	parUsage['ddfIce'] = []
	parUsage['lapse'] = []
	parUsage['BsR2'] = []
	parUsage['BnR2'] = []
	counter = 0
	for Snw in Snow:
		for Spr in Super:
			for Frn in Firn:
				for Eci in Ice:
					for Lps in Lapse:
						# Parameters read into dictionary to be passed to model
						paramDict = {}
						paramDict['ddfSnow'] = ddfSnow[Snw]
						paramDict['ddfSi'] = ddfSi[Spr]
						paramDict['ddfFirn'] = ddfFirn[Frn]
						paramDict['ddfIce'] = ddfIce[Eci]
						paramDict['lapse'] = lapse[Lps]
						paramDict['elevLapse'] = elevLapse
						paramDict['sfe'] = sfe
						paramDict['ELA'] = ELA
						paramDict['refElev'] = refElev
						#
						# Create vector for shade values
						vals = []
						for d in range(365):
							vals.append(1)
						#
						# Read stake data
						SinFile = open(SfileName,'rb')
						points = {}
						# Read point data file with position and winter balance as of last winter probing (jdayBw) and send to model
						for line in csv.DictReader(SinFile, delimiter=','):
							data = {}
							#data['name'] = line['Stake'].strip()
							# Coordinates
							data['easting'] = float(line['Easting'].strip())
							data['northing'] = float(line['Northing'].strip())
							data['elevation'] = float(line['Elevation'].strip())
							# Get shading factor for location
							try:
								data['Shadevals'] = GetShadeVals(data['easting'], data['northing'], raster, transf, bandcount, vals, startday)
								#print line['Stake'].strip()
							except:
								data['Shadevals'] = vals
							# Get the winter balance
							data['Bw'] = float(line['Bw'].strip())
							# Send input data to Degree Day Model object
							data['MeltModel'] = DdfCell(data['easting'], data['northing'], data['elevation'], data['Bw'], jdayBw, data['Shadevals'], paramDict)
							points[line['Stake'].strip()] = data
						pointKeys = points.keys()
						pointKeys.sort()
						#
						# For each julian day in the "times" vector call the meltInst method for each point object, passing the temperature and the day number.
						# This is what runs the model at each time step in the temperature time series file
						for i in range(len(temps)):
							for j in pointKeys:
								points[j]['MeltModel'].meltInst(temps[i],times[i])
						#
						# Only write out results now if not running test, otherwise wait until assessed
						if test == 0:
							# Output model data to file
							# Get name of mb file
							inname,inext,inpath,innamefull = namer(SfileName)
							outDir = makeDir(outputDir, inname)
							meltDataOut(pointKeys, points, outDir)
							# Write parameters used to text file
							paramFile = os.path.join(outDir,'Parameters.txt')
							with open (paramFile, 'w') as fp:
								for p in paramDict.items():
									fp.write("%s:%2.6f\n" % p)
						#
						# MARKER - end of modelling
						#
						if test != 0:
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
								# Calculate 'errors'
								modBs = dbcomp[bsn]
								modBn = dbcomp[bnn]
								obsBsmean = nanmean(obsBs)
								obsBnmean = nanmean(obsBn)
								obsBsMinModBs = obsBs - modBs
								obsBnMinModBn = obsBn - modBn
								obsBsMinMean = obsBs - obsBsmean
								obsBnMinMean = obsBn - obsBnmean
								BsR2 = 1 - ((np.nansum(obsBsMinModBs**2)) / (np.nansum(obsBsMinMean**2)))
								BnR2 = 1 - ((np.nansum(obsBnMinModBn**2)) / (np.nansum(obsBnMinMean**2)))
								paramDict[(bsn+'_R2')] = BsR2
								paramDict[(bnn+'_R2')] = BnR2
								bsR2list.append(BsR2)
								bnR2list.append(BnR2)
								# Add parameters used to results storage
								parUsage['ddfSnow'].append(paramDict['ddfSnow'])
								parUsage['ddfSi'].append(paramDict['ddfSi'])
								parUsage['ddfFirn'].append(paramDict['ddfFirn'])
								parUsage['ddfIce'].append(paramDict['ddfIce'])
								parUsage['lapse'].append(paramDict['lapse'])
								parUsage['BsR2'].append(BsR2)
								parUsage['BnR2'].append(BnR2)
						#if max(bsR2list) >= bestBsR2 and max(bnR2list) >= bestBnR2:
						if test == 0:
							continue
						else:
							if test != 0 and ((max(bsR2list) < bestBsR2) or np.isnan(bestBsR2)):
								time_i = datetime.now()
								print "\n", time_i
								if (max(bsR2list) > bestBsR2) or np.isnan(bestBsR2):
									bestBsR2 = max(bsR2list)
									print "new best Bs R2: %2.3f" % bsR2list[-1]
								if (max(bnR2list) > bestBnR2)or np.isnan(bestBnR2):
									bestBnR2 = max(bnR2list)
									print "new best Bn R2: %2.3f" % bnR2list[-1]
								#
								# Output model data to file
								# Get name of mb file
								# inname,inext,inpath,innamefull = namer(SfileName)
								flnm = str(counter)
								outDir = makeDir(outputDir, flnm)
								meltDataOut(pointKeys, points, outDir)
								# Write parameters used to text file
								paramFile = os.path.join(outDir,'Parameters.txt')
								with open (paramFile, 'w') as fp:
									for p in paramDict.items():
										fp.write("%s:%2.6f\n" % p)
								#
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
									pltnmBs = 'Bs_diff_' + str(d)
									pltnmBn = 'Bn_diff_' + str(d)
									plotDifElev(pltnmBs,outDir,x,bsDiff,'r')
									plotDifElev(pltnmBn,outDir,x,bnDiff,'k')
							pltot = pltot -1
							counter = counter +1
	# For test runs, write out csv list of R squared for parameter values used
	if test != 0:
		paramFile = os.path.join(outputDir,'UsedParameters.csv')
		keys = ['BsR2','BnR2','ddfSnow','ddfSi','ddfFirn','ddfIce','lapse']
		with open (paramFile, 'w') as fp:
			writer = csv.writer(fp, delimiter = ",")
			writer.writerow(keys)
			writer.writerows(zip(*[parUsage[key] for key in keys]))
		ParamDir = makeDir('../Edit/', 'Parameters')
		shutil.copyfile(paramFile,os.path.join(ParamDir,(str(year)+'UsedParameters.csv')))
# If test run, do analysis of the best values
if test != 0:
	ParameterAnalysis.main()
