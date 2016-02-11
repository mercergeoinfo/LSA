#!/Users/andrew/anaconda/bin/python
#from __future__ import division
"""Late Season Accumulation Model"""
__author__ = 'Andrew Mercer'
__contact__ = 'mercergeoinfo@gmail.com'
__maintainer__ = 'Andrew Mercer'
__version__ = '1.0'
__date__ = '29/01/2016'
#
# Created: 29/01/2016
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
from scipy import stats
from scipy.optimize import curve_fit
from scipy.stats.stats import nanmean
from scipy.optimize import minimize
from osgeo import gdal
#from osgeo import ogr
from osgeo import osr
from osgeo.gdalconst import *
## PLOT TOOLS
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
## OWN MODULES
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
		# Elevation difference to AWS.
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
		elif self.currSnowMass <= 0.0 and self.elevation > self.ELA:
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
	paramDict = {}
	try:
		paramDict['ddfSnow'] = float(settings['Snow'])
		paramDict['ddfSi'] = float(settings['Si'])
		paramDict['ddfFirn'] = float(settings['Firn'])
		paramDict['ddfIce'] = float(settings['Ice'])
		paramDict['lapse'] = float(settings['lapse'])
		paramDict['elevLapse'] = float(settings['elevLapse'])
		paramDict['sfe'] = float(settings['sfe'])
		paramDict['ELA'] = float(settings['ELA'])
	except:
		pass
	return refElev, jdayBw, jdatelist, startday, paramDict
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
def runModel(inData, stakeNames, temps, times, jdatelist, jdayBw, paramDict, trudbKeys, trudb, counter):
	data = copy.deepcopy(inData)
	data['DataSets'] = {}
	for stake in stakeNames:
		# For ordered headers/keys
		data[stake]['Headers'] = ['MeltModel', 'Shadevals', 'Easting', 'Northing', 'Elevation', 'Org_Bw']
		if 'Org_Bs' in data[stake].keys():
			data[stake]['Headers'].append('Org_Bs')
		if 'Org_Bn' in data[stake].keys():
			data[stake]['Headers'].append('Org_Bn')
		# Send input data to Degree Day Model object
		if	'Mod' + str(counter - 1) + '_Bw_' + str(jdatelist[0]) in data[stake].keys():
			BwColumn =	'Mod' + str(counter - 1) + '_Bw_' + str(jdatelist[0])
		else:
			BwColumn =	'Org_Bw'
		# print "Using column ", BwColumn
		data[stake]['MeltModel'] = 0
		data[stake]['MeltModel'] = DdfCell(data[stake]['Easting'], data[stake]['Northing'], data[stake]['Elevation'], data[stake][BwColumn], jdayBw, data[stake]['Shadevals'], paramDict)
		# For each julian day in the "times" vector call the meltInst method for each point object, passing the temperature and the day number.
		# This is what runs the model at each time step in the temperature time series file
		for i in range(len(temps)):
			data[stake]['MeltModel'].meltInst(temps[i],times[i])
		for day in jdatelist:
			# Fetch modelled melt and net balance for each julian day specific in settings and create new entry for each
			loc = data[stake]['MeltModel'].jTimeSeries.index(day)
			if 'Org_Bn_' + str(day) in data[stake].keys():
				orgBnColumn = 'Org_Bn_' + str(day)
			else:
				orgBnColumn = 'Org_Bn'
			data[stake]['Mod' + str(counter) + '_Bw_' + str(day)] = round((data[stake][orgBnColumn] + data[stake]['MeltModel'].meltSumSeries[loc]), 3)
			data[stake]['Diff' + str(counter) + '_Bw_' + str(day)] = round((data[stake][orgBnColumn] + data[stake]['MeltModel'].meltSumSeries[loc] - data[stake]['Org_Bw']), 3)
			data[stake]['Mod' + str(counter) + '_Bs_' + str(day)] =	 round(data[stake]['MeltModel'].meltSumSeries[loc],3)
			data[stake]['Mod' + str(counter) + '_Bn_' + str(day)] =	 round(data[stake]['MeltModel'].BnSeries[loc],3)
			data[stake]['Headers'].append('Mod' + str(counter) + '_Bw_' + str(day))
			data[stake]['Headers'].append('Diff' + str(counter) + '_Bw_' + str(day))
			data[stake]['Headers'].append('Mod' + str(counter) + '_Bs_' + str(day))
			data[stake]['Headers'].append('Mod' + str(counter) + '_Bn_' + str(day))
			# Fetch any truthing data available
			if str(day) in trudbKeys:
				try:
					loc = np.where(trudb[str(day)]['Stake']==stake)[0][0]
					data[stake]['Org_Bs_' + str(day)] = round(trudb[str(day)]['Bs'][loc],3)
					data[stake]['Org_Bn_' + str(day)] = round(trudb[str(day)]['Bn'][loc],3)
					# data[stake]['Mod' + str(counter) + '_Bw_' + str(day)] = round((data[stake]['Org_Bn_' + str(day)] + data[stake]['Mod' + str(counter) + '_Bs_' + str(day)]), 3)
				except:
					data[stake]['Org_Bs_' + str(day)] = np.nan
					data[stake]['Org_Bn_' + str(day)] = np.nan
					# data[stake]['Mod' + str(counter) + '_Bw_' + str(day)] = np.nan
				data[stake]['Headers'].insert(-4, 'Org_Bs_' + str(day))
				data[stake]['Headers'].insert(-4, 'Org_Bn_' + str(day))
				# Add values to lists for calculating score later
				if 'Mod' + str(counter) + '_Bs_' + str(day) not in data['DataSets'].keys():
					data['DataSets']['Mod' + str(counter) + '_Bs_' + str(day)] = []
				data['DataSets']['Mod' + str(counter) + '_Bs_' + str(day)].append(data[stake]['Mod' + str(counter) + '_Bs_' + str(day)])
				if 'Org_Bs_' + str(day) not in data['DataSets'].keys():
					data['DataSets']['Org_Bs_' + str(day)] = []
				data['DataSets']['Org_Bs_' + str(day)].append(data[stake]['Org_Bs_' + str(day)])

	return data
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
def reportCreate(data, paramDict):
	report = copy.deepcopy(paramDict)
	setKeys = data['DataSets'].keys()
	# Order all Mod first, then all Org
	setKeys.sort()
	start = 0
	end = len(setKeys)
	middle = end/2
	i = start
	while i < end/2:
		# Calculate Score
		modBs = np.array(data['DataSets'][setKeys[i]])
		obsBs = np.array(data['DataSets'][setKeys[middle]])
		modBsmean = nanmean(modBs)
		obsBsmean = nanmean(obsBs)
		obsBsMinModBs = obsBs - modBs
		obsBsMinMean = obsBs - obsBsmean
		SSres = (np.nansum(obsBsMinModBs**2))
		SStot = (np.nansum(obsBsMinMean**2))
		ResNorm = SSres**0.5
		report[(setKeys[i]+'_RN')] = ResNorm # Norm of residuals version
		i = i+1
		middle = middle+1
	return report, ResNorm
#
def reportWrite(report, outDir):
	reportFile = os.path.join(outDir,'Report.txt')
	with open (reportFile, 'w') as fp:
		for p in sorted(report.keys()):
			fp.write("%s,%2.4f\n" % (p, report[p]))
	return 0
#
def plotDifElev(outputname,outDir, title, x, y, colour, xlabel, ylabel, plttitle):
	'''Plot modelled data.
	To use: plotDifElev(outputname,outDir, title, x, y, colour)'''
	#
	matplotlib.rcParams['axes.grid'] = True
	matplotlib.rcParams['legend.fancybox'] = True
	matplotlib.rcParams['figure.figsize'] = 11.69, 8.27 # A4
	matplotlib.rcParams['savefig.dpi'] = 300
	plotName = outputname + '.pdf'
	pp1 = PdfPages(os.path.join(outDir,plotName))
	fig1 = plt.figure(1)
	ax1 = fig1.add_subplot(111)
	xmax = max(x)
	xmin = min(x)
	labelString = title
	ax1.plot(x,y, color=colour, marker='o', linestyle='None', label=labelString)
	matplotlib.pyplot.axes().set_position([0.04, 0.065, 0.8, 0.9])
	ax1.plot([xmin/1.01,xmax*1.01],[0,0],'-k')
	plt.axis([xmin/1.01, xmax*1.01, -2, 2])
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(plttitle)
	pp1.savefig(bbox_inches='tight')
	pp1.close()
	plt.close()
	return 0
#
def plotDif(outputname,outDir, title, x, y, colour):
	'''Plot modelled data.
	To use: plotDif(outputname,outDir, title, x, y, colour)'''
	#
	matplotlib.rcParams['axes.grid'] = True
	matplotlib.rcParams['legend.fancybox'] = True
	matplotlib.rcParams['figure.figsize'] = 11.69, 8.27 # A4
	matplotlib.rcParams['savefig.dpi'] = 300
	plotName = outputname + '.pdf'
	pp1 = PdfPages(os.path.join(outDir,plotName))
	fig1 = plt.figure(1)
	ax1 = fig1.add_subplot(111)
	xmax = max(x)
	xmin = min(x)
	ymax = max(y)
	ymin = min(y)
	labelString = title
	ax1.plot(x,y, color=colour, marker='o', linestyle='None', label=labelString)
	matplotlib.pyplot.axes().set_position([0.04, 0.065, 0.8, 0.9])
	ax1.legend(bbox_to_anchor=(0.0, 1), loc=2, borderaxespad=0.1, ncol=3, title = "Julian Day")
	ax1.plot([xmin,xmax],[xmin,xmax],'-k')
	plt.axis([xmin+0.1, xmax+0.1, ymin+0.1, ymax+0.1])
	plt.xlabel('Measured Melt (m.w.e.)')
	plt.ylabel('Modelled Melt (m.w.e.)')
	plt.title('Modelled against Measured')
	# plt.show()
	pp1.savefig(bbox_inches='tight')
	pp1.close()
	plt.close()
	return 0
#
def scoring(scores, paramDict, BsScore):
	'''Set scores for each parameter when iterating over multiple possible combinations'''
	parameterNames = ['ddfSnow', 'ddfSi', 'ddfFirn', 'ddfIce', 'lapse', 'elevLapse', 'sfe', 'ELA', 'refElev']
	for p in parameterNames:
		value = paramDict[p]
		if value not in scores[p].keys():
			scores[p][value] = BsScore
		elif value in scores[p].keys():
			scores[p][value] = scores[p][value] + BsScore
		else:
			print 'Scoring error'
	return scores
#
def usageUpdate(parUsage, paramDict, BsScore):
	parameterNames = ['ddfSnow', 'ddfSi', 'ddfFirn', 'ddfIce', 'lapse', 'elevLapse', 'sfe', 'ELA', 'refElev']
	for p in parameterNames:
		value = paramDict[p]
		parUsage[p].append(value)
	parUsage['BsScore'].append(BsScore)
	return parUsage
#
def parameterCheckWrite(outputDir, year, parUsage):
	paramFile = os.path.join(outputDir,'UsedParameters.csv')
	keys = ['BsScore', 'ddfSnow', 'ddfSi', 'ddfFirn', 'ddfIce', 'lapse', 'elevLapse', 'sfe', 'ELA', 'refElev']
	with open (paramFile, 'w') as fp:
		writer = csv.writer(fp, delimiter = ",")
		writer.writerow(keys)
		writer.writerows(zip(*[parUsage[key] for key in keys]))
	return 0
#
def snowCalc(data, stakeNames, temps, times, jdatelist, jdayBw, paramDict, trudbKeys, trudb, outputDir, outputname):
	# Settings and counters for multiple parameter values
	counter = 0
	#
	previousModBs = []
	modelDiffScore = 1
	while modelDiffScore > 0.01 and counter < 20:
		# RUN MODEL
		data = runModel(data, stakeNames, temps, times, jdatelist, jdayBw, paramDict, trudbKeys, trudb, counter)
		dataKeys = data.keys()
		dataKeys.sort()
		#
		# REPORT ON MODEL
		if len(data['DataSets']) > 0:
			report, resNorm = reportCreate(data, paramDict)
		#if resNorm < best:
		print "Snow back calculation run {0:d}".format(counter)
		print "Score: {0:2.3f}\n".format(resNorm)
		# Output model data to file
		pfnm = str(counter)
		outDir = makeDir(outputDir, pfnm)
		meltDataWrite(data, outDir)
		#	Write report to text file
		reportWrite(report, outDir)
		#	Plot model results
		x = []
		for stake in stakeNames:
			x.append(data[stake]['Elevation'])
		if len(data['DataSets']) > 0:
			setKeys = data['DataSets'].keys()
			# Order all Mod first, then all Org
			setKeys.sort()
			# Plot differences between measured and modelled and test difference to previous run
			start = 0
			end = len(setKeys)
			middle = end/2
			i = start
			while i < end/2:
				modBs = np.array(data['DataSets'][setKeys[i]])
				if len(previousModBs) == 0:
					previousModBs = copy.copy(modBs)
				else:
					if i == 0:
						diffModBs = modBs - previousModBs
						modelDiffScore = np.mean(diffModBs)
						# print "Average difference of modelled ablation to previous estimate: {0:2.3f}".format(modelDiffScore)
						# print diffModBs
						previousModBs = copy.copy(modBs)
						pltnmMosDiff = outputname + setKeys[i] + '_moddiff'
						xlabel = 'Elevation (m.a.s.l.)'
						ylabel = 'Model improvement (m w.e.)'
						plttitle = 'Model improvement on previous run against Elevation'
						plotDifElev(pltnmMosDiff, outDir, setKeys[i], x, diffModBs, 'c',xlabel, ylabel, plttitle)
				obsBs = np.array(data['DataSets'][setKeys[middle]])
				bsDiff = obsBs - modBs
				pltnmBs = outputname + setKeys[i] + '_measured'
				xlabel = 'Elevation (m.a.s.l.)'
				ylabel = 'Measured - Modelled Melt (m w.e.)'
				plttitle = 'Measured - Modelled against Elevation'
				plotDifElev(pltnmBs, outDir, setKeys[i], x, bsDiff, 'r',xlabel, ylabel, plttitle)
				# pltnmBsmm = outputname + setKeys[i] + '_modmeas'
				# plotDif(pltnmBsmm, outDir, setKeys[i], obsBs, modBs, 'b')
				i = i+1
				middle = middle+1
		counter = counter+1
	return resNorm
#
def  modelOptimiser(data, stakeNames, temps, times, jdatelist, jdayBw, param, refElev, trudbKeys, trudb):
	paramDict = {}
	paramDict['ddfSnow'] = param[0]
	paramDict['ddfSi'] = param[1]
	paramDict['ddfFirn'] = param[2]
	paramDict['ddfIce'] = param[3]
	paramDict['lapse'] = param[4]
	paramDict['elevLapse'] = param[5]
	paramDict['sfe'] = param[6]
	paramDict['ELA'] = param[7]
	paramDict['refElev'] = refElev
	data = runModel(data, stakeNames, temps, times, jdatelist, jdayBw, paramDict, trudbKeys, trudb, 0)
	report, resNorm = reportCreate(data, paramDict)
	return resNorm
#
#
################################### MAIN ###################################
#
# Try calling with e.g "echo 2005 2006 2007 2008 2009 2010 2011 2013 2014 | xargs ./lsa.py"
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
# /Output/Shades/SG_shade.tif		This file is here as it is the output of another script
# /Output/yyyy/		These are created by this script as needed
# /Scripts/
#
# Format examples:
	# settings.txt (note that if 'ExportDate' omitted then all dates between bW and final temperature reading exported):
		# Elevation=1150
		# jdayBw=115
		# ExportDate=236,256
		# ShadeStart=100
		# These parameter settings must ALL be present or those present will be ignored
		# Snow=0.0046
		# Si=0.0054
		# Firn=0.0058
		# Ice=:0.0064
		# lapse=0.0044
		# elevLapse=950
		# sfe=1.5
		# ELA=1500
	# weatheryyyy.csv:
		# Date,Temp
		# 2010-01-25,-8.3
	# StakeDatayyyy.csv:
		# Stake,Easting,Northing,Elevation,Bw,Bs,Bn,Surface
		# 04C,651103.586397,7536381.86553,1219,0.334,2.53,-2.196,ice
#
#
# Get shading data
# Location of file containing a multiband raster, each band represents the shading on one day. 1 = no shade, 0 = really quite dark
shadefile = '../InData/Shades/SG_shade.tif'
# shadefile = '../InData/Shades/Reduced/SG_shade.tif'
# Read the shade factor raster in to memory
raster, transf, bandcount = getShadeFile(shadefile)
#
#
def main():
# Set up list of years from command line arguments
	paramChoice = ''
	if len(sys.argv) > 1:
		argstart = 1
		if sys.argv[1] == 'p':
			paramChoice = 'p'
			argstart = 2
		years = []
		for arg in sys.argv[argstart:]:
			try:
				years.append(int(arg))
			except:
				print sys.argv
				sys.exit("Argument Error")
		print years
	else:
		years = []
		while len(years) <1:
			try:
				yearsFromUser = raw_input("Give years as space separated list: ")
				yearsAsString = yearsFromUser.strip().split(' ')
				for yearString in yearsAsString:
					years.append(int(yearString))
			except:
				years = []
#
# Run main function
	for year in years:
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
		refElev, jdayBw, jdatelist, startday, paramDict = getSettings(dataLoc, times[-1])
		print "For year %s following settings used: " %(strYear)
		print "refElev set to %s" %(refElev)
		print "jdayBw set to %s" %(jdayBw)
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
				trudb[file.split('.')[0]] = import2vector(os.path.join(truthDir,file),	importError = 'no')
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
		paramDict['refElev'] = refElev

		# Prepare database for passing into loop and model run
		stakeNames = stakeData.keys()
		stakeNames.sort()
		data = copy.deepcopy(stakeData)
		best = 9999
		#

		if paramDict.keys() > 0 and paramChoice == 'p':
			resNorm = snowCalc(data, stakeNames, temps, times, jdatelist, jdayBw, paramDict, trudbKeys, trudb, outputDir, outputname)
		else:
			x0 = [0.0050, 0.0050, 0.0050, 0.0050, 0.0050, 1000, 1.0, 1500]
			def rosen(x):
				resNorm = modelOptimiser(data, stakeNames, temps, times, jdatelist, jdayBw, x, refElev, trudbKeys, trudb)
				return resNorm
			res = minimize(rosen, x0, method='nelder-mead', options={'xtol': 1e-8, 'disp': True})
			print (res.x)
			paramDict['ddfSnow'] = res.x[0]
			paramDict['ddfSi'] = res.x[1]
			paramDict['ddfFirn'] = res.x[2]
			paramDict['ddfIce'] = res.x[3]
			paramDict['lapse'] = res.x[4]
			paramDict['elevLapse'] = res.x[5]
			paramDict['sfe'] = res.x[6]
			paramDict['ELA'] = res.x[7]
			resNorm = snowCalc(data, stakeNames, temps, times, jdatelist, jdayBw, paramDict, trudbKeys, trudb, outputDir, outputname)
		print "Final score: {:2.4f}".format(resNorm)

#
if __name__ == "__main__":
	main()