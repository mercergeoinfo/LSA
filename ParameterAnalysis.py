#!/Users/andrew/anaconda/bin/python
"""Description of file contents"""
__author__ = 'Andrew Mercer'
__contact__ = 'mercergeoinfo@gmail.com'
__maintainer__ = 'Andrew Mercer'
__version__ = '0.1.1'
__date__ = '15/07/2014'
#
## IMPORTS
import sys
import os
# import glob
import fnmatch
# import gc
# import csv
# import pickle
# import time
# from datetime import datetime, timedelta
#
## SCIENTIFIC
import math
import numpy as np
# import numpy.ma as ma
import scipy.stats as stats
from scipy.stats.stats import nanmean
# from pandas import Series, DataFrame
# import pandas as pd
#
## GIS SUPPORT
# import gdal
# import gdalconst
#
## PLOT TOOLS
# import matplotlib
# import matplotlib.pyplot as plt
# import matplotlib.dates as dates
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib.backends.backend_pdf import PdfPages
# import pylab
# pylab.ion()
# pylab.ioff()
#
## OWN
#import AM_Func as AM
#
## FUNCTIONS
#
def namer(pathstring):
	'''Get name and extension of a file
	To use: name,ext,path,namefull = namer(pathstring) '''
	## Called by: kriging
	namefull = os.path.basename(pathstring)
	path = os.path.dirname(pathstring)
	ext = namefull.split('.')[1]
	name = namefull.split('.')[0]
	return name,ext,path,namefull
#
def filelist(folder,ending):
	'''Return list of specific file type in folder other than current folder.'''
	matchstring = '*.' + ending
	filelist = fnmatch.filter(os.listdir(folder),matchstring)
	print 'From ', folder, ' the following matched ', matchstring, '\n', filelist
	return filelist
#
def import2vector(fileName, dateString = '%d/%m/%y %H:%M:%S'):
	'''Imports the data as vectors in a dictionary.'''
	# Open file
	InFile = open(fileName,'rb')
	line = InFile.next()
	Name = os.path.basename(fileName).split('.')[0]
	# Get headers
	Headers = line.strip().split(',')
	# Create dictionary for data
	data = {}
	data['Source'] = Name
	i=0
	# Set up list of Nan values
	nanlist = ['NAN','NaN','nan','NULL','null','-9999',-9999]
	# Read through data file
	for i in Headers:
		data[i] = []
	for row in InFile:
		if row != "\n":
			# Split read line of data into list
			dataIn = row.strip().split(',')
			for i in range(len(dataIn)):
				# Check for NaN and empty values
				if dataIn[i] in nanlist:
					 dataIn[i] = np.nan
				elif dataIn[i] == "":
					dataIn[i] = np.nan
				else:
					# Try date formatted data conversion
					try:
						dataIn[i] = datetime.strptime(dataIn[i],dateString)
					except:
						# Try converting to float
						try:
							dataIn[i] = float(dataIn[i])
						except:
							# Leave as string
							dataIn[i] = dataIn[i]
				# Add to vector
				data[Headers[i]].append(dataIn[i])
	for i in Headers:
		# Convert to numpy arrays
		data[i] = np.array(data[i])
		try:
			# Create posts containing basic statistics for each numerical column (vector)
			data[(str(i)+'_min')] = np.nanmin(data[i])
			data[(str(i)+'_max')] = np.nanmax(data[i])
			data[(str(i)+'_mean')] = np.nanmean(data[i])
			data[(str(i)+'_stdDev')] = np.nanstd(data[i])
			data[(str(i)+'_median')] = np.median(data[i])
		except:
			print "\nStatistics not computable for %s\n" % str(i)
	return data
#
## MAIN
#
def main():
	# Read in parameter usage files (csv) from folder into dictionary of vectors
	inFolder = '../Edit/Parameters'
	pFileList = filelist(inFolder,'csv')
	pDict = {}
	for file in pFileList:
		fileIn = os.path.join(inFolder,file)
		name,ext,path,namefull = namer(fileIn)
		pDict[name] = import2vector(fileIn)
	keys1 = sorted(pDict.keys())
	# From the Bs R^2 value, calculate weights
	parKeys1 = ['ddfSnow', 'ddfSi', 'ddfFirn','ddfIce', 'lapse']
	prtxt = open('../Edit/Parameters.txt', "wb")
	for key in keys1:
		prtxt.write("Year,BsR2,BnR2,Parameter,Mean,s,WeightMean,Best\n")
		pDict[key]['Weight'] = pDict[key]['BsR2'] / pDict[key]['BsR2_mean']
	# Calculate weighted mean for each parameter
		for parKey in parKeys1:
			pDict[key][('Weighted_'+parKey)] = pDict[key][parKey] * pDict[key]['Weight']
			pDict[key][('Weighted_'+parKey+'_mean')] = np.nanmean(pDict[key][('Weighted_'+parKey)])
			ind = np.where(pDict[key]['BsR2'] == pDict[key][('BsR2_max')])
			for i in ind[0]:
				prtxt.write("%s,%2.3f,%2.3f,%s,%2.4f,%2.4f,%2.4f,%2.4f\n" % (key[:4], pDict[key][('BsR2_max')],pDict[key][('BnR2_max')], parKey, pDict[key][(parKey+'_mean')], pDict[key][(parKey+'_stdDev')], pDict[key][('Weighted_'+parKey+'_mean')], pDict[key][parKey][i]))
	prtxt.close()
#
#
if __name__ == "__main__":
	 main() # Calls first function, named "main" which contains main body of programme. For use as library file.