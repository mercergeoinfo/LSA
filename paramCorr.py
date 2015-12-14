#!/Users/andrew/anaconda/bin/python
# -*- coding: utf-8 -*-
#from __future__ import division
"""Calculate correlation coefficient and covariance matrices"""
__author__ = 'Andrew Mercer'
__contact__ = 'mercergeoinfo@gmail.com'
__maintainer__ = 'Andrew Mercer'
__version__ = '1.0.0'
__date__ = '03/11/2015'
#
## IMPORTS
import sys
import os
# import glob
# import fnmatch
# import gc
import csv
# import pickle
# import time
# from datetime import datetime, timedelta
# import subprocess
# import urllib2
# import os.path
# import zipfile
# import StringIO
#
## SCIENTIFIC
import math
from math import *
import numpy as np
import numpy.ma as ma
# import scipy.linalg
import scipy.stats as stats
from scipy.stats.stats import nanmean
# from scipy.stats import norm
# from pandas import Series, DataFrame
# import pandas as pd
#
## GIS SUPPORT
# import gdal
# import gdalconst
# from geostatsmodels import utilities, kriging, variograms, model, geoplot
# from qgis.core import QgsProject
# from PyQt4.QtCore import QFileInfo
# To run in QGIS enter following in QGIS console:
# execfile(u'/path/to/script.py'.encode('utf-8'))
#
## PLOT TOOLS
# import matplotlib
# import matplotlib.pyplot as plt
# import matplotlib.dates as dates
# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib.backends.backend_pdf import PdfPages
# import pylab
# pylab.ion()
# pylab.ioff()
#
## OWN
# import AM_Func as AM
#
## FUNCTIONS
#
def getClipboardData():
	'''Get contents of system clipboard
	To use: data = getClipboardData()'''
	## Called by:
	p = subprocess.Popen(['pbpaste'], stdout=subprocess.PIPE)
	retcode = p.wait()
	data = p.stdout.read()
	return data
#
def setClipboardData(data):
	'''Set contents of system clipboard
	To use: setClipboardData(data)'''
	## Called by:
	p = subprocess.Popen(['pbcopy'], stdin=subprocess.PIPE)
	p.stdin.write(data)
	p.stdin.close()
	retcode = p.wait()
	return 0
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
def import2vector(fileName, dateString = '%d/%m/%y %H:%M:%S'):
	'''Imports the data as vectors in a dictionary. dateString is optional and can be set to match datetime format
	To use: db = import2vector(filename) or db = import2vector(filename, dateString = '%d/%m/%y %H:%M:%S')'''
	# Called by:
	# Open file
	InFile = open(fileName,'rb')
	line = InFile.next()
	Name = os.path.basename(fileName).split('.')[0]
	# Get headers
	Headers = line.strip().split(',')
	# Create dictionary for data
	data = {}
	data['Data'] = {}
	data['Description'] = {}
	data['Description']['Source'] = Name
	data['Description']['Numerical'] = []
	i=0
	# Set up list of Nan values
	nanlist = ['NAN','NaN','nan','NULL','null','-9999',-9999,'']
	# Read through data file
	for i in Headers:
		data['Data'][i] = []
	for row in InFile:
		if row != "\n":
			# Split read line of data into list
			dataIn = row.strip().split(',')
			for i in range(len(dataIn)):
				# Check for NaN and empty values
				if dataIn[i] in nanlist:
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
				data['Data'][Headers[i]].append(dataIn[i])
	for i in Headers:
		# Convert to numpy arrays
		data['Data'][i] = np.array(data['Data'][i])
		try:
			# Create posts containing basic statistics for each numerical column (vector)
			data['Description'][(str(i)+'_min')] = np.nanmin(data['Data'][i])
			data['Description'][(str(i)+'_max')] = np.nanmax(data['Data'][i])
			data['Description'][(str(i)+'_mean')] = np.nanmean(data['Data'][i])
			data['Description'][(str(i)+'_stdDev')] = np.nanstd(data['Data'][i])
			data['Description'][(str(i)+'_median')] = np.median(data['Data'][i])
			data['Description']['Numerical'].append(i)
		except:
			print "\nStatistics not computable for %s\n" % str(i)
	return data
#
## MAIN
#
file = "UsedParameters.csv"
data = import2vector(file)
x1 = data['Data']['BsR2']
x2 = data['Data']['BnR2']
y1 = data['Data']['ddfSnow']
y2 = data['Data']['ddfSi']
y3 = data['Data']['ddfFirn']
y4 = data['Data']['ddfIce']
y5 = data['Data']['lapse']
X = np.vstack((x1, x2, y1, y2, y3, y4, y5))
covX = np.cov(X)
corX = np.corrcoef(X)
#
reportFile = 'CovarianceReport.csv'
keys = [' ','BsR2','BnR2','ddfSnow','ddfSi','ddfFirn','ddfIce','lapse']
with open (reportFile, 'w') as fp:
	writer = csv.writer(fp, delimiter = ",")
	writer.writerow(keys)
	m, n = np.shape(covX)
	for i in range(m):
		row = [keys[i+1]]
		for j in range(n):
			#row.append(covX[i,j])
			row.append('{:7.3f}'.format(covX[i,j]))
		writer.writerow(row)
#
reportFile = 'CorrelationReport.csv'
keys = [' ','BsR2','BnR2','ddfSnow','ddfSi','ddfFirn','ddfIce','lapse']
with open (reportFile, 'w') as fp:
	writer = csv.writer(fp, delimiter = ",")
	writer.writerow(keys)
	m, n = np.shape(corX)
	for i in range(m):
		row = [keys[i+1]]
		for j in range(n):
			#row.append(corX[i,j])
			row.append('{:7.3f}'.format(corX[i,j]))
		writer.writerow(row)

# main():
# if __name__ == "__main__":
#	 main() # Calls first function, named "main" which contains main body of programme.

