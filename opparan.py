#!/Users/andrew/anaconda/bin/python
# -*- coding: utf-8 -*-
#from __future__ import division
"""Examine parameters from optimising model"""
__author__ = 'Andrew Mercer'
__contact__ = 'mercergeoinfo@gmail.com'
__maintainer__ = 'Andrew Mercer'
__version__ = '1.0.0'
__date__ = '12/02/2016'
#
### IMPORTS
import sys
import os
import copy
# import glob
# import fnmatch
# import gc
# import csv
# import pickle
# import time
# from datetime import datetime, timedelta
# import subprocess
# import urllib2
# import os.path
# import zipfile
# import StringIO
#
## IMAGE MANIPULATION
# import cv2
# import PIL
# from PIL import Image
# from PIL.GifImagePlugin import getheader, getdata
#
## SCIENTIFIC
# import math
from math import *
import numpy as np
# import numpy.ma as ma
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
# import osr
# from geostatsmodels import utilities, kriging, variograms, model, geoplot
# from qgis.core import QgsProject
# from PyQt4.QtCore import QFileInfo
# #To run in QGIS enter following in QGIS console:
# #execfile(u'/path/to/script.py'.encode('utf-8'))
#
# #To create python script using QGIS tools externally
# from qgis.core import *
# # supply path to where is your qgis installed
# QgsApplication.setPrefixPath("/Applications/QGIS.app/Contents/Resources/python/", True)
# # load providers
# QgsApplication.initQgis()
# # At end of script
# QgsApplication.exitQgis()
#
## PLOT TOOLS
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as dates
# from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
# import pylab
# pylab.ion()
# pylab.ioff()
#
## OWN
from standard import *
# from spatialfiles import *
# from kriging import *
#
### CLASSES
#
### FUNCTIONS

def readData(year, mainPath, db):
	db[year] = {}
	yearPath = os.path.join(mainPath, str(year))
	# Find highest numbered subdirectory in each year (last iteration of snow back calculation)
	for folder in range(20,0,-1):
		if  os.path.exists(os.path.join(yearPath, str(folder))):
			folderPath = os.path.join(yearPath, str(folder))
			break
		else:
			continue

	# Open report file
	file = os.path.join(folderPath, 'Report.txt')
	inFile = open(file,'rb')

	# Read data line by line into dictionary
	for row in inFile:
		if row != "\n":
			# Split read line of data into list
			dataIn = row.strip().split(',')
			key = dataIn[0]
			value = float(dataIn[1])
			if  'RN' in key:
				key = 'Score'
			if key in db[year].keys():
				continue
			else:
				db[year][key] = value
	return db

def reHash(db, year):
	db2 = {}
	for param in db[year].keys():
		db2[param] = {}
		for year in db.keys():
			db2[param][year] = db[year][param]
	return db2

def vectorBase(db2, years):
	db3 = {}
	for param in db2.keys():
		db3[param] = []
		for year in years:
			db3[param].append(db2[param][year])
	return db3

def printRes(txt, headers, years, span, m, db3):
	print "\n{} years".format(txt)
	print "\n\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(*headers)
	for i in span:
		print "{}\t{:2.3f}\t{:2.4f}\t{:2.4f}\t{:2.4f}\t{:2.4f}\t{:2.0f}\t{:2.4f}\t{:2.2f}\t{:2.0f}".format(years[i],*m[:,i])

	print "\n"
	for param in db3.keys():
		x = []
		for i in span:
			x.append(db3[param][i])
		xmean = np.mean(x)
		print "{}={:2.4f}".format(param, xmean)
	return 0

def plotetc(x,y,stat,season):
	cc_all = np.corrcoef(x, y['All'])[0][1]
	cc_opt = np.corrcoef(x, y['Optimal'])[0][1]
	cc_b1 = np.corrcoef(x, y['b1'])[0][1]
	cc_b2 = np.corrcoef(x, y['b2'])[0][1]
	print "Correlation coefficients for scores with {0} NAO during {1}".format(stat, season)
	print "Optimal\tb1\tb2\tAll"
	print "{0:.3f}\t{1:.3f}\t{2:.3f}\t{3:.3f}\n".format(cc_opt, cc_b1, cc_b2, cc_all)
	return 0

def plotNAO(optimised, all_, b1, b2):
	pltnao_path = '/Users/andrew/Google Drive/Work/MeltModel/InData/nao_2004to2014.csv'
	pltnao_data = import2vector(pltnao_path, dateString = '%Y/%m')

	averageNAO = np.zeros_like(pltnao_data['NAO'])
	averageNAO[averageNAO == 0] = np.nan
	minNAO = np.zeros_like(pltnao_data['NAO'])
	minNAO[minNAO == 0] = np.nan
	maxNAO = np.zeros_like(pltnao_data['NAO'])
	maxNAO[maxNAO == 0] = np.nan
	ispan = [-5, -4, -3, -2, -1, 0]
	i = 5
	while i < len(pltnao_data['NAO']):
		part = []
		for j in ispan:
			part.append(pltnao_data['NAO'][i + j])
		averageNAO[i] =  np.mean(part)
		minNAO[i] = np.min(part)
		maxNAO[i] = np.max(part)
		i = i+1

	matplotlib.rcParams['axes.grid'] = True
	matplotlib.rcParams['legend.fancybox'] = True
	matplotlib.rcParams['figure.figsize'] = 18, 9
	matplotlib.rcParams['savefig.dpi'] = 300
	# Set figure name and number for pdf ploting
	pdfName = 'NAO.pdf'
	pp1 = PdfPages(os.path.join('/Users/andrew/Google Drive/Work/MeltModel/Output/',pdfName))
	fig1 = plt.figure(1)
	ax1 = fig1.add_subplot(111)
	ax1.plot(pltnao_data['DATE'], pltnao_data['NAO'], '-k', label='NAO')
	#ax1.plot(pltnao_data['DATE'], minNAO, '-', color='0.75', label='6 month min NAO')
	ax1.fill_between(pltnao_data['DATE'], minNAO, maxNAO, color='0.85', facecolor='0.85', label='6 month min-max')
	#ax1.plot(pltnao_data['DATE'], maxNAO, '-', color='0.75', label='6 month max NAO')
	ax1.plot(pltnao_data['DATE'], averageNAO, '-', color='0.95', label='NAO 6 month mean')
	ax1.set_ylabel('NAO')
	ax1.grid(b=True, which='major', color='0.75', linestyle='-')
	plt.legend(loc='lower left')
	#ax1.grid(False)
	#
	pltyear = []
	for yr in years:
		yr = str(yr) + '/06'
		pltyear.append(datetime.strptime(str(yr), "%Y/%m"))

	ax2 = ax1.twinx()
	ax2.plot(pltyear, optimised, '-ok', label='Optimised')
	ax2.plot(pltyear, b1, '-or', label='b1')
	ax2.plot(pltyear, b2, '-og', label='b2')
	ax2.plot(pltyear, all_, '-ob', label='All years')
	ax2.set_ylabel('Score')
	ax2.grid(False)
	fig1.autofmt_xdate()
	plt.xlabel('Year')
	ax1.set_xlabel('Year')
	plt.legend(loc='upper left')
	pp1.savefig()
	pp1.close()
	plt.show()
	return averageNAO, minNAO, maxNAO
#
### MAIN
#
# def main():
# Set path to containing directory
mainPath = '/Users/andrew/Google Drive/Work/MeltModel/Output/Param_notELA'
allPath = '/Users/andrew/Google Drive/Work/MeltModel/Output/Years_All'
b1Path = '/Users/andrew/Google Drive/Work/MeltModel/Output/Years_b1'
b2Path = '/Users/andrew/Google Drive/Work/MeltModel/Output/Years_b2'
# Create list of years, subdirectories
years = range(2005, 2014)

## Raw, optimised data
# Create dictionary for data storage
db = {}

# Read data year by year
for year in years:
	db = readData(year, mainPath, db)

# Rehash db to order by parameter first
db2 = reHash(db, year)

# This time make vectors cor calculating correlation
db3 = vectorBase(db2, years)

# Correlation of scores with parameter values
m = np.array([db3['Score'], db3['ddfSnow'], db3['ddfFirn'], db3['ddfSi'], db3['ddfIce'], db3['elevLapse'], db3['lapse'], db3['sfe'], db3['ELA']])
corArray = np.corrcoef(m)

print "\nCorreletion Coefficient matrix for optimised parameters and score at first available summer survey\n"
headers = ['Score', 'ddfSnow', 'ddfFirn',  'ddfSi', 'ddfIce' , 'eLps',  'lapse', 'sfe', 'ELA']
print "\n\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(*headers)
for i in range(np.shape(corArray)[0]):
	print "{}\t{:2.3f}\t{:2.3f}\t{:2.3f}\t{:2.3f}\t{:2.3f}\t{:2.3f}\t{:2.3f}\t{:2.3f}\t{:2.3f}".format(headers[i],*corArray[i,])

# Calculate average values for entire period (ex 2012), b1 and b2
entire = range(len(years))
yall = copy.copy(entire)
del yall[-2]
b1 = yall[0:4]
b2 = yall[4:]

print "\nOptimised parameter values and scores at first summer survey\n"
headers = ['Score', 'ddfSnow', 'ddfFirn',  'ddfSi', 'ddfIce' , 'eLps',  'lapse', 'sfe', 'ELA']
print "\n\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(*headers)
for i in entire:
	print "{}\t{:2.3f}\t{:2.4f}\t{:2.4f}\t{:2.4f}\t{:2.4f}\t{:2.0f}\t{:2.4f}\t{:2.2f}\t{:2.0f}".format(years[i],*m[:,i])

print "\n", "*"*60

printRes('All', headers, years, yall, m, db3)
printRes('b1', headers, years, b1, m, db3)
printRes('b2', headers, years, b2, m, db3)

## Data from period average values (All, b1, b2)

print "\n", "*"*60
print "\nModelled using parameter mean of all years\n"
# All
# Create dictionary for data storage
dbA = {}

# Read data year by year
for year in years:
	dbA = readData(year, allPath, dbA)

# Rehash db to order by parameter first
dbA2 = reHash(dbA, year)

# Vectors
dbA3 = vectorBase(dbA2, years)

mA = np.array([dbA3['Score'], dbA3['ddfSnow'], dbA3['ddfFirn'], dbA3['ddfSi'], dbA3['ddfIce'], dbA3['elevLapse'], dbA3['lapse'], dbA3['sfe'], dbA3['ELA']])

printRes('All', headers, years, entire, mA, dbA3)

print "\n", "*"*60
print "\nModelled using parameter mean of b1 years\n"
# All
# Create dictionary for data storage
dbb1 = {}

# Read data year by year
for year in years:
	dbb1 = readData(year, b1Path, dbb1)

# Rehash db to order by parameter first
dbb12 = reHash(dbb1, year)

# Vectors
dbb13 = vectorBase(dbb12, years)

mb1 = np.array([dbb13['Score'], dbb13['ddfSnow'], dbb13['ddfFirn'], dbb13['ddfSi'], dbb13['ddfIce'], dbb13['elevLapse'], dbb13['lapse'], dbb13['sfe'], dbb13['ELA']])

printRes('b1', headers, years, entire, mb1, dbb13)

print "\n", "*"*60
print "\nModelled using parameter mean of b2 years\n"
# All
# Create dictionary for data storage
dbb2 = {}

# Read data year by year
for year in years:
	dbb2 = readData(year, b2Path, dbb2)

# Rehash db to order by parameter first
dbb22 = reHash(dbb2, year)

# Vectors
dbb23 = vectorBase(dbb22, years)

mb2 = np.array([dbb23['Score'], dbb23['ddfSnow'], dbb23['ddfFirn'], dbb23['ddfSi'], dbb23['ddfIce'], dbb23['elevLapse'], dbb23['lapse'], dbb23['sfe'], dbb23['ELA']])

printRes('b2', headers, years, entire, mb2, dbb23)

# Import NAO data
nao_path = '/Users/andrew/Google Drive/Work/MeltModel/InData/nao_monthly.csv'
nao_data = import2vector(nao_path)

#years = [2005, 2006, 2007, 2008, 2009, 2010, 2011, 2013]
accum_season = [-3, -2, -1, 0, 1, 2, 3]
ablat_season = [4, 5, 6, 7, 8]
stats = ['Maximum', 'Minimum', 'Mean']
naodate = []
naoval = []
for stat in stats:
	acc_nao = []
	abl_nao = []
	all_nao = []
	nao = {}
	for year in years:
		accum = []
		ablat = []
		location_year = np.where(nao_data['YEAR'] == year)
		location_jan = location_year[0][0]
		for i in accum_season:
			accum.append(nao_data['NAO'][location_jan + i])
		for j in ablat_season:
			ablat.append(nao_data['NAO'][location_jan + j])
		bothparts = accum + ablat
		nao[year] = {'Acc':accum, 'Abl':ablat, 'All':bothparts}
		nao[year]['AccStatistics'] = quickStats(accum)
		nao[year]['AblStatistics'] = quickStats(ablat)
		nao[year]['AllStatistics'] = quickStats(bothparts)
		acc_nao.append(nao[year]['AccStatistics'][stat])
		abl_nao.append(nao[year]['AblStatistics'][stat])
		all_nao.append(nao[year]['AllStatistics'][stat])
	#plotetc(acc_nao, scores, stat, 'Accumulation')
	#plotetc(abl_nao, scores, stat, 'Ablation')
	#plotetc(all_nao, scores, stat, 'Entire')

averageNAO, minNAO, maxNAO = plotNAO(db3['Score'], dbA3['Score'], dbb13['Score'], dbb23['Score'])

mC = np.array([acc_nao, abl_nao, all_nao, db3['Score'], dbA3['Score'], dbb13['Score'], dbb23['Score']])
NAOcorArray = np.corrcoef(mC)

print "\nCorreletion Coefficient matrix for scores with NAO\n"
NAOheaders = ['NAO_acc', 'NAO_abl', 'NAO_all', 'Optm', 'All', 'b1', 'b2']
print "\n\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(*NAOheaders)
for i in range(np.shape(NAOcorArray)[0]):
	print "{}\t{:2.3f}\t{:2.3f}\t{:2.3f}\t{:2.3f}\t{:2.3f}\t{:2.3f}\t{:2.3f}".format(NAOheaders[i],*NAOcorArray[i,])




# if __name__ == "__main__":
#	 main() # Calls first function, named "main" which contains main body of programme.

