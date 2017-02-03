#!/Users/andrew/anaconda/bin/python
# -*- coding: utf-8 -*-
#from __future__ import division
"""NAO Correlation with Temperature at Tarfala"""
__author__ = 'Andrew Mercer'
__contact__ = 'mercergeoinfo@gmail.com'
__maintainer__ = 'Andrew Mercer'
__version__ = '1.0.0'
__date__ = 'dd/mm/yyyy'
#
### IMPORTS
import sys
import os
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
# import scipy.stats as stats
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
from standard import *
# from spatialfiles import *
# from kriging import *
#
### CLASSES
#
### FUNCTIONS
#
### MAIN
#
# Import Temperature data
temperature_path = '/Users/andrew/Google Drive/Work/MeltModel/Edit/MonthlyAverageTemp.csv'
temperature_data = import2vector(temperature_path)
#Import NAO data
nao_path = '/Users/andrew/Google Drive/Work/MeltModel/Edit/nao_monthly_e.csv'
nao_data = import2vector(nao_path)
#
print "Length of temperature data: {0:.0f}, length of NAO data: {1:.0f}".format(len(temperature_data['YEAR']), len(nao_data['YEAR']))
# Due to a great deal of missing monthly averages in SMHI data (WTF!) we need to sort through it all.
years = [2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015]
months = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
dates_filtered = []
temperature_filtered = []
nao_filtered = []
for year in years:
	for month in months:
		nao_loc = np.where(np.logical_and(nao_data['YEAR'] == year, nao_data['MONTH'] == month))
		temperature_loc = np.where(np.logical_and(temperature_data['YEAR'] == year, temperature_data['MONTH'] == month))
		if (len(nao_loc[0]) == 1 and len(temperature_loc[0]) == 1):
			dates_filtered.append((year, month))
			temperature_filtered.append(temperature_data['TEMP'][temperature_loc][0])
			nao_filtered.append(nao_data['NAO'][nao_loc][0])
			# print "Year: {0:.0f} Month: {1:.0f}".format(year, month)
			# print "Temperature at {0:s}: {1:s}".format(temperature_loc, temperature_data['TEMP'][temperature_loc][0])
			# print "NAO at {0:s}: {1:s}".format(nao_loc, nao_data['NAO'][nao_loc][0])

# Calculate and report the month by month correlation
monthly_cc = np.corrcoef(temperature_filtered, nao_filtered)
print "\n\nCorrelation month by month"
print monthly_cc

# Get seasonal averages for each data set
accum_season = [10, 11,12, 1, 2, 3, 4]
ablat_season = [5, 6, 7, 8, 9]
# Create dictionary for data
data = {}
# Loop through filtered data
for i in range(len(dates_filtered)):
	date = dates_filtered[i]
	temperature = temperature_filtered[i]
	nao = nao_filtered[i]
	# If the date belongs to the following year's accumulation season shift year up one
	if date[1] > 9:
		year = date[0] + 1
	else:
		year = date[0]
	# If year not already encountered create a new post in dictionary for it
	if year not in data.keys():
		data[year] = {}
		data[year]['AccSeasonTEMP'] = []
		data[year]['AblSeasonTEMP'] = []
		data[year]['AccSeasonNAO'] = []
		data[year]['AblSeasonNAO'] = []
	if date[1] in accum_season:
		data[year]['AccSeasonTEMP'].append(temperature)
		data[year]['AccSeasonNAO'].append(nao)
	if date[1] in ablat_season:
		data[year]['AblSeasonTEMP'].append(temperature)
		data[year]['AblSeasonNAO'].append(nao)
# Get array of years in data dictionary into correct order
data_years = data.keys()
data_years.sort()
# Create new dictionary posts for seasonal averages
data['acc_temp'] = []
data['abl_temp'] = []
data['acc_nao'] = []
data['abl_nao'] = []
# Add seasonal averages to respective post
for year in data_years[:-1]:
		data['acc_temp'].append(np.mean(data[year]['AccSeasonTEMP']))
		data['abl_temp'].append(np.mean(data[year]['AblSeasonTEMP']))
		data['acc_nao'].append(np.mean(data[year]['AccSeasonNAO']))
		data['abl_nao'].append(np.mean(data[year]['AblSeasonNAO']))

accacc_cc = np.corrcoef(data['acc_temp'], data['acc_nao'])
print "\n\nCorrelation accumulation season temperature with accumulation season NAO"
print accacc_cc

ablabl_cc = np.corrcoef(data['abl_temp'], data['abl_nao'])
print "\n\nCorrelation ablation season temperature with ablation season NAO"
print ablabl_cc

ablacc_cc = np.corrcoef(data['abl_temp'], data['acc_nao'])
print "\n\nCorrelation ablation season temperature with accumulation season NAO"
print ablacc_cc

accabl_cc = np.corrcoef(data['acc_temp'], data['abl_nao'])
print "\n\nCorrelation accumulation season temperature with ablation season NAO"
print accabl_cc



