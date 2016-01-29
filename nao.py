#!/Users/andrew/anaconda/bin/python
# -*- coding: utf-8 -*-
#from __future__ import division
"""Template file for Python 2.7"""
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
import matplotlib
import matplotlib.pyplot as plt
# import matplotlib.dates as dates
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
#
def plotetc(x,y,stat,season):
	cc_all = np.corrcoef(x, y['All'])[0][1]
	cc_opt = np.corrcoef(x, y['Optimal'])[0][1]
	cc_b1 = np.corrcoef(x, y['b1'])[0][1]
	cc_b2 = np.corrcoef(x, y['b2'])[0][1]
	print "Correlation coefficients for scores with {0} NAO during {1}".format(stat, season)
	print "Optimal\tb1\tb2\tAll"
	print "{0:.3f}\t{1:.3f}\t{2:.3f}\t{3:.3f}\n".format(cc_opt, cc_b1, cc_b2, cc_all)
	matplotlib.rcParams['axes.grid'] = True
	matplotlib.rcParams['legend.fancybox'] = True
	matplotlib.rcParams['figure.figsize'] = 18, 9
	matplotlib.rcParams['savefig.dpi'] = 300
	# Set figure name and number for pdf ploting
	pdfName = '{0}_{1}.pdf'.format(stat, season)
	pp1 = PdfPages(os.path.join('/Users/andrew/Google Drive/Work/MeltModel/Output/',pdfName))
	fig1 = plt.figure(1)
	ax1 = fig1.add_subplot(111)
	ax1.plot(x, y['Optimal'], 'ok', label='Optimum')
	ax1.plot(x, y['All'], 'or', label='All')
	ax1.plot(x, y['b1'], 'og', label='b1')
	ax1.plot(x, y['b2'], 'ob', label='b2')
	ax1.set_xlabel("NAO")
	ax1.set_xlim((-3,3))
	ax1.set_ylabel("Score")
	#
	#ax2 = ax1.twinx()
	#ax2.plot(x, y['AdjOptimal'], 'ok', label='Adjusted')
	#ax2.set_ylabel("Adjusted Score")
	plt.title(stat)
	plt.legend(loc='upper left')
	pp1.savefig(bbox_inches='tight')
	pp1.close()
	plt.close()
	return 0
### MAIN
#
score_path = '/Users/andrew/Google Drive/Work/MeltModel/Output/Performance/Adjusted-Absolute.csv'
# score_data = np.genfromtxt(score_path, delimiter=',')
scores = {}
infile = open(score_path, 'rb')
line = infile.next()
headers = line.strip().split(',')
for row in infile:
	score_values = []
	if row != "\n":
		datain = row.strip().split(',')
		group = datain[0]
		for k in range(len(datain) -1):
			score_values.append(float(datain[k+1]))
		scores[group] = score_values
		scores['Year'] = headers[1:]
#
nao_path = '/Users/andrew/Google Drive/Work/MeltModel/InData/nao_monthly.csv'
nao_data = import2vector(nao_path)
years = [2005, 2006, 2007, 2008, 2009, 2010, 2011, 2013]
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
	plotetc(acc_nao, scores, stat, 'Accumulation')
	plotetc(abl_nao, scores, stat, 'Ablation')
	plotetc(all_nao, scores, stat, 'Entire')



pltyear = []

Optimal = [0.018,0.044,0.025,0.019,0.028,0.019,0.030,0.020]
b1 = [0.018,0.049,0.061,0.039,0.040,0.030,0.032,0.053]
b2 = [0.039,0.072,0.032,0.021,0.033,0.026,0.060,0.022]
All = [0.025,0.059,0.044,0.027,0.031,0.023,0.045,0.035]

for yr in years:
	yr = str(yr) + '/06'
	pltyear.append(datetime.strptime(str(yr), "%Y/%m"))


# Plot NAO

pltnao_path = '/Users/andrew/Google Drive/Work/MeltModel/InData/nao_2004to2014.csv'
pltnao_data = import2vector(pltnao_path, dateString = '%Y/%m')

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
ax1.set_ylabel('NAO')
ax1.grid(b=True, which='major', color='0.75', linestyle='-')
#ax1.grid(False)
#
ax2 = ax1.twinx()
ax2.plot(pltyear, Optimal, ':ok', label='Optimal')
ax2.plot(pltyear, b1, ':or', label='b1')
ax2.plot(pltyear, b2, ':og', label='b2')
ax2.plot(pltyear, All, ':ob', label='All years')
ax2.set_ylabel('Score')
ax2.grid(False)
fig1.autofmt_xdate()
plt.xlabel('Year')
ax1.set_xlabel('Year')
plt.legend(loc='upper left')
pp1.savefig()
pp1.close()
plt.show()


# def main():
	# print command line arguments
	# for arg in sys.argv[1:]:
	# print arg
#
# if __name__ == "__main__":
#	 main() # Calls first function, named "main" which contains main body of programme.

