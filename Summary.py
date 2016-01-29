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
data = import2vector('../Output/Iterations/Summary.csv')
with open('../Output/Iterations/Summary.csv') as f:
	data2 = list(csv.DictReader(f))
years = ['2005', '2006', '2007', '2008', '2009', '2010', '2011', '2013']
data3 = {}
for i in range(len(data2)):
	dataset = data2[i]
	name = dataset['Parameter']
	values = []
	for year in years:
		if dataset[year] == '':
			values.append(np.nan)
		else:
			values.append(float(dataset[year]))
	data3[name] = values
keys = data3.keys()
keys.sort()
X = np.vstack((data3[keys[0]],data3[keys[1]],data3[keys[2]],data3[keys[3]],data3[keys[4]]))
covX =  np.cov(X)
Xmask = np.ma.masked_where(X == np.nan, X)
corrX =  np.ma.corrcoef(Xmask)
header = ""
for j in keys:
	header = header + j + ','
header = header[:-1]
np.savetxt('../Output/Iterations/SummaryCorr.csv', corrX, fmt='%.2f', delimiter=',', newline='\n',  footer='', comments='', header=header)
# def main():
	# print command line arguments
	# for arg in sys.argv[1:]:
	# print arg
#
# if __name__ == "__main__":
#	 main() # Calls first function, named "main" which contains main body of programme.

