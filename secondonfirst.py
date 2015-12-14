#!/Users/andrew/anaconda/bin/python
# -*- coding: utf-8 -*-
#from __future__ import division
"""Read through modelled melt results and compare year by year"""
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
# def main():
    # print command line arguments
    # for arg in sys.argv[1:]:
    # print arg
#
# 1.0 Get list of files to read
# 1.1 Names of top folder
topFolder1 = '/Users/andrew/Google Drive/Work/MeltModel/Output/ParamFrom2005to2008'
topFolder2 = '/Users/andrew/Google Drive/Work/MeltModel/Output/ParamFrom2009to2013'
# 1.2 Path to files within top folders
yearDirectory = ['2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013']
subDirectory = '0'
dataFile = 'melt.csv'
# 2.0 Stor data as two lists for each year
allData = {}
# 2.1 Loop through year by, picking data file from each top folder subdirectory
print "year&mean difference&std. dev.&Pearson's r \\\\"
for year in yearDirectory:
    # Read in data for each data set
    db1 = import2vector(os.path.join(topFolder1, year, subDirectory, dataFile), importError = 'no')
    db2 = import2vector(os.path.join(topFolder2, year, subDirectory, dataFile), importError = 'no')
    # Transfer relevant data to new database
    # Start with common information: name and position
    db = {}
    db['Stake'] = db1['Stake']
    db['Easting'] = db1['Easting']
    db['Northing'] = db1['Northing']
    db['Elevation'] = db1['Elevation']
    # Fetch modelled ablation
    keys1 = db1.keys()
    indices1 = [i for i, s in enumerate(keys1) if '_Bs_' in s]
    db['Bs1'] = db1[ keys1[ indices1[0] ] ]
    keys2 = db2.keys()
    indices2 = [i for i, s in enumerate(keys2) if '_Bs_' in s]
    db['Bs2'] = db2[ keys2[ indices2[0] ] ]
    db['Difference'] = db['Bs1'] - db['Bs2']
    db['MeanDifference'] = np.mean(db['Difference'])
    db['StDevDifference'] = np.std(db['Difference'])
    db['CorrelationCoefficient'] = np.corrcoef(db['Bs1'], db['Bs2'])
    print "{0}&{1:.2f}&{2:.2f}&{3:.3f} \\\\".format(year, db['MeanDifference'], db['StDevDifference'], db['CorrelationCoefficient'][0,1])
    # Transfer year database to super database
    allData[year] = db

# Compare year by year the difference at each stake of modelled accumulation
# Summarise differences
# if __name__ == "__main__":
#    main() # Calls first function, named "main" which contains main body of programme.

