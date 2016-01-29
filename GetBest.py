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
years = ['2005', '2006', '2007', '2008', '2009', '2010', '2011', '2013']
data= {}
info = {}
for year in years:
	file = os.path.join('../Output/Iterations/',year)
	file = os.path.join(file,'UsedParameters.csv')
	data[year] = import2vector(file)
	bestScore = data[year]['Description']['BsR2_min']
	location = np.where(data[year]['BsR2'] == data[year]['Description']['BsR2_min'])
	Snow = data[year]['ddfSnow'][location]
	Ice = data[year]['ddfIce'][location]
	Si = data[year]['ddfSi'][location]
	Firn = data[year]['ddfFirn'][location]
	Lapse = data[year]['lapse'][location]
	info[year] = {'Score':bestScore, 'Snow':Snow, 'Ice':Ice, 'Si':Si, 'Firn':Firn, 'Lapse':Lapse}
# def main():
	# print command line arguments
	# for arg in sys.argv[1:]:
	# print arg
#
# if __name__ == "__main__":
#	 main() # Calls first function, named "main" which contains main body of programme.

