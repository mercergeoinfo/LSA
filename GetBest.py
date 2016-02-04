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
def checkVal(x):
	if min(x) != max(x):
		x = np.nan
	else:
		x = min(x)
	return x
#
def printArray(args):
	print "\t".join(args)
#
### MAIN
#

dates = {2005:257-114, 2006:193-120, 2007:239-111, 2008:198-120, 2009:209-115, 2010:209-105, 2011:257-115, 2013:194-114}
years = dates.keys()
years.sort()
b1years = [2005, 2006, 2007, 2008]
b2years = [2009, 2010, 2011, 2013]
#
data= {}
info = {}
info['All_wb'] = {}
info['b1_wb'] = {}
info['b2_wb'] = {}
info['All_sb'] = {}
info['b1_sb'] = {}
info['b2_sb'] = {}

for year in years:
	file = os.path.join('../Output/Iterations/',str(year))
	file = os.path.join(file,'UsedParameters.csv')
	data[year] = import2vector(file)
	correlationHeaders = ['BsScore', 'ddfSnow', 'ddfSi', 'ddfFirn', 'ddfIce', 'lapse']
	X = np.vstack((data[year][correlationHeaders[0]], data[year][correlationHeaders[1]], data[year][correlationHeaders[2]], data[year][correlationHeaders[3]], data[year][correlationHeaders[4]], data[year][correlationHeaders[5]]))
	Xmask = np.ma.masked_where(X == np.nan, X)
	correl =  np.ma.corrcoef(Xmask)
	bestScore = data[year]['Description']['BsScore_min']
	location = np.where(data[year]['BsScore'] == data[year]['Description']['BsScore_min'])
	Snow = checkVal(data[year]['ddfSnow'][location])
	Ice = checkVal(data[year]['ddfIce'][location])
	Si = checkVal(data[year]['ddfSi'][location])
	Firn = checkVal(data[year]['ddfFirn'][location])
	Lapse = checkVal(data[year]['lapse'][location])
	info[year] = {'RawScore':bestScore, 'Score':bestScore/dates[year], 'Snow':Snow, 'Ice':Ice, 'Si':Si, 'Firn':Firn, 'Lapse':Lapse, 'Correlation':correl}
	print year
	printArray([x for x in correlationHeaders])
	for row in info[year]['Correlation']:
		printArray(["{:1.3f}".format(x) for x in row])

#
## SNOW
SnowYear = {}
for year in years:
	SnowYear[year] = {}
	file = os.path.join('../Output/Iterations/',str(year))
	file = os.path.join(file,'ddfSnow.csv')
	inFile = open(file, 'rb')
	for row in inFile:
		if row != "\n":
			dataIn = row.strip().split(',')
			SnowYear[year][dataIn[0]] = float(dataIn[1])
#
SnowAll = {}
vals = SnowYear[SnowYear.keys()[0]].keys()
vals.sort()
SnowAll['Total'] = 0
SnowAll['ScoreSum'] = 0
SnowAll['WeightedSum'] = 0
bestScore = 0
for val in vals:
	SnowAll[val] = []
	for year in years:
		SnowAll[val].append(SnowYear[year][val])
	sumVal = sum(SnowAll[val])
	SnowAll[val].append(sumVal)
	SnowAll['Total'] = SnowAll['Total'] + sumVal
for val in vals:
	score = 1-(SnowAll[val][-1] / SnowAll['Total'])
	SnowAll[val].append(score)
	if score > bestScore:
		bestScore = score
		SnowAll['ScoreBest'] = float(val)
	SnowAll[val].append(float(val) * SnowAll[val][-1])
for val in vals:
	SnowAll['ScoreSum'] = SnowAll['ScoreSum'] + SnowAll[val][-2]
	SnowAll['WeightedSum'] = SnowAll['WeightedSum'] + SnowAll[val][-1]
SnowAll['WeightedBest'] = SnowAll['WeightedSum'] / SnowAll['ScoreSum']
info['All_wb']['Snow'] = SnowAll['WeightedBest']
info['All_sb']['Snow'] = SnowAll['ScoreBest']
#
Snowb1 = {}
vals = SnowYear[SnowYear.keys()[0]].keys()
vals.sort()
Snowb1['Total'] = 0
Snowb1['ScoreSum'] = 0
Snowb1['WeightedSum'] = 0
bestScore = 0
for val in vals:
	Snowb1[val] = []
	for year in b1years:
		Snowb1[val].append(SnowYear[year][val])
	sumVal = sum(Snowb1[val])
	Snowb1[val].append(sumVal)
	Snowb1['Total'] = Snowb1['Total'] + sumVal
for val in vals:
	score = 1-(Snowb1[val][-1] / Snowb1['Total'])
	Snowb1[val].append(score)
	if score > bestScore:
		bestScore = score
		Snowb1['ScoreBest'] = float(val)
	Snowb1[val].append(float(val) * Snowb1[val][-1])
for val in vals:
	Snowb1['ScoreSum'] = Snowb1['ScoreSum'] + Snowb1[val][-2]
	Snowb1['WeightedSum'] = Snowb1['WeightedSum'] + Snowb1[val][-1]
Snowb1['WeightedBest'] = Snowb1['WeightedSum'] / Snowb1['ScoreSum']
info['b1_wb']['Snow'] = Snowb1['WeightedBest']
info['b1_sb']['Snow'] = Snowb1['ScoreBest']
#
Snowb2 = {}
vals = SnowYear[SnowYear.keys()[0]].keys()
vals.sort()
Snowb2['Total'] = 0
Snowb2['ScoreSum'] = 0
Snowb2['WeightedSum'] = 0
bestScore = 0
for val in vals:
	Snowb2[val] = []
	for year in b2years:
		Snowb2[val].append(SnowYear[year][val])
	sumVal = sum(Snowb2[val])
	Snowb2[val].append(sumVal)
	Snowb2['Total'] = Snowb2['Total'] + sumVal
for val in vals:
	score = 1-(Snowb2[val][-1] / Snowb2['Total'])
	Snowb2[val].append(score)
	if score > bestScore:
		bestScore = score
		Snowb2['ScoreBest'] = float(val)
	Snowb2[val].append(float(val) * Snowb2[val][-1])
for val in vals:
	Snowb2['ScoreSum'] = Snowb2['ScoreSum'] + Snowb2[val][-2]
	Snowb2['WeightedSum'] = Snowb2['WeightedSum'] + Snowb2[val][-1]
Snowb2['WeightedBest'] = Snowb2['WeightedSum'] / Snowb2['ScoreSum']
info['b2_wb']['Snow'] = Snowb2['WeightedBest']
info['b2_sb']['Snow'] = Snowb2['ScoreBest']
#
## Ice
IceYear = {}
for year in years:
	IceYear[year] = {}
	file = os.path.join('../Output/Iterations/',str(year))
	file = os.path.join(file,'ddfIce.csv')
	inFile = open(file, 'rb')
	for row in inFile:
		if row != "\n":
			dataIn = row.strip().split(',')
			IceYear[year][dataIn[0]] = float(dataIn[1])
#
IceAll = {}
vals = IceYear[IceYear.keys()[0]].keys()
vals.sort()
IceAll['Total'] = 0
IceAll['ScoreSum'] = 0
IceAll['WeightedSum'] = 0
bestScore = 0
for val in vals:
	IceAll[val] = []
	for year in years:
		IceAll[val].append(IceYear[year][val])
	sumVal = sum(IceAll[val])
	IceAll[val].append(sumVal)
	IceAll['Total'] = IceAll['Total'] + sumVal
for val in vals:
	score = 1-(IceAll[val][-1] / IceAll['Total'])
	IceAll[val].append(score)
	if score > bestScore:
		bestScore = score
		IceAll['ScoreBest'] = float(val)
	IceAll[val].append(float(val) * IceAll[val][-1])
for val in vals:
	IceAll['ScoreSum'] = IceAll['ScoreSum'] + IceAll[val][-2]
	IceAll['WeightedSum'] = IceAll['WeightedSum'] + IceAll[val][-1]
IceAll['WeightedBest'] = IceAll['WeightedSum'] / IceAll['ScoreSum']
info['All_wb']['Ice'] = IceAll['WeightedBest']
info['All_sb']['Ice'] = IceAll['ScoreBest']
#
Iceb1 = {}
vals = IceYear[IceYear.keys()[0]].keys()
vals.sort()
Iceb1['Total'] = 0
Iceb1['ScoreSum'] = 0
Iceb1['WeightedSum'] = 0
bestScore = 0
for val in vals:
	Iceb1[val] = []
	for year in b1years:
		Iceb1[val].append(IceYear[year][val])
	sumVal = sum(Iceb1[val])
	Iceb1[val].append(sumVal)
	Iceb1['Total'] = Iceb1['Total'] + sumVal
for val in vals:
	score = 1-(Iceb1[val][-1] / Iceb1['Total'])
	Iceb1[val].append(score)
	if score > bestScore:
		bestScore = score
		Iceb1['ScoreBest'] = float(val)
	Iceb1[val].append(float(val) * Iceb1[val][-1])
for val in vals:
	Iceb1['ScoreSum'] = Iceb1['ScoreSum'] + Iceb1[val][-2]
	Iceb1['WeightedSum'] = Iceb1['WeightedSum'] + Iceb1[val][-1]
Iceb1['WeightedBest'] = Iceb1['WeightedSum'] / Iceb1['ScoreSum']
info['b1_wb']['Ice'] = Iceb1['WeightedBest']
info['b1_sb']['Ice'] = Iceb1['ScoreBest']
#
Iceb2 = {}
vals = IceYear[IceYear.keys()[0]].keys()
vals.sort()
Iceb2['Total'] = 0
Iceb2['ScoreSum'] = 0
Iceb2['WeightedSum'] = 0
bestScore = 0
for val in vals:
	Iceb2[val] = []
	for year in b2years:
		Iceb2[val].append(IceYear[year][val])
	sumVal = sum(Iceb2[val])
	Iceb2[val].append(sumVal)
	Iceb2['Total'] = Iceb2['Total'] + sumVal
for val in vals:
	score = 1-(Iceb2[val][-1] / Iceb2['Total'])
	Iceb2[val].append(score)
	if score > bestScore:
		bestScore = score
		Iceb2['ScoreBest'] = float(val)
	Iceb2[val].append(float(val) * Iceb2[val][-1])
for val in vals:
	Iceb2['ScoreSum'] = Iceb2['ScoreSum'] + Iceb2[val][-2]
	Iceb2['WeightedSum'] = Iceb2['WeightedSum'] + Iceb2[val][-1]
Iceb2['WeightedBest'] = Iceb2['WeightedSum'] / Iceb2['ScoreSum']
info['b2_wb']['Ice'] = Iceb2['WeightedBest']
info['b2_sb']['Ice'] = Iceb2['ScoreBest']
#
## Firn
FirnYear = {}
for year in years:
	FirnYear[year] = {}
	file = os.path.join('../Output/Iterations/',str(year))
	file = os.path.join(file,'ddfFirn.csv')
	inFile = open(file, 'rb')
	for row in inFile:
		if row != "\n":
			dataIn = row.strip().split(',')
			FirnYear[year][dataIn[0]] = float(dataIn[1])
#
FirnAll = {}
vals = FirnYear[FirnYear.keys()[0]].keys()
vals.sort()
FirnAll['Total'] = 0
FirnAll['ScoreSum'] = 0
FirnAll['WeightedSum'] = 0
bestScore = 0
for val in vals:
	FirnAll[val] = []
	for year in years:
		FirnAll[val].append(FirnYear[year][val])
	sumVal = sum(FirnAll[val])
	FirnAll[val].append(sumVal)
	FirnAll['Total'] = FirnAll['Total'] + sumVal
for val in vals:
	score = 1-(FirnAll[val][-1] / FirnAll['Total'])
	FirnAll[val].append(score)
	if score > bestScore:
		bestScore = score
		FirnAll['ScoreBest'] = float(val)
	FirnAll[val].append(float(val) * FirnAll[val][-1])
for val in vals:
	FirnAll['ScoreSum'] = FirnAll['ScoreSum'] + FirnAll[val][-2]
	FirnAll['WeightedSum'] = FirnAll['WeightedSum'] + FirnAll[val][-1]
FirnAll['WeightedBest'] = FirnAll['WeightedSum'] / FirnAll['ScoreSum']
info['All_wb']['Firn'] = FirnAll['WeightedBest']
info['All_sb']['Firn'] = FirnAll['ScoreBest']
#
Firnb1 = {}
vals = FirnYear[FirnYear.keys()[0]].keys()
vals.sort()
Firnb1['Total'] = 0
Firnb1['ScoreSum'] = 0
Firnb1['WeightedSum'] = 0
bestScore = 0
for val in vals:
	Firnb1[val] = []
	for year in b1years:
		Firnb1[val].append(FirnYear[year][val])
	sumVal = sum(Firnb1[val])
	Firnb1[val].append(sumVal)
	Firnb1['Total'] = Firnb1['Total'] + sumVal
for val in vals:
	score = 1-(Firnb1[val][-1] / Firnb1['Total'])
	Firnb1[val].append(score)
	if score > bestScore:
		bestScore = score
		Firnb1['ScoreBest'] = float(val)
	Firnb1[val].append(float(val) * Firnb1[val][-1])
for val in vals:
	Firnb1['ScoreSum'] = Firnb1['ScoreSum'] + Firnb1[val][-2]
	Firnb1['WeightedSum'] = Firnb1['WeightedSum'] + Firnb1[val][-1]
Firnb1['WeightedBest'] = Firnb1['WeightedSum'] / Firnb1['ScoreSum']
info['b1_wb']['Firn'] = Firnb1['WeightedBest']
info['b1_sb']['Firn'] = Firnb1['ScoreBest']
#
Firnb2 = {}
vals = FirnYear[FirnYear.keys()[0]].keys()
vals.sort()
Firnb2['Total'] = 0
Firnb2['ScoreSum'] = 0
Firnb2['WeightedSum'] = 0
bestScore = 0
for val in vals:
	Firnb2[val] = []
	for year in b2years:
		Firnb2[val].append(FirnYear[year][val])
	sumVal = sum(Firnb2[val])
	Firnb2[val].append(sumVal)
	Firnb2['Total'] = Firnb2['Total'] + sumVal
for val in vals:
	score = 1-(Firnb2[val][-1] / Firnb2['Total'])
	Firnb2[val].append(score)
	if score > bestScore:
		bestScore = score
		Firnb2['ScoreBest'] = float(val)
	Firnb2[val].append(float(val) * Firnb2[val][-1])
for val in vals:
	Firnb2['ScoreSum'] = Firnb2['ScoreSum'] + Firnb2[val][-2]
	Firnb2['WeightedSum'] = Firnb2['WeightedSum'] + Firnb2[val][-1]
Firnb2['WeightedBest'] = Firnb2['WeightedSum'] / Firnb2['ScoreSum']
info['b2_wb']['Firn'] = Firnb2['WeightedBest']
info['b2_sb']['Firn'] = Firnb2['ScoreBest']
#
## Si
SiYear = {}
for year in years:
	SiYear[year] = {}
	file = os.path.join('../Output/Iterations/',str(year))
	file = os.path.join(file,'ddfSi.csv')
	inFile = open(file, 'rb')
	for row in inFile:
		if row != "\n":
			dataIn = row.strip().split(',')
			SiYear[year][dataIn[0]] = float(dataIn[1])
#
SiAll = {}
vals = SiYear[SiYear.keys()[0]].keys()
vals.sort()
SiAll['Total'] = 0
SiAll['ScoreSum'] = 0
SiAll['WeightedSum'] = 0
bestScore = 0
for val in vals:
	SiAll[val] = []
	for year in years:
		SiAll[val].append(SiYear[year][val])
	sumVal = sum(SiAll[val])
	SiAll[val].append(sumVal)
	SiAll['Total'] = SiAll['Total'] + sumVal
for val in vals:
	score = 1-(SiAll[val][-1] / SiAll['Total'])
	SiAll[val].append(score)
	if score > bestScore:
		bestScore = score
		SiAll['ScoreBest'] = float(val)
	SiAll[val].append(float(val) * SiAll[val][-1])
for val in vals:
	SiAll['ScoreSum'] = SiAll['ScoreSum'] + SiAll[val][-2]
	SiAll['WeightedSum'] = SiAll['WeightedSum'] + SiAll[val][-1]
SiAll['WeightedBest'] = SiAll['WeightedSum'] / SiAll['ScoreSum']
info['All_wb']['Si'] = SiAll['WeightedBest']
info['All_sb']['Si'] = SiAll['ScoreBest']
#
Sib1 = {}
vals = SiYear[SiYear.keys()[0]].keys()
vals.sort()
Sib1['Total'] = 0
Sib1['ScoreSum'] = 0
Sib1['WeightedSum'] = 0
bestScore = 0
for val in vals:
	Sib1[val] = []
	for year in b1years:
		Sib1[val].append(SiYear[year][val])
	sumVal = sum(Sib1[val])
	Sib1[val].append(sumVal)
	Sib1['Total'] = Sib1['Total'] + sumVal
for val in vals:
	score = 1-(Sib1[val][-1] / Sib1['Total'])
	Sib1[val].append(score)
	if score > bestScore:
		bestScore = score
		Sib1['ScoreBest'] = float(val)
	Sib1[val].append(float(val) * Sib1[val][-1])
for val in vals:
	Sib1['ScoreSum'] = Sib1['ScoreSum'] + Sib1[val][-2]
	Sib1['WeightedSum'] = Sib1['WeightedSum'] + Sib1[val][-1]
Sib1['WeightedBest'] = Sib1['WeightedSum'] / Sib1['ScoreSum']
info['b1_wb']['Si'] = Sib1['WeightedBest']
info['b1_sb']['Si'] = Sib1['ScoreBest']
#
Sib2 = {}
vals = SiYear[SiYear.keys()[0]].keys()
vals.sort()
Sib2['Total'] = 0
Sib2['ScoreSum'] = 0
Sib2['WeightedSum'] = 0
bestScore = 0
for val in vals:
	Sib2[val] = []
	for year in b2years:
		Sib2[val].append(SiYear[year][val])
	sumVal = sum(Sib2[val])
	Sib2[val].append(sumVal)
	Sib2['Total'] = Sib2['Total'] + sumVal
for val in vals:
	score = 1-(Sib2[val][-1] / Sib2['Total'])
	Sib2[val].append(score)
	if score > bestScore:
		bestScore = score
		Sib2['ScoreBest'] = float(val)
	Sib2[val].append(float(val) * Sib2[val][-1])
for val in vals:
	Sib2['ScoreSum'] = Sib2['ScoreSum'] + Sib2[val][-2]
	Sib2['WeightedSum'] = Sib2['WeightedSum'] + Sib2[val][-1]
Sib2['WeightedBest'] = Sib2['WeightedSum'] / Sib2['ScoreSum']
info['b2_wb']['Si'] = Sib2['WeightedBest']
info['b2_sb']['Si'] = Sib2['ScoreBest']
#
## Lapse
LapseYear = {}
for year in years:
	LapseYear[year] = {}
	file = os.path.join('../Output/Iterations/',str(year))
	file = os.path.join(file,'lapse.csv')
	inFile = open(file, 'rb')
	for row in inFile:
		if row != "\n":
			dataIn = row.strip().split(',')
			LapseYear[year][dataIn[0]] = float(dataIn[1])
#
LapseAll = {}
vals = LapseYear[LapseYear.keys()[0]].keys()
vals.sort()
LapseAll['Total'] = 0
LapseAll['ScoreSum'] = 0
LapseAll['WeightedSum'] = 0
bestScore = 0
for val in vals:
	LapseAll[val] = []
	for year in years:
		LapseAll[val].append(LapseYear[year][val])
	sumVal = sum(LapseAll[val])
	LapseAll[val].append(sumVal)
	LapseAll['Total'] = LapseAll['Total'] + sumVal
for val in vals:
	score = 1-(LapseAll[val][-1] / LapseAll['Total'])
	LapseAll[val].append(score)
	if score > bestScore:
		bestScore = score
		LapseAll['ScoreBest'] = float(val)
	LapseAll[val].append(float(val) * LapseAll[val][-1])
for val in vals:
	LapseAll['ScoreSum'] = LapseAll['ScoreSum'] + LapseAll[val][-2]
	LapseAll['WeightedSum'] = LapseAll['WeightedSum'] + LapseAll[val][-1]
LapseAll['WeightedBest'] = LapseAll['WeightedSum'] / LapseAll['ScoreSum']
info['All_wb']['Lapse'] = LapseAll['WeightedBest']
info['All_sb']['Lapse'] = LapseAll['ScoreBest']
#
Lapseb1 = {}
vals = LapseYear[LapseYear.keys()[0]].keys()
vals.sort()
Lapseb1['Total'] = 0
Lapseb1['ScoreSum'] = 0
Lapseb1['WeightedSum'] = 0
bestScore = 0
for val in vals:
	Lapseb1[val] = []
	for year in b1years:
		Lapseb1[val].append(LapseYear[year][val])
	sumVal = sum(Lapseb1[val])
	Lapseb1[val].append(sumVal)
	Lapseb1['Total'] = Lapseb1['Total'] + sumVal
for val in vals:
	score = 1-(Lapseb1[val][-1] / Lapseb1['Total'])
	Lapseb1[val].append(score)
	if score > bestScore:
		bestScore = score
		Lapseb1['ScoreBest'] = float(val)
	Lapseb1[val].append(float(val) * Lapseb1[val][-1])
for val in vals:
	Lapseb1['ScoreSum'] = Lapseb1['ScoreSum'] + Lapseb1[val][-2]
	Lapseb1['WeightedSum'] = Lapseb1['WeightedSum'] + Lapseb1[val][-1]
Lapseb1['WeightedBest'] = Lapseb1['WeightedSum'] / Lapseb1['ScoreSum']
info['b1_wb']['Lapse'] = Lapseb1['WeightedBest']
info['b1_sb']['Lapse'] = Lapseb1['ScoreBest']
#
Lapseb2 = {}
vals = LapseYear[LapseYear.keys()[0]].keys()
vals.sort()
Lapseb2['Total'] = 0
Lapseb2['ScoreSum'] = 0
Lapseb2['WeightedSum'] = 0
bestScore = 0
for val in vals:
	Lapseb2[val] = []
	for year in b2years:
		Lapseb2[val].append(LapseYear[year][val])
	sumVal = sum(Lapseb2[val])
	Lapseb2[val].append(sumVal)
	Lapseb2['Total'] = Lapseb2['Total'] + sumVal
for val in vals:
	score = 1-(Lapseb2[val][-1] / Lapseb2['Total'])
	Lapseb2[val].append(score)
	if score > bestScore:
		bestScore = score
		Lapseb2['ScoreBest'] = float(val)
	Lapseb2[val].append(float(val) * Lapseb2[val][-1])
for val in vals:
	Lapseb2['ScoreSum'] = Lapseb2['ScoreSum'] + Lapseb2[val][-2]
	Lapseb2['WeightedSum'] = Lapseb2['WeightedSum'] + Lapseb2[val][-1]
Lapseb2['WeightedBest'] = Lapseb2['WeightedSum'] / Lapseb2['ScoreSum']
info['b2_wb']['Lapse'] = Lapseb2['WeightedBest']
info['b2_sb']['Lapse'] = Lapseb2['ScoreBest']
#
parameters = ['Snow', 'Si', 'Firn', 'Ice', 'Lapse']
infokeys = info.keys()
infokeys.sort()
infok = [' '] + infokeys
printArray([str(x) for x in infok])

for p in parameters:
	row = []
	row.append(p)
	for i in infokeys:
		row.append("{:1.4f}".format(info[i][p]))
	printArray([x for x in row])

