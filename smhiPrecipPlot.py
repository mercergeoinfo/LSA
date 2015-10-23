#!/Users/andrew/anaconda/bin/python
# -*- coding: utf-8 -*-
"""Plot SMHI Data"""
__author__ = 'Andrew Mercer'
__contact__ = 'mercergeoinfo@gmail.com'
__maintainer__ = 'Andrew Mercer'
__version__ = '1.0'
__date__ = '13/02/2015'
# Plot SMHI data
#
from datetime import datetime
#import sys
import os
import csv
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.dates as mdates
#import AM_Func as AM

def plotAllone(dataDict):
    """Plot data in the vector dictionary on one single plot"""
    matplotlib.rcParams['axes.grid'] = True
    matplotlib.rcParams['legend.fancybox'] = True
    matplotlib.rcParams['figure.figsize'] = 18, 9
    matplotlib.rcParams['savefig.dpi'] = 300
    # Set figure name and number for pdf ploting
    plotName = 'CumulativePrecipitation.pdf'
    destDir = '../Images'
    pp1 = PdfPages(os.path.join(destDir,plotName))
    #
    lnclr = ['k','r','g','b','y','c','m']
    lnsty = ['-','--','-.',':']
    #
    days   = mdates.MonthLocator()
    daysFmt = mdates.DateFormatter('%b')
    #
    ata=0
    cnt = 0
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111)#new
    sets = sorted(dataDict.keys())
    for nm in sets:
        lab = str(nm)
        ax1.plot(dataDict[nm]['jdate'], dataDict[nm]['cumulative'], color = lnclr[cnt-7*(cnt/7)], linestyle = lnsty[cnt/7], label = lab)
        matplotlib.pyplot.axes().set_position([0.05, 0.05, 0.70, 0.85])
        ax1.legend(bbox_to_anchor=(0.01, 0.99), loc=2, borderaxespad=0.)
        plt.xlabel('2012')
        plt.ylabel('mm')
        #ax1.set_xlim(xlimit)
        plt.title('Precipitation at nearby stations')
        cnt = cnt + 1
        ax1.xaxis.set_major_locator(days)
        ax1.xaxis.set_major_formatter(daysFmt)
        ax1.format_xdata = mdates.DateFormatter('%b')
        #ax1.grid(True)
        fig1.autofmt_xdate(rotation=30)
    #plt.suptitle(PlotTable,fontsize=16, fontweight='bold')
    pp1.savefig()
    pp1.close()
    plt.close()
    ata = ata + 1
    return ata

dataDict = {}
flist = ['../InData/SMHI/2012_5_177930.csv', '../InData/SMHI/2012_5_179960.csv', '../InData/SMHI/2012_5_188800.csv', '../InData/SMHI/2012_5_178820.csv', '../InData/SMHI/2012_5_189720.csv']
nlist = ['Ritsem', 'Nikkaluokta', 'Abisko', 'Bjorkudden', 'Rensjon']
print range(len(flist))
for i in range(len(flist)):
    TfileName = flist[i]
    name = nlist[i]
    dataDict[name] = {}
    print TfileName, name, dataDict.keys()
    # Read Temperature data from csv file and convert dates to julian days. Date format '%Y-%m-%d' is SMHI's
    TinFile = open(TfileName,'rb')
    dates = []
    times = []
    temps = []
    cumul = []
    cumsum = 0
    for line in csv.DictReader(TinFile, delimiter=','):
        dates.append(line['Date'].strip())
        date = datetime.strptime(line['Date'].strip(),'%Y-%m-%d')
        jdate = datetime.strftime(date,'%j')
        times.append(int(jdate))
        temps.append(float(line['mm'].strip()))
        cumsum = cumsum + float(line['mm'].strip())
        cumul.append(cumsum)
    TinFile.close()
    dataDict[name]['date'] = dates
    dataDict[name]['jdate'] = times
    dataDict[name]['precipitation'] = temps
    dataDict[name]['cumulative'] = cumul
    print "Done: ", name
ata = plotAllone(dataDict)