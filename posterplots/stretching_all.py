#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
from obspy.signal.cross_correlation import xcorr
import os
import glob
import cPickle as pickle
from obspy.core import UTCDateTime
from obspy.signal.filter import bandpass
import datetime
from obspy.core import read, UTCDateTime
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, DateFormatter, date2num
from stretching_fct import stretching
import time


# PERFORM STRETCHING

dataDir = '/import/sun-data/hable//STACK/'

starttime = time.time()

frequencies = glob.glob(os.path.join(dataDir, 'W1*'))
frequencies = [x.split(os.path.sep)[-1] for x in list(frequencies)]
print(frequencies)
frequencies.sort()

for frequency in frequencies:
  
    print(frequency)
    
    f = frequency.split('_')
    if f[0][1] == '0':
      freq1 = '0.' + f[0][2:]
    else:
      freq1 = f[0][1:] + '.0'
      
    if f[1][0] == '0':
      freq2 = '0.' + f[1][1:]
    else: 
      freq2 = f[1] + '.0'
      
    freqband = freq1 + '-' + freq2

    # select only directories which start with a number (excluding READMEs etc.)
    components = glob.glob(os.path.join(dataDir, frequency, '[0-9]*'))
    components.sort()

    # loop through directory structure, load sliding stack file 
    for component in components:

	print(component)
		  
	if component[-9:-7] == '11':
	    comp = 'ZZ'
	elif component[-9:-7] == '12':
	    comp = 'ZN'
	elif component[-9:-7] == '13':
	    comp = 'ZE'
	elif component[-9:-7] == '21':
	    comp = 'NZ'
	elif component[-9:-7] == '22':
	    comp = 'NN'
	elif component[-9:-7] == '23':
	    comp = 'NE'
	elif component[-9:-7] == '31':
	    comp = 'EZ'
	elif component[-9:-7] == '32':
	    comp = 'EN'
	else:
	    comp = 'EE'
	    
	filenames = glob.glob(os.path.join(dataDir, frequency, component, '*.pck'))
	filenames.sort()
	for filename in filenames:
	  
	    # discard station pairs already processed
	    pair = [filename[-20:-16], filename[-8:-4]]
	    pair.sort()
	    if filename[-20:-16] != filename[-8:-4] and pair[1] == filename[-20:-16]:
		continue
	    
	    print(filename)
	    try:
		# open pickle files
		[ref, datevec, availability, slidestack, startdate, enddate, counter, SNRvec, timevec, stacklength, stackshift, pairname, fs ] =\
		    pickle.load(open(filename, "rb"))

	    except IOError:
		print('opening %s failed' % (filename))
		
		    
	    # define epsilons to try, starttime and endtime of time window
	    epsilons = np.linspace(-0.005,0.005,101)
	    starttime = 20.0
	    endtime = 50.0

	    # determine epsilon and corresponding correlation coefficient for different frequency band
	    epsCC = np.zeros((slidestack.shape[0],6))
	    for i in xrange(0,slidestack.shape[0]):
		poseps, posmaxCC, negeps, negmaxCC, botheps, bothmaxCC = stretching(ref, slidestack[i], epsilons, timevec, starttime, endtime)
		res = np.array([poseps, posmaxCC, negeps, negmaxCC, botheps, bothmaxCC])
		epsCC[i,:] = res
	    np.save('/import/sun-data/hable/epsCC/F%s/%s/arrays/%s_F%s.npy' % (frequency[1:], component[-9:], pairname, frequency[1:]), epsCC)



	    # determine data quality: counter = 0 if SNR<5000 for ZZ or SNR<2500 for other comp, counter = 1 if SNR>5000
	    #                         --> data_quality = 30.0/30.0 = 1, if counter is 1 for all daily CCFs in the stack
	    #                         --> data_quality = 0.0/30.0 = 0, if counter is 0 for all daily CCFs in the stack
	    quality = np.zeros((slidestack.shape[0],stacklength))
	    c = 0
	    for i in xrange(0,slidestack.shape[0]):
		count = counter[c:(c+stacklength)]
		quality[i,:] = count
		c = c + stackshift

	    ssum = []
	    for  i2 in xrange(0,slidestack.shape[0]):
		s = np.sum(quality[i2])
		ssum.append(s)
	    data_quality = (np.divide(ssum,float(stacklength)))
	    
	    # get lists of epsilon and CC values
	    poseps = []
	    posCC = []
	    negeps = []
	    negCC = []
	    botheps = []
	    bothCC = []
	    
	    for i in xrange(0, epsCC.shape[0]):
		    poseps.append(epsCC[i][0])
		    posCC.append(epsCC[i][1])
		    negeps.append(epsCC[i][2])
		    negCC.append(epsCC[i][3])
		    botheps.append(epsCC[i][4])
		    bothCC.append(epsCC[i][5])
		    
	    
	    # convert epsilons to rel. velocity change in %: eps = dt/t = -dv/v
	    poseps = np.asarray(poseps)
	    negeps = np.asarray(negeps)
	    botheps = np.asarray(botheps)
	    
	    poseps = poseps * (-100.0)
	    negeps = negeps * (-100.0)
	    botheps = botheps * (-100.0)
	    

	    # add 15 days to startdate and convert UTCDateTime list to python datetime objects:
	    datelist = []
	    for i in xrange(0,len(startdate)):
		    datum = startdate[i] + 15*24*3600
		    datelist.append(datum)
	    date = [x.datetime for x in datelist]
	    
	    # define date of earthquakes
	    earthquake1 = datetime.datetime(2004,9,28)
	    earthquake2 = datetime.datetime(2003,12,22)
	    parkfield = date2num(earthquake1)
	    sansimeon = date2num(earthquake2)

	    # format the ticks
	    years = YearLocator()   # every year
	    months = MonthLocator()  # every month

	    # plot of epsilon values and corresponding CC
	    fig = plt.figure(figsize=(11.69,8.27))
	    fig.suptitle('\n Relative velocity change and correlation coefficients for %s-%s \n \
			%s component, frequency band: %s Hz' % (pairname[0:4], pairname[12:16], comp, freqband), fontsize = 15, fontweight='bold')
	    
	    fig.add_subplot(2,3,1)
	    ax1 = plt.gca()
	    im = ax1.scatter(date, poseps, c=data_quality, vmin=0.0, vmax=1.0, cmap='YlOrRd')
	    ax1.set_title('positive time window \n', fontsize = 12)
	    ax1.set_ylabel('$\Delta$v/v (%) \n', fontsize = 12)
	    fig.add_subplot(2,3,4)
	    ax4 = plt.gca()
	    ax4.scatter(date, posCC, c=data_quality, vmin=0.0, vmax=1.0, cmap='YlOrRd')
	    ax4.set_ylabel('correlation coefficient \n', fontsize = 12)

	    fig.add_subplot(2,3,2)
	    ax2 = plt.gca()
	    ax2.scatter(date,negeps, c=data_quality, vmin=0.0, vmax=1.0, cmap='YlOrRd')
	    ax2.set_title('negative time window \n', fontsize = 12)
	    fig.add_subplot(2,3,5)
	    ax5 = plt.gca()
	    ax5.scatter(date, negCC, c=data_quality, vmin=0.0, vmax=1.0, cmap='YlOrRd')

	    fig.add_subplot(2,3,3)
	    ax3 = plt.gca()
	    ax3.scatter(date, botheps, c=data_quality, vmin=0.0, vmax=1.0, cmap='YlOrRd')
	    ax3.set_title('both time windows \n', fontsize = 12)
	    fig.add_subplot(2,3,6)
	    ax6 = plt.gca()
	    ax6.scatter(date,bothCC, c=data_quality, vmin=0.0, vmax=1.0, cmap='YlOrRd')


	    for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
		    ax.xaxis.set_major_locator(years)
		    ax.xaxis.set_major_formatter(DateFormatter('%Y'))
		    ax.xaxis.set_minor_locator(months)
		    ax.grid()
		    ax.axvline(x=parkfield, color='k', ls='dashed', lw=1.4)
		    ax.axvline(x=sansimeon, color='k', ls='dashed', lw=1.4)

	    for axx in [ax1, ax2, ax3]:
		    axx.set_ylim(-0.55, 0.55)
		    axx.yaxis.set_ticks([-0.5, -0.3, 0.0, 0.3, 0.5])
	    for axxx in [ax4, ax5, ax6]:
		    axxx.set_ylim(-0.3, 1.1)
		    axxx.yaxis.set_ticks([-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])


	    fig.subplots_adjust(bottom=0.2, wspace=0.45, top=0.81, hspace=0.3) 
	    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.02])
	    cb = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
	    cb.set_label('data quality')
	    
	    plt.savefig('/import/sun-data/hable/epsCC/F%s/%s/plots/%s_F%s.pdf' % (frequency[1:], component[-9:], pairname, frequency[1:]))
	    plt.close()
	    
	    
endtime = time.time()
print("elapsed time: %s s" % (endtime - starttime))