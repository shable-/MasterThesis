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
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, DateFormatter, date2num, num2date
from matplotlib import ticker
import time

# PLOTPOSTER ZNE
dates = np.load('/import/walchen-data/hable/STACK/W025_1/dates.npy')


startdate = dates[0].tolist()
startdate = [num2date(x) for x in startdate]

# add 15 days to startdate and convert UTCDateTime list to python datetime objects:
date = []
for i in xrange(0,len(startdate)):
    datum = startdate[i] + datetime.timedelta(days=15)
    date.append(datum)


ar_f1 = np.load('/import/sun-data/hable/epsCC_30/F025_1/averages/arrays/average_F025_1_ZNE_west.npy')
ar_f2 = np.load('/import/sun-data/hable/epsCC_30/F05_2/averages/arrays/average_F05_2_ZNE_west.npy')
ar_f3 = np.load('/import/sun-data/hable/epsCC_30/F1_4/averages/arrays/average_F1_4_ZNE_west.npy')

botheps_f1 = []
botheps_f2 = []
botheps_f3 = []
bothCC_f1 = []
bothCC_f2 = []
bothCC_f3 = []
#bothquality_f1 = []
#bothquality_f2 = []
#bothquality_f3 = []

for i in xrange(0, ar_f1.shape[0]):
    botheps_f1.append(ar_f1[i][6])
    bothCC_f1.append(ar_f1[i][7])
    #bothquality_f1.append(ar_f1[i][8])
    
for i in xrange(0, ar_f2.shape[0]):
    botheps_f2.append(ar_f2[i][6])
    bothCC_f2.append(ar_f2[i][7])
    #bothquality_f2.append(ar_f2[i][8])
    
for i in xrange(0, ar_f3.shape[0]):
    botheps_f3.append(ar_f3[i][6])
    bothCC_f3.append(ar_f3[i][7])
    #bothquality_f3.append(ar_f3[i][8])
    
print(len(botheps_f3), len(date[0:-1]))

botheps_f1 = np.asarray(botheps_f1)
botheps_f2 = np.asarray(botheps_f2)
botheps_f3 = np.asarray(botheps_f3)

lineup_f1 = 0.025
lineup_f2 = 0.015
lineup_f3 = 0.0

botheps_f1 = botheps_f1 * (-100) + lineup_f1
botheps_f2 = botheps_f2 * (-100) + lineup_f2
botheps_f3 = botheps_f3 * (-100) + lineup_f3

# format the ticks
years = YearLocator()   # every year
months = MonthLocator()  # every month

# define date of earthquakes
parkfield = datetime.datetime(2004,9,28)
#parkfield = date2num(earthquake1)
sansimeon = datetime.datetime(2003,12,22)
#sansimeon = date2num(earthquake2)

start = 20.0
end = 50.0

omega1 = 0.625
T1 = 0.375
fac1 = np.sqrt((6*np.sqrt(0.5*np.pi)*T1)/(omega1**2*(end**3-start**3)))
print(fac1)
rms1 = np.zeros(len(bothCC_f1))
for i in xrange(0,len(bothCC_f1)):
    rms1[i] = fac1 * (np.sqrt(1-np.power(bothCC_f1[i] ,2))/(2*bothCC_f1[i]))
rms1 = rms1 * 100 / np.sqrt(66) / np.sqrt(9)
print(bothCC_f1[0])
print(rms1[0])

omega2 = 1.25
T2 = 0.75
fac2 = np.sqrt((6*np.sqrt(0.5*np.pi)*T2)/(omega2**2*(end**3-start**3)))
print(fac2)
rms2 = np.zeros(len(bothCC_f2))
for i in xrange(0,len(bothCC_f2)):
    rms2[i] = fac2 * (np.sqrt(1-np.power(bothCC_f2[i] ,2))/(2*bothCC_f2[i]))
rms2 = rms2 * 100 / np.sqrt(78) / np.sqrt(9)

omega3 = 2.5
T3 = 1.5
fac3 = np.sqrt((6*np.sqrt(0.5*np.pi)*T3)/(omega3**2*(end**3-start**3)))
print(fac3)
rms3 = np.zeros(len(bothCC_f3))
for i in xrange(0,len(bothCC_f3)):
    rms3[i] = fac3 * (np.sqrt(1-np.power(bothCC_f3[i] ,2))/(2*bothCC_f3[i]))
rms3 = rms3 * 100 / np.sqrt(78) / np.sqrt(9)

# define xlim:
startdatum = datetime.datetime(2002, 11, 1)
enddatum = datetime.datetime(2007, 3, 1)

# plot epsilons with corresponding CC
fig = plt.figure(figsize=(11.69,8.27))
fig.suptitle('\n Relative velocity change and correlation coefficients for western stations, \n' \
	     'all components (NN, NE, NZ, EN, EE, EZ, ZN, ZE, ZZ)', fontsize=15, fontweight='bold')

fig.add_subplot(2,1,1)
ax1 = plt.gca()
ax1.plot(date, botheps_f1, c='YellowGreen', lw=1.8, label='0.25-1.0 Hz', zorder=200)
ax1.fill_between(date, botheps_f1-rms1, botheps_f1+rms1, alpha=0.3, facecolor='Green', edgecolor='YellowGreen')

ax1.plot(date, botheps_f2, c='Green', lw=1.8, label='0.5-2.0 Hz', zorder=200)
ax1.fill_between(date, botheps_f2-rms2, botheps_f2+rms2, alpha=0.2, facecolor='Green', edgecolor='Green')

ax1.plot(date[0:-1], botheps_f3, c='DarkGreen', lw=1.8, label='1.0-4.0 Hz', zorder=200)
ax1.fill_between(date[0:-1], botheps_f3-rms3, botheps_f3+rms3, alpha=0.2, facecolor='DarkGreen', edgecolor='DarkGreen')
ax1.legend(loc=4, prop={'size':12})
#ax1.set_title('0.25-1.0 Hz \n', fontsize = 13)
ax1.set_ylabel('$\Delta$v/v (%) \n', fontsize = 12)
fig.add_subplot(2,1,2)
ax2 = plt.gca()
ax2.plot(date, bothCC_f1, c='YellowGreen', lw=1.8, label='0.25-1.0 Hz')
ax2.plot(date, bothCC_f2, c='Green', lw=1.8, label='0.5-2.0 Hz')
ax2.plot(date[0:-1], bothCC_f3, c='DarkGreen', lw=1.8, label='1.0-4.0 Hz')
ax2.legend(loc=4, prop={'size':12})
ax2.set_ylabel('correlation coefficient \n \n', fontsize = 12)

for axx in [ax1]:
    axx.set_ylim(-0.165, 0.165)
    axx.yaxis.set_ticks([-0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15])
    axx.set_xlim(startdatum, enddatum)
for axxx in [ax2]:
    axxx.set_ylim(-0.05, 1.05)
    axxx.yaxis.set_ticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    axxx.set_xlim(startdatum, enddatum)
    
for ax in [ax1, ax2]:
    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_major_formatter(DateFormatter('%Y'))
    ax.xaxis.set_minor_locator(months)
    ax.grid()
    ax.axvline(x=parkfield, color='k',ls='dashed', lw=1.4)
    ax.axvline(x=sansimeon, color='k',ls='dashed', lw=1.4)
    
fig.subplots_adjust(bottom=0.1, top=0.85, hspace=0.3)

plt.savefig('/import/sun-data/hable/posterplots_new/ZNEwest_in_one.pdf')
plt.close()
