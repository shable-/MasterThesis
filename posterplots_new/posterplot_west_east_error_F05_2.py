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
dates = np.load('/import/walchen-data/hable/STACK_30/W05_2/dates.npy')


startdate = dates[0].tolist()
startdate = [num2date(x) for x in startdate]

# add 15 days to startdate and convert UTCDateTime list to python datetime objects:
date = []
for i in xrange(0,len(startdate)):
    datum = startdate[i] + datetime.timedelta(days=15)
    date.append(datum)
    


west = np.load('/import/sun-data/hable/epsCC_30/F05_2/averages/arrays/average_F05_2_ZNE_west.npy')
east = np.load('/import/sun-data/hable/epsCC_30/F05_2/averages/arrays/average_F05_2_ZNE_east.npy')

print(west.shape[0])
print(east.shape[0])
print(west.shape[1])
print(east.shape[1])
print(len(date))

westeps = []
easteps = []

westCC = []
eastCC = []

#bothquality_f1 = []
#bothquality_f2 = []
#bothquality_f3 = []

for i in xrange(0, west.shape[0]):
    westeps.append(west[i][6])
    westCC.append(west[i][7])
    #bothquality_f1.append(ar_f1[i][8])
    
for i in xrange(0, east.shape[0]):
    easteps.append(east[i][6])
    eastCC.append(east[i][7])
    #bothquality_f2.append(ar_f2[i][8])
    

westeps = np.asarray(westeps)
easteps = np.asarray(easteps)


lineup_west = 0.022
lineup_east = 0.0

westeps = westeps * (-100) + lineup_west
easteps = easteps * (-100) + lineup_east

print(len(westeps))
print(len(easteps))


start = 20.
end = 50.


omega1 = 1.25
T1 = 0.75
fac1 = np.sqrt((6*np.sqrt(0.5*np.pi)*T1)/(omega1**2*(end**3-start**3)))


rms1 = np.zeros(len(eastCC))
for i in xrange(0,len(eastCC)):
    rms1[i] = fac1 * (np.sqrt(1-np.power(eastCC[i] ,2))/(2*eastCC[i]))
rms_east = rms1 * 100 / np.sqrt(66) / np.sqrt(9)


rms2 = np.zeros(len(westCC))
for i in xrange(0,len(westCC)):
    rms2[i] = fac1 * (np.sqrt(1-np.power(westCC[i] ,2))/(2*westCC[i]))
rms_west = rms2 * 100 / np.sqrt(78) / np.sqrt(9)


# format the ticks
years = YearLocator()   # every year
months = MonthLocator()  # every month

# define date of earthquakes
parkfield = datetime.datetime(2004,9,28)
#parkfield = date2num(earthquake1)
sansimeon = datetime.datetime(2003,12,22)
#sansimeon = date2num(earthquake2)

# define xlim:
startdatum = datetime.datetime(2002, 11, 1)
enddatum = datetime.datetime(2007, 3, 1)

# plot epsilons with corresponding CC
fig = plt.figure(figsize=(11.69,8.27))
fig.suptitle('\n Relative velocity change and correlation coefficients for east- and westside of the SAF, \n' \
	     'all components (NN, NE, NZ, EN, EE, EZ, ZN, ZE, ZZ), frequency band: 0.5-2.0 Hz', fontsize=15, fontweight='bold')

fig.add_subplot(2,1,1)
ax1 = plt.gca()
ax1.fill_between(date, easteps-rms_east, easteps+rms_east, facecolor='SteelBlue', edgecolor='SteelBlue', zorder=1, alpha=0.2)
ax1.fill_between(date, westeps-rms_west, westeps+rms_west, facecolor='MidnightBlue', edgecolor='MidnightBlue', zorder=1, alpha=0.2)
ax1.plot(date, easteps, c='SteelBlue', lw=1.8, label='eastside')
ax1.plot(date, westeps, c='MidnightBlue', lw=1.8, label='westside')
ax1.legend(loc=4, prop={'size':12})
#ax1.set_title('0.25-1.0 Hz \n', fontsize = 13)
ax1.set_ylabel('$\Delta$v/v (%) \n', fontsize = 12)
fig.add_subplot(2,1,2)
ax2 = plt.gca()
ax2.plot(date, eastCC, c='SteelBlue', lw=1.8, label='eastside')
ax2.plot(date, westCC, c='MidnightBlue', lw=1.8, label='westside')
ax2.legend(loc=4, prop={'size':12})
ax2.set_ylabel('correlation coefficient \n \n', fontsize = 12)

for axx in [ax1]:
    axx.set_ylim(-0.125, 0.125)
    axx.yaxis.set_ticks([-0.1, -0.05, 0.0, 0.05, 0.1])
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

plt.savefig('/import/sun-data/hable/posterplots_new/east_west_F05_2.pdf')
plt.close()