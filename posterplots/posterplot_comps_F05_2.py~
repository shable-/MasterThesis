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
dates = np.load('/import/sun-data/hable/STACK/W025_1/dates.npy')


startdate = dates[0].tolist()
startdate = [num2date(x) for x in startdate]

# add 15 days to startdate and convert UTCDateTime list to python datetime objects:
date = []
for i in xrange(0,len(startdate)):
    datum = startdate[i] + datetime.timedelta(days=15)
    date.append(datum)


NN = np.load('/import/sun-data/hable/epsCC/F025_1/averages2/arrays/average_F025_1_NN.npy')
NE = np.load('/import/sun-data/hable/epsCC/F025_1/averages2/arrays/average_F025_1_NE.npy')
NZ = np.load('/import/sun-data/hable/epsCC/F025_1/averages2/arrays/average_F025_1_NZ.npy')
EN = np.load('/import/sun-data/hable/epsCC/F025_1/averages2/arrays/average_F025_1_EN.npy')
EE = np.load('/import/sun-data/hable/epsCC/F025_1/averages2/arrays/average_F025_1_EE.npy')
EZ = np.load('/import/sun-data/hable/epsCC/F025_1/averages2/arrays/average_F025_1_EZ.npy')
ZN = np.load('/import/sun-data/hable/epsCC/F025_1/averages2/arrays/average_F025_1_ZN.npy')
ZE = np.load('/import/sun-data/hable/epsCC/F025_1/averages2/arrays/average_F025_1_ZE.npy')
ZZ = np.load('/import/sun-data/hable/epsCC/F025_1/averages2/arrays/average_F025_1_ZZ.npy')

print(NN.shape[0])
print(NN.shape[1])


NNeps = []
NEeps = []
NZeps = []
ENeps = []
EEeps = []
EZeps = []
ZNeps = []
ZEeps = []
ZZeps = []

NNCC = []
NECC = []
NZCC = []
ENCC = []
EECC = []
EZCC = []
ZNCC = []
ZECC = []
ZZCC = []

NNquality = []
NEquality = []
NZquality = []
ENquality = []
EEquality = []
EZquality = []
ZNquality = []
ZEquality = []
ZZquality = []

#bothquality_f1 = []
#bothquality_f2 = []
#bothquality_f3 = []

for i in xrange(0, 144):
    NNeps.append(NN[i][6])
    NNCC.append(NN[i][7])
    NNquality.append(NN[i][8])
    NEeps.append(NE[i][6])
    NECC.append(NE[i][7])
    NEquality.append(NE[i][8])
    NZeps.append(NZ[i][6])
    NZCC.append(NZ[i][7])
    NZquality.append(NZ[i][8])
    ENeps.append(EN[i][6])
    ENCC.append(EN[i][7])
    ENquality.append(EN[i][8])
    EEeps.append(EE[i][6])
    EECC.append(EE[i][7])
    EEquality.append(EE[i][8])
    EZeps.append(EZ[i][6])
    EZCC.append(EZ[i][7])
    EZquality.append(EZ[i][8])
    ZNeps.append(ZN[i][6])
    ZNCC.append(ZN[i][7])
    ZNquality.append(ZN[i][8])
    ZEeps.append(ZE[i][6])
    ZECC.append(ZE[i][7])
    ZEquality.append(ZE[i][8])
    ZZeps.append(ZZ[i][6])
    ZZCC.append(ZZ[i][7])
    ZZquality.append(ZZ[i][8])
  
     

NNeps = np.asarray(NNeps)
NEeps = np.asarray(NEeps)
NZeps = np.asarray(NZeps)
ENeps = np.asarray(ENeps)
EEeps = np.asarray(EEeps)
EZeps = np.asarray(EZeps)
ZNeps = np.asarray(ZNeps)
ZEeps = np.asarray(ZEeps)
ZZeps = np.asarray(ZZeps)


NNeps = NNeps * (-100)
NEeps = NEeps * (-100)
NZeps = NZeps * (-100)
ENeps = ENeps * (-100)
EEeps = EEeps * (-100)
EZeps = EZeps * (-100)
ZNeps = ZNeps * (-100)
ZEeps = ZEeps * (-100)
ZZeps = ZZeps * (-100)


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

# plot epsilons
fig = plt.figure(figsize=(11.69,8.27))
fig.suptitle('\n Relative velocity change for all stations, ' \
	     'frequency band: 0.25-1.0 Hz', fontsize=15, fontweight='bold')

fig.add_subplot(3,3,1)
ax1 = plt.gca()
im = ax1.scatter(date, NNeps, c=NNquality, vmin=0.0, vmax=1.0, cmap='YlOrRd', edgecolor='None', s=15)
ax1.set_title('NN component', fontsize = 13)
ax1.set_ylabel('$\Delta$v/v (%)', fontsize = 12)

fig.add_subplot(3,3,2)
ax2 = plt.gca()
ax2.scatter(date, NEeps, c=NEquality, vmin=0.0, vmax=1.0, cmap='YlOrRd', edgecolor='None', s=15)
ax2.set_title('NE component', fontsize = 13)

fig.add_subplot(3,3,3)
ax3 = plt.gca()
ax3.scatter(date, NZeps, c=NZquality, vmin=0.0, vmax=1.0, cmap='YlOrRd', edgecolor='None', s=15)
ax3.set_title('NZ component', fontsize = 13)

fig.add_subplot(3,3,4)
ax4 = plt.gca()
ax4.scatter(date, ENeps, c=ENquality, vmin=0.0, vmax=1.0, cmap='YlOrRd', edgecolor='None', s=15)
ax4.set_title('EN component', fontsize = 13)
ax4.set_ylabel('$\Delta$v/v (%)', fontsize = 12)

fig.add_subplot(3,3,5)
ax5 = plt.gca()
ax5.scatter(date, EEeps, c=EEquality, vmin=0.0, vmax=1.0, cmap='YlOrRd', edgecolor='None', s=15)
ax5.set_title('EE component', fontsize = 13)

fig.add_subplot(3,3,6)
ax6 = plt.gca()
ax6.scatter(date, EZeps, c=EZquality, vmin=0.0, vmax=1.0, cmap='YlOrRd', edgecolor='None', s=15)
ax6.set_title('EZ component', fontsize = 13)

fig.add_subplot(3,3,7)
ax7 = plt.gca()
ax7.scatter(date, ZNeps, c=ZNquality, vmin=0.0, vmax=1.0, cmap='YlOrRd', edgecolor='None', s=15)
ax7.set_title('ZN component', fontsize = 13)
ax7.set_ylabel('$\Delta$v/v (%)', fontsize = 12)

fig.add_subplot(3,3,8)
ax8 = plt.gca()
ax8.scatter(date, ZEeps, c=ZEquality, vmin=0.0, vmax=1.0, cmap='YlOrRd', edgecolor='None', s=15)
ax8.set_title('ZE component', fontsize = 13)

fig.add_subplot(3,3,9)
ax9 = plt.gca()
ax9.scatter(date, ZZeps, c=ZZquality, vmin=0.0, vmax=1.0, cmap='YlOrRd', edgecolor='None', s=15)
ax9.set_title('ZZ component', fontsize = 13)

for ax in [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]:
    ax.set_ylim(-0.125, 0.125)
    ax.yaxis.set_ticks([-0.1, -0.05, 0.0, 0.05, 0.1])
    ax.set_xlim(startdatum, enddatum)
    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_major_formatter(DateFormatter('%Y'))
    ax.xaxis.set_minor_locator(months)
    ax.grid()
    ax.axvline(x=parkfield, color='k',ls='dashed', lw=1.4)
    ax.axvline(x=sansimeon, color='k',ls='dashed', lw=1.4)
    for tick in ax.xaxis.get_major_ticks():
	tick.label.set_fontsize(10)
    for tick in ax.yaxis.get_major_ticks():
	tick.label.set_fontsize(10)
    
fig.subplots_adjust(left=0.08, bottom=0.05, top=0.87, wspace=0.35, hspace=0.4)
cbar_ax = fig.add_axes([0.93, 0.1, 0.01, 0.7])
cb = fig.colorbar(im, cax=cbar_ax, orientation='vertical', ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
cb.set_label('data quality')

plt.savefig('/import/sun-data/hable/posterplots/comps_eps_F025_1.pdf')
plt.close()


# plot CC
fig = plt.figure(figsize=(11.69,8.27))
fig.suptitle('\n Correlation coefficients for all stations, ' \
	     'frequency band: 0.25-1.0 Hz', fontsize=15, fontweight='bold')

fig.add_subplot(3,3,1)
ax1 = plt.gca()
im = ax1.scatter(date, NNCC, c=NNquality, vmin=0.0, vmax=1.0, cmap='YlOrRd', edgecolor='None', s=15)
ax1.set_title('NN component', fontsize = 13)
ax1.set_ylabel('correlation coefficient \n', fontsize = 12)

fig.add_subplot(3,3,2)
ax2 = plt.gca()
ax2.scatter(date, NECC, c=NEquality, vmin=0.0, vmax=1.0, cmap='YlOrRd', edgecolor='None', s=15)
ax2.set_title('NE component', fontsize = 13)

fig.add_subplot(3,3,3)
ax3 = plt.gca()
ax3.scatter(date, NZCC, c=NZquality, vmin=0.0, vmax=1.0, cmap='YlOrRd', edgecolor='None', s=15)
ax3.set_title('NZ component', fontsize = 13)

fig.add_subplot(3,3,4)
ax4 = plt.gca()
ax4.scatter(date, ENCC, c=ENquality, vmin=0.0, vmax=1.0, cmap='YlOrRd', edgecolor='None', s=15)
ax4.set_title('EN component', fontsize = 13)
ax4.set_ylabel('correlation coefficient \n', fontsize = 12)

fig.add_subplot(3,3,5)
ax5 = plt.gca()
ax5.scatter(date, EECC, c=EEquality, vmin=0.0, vmax=1.0, cmap='YlOrRd', edgecolor='None', s=15)
ax5.set_title('EE component', fontsize = 13)

fig.add_subplot(3,3,6)
ax6 = plt.gca()
ax6.scatter(date, EZCC, c=EZquality, vmin=0.0, vmax=1.0, cmap='YlOrRd', edgecolor='None', s=15)
ax6.set_title('EZ component', fontsize = 13)

fig.add_subplot(3,3,7)
ax7 = plt.gca()
ax7.scatter(date, ZNCC, c=ZNquality, vmin=0.0, vmax=1.0, cmap='YlOrRd', edgecolor='None', s=15)
ax7.set_title('ZN component', fontsize = 13)
ax7.set_ylabel('correlation coefficient \n', fontsize = 12)

fig.add_subplot(3,3,8)
ax8 = plt.gca()
ax8.scatter(date, ZECC, c=ZEquality, vmin=0.0, vmax=1.0, cmap='YlOrRd', edgecolor='None', s=15)
ax8.set_title('ZE component', fontsize = 13)

fig.add_subplot(3,3,9)
ax9 = plt.gca()
ax9.scatter(date, ZZCC, c=ZZquality, vmin=0.0, vmax=1.0, cmap='YlOrRd', edgecolor='None', s=15)
ax9.set_title('ZZ component', fontsize = 13)

for ax in [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]:
    ax.set_ylim(-0.05, 1.05)
    ax.yaxis.set_ticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlim(startdatum, enddatum)
    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_major_formatter(DateFormatter('%Y'))
    ax.xaxis.set_minor_locator(months)
    ax.grid()
    ax.axvline(x=parkfield, color='k',ls='dashed', lw=1.4)
    ax.axvline(x=sansimeon, color='k',ls='dashed', lw=1.4)
    for tick in ax.xaxis.get_major_ticks():
	tick.label.set_fontsize(10)
    for tick in ax.yaxis.get_major_ticks():
	tick.label.set_fontsize(10)
    
fig.subplots_adjust(left=0.08, bottom=0.05, top=0.87, wspace=0.35, hspace=0.4)
cbar_ax = fig.add_axes([0.93, 0.1, 0.01, 0.7])
cb = fig.colorbar(im, cax=cbar_ax, orientation='vertical', ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
cb.set_label('data quality')

plt.savefig('/import/sun-data/hable/posterplots/comps_CC_F025_1.pdf')
plt.close()
