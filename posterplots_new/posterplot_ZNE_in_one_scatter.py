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


ar_f1 = np.load('/import/sun-data/hable/epsCC/F025_1/averages2/arrays/average_F025_1_ZNE.npy')
ar_f2 = np.load('/import/sun-data/hable/epsCC/F05_2/averages2/arrays/average_F05_2_ZNE.npy')
ar_f3 = np.load('/import/sun-data/hable/epsCC/F1_4/averages2/arrays/average_F1_4_ZNE.npy')

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

#start = 20.0
#end = 50.0

#omega1 = 0.625
#T1 = 0.375
#fac1 = np.sqrt((6*np.sqrt(0.5*np.pi)*T1)/(omega1**2*(end**3-start**3)))
#print(fac1)
#rms1 = np.zeros(len(bothCC_f1))
#for i in xrange(0,len(bothCC_f1)):
    #rms1[i] = fac1 * (np.sqrt(1-np.power(bothCC_f1[i] ,2))/(2*bothCC_f1[i]))
#rms1 = rms1 * 100 / np.sqrt(66) / np.sqrt(9)
#print(bothCC_f1[0])
#print(rms1[0])

#omega2 = 1.25
#T2 = 0.75
#fac2 = np.sqrt((6*np.sqrt(0.5*np.pi)*T2)/(omega2**2*(end**3-start**3)))
#print(fac2)
#rms2 = np.zeros(len(bothCC_f2))
#for i in xrange(0,len(bothCC_f2)):
    #rms2[i] = fac2 * (np.sqrt(1-np.power(bothCC_f2[i] ,2))/(2*bothCC_f2[i]))
#rms2 = rms2 * 100 / np.sqrt(78) / np.sqrt(9)

#omega3 = 2.5
#T3 = 1.5
#fac3 = np.sqrt((6*np.sqrt(0.5*np.pi)*T3)/(omega3**2*(end**3-start**3)))
#print(fac3)
#rms3 = np.zeros(len(bothCC_f3))
#for i in xrange(0,len(bothCC_f3)):
    #rms3[i] = fac3 * (np.sqrt(1-np.power(bothCC_f3[i] ,2))/(2*bothCC_f3[i]))
#rms3 = rms3 * 100 / np.sqrt(78) / np.sqrt(9)

# define xlim:
startdatum = datetime.datetime(2002, 11, 1)
enddatum = datetime.datetime(2007, 3, 1)

# plot epsilons with corresponding CC
fig = plt.figure(figsize=(11.69,8.27))
fig.suptitle('\n Relative velocity change and correlation coefficients for all stations, \n' \
	     'all components (NN, NE, NZ, EN, EE, EZ, ZN, ZE, ZZ)', fontsize=15, fontweight='bold')

fig.add_subplot(2,1,1)
ax1 = plt.gca()
ax1.scatter(date, botheps_f1, c='Orange', label='0.25-1.0 Hz', zorder=200)

ax1.scatter(date, botheps_f2, c='OrangeRed', label='0.5-2.0 Hz', zorder=200)

ax1.scatter(date[0:-1], botheps_f3, c='FireBrick', label='1.0-4.0 Hz', zorder=200)
ax1.legend(loc=4, prop={'size':12})
#ax1.set_title('0.25-1.0 Hz \n', fontsize = 13)
ax1.set_ylabel('$\Delta$v/v (%) \n', fontsize = 12)
fig.add_subplot(2,1,2)
ax2 = plt.gca()
ax2.scatter(date, bothCC_f1, c='Orange', label='0.25-1.0 Hz')
ax2.scatter(date, bothCC_f2, c='OrangeRed', label='0.5-2.0 Hz')
ax2.scatter(date[0:-1], bothCC_f3, c='FireBrick', label='1.0-4.0 Hz')
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

plt.savefig('/import/sun-data/hable/posterplots_new/ZNE_in_one_scatter.pdf')
plt.close()
	
# 	compcounter = compcounter + 1
	
	
#     print(compcounter)
    
#     poseps_all_comps_av = []
#     posCC_all_comps_av = []
#     pos_data_quality_all_comps_av = []
#     for i1 in xrange(0, stacknumber):
# 	add1 = 0.0
# 	add2 = 0.0
# 	add3 = 0.0
# 	add4 = 0.0
# 	for i2 in xrange(0,compcounter,1):
# 	    add1 = add1 + (np.power(posCC_all_comps[i2][i1], 2) * poseps_all_comps[i2][i1])
# 	    add2 = add2 + np.power(posCC_all_comps[i2][i1], 2)
# 	    add3 = add3 + np.power(posCC_all_comps[i2][i1], 3)
# 	    add4 = add4 + (np.power(posCC_all_comps[i2][i1], 2) * pos_data_quality_all_comps[i2][i1])
#         poseps_all_comps_av.append(add1/add2)
# 	posCC_all_comps_av.append(add3/add2)
# 	pos_data_quality_all_comps_av.append(add4/add2)

	
#     negeps_all_comps_av = []
#     negCC_all_comps_av = []
#     neg_data_quality_all_comps_av = []
#     for i1 in xrange(0,stacknumber):
# 	add1 = 0.0
# 	add2 = 0.0
# 	add3 = 0.0
# 	add4 = 0.0
# 	for i2 in xrange(0,compcounter,1):
# 	    add1 = add1 + (np.power(negCC_all_comps[i2][i1], 2) * negeps_all_comps[i2][i1])
# 	    add2 = add2 + np.power(negCC_all_comps[i2][i1], 2)
# 	    add3 = add3 + np.power(negCC_all_comps[i2][i1], 3)
# 	    add4 = add4 + (np.power(negCC_all_comps[i2][i1], 2) * neg_data_quality_all_comps[i2][i1])
# 	negeps_all_comps_av.append(add1/add2)
# 	negCC_all_comps_av.append(add3/add2)
# 	neg_data_quality_all_comps_av.append(add4/add2)
	    
#     botheps_all_comps_av = []
#     bothCC_all_comps_av = []
#     both_data_quality_all_comps_av = []
#     for i1 in xrange(0,stacknumber):
# 	add1 = 0.0
# 	add2 = 0.0
# 	add3 = 0.0
# 	add4 = 0.0
# 	for i2 in xrange(0,compcounter,1):
# 	    add1 = add1 + (np.power(bothCC_all_comps[i2][i1], 2) * botheps_all_comps[i2][i1])
# 	    add2 = add2 + np.power(bothCC_all_comps[i2][i1], 2)
# 	    add3 = add3 + np.power(bothCC_all_comps[i2][i1], 3)
# 	    add4 = add4 + (np.power(bothCC_all_comps[i2][i1], 2) * both_data_quality_all_comps[i2][i1])
# 	botheps_all_comps_av.append(add1/add2)
# 	bothCC_all_comps_av.append(add3/add2)
# 	both_data_quality_all_comps_av.append(add4/add2)
	
	
#     poseps_all_comps_av = np.asarray(poseps_all_comps_av)
#     posCC_all_comps_av = np.asarray(posCC_all_comps_av)	
#     pos_data_quality_all_comps_av = np.asarray(pos_data_quality_all_comps_av)
#     negeps_all_comps_av = np.asarray(negeps_all_comps_av)
#     negCC_all_comps_av = np.asarray(negCC_all_comps_av)
#     neg_data_quality_all_comps_av = np.asarray(neg_data_quality_all_comps_av)
#     botheps_all_comps_av = np.asarray(botheps_all_comps_av)
#     bothCC_all_comps_av = np.asarray(bothCC_all_comps_av)
#     both_data_quality_all_comps_av = np.asarray(both_data_quality_all_comps_av) 

#     epsCCquality_all_comps = np.zeros((stacknumber, compnumber))
#     for i in xrange(0,stacknumber):
# 	    epsCCquality_all_comps[i, 0] = poseps_all_comps_av[i]
# 	    epsCCquality_all_comps[i, 1] = posCC_all_comps_av[i]
# 	    epsCCquality_all_comps[i, 2] = pos_data_quality_all_comps_av[i]
# 	    epsCCquality_all_comps[i, 3] = negeps_all_comps_av[i]
# 	    epsCCquality_all_comps[i, 4] = negCC_all_comps_av[i]
# 	    epsCCquality_all_comps[i, 5] = neg_data_quality_all_comps_av[i]
# 	    epsCCquality_all_comps[i, 6] = botheps_all_comps_av[i]
# 	    epsCCquality_all_comps[i, 7] = bothCC_all_comps_av[i]
# 	    epsCCquality_all_comps[i, 8] = both_data_quality_all_comps_av[i]    
	    
#     np.save('/import/sun-data/hable/epsCC/%s/averages2/arrays/average_%s_ZNE.npy' %(frequency, frequency), epsCCquality_all_comps)

#     poseps_all_comps_av = poseps_all_comps_av * (-100.0)
#     negeps_all_comps_av = negeps_all_comps_av * (-100.0)
#     botheps_all_comps_av = botheps_all_comps_av * (-100.0)


		    
#     # plot epsilons with corresponding CC
#     fig = plt.figure(figsize=(11.69,8.27))
#     fig.suptitle('\n Relative velocity change and correlation coefficients for all stations \n' \
# 		  'all components, frequency band: %s Hz' % (freqband), fontsize=15, fontweight='bold')

#     fig.add_subplot(2,3,1)
#     ax1 = plt.gca()
#     im1 = ax1.scatter(date, poseps_all_comps_av, c=pos_data_quality_all_comps_av, vmin=0.0, vmax=1.0, cmap='YlOrRd')
#     ax1.set_title('positive time window \n', fontsize = 12)
#     ax1.set_ylabel('$\Delta$v/v (%) \n', fontsize = 12)
#     fig.add_subplot(2,3,4)
#     ax4 = plt.gca()
#     ax4.scatter(date, posCC_all_comps_av, c=pos_data_quality_all_comps_av, vmin=0.0, vmax=1.0, cmap='YlOrRd')
#     ax4.set_ylabel('correlation coefficient \n \n', fontsize = 12)

#     fig.add_subplot(2,3,2)
#     ax2 = plt.gca()
#     im2 = ax2.scatter(date, negeps_all_comps_av, c=neg_data_quality_all_comps_av, vmin=0.0, vmax=1.0, cmap='YlOrRd')
#     ax2.set_title('negative time window \n', fontsize = 12)
#     fig.add_subplot(2,3,5)
#     ax5 = plt.gca()
#     ax5.scatter(date, negCC_all_comps_av, c=neg_data_quality_all_comps_av, vmin=0.0, vmax=1.0, cmap='YlOrRd')

#     fig.add_subplot(2,3,3)
#     ax3 = plt.gca()
#     im3 = ax3.scatter(date, botheps_all_comps_av, c=both_data_quality_all_comps_av, vmin=0.0, vmax=1.0, cmap='YlOrRd')
#     ax3.set_title('both time windows \n', fontsize = 12)
#     fig.add_subplot(2,3,6)
#     ax6 = plt.gca()
#     ax6.scatter(date, bothCC_all_comps_av, c=both_data_quality_all_comps_av, vmin=0.0, vmax=1.0, cmap='YlOrRd')


#     for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
# 	    ax.xaxis.set_major_locator(years)
# 	    ax.xaxis.set_major_formatter(DateFormatter('%Y'))
# 	    ax.xaxis.set_minor_locator(months)
# 	    ax.grid()
# 	    ax.axvline(x=parkfield, color='k',ls='dashed', lw=1.4)
# 	    ax.axvline(x=sansimeon, color='k',ls='dashed', lw=1.4)

#     for axx in [ax1, ax2, ax3]:
# 	    axx.set_ylim(-0.125, 0.125)
# 	    axx.yaxis.set_ticks([-0.1, -0.05, 0.0, 0.05, 0.1])
#     for axxx in [ax4, ax5, ax6]:
# 	    axxx.set_ylim(-0.05, 1.05)
# 	    axxx.yaxis.set_ticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])


#     fig.subplots_adjust(bottom=0.2, wspace=0.45, top=0.81, hspace=0.3) 
#     cbar_ax1 = fig.add_axes([0.12, 0.10, 0.20, 0.02])
#     cbar_ax2 = fig.add_axes([0.415, 0.10, 0.20, 0.02])
#     cbar_ax3 = fig.add_axes([0.70, 0.10, 0.20, 0.02])
#     cb1 = fig.colorbar(im1, cax=cbar_ax1, orientation='horizontal', ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
#     cb2 = fig.colorbar(im2, cax=cbar_ax2, orientation='horizontal', ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
#     cb3 = fig.colorbar(im3, cax=cbar_ax3, orientation='horizontal', ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

#     for cb in [cb1, cb2, cb3]:
# 	cb.set_label('data quality')
# 	axcb = plt.gca() 
	
	
#     plt.savefig('/import/sun-data/hable/epsCC/%s/averages2/plots/average_%s_ZNE.pdf' %(frequency, frequency))
#     plt.close()

    
# endtime = time.time()
# print("elapsed time: %s s" % (endtime - starttime))