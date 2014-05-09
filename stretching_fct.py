#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from obspy.signal.cross_correlation import xcorr


# highest CC is used, no discrimination between positive, negative

# STRETCHING CODE

def stretching(signalRef, signalStr, epsilons, timevec, starttime=None, endtime=None):
	"""
	Calculates the stretching factor eps. This is the factor with which a signal (signalStr)
	must be stretched to get the highest correlation with a reference signal (signalRef).
	The factor eps is chosen from an array epsilons. The time vector of both signals is timevec.
	If starttime and endtime for a time window are provided, eps is calcutlated for the positive
	and negative time window as well for both time windows together. Starttime and endtime refer
	to the positive time window; from that a negative time window is calculated
	(e.g. starttime = 20.0, endtime = 50.0 --> -20.0 and -50.0 for the negative window).
	If no starttime and endtime are given, eps is computed for the whole data.
	"""

   
	if starttime!=None and endtime!=None: # eps for time windows

		if endtime > timevec[-1]:
			raise ValueError('Time window exceeds bound of time vector!')
		if starttime < 0.0:
			raise ValueError('Positive and negative time window are overlapping!')
		if starttime > endtime:
			raise ValueError('Starttime must be smaller than endtime!')	

		# indices of starttime and endtime of the time windows
		pos_t1 = np.abs(timevec-starttime).argmin()
		pos_t2 = np.abs(timevec-endtime).argmin()
		neg_t1 = np.abs(timevec-(-endtime)).argmin()
		neg_t2 = np.abs(timevec-(-starttime)).argmin()
		
		# taper the time windows
		pos_time = timevec[pos_t1:(pos_t2+1)]
		pos_taper_percentage = 0.1
		pos_taper = np.blackman(int(len(pos_time) * pos_taper_percentage))
		pos_taper_left, pos_taper_right = np.array_split(pos_taper, 2)
		pos_taper = np.concatenate([pos_taper_left, np.ones(len(pos_time)-len(pos_taper)), pos_taper_right])

		neg_time = timevec[neg_t1:(neg_t2+1)]
		neg_taper_percentage = 0.1
		neg_taper = np.blackman(int(len(neg_time) * neg_taper_percentage))
		neg_taper_left, neg_taper_right = np.array_split(neg_taper, 2)
		neg_taper = np.concatenate([neg_taper_left, np.ones(len(neg_time)-len(neg_taper)), neg_taper_right])

		both_time = np.concatenate([neg_time, pos_time])
		both_taper_percentage = 0.1
		both_taper = np.blackman(int(len(both_time) * both_taper_percentage))
		both_taper_left, both_taper_right = np.array_split(both_taper, 2)
		both_taper = np.concatenate([both_taper_left, np.ones(len(both_time)-len(both_taper)), both_taper_right])

		pos_signalRef = pos_taper * signalRef[pos_t1:(pos_t2+1)]
		pos_signalStr = pos_taper * signalStr[pos_t1:(pos_t2+1)]
		neg_signalRef = neg_taper * signalRef[neg_t1:(neg_t2+1)]
		neg_signalStr = neg_taper * signalStr[neg_t1:(neg_t2+1)]
		both_signalRef = both_taper * np.concatenate([signalRef[neg_t1:(neg_t2+1)], signalRef[pos_t1:(pos_t2+1)]])
		both_signalStr = both_taper * np.concatenate([signalStr[neg_t1:(neg_t2+1)], signalStr[pos_t1:(pos_t2+1)]])
		
		# calculate the correlation coefficient CC for each epsilon
		posCC = []
		negCC = []
		bothCC = []
		for i in xrange(0,len(epsilons),1):
			# positive time window
			pos_time_new = (1.0-epsilons[i])*pos_time
			pos_s = InterpolatedUnivariateSpline(pos_time_new, pos_signalStr)
			pos_stretch = pos_s(pos_time)
			pos_coeffs = xcorr(pos_stretch,pos_signalRef,0)
			posCC.append(pos_coeffs[1])
			
			# negative time window
			neg_time_new = (1.0-epsilons[i])*neg_time
			neg_s = InterpolatedUnivariateSpline(neg_time_new, neg_signalStr)
			neg_stretch = neg_s(neg_time)
			neg_coeffs = xcorr(neg_stretch,neg_signalRef,0)
			negCC.append(neg_coeffs[1])
			
			# both time windows
			both_time_new = (1.0-epsilons[i])*both_time
			both_s = InterpolatedUnivariateSpline(both_time_new, both_signalStr)
			both_stretch = both_s(both_time)
			both_coeffs = xcorr(both_stretch,both_signalRef,0)
			bothCC.append(both_coeffs[1])
			
			
		# determine the max. CC and corresponding epsilon
		posmaxCC = max(posCC)
		posindex = posCC.index(posmaxCC)
		poseps = epsilons[posindex]

		negmaxCC = max(negCC)
		negindex = negCC.index(negmaxCC)
		negeps = epsilons[negindex]	

		bothmaxCC = max(bothCC)
		bothindex = bothCC.index(bothmaxCC)
		botheps = epsilons[bothindex]

		return poseps, posmaxCC, negeps, negmaxCC, botheps, bothmaxCC
	
	elif (starttime == None and endtime != None) or (starttime != None and endtime == None):
		raise SyntaxError('Both starttime and endtime must be given!')


	else: # eps for whole data

		# taper the signal and the reference
		taper_percentage = 0.1
		taper = np.blackman(int(len(timevec) * taper_percentage))
		taper_left, taper_right = np.array_split(taper, 2)
		taper = np.concatenate([taper_left, np.ones(len(timevec)-len(taper)), taper_right])

		signalStr = signalStr * taper
		signalRef = signalRef * taper
		
		# calculate the correlation coefficient CC for each epsilon
		CC = []
		for i in xrange(0,len(epsilons),1):
			time_new = (1.0-epsilons[i])*timevec
			s = InterpolatedUnivariateSpline(time_new, signalStr)
			stretch = s(timevec)
			coeffs = xcorr(stretch,signalRef,0)
			CC.append(coeffs[1])
		
		# determine the max. CC and corresponding epsilon
		maxCC = max(CC)
		index = CC.index(maxCC)
		eps = epsilons[index]
	
		return eps, maxCC
