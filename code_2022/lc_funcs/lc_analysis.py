#!/bin/bash
# -*- coding: utf-8 -*-
# written by Tong
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.optimize import curve_fit
import pandas as pd
from astropy.timeseries import LombScargle
import stingray
from stingray import Lightcurve, Powerspectrum, AveragedPowerspectrum
from stingray.events import EventList
import lc_funcs as funcs

font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 16, }

def get_hist(t, len_bin,tstart=0,tstop=0):
    ###将输入的time信息，按照len_bin的长度输出为lc
    if tstart==0 and tstop==0:
        tstart=t[0]
        tstop=t[-1]
        tseg=tstop-tstart
    else:tseg=tstop-tstart

    t_test = t;dt=len_bin
    ev = EventList()
    ev.time = t_test
    lc_new = ev.to_lc(dt=dt, tstart=tstart, tseg=tseg)
    return lc_new

def plot_lc(time,flux):
    ## 画光变曲线，即光子流量随时间的变化 ##
    lc = Lightcurve(time, flux)
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.plot(lc.time, lc.counts, lw=2, color='blue')
    ax.set_xlabel("Time (s)", fontproperties=font1)
    ax.set_ylabel("Counts/s›", fontproperties=font1)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.tick_params(which='major', width=1.5, length=7)
    ax.tick_params(which='minor', width=1.5, length=4)
    plt.show()
def plot_pds(time,flux):
    ## 画功率谱图 ##
    ps = Powerspectrum(lc,norm='leahy')
    ps = ps.rebin(df=2e-5)
    fig, ax1 = plt.subplots(1, 1, figsize=(9, 6), sharex=True)
    ax1.loglog()
    ax1.step(ps.freq, ps.power, lw=2, color='blue')
    ax1.set_xlabel("Frequency (Hz)", fontproperties=font1)
    ax1.set_ylabel("Power ", fontproperties=font1)
    ax1.set_yscale('log')
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax1.tick_params(which='major', width=1.5, length=7)
    ax1.tick_params(which='minor', width=1.5, length=4)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(1.5)
    plt.show()


    ## 画average的功率谱，对于绝大多数的源不需要做这张图，除非所给的参考文献中有做average power density spectrum ###
    avg_ps = AveragedPowerspectrum(lc, 500,dt=lc.time[1]-lc.time[0],norm='leahy')
    print("Number of segments: %d" % avg_ps.m)
    fig, ax1 = plt.subplots(1, 1, figsize=(9, 6))
    ax1.loglog()
    ax1.step(avg_ps.freq, avg_ps.power, lw=2, color='blue')
    ax1.set_xlabel("Frequency (Hz)", fontproperties=font1)
    ax1.set_ylabel("Power ", fontproperties=font1)
    ax1.set_yscale('log')
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax1.tick_params(which='major', width=1.5, length=7)
    ax1.tick_params(which='minor', width=1.5, length=4)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(1.5)
    plt.show()

def plot_SAS_lc(path,mode='pn',dt=5,freq=None):
    # 1,2,3分别代表mos1,mos2,pn的light curve，也可以加起来用，记为_all;
    # 根据实际情况来决定lomb-scargle的输入
    # 要注意的是，如果在run_XMMproducts_spectra.py中没有自行指定tmin和tmax，这里就不能直接对三个detector的lc进行加减
    os.chdir(path)
    filename3=mode+'_lccorr_bin{0}.lc'.format(int(dt))
    lc3=fits.open(filename3)
    time3=lc3[1].data['TIME'];rate3=lc3[1].data['RATE']
    rate3=np.nan_to_num(rate3)
    rate3[np.where(rate3<0)]=0
    print('counts=',np.sum(rate3))
    time_all=time3;rate_all=rate3
    index_gti=np.where(rate_all>0)
    time_all=time_all[index_gti];rate_all=rate_all[index_gti]
    T_all=time_all[-1]-time_all[0]
    if not freq.any(): freq=np.arange(1/T_all,0.5/dt,1/(5*T_all))
    plot_lc(time_all,rate_all)
    funcs.get_LS(time_all,rate_all,freq)
    # plot_pds(time_all,rate_all)


if __name__ == "__main__":
    print('How you doing?')