#!/bin/bash
# -*- coding: utf-8 -*-
# written by Tong
# plot the phase-folded light curve from txt file (version for xmm)
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
#import correct as correct
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.optimize import curve_fit
import pandas as pd
from astropy.timeseries import LombScargle

def filter_energy(time,energy,band):
    T=time
    E=energy
    i=0
    if len(time)!=len(energy):
        print('error')
        return None
    else:
        while i <len(E):
            if E[i]<band[0] or E[i]>band[1]:
                E=np.delete(E,i)
                T=np.delete(T,i)
            else:
                i+=1
    return T

def get_LS(time, flux,freq):
    path='/Users/baotong/xmm/0201290301/txt/'
    x = time
    y = flux
    # dy=np.sqrt(y)
    # plt.scatter(x,y)
    # plt.show()

    # LS = LombScargle(x, y, dy = 1, normalization = 'standard', fit_mean = True,
    #                  center_data = True).power(freq, method = 'cython')
    LS = LombScargle(x, y,normalization = 'standard')
    # LS = LombScargle(x, y, dy, normalization='psd')
    power = LS.power(freq)

    # print('freq_num={0}'.format(len(freq)))
    FP=LS.false_alarm_probability(power.max(),minimum_frequency = freq[0], maximum_frequency = freq[-1],method='baluev')
    FP_99 = LS.false_alarm_level(0.0027, minimum_frequency = freq[0], maximum_frequency = freq[-1],method='baluev')
    FP_90 = LS.false_alarm_level(0.05,  minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')
    FP_68 = LS.false_alarm_level(0.32, minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')
    plt.title('FP={0}'.format(FP))
    plt.semilogx()
    # plt.xlim(1000.,1500.)
    plt.plot(1/freq, power)
    print(1./freq[np.where(power==np.max(power))])
    # plt.plot([freq[0], freq[-1]], [FP_99, FP_99], '--')
    # plt.plot([freq[0], freq[-1]], [FP_90, FP_90], '--')
    # plt.plot([freq[0], freq[-1]], [FP_68, FP_68], '--')
    # plt.plot([1/(20.3*60),1/(20.3*60)],[0,0.08],'--',color='green')
    # plt.plot([1 / (22.5 * 60), 1 / (22.5 * 60)], [0, 0.08], '--', color='yellow')
    # plt.plot([1 / (18.6 * 60), 1 / (18.6 * 60)], [0, 0.08], '--', color='blue')

    plt.plot([20.3*60,20.3*60],[0,0.08],'--',color='green')
    plt.plot([22.5 * 60, 22.5 * 60], [0, 0.08], '--', color='yellow')
    plt.plot([18.6 * 60, 18.6 * 60], [0, 0.08], '--', color='blue')
    plt.show()
    res=1e5*power
    res=np.round(res,2)
    # res=res[::smooth]
    # np.savetxt(path+'LS_simres_{0}.txt'.format(trial+1),res,fmt='%10.5f')

    return [FP, 1. / freq[np.where(power == np.max(power))],np.max(power),res]

def read_SAS_lc():
    path='/Users/baotong/xmm/0201290301/cal/'
    os.chdir(path)
    mode=['MOS1','MOS2','PN']
    filename1=mode[0]+'_lccorr.lc';filename2=mode[1]+'_lccorr.lc';filename3=mode[2]+'_lccorr.lc'
    lc1=fits.open(filename1);lc2=fits.open(filename2);lc3=fits.open(filename3)
    time1=lc1[1].data['TIME'][0:1764];rate1=lc1[1].data['RATE'][0:1764]
    time2=lc2[1].data['TIME'][0:1764];rate2=lc2[1].data['RATE'][0:1764]
    time3=lc3[1].data['TIME'][0:1764];rate3=lc3[1].data['RATE'][0:1764]
    dt=10.
    # print(len(time1))
    freq=np.arange(1./20000,0.5/dt,1./(10*20000))

    time_all=time1;rate_all=rate3+rate2+rate1
    get_LS(time_all,rate_all,freq)
# read_SAS_lc()



def get_hist(t, len_bin):
    ###将输入的time信息，按照len_bin的长度输出为lc
    t_test = t-t[0]
    a = [0 for i in range(int(t_test[-1] / len_bin) + 1)]
    for i in range(len(t_test)):
        a[int(t_test[i] / len_bin)] += 1
    a = np.array(a)
    return a

def plot_LS(srcname):
    path = '/Users/baotong/xmm/0201290301/txt/'
    mode='pn';dt=10;
    time=np.loadtxt(path+srcname+'_'+mode+'_cut.txt')[:,0]
    energy=np.loadtxt(path+srcname+'_'+mode+'_cut.txt')[:,1]
    time=filter_energy(time,energy,[200,10000])
    flux=get_hist(time,dt)
    x=np.arange(time[0],time[-1],dt)
    freq=np.arange(1./17000,0.5/dt,1./(5*17000))
    # plt.scatter(x,flux)
    # plt.show()
    get_LS(x,flux,freq)
# plot_LS(srcname='VZ_Sex')

def filter_random_photon(time):
    i=0
    a=len(time)/1000.
    print(a)
    while i <len(time):
        if np.random.randint(0,a)==0:
            i += 1
            continue
        else:
            time=np.delete(time,i)
    print(len(time))
    return time

def get_T_in_mbins(epoch_file,w,m,fi):
    T=2*np.pi/w
    T_in_perbin = np.zeros(m)
    # 每个bin的总积分时间
    tbin = T/m
    # 每个bin的时间长度
    epoch_info = np.loadtxt(epoch_file)
    t_start = np.array([epoch_info[0]])
    t_end = np.array([epoch_info[1]])

    # t_start = epoch_info[:, 0]
    # t_end = epoch_info[:, 1]
    exptime=np.sum(t_end-t_start)
    #print(exptime)
    print('exptime={0}'.format(exptime))
    # plt.scatter(t_start,np.zeros(len(t_start))+1.)
    # plt.scatter(t_end, np.zeros(len(t_end)) + 1.)
    # plt.show()


    N_bin_t_start=t_start/tbin+m*fi/(2*np.pi)
    N_bin_t_end=t_end/tbin+m*fi/(2*np.pi)
    intN_bin_t_start=np.floor(N_bin_t_start)+1
    intN_bin_t_end=np.floor(N_bin_t_end)
    intN_bin_t_start=intN_bin_t_start.astype(int)
    intN_bin_t_end=intN_bin_t_end.astype(int)
    for i in range(len(N_bin_t_start)):
        if intN_bin_t_end[i]>=intN_bin_t_start[i]:
            T_in_perbin+=int((intN_bin_t_end[i]-intN_bin_t_start[i])/m)*tbin
            #print(intN_bin_t_start[i]-1)
            T_in_perbin[np.mod(intN_bin_t_start[i],m)-1]+=(intN_bin_t_start[i]-N_bin_t_start[i])*tbin
            T_in_perbin[np.mod(intN_bin_t_end[i],m)]+=(N_bin_t_end[i]-intN_bin_t_end[i])*tbin
            rest=np.mod(intN_bin_t_end[i]-intN_bin_t_start[i],m)
            for k in range(rest):
                T_in_perbin[int(np.mod((intN_bin_t_start[i] + k), m))] += tbin
            #print(rest)
        else:
            T_in_perbin[np.mod(intN_bin_t_start[i],m)-1]+=(N_bin_t_end[i]-N_bin_t_start[i])*tbin
    return T_in_perbin+1e-5

def phase_fold(path,data_file,p_test,bin,net_percent,shift,label):
    time=np.loadtxt(data_file)[:,0]
    energy=np.loadtxt(data_file)[:,1]
    time = filter_energy(time, energy, [200, 50000])
    plt.hist(time,bins=100,histtype = 'step')
    plt.show()

    epoch_file=path+'txt/' +'epoch_'+dataname+'.txt'
    T_in_perbin = get_T_in_mbins(epoch_file, 2 * np.pi / p_test, bin, shift*2*np.pi)

    def trans(t,p_test,shift):
        ti =t
        v = 1.0 /p_test
        turns = v * ti
        turns += shift
        # 初始相位
        for i in range(len(turns)):
            turns[i] = turns[i] - int(turns[i])
        return turns

    turns=trans(time,p_test,shift)
    loc=np.zeros(bin)
    for index in turns:
        loc[int(index*bin)] += 1

    x = np.array([(i / bin + 0.5 / bin) for i in range(bin)])

    src_bkg = 1 - net_percent
    bkg_y = len(time) * src_bkg
    bkg_y_low=bkg_y-bkg_y**0.5
    bkg_y_high = bkg_y+bkg_y**0.5

    bkg_y /= bin
    bkg_y_low/=bin
    bkg_y_high/=bin

    fig=plt.figure()
    ax1 = fig.add_subplot(111)

    bkg_x = [0, 2]
    plt.fill_between(bkg_x, bkg_y_low, bkg_y_high,facecolor = 'blue', alpha = 0.9)

    x2=np.concatenate((x,x+1))
    y2=np.concatenate((loc,loc))
    correct_gap = T_in_perbin / (sum(T_in_perbin) / len(T_in_perbin))
    print(correct_gap)
    y2 /= np.concatenate((correct_gap, correct_gap))

    #plt.figure(1,(8,8))
    plt.title("{0} P={1},cts={2}".format(label[0:-3],str(p_test),str(len(time))), fontsize = 18)
    plt.xlabel('phase')
    plt.ylabel('counts/bin')
    plt.ylim(0,(np.max(y2)+np.max(y2)**0.5)*1.05)
    plt.step(x2,y2,color='red')
    plt.errorbar(x2 - 0.5 / bin, y2, yerr = y2 ** 0.5, fmt = '.', capsize = 1, elinewidth = 1, ecolor = 'red')

    ax2 = ax1.twinx()
    yhigh=(np.max(y2)+np.max(y2)**0.5)*1.05/np.mean(y2)
    ax2.set_ylabel('Normalized flux')
    ax2.plot([0,2],[1.0,1.0],'--',color='green')
    ax2.set_ylim([0,yhigh])
    # plt.figure(2)
    # plt.hist(time,bins=100,histtype = 'step')
    # plt.show()
    #print('exptime={0}'.format(time[-1]-time[0]))

    plt.savefig(path + 'pfold_lc_spin_{0}.eps'.format(dataname[0:-3]))
    plt.show()

DN_name=['HT_CAS','OY_CAR','QZ_VIR','RU_PEG','SS_AUR',
      'V893_SCO','VW_HYI','YZ_CNC','V405_PEG','WX_HYI']

IP_name=['AO_PSC','DW_CNC','FO_AQR','HT_CAM','J1509_6649',
         'J1649_3307','J1719_4100','J1817_2508','J1830_1232','XY_ARI','V1025_Cen']
Polar_name=['EF_Eri','BM_CrB','FL_Cet','V379_Vir','CP_Tuc',
            'PW Aqr','WW Hor','HU Aqr','EU Cnc','AM_Her']
period_DN=[6363.1008,5453.65,5218.56,32365.44,15793.91,
        6563.0304,6417.0144,7499.52,15348.7008,6704.64]
period_IP=[14325.12,5166.1152,17457.984,5159.1168,21202.56,
           13020.48,14420.16,5514.048,19344.96,21833.0208,5076.864]
period_Polar=[4861.3824,5054.4,5228.5824,5305.9968,5342.2848,
              5652.72,6929.1936,7501.248,7525.44,11139.2928]
spin_IP=[805.2,2315.026,1254.284,514.6,809.42,
         597.920,1139.550,1660.8,1820,206.298,2146.6]

#path='/Users/baotong/xmm/0802410101/'
path='/Volumes/pulsar/WR/0109110101/'
net_p=0.985  ##get from your spectra
i=-1
# dataname=IP_name[i]+'_pn'
# period=spin_IP[i]
dataname='WR46_pn'
period=16123.831022250884*2

#period=14549.22
label=dataname
# phase_fold(path,path +'txt/'+dataname+'.txt', period, bin = 20, net_percent = net_p, shift = 0., label = label)
