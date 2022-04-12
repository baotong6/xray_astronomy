#!/bin/bash
# -*- coding: utf-8 -*-
# written by Tong
# extract the arrival time of photons in xmm data
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
#import correct as correct
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.fftpack import fft,ifft
import scipy.signal as ss
import random
import pandas as pd
#path='/Volumes/pulsar/xmm_obs_16/'
# path='/Volumes/pulsar/period/'
# path="/Volumes/pulsar/J1302/"
path='/Users/baotong/xmm/'
def get_txt(obsID,mode,src):
    hdul_evt= fits.open(path+obsID+'/cal/'+mode+'_bary.fits')
    # # Be careful! Different column if used for other satellite-----
    x=hdul_evt[1].data['X']
    y=hdul_evt[1].data['Y']
    energy=hdul_evt[1].data['PI']
    time=hdul_evt[1].data['TIME']
    tstart=hdul_evt[1].header['TSTART']
    tstop=hdul_evt[1].header['TSTOP']
    #obs_ID=hdul_evt[0].data.field(11)
    # #--------------------------------------------------------------
    #-----------read region from filename-----------------
    # def read_region(regname):
    #     #physics坐标
    #     reg_file = []
    #     #with open(path+regname+'_s'+'.reg', 'r') as file_to_read:
    #     with open(path +obsID+ '/cal/'+regname + '.reg', 'r') as file_to_read:
    #         while True:
    #             lines = file_to_read.readline() # 整行读取数据
    #             reg_file.append(lines)
    #             if not lines:
    #                 break
    #                 pass
    #     region=reg_file[-2][7:-2]
    #     reg_x,reg_y,reg_r=[float(i) for i in region.split(',')]
    #     return [reg_x,reg_y,reg_r]
    # reg=read_region(src)

    # #-----------define region yourself-----------------
    # #----------(x,y,r) in physical coordinate----------
    reg=[26389.193,27940.813,800.00028]
    # #--------------------------------------------------


    def where_region(x,y,reg):
        #输入所有光子x,y，输出region之中光子的index
        r=np.array((x-reg[0],y-reg[1]))
        len_r=np.sqrt(r[0]**2+r[1]**2)
        temp=len_r-reg[2]
        return np.where(temp<=0)

    src_index=where_region(x,y,reg)
    src_x=x[src_index]
    src_y=y[src_index]
    src_t=time[src_index]
    src_E=energy[src_index]
    #src_ID=obs_ID[src_index]
    #src_***就是该源的所有光子的信息

    def delete_photon_ID(time,energy):
        i=0
        while i < len(energy):
            if energy[i]>10000 or energy[i]<100:
                energy=np.delete(energy,i)
                time=np.delete(time,i)
                i=i-1
            i=i+1
        return [time,energy]

    [src_t,src_E]=delete_photon_ID(src_t,src_E)

    src_txt=np.column_stack((src_t,src_E))
    src_txt = src_txt[src_txt[:,0].argsort()]
    epoch=np.array([tstart,tstop])
    #print src_txt
    os.chdir(path+obsID)
    os.system('mkdir txt')
    #os.system('rm ./txt/*.txt')
    np.savetxt(path +obsID+ '/txt/' +src+'_'+mode+ '.txt', src_txt, fmt="%.7f  %.3f ")
    print(epoch[0])
    np.savetxt(path + obsID + '/txt/' +'epoch_'+ src + '_' + mode + '.txt', [epoch],fmt="%.2f  %.2f ")
get_txt('0201290301','pn','VZ_Sex')
# obsID=['0111310101','0099020301','0111970701','0551920101','0502640201',
#        '0553720101','0111970301','0152530101','0604060101','0111970401',
#        '0009650101','0673140101','0009650201','0144840101','0551430301',
#        '0601270401','0601270201','0601270301','0601270501','0652480201',
#        '0652290101','0300870501','0145450201','0760440101','0079940101',
#        '0203050501','0098810101','0724431101','0212080601','0744180801']
#
# mode=['mos1','mos2','pn']
# name=['HT_CAS','OY_CAR','QZ_VIR','RU_PEG','SS_AUR',
#       'V893_SCO','VW_HYI','YZ_CNC','V405_PEG','WX_HYI',
#       'AO_PSC','DW_CNC','FO_AQR','HT_CAM','J1509_6649',
#       'J1649_3307','J1719_4100','J1817_2508','J1830_1232','XY_ARI',
#       'EF_Eri','BM_CrB','FL_Cet','V379_Vir','CP_Tuc',
#       'PW Aqr','WW Hor','HU Aqr','EU Cnc','AM_Her']
#get_txt(obsID[24], mode[2], name[24])
# get_txt('0851180501',mode[2],'J1302')


t10=92870637+4000;t11=92870637+12100
t20=201279094.5068995;t21=201279094.5068995+18000
def cut_txt(obsID,mode,t0,t1):
    path_out='/Users/baotong/xmm/0201290301/'
    src='VZ_Sex'
    evt_name=path +obsID+ '/txt/' +src+'_'+mode+ '.txt'
    print(evt_name)
    epoch_name=path + obsID + '/txt/' +'epoch_'+ src + '_' + mode + '.txt'
    time=np.loadtxt(evt_name)[:,0]
    energy=np.loadtxt(evt_name)[:,1]
    epoch = np.array([t0, t1])
    def filter_time(time, energy,t0,t1):
        i = 0
        if len(time) != len(energy):
            print('error')
            return None
        else:
            while i < len(time):
                if time[i] < t0 or time[i] > t1:
                    print(time[i])
                    energy = np.delete(energy, i)
                    time = np.delete(time, i)
                else:
                    i += 1
        return [time,energy]
    [time,energy]=filter_time(time,energy,t0,t1)
    np.savetxt(path_out+src+'_'+mode+'_cut.txt',np.column_stack((time,energy)),fmt="%.7f  %.3f ")
    np.savetxt(path_out+'epoch_'+ src +'_'+mode+ '_cut.txt',[epoch],fmt="%.2f  %.2f ")

cut_txt('0201290301','pn',t20,t21)