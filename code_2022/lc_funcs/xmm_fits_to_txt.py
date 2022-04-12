#!/bin/bash
# -*- coding: utf-8 -*-
# written by Tong
# extract the arrival time of photons in xmm data
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import random
import pandas as pd

def where_region(x, y, reg):
    # 输入所有光子x,y，输出region之中光子的index
    r = np.array((x - reg[0], y - reg[1]))
    len_r = np.sqrt(r[0] ** 2 + r[1] ** 2)
    temp = len_r - reg[2]
    return np.where(temp <= 0)

def delete_photon_ID(time,energy,band):
    i=0
    while i < len(energy):
        if energy[i]>band[1] or energy[i]<band[0]:
            energy=np.delete(energy,i)
            time=np.delete(time,i)
            i=i-1
        i=i+1
    return [time,energy]

def check_gti(epoch,src_t):
    ##check for ##
    if len(epoch)==1:
        print(np.where((src_t<epoch[0][0])|(src_t>epoch[0][1])))
    for i in range(len(epoch)-1):
        gap_start=epoch[i][1]
        gap_stop=epoch[i+1][0]
        t1=gap_stop-src_t
        t2=src_t-gap_start
        print(np.where((t1>0)&(t2>0)))

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
def get_txt(path,obsID,mode,srcname,reg_cicrle,band=[200,10000]):
    ## reg_circle should be a list[reg_x,reg_y,reg_radius]
    hdul_evt= fits.open(path+obsID+'/cal/'+mode+'_filt_time_bary.fits')
    gti_file=fits.open(path+obsID+'/cal/'+f'gti_{mode}.fits')
    epoch=[]
    gti_info=gti_file[1].data
    for i in range(len(gti_info)):
        epoch.append(list(gti_info[i]))
    epoch=np.array(epoch)
    ## Be careful! Different column if used for other satellite-----
    x=hdul_evt[1].data['X']
    y=hdul_evt[1].data['Y']
    energy=hdul_evt[1].data['PI']
    time=hdul_evt[1].data['TIME']
    tstart=hdul_evt[1].header['TSTART']
    tstop=hdul_evt[1].header['TSTOP']
    src_index=where_region(x,y,reg_cicrle)
    src_x=x[src_index]
    src_y=y[src_index]
    src_t=time[src_index]
    src_E=energy[src_index]
    #src_***就是该源的所有光子的信息

    # check_gti(epoch,src_t)

    [src_t,src_E]=delete_photon_ID(src_t,src_E,band=band)
    src_txt=np.column_stack((src_t,src_E))
    src_txt = src_txt[src_txt[:,0].argsort()]
    os.chdir(path+obsID)
    os.system('mkdir txt')
    os.system('rm ./txt/*.txt')
    np.savetxt(path+obsID+ '/txt/' +srcname+'_'+mode+ '.txt', src_txt, fmt="%.7f  %.3f ")
    np.savetxt(path+obsID + '/txt/' +'epoch_'+ srcname + '_' + mode + '.txt', epoch,fmt="%.2f  %.2f ")

if __name__=="__main__":
    print('How you doing?')