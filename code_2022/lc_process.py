#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys,os
import lc_funcs as funcs

path='/Volumes/pulsar/xmm_CV/'
reg_circle=[26719.166,27810.498,400]
obsID='0784670101'
srcname='HM_CnC'
mode='pn'

###------ 1. get event list -----###
# funcs.get_txt(path,obsID=obsID,mode=mode,srcname=srcname,reg_cicrle=reg_circle,band=[200,10000])

###------ 2. load data -----###
src_evt=np.loadtxt(path+obsID+'/txt/'+srcname+'_'+mode+'.txt')
epoch_info=np.loadtxt(path+obsID+'/txt/epoch_'+srcname+'_'+mode+'.txt')
if epoch_info.ndim == 1: epoch_info = np.array([epoch_info])
time=src_evt[:,0];energy=src_evt[:,1]

###----- 3.plot phase folding light curve ------###
time = funcs.filter_energy(src_evt[:, 0], src_evt[:, 1], [200, 1000])
funcs.phase_fold(time,epoch_info,p_test= 321.54,outpath=path+obsID+'/txt/',bin=20,net_percent=0.9,shift=0.0,label=srcname,save=False,show=True)

##------ 4. Lomb-Scargle periodogram ------###
funcs.read_SAS_lc(path=path+obsID+'/cal/',mode='pn',dt=5,freq=np.arange(1/5000.,0.5/5,1/500000))