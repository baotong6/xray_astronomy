#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import poisson_conf_interval

font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }

def trans(t, p_test, shift):
    ti = t
    v = 1.0 / p_test
    turns = v * ti
    turns += shift
    # 初始相位
    for i in range(len(turns)):
        turns[i] = turns[i] - int(turns[i])
    return turns

def get_T_in_mbins(epoch_info,w,m,fi):
    T=2*np.pi/w
    T_in_perbin = np.zeros(m)
    # 每个bin的总积分时间
    tbin = T/m
    # 每个bin的时间长度
    if epoch_info.ndim==1:epoch_info=np.array([epoch_info])
    t_start=epoch_info[:,0];t_end = epoch_info[:, 1]

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
        else:
            T_in_perbin[np.mod(intN_bin_t_start[i],m)-1]+=(N_bin_t_end[i]-N_bin_t_start[i])*tbin
    return T_in_perbin

def phase_fold(time,epoch_info,p_test,outpath,bin=20,net_percent=0.9,shift=0.0,label='test',text=None,save=False,show=True):
    turns=trans(time,p_test,shift)
    loc=np.zeros(bin)
    for index in turns:
        loc[int(index*bin)] += 1
    AM=1-min(loc)/max(loc)
    A0=AM/(2-AM)
    print('A0={0}'.format(A0))
    x = np.array([(i / bin + 0.5 / bin) for i in range(bin)])
    src_bkg = 1 - net_percent
    bkg_y = len(time) * src_bkg/bin
    b_1sigma = poisson_conf_interval(bkg_y, interval='frequentist-confidence').T
    bkg_y_low=b_1sigma[0];bkg_y_high=b_1sigma[1]
    fig=plt.figure(1,(10,7.5))
    ax1 = fig.add_subplot(111)
    bkg_x = [0, 2]
    plt.fill_between(bkg_x, bkg_y_low, bkg_y_high,facecolor = 'green', alpha = 0.5)
    x2 = np.concatenate((x, x + 1))
    y2 = np.concatenate((loc, loc))
    T_in_perbin = get_T_in_mbins(epoch_info, 2 * np.pi / p_test, bin, shift * 2 * np.pi)

    correct_gap = T_in_perbin / (sum(T_in_perbin) / len(T_in_perbin))
    # print('correct_gap=',correct_gap)
    y2 /= np.concatenate((correct_gap, correct_gap))
    y2_err = np.array(poisson_conf_interval(y2, interval='frequentist-confidence'))
    y2_err[0] = y2 - y2_err[0]
    y2_err[1] = y2_err[1] - y2

    # plt.title("#{0} P={1:.2f},C={2}".format(label, p_test, str(len(time))), fontsize=18)
    plt.xlabel('phase', font1)
    plt.ylabel('counts/bin', font1)
    plt.tick_params(labelsize=18)
    plt.ylim(0, (np.max(y2) + np.max(y2) ** 0.5) * 1.05)
    plt.step(np.concatenate(([0], x2)), np.concatenate(([y2[0]], y2)), color='red',linewidth=1.0)
    plt.errorbar(x2 - 0.5 / bin, y2, yerr=y2_err, fmt='.', capsize=1, elinewidth=1.0, ecolor='red',linewidth=1.0)
    if text:plt.text(0.,0.95, '{0}, P={1:.2f}s'.format(text,p_test), fontsize=18,fontweight='semibold',transform=ax1.transAxes)
    plt.text(1.7,0.03*np.max(y2),'C={0}'.format(str(len(time))),fontsize=18)

    ax2 = ax1.twinx()
    yhigh = (np.max(y2) + np.max(y2) ** 0.5) * 1.05 / np.mean(y2)
    ax2.set_ylabel('Normalized flux', font1)
    ax2.plot([0, 2], [1.0, 1.0], '--', color='green')
    ax2.set_ylim([0, yhigh])
    ax2.tick_params(labelsize=18)
    if save:plt.savefig(outpath + 'pfold_lc_{0}.pdf'.format(label),bbox_inches='tight', pad_inches=0.1)
    if show:plt.show()
    else:plt.close()

