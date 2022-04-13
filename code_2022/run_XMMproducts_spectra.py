#!/usr/bin/env python3
#written by Ben; modified by Tong Bao
# written on 2018-6-13 to run XMM products on SgrA PWN for 5 XMM obsID
# modified on 2019-3-9 to be the pipeline for point source (dummy users version,not recommended for research)
# modified on 2020-4-7 to add the extraction of light curve
# modified on 2020-4-8 to add the particle background filtering
# point source or extended source: change flag in afrgen: extendedsource=yes/no; currently set to no for point source

import sys
import os
import string
from pathlib import Path
from astropy.io import fits

# ------obsID List---------------
path="/Volumes/pulsar/xmm_CV/"
obsID1 = "0784670101"
## 改成你的路径及观测号 ##
# -------------------------------
# ------choose obsID-------------
obsList = [obsID1]
# -------------------------------
# ------detName List-------------
det1 = "mos1"
det2 = "mos2"
det3 = "pn"
# -------------------------------
#
# ------choose det --------------
detList = [det1,det2,det3]
# -------------------------------
process=0
filt_particle_bkg=0
lc=0
spectra=1
combine_spec=0

ra=obsID=121.5956395511654
dec=15.4586056893463
##改成你所做的源的坐标

for obsID in obsList:
   os.chdir(path+obsID)
   mypath=Path("./cal")
   mypath=Path("./cal")

   if process:
      if mypath.is_dir():
         print("continued")
      else:
         os.system("mkdir cal")
         os.system("gunzip ODF/*.gz")
         #------------ set environment -----------------
         os.chdir("./cal")
         os.environ['SAS_ODF'] = path + obsID + "/ODF"
         os.system("cifbuild")
         os.environ['SAS_CCF'] = path+obsID+"/cal/ccf.cif"
         os.system("odfingest")
         os.system("rm ../ODF/._*.ASC")
         os.system("rm SAS.txt")
         os.system("ls *.SAS >SAS.txt")
         with open("SAS.txt",'r') as f:
            sasname=f.readline()
         print(sasname)
         os.environ['SAS_ODF'] = path + obsID + "/cal/"+sasname[0:-1]
         # # ---------------------------------------------
         # -----------processing for evt2 file----------
         os.system("emchain")
         os.system("epchain")
         # # ---------------------------------------------
         #------------rename file-----------------------
         os.system("rm *FIT.txt")
         v1=os.popen("ls *M1*EVLI*.FIT")
         mos1name=v1.read()[0:-1]
         v2=os.popen("ls *M2*EVLI*.FIT")
         mos2name=v2.read()[0:-1]
         v3=os.popen("ls *PI*EVLI*.FIT")
         pnname=v3.read()[0:-1]

         cmd="cp "+mos1name+" mos1.fits"
         os.system(cmd)
         cmd="cp "+mos2name+" mos2.fits"
         os.system(cmd)
         cmd="cp "+pnname+" pn.fits"
         os.system(cmd)
         # # ---------------------------------------------
         # #-----------barycen----------------------------
         for det in detList:
            cmd="cp "+det+".fits "+det+"_bary.fits"
            print(cmd)
            os.system(cmd)
            cmd="barycen withtable=yes table="+det+"_bary.fits"+": withsrccoordinates=yes srcra="+str(ra)+" srcdec="+str(dec)
            print(cmd)
            os.system(cmd)
         # # ---------------------------------------------
   #---------choose region-------------------
   #You should open pn.fits and choose the region, saved in physical coordinate.
   # # ---------------------------------------------
   #---------reg def------------------------------------
   # ---------------------------------------------
   srcName = "HM_CnC"
   srcReg = "circle(26710.9,27840.9,400.00)"
   bkgReg = "circle(28029.7,25976.26,800.00026)"

   ## 改成你所选取的region ##

   if filt_particle_bkg:
      pn_threshold=0.5;mos_threshold=0.4
      for det in detList:
         print("running obsservation"+" "+obsID)
         print("running detector"+" "+det)

         datapath = path+obsID+"/cal/"
         print(datapath)
         os.environ['SAS_CCF'] = path + obsID + "/cal/ccf.cif"
         os.chdir(datapath)

         if det == "pn":
             cmd = "evselect table=pn_bary.fits withrateset=Y rateset=rateEPIC_pn.fits maketimecolumn=Y " \
                   "timebinsize=100 makeratecolumn=Y expression='#XMMEA_EP && (PI>10000&&PI<12000) && (PATTERN==0)'"

         else:
            cmd = "evselect table={0}_bary.fits withrateset=Y rateset=rateEPIC_{0}.fits maketimecolumn=Y " \
                   "timebinsize=100 makeratecolumn=Y expression='#XMMEA_EM && (PI>10000) && (PATTERN==0)'".format(det)
         print(" ")
         print("1 make lc")
         print(cmd)
         os.system(cmd)

         if det == "pn":
             cmd = "tabgtigen table=rateEPIC_pn.fits expression='RATE<={0}' gtiset=gti_pn.fits".format(pn_threshold)

         else:
            cmd = "tabgtigen table=rateEPIC_{0}.fits expression='RATE<={1}' gtiset=gti_{0}.fits".format(det,mos_threshold)
         print(" ")
         print("2 make GTI file")
         print(cmd)
         os.system(cmd)


         if det == "pn":
             cmd = "evselect table=pn_bary.fits withfilteredset=yes expression='#XMMEA_EP && gti(gti_pn.fits,TIME) && (PI>150)' " \
                   "filteredset=pn_filt_time_bary.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes"
         else:
            cmd = "evselect table={0}_bary.fits withfilteredset=yes expression=" \
                  "'#XMMEA_EM && gti(gti_{0}.fits,TIME) && (PI>150)' filteredset={0}_filt_time_bary.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes".format(det)

         print(" ")
         print("3 filter GTI fits")
         print(cmd)
         os.system(cmd)

   if lc:
      # you should run this step multiple times to determine the best tmin and tmax
      lenbin=5;tmin=0;tmax=0
      ## you can also specify the tmin and tmax (read from lc)
      ## otherwise they would be total duration
      for det in detList:

         print("running obsservation"+" "+obsID)
         print("running detector"+" "+det)

         datapath = path+obsID+"/cal/"
         print(datapath)
         os.environ['SAS_CCF'] = path + obsID + "/cal/ccf.cif"
         os.chdir(datapath)
         ###----define SAS_ODF----------------###
         with open("SAS.txt",'r') as f:
            sasname=f.readline()
         print(sasname)
         os.environ['SAS_ODF'] = path + obsID + "/cal/"+sasname[0:-1]
         ##------------------##
         if tmin==0 and tmax==0:
            event=fits.open(det+"_filt_time_bary.fits")
            tmin=event[1].header['TSTART']
            tmax=event[1].header['TSTOP']

         if det == "pn":
             cmd = "evselect table="+det+"_filt_time_bary.fits energycolumn=PI " \
                                                  "expression='#XMMEA_EP && (PATTERN<=4) && ((X,Y) IN "+srcReg+ ")"+" &&(PI in [200:10000])' withrateset=yes rateset="\
                   +det+"_src_lc_bin{0}.lc timebinsize={0} maketimecolumn=yes makeratecolumn=yes timemin={1} timemax={2}".format(lenbin,tmin,tmax)
         else:
            cmd = "evselect table=" + det + "_filt_time_bary.fits energycolumn=PI " \
                                                       "expression='#XMMEA_EM && (PATTERN<=12) && ((X,Y) IN " + srcReg + ")"+" &&(PI in [200:10000])' withrateset=yes rateset=" \
                  + det + "_src_lc_bin{0}.lc timebinsize={0} maketimecolumn=yes makeratecolumn=yes timemin={1} timemax={2}".format(lenbin, tmin, tmax)
         print(" ")
         print("1 extract src light curve")
         print(cmd)
         os.system(cmd)

         if det == "pn":
             cmd = "evselect table="+det+"_filt_time_bary.fits energycolumn=PI " \
                                                  "expression='#XMMEA_EP && (PATTERN<=4) && ((X,Y) IN "+bkgReg+")"+" &&(PI in [200:10000])' withrateset=yes rateset="\
                   +det+"_bkg_lc_bin{0}.lc timebinsize={0} maketimecolumn=yes makeratecolumn=yes timemin={1} timemax={2}".format(lenbin,tmin,tmax)
         else:
            cmd = "evselect table=" + det + "_filt_time_bary.fits energycolumn=PI " \
                                                       "expression='#XMMEA_EM && (PATTERN<=12) && ((X,Y) IN " + bkgReg + ")"+" &&(PI in [200:10000])' withrateset=yes rateset=" \
                  + det + "_bkg_lc_bin{0}.lc timebinsize={0} maketimecolumn=yes makeratecolumn=yes timemin={1} timemax={2}".format(lenbin, tmin, tmax)
         print(" ")
         print("2 extract bkg light curve")
         print(cmd)
         os.system(cmd)


         cmd ="epiclccorr srctslist={0}_src_lc_bin{1}.lc eventlist={0}_filt_time_bary.fits outset={0}_lccorr_bin{1}.lc bkgtslist={0}_bkg_lc_bin{1}.lc withbkgset=yes applyabsolutecorrections=yes".format(det,lenbin)

         print(" ")
         print("3 extract corrected light curve")
         print(cmd)
         os.system(cmd)

   if spectra:
      for det in detList:
         print("running obsservation"+" "+obsID)
         print("running detector"+" "+det)
         filt_label='_filt_time_bary'

         datapath = path+obsID+"/cal/"
         print(datapath)
         os.environ['SAS_CCF'] = path + obsID + "/cal/ccf.cif"

         if det == "pn":
             cmd = "evselect table="+datapath+det+filt_label+".fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [200:15000])&&#XMMEA_EP' filteredset="+datapath+det+filt_label+"_filt.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes"
         else:
             cmd = "evselect table="+datapath+det+filt_label+".fits withfilteredset=yes expression='(PATTERN <= 12)&&(PI in [200:12000])&&#XMMEA_EM' filteredset="+datapath+det+filt_label+"_filt.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes"
         print(" ")
         print("1 filter by energy")
         print(cmd)
         os.system(cmd)
         cmd = "evselect table="+datapath+det+filt_label+"_filt.fits withimageset=yes imageset="+datapath+det+filt_label+"_filt.img xcolumn=X ycolumn=Y imagebinning=imageSize ximagesize=600 yimagesize=600"
         print(" ")
         print("2 create det image")
         print(cmd)
         os.system(cmd)
         cmd = "evselect table='"+datapath+det+filt_label+"_filt.fits:EVENTS' destruct=false withfilteredset=true withimageset=true imageset="+datapath+det+filt_label+"_detmap.ds xcolumn=DETX ycolumn=DETY #withxranges=true ximagemin=-1500 ximagemax=1500 withyranges=true #yimagemin=-1500 yimagemax=1500 imagebinning='imageSize' #ximagesize=20 yimagesize=20 #writedss=true updateexposure=true"
         print(" ")
         print("3 create det map")
         print(cmd)
         os.system(cmd)

         if det == "pn":
             cmd = "evselect table='"+datapath+det+filt_label+"_filt.fits:EVENTS' withspectrumset=yes " \
                                                              "spectrumset="+datapath+det+"_"+srcName+".pha energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479 expression='#XMMEA_EP && (PATTERN<=4) && ((X,Y) IN "+srcReg+")'"
         else:
             cmd = "evselect table='"+datapath+det+filt_label+"_filt.fits:EVENTS' withspectrumset=yes " \
                                                              "spectrumset="+datapath+det+"_"+srcName+".pha energycolumn=PI spectralbinsize=15 withspecranges=yes specchannelmin=0 specchannelmax=11999 expression='#XMMEA_EM && (PATTERN<=12) && ((X,Y) IN "+srcReg+")'"
         print(" ")
         print("4 extract source spectrum")
         print(cmd)
         os.system(cmd)
         if det == "pn":
             cmd = "evselect table='"+datapath+det+filt_label+"_filt.fits:EVENTS' withspectrumset=yes " \
                                                              "spectrumset="+datapath+det+"_BKG_for"+srcName+".pha energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479 expression='#XMMEA_EP && (PATTERN<=4) && ((X,Y) IN "+bkgReg+")'"
         else:
             cmd = "evselect table='"+datapath+det+filt_label+"_filt.fits:EVENTS' withspectrumset=yes " \
                                                              "spectrumset="+datapath+det+"_BKG_for"+srcName+".pha energycolumn=PI spectralbinsize=15 withspecranges=yes specchannelmin=0 specchannelmax=11999 expression='#XMMEA_EM && (PATTERN<=12) && ((X,Y) IN "+bkgReg+")'"
         print(" ")
         print("5 extract background spectrum")
         print(cmd)
         os.system(cmd)

         cmd = "backscale spectrumset="+datapath+det+"_"+srcName+".pha badpixlocation="+datapath+det+filt_label+"_filt.fits"
         print(" ")
         print("6 create source backscale keyword")
         print(cmd)
         os.system(cmd)
         cmd = "backscale spectrumset="+datapath+det+"_BKG_for"+srcName+".pha badpixlocation="+datapath+det+filt_label+"_filt.fits"
         print(" ")
         print("7 create background backscale keyword")
         print(cmd)
         os.system(cmd)

         cmd = "rmfgen spectrumset="+datapath+det+"_"+srcName+".pha rmfset="+datapath+det+"_"+srcName+".rmf"
         print(" ")
         print("8 create rmf")
         print(cmd)
         os.system(cmd)
         cmd = "arfgen spectrumset="+datapath+det+"_"+srcName+".pha arfset="+datapath+det+"_"+srcName+".arf withrmfset=yes " \
                                                                                                      "rmfset="+datapath+det+"_"+srcName+".rmf badpixlocation="+datapath+det+filt_label+"_filt.fits extendedsource=yes detmaptype=dataset detmaparray="+datapath+det+filt_label+"_detmap.ds"
         print(" ")
         print("9 create arf")
         print(cmd)
         os.system(cmd)
         cmd="fparkey " +det+"_BKG_for"+srcName+".pha " +datapath+det+"_"+srcName+".pha " +"BACKFILE add=yes"
         print("10 add key")
         print(" ")
         print(cmd)
         os.system(cmd)
         cmd="fparkey "+det+"_"+srcName+".rmf " +datapath+det+"_"+srcName+".pha " +"RESPFILE add=yes"
         print(" ")
         print(cmd)
         os.system(cmd)
         cmd="fparkey " +det+"_"+srcName+".arf " + datapath + det + "_" + srcName + ".pha " + "ANCRFILE add=yes"
         print(" ")
         print(cmd)
         os.system(cmd)
         print(" ")
         print(" ")


   ## 合并多个光谱，一般来说你不需要做这一步 ##
   if combine_spec:
      datapath = path + obsID + "/cal/"
      os.chdir(datapath+'add_2obs_spec/')
      srcName='vcc1154_polygon'
      # cmd="epicspeccombine pha="+'"mos1_{0}.pha mos2_{0}.pha mos1_{0}_02.pha mos2_{0}_02.pha"'.format(srcName) \
      #         +" bkg="+'"mos2_BKG_for{0}.pha mos2_BKG_for{0}.pha mos1_BKG_for{0}_02.pha mos2_BKG_for{0}_02.pha"'.format(srcName) \
      #         +" rmf="+'"mos1_{0}.rmf mos2_{0}.rmf mos1_{0}_02.rmf mos2_{0}_02.rmf"'.format(srcName) \
      #         +" arf="+'"mos1_{0}.arf mos2_{0}.arf mos1_{0}_02.arf mos2_{0}_02.arf"'.format(srcName) \
      #         +" filepha="+'"src_spec_grp_mos_2obs.pha"'\
      #         +" filebkg="+'"bkg_spec_grp_mos_2obs.pha"'\
      #         +" filersp="+'"response_grp_mos_2obs.rmf"'

      cmd="epicspeccombine pha="+'"mos1_{0}.pha mos2_{0}.pha mos1_{0}_obs2.pha mos2_{0}_obs2.pha "'.format(srcName) \
              +" bkg="+'"mos1_BKG_for{0}.pha mos2_BKG_for{0}.pha mos2_BKG_for{0}_obs2.pha mos2_BKG_for{0}_obs2.pha "'.format(srcName) \
              +" rmf="+'"mos1_{0}.rmf mos2_{0}.rmf mos1_{0}_obs2.rmf mos2_{0}_obs2.rmf"'.format(srcName) \
              +" arf="+'"mos1_{0}.arf mos2_{0}.arf mos1_{0}_obs2.arf mos2_{0}_obs2.arf"'.format(srcName) \
              +" filepha="+'"merge_2obs_mos_src.pha"'\
              +" filebkg="+'"merge_2obs_mos_bkg.pha"'\
              +" filersp="+'"response_merge_2obs_mos.rmf"'

      # cmd="epicspeccombine pha="+'"pn_{0}.pha pn_{0}_obs2.pha "'.format(srcName) \
      #         +" bkg="+'"pn_BKG_for{0}.pha pn_BKG_for{0}_obs2.pha "'.format(srcName) \
      #         +" rmf="+'"pn_{0}.rmf pn_{0}_obs2.rmf"'.format(srcName) \
      #         +" arf="+'"pn_{0}.arf pn_{0}_obs2.arf"'.format(srcName) \
      #         +" filepha="+'"merge_2obs_pn_src.pha"'\
      #         +" filebkg="+'"merge_2obs_pn_bkg.pha"'\
      #         +" filersp="+'"response_merge_2obs_pn.rmf"'

      print(cmd)
      os.system(cmd)
      cmd = "fparkey " + "response_merge_2obs_mos.rmf " + "merge_2obs_mos_src.pha" + " RESPFILE add=yes"
      print(cmd)
      os.system(cmd)
      cmd = "fparkey " + "merge_2obs_mos_bkg.pha " + "merge_2obs_mos_src.pha" + " BACKFILE add=yes"
      print(cmd)
      os.system(cmd)
      cmd = "fparkey " + "None " + "merge_2obs_mos_src.pha" + " ANCRFILE add=yes"
      print(cmd)
      os.system(cmd)