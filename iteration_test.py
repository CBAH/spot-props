# -*- coding: utf-8 -*-
"""
Created on Thu May  2 17:40:38 2024

@author: TheLa

For Section 'Iteration Relibility'

"""



import os
import numpy
#import all_hc_functions as hc
import matplotlib.pyplot as plt
import matplotlib.tri as tri

import pysynphot as S
import speclite.filters
import math
from scipy import stats
from scipy.ndimage.filters import gaussian_filter
import random
import csv
import datetime
import sys
from sklearn.metrics import mean_squared_error
from scipy.interpolate import griddata
from scipy.spatial.distance import pdist
import matplotlib.pylab as pylab
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.transforms import blended_transform_factory
import matplotlib.colors as mcol
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap

now = datetime.datetime.now()
print ("Current date and time : ")
print (now.strftime("%Y-%m-%d %H:%M:%S"))

#Plotting info
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

#original object from stageii of project, files are in original table type text files. 	    
rms_tabletypeog = numpy.dtype([('id',int),		#running integer ID of star []
                                 ('og_amp_u', float),           
                                 ('og_amp_b', float),
                                 ('og_amp_v', float),
                                 ('og_amp_r', float),
                                 ('og_amp_i', float),
                                 ('og_ampe_u', float),
                                 ('og_ampe_b', float),
                                 ('og_ampe_v', float),
                                 ('og_ampe_r', float),
                                 ('og_ampe_i', float),
                                 ('n_amp_u', float),
                                 ('n_amp_b', float),
                                 ('n_amp_v', float),
                                 ('n_amp_r', float),
                                 ('n_amp_i', float),
                                 ('Teff', float),
                                 ('st_model', numpy.unicode_, 16),
                                 ('sp_model', numpy.unicode_, 16),#number of column
                                 ('min_rms_irv', float),
                                 ('min_rms_irvb', float),
                                 ('min_rms_irvbu', float),
                                 ('sp_size_irv', float),
                                 ('sp_temp_irv', float),
                                 ('sp_size_irvb', float),
                                 ('sp_temp_irvb', float),
                                 ('sp_size_irvbu', float),
                                 ('sp_temp_irvbu', float),
                                 ])
				 

mypath = 'D:/PhD/Codes/stageii/'

def errorprop(p, n): #error of proportion --- p = x/n where x number of instances
    return(numpy.sqrt((p*(1-p))/(n)))

#print(errorprop(6, 27))



def ittest(mypath):
    
    plt.clf()
    plt.close()
    
    
    fig, (ax1,ax2, ax3) = plt.subplots(3, 1, sharex = True, figsize = (8,6), gridspec_kw={'hspace': 0})
            
    
    ax3.set_xlabel('Iteration Count')
    ax1.set_ylabel('Spot Temperature [K]')
    ax2.set_ylabel('Spot Coverage [K]')
    ax3.set_ylabel('HS:CS$_{\{V\}}$')
    
    data = numpy.genfromtxt(mypath+'error_files_3/8038.txt', dtype = rms_tabletypeog, delimiter = ',') #uses object 8038 as example object
    Teff = data['Teff'][0]
    
    
    print(data[0])
    
    #numbers to be plotted
    medl = []
    medlmad = []
    meds = []
    medsmad = []
    hsr = []
    stds = []
    
    rang = numpy.arange(0, len(data), 20)
   
    #to see the effect of re-organising the data - changes the shape but not the statistics, why you can't use 'first 200' or similar
  #  sortind = numpy.argsort(data['min_rms_irv'])    
#    data = data[sortind][::-1]
  #  randind = random.shuffle(data)
   
    print('teff', Teff)
    for i in rang:
        data2 = data[0:i+20]
        hs_d = data2['sp_temp_irv'][numpy.where(data2['sp_temp_irv'] > Teff)] #hot spot temp
        hs_s = data2['sp_size_irv'][numpy.where(data2['sp_size_irv'] > Teff)] #hot size
        
        cs_d = data2['sp_temp_irv'][numpy.where(data2['sp_temp_irv'] < Teff)] #cold spot temp
        cs_s = data2['sp_size_irv'][numpy.where(data2['sp_size_irv'] < Teff)] # hot spot size
        
    
        if len(hs_d) > len(cs_d):
            medl.append(numpy.median(cs_d))
            meds.append(numpy.median(cs_s))
            medlmad.append(stats.median_abs_deviation(cs_d))
            medsmad.append(stats.median_abs_deviation(cs_s))
            hsr.append(len(hs_d)/(len(cs_d)+len(hs_d)))
            stds.append(numpy.std(hsr))
            
        if len(hs_d) < len(cs_d):
            medl.append(numpy.median(cs_d))
            meds.append(numpy.median(cs_s))
            medlmad.append(stats.median_abs_deviation(cs_d))
            medsmad.append(stats.median_abs_deviation(cs_s))
            hsr.append(len(hs_d)/(len(cs_d)+len(hs_d)))
            stds.append(numpy.std(hsr))
            
     
    #print('std first 250 iterations', numpy.std(hsr[0:250]), 'std all', numpy.std(hsr))
    
    medl, meds, medlmad, medsmad, stds = numpy.array(medl), numpy.array(meds), numpy.array(medlmad)*0.1, numpy.array(medsmad)*0.1, numpy.array(stds)


    ### !! change limits if not using object 8038
    ax3.set_ylim(0.101, 0.148)
    ax1.set_ylim(3320, 3430)
    ax2.set_ylim(0.139, 0.154)
    ax3.set_xlim(0, 10200)
    
    fhsr = hsr[-1:][0]
    print(fhsr)
    err = errorprop(fhsr, len(data))
    print('hotspot ratio error', err)
    
    ax1.plot(rang, medl, c = 'black', rasterized = True)
    ax1.fill_between(rang, medl+medlmad, medl-medlmad, alpha = 0.6, rasterized = True)
    
    ax2.plot(rang, meds, c = 'black', rasterized = True)
    ax2.fill_between(rang, meds+medsmad, meds-medsmad, alpha = 0.6, rasterized = True)
    
    ax3.plot(rang, hsr, c = 'black', rasterized = True)
    #ax3.fill_between(rang, hsr+stds, hsr-stds, alpha = 0.6, rasterized = True) 
    ax3.errorbar(rang[-1:][0], fhsr, yerr = err, dash_capstyle = 'projecting', lw = 2, alpha = 1, capsize = 6, capthick = 1) #error bar on final point
    
    plt.savefig(mypath+'iteration_test_8038.pdf', format = 'pdf', dpi = 600)
    return()

#ittest(mypath)
