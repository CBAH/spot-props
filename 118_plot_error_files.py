"""
Created on Tue Nov 23 11:23:15 2021

@author: TheLa
"""


import os
import numpy
#import all_hc_functions as hc
import matplotlib.pyplot as plt
import matplotlib.tri as tri

#import pysynphot as S
#import speclite.filters
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
import mk_lc_fmt_5 as lcfmt


now = datetime.datetime.now()
print ("Current date and time : ")
print (now.strftime("%Y-%m-%d %H:%M:%S"))

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

mypath = 'D:/PhD/Codes/stageiii/lc_118_short/'
#mypath = '/data3/cbah1/stageiii/lc_118_short/'


fam = numpy.load('118_finamp.npy')

amptabletype = numpy.dtype([('number',int),
        ('st_model', numpy.unicode_, 16),
        ('sp_model', numpy.unicode_, 16),#number of column
    	('st_temp',float),               #temperatureof star [K]
    	('sp_temp',float),              #temperature of spot [K]
    	('sp_size',float), #size of spot in fraction
        ('met', float),       #metalicity [M/H]
        ('logg', float),    #log g      
        ('amp_U',float),               #mag in U
        ('amp_B',float),               #mag in B
        ('amp_V',float),               #mag in V
        ('amp_R',float),               #mag in R
        ('amp_I',float),   
        ])#
	    
rms_tabletype = numpy.dtype([('id',int),		#running integer ID of star []
                                 ('og_amp_v', float),
                                 ('og_amp_r', float),
                                 ('og_amp_i', float),
                                 ('og_ampe_v', float),
                                 ('og_ampe_r', float),
                                 ('og_ampe_i', float),
                                 ('n_amp_v', float),
                                 ('n_amp_r', float),
                                 ('n_amp_i', float),
                                 ('Teff', float),
                                 ('st_model', numpy.unicode_, 16),
                                 ('sp_model', numpy.unicode_, 16),
                                 ('min_rms_irv', float),
                                 ('sp_size_irv', float),
                                 ('sp_temp_irv', float),
                                 ])
				 

#matplotlib.rcParams['text.usetex'] = True

def find_nearest(array, value):
    array = numpy.asarray(array)
    idx = (numpy.abs(array - value)).argmin()
    return array[idx]
                    	

#plots the files from 118_make_error file
def plt_error_table(mypath, st_id, sp_model, s):
    
    plt.clf()
    plt.close()
    print('plotting error plot for', st_id)
   
    data = numpy.genfromtxt(mypath+'error_files/'+str(int(st_id))+'/s'+str(s)+'.txt', dtype = rms_tabletype, delimiter = ',')
    
    select = numpy.where(data['sp_model'] == sp_model)
    
    data = data[select[0]]
    print((data[0]))

    Teff = data['Teff'][0]
    
    rms_list = []
    
    data = data[numpy.where((abs(data['og_amp_i'] - data['n_amp_i'] ) < data['og_ampe_i']))]
    data = data[numpy.where((abs(data['og_amp_r'] - data['n_amp_r'] ) < data['og_ampe_r']))]
    data = data[numpy.where((abs(data['og_amp_v'] - data['n_amp_v'] ) < data['og_ampe_v']))]
    
    og_list = [data['og_amp_i'][0],data['og_amp_r'][0],data['og_amp_v'][0]]
    oge_list = [data['og_ampe_i'][0],data['og_ampe_r'][0],data['og_ampe_v'][0]] 
	
    print('final len of data',len(data))
    		
    for i in range(len(data)):
    
        calc = [data['n_amp_i'][i], data['n_amp_r'][i], data['n_amp_v'][i]]
        denom, num = (((og_list[0] - calc[0])**2)/(oge_list[0]**2))+(((og_list[1] - calc[1])**2)/(oge_list[1]**2))+(((og_list[2] - calc[2])**2)/(oge_list[2]**2)), 3.0
        mse = denom/num
        rms = numpy.sqrt(mse)
        rms_list.append(rms)
       # print(rms)
       
    rectangle1 = [0.1,0.1, 0.8,0.8]
    ax1 = plt.axes(rectangle1)
    ax1.set_xlabel('Spot Coverage')
    ax1.set_xlim(0, 0.5)
    ax1.set_ylabel('Spot Temperature [K]')

    p_s = data['min_rms_irv']*1000. #- point size
    rms_arr = numpy.array(rms_list)
    
    
    sp_size, sp_temp = data['sp_size_irv'], data['sp_temp_irv']
    
    hot_select = numpy.where(sp_temp > Teff)
    cold_select = numpy.where(sp_temp < Teff)
    
    hotspot_ratio = len(sp_temp[hot_select[0]]) / len(sp_temp)
    
    
    #medians of the dominant side - hot or cold 
    if len(sp_temp[hot_select[0]]) < len(sp_temp[cold_select[0]]):
        sp_size_select, sp_temp_select, rms_arr_select, p_s_select = sp_size[cold_select[0]], sp_temp[cold_select[0]], rms_arr[cold_select[0]], p_s[cold_select[0]]

        A1 = [0,0.5]
        B1 = [2000, 2000]
        ax1.plot(A1, B1, linestyle = 'dashed', color = 'gray', linewidth = 0.5 )
            
    else: 
        sp_size_select, sp_temp_select, rms_arr_select, p_s_select = sp_size[hot_select[0]], sp_temp[hot_select[0]], rms_arr[hot_select[0]], p_s[hot_select[0]]
 
    
	
    sort = numpy.argsort(rms_arr)[::-1]
	
    sp_size, sp_temp, rms_arr, p_s = sp_size[sort], sp_temp[sort], rms_arr[sort], p_s[sort]
    
    A = [0,0.5]
    B = [Teff, Teff]
    ax1.plot(A, B, linestyle = 'dotted', color = 'black', linewidth = 0.5 )
    
    sp_size_med, sp_temp_med = numpy.median(sp_size_select), numpy.median(sp_temp_select)
    sp_size_mad, sp_temp_mad = stats.median_abs_deviation(sp_size_select), stats.median_abs_deviation(sp_temp_select)
    
    plt.scatter(sp_size, sp_temp,  c = rms_arr, s = p_s*6, cmap = plt.cm.viridis, vmin = 0, vmax = 1.0,  zorder = 0, rasterized = True)    
    
    ax1.errorbar(sp_size_med, sp_temp_med, yerr = sp_temp_mad, xerr = sp_size_mad, color = 'red', zorder = 50, rasterized = True)
   
    #cbar = plt.colorbar()
    
    plt.savefig(mypath+'error_plots/'+str(st_id)+'/s'+str(s)+'_error_plot.png', format='png', bbox_inches='tight', dpi=600) 
    plt.savefig(mypath+'error_plots/'+str(st_id)+'/s'+str(s)+'_error_plot.pdf', format='pdf', bbox_inches='tight', dpi=600)
   # plt.show()
    print('sp_size_med, sp_temp_med, sp_size_mad, sp_temp_mad, hotspot_ratio',sp_size_med, sp_temp_med,sp_size_mad, sp_temp_mad, hotspot_ratio)
    
    return(sp_size_med, sp_temp_med, sp_size_mad, sp_temp_mad, hotspot_ratio)


## saves the median temp/size and their median absolute uncertainties
def loop_id_error_tbl(mypath, fam):
    ids = numpy.unique(fam['sid'][numpy.where(fam['hid'] > 99)])

   
    for i in range(len(ids)):
        try:
            st_data = fam[numpy.where(fam['sid'] == ids[i])] #st_data
	
            slcs = numpy.unique(st_data['slice'])
	
            for s in slcs:
                sl_data = st_data[numpy.where(st_data['slice'] == s)] #slice_data
		
                sp_size_med, sp_temp_med, sp_size_mad, sp_temp_mad, hotspot_ratio  = plt_error_table(mypath, ids[i], 'ph', s)
                sl_data['sp_size_med'], sl_data['sp_temp_med'], sl_data['sp_size_mad'], sl_data['sp_temp_mad'], sl_data['hs_ratio'] = sp_size_med, sp_temp_med, sp_size_mad, sp_temp_mad, hotspot_ratio

                st_data[numpy.where(st_data['slice'] == s)] = sl_data
            print(st_data)        
            fam[numpy.where(fam['sid'] == ids[i])] = st_data 
        except: continue
    
    numpy.save('118_finamp.npy', fam)
    
    return(fam)


#loop_id_error_tbl(mypath, fam)
