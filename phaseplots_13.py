# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 11:18:24 2022

@author: TheLa
"""

import numpy

import all_hc_functions as hc
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
import mk_lc_fmt_5 as lcfmt
from scipy import stats
import csv

#print((lk_tbl[numpy.where(lk_tbl['slice'] == 18)]))


SMALL_SIZE = 8
MEDIUM_SIZE = 12
BIGGER_SIZE = 12

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

mypath = 'D:/PhD/Codes/stageiii/lc_118_short/'
ids = numpy.genfromtxt('lc_118_pids.txt', dtype = str)
lids = numpy.genfromtxt('lc_118_lpids.txt', dtype = 'int64')
fam = numpy.load('118_finamp.npy')


lk_tbl = numpy.load('D:/PhD/Codes/stageiii/lc_118_short/lc_sliced/lookup/lk_table.npy')
bdates  = lcfmt.make_boundary_slices()

#takes in amps
def circ_mean_phase(ampp):
    slcs = numpy.unique(ampp['slice'])
    farr = numpy.full(len(slcs), 99.9)
    std = numpy.full(len(slcs), 99.9)
   # print((slcs))
    for s in range(len(slcs)):
        amps = ampp['phmax'][numpy.where(ampp['slice'] == slcs[s])]
        farr[s] = stats.circmean(amps, high = 1, low = 0)
        std[s] = stats.circstd(amps, high = 1, low = 0)
      #  print(amps)
    #print(len(farr))
    return(slcs, farr, std)


def circ_mean_phase2(st_data):
	slcs = numpy.unique(st_data['slice'])
	farr = numpy.full(len(slcs), 99.9)
	std = numpy.full(len(slcs), 99.9)
    
	for s in range(len(slcs)):
		sl_data = st_data[numpy.where(st_data['slice'] == slcs[s])]
		ph_data = [sl_data['phmaxv'],sl_data['phmaxr'],sl_data['phmaxi']]
		farr[s] = stats.circmean(ph_data, high = 1, low = 0)
		std[s] = stats.circstd(ph_data, high = 1, low = 0)
	return(slcs, farr, std)

#works in deg
def circ_mean_phase3(phv, phr, phi):
    pharr = numpy.full(len(phv), 99.9)
    std = numpy.full(len(phv), 99.9)
    
    for p in range(len(phv)):

        ph_data = [phv[p], phr[p], phi[p]]
        pharr[p] = stats.circmean(ph_data, high = 180, low = -180)
        std[p] = stats.circstd(ph_data, high = 180, low = -180)
    return(pharr, std)

#indexes every 4th slice
def slices4(arr):
    min_value = min(arr)
    max_value = max(arr)
    
    result = [num for num in range(min_value, max_value + 1) if num % 4 == 0]
    return result
    
#MAKE NEW FAM WITH SNR AND PHASE CUTS!

def snr_phase_fam(mypath, fam):
    
    gids = numpy.unique(fam['gid'])
    refam = fam[0]
    
    #make an initial array to append to
    for i in range(len(refam)):
        refam[i] = 99
        
    print(refam)
        
    for s in range(len(gids)):
        gid = gids[s]
        st_data = fam[numpy.where(fam['gid'] == gid)]
        Teff = st_data['Teff'][0]
        per = st_data['finper'][0]
        a = (len(st_data))
        
#### remove any slices that have less than snr = 3 in all filters. Slow but works. 
        slcs2 = [] 	
        
        for slc in numpy.unique(st_data['slice']):
            sl_data = st_data[numpy.where(st_data['slice'] == slc)]
            
            if min(sl_data['snrv'],sl_data['snrr'],sl_data['snri']) >= 3.:
                slcs2.append(slc)
                
        print(slcs2)
        
        slc2mask = numpy.isin(st_data['slice'], slcs2)
        st_data = st_data[slc2mask]
        
        slcs, pharr, std = circ_mean_phase2(st_data)            
             
        chslc = slcs[numpy.where(std <= 0.125) ]
        print('chslc', chslc)
        slmask = numpy.isin(st_data['slice'], chslc)
        
        st_data = st_data[slmask]
        
        b = (len(st_data))
        print('object', s+1, 'precut', a, 'post', b)
        
        print(gid, st_data['slice'])
        refam = numpy.append(refam, st_data)
        
    refam = refam[numpy.where(refam['gid'] > 100)] #remove the intial 99 array
  #  print(refam['slice'][numpy.where(refam['gid'] == 2163156056685634944)])   
    return(refam)

#find the slices with bad snr
def snr_phase_fam_snrtrash(mypath, fam):
    
    gids = numpy.unique(fam['gid'])
    
    snrtfam = fam[0]
    #make an initial array to append to
    for i in range(len(snrtfam)):
        snrtfam[i] = 99
        
    print(snrtfam)
        
    for s in range(len(gids)):
        gid = gids[s]
        st_data = fam[numpy.where(fam['gid'] == gid)]
        #print(st_data)
		#print(st_data)
        Teff = st_data['Teff'][0]
        per = st_data['finper'][0]
        a = (len(st_data))
        
        	#### remove any slices that have less than snr = 3 in all filters. Slow but works. 
        slcs2 = [] 	
        
        for slc in numpy.unique(st_data['slice']):
            sl_data = st_data[numpy.where(st_data['slice'] == slc)]
            
            if min(sl_data['snrv'],sl_data['snrr'],sl_data['snri']) <= 3.:
                slcs2.append(slc)
                
        #print(slcs2)
        
        slc2mask = numpy.isin(st_data['slice'], slcs2)
        st_data = st_data[slc2mask]
        
        slcs, pharr, std = circ_mean_phase2(st_data)            
             
        
        print(len(st_data))
        b = (len(st_data))
        print('ob', s+1, 'precut', a, 'post', b)
        
        print(gid, st_data['slice'])
        snrtfam = numpy.append(snrtfam, st_data)
        
    snrtfam = snrtfam[numpy.where(snrtfam['gid'] > 100)] #remove the intial 99 array
  #  print(snrtfam['slice'][numpy.where(snrtfam['gid'] == 2163156056685634944)])   
    return(snrtfam)

#slices with good snr but bad phase
def snr_phase_fam_phasetrash(mypath, fam):
    
    gids = numpy.unique(fam['gid'])
    
    phtfam = fam[0]
    #make an initial array to append to
    for i in range(len(phtfam)):
        phtfam[i] = 99
        
    print(phtfam)
        
    for s in range(len(gids)):
        gid = gids[s]
        st_data = fam[numpy.where(fam['gid'] == gid)]
        #print(st_data)
		#print(st_data)
        Teff = st_data['Teff'][0]
        per = st_data['finper'][0]
        a = (len(st_data))
        
        #### remove any slices that have less than snr = 3 in all filters. Slow but works. 
        slcs2 = [] 	
        
        for slc in numpy.unique(st_data['slice']):
            sl_data = st_data[numpy.where(st_data['slice'] == slc)]
            
            if min(sl_data['snrv'],sl_data['snrr'],sl_data['snri']) >= 3.:
                slcs2.append(slc)
                
        
        slc2mask = numpy.isin(st_data['slice'], slcs2)
        st_data = st_data[slc2mask]
        
        
        slcs, pharr, std = circ_mean_phase2(st_data)    
        chslc = slcs[numpy.where(std >= 0.125) ] #max seperation between filters!
        print(slcs, pharr, std)
        print('slices where over limit', chslc)
        
        slmask = numpy.isin(st_data['slice'], chslc)
        
        st_data = st_data[slmask]
        
        print(len(st_data))
        b = (len(st_data))
        
        print('ob', s+1, 'precut', a, 'post', b)
        
        print(gid, st_data['slice'])
        phtfam = numpy.append(phtfam, st_data)
        
    phtfam = phtfam[numpy.where(phtfam['gid'] > 100)] #remove the intial 99 array
    print(phtfam['slice'])#[numpy.where(phtfam['gid'] == 2163156056685634944)])   
    return(phtfam)

#
#refam = snr_phase_fam(mypath, fam) #good snr and phase

#snrt = snr_phase_fam_snrtrash(mypath, fam)
#numpy.save('118_finamp_snrtrash.npy', snrt)

refam = numpy.load('118_finamp_snrph.npy') #good snr and phase
snrt = numpy.load('118_finamp_snrtrash.npy') #bad snr
pht = numpy.load('118_finamp_pht.npy') #bad phase


def save_as_file(mypath, w):
	with open(mypath+'phasevalues.txt', 'a+') as file:
		file.write(str(w))
		file.write('\n')
		file.close()
	return()

#to create 4 panel plots with for each object  data points, phase, spot temp and spot size. 

def plot_spot_prop_amp_4_4(mypath, fam, snrt, pht): #versions_4 and _3 are in phaseplots_11  in archive. 
	gids = numpy.unique(fam['gid'])
	#gids = gids[5:6]
	cols = ['green', 'red', 'black']
	mas = [5,4,3]
    
	for s in range(len(gids)):
		plt.clf()
		plt.close()
		gid = gids[s]
        
		st_data = fam[numpy.where(fam['gid'] == gid)]
		snrt_st_data = snrt[numpy.where(snrt['gid'] == gid)]
		pht_st_data = pht[numpy.where(pht['gid'] == gid)]
        
        
		Teff = st_data['Teff'][0]
		per = st_data['finper'][0]
	    
		fig, (ax1,ax2,ax3, ax4) = plt.subplots(4, 1, sharex = False, figsize = (6, 8), gridspec_kw={'hspace': 0})
            
		ax1.set_title('P = '+str("%.4f"%per)+' days')
		ax1.set_ylabel('Amplitude [mag]')
		ax4.set_xlabel('MJD')
		ax2.set_ylabel('Phase [deg]')
		ax3.set_ylabel('Spot Temperature (K)')
		ax4.set_ylabel('Spot Coverage')
    

		aa = numpy.concatenate((st_data['slice'], snrt_st_data['slice'], pht_st_data['slice']))
		aa = numpy.unique(aa)
		yrs = slices4(aa)
		aa = lcfmt.slice_to_mjd(aa)
		print(yrs)
		yrsmjd = lcfmt.slice_to_mjd(yrs)
        
		ax1.errorbar(lcfmt.slice_to_mjd(st_data['slice']), st_data['finampv'], st_data['finampve'], color = 'green', ls = 'none', marker = 'o', markersize = 4, rasterized = True) 
		ax1.errorbar(lcfmt.slice_to_mjd(st_data['slice']), st_data['finampr'], st_data['finampre'], color = 'red', ls = 'none', marker = 'o', markersize = 4, rasterized = True)  
		ax1.errorbar(lcfmt.slice_to_mjd(st_data['slice']), st_data['finampi'], st_data['finampie'], color = 'black', ls = 'none', marker = 'o', markersize = 4, rasterized = True) 
		
		ax1.errorbar(lcfmt.slice_to_mjd(snrt_st_data['slice']), 3*snrt_st_data['finampve'] , snrt_st_data['finampve'], uplims = True, color = 'green', alpha = 0.4, ls = 'none', marker = '_',  rasterized = True) 
		ax1.errorbar(lcfmt.slice_to_mjd(snrt_st_data['slice']), 3*snrt_st_data['finampre'] , snrt_st_data['finampre'], uplims = True, color = 'red', alpha = 0.4, ls = 'none', marker = '_',  rasterized = True) 
		ax1.errorbar(lcfmt.slice_to_mjd(snrt_st_data['slice']), 3*snrt_st_data['finampie'] , snrt_st_data['finampie'], uplims = True, color = 'black', alpha = 0.4, ls = 'none', marker = '_',  rasterized = True) 
		
		ax1.plot(lcfmt.slice_to_mjd(pht_st_data['slice']), 3*pht_st_data['finampve'], linestyle = ' ', marker  = '+', color = 'green', alpha = 0.4,markersize = 10)
		ax1.plot(lcfmt.slice_to_mjd(pht_st_data['slice']), 3*pht_st_data['finampre'], linestyle = ' ', marker  = '+', color = 'red', alpha = 0.4,markersize = 10)
		ax1.plot(lcfmt.slice_to_mjd(pht_st_data['slice']), 3*pht_st_data['finampie'], linestyle = ' ', marker  = '+', color = 'black', alpha = 0.4,markersize = 10)
        
		ax1.set_xlim(min(aa)-91, max(aa)+91)
		ax2.set_xlim(min(aa)-91, max(aa)+91)
		ax3.set_xlim(min(aa)-91, max(aa)+91)
		ax4.set_xlim(min(aa)-91, max(aa)+91)
        
		ax2.set_xticks([])
		ax3.set_xticks([])
        
        #set limits on maximum/min ampliotudes 
		allam = numpy.concatenate((st_data['finampv'], st_data['finampr'], st_data['finampi']))
		allame = numpy.concatenate((st_data['finampve'], st_data['finampre'], st_data['finampie']))       
		argall = numpy.argmax(allam)
        
		ax1.set_ylim(0, (allam[argall]+allame[argall])*1.2)
    

		allph = numpy.concatenate((st_data['phmaxv'],st_data['phmaxr'],st_data['phmaxi']))
		medph = numpy.median(allph)*360
        
		rvphase = (st_data['phmaxv']*360 - medph + 180) % 360 - 180	
		rrphase = (st_data['phmaxr']*360 - medph + 180) % 360 - 180
		riphase = (st_data['phmaxi']*360 - medph + 180) % 360 - 180
        
		ax2.plot(lcfmt.slice_to_mjd(st_data['slice']), rvphase, color = 'green', ls = 'none', marker = 'o', markersize = 4, rasterized = True) 
		ax2.plot(lcfmt.slice_to_mjd(st_data['slice']), rrphase, color = 'red', ls = 'none', marker = 'o', markersize = 4, rasterized = True) 
		ax2.plot(lcfmt.slice_to_mjd(st_data['slice']), riphase, color = 'black', ls = 'none', marker = 'o', markersize = 4, rasterized = True) 
        
		pharr, std = circ_mean_phase3(rvphase, rrphase, riphase)
        
		ax2.fill_between(lcfmt.slice_to_mjd(st_data['slice']),  pharr - std, pharr + std, color = 'b', alpha = 0.3, rasterized = True)   
		ax2.plot(lcfmt.slice_to_mjd(st_data['slice']),  pharr,  color = 'b', linewidth = 0.5, rasterized = True)
        
		for m in range(len(std)):
				save_as_file(mypath, std[m])
            
		ax2.set_ylim(-190, 190)
        
		hsch = numpy.where(((st_data['hs_ratio']  < 0.4) | (st_data['hs_ratio'] > 0.6)) & (st_data['sp_temp_med'] > 2200) & (st_data['sp_size_med'] < 0.45))
		#print(st_data['sp_temp_med'][hsch[0]])#
        
		if len(st_data['sp_size_med'][hsch[0]]) > 0 :
    		
				ax3.errorbar(lcfmt.slice_to_mjd(st_data['slice'][hsch[0]]), st_data['sp_temp_med'][hsch[0]],st_data['sp_temp_mad'][hsch[0]], color = 'black', ls = 'none', marker = 'o', markersize = 1, rasterized = True, zorder = 10)       
				ax4.errorbar(lcfmt.slice_to_mjd(st_data['slice'][hsch[0]]), st_data['sp_size_med'][hsch[0]],st_data['sp_size_mad'][hsch[0]], color = 'black', ls = 'none', marker = 'o', markersize = 1, rasterized = True, zorder = 10)       
             
				ax3.scatter(lcfmt.slice_to_mjd(st_data['slice'][hsch[0]]), st_data['sp_temp_med'][hsch[0]], c = st_data['hs_ratio'][hsch[0]], cmap = 'seismic_r',  s = 26, rasterized = True, zorder = 40, vmin = 0.0, vmax = 1.0)       
				ax4.scatter(lcfmt.slice_to_mjd(st_data['slice'][hsch[0]]), st_data['sp_size_med'][hsch[0]], c = st_data['hs_ratio'][hsch[0]], cmap = 'seismic_r', s = 26, rasterized = True, zorder = 40, vmin = 0.0, vmax = 1.0)   
      
				cbar_ax = fig.add_axes([0.92, 0.11, 0.02, 0.77])  # Adjust position and size as needed #cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.02]) horizontal
				cmap = plt.cm.get_cmap('seismic_r')
				sm = plt.cm.ScalarMappable(cmap=cmap)
				sm.set_array([])
				colorbar3 = plt.colorbar(mappable = sm, cax=cbar_ax, orientation='vertical', cmap= cmap)
				colorbar3.set_label('HS:CS$_{\{V\}}$')
    
    		
				A = [min(aa)-1000, max(aa)+1000]
				B = [Teff, Teff]
				ax3.plot(A, B, linestyle = 'dotted', color = 'black', linewidth = 0.5 , rasterized = True)    
        
				maxsize = numpy.argmax(st_data['sp_size_med'][hsch[0]])
				maxt = numpy.argmax(st_data['sp_temp_med'][hsch[0]])
				mint = numpy.argmin(st_data['sp_temp_med'][hsch[0]])
        
				ax3.set_ylim((st_data['sp_temp_med'][hsch[0]][mint] - st_data['sp_temp_mad'][hsch[0]][mint])*0.8 , (st_data['sp_temp_med'][hsch[0]][maxt] + st_data['sp_temp_mad'][hsch[0]][maxt])*1.2)
        
        
				ax4.set_ylim(0, (st_data['sp_size_med'][hsch[0]][maxsize]+st_data['sp_size_mad'][hsch[0]][maxsize])*1.2)
		else: 
				ax3.set_ylim(0,10000)				
				ax4.set_ylim(0,0.5)

		for yr in yrs:
						ax1.plot([lcfmt.slice_to_mjd(yr), lcfmt.slice_to_mjd(yr)], [0, 3], color = 'gray', linestyle = 'dashed', linewidth = 0.5)
						ax2.plot([lcfmt.slice_to_mjd(yr), lcfmt.slice_to_mjd(yr)], [-190, 190], color = 'gray', linestyle = 'dashed', linewidth = 0.5)
						ax3.plot([lcfmt.slice_to_mjd(yr), lcfmt.slice_to_mjd(yr)], [100, 14000], color = 'gray', linestyle = 'dashed', linewidth = 0.5)
						ax4.plot([lcfmt.slice_to_mjd(yr), lcfmt.slice_to_mjd(yr)], [0, 0.8], color = 'gray', linestyle = 'dashed', linewidth = 0.5)
                               
        
        
		
        
		for slc in range(len(st_data['slice'])):
				print(s+1, st_data['slice'][slc], 'temp',st_data['sp_temp_med'][slc],'size',st_data['sp_size_med'][slc],'hscs', st_data['hs_ratio'][slc])
		print(yrs)
        
		ax1.xaxis.set_ticks_position('top')
		ax1.set_xticks(lcfmt.slice_to_mjd(yrs))
		ax1.set_xticklabels(yrs)
        
		ax4.xaxis.set_ticks_position('bottom')
		ax4.set_xticks(lcfmt.slice_to_mjd(yrs))
		ax4.set_xticklabels(lcfmt.slice_to_mjd(yrs))
        
		plt.savefig(mypath+'analysis_plots/ampprop_plots/'+str(s+1)+'_spot_properties_amp.png', format = 'png', bbox_inches='tight', dpi=600)
		plt.savefig(mypath+'analysis_plots/ampprop_plots/'+str(s+1)+'_spot_properties_amp.pdf', format = 'pdf', bbox_inches='tight', dpi=600)
		#else: pass
	return()

#plot_spot_prop_amp_4_4(mypath, refam, snrt, pht)
