# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 19:35:27 2024

@author: TheLa

"""
#from __future__ import print_function
from astropy.coordinates import match_coordinates_sky , SkyCoord
from astropy import stats
from astropy import units as u
from astropy.timeseries import LombScargle
from matplotlib.font_manager import FontProperties
from operator import itemgetter
from scipy.cluster.hierarchy import linkage, fcluster
from scipy import optimize #Leastsq Levenberg-Marquadt Algorithm
from scipy.optimize import curve_fit
from scipy.spatial import distance_matrix
import matplotlib.pylab as plt
import numpy
from astroquery.vizier import Vizier
import astropy.coordinates as coord
from sklearn.linear_model import LinearRegression

mypath = 'D:/PhD/Codes/stageiv/'
psl = numpy.load(mypath+'mastersl_periodic.npy') #source list all periodic
gin = numpy.load('D:/PhD/Codes/stageiii/118_allgaia_com.npy') #sourcelist ic5070 all objects
fam118 = numpy.load('D:/PhD/Codes/stageiii/118_finamp_snrph.npy') #ic5070 only
fam = numpy.load(mypath+'master_famhi8.npy') #most recent finamps of all regions

#fitting function
def test_T5(X, x0, a, b, c, d, e):
    G, col, c1, c2 = X
    return(x0 + a*G + b*col + c*G*col + d*c1 + e*c2)
                    	
#apparent magitude to abs	
def app_to_abs(m, d):
    return(m - 5*numpy.log10(d/10))

#straight line
def func_line(x, m, c):
    return(m*x + c)

#match Fang catalogue members with source list params - sourcelist, takes four colours which are subtracted thus c11 - c12, c21 - c22 to create 'col1', 'col2'
def match_A2(gin, c11, c12, c21, c22):
    cat_name = "J/ApJ/904/146/table4"
    v = Vizier(columns= ["*", "+_r"], catalog = cat_name)
    #print(v)
    gids = numpy.unique(gin['gid'])
    arrdtype = numpy.dtype([('gid', 'int64'),
                      ('bprp', float),
                      ('gmag', float),
                      ('plx', float),
                      ('FTeff', float),
                      ('col1', float),
                      ('col2', float),
                      ])
    arr = numpy.full(len(gin), 99., dtype = arrdtype)
    
    badg = []
    for g in range(len(gids)):
        st_fin = gin[numpy.where(gin['gid'] == gids[g])]
        
        result = v.query_region(coord.SkyCoord(ra = st_fin['ra'], dec = st_fin['de'], unit = (u.deg, u.deg), frame = 'icrs'), radius = "3s", catalog = cat_name)
      #  print(len(result))
        if len(result) == 1: 
            arr['gid'][g] = st_fin['gid']
            arr['bprp'][g] = st_fin['bprp']
            arr['gmag'][g] = st_fin['gmag']
            arr['plx'][g] = st_fin['plx']
            arr['FTeff'][g] =  result[0]['Teff'][0] 
            
            arr['col1'][g] = st_fin[c11] - st_fin[c12]
            arr['col2'][g] = st_fin[c21] - st_fin[c22]
        
            if (abs(arr['col1'][g]) > 100) | (abs(arr['col2'][g])> 100):
                badg.append(arr['gid'][g])
        
    cdd = 1/(arr['plx']*1e-3)   
    arr['gmag'] = app_to_abs(arr['gmag'], cdd)        
    
    arr = arr[numpy.where(arr['FTeff'] > 100)]   
    arr = arr[numpy.where(arr['FTeff'] < 8000)]   
    arr = arr[numpy.where(arr['col1'] < 80)]    # takes care of any missing colours. 
    arr = arr[numpy.where(arr['col2'] < 80)]   
    return(arr)

#fits parameters excluding the group between 3900 - 4000, returns parameters
def fit_Fang_noplot(arr): 
    arr = arr[numpy.where(((arr['FTeff'] < 3900) | (arr['FTeff'] > 4000)))]
    param, cov = optimize.curve_fit(test_T5, (arr['gmag'], arr['bprp'], arr['col1'], arr['col2']), arr['FTeff'], sigma = 1/arr['FTeff'])    #sigma deweights high tempreatures
    return(param) 

# Recreates green and blue histogram plot of the Fang and fitted temperatures 
def FangFit_hist(gin):
    plt.clf()
    plt.close()
    
    rectangle1 = [0.1,0.1, 0.8,0.8]
    ax1 = plt.axes(rectangle1)
    ymax = 54
    ax1.set_ylim(0, ymax)
    ax1.set_xlim(2800, 8000)
    ax1.set_ylabel('Number of objects')
    ax1.set_xlabel('T$_{eff}$ [K]')
    
    arr = match_A2(gin, 'gmag', 'Jmag', 'W1mag', 'W2mag')
    bns = numpy.arange(2800, 8100, 100)
    bns2 = numpy.arange(2800, 8020, 20)
    
    hist, bins, _ = ax1.hist(arr['FTeff'], bins = bns, color = 'green', alpha = 0.6)
    hist2, bins2 = numpy.histogram(arr['FTeff'], bins = bns2)
    cdf = numpy.cumsum(hist2)
    
    sf = ymax/max(cdf)
    cdf = cdf *sf
    ax2 = plt.twinx()
    ax2.set_ylim(0, 1)
    ax2.set_ylabel('Fraction of objects')
    
    ax1.plot(bins2[:-1], cdf, marker= None, linestyle='-', color='black',linewidth = 1, label='CDF', rasterized = True)
    
    
    gids = numpy.unique(arr['gid'])
    
    params1 = fit_Fang_noplot(arr)
    
    nTs = []
    c11, c12, c21, c22 = 'bpmag', 'Jmag', 'W1mag', 'W2mag'
    
    for i in range(len(gids)):
        gid = gids[i]
        
        narr = arr[numpy.where(arr['gid'] == gid)]
        st_fin = gin[numpy.where(gin['gid'] == gid)]#
        
        cdd = 1/(st_fin['plx']*1e-3)       
        
        nT = test_T5((app_to_abs(st_fin['gmag'], cdd), st_fin['bprp'], (st_fin[c11] - st_fin[c12]), (st_fin[c21] - st_fin[c22])), *params1)#
        
        nTs.append(nT)
        
    nTs = numpy.array(nTs)
    
    hist, bins, _ = ax1.hist(nTs, bins = bns, color = 'blue', alpha = 0.6)
    hist2, bins2 = numpy.histogram(nTs, bins = bns2)
    cdf = numpy.cumsum(hist2)
    
    sf = ymax/max(cdf)
    cdf = cdf *sf
    
    ax1.plot(bins2[:-1], cdf, marker= None, linestyle='dashed', color='black',linewidth = 1, label='CDF', rasterized = True)
    

    plt.savefig('FangFittemp_hist.png', format ='png', dpi = 600)
    plt.savefig('FangFittemp_hist.pdf', format ='pdf', dpi = 600)
    return()

#FangFit_hist(gin)
    
#plots histogram results as a scatter plot. 
def FangFit_scatter(gin):
    plt.clf()
    plt.close()
    
    
    rectangle1 = [0.1,0.1, 0.8,0.8]
    ax1 = plt.axes(rectangle1)
   # ymax = 54
    ax1.set_ylim(2800, 8000)
    ax1.set_xlim(2800, 8000)
    ax1.set_ylabel('T$_{fit}$ [K]')
    ax1.set_xlabel('T$_{eff}$ [K]')
    
    arr = match_A2(gin, 'gmag', 'Jmag', 'W1mag', 'W2mag')
    
    gids = numpy.unique(arr['gid'])
    
    params1 = fit_Fang_noplot(arr)
    
    nTs = []
    c11, c12, c21, c22 = 'bpmag', 'Jmag', 'W1mag', 'W2mag'
    
    for i in range(len(gids)):
        gid = gids[i]
        st_fin = gin[numpy.where(gin['gid'] == gid)]#
        
        cdd = 1/(st_fin['plx']*1e-3)       
        
        nT = test_T5((app_to_abs(st_fin['gmag'], cdd), st_fin['bprp'], (st_fin[c11] - st_fin[c12]), (st_fin[c21] - st_fin[c22])), *params1)#new temperature
         
        nTs.append(nT)
        
        
    nTs = numpy.array(nTs)
    Ts = arr['FTeff']

    nTs = nTs.reshape(-1,1)
    Ts = Ts.reshape(-1, 1)
    
    model = LinearRegression().fit(Ts, nTs)
    AA = [1000, 10000]
    BB = [func_line(1000, model.coef_, model.intercept_)[0], func_line(10000, model.coef_, model.intercept_)[0]] 
    
    ax1.plot(AA, BB, 'black')
    
    
    plt.scatter(Ts, nTs, s = 12, c = 'black')
    
        
    plt.savefig('FangFittemp_scatter.png', format ='png', dpi = 600)
    plt.savefig('FangFittemp_scatter.pdf', format ='pdf', dpi = 600)
    return()

#FangFit_scatter(gin)

def new_CMD(gin, fam): #remakes CMD with fitted temps. 
    
    rectangle1 = [0.1,0.1, 0.8,0.8]
    ax1 = plt.axes(rectangle1)
    ax1.set_xlim(0.5,4.5)
    ax1.set_ylim(18.2,10.6)
    
    oids = numpy.unique(fam['oid'])
    
    
 #   arr = match_Av2(gin, c1, c2)
    arr = match_A2(gin, 'gmag', 'Jmag', 'W1mag', 'W2mag') #c11 - c12, c21 - c22
    params1 = fit_Fang_noplot(arr)
    #2params
    
    gids = numpy.unique(gin['gid'])
    
    t1 = []
    t2 = []
    G = []
    bprp = []
    
    for gid in gids:
        
        narr = arr[numpy.where(arr['gid'] == gid)]
        st_fin = gin[numpy.where(gin['gid'] == gid)]
        
        
        narr = arr[numpy.where(arr['gid'] == gid)]
        st_fin = gin[numpy.where(gin['gid'] == gid)]
        
        if len(narr) >= 1:
            cdd = 1/(st_fin['plx']*1e-3)      
            nT = test_T5(( app_to_abs(st_fin['gmag'], cdd), st_fin['bprp'], (st_fin['gmag'] - st_fin['Jmag']), (st_fin['W1mag'] - st_fin['W2mag'])), *params1)
            t1.append(narr['FTeff'][0])
            G.append(st_fin['gmag'][0])
            bprp.append(st_fin['bprp'][0])
            t2.append(nT[0])
            
            
    plt.scatter(bprp, G, s = 22, c = t2, cmap = 'jet_r')
    #plt.scatter(t1, t2, s = 22, c = G,  vmin = 17, vmax = 12, cmap = 'plasma')
  #  plt.colorbar()
    
    ax1.set_xlabel('BP - RP [mag]')
    ax1.set_ylabel('G [mag]')
    
    
    cbar = plt.colorbar()#shrink = 0.78, pad = 0.075)
    cbarlabel = 'T$_{fit}$ [K]'
    cbar.ax.set_ylabel(cbarlabel,  labelpad = +10)
    #cbar.ax.invert_yaxis()
    
    plt.savefig('Tfitplots/118_CMD.png', format = 'png', dpi = 600)
    plt.savefig('Tfitplots/118_CMD.pdf', format = 'pdf', dpi = 600)
    return()

#new_CMD(gin, refam)        
 

#Up to here, it's been only concerned with IC5070. Now we go to the rest of the objects 

#fits new temperatures for objects in the periodic source list in all regions
def return_idT2(psl, fam, c11,c12,c21,c22): 
    
    oids = numpy.unique(fam['oid'])
    arr = match_A2(gin, c11, c12, c21, c22) #c11 - c12, c21 - c22
    params1 = fit_Fang_noplot(arr)
    
    for i in range(len(oids)):
        oid = oids[i]
        gid = numpy.unique(fam['gid'][numpy.where(fam['oid'] == oid)])
        
        st_fin = psl[numpy.where(psl['gid'] == gid)]
        
        cdd = 1/(st_fin['plx']*1e-3)       
        
        nT = test_T5((app_to_abs(st_fin['gmag'], cdd), st_fin['bprp'], (st_fin[c11] - st_fin[c12]), (st_fin[c21] - st_fin[c22])), *params1)
        
        fam['Tfit'][numpy.where(fam['oid'] == oid)] = nT
         
    return(fam)



#fam2  = return_idT2(psl, fam, 'bpmag', 'Jmag', 'W1mag', 'W2mag') 
#numpy.save('mypath+'master_famhi8_test.npy', fam2)
 

#plots new temps as a CMD
def plot_cmd_newT(fam, psl):
    
    plt.clf()
    
    rectangle1 = [0.1,0.1, 0.8,0.8]
    ax1 = plt.axes(rectangle1)
    ax1.set_xlim(0.5,3.2)
    ax1.set_ylim(18.2,9.6)
    
    ax1.set_xlabel('BP - RP [mag]')
    ax1.set_ylabel('G [mag]')
        
    Gmag = []
    col = []
    Tfit = []
    
    per = []
    
    oids = numpy.unique(fam['oid'])
 #   arr = match_Av2(gin, c1, c2)
    for i in range(len(oids)):
        oid = oids[i]
        st_fin = psl[numpy.where(psl['oid'] == oid)] 
        st_data = fam[numpy.where(fam['oid'] == oid)]
        
        cdd = 1/(st_fin['plx']*1e-3)   
        Gmag.append(st_fin['gmag'])
        col.append(st_fin['bprp'])
        per.append(st_data['finper'][0])
        Tfit.append((st_data['Tfit'][0]))
   
    Gmag, col, Tfit, per = numpy.array(Gmag), numpy.array(col), numpy.array(Tfit), numpy.array(per)
        

    sort = numpy.argsort(per)[::-1]
    Gmag, col, Tfit, per = Gmag[sort], col[sort], Tfit[sort], per[sort]
    plt.scatter(col, Gmag, s = 10+per**2,  c = Tfit, cmap = 'jet_r') 
    cbar = plt.colorbar()
    cbarlabel = ' $T_{fit}$ '
    cbar.ax.set_ylabel(cbarlabel,  labelpad = +10)
    
    testpers = [2, 4, 8]
    for t in range(len(testpers)):
        plt.scatter(100, 100, s = 10+testpers[t]**2, c = 'black', label = 'P = '+str(testpers[t])+' d')
    plt.legend()
    
    plt.savefig('Tfitplots/fullsample_Tfit.png', format = 'png', dpi = 600)
    plt.savefig('Tfitplots/fullsample_Tfit.pdf', format = 'pdf', dpi = 600)
    return()

#plot_cmd_newT(fam2, psl)

#just an admin function to put known temperatures in the final amp table
def put_Fang_in_Fam(gin, psl, c11, c12, c21, c22):
    
    arr = match_A2(gin, c11, c12, c21, c22) #c11 - c12, c21 - c22
    oids = numpy.unique(fam['oid'])
 #   arr = match_Av2(gin, c1, c2)
    for i in range(len(oids)):
        oid = oids[i]
        st_fin = psl[numpy.where(psl['oid'] == oid)]#
        st_data = fam[numpy.where(fam['oid'] == oid)]
        st_T = arr['FTeff'][numpy.where(arr['gid'] == st_fin['gid'])]
    #    print(st_T)
        if len(st_T) > 0:
            st_data['Teff'] = st_T
            fam[numpy.where(fam['oid'] == oid)] = st_data
            
    print(fam[numpy.where(fam['Teff'] > 100.)])
#    numpy.save('D:/PhD/Codes/stageiii/118_fam_Tfit2.npy', fam)
    
    return()

#put_Fang_in_Fam(gin, psl, c11, c12, c21, c22)

