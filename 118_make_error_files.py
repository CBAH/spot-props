"""
Created on Wed Nov 23 11:23:15 2022

@author: TheLa

Generates files with 10000 variations of an objects amplitudes and each variations is compared with the spotgrid to find best fitting spot model, the new amplitudes and the results are put to file

"""


#import os
import numpy
import matplotlib.pyplot as plt
import csv
import datetime
now = datetime.datetime.now()


mypath = 'D:/PhD/Codes/stageiii/lc_118_short/' 
#mypath = '/data3/cbah1/stageiii/lc_118_short/'

print ("Current date and time : ")
print (now.strftime("%Y-%m-%d %H:%M:%S"))


fam = numpy.load('118_finamp.npy') #final amplitude table for region 118. Finamp becomes fam. 

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
                                 ('og_amp_v', float), #original amp
                                 ('og_amp_r', float),
                                 ('og_amp_i', float),
                                 ('og_ampe_v', float),#original amp error
                                 ('og_ampe_r', float),
                                 ('og_ampe_i', float),
                                 ('n_amp_v', float), #new amp error
                                 ('n_amp_r', float),
                                 ('n_amp_i', float),
                                 ('Teff', float),
                                 ('st_model', numpy.unicode_, 16), 
                                 ('sp_model', numpy.unicode_, 16),
                                 ('min_rms_irv', float),
                                 ('sp_size_irv', float),
                                 ('sp_temp_irv', float),
                                 ])
				 


def find_nearest(array, value):
    array = numpy.asarray(array)
    idx = (numpy.abs(array - value)).argmin()
    return array[idx]
                    		


def new_error(mypath, st_data, N):
	st_temp_list = numpy.arange(3500, 5550, 50)
    
	step_temp = find_nearest(st_temp_list, st_data['Teff'])
	g_data = numpy.genfromtxt(mypath+'/spotgrid/'+str(int(step_temp))+'.txt', dtype = amptabletype, delimiter = ',')
	
	g_data = g_data[numpy.where(g_data['sp_model'] == 'ph')]
	g_data = g_data[numpy.where((g_data['sp_temp'] >= 2000) & (g_data['sp_temp'] <= 10000) & (g_data['sp_size'] >= 0.00) & (g_data['sp_size'] <= 0.5))]

	n_amp_v = numpy.random.normal(st_data['finampv'], st_data['finampve'], size = [N])
	n_amp_r = numpy.random.normal(st_data['finampr'], st_data['finampre'], size = [N])
	n_amp_i = numpy.random.normal(st_data['finampi'], st_data['finampie'], size = [N])
	
	n_amp_v = n_amp_v[numpy.where(abs(st_data['finampv'] - n_amp_v) < st_data['finampve'])]
	n_amp_r = n_amp_r[numpy.where(abs(st_data['finampr'] - n_amp_r) < st_data['finampre'])]
	n_amp_i = n_amp_i[numpy.where(abs(st_data['finampi'] - n_amp_i) < st_data['finampie'])]
	
	Nm = numpy.min([len(n_amp_v),len(n_amp_r),len(n_amp_i)])
	if Nm > 10000: Nm = 10000 #ensures exactly 10000 iterations - N is generated, cut to 10000
	
	n_amp_v = n_amp_v[0:Nm]	
	n_amp_r = n_amp_r[0:Nm]	
	n_amp_i = n_amp_i[0:Nm]
	
	#new data	
	n_data = numpy.empty(Nm, rms_tabletype)
		
	n_data['og_amp_v'] = st_data['finampv']
	n_data['og_amp_r'] = st_data['finampr']
	n_data['og_amp_i'] = st_data['finampi']
    
	n_data['og_ampe_v'] = st_data['finampve']
	n_data['og_ampe_r'] = st_data['finampre']
	n_data['og_ampe_i'] = st_data['finampie']
    
	n_data['n_amp_v'] = n_amp_v
	n_data['n_amp_r'] = n_amp_r
	n_data['n_amp_i'] = n_amp_i
	n_data['Teff'] = st_data['Teff']
	n_data['st_model'] = g_data['st_model'][0:Nm]
	n_data['sp_model'] = g_data['sp_model'][0:Nm]
	

	
	g_array_i = numpy.array(g_data['amp_I'])
	g_array = numpy.array([g_data['amp_I'],g_data['amp_R'],g_data['amp_V'] ])
	#print(g_array.shape)
	for j in range(0,Nm):
		found_array = numpy.array( [ g_array_i*0.0 + n_amp_i[j], g_array_i*0.0 + n_amp_r[j], g_array_i*0.0 + n_amp_v[j]]) #get array with new amps, shape of old amps
		rms_list2 = numpy.sqrt(numpy.square(numpy.subtract(found_array, g_array)).mean(axis=0)) 
		sort2 = numpy.argsort(rms_list2)

		data_sort2 = g_data[sort2]
        #single best fit
		n_data['min_rms_irv'][j] = rms_list2[sort2[0]]
		n_data['sp_size_irv'][j] = data_sort2['sp_size'][0]
		n_data['sp_temp_irv'][j] = data_sort2['sp_temp'][0]
       				
	return(n_data)

#put all the iterations in a file and save
def save_as_file(mypath, st_data, n_data, s):
	st_id = st_data['sid'][0]
	with open(mypath+'error_files/'+str(st_id)+'/s'+str(s)+'.txt', 'a+') as file:
		writer = csv.writer(file)
		writer.writerows(n_data)
		file.close()
	return()
	
#loop over all the objects
def loop_new_error(mypath, fam):
    ids = numpy.unique(fam['sid'][numpy.where(fam['hid'] > 99)])
    
    for sid in range(len(ids)):        
        slcs = fam['slice'][numpy.where(fam['sid'] == ids[sid])]
        st_data = fam[numpy.where(fam['sid'] == ids[sid])] 
   
        print(ids[sid])
        for s in (slcs):
             sl_data = st_data[numpy.where(st_data['slice'] == s)]   
             n_data = new_error(sl_data, 20000)
             save_as_file(mypath, st_data, n_data, s)
             
             now = datetime.datetime.now()
             print(now.strftime("%Y-%m-%d %H:%M:%S"))

    return()
	
#loop_new_error(mypath, fam) #To generate error files - the 10000 iterations
