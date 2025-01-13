import sys
import numpy
import math
import csv
import JDmod
import all_hc_functions as hc

mypath = 'D:/PhD/Codes/stageiii/lc_118/'
numpy.set_printoptions(threshold=sys.maxsize)

'''
File to slice light curves in HOYS format into 6month slices. New lightcurves are saves in subfolders, numbered s3 as 's3'. Folders and subfolders structure are made in 'generate_lists_folders.py' 
Lookup folders/files are made as light curves are processed.

Iteration -- 5
'''

#using 02/14/2014 as start date - unfavourable position in the sky -less likely to cut data

def make_boundary_slices():
    deltat = 91.3105	 #1/4 year
    nslices = math.floor((JDmod.get_current_JD() - JDmod.turn_date_JD('2014-02-14 12:00:00'))/deltat) +1 
    bdates = [JDmod.turn_date_JD('2014-02-14 12:00:00')] #boundarydates! boundaries starting at 2456703.83333. 2014-02-15T0800
    for n in range(nslices):
        bdates.append(int((JDmod.turn_date_JD('2014-02-14 12:00:00')+deltat*(n+1)))) #integer part to make changeover midday 
    return(bdates) 


def slice_to_mjd(sls):
    dates = []
    if type(sls) == int: dates = bdates[sls] - 2400000.5 #turns julian date to modified jd
    else:
        for i in sls:
            dates.append(bdates[i] -2400000.5)
    return(dates) 

bdates = make_boundary_slices()

'''
for i in range(len(bdates)):
    print(i, bdates[i], JDmod.turn_date_str(bdates[i]))
'''

ids = numpy.genfromtxt('lc_118_ids.txt', dtype = 'int64') #gaia ids too long for int - use int64


def slice_lc_period(bdates, mypath, ids):
    lk_table = numpy.load(mypath+'s_slc/lookup/'+'lk_table.npy') #lookup table is generated, but empty (full of 99s). 
    
    for i in range(0,len(ids)):
        sid = ids[i]
        lk_sid = lk_table[numpy.where(lk_table['sid'] == sid)]
       
        lc = numpy.load(mypath+'lc_'+str(ids[i])+'.npy') #lightcurve
        
        lkl = ['Object '+str(sid)] #lookup_list goes to write in the file. 
        zerofit_counter = 0
        onefit_counter = 0
        twofit_counter = 0
        threefit_counter = 0
        for j in range(len(bdates) - 2):
            clc = lc[numpy.where((lc['date'] >= bdates[j]) & (lc['date'] <= bdates[j+2]))] #between start date and up to date +2. if deltat is 1/4 year, slice is 1/2 year

            fff = ['V','R','I']
            filt_num = hc.check_filter_single(clc,filt_list=fff) #checks how many points of each filter there are
            len_check = filt_num[numpy.where(filt_num > 50.)]
            
            
            lk_sid['s_count'][numpy.where(lk_sid['slice'] == j)] = numpy.max(filt_num)
            lk_sid['f_count'][numpy.where(lk_sid['slice'] == j)] = len(len_check)
            if len(len_check) < 2:
                
                if len(len_check) == 0: zerofit_counter = zerofit_counter + 1
                if len(len_check) == 1: onefit_counter = onefit_counter + 1
                
                fff = numpy.array(fff)
                fff = fff[numpy.where((filt_num > 50.))]
                lk_sid['s_filters'][numpy.where(lk_sid['slice'] == j)] = str(fff)
                
            #    print('there are less than 50 points in at least 2 filters') 
                lkl.append('Slice '+str( j )+' with int boundary start date '+ str(int(bdates[j]))+': < 50 points in at least 2 filters')
            
            else:
                fff = numpy.array(fff)
                fff = fff[numpy.where((filt_num > 50.))]
                
                clc = clc[numpy.where(clc['filter'] )]###
                
                lk_sid['s_filters'][numpy.where(lk_sid['slice'] == j)] = str(fff)
                
                if len(len_check) == 2: twofit_counter = twofit_counter + 1
                if len(len_check) == 3: threefit_counter = threefit_counter + 1
                
                lkl.append('slice '+str( j )+' with int boundary start date '+ str(int(bdates[j]))+': > 50 points in filters' + str(fff))                
                numpy.save(mypath+'s_slc/s'+str(j)+'/'+str(sid)+'.npy', clc) #saves shorter light curve with slice s1, with the object id.
			
            
            #### Uncomment this to generate text files to use with Fenia's code, generates text files to input into L2Beta_FEnia.R.
            
            
               # for f in fff:
               #     clcf = clc[numpy.where(clc['filter'] == f)]
               #     Lf = numpy.empty([len(clcf), 3])
		
				#Fenia's code uses format mag, mage, date.
				
                   # Lf[:,2] = clcf['date']
                   # Lf[:,0] = clcf['calibrated_magnitude']
                   # Lf[:,1] = clcf['calibrated_error']
				
#                    with open(mypath+'s_slc/s'+str(j)+'/'+str(sid)+'_'+str(f)+'.txt', 'w+') as file:  
					
                       # file.write(str(Lf[0:2]))
#                        writer = csv.writer(file, delimiter = ' ')
#                        writer.writerows(Lf)
#                        file.close()
                        
        #lookup folder - for each object there will be a text summarising each slice
        with open(mypath+'s_slc/lookup/object_lookup/'+str(sid)+'.txt', 'w+') as file:
            for i in range(len(lkl)):
                file.write(lkl[i] + '\n')
            file.close()
            
        with open(mypath+'s_slc/lookup/object_lookup/'+str(sid)+'.txt', 'a+') as file: #generates text files to put in lookup folder.

            file.write('There are '+str(zerofit_counter)+' slices with 0 filters'+'\n')
            file.write('There are '+str(onefit_counter)+' slices with 1 filters'+'\n')
            file.write('There are '+str(twofit_counter)+' slices with 2 filters'+'\n')
            file.write('There are '+str(threefit_counter)+' slices with 3 filters'+'\n')
            file.close()
    
        lk_table[numpy.where(lk_table['sid'] == sid)]  = lk_sid
        
        numpy.save(mypath+'s_slc/lookup/'+'lk_table.npy', lk_table)
        
    return()

#slice_lc_period(bdates, mypath, ids)
