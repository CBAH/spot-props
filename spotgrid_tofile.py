
"""

Created on Tue Jan 26 16:34:17 2021



@author: TheLa

THIS IS THE SCRIPT THAT GENERATES A MAGNITUDE DIFFERENCE FOR RANDOM SPOT TEMP, SPOT TEMP AND STAR TEMP.

2024 edit

"""



import os

import numpy

#import all_hc_functions as hc

import matplotlib.pyplot as plt

import pysynphot as S

import speclite.filters

import math

from scipy import stats



import random

import csv

import datetime

import sys




now = datetime.datetime.now()

print ("Current date and time : ")

print (now.strftime("%Y-%m-%d %H:%M:%S"))



def tablemakery(st_temp, met, logg, st_model, sp_model):


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

    ('amp_I',float),   #mag in I


    ])


    calitable = numpy.load('vegamag_cali.npy') #calibration table with constants to calibrate to vega - made in model_differences/'modelfunctions.py'


    counter = 0

    bessell = speclite.filters.load_filters('bessell-*')

    amptable= numpy.empty(1, dtype = amptabletype)

    fa = random.uniform(0.00, 0.500) 
    
    st_spec = S.Icat('phoenix', st_temp, met, logg)

    sp_temp = random.uniform(2000, 12000)

    sp_spec = S.Icat('phoenix', sp_temp, met, logg)
	
    calitable = calitable[2] #[0] for castelli, [2] for phoenix, [1] is black body

    sys_spec = (st_spec*(1.0-fa))+(sp_spec*fa)


    st_mag = bessell.get_ab_magnitudes(st_spec.flux,st_spec.wave).as_array()  #get PEAK-TO-PEAK amplitudes

    sys_mag = bessell.get_ab_magnitudes(sys_spec.flux, sys_spec.wave).as_array()

    

    amp_U = abs((st_mag['bessell-U']*calitable['U'])-(sys_mag['bessell-U']*calitable['U']))

    amp_B = abs((st_mag['bessell-B']*calitable['B'])-(sys_mag['bessell-B']*calitable['B']))
     
    amp_V = abs((st_mag['bessell-V']*calitable['V'])-(sys_mag['bessell-V']*calitable['V']))

    amp_R = abs((st_mag['bessell-R']*calitable['R'])-(sys_mag['bessell-R']*calitable['R']))

    amp_I = abs((st_mag['bessell-I']*calitable['I'])-(sys_mag['bessell-I']*calitable['I']))

               
    amptable['number'][counter] = counter

    amptable['st_model'][counter] = st_model

    amptable['sp_model'][counter] = sp_model

    amptable['st_temp'][counter] = st_temp

    amptable['sp_temp'][counter] = sp_temp

    amptable['sp_size'][counter] = fa

    amptable['met'][counter] = met

    amptable['logg'][counter] = logg

    amptable['amp_U'][counter] = amp_U

    amptable['amp_B'][counter] = amp_B

    amptable['amp_V'][counter] = amp_V

    amptable['amp_R'][counter] = amp_R

    amptable['amp_I'][counter] = amp_I

    counter = counter + 1 


    return(amptable)

    
#tablemakery(0,4)


def loopcount():

    loopcount = 0 

    while loopcount < 1e20:

#
        st_temp_list = numpy.arange(3500, 5550, 50)
        
        #st_temp_list = numpy.arange(5550, 8050, 50) #for hot

        st_temp = st_temp_list[random.randint(0, (len(st_temp_list) - 1))]

        data2 = tablemakery(st_temp, 0, 4, 'ph', 'ph')

        with open('step_temp_files_2024/'+str(int(st_temp))+'.txt', 'a+') as file:


           writer = csv.writer(file)
##
           writer.writerows(data2)

           if loopcount % 1000 == 0 :

               print(loopcount, 'cool')

               now = datetime.datetime.now()

               print (now.strftime('%H:%M:%S'))

           file.close()

           loopcount = loopcount + 1

    return()


loopcount()


def test_length():
    
    gdata = numpy.genfromtxt('step_temp_files_2024/4000.txt')
    print(len(gdata))
    return()


#test_length()

#to test format
def single():

    st_temp, data1 = tablemakery( 0, 4, 'ph', 'ph')
    

    print(data1)



    return()



#single()



now = datetime.datetime.now()

print ("Current date and time : ")

print (now.strftime("%Y-%m-%d %H:%M:%S"))
