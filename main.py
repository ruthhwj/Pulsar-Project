# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 09:57:02 2020

@author: ruthw
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
from astropy.io import ascii
import glob, os
import sklearn 
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import minmax_scale

def generate_pulsars(i, arg1, *argv): 
#pulsar arguments list: fakepulsar -cone "1-9" -ellipse "10-12" -cone'13-21' -ellipse '22-24' -a '25' -b '26' -gg'27'
pulsar_arg=["./pulsar-getter.sh", "1", "10.5", "1","15","750","12","-1","300","0","0.85","45","1","0.5","7.7","1","15","750","12","-1","300","0","0.8","45","1","20","-2.4", "w6test3model.gg"]
#[script name, I1, rho1, w1, n1, p4_1, phi0_1, fP4_1, t4_1, phi4_1 e1, or1, s2, I2, rho2, w2, n2, p4_2, phi0_2, fP4_2, t4_2, phi4_2 e2, or2, s2, a, b, filename (include .gg)]
test = subprocess.Popen(pulsar_arg)


#    for arg in argv:
#        int_arg = int(pulsar_arg[i]
#    subprocess.check_call(pulsar_arg)
#call bash script looping over pulsars
#def read_pulsars():
#    global data=[]
#    for file in glob.glob("*.ASCII"):
#        data.append(file)
#read in pulsars generated by generate_pulsars
#store as an array

#def bi_test():
#def plot_pulsars():
#loop over array formed by read_pulsars()
#def constrain_parameter(): #will likely extend this to individual functions for each paramter or a class

#Read in experimental data
data_exp = ascii.read("weak.all37.p3fold.ASCII", data_start=1)

df_exp = data_exp.to_pandas()

intensities_exp = np.array(df_exp.col4)

df_pixelarray_exp = pd.DataFrame(np.array(intensities_exp).reshape(50,2246)) 

#plt.imshow(df_pixelarray_exp, 'twilight', origin='lower', interpolation='none', aspect='auto')

#read in simulated data
data_model = ascii.read("W5testmodel.p3fold.ASCII", data_start=1)

df_model = data_model.to_pandas()

intensities_model = np.array(df_model.col4)

df_pixelarray_model = pd.DataFrame(np.array(intensities_model).reshape(50,2246)) 

#plt.imshow(df_pixelarray_model, 'twilight', origin='lower', interpolation='none', aspect='auto')





#replace experimental intensities with the smallest non zero value
ref = minmax_scale(intensities_exp)
img = minmax_scale(intensities_model)

DoF = (len(ref) -1)
chi=0
chi_0=0

ones = np.ones(len(ref))

for i in range(len(ref)):
    x1 = (ones[i]-ref[i])
    if img[i] != 0:
        chi_0 += x1*x1/img[i]

print(chi_0/DoF) 


for i in range(len(ref)):
    x1 = (img[i]-ref[i])
    if img[i] != 0:
        chi_i = x1*x1/img[i]
        #print(chi_i)
        chi += chi_i

print(chi/DoF)