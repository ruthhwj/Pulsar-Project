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

#pulsar arguments: fakepulsar -cone "1-9" -ellipse "10-12" -cone'13-21' -ellipse '22-24' -a '25' -b '26' -gg'27'

pulsar_arg=["pulsar_getter.sh", "1", "10.5", "1","15","750","12","-1","300","0","0.85","45","1","0.5","7.7","1","15","750","12","-1","300","0","0.8","45","1","20","-2.4", "w6test3model.gg"]
test = subprocess.call(pulsar_arg)
output = test.communicate()[0]
print(output)
#[script name, I1, rho1, w1, n1, p4_1, phi0_1, fP4_1, t4_1, phi4_1 e1, or1, s2, I2, rho2, w2, n2, p4_2, phi0_2, fP4_2, t4_2, phi4_2 e2, or2, s2, a, b, filename (include .gg)]
# def

#def generate_pulsars(i, arg1, *argv): #arg1 is the starting parameter, argv is the total list of parameters
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
"""
def plot_pulsars():
#loop over array formed by read_pulsars()
#def constrain_parameter(): #will likely extend this to individual functions for each paramter or a class


#Read in experimental data
data_exp = ascii.read("weak.all37.p3fold.ASCII", data_start=1)
df_exp = data_exp.to_pandas()

#calculate baseline and add it so that all intensity values are >0

intensities_exp = np.array(df_exp.col4)

df_pixelarray_exp = pd.DataFrame(np.array(intensities_exp).reshape(50,2246)) 

print(df_pixelarray_exp)

plt.imshow(df_pixelarray_exp, 'twilight', origin='lower', interpolation='none', aspect='auto')

print(df_pixelarray_exp)
#read in simulated data
#data_model = ascii.read("", data_start=1)

#df_model = data_model.to_pandas()

#intensities_model = np.array(df_model.col4)

#df_pixelarray_model = pd.DataFrame(np.array(intensities_model).reshape(50,2246)) 

#plt.imshow(df_pixelarray_model, 'gray', origin='lower', interpolation='none', aspect='auto')


#calculate chi squared

testarray=np.array(50*2246)
print(test)

#calculate chi squared 

x1 = intensities_exp - testarray
x2 = np.divide(x1*x1, intensities_exp)
chi = np.sum(x2)

"""
