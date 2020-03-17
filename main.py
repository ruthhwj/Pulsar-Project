# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 09:57:02 2020

@author: ruthw
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import string
from astropy.io import ascii
import glob, os
import sklearn
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import minmax_scale

#pulsar arguments list: fakepulsar -cone "1-9" -ellipse "10-12" -cone'13-21' -ellipse '22-24' -a '25' -b '26' -gg'27'
pulsar_arg=["./pulsar-getter.sh", "0.5", "10.5", "1","15","0.85","45","0.5","7.7","1", "15", "20", "-2.4", "refpulsar.gg"]
#[script name, I1, rho1, w1, n1, e1, or1, s2, I2, rho2, w2, n2,, e2, or2, s2, a, b, filename (include .gg)]
pulsars = {}
def generate_pulsars(inc):
  pulsars["RefPulsar"] = subprocess.check_output(pulsar_arg)
  for j in range(5):
    pulsar_arg[1]=str((float(pulsar_arg[1])+inc))
    pulsar_arg[13]="SimPulseI1{}.gg".format(str(j))
    pulsars["SimPulseI1{}.format(str(j))"] = subprocess.check_output(pulsar_arg)

    
print("attempting to generate model pulsars")
generate_pulsars(0.1)
print("Finished Generating Pulsars")
for x in pulsars:
  print(x.values)



#loop over array formed by read_pulsars()
#def constrain_parameter(): #will likely extend this to individual functions for each paramter or a class

def read_pulsar(string):  # Reads ASCII, returns dataframe  #"weak.all37.p3fold.ASCII" "W5testmodel.p3fold.ASCII"
    data = ascii.read(string, data_start=1)
    df = data.to_pandas()
    return df


def get_intensities(df, flag):  # Reads dataframe, returns 50x2246 array for plotting OR as a list
    intensities = np.array(df.col4)  # extract intensities column    
    if flag == 0: 
	df_pixelarray = pd.DataFrame(np.array(intensities).reshape(50, 2246))# shape into array with dimensions of image
        return df_pixelarray  # want this for plotting
    if flag != 0:
        return intensities  # want this for analysis

print(df_pixelarray_exp)
#read in simulated data
#data_model = ascii.read("", data_start=1)

def plot_pulsar(df_pixelarray):
    plt.imshow(df_pixelarray, 'twilight', origin='lower', interpolation='none', aspect='auto')

#intensities_model = np.array(df_model.col4)

def fit_measure(intensities_ref, intensities_img):
    ref = minmax_scale(intensities_ref)
    img = minmax_scale(intensities_img)

    DoF, chi = (len(ref) - 1), 0

    for i in range(len(ref)):
        x1 = (img[i] - ref[i])
        if img[i] != 0:
            chi += x1 * x1 / img[i]
    return (chi/DoF)


""" 
    #produce chi squared for all pixel intensities = 1
    
    chi_0 = 0
    ones = np.ones(len(ref))
    
    for i in range(len(ref)):
        x1 = (ones[i]-ref[i])
        if img[i] != 0:
        chi_0 += x1*x1/img[i]
 """

# Main code starts here

df_exp = read_pulsar("weak.all37.p3fold.ASCII")
intensities_exp = get_intensities(df_exp, 1)

# loop for each pulsar starts here

for i in whatever:
    df_model = read_pulsar("W5testmodel.p3fold.ASCII")
    intensities_model = get_intensities(df_model, 1)
    chi = fit_measure(intensities_exp, intensities_model)
    #append chi to i of the library or whatever 


