# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 09:57:02 2020

@author: ruthw
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from astropy.io import ascii

#def generate_pulsars(): 

#def read_pulsars():

#def bi_test():

#def plot_pulsars(): 

#def constrain_parameter(): #will likely extend this to individual functions for each paramter or a class


#Read in experimental data
data_exp = ascii.read("weak.all37.p3fold.ASCII", data_start=1)
df_exp = data_exp.to_pandas()

#calculate baseline and add it so that all intensity values are >0

intensities_exp = np.array(df_exp.col4)

df_pixelarray_exp = pd.DataFrame(np.array(intensities_exp).reshape(50,2246)) 

print(df_pixelarray_exp)

plt.imshow(df_pixelarray_exp, 'twilight', origin='lower', interpolation='none', aspect='auto')


#read in simulated data
#data_model = ascii.read("", data_start=1)

#df_model = data_model.to_pandas()

#intensities_model = np.array(df_model.col4)

#df_pixelarray_model = pd.DataFrame(np.array(intensities_model).reshape(50,2246)) 

#plt.imshow(df_pixelarray_model, 'gray', origin='lower', interpolation='none', aspect='auto')


#calculate chi squared

#array = np.array(((df_exp.col4-df_model.col4)^2)/df_model.col4 )
#chi = np.sum(array)
