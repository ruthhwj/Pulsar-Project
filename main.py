# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 09:57:02 2020

@author: ruthw
"""

import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
import subprocess
import string
from astropy.io import ascii
import glob, os
import sklearn
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import minmax_scale

#push comment
pulsar_arg=["./pulsar-getter.sh", "0.5", "10.5", "1","15","0.85","45","0.5","7.7","1", "15", "20", "-2.4", "refpulsar.gg"]

pulsars_args = {}  # pulsar_number : pulsar_arg list
results = []

param_dict = {
  1 : list(np.arange(0.3, 0.8, 0.1)),  # intensity
  2 : list(np.arange(10, 11.5, 0.5)),  # half opening angle of beam
  3 : list(np.arange(0.5, 2, 0.5)),  # half opening angle of beamlets
  4 : list(np.arange(14, 17, 1)),  # number of sparks
  5 : list(np.arange(0.8, 0.95, 0.05)),  # eccentricity
  6 : list(np.arange(40, 55, 5)),  # orientation of semi major axis
  7 : list(np.arange(0.4, 0.7, 0.1)),  # intensity
  8 : list(np.arange(7.5, 8.1, 0.2)),  # half opening angle of beam
  9 : list(np.arange(0.5, 2, 0.5))  # half opening angle of beamlets
}

def read_pulsar(string): # Reads ASCII, returns dataframe  #"weak.all37.p3fold.ASCII" "W5testmodel.p3fold.ASCII"
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



def plot_pulsar(df_pixelarray):
  plt.imshow(df_pixelarray, 'twilight', origin='lower', interpolation='none', aspect='auto')



def fit_measure(intensities_ref, intensities_img):

  DoF, chi = (len(intensities_ref) - 1), 0

  for i in range(len(intensities_ref)):
    x1 = (intensities_img[i] - intensities_ref[i])
    if intensities_img[i] != 0:
      chi += x1 * x1 / intensities_img[i]
      return (chi/(DoF))


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



#noise_array = np.array(intensities_exp).reshape(50, 2246)
# find average of white noise in off pulse region, col 0->1450 and 1950->2246
#x1 = (noise_array[:, 0:1450])
#x2 = (noise_array[:, 1950:2246])

#noise = (1450 / (1450 + 296)) * np.mean(x1) + (296 / (1450 + 296)) * np.mean(x2)
#print(noise)

c = 1

df_exp = read_pulsar("norm_exp.ASCII")
intensities_exp = get_intensities(df_exp, 1)


for i in range(len(param_dict[1])):
 # set arguments
 pulsar_number = str(c)
 c+=1

 pulsar_arg[1] = str((param_dict[1][i]))
 pulsar_arg[13] = "SimPulse{}.gg".format(str(pulsar_number))

 #pulsar_arg.pop(14) # weird 14th argument showing up, idk why just get rid


 subprocess.check_output(pulsar_arg)

 df_sim = read_pulsar("SimPulse"+pulsar_number+".gg.final.ASCII")
 intensities_sim = get_intensities(df_sim, 1)

 chi = fit_measure(intensities_exp, intensities_sim)

 print( "Pulsar "+ pulsar_number + " has a chi squared of " + str(chi))
 results.append([(param_dict[1][i]), chi])

#clean up
 os.remove("SimPulse" + pulsar_number + ".gg")
 os.remove("SimPulse"+pulsar_number+".gg.D.normalised")
 os.remove("SimPulse"+pulsar_number+".gg.noise")
 os.remove("SimPulse"+pulsar_number+".gg.final.ASCII")



np.savetxt('results.txt', results, delimiter=',')
np.save("results", results)