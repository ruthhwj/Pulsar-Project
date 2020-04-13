# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 09:57:02 2020

@author: ruthw
"""

import numpy as np
import pandas as pd
import random as rd
import matplotlib.pyplot as plt
import subprocess
import string
from astropy.io import ascii
import glob, os


#push comment
pulsar_arg=["./pulsar-getter.sh", "0.5", "8.773202", "1","15","0.85","45","0.5","7.7","1", "15", "20", "-2.4", "refpulsar.gg"]

pulsars_args = {}
results = []



def read_pulsar(string): # Reads ASCII, returns dataframe  #"weak.all37.p3fold.ASCII" "W5testmodel.p3fold.ASCII"
  data = ascii.read(string, data_start=1)
  df = data.to_pandas()
  return df



def get_intensities(df, flag):  # Reads dataframe, returns 50x2246 array for plotting OR as a list
  intensities = np.array(df.col4)  # extract intensities column
  pixelarray = np.array(intensities).reshape(50, 2246)  # shape into array with dimensions of image
  croppedarray = pixelarray[:, 1350:2000] #rough onpulse region of exp data

  if flag == 0:
    return croppedarray # want this for plotting
  if flag != 0:
    return croppedarray.flatten()  # want this for analysis



def plot_pulsar(df_pixelarray):
  plt.imshow(df_pixelarray, 'twilight', origin='lower', interpolation='none', aspect='auto')



def fit_measure(intensities_ref, intensities_img):

  chi = 0

  for i in range(len(intensities_ref)):
    x1 = (intensities_img[i] - intensities_ref[i])
    if intensities_img[i] != 0:
      chi += abs(x1 * x1 / intensities_ref[i])
      return (chi)


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

c = 1

df_exp = read_pulsar("norm_exp.ASCII")
intensities_exp = get_intensities(df_exp, 1)

subprocess.check_output(pulsar_arg)

df_ref = read_pulsar("refpulsar.gg.final.ASCII")
intensities_ref = get_intensities(df_ref, 1)


os.remove("refpulsar.gg")
os.remove("refpulsar.gg.D.normalised")
os.remove("refpulsar.gg.final.ASCII")

chi = fit_measure(intensities_exp, intensities_ref)

print( "Reference pulsar has a fit measure of " + str(chi))

results.append(["intensity 1","half opening beam angle 1","beamlets half opening angle 1","eccentricity","orientation of semi major axis","intensity 2","half opening beam angle 2","beamlets half opening angle 2"])
while c < 10001:
 # set arguments
 pulsar_number = str(c)
 c+=1

 pulsar_arg[13] = "SimPulse{}.gg".format(str(pulsar_number))

 a1 = rd.uniform(20, 600)  # 1
 b1 = rd.uniform(7, 11)  # 2
 c1 = rd.uniform(1, 6)  # 3
 E = rd.uniform(0.5, 0.95)  # 5
 osm = rd.uniform(40, 60)  # 6
 a2 = rd.uniform(20, 600)  # 7
 b2 = rd.uniform(4, 11)  # 8
 c2 = rd.uniform(1, 6)  # 9

 pulsar_arg[1] = str(a1)
 pulsar_arg[2] = str(b1)
 pulsar_arg[3] = str(c1)
 pulsar_arg[5] = str(E)
 pulsar_arg[6] = str(osm)
 pulsar_arg[7] = str(a2)
 pulsar_arg[8] = str(b2)
 pulsar_arg[9] = str(c2)


 subprocess.check_output(pulsar_arg)

 df_sim = read_pulsar("SimPulse"+pulsar_number+".gg.final.ASCII")
 intensities_sim = get_intensities(df_sim, 1)

 chi = fit_measure(intensities_exp, intensities_sim)

 print( "Pulsar "+ pulsar_number + " has a fit measure of " + str(chi))

 results.append([a1, b1, c1, E, osm, a2, b2, c2])

 #clean up
 os.remove("SimPulse" + pulsar_number + ".gg")
 os.remove("SimPulse"+pulsar_number+".gg.ASCII")


np.savetxt('results_fullmonte.txt', results, delimiter=',')
