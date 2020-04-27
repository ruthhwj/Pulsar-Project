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
pulsar_arg=["./pulsar-getter.sh", "233.940149", "10.5", "1","15","0.8", "45", "60.108976","7.7","0.85", "15", "20", "-2.4", "refpulsar.gg", "17", "29"]

results = []



def read_pulsar(string): # Reads ASCII, returns dataframe  #"weak.all37.p3fold.ASCII" "W5testmodel.p3fold.ASCII"
  data = ascii.read(string, data_start=1)
  df = data.to_pandas()
  return df



def get_intensities(df, flag):  # Reads dataframe, returns 50x2246 array for plotting OR as a list
  intensities = np.array(df.col4)  # extract intensities column
  pixelarray = np.array(intensities).reshape(50, 1123)  # shape into array with dimensions of image
  croppedarray = pixelarray[:, 700:1000] #rough onpulse region of exp data

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
    chi += abs(x1*x1)

  return (chi/(RMS_noise*(50*300-2)))  #x1*x1/noise*DoF

def gaussian(x, mu, sig):
  return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))



c = 1
min_chi = 10000
min_a1 = 0
min_a2 = 0
i=0
j=0
param = 1

# Main code starts here #

df_exp = read_pulsar("weak.all37.p3fold.rebinned.ASCII")
intensities_exp = get_intensities(df_exp, 0)

# artifically brighten left component using a gaussian centred on the brightest column #

while i < 50:
 while j < 300:
  gauss_j = (gaussian(float(j), 85, 25)) #centred on bin 85, width 25
  intensities_exp[i][j] = ((param*gauss_j+1)*intensities_exp[i][j])
  j += 1
 i+=1
 j=0

intensities_exp_flat = intensities_exp.flatten() #1d list used for analysis

# find RMS noise for reduced chi squared #

intensities_RMS = np.array(df_exp.col4)
exp_croppedlist = ((intensities_RMS.reshape(50, 1123))[:, 600:700]).flatten()  #off pulse RMS noise
RMS_noise = np.var(exp_croppedlist)

# process reference pulsar #

subprocess.check_output(pulsar_arg)

df_ref = read_pulsar("refpulsar.gg.ASCII")
intensities_ref = get_intensities(df_ref, 1)

chi = fit_measure(intensities_exp_flat, intensities_ref)

print( "Reference pulsar has a fit measure of " + str(chi))


while c < 10001:
 c+=1

 # set arguments

 pulsar_number = str(c)
 a1 = rd.uniform(230, 250)  # 1
 b1 = rd.uniform(9, 12)  # 2
 c1 = rd.uniform(1, 2)  # 3
 E = [0.7, 0.72, 0.74, 0.76, 0.78, 0.80, 0.82, 0.84, 0.86, 0.88, 0.90] # avoid weird floating point error
 osm = [43, 44, 45, 46, 47]
 a2 = rd.uniform(40, 80)  # 7
 b2 = rd.uniform(4, 8)  # 8
 c2 = rd.uniform(0.5, 1.5)  # 9
 d1 = rd.uniform(14, 20)
 d2 = rd.uniform(26, 32)

 E_choice = round(rd.choice(E), 2)
 osm_choice = round(rd.choice(osm), 0)



 pulsar_arg[1] = str(a1)
 pulsar_arg[2] = str(b1)
 pulsar_arg[3] = str(c1)
 pulsar_arg[5] = str(E_choice)
 pulsar_arg[6] = str(osm_choice)
 pulsar_arg[7] = str(a2)
 pulsar_arg[8] = str(b2)
 pulsar_arg[9] = str(c2)
 pulsar_arg[13] = "SimPulse{}.gg".format(str(pulsar_number))
 pulsar_arg[14] = str(d1)
 pulsar_arg[15] = str(d2)


 subprocess.check_output(pulsar_arg)

 df_sim = read_pulsar("SimPulse"+pulsar_number+".gg.ASCII")
 intensities_sim = get_intensities(df_sim, 1)

 chi = fit_measure(intensities_exp_flat, intensities_sim)
 results.append([a1, a2, b1, b2, c1, c2, d1, d2, E_choice, osm_choice, chi])


 print( "Pulsar "+ pulsar_number + " has a fit measure of " + str(chi))
 #print("(b1,b2) = ("+str(b1)+", "+str(b2)+")")
 #print("(a1,a2) = ("+ str(a1)+","+str(a2)+")")


 if chi < min_chi:
   min_chi = chi


 print("current minimum reduced chi squared = " + str(min_chi))
 #print("for (a1,a2) = (" + str(min_a1) + ", " + str(min_a2) + ")")
 #print("for b2 = "+str(min_b2))



 #clean up
 os.remove("SimPulse" + pulsar_number + ".gg")
 os.remove("SimPulse"+pulsar_number+".gg.ASCII")

 #save results file every 500 runs

 if (c/500).is_integer():
   np.savetxt('results_full_270420_2.txt', results, delimiter=',')



