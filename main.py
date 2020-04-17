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
pulsar_arg=["./pulsar-getter.sh", "233.940149", "10.621519", "1","15","0.8", "45", "60.108976","6.917815","0.85", "15", "20", "-2.4", "refpulsar.gg"]

results = []



def read_pulsar(string): # Reads ASCII, returns dataframe  #"weak.all37.p3fold.ASCII" "W5testmodel.p3fold.ASCII"
  data = ascii.read(string, data_start=1)
  df = data.to_pandas()
  return df



def get_intensities(df, flag):  # Reads dataframe, returns 50x2246 array for plotting OR as a list
  intensities = np.array(df.col4)  # extract intensities column
  pixelarray = np.array(intensities).reshape(50, 2246)  # shape into array with dimensions of image
  croppedarray = pixelarray[:, 1400:2000] #rough onpulse region of exp data

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

  return (chi/(np.var(intensities_img)*(50*600-2)))  #x1*x1/DoF



# Main code starts here

c = 1
min_chi = 10
min_b1 = 0
min_b2 = 0

df_exp = read_pulsar("weak.all37.p3fold.ASCII")
intensities_exp = get_intensities(df_exp, 1)

subprocess.check_output(pulsar_arg)

df_ref = read_pulsar("refpulsar.gg.ASCII")
intensities_ref = get_intensities(df_ref, 1)


chi = fit_measure(intensities_exp, intensities_ref)

print( "Reference pulsar has a fit measure of " + str(chi))


while c < 101:
 # set arguments
 pulsar_number = str(c)
 c+=1

 #a1 = rd.uniform(220, 250)  # 1
 #b1 = rd.uniform(10,11.5)  # 2 10.5
 #c1 = rd.uniform(1, 6)  # 3
 #E = rd.uniform(0.5, 0.95)  # 5
 #osm = rd.uniform(40, 60)  # 6
 #a2 = rd.uniform(50, 80)  # 7
 b2 = rd.uniform(6, 8)  # 8 7.7
 #c2 = rd.uniform(1, 6)  # 9

 #pulsar_arg[1] = str(a1)
 #pulsar_arg[2] = str(b1)
 #pulsar_arg[3] = str(c1)
 #pulsar_arg[5] = str(E)
 #pulsar_arg[6] = str(osm)
 #pulsar_arg[7] = str(a2)
 pulsar_arg[8] = str(b2)
 #pulsar_arg[9] = str(c2)
 pulsar_arg[13] = "SimPulse{}.gg".format(str(pulsar_number))

 subprocess.check_output(pulsar_arg)

 df_sim = read_pulsar("SimPulse"+pulsar_number+".gg.ASCII")
 intensities_sim = get_intensities(df_sim, 1)

 chi = fit_measure(intensities_exp, intensities_sim)
 results.append([b2, chi])


 print( "Pulsar "+ pulsar_number + " has a fit measure of " + str(chi))
 #print("(b1,b2) = ("+str(b1)+", "+str(b2)+")")
 print("b2 = "+ str(b2))


 if chi < min_chi:
   min_chi = chi
   #min_b1 = b1
   min_b2 = b2


 print("current minimum reduced chi squared = " + str(min_chi))
 #print("for (b1,b2) = (" + str(min_b1) + ", " + str(min_b2) + ")")
 print("for b2 = "+str(min_b2))



 #clean up
 os.remove("SimPulse" + pulsar_number + ".gg")
 os.remove("SimPulse"+pulsar_number+".gg.ASCII")


np.savetxt('results_b1b2_b2refine.txt', results, delimiter=',')
