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

pulsars_args = {}  # pulsar_number : pulsar_arg list

param_dict = {
  1 : list(np.arange(0.4, 0.7, 0.1)),  # intensity
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
  ref = minmax_scale(intensities_ref) # concerned about validity of doing this
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

subprocess.check_output(pulsar_arg)  # reference pulsar
pulsars_args[0] = pulsar_arg

df_sim = read_pulsar("refpulsar.p3fold.ASCII")
intensities_sim = get_intensities(df_sim, 1)
chi = fit_measure(intensities_exp, intensities_sim)
pulsars_args[0][14] = chi

for i in range(len(param_dict[1])):
  for j in range(len(param_dict[2])):
    for k in range(len(param_dict[3])):
      for l in range(len(param_dict[4])):
        for m in range(len(param_dict[5])):
          for n in range(len(param_dict[6])):
            for o in range(len(param_dict[7])):
              for p in range(len(param_dict[8])):
                for q in range(len(param_dict[9])):

                  # set arguments

                  pulsar_number = str(i + j + k + l + m + n + o + p + q + 1)

                  pulsar_arg[1] = str(param_dict[1][i])
                  pulsar_arg[2] = str(param_dict[2][j])
                  pulsar_arg[3] = str(param_dict[3][k])
                  pulsar_arg[4] = str(param_dict[4][l])
                  pulsar_arg[5] = str(param_dict[5][m])
                  pulsar_arg[6] = str(param_dict[6][n])
                  pulsar_arg[7] = str(param_dict[7][o])
                  pulsar_arg[8] = str(param_dict[8][p])
                  pulsar_arg[9] = str(param_dict[9][q])

                  pulsar_arg[13] = "SimPulse{}.gg".format(pulsar_number)

                  subprocess.check_output(pulsar_arg)  # run fake Pulsar

                  pulsars_args[int(pulsar_number)] = pulsar_arg # dump arg list into dictionary

                  df_sim = read_pulsar("SimPulse"+pulsar_number+".p3fold.ASCII")
                  intensities_sim = get_intensities(df_sim, 1)
                  chi = fit_measure(intensities_exp, intensities_sim)

                  pulsars_args[int(pulsar_number)].append(chi) # append chi onto list in dictionary


                  print( "Pulsar "+ pulsar_number + " has a chi squared of " + chi)

                  #clean up

                  os.remove("SimPulse"+pulsar_number+".p3fold.ASCII")
                  os.remove("SimPulse"+pulsar_number+".gg")
