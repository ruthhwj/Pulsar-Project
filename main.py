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
pulsar_arg=["./pulsar-getter.sh", "0.5", "10.5", "1","15","0.85","45","0.5","7.7","1", "15", "20", "-2.4", "refpulsar.gg"]
pulsar_arg_names = ["scriptname", "Cone1Intensity", "Cone1BeamAngle", "Cone1BeamletAngle","Cone1NumberOfSparks", "Eccentricity", "Orientation", "Cone2Intensity",
                    "Cone2BeamAngle", "Cone2BeamletAngle","Cone2NumberOfSparks", "Alpha", "Beta", "Filename"]
pulsar_arg_ranges = [[0.2, 1], [10, 11.5], [0.5, 2], [10, 15], [0, 0.99], [40, 55], [0.2, 1], [7.5,8.1], [0.5,2],[5,15],[0, 45],[0, -45]]
# [cone1intensity, cone1beamangle, cone1beamletangle, cone1numberofsparks, eccentricity, orientation, cone2intensity, cone2beamangle, cone2beamletangle, cone2numberofsparks, alpha, beta
pulsars_args = {}  # pulsar_number : pulsar_arg list


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

def compare_pulsars(pulsar_number):
    results = []
    try:
        df_sim = read_pulsar("SimPulse" + pulsar_number + ".gg.final.ASCII")
    except:
        raise
    intensities_sim = get_intensities(df_sim, 1)
    chi = fit_measure(intensities_exp, intensities_sim)
    print("Pulsar " + pulsar_number + "has a fit measurement of " + str(chi))
    results.append([b1,chi])
    np.savetxt('results_{}.txt'.format(pulsar_arg_names[i]), results, delimiter=',')
    os.remove("SimPulse" + pulsar_number + ".gg")
    os.remove("SimPulse" + pulsar_number + ".gg.D.normalised")
    os.remove("SimPulse" + pulsar_number + ".gg.final.ASCII")

# Main code starts here
def main():
    df_exp = read_pulsar("norm_exp.ASCII")
    intensities_exp = get_intensities(df_exp, 1)
    for i in range(4,13):
        c=1
        while c <= 20:
            try:
                pulsar_number = str(c)
                c+=1
             #pulsar_arg[1] = str((param_dict[1][i]))
                pulsar_arg[13] = "SimPulse{}.gg".format(str(pulsar_number))
                b1 = np.random.uniform(pulsar_arg_ranges[i-1][0], pulsar_arg_ranges[i-1][1])
                print(b1)
                pulsar_arg[4] = str(14)
                pulsar_arg[i] ='{0:.4f}'.format(float(str(b1)))
                x = pulsar_arg[i]
                print(pulsar_arg[i])
                subprocess.check_output(pulsar_arg)
                df_sim = read_pulsar("SimPulse" + pulsar_number + ".gg.final.ASCII")
            except:
                try:
                    print("Exception Called, testing for limits error")
                    pulsar_arg[i]=str(float(pulsar_arg[i])+0.0001)
                    subprocess.check_output(pulsar_arg)
                    df_sim = read_pulsar("SimPulse" + pulsar_number + ".gg.final.ASCII")
                except:
                    print("Type conversion error, rounding to nearest integer")
                    pulsar_arg[i] = int(float(pulsar_arg[i]))
                    x=pulsar_arg[i]
                    subprocess.check_output(pulsar_arg[i])
                    compare_pulsars(pulsar_number)

         # set arguments


if __name__ == '__main__':
    main()
