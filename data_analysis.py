import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import string
from astropy.io import ascii
import glob, os
import sklearn
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cmx
import scipy as sp
from mpl_toolkits import mplot3d
import numpy as np


pulsar=["./pulsar-getter.sh", "233.940149", "10.5", "1", "15" , "17", "0.85", "45", "60", "7.7", "1", "15", "29", "refpulsar.gg"]
pulsar_arg_names = ["scriptname", "Cone1Intensity", "Cone1BeamAngle", "Cone1BeamletAngle","Cone1NumberOfSparks", "Cone1phi0", "Eccentricity", "Orientation", "Cone2Intensity",
                    "Cone2BeamAngle", "Cone2BeamletAngle","Cone2NumberOfSparks", "Cone2phi0", "Filename"]
pulsar_arg_ranges = [[230, 260], [8, 12], [1, 2], [15, 15], [10,20] , [0.5, 0.9], [40, 50], [40, 100], [6,10], [0.5,1.5], [15,15], [22,32]] #ranges over which to search for each variable

"""
1 cone 1 intensity
2 cone 1 beam angle
3 cone 1 beamlet angle
4 cone 1 number of sparks
5 cone 1 offset
6 eccentricity
7 orientation of semi major
8 cone 2 intensity
9 cone 2 beam angle
10 cone 2 beamlet angle
11 cone 2 number of sparks
12 cone 2 offset
"""

def read_pulsar(string):  # Reads ASCII, returns dataframe  #"weak.all37.p3fold.ASCII"
  data = ascii.read(string, data_start=1)
  df = data.to_pandas()
  return df

def get_intensities(df, flag):  # Reads dataframe, returns 50x2246 array for plotting OR as a list
  intensities = np.array(df.col4)  # extract intensities column
  pixelarray = np.array(intensities).reshape(50, 1123)  # shape into array with dimensions of image
  croppedarray = pixelarray[:, 700:1000]  # rough onpulse region of exp data

  if flag == 0:
    return croppedarray  # want this for plotting
  if flag != 0:
    return croppedarray.flatten()  # want this for analysis

def plot_pulsar(df_pixelarray):
 print("Plotting pulsar...")
 plt.imshow(df_pixelarray, 'afmhot', origin='lower', interpolation='none', aspect='auto',extent=[700,1000,0,50])
 plt.colorbar()
 plt.show()

def gaussian(x, mu, sig):
  return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def brighten(exp_data):
  i = 0
  j = 0
  param = 1
  while i < 50:
    while j < 300:
      gauss_j = (gaussian(float(j), 85, 25))  # centred on bin 85, width 25
      exp_data[i][j] = ((param * gauss_j + 1) * exp_data[i][j])
      j += 1
    i += 1
    j = 0
  return exp_data


pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


      #main code#

#data preprocessing

#df_exp = read_pulsar("weak.all37.p3fold.rebinned.ASCII")  # experimental p3fold here
#intensities_exp = brighten(get_intensities(df_exp, 0)).flatten()

#

# N=1000

for i in [x for x in range(1,13) if (x!=4 and x!=11 and x!=1 and x!=8)]:
  data1d = pd.read_csv(r'Results\1k 1D\results{}.txt'.format(pulsar_arg_names[i]), sep=",", header=None)
  data1d.columns = ["col1", "chi"]
  print(data1d.nsmallest(5, 'chi'))
  plt.scatter(data1d.col1, data1d.chi, linewidth=1)
  plt.xlabel(pulsar_arg_names[i])
  plt.ylabel('Reduced Chi Squared')
  plt.title("N=1000")
  plt.show()

# N=5000

for i in [x for x in range(1,13) if (x!=4 and x!=11 and x!=1 and x!=8)]:
  data1d = pd.read_csv(r'Results\5k 1D\results{}.txt'.format(pulsar_arg_names[i]), sep=",", header=None)
  data1d.columns = ["col1", "chi"]
  print(data1d.nsmallest(5, 'chi'))
  plt.scatter(data1d.col1, data1d.chi, linewidth=1)
  plt.xlabel(pulsar_arg_names[i])
  plt.ylabel('Reduced Chi Squared')
  plt.title("N=5000")
  plt.show()



# N=10000

for i in [x for x in range(1,13) if (x!=4 and x!=11 and x!=1 and x!=8)]:
  data1d = pd.read_csv(r'Results\10k 1D\results{}.txt'.format(pulsar_arg_names[i]), sep=",", header=None)
  data1d.columns = ["col1", "chi"]
  print(data1d.nsmallest(5, 'chi'))
  plt.scatter(data1d.col1, data1d.chi, linewidth=1)
  plt.xlabel(pulsar_arg_names[i])
  plt.ylabel('Reduced Chi Squared')
  plt.title("N=10000")
  plt.show()


# N=50000

for i in [x for x in range(1,13) if (x!=4 and x!=11 and x!=1 and x!=8)]:
  data1d = pd.read_csv(r'Results\500 1D\results{}.txt'.format(pulsar_arg_names[i]), sep=",", header=None)
  data1d.columns = ["col1", "chi"]
  print(data1d.nsmallest(5, 'chi'))
  plt.scatter(data1d.col1, data1d.chi, linewidth=1)
  plt.xlabel(pulsar_arg_names[i])
  plt.ylabel('Reduced Chi Squared')
  plt.title("N=50000")
  plt.show()







