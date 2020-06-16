import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import ascii
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
 plt.xlabel("Pulse bin")
 plt.ylabel('Pulse Subint')
 #plt.colorbar()
 plt.show()

def gaussian(x, mu, sig):
  return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def fit_measure(intensities_ref, intensities_img):
    chi = 0
    for i in range(len(intensities_ref)):
        x1 = (intensities_img[i] - intensities_ref[i])
        chi += abs(x1 * x1)/ (RMS_noise * (50 * 300 - 2))
    return chi   #


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



df_exp = read_pulsar("weak.all37.p3fold.rebinned.ASCII")  # experimental p3fold here
intensities_exp = brighten(get_intensities(df_exp, 0)).flatten()
intensities_RMS = np.array(df_exp.col4)
exp_croppedlist = ((intensities_RMS.reshape(50, 1123))[:, 0:600]).flatten()  # off pulse RMS noise
RMS_noise = np.var(exp_croppedlist)
print(RMS_noise)


#attempt to get individual chis out for chosen ASCII files

"""
df_ref = read_pulsar("realreferencepulsar.ASCII")
intensities_ref = get_intensities(df_exp,1)
fitmeasure = fit_measure(intensities_ref, intensities_exp)
print(fitmeasure)
plot_pulsar(get_intensities(df_ref,0))
"""




# 1D plots
N=[100,500,1000]

for j in N:
  for i in [x for x in range(1,13) if (x!=4 and x!=11)]:
    data1d = pd.read_csv(r'Results\N = {} 1D\results{}.txt'.format(j,pulsar_arg_names[i]), sep=",", header=None)
    data1d.columns = ["col1", "chi"]
    print(data1d.nsmallest(5, 'chi'))
    plt.scatter(data1d.col1, data1d.chi, linewidth=1)
    plt.xlabel(pulsar_arg_names[i])
    plt.ylabel('Reduced Chi Squared')
    plt.title("N="+str(j))
    plt.savefig('1D_graphs\AllResults_N{}_{}.png'.format(j, pulsar_arg_names[i]))
    plt.show()

"""

# ND plots

N = 10000  #change this depending on which dataset you want to look at
col = []
dataND = pd.read_csv('AllResults_N{}.csv'.format(N), sep=",", header=None) #has 11 columns

for i in [x for x in range(1,13) if (x!=4 and x!=11)]:
  col.append(pulsar_arg_names[i])
col.append("chi")
print(col)
dataND.columns = col
print(dataND.nsmallest(20, 'chi')) # outputs args of the 20 minimum chi pulsars
for column in dataND:
  contents = dataND[column]
  plt.scatter(contents, dataND.chi, linewidth=1)
  plt.xlabel(column)
  plt.ylabel('Reduced Chi Squared')
  plt.title("N="+str(N))
  plt.savefig('NDResults\AllResults_N{}_{}.png'.format(N, column))
  plt.show()
  
  
# Plot minimum points of a data set against simulation number

array = np.array(dataND.chi)
min_chi_results = []
min_chi =10000
for j in range(len(array)):
  if array[j] < min_chi:
    min_chi = array[j]
    min_chi_results.append([int(j), min_chi])
    
min_chi_results = np.array(min_chi_results)
plt.scatter(min_chi_results[:,0],min_chi_results[:,1], linewidth=1)
plt.xlabel("Simulation number")
plt.ylabel("New minimum reduced Chi squared")
plt.show()

"""
