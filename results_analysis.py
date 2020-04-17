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



def read_pulsar(string):  # Reads ASCII, returns dataframe  #"weak.all37.p3fold.ASCII" "W5testmodel.p3fold.ASCII"
  data = ascii.read(string, data_start=1)
  df = data.to_pandas()
  return df

def get_intensities(df, flag):  # Reads dataframe, returns 50x2246 array for plotting OR as a list
  intensities = np.array(df.col4)  # extract intensities column
  pixelarray = np.array(intensities).reshape(50, 2246)  # shape into array with dimensions of image
  croppedarray = pixelarray[:, 1400:2000]  # rough onpulse region of exp data

  if flag == 0:
    return croppedarray  # want this for plotting
  if flag != 0:
    return croppedarray.flatten()  # want this for analysis

def plot_pulsar(df_pixelarray, title):
 print("Plotting pulsar...")
 plt.imshow(df_pixelarray, 'afmhot', origin='lower', interpolation='none', aspect='auto',extent=[1400,2000,0,50])
 plt.set_title(title)
 plt.colorbar()
 plt.show()

      #main code

        #######scatter plots of the 1D Monte Carlo results






data2d = pd.read_csv('results_b1b2.txt', sep=",", header=None)

# 16/04/20 results of N=500 monte carlo simulation for cone intensities
#
""""
data2d.columns = ["b1", "b2", "fmeasure"]

print(data2d.nsmallest(10, 'fmeasure')) #print parameters which give the minimum fmeasures
print(data2d.nlargest(10, 'fmeasure'))
#PLOT 2: scatter plot
plt.scatter(data2d.b2, data2d.fmeasure, linewidth=1)#cool,BrBg, twilight_shifted
plt.xlabel('Cone 2 Half Opening Beam Angle')
plt.ylabel('Reduced Chi Squared')
plt.show()




        #######scatter plots of the 2D Monte Carlo results

#data = pd.read_csv('results_p1p2.txt', sep=",", header=None)

#08/04/20 results_p1p2 is the N=50 monte carlo simulation
# p1 and p2 were randomly generated between 0.1 and 0.8

data = pd.read_csv('results_p1p2N500.txt', sep=",", header=None)

#08/04/20 results_p1p2 is the N=500 monte carlo simulation
# p1 and p2 were randomly generated between 0 and 1

data = pd.read_csv('results_b1b2.txt', sep=",", header=None)

#09/04/20 results_b1b2 is the N=500 monte carlo simulation
# b1 and b2 were randomly generated between 5 and 15 (half opening angle of the beam)

data.columns = ["p1", "p2", "fmeasure"] #cone 1 intensity, cone 2 intensity, fit measure
#rdata = data[data['chi'] > 10] # can crop out crazy outliers that make the plots a bit hard to interpret


#PLOT 1: 3d scatter plot
dplot = plt.figure().gca(projection='3d')
dplot.scatter(data2d.a1, data2d.a2, data2d.fmeasure,  c=data2d.fmeasure, cmap='BrBG', linewidth=1 )
dplot.set_xlabel('Cone 1 Intensity')
dplot.set_ylabel('Cone 2 Intensity')
dplot.set_zlabel('Reduced Chi Squared')
plt.show()
"""

#Looking at the N=500 data shows that the fit measures for (p1,p2) = (1,0) and (0.4,0) are very similar
#this is due to global_norm normalising the peak intensity value to 1 regardless of intensity.
#so we need to analyse the ratio of p1/p2.
"""

#PLOT 2: present as 2d with colourmap as fitmeasure
plt.scatter(data2d.b1, data2d.b2, c=data2d.fmeasure, cmap='BrBG', linewidth=1)#cool,BrBg, twilight_shifted
plt.xlabel('Cone 1 Half Opening Beam Angle')
plt.ylabel('Cone 2 Half Opening Beam Angle')
cbar = plt.colorbar()
cbar.set_label('Reduced Chi Squared')
plt.show()


#PLOT 3: make it a 3d surface
ax = plt.axes(projection='3d')
ax.plot_trisurf(data.p1, data.p2, data.fmeasure, cmap='cool') #the sickest plot you've ever seen
plt.show()
"""
      #Histograms of the pixel intensity ranges

df_exp = read_pulsar("weak.all37.p3fold.ASCII")
df_model = read_pulsar("refpulsar.gg.ASCII")
"""
exp = get_intensities(df_exp, 1)
model = get_intensities(df_model, 1)

print("max intensity of model is " + str(max(model)))
print("min intensity of model is " + str(min(model)))
print("max intensity of exp is " + str(max(exp)))
print("min intensity of exp is " + str(min(exp)))


A = [(exp),(model)]

plt.boxplot(A, labels=["Experimental Intensities", "Model Intensities"], sym="")
plt.show()

#without outliers
#plt.boxplot(df_exp.col4)


df = pd.DataFrame([exp, model], columns=["Experimental Intensities", "Model Intensities"])
boxplot = df.boxplot(column=["Experimental Intensities", "Model Intensities"], sym="")
print("here")
boxplot.show()"""



p1_0 = read_pulsar("a1a2_3_pulsar.ASCII")
p2_0 = read_pulsar("b1b2_pulsar.ASCII")

p1 = get_intensities(p1_0,0)
p2 = get_intensities(p2_0,0)
exp = get_intensities(df_exp,0)

dif_1 = abs(exp - p1)
dif_2 = abs(exp - p2)


plt.imshow(dif_1, cmap='afmhot', aspect='auto', extent=[1400,2000,0,50])
plt.colorbar()
plt.show()

plt.imshow(dif_2, cmap='afmhot', aspect='auto', extent=[1400,2000,0,50])
plt.colorbar()
plt.show()
