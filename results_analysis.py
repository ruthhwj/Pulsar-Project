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


        #main code

        #######scatter plots of the 1D Monte Carlo results

#data1d = pd.read_csv('results_b1.txt', sep=",", header=None)

# 08/04/20 results_b1 is the N=500 monte carlo simulation
# the half opening angle of the first cone was randomly generated between 5 and 15 degrees
# min fmeasure given by b1 = 8.773202

#data1d = pd.read_csv('results_b2.txt', sep=",", header=None)

# 09/04/20 results_b2 is the N=500 monte carlo simulation
# the half opening angle of the second cone was randomly generated between 5 and 15 degrees
# min feasure of 0.000123 given by b2 = 9.797606 (for usual ref pulsar para but b1 = 8.773202)


data4d = pd.read_csv('results_normdisable_a1a2b1b2.txt', sep=",", header=None)

# 14/04/20 results of N=1000 monte carlo simulation with -normdisable
#

data4d.columns = ["a1", "a2", "b1", "b2", "fmeasure"]

print(data4d.nsmallest(10, 'fmeasure')) #print parameters which give the minimum fmeasures

#PLOT 2: scatter plot
#plt.scatter(data1d.b2, data1d.fmeasure, linewidth=1)#cool,BrBg, twilight_shifted
#plt.xlabel('Cone 2 half opening beam angle')
#plt.ylabel('Fit measure')
#plt.show()
""""
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
dplot.scatter(data.p1, data.p2, data.fmeasure,  c=data.fmeasure, cmap='BrBG', linewidth=1 )
dplot.set_xlabel('Cone 1 half opening angle')
dplot.set_ylabel('Cone 2 half opening angle')
dplot.set_zlabel('Fit measure')
plt.show()

"""
#Looking at the N=500 data shows that the fit measures for (p1,p2) = (1,0) and (0.4,0) are very similar
#this is due to global_norm normalising the peak intensity value to 1 regardless of intensity.
#so we need to analyse the ratio of p1/p2.

"""
#PLOT 2: present as 2d with colourmap as fitmeasure
plt.scatter(data.p1, data.p2, c=data.fmeasure, cmap='BrBG', linewidth=1)#cool,BrBg, twilight_shifted
plt.xlabel('Cone 1 half opening angle')
plt.ylabel('Cone 2 half opening angle')
cbar = plt.colorbar()
cbar.set_label('Fit measure')

plt.show()
#PLOT 3: make it a 3d surface
ax = plt.axes(projection='3d')
ax.plot_trisurf(data.p1, data.p2, data.fmeasure, cmap='cool') #the sickest plot you've ever seen
plt.show()

      #Histograms of the pixel intensity ranges

df_exp = read_pulsar("weak.all37.p3fold.ASCII.normalised.ASCII")
df_model_D = read_pulsar("testmodel.gg.D.p3fold.ASCII")
df_model = read_pulsar("testmodel.noise.normalised.ASCII")

print("max intensity of model is " + str(max(df_model.col4)))
print("min intensity of model is " + str(min(df_model.col4)))
print("max intensity of exp is " + str(max(df_exp.col4)))
print("min intensity of exp is " + str(min(df_exp.col4)))


#A = [(df_model.col4), (df_exp.col4), (df_model_D.col4)]

#plt.boxplot(A)

#without outliers
plt.boxplot(df_model.col4, sym='')
plt.boxplot(df_exp.col4, sym='')

plt.show()
"""