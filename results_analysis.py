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


def read_pulsar(string):  # Reads ASCII, returns dataframe  #"weak.all37.p3fold.ASCII" "W5testmodel.p3fold.ASCII"
  data = ascii.read(string, data_start=1)
  df = data.to_pandas()
  return df


        #main code

        #scatter plots of the 2D Monte Carlo results

data = pd.read_csv('results_p1p2.txt', sep=",", header=None)

#08/04/20 results_p1p2 is the N=50 monte carlo simulation
# p1 and p2 were randomly generated between 0.1 and 0.8

data.columns = ["p1", "p2", "fmeasure"] #cone 1 intensity, cone 2 intensity, fit measure

#rdata = data[data['chi'] > 10] # can crop out crazy outliers that make the plots a bit hard to interpret

dplot = plt.figure().gca(projection='3d')
dplot.scatter(data.p1, data.p2, data.fmeasure )
dplot.set_xlabel('Cone 1 intensity')
dplot.set_ylabel('Cone 2 intensity')
dplot.set_zlabel('Fit measure')
plt.show()

plt.scatter(data.p1, data.p2, c=data.fmeasure, cmap='cool')
plt.xlabel('Cone 1 intensity')
plt.ylabel('Cone 2 intensity')
cbar = plt.colorbar()
cbar.set_label('fit measure')

plt.show()


      #Histograms of the pixel intensity ranges
""""
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