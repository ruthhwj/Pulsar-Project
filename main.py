import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from astropy.io import ascii

#def generate_pulsars(): 

#def read_pulsars():

#def bi_test():

#def plot_pulsars(): 

#def constrain_parameter(): #will likely extend this to individual functions for each paramter or a class

#test

data = ascii.read("weak.all37.p3fold.ASCII", data_start=1)
print(data)
df = data.to_pandas()
print(df)

plt.hist2d(df.col4, df.col3, bins=(2245, 2245), cmap=plt.cm.jet)
plt.show()

 

