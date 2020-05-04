from multiprocessing import Pool
import subprocess
import glob, os
from astropy.io import ascii
import numpy as np

pulsar_arg=["./pulsar-getter.sh", "0.5", "10.5", "1", "15" , "17", "0.85", "45", "0.5", "7.7", "1", "15", "4", "9", "-2.4", "refpulsar.gg"]
pulsar_arg_names = ["scriptname", "Cone1Intensity", "Cone1BeamAngle", "Cone1BeamletAngle","Cone1NumberOfSparks", "Cone1phi0", "Eccentricity", "Orientation", "Cone2Intensity",
                    "Cone2BeamAngle", "Cone2BeamletAngle","Cone2NumberOfSparks", "Cone2phi0", "Alpha", "Beta", "Filename"]
pulsar_arg_ranges = [[0.2, 1], [10, 11.5], [0.5, 2], [10, 15], [14,20] , [0, 0.99], [40, 55], [0.2, 1], [7.5,8.1], [0.5,2], [5,15], [1 ,7] ,[0, 45],[0, -45]]

def read_pulsar(string): # Reads ASCII, returns dataframe  #"weak.all37.p3fold.ASCII" "W5testmodel.p3fold.ASCII"
  data = ascii.read(string, data_start=1)
  df = data.to_pandas()
  return df

def get_intensities(df, flag):  # Reads dataframe, returns 50x2246 array for plotting OR as a list
  intensities = np.array(df.col4)  # extract intensities column
  pixelarray = np.array(intensities).reshape(50, 2246)  # shape into array with dimensions of image
  croppedarray = pixelarray[:, 700:1000] #rough onpulse region of exp data

  if flag == 0:
    return croppedarray # want this for plotting
  if flag != 0:
    return croppedarray.flatten()  # want this for analysis

def fit_measure(intensities_ref, intensities_img):
  chi = 0
  for i in range(len(intensities_ref)):
    x1 = (intensities_img[i] - intensities_ref[i])
    chi += abs(x1 * x1 / intensities_ref[i])
  return chi

def compare_pulsars(pulsar_number, pulsar_arg, intensities_exp, b1):
    results = []
    df_sim = read_pulsar("SimPulse" + pulsar_arg + pulsar_number + ".gg.ASCII")
    intensities_sim = get_intensities(df_sim, 1)
    chi = fit_measure(intensities_exp, intensities_sim)
    print("Pulsar " + pulsar_number + "has a fit measurement of " + str(chi))
    results.append([b1,chi])
    os.remove("SimPulse{}{}" + pulsar_arg + pulsar_number + ".gg")
    os.remove("SimPulse{}{}gg.ASCII")
    return results


def main(i):
    results = []
    c = 1
    while c<=25:
        pulsar_number=str(c)
        b1 = np.random.uniform(pulsar_arg_ranges[i-1][0], pulsar_arg_ranges[i-1][1])
        pulsar_arg[15] ="SimPulse{}{}.gg".format(pulsar_arg_names[i],str(pulsar_number))
        pulsar_arg[i]='{0:.2f}'.format(float(str(b1)))
        print('Argument {} = {}'.format(pulsar_arg_names[i],pulsar_arg[i]))
        proc = subprocess.run(pulsar_arg)
        try:
            results=results+compare_pulsars(pulsar_number, pulsar_arg_names[i],intensities_exp, b1)
        except Exception:
            continue
        finally:
            c+=1
    return results



if __name__=='__main__':
    main()