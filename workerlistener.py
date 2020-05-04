import multiprocessing as mp
import time
import subprocess
import glob, os
from astropy.io import ascii
import numpy as np
#REBINNED WEAK_ALL37 REBINNED RUTH FILE
pulsar=["./pulsar-getter.sh", "0.5", "10.5", "1", "15" , "17", "0.85", "45", "0.5", "7.7", "1", "15", "4", "refpulsar.gg"]
pulsar_arg_names = ["scriptname", "Cone1Intensity", "Cone1BeamAngle", "Cone1BeamletAngle","Cone1NumberOfSparks", "Cone1phi0", "Eccentricity", "Orientation", "Cone2Intensity",
                    "Cone2BeamAngle", "Cone2BeamletAngle","Cone2NumberOfSparks", "Cone2phi0", "Filename"]
pulsar_arg_ranges = [[0.2, 1], [10, 11.5], [0.5, 2], [10, 15], [14,20] , [0, 0.99], [40, 55], [0.2, 1], [7.5,8.1], [0.5,2], [5,15], [1 ,7]]

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

# CHANGE FIT MEASURE USING RMS NOISE FROM RUTHS BRANCH
# CHANGE FIT MEASURE USING RMS NOISE FROM RUTHS BRANCH
def fit_measure(intensities_ref, intensities_img):
  chi = 0
  for i in range(len(intensities_ref)):
    x1 = (intensities_img[i] - intensities_ref[i])
    chi += abs(x1 * x1 / intensities_ref[i])
  return chi
# CHANGE FIT MEASURE USING RMS NOISE FROM RUTHS BRANCH
# CHANGE FIT MEASURE USING RMS NOISE FROM RUTHS BRANCH

def compare_pulsars(pulsar_number, pulsar_variable, intensities_exp):
    df_sim = read_pulsar("SimPulse{}{}.gg.ASCII".format(pulsar_variable,pulsar_number))
    intensities_sim = get_intensities(df_sim, 1)
    chi = fit_measure(intensities_exp, intensities_sim)
    print("returning chi")
    return chi

def pulsar_worker(arg, exp):
    N = 1
    res = []
    while N<=10000:
        pulsar_number=str(N)
        b1=np.random.uniform(pulsar_arg_ranges[arg-1][0], pulsar_arg_ranges[arg-1][1])
        pulsar[13]="SimPulse{}{}.gg".format(pulsar_arg_names[arg], str(pulsar_number))
        pulsar[arg]='{0:.2f}'.format(float(str(b1)))
        subprocess.run(pulsar)
        try:
            x = compare_pulsars(pulsar_number, pulsar_arg_names[arg] ,exp)
            result = [b1,x]
            res.append(result)
        except Exception:
            print("Skipping")
            continue
        finally:
            print("cleaning up")
            try:
                os.remove("SimPulse{}{}.gg".format(pulsar_arg_names[arg], str(pulsar_number)))
                os.remove("SimPulse{}{}.gg.ASCII".format(pulsar_arg_names[arg], str(pulsar_number)))
            except FileNotFoundError as e:
                print("Value for {} at number {} skipped.".format(pulsar_arg_names[arg], str(pulsar_number)))

            N+=1
    print("writing Results{} to file".format(pulsar_arg_names[arg]))
    np.savetxt('Results/results{}.txt'.format(pulsar_arg_names[arg]), res, delimiter=',')

# def pulsar_listener(q, fn):
#     '''listens for messages on the q, writes to file. '''
#
#     with open(fn, 'w') as f:
#         while 1:
#             m = q.get()
#             if m == 'kill':
#                 f.write('killed')
#                 break
#             f.write(str(m) + '\n')
#             f.flush()

def main():
    df_exp = read_pulsar("norm_exp.ASCII")
    intensities_exp = get_intensities(df_exp, 1)

    #must use Manager queue here, or will not work
    # manager = mp.Manager()
    # q = manager.Queue()
    pool = mp.Pool(mp.cpu_count() + 2)

    #fire off workers
    for i in [x for x in range(1,13) if (x!=4 and x!=11)]:
        job = pool.apply_async(pulsar_worker, (i, intensities_exp))


    # collect results from the workers through the pool result queue

    #now we are done, kill the listener
    pool.close()
    pool.join()

if __name__ == "__main__":
   main()