import multiprocessing as mp
import time
import subprocess
import glob, os
from astropy.io import ascii
import numpy as np
#REBINNED WEAK_ALL37 REBINNED RUTH FILE
pulsar=["./pulsar-getter.sh", "233.940149", "10.5", "1", "15" , "17", "0.85", "45", "0.5", "7.7", "1", "15", "4", "refpulsar.gg"]
pulsar_arg_names = ["scriptname", "Cone1Intensity", "Cone1BeamAngle", "Cone1BeamletAngle","Cone1NumberOfSparks", "Cone1phi0", "Eccentricity", "Orientation", "Cone2Intensity",
                    "Cone2BeamAngle", "Cone2BeamletAngle","Cone2NumberOfSparks", "Cone2phi0", "Filename"]
pulsar_arg_ranges = [[230, 250], [9, 12], [1, 2], [10, 15], [14,20] , [0, 0.99], [40, 55], [0.2, 1], [7.5,8.1], [0.5,2], [5,15], [1 ,7]] #ranges over which to search for each variable


def read_pulsar(string): # Reads ASCII, returns dataframe  #"weak.all37.p3fold.ASCII" "W5testmodel.p3fold.ASCII"
  data = ascii.read(string, data_start=1)
  df = data.to_pandas()
  return df


def get_intensities(df, flag):  # Reads dataframe, returns 50x2246 array for plotting OR as a list
  intensities = np.array(df.col4)  # extract intensities column
  pixelarray = np.array(intensities).reshape(50, 1123)  # shape into array with dimensions of image
  croppedarray = pixelarray[:, 700:1000] #rough onpulse region of exp data

  if flag == 0:
    return croppedarray # want this for plotting
  if flag != 0:
    return croppedarray.flatten()  # want this for analysis


def fit_measure(intensities_ref, intensities_img):
    min_chi = 10000
    chi = 0
    for i in range(len(intensities_ref)):
        x1 = (intensities_img[i] - intensities_ref[i])
        chi += abs(x1 * x1 / intensities_ref[i])
    if chi<min_chi:
        chi=min_chi
    return chi / (RMS_noise * (50 * 300 - 2))  #


def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


def brighten(exp_data):
    i = 0
    j = 0
    param = 1
    while i < 50:
        while j < 300:
            gauss_j = (gaussian(float(j), 85, 25)) #centred on bin 85, width 25
            exp_data[i][j] = ((param*gauss_j+1)*exp_data[i][j])
            j+=1
        i+=1
        j=0
    return exp_data


def compare_pulsars_1d(pulsar_number, pulsar_variable, intensities_exp):
    df_sim = read_pulsar("SimPulse{}{}.gg.ASCII".format(pulsar_variable, pulsar_number))
    intensities_sim = get_intensities(df_sim, 1)
    chi = fit_measure(intensities_exp, intensities_sim)
    print("returning chi")
    return chi


def compare_pulsars_all(pulsar_number, N, intensities_exp_flat):
    df_sim = read_pulsar("SimPulse{}{}.gg.ASCII".format(N, pulsar_number))
    intensities_sim = get_intensities(df_sim, 1)
    chi = fit_measure(intensities_exp_flat, intensities_sim)
    print("returning chi")
    return chi


def pulsar_worker_1d(arg, exp):
    n = 1
    res = []
    while n<=10:
        pulsar_number=str(n)
        b1=np.random.uniform(pulsar_arg_ranges[arg-1][0], pulsar_arg_ranges[arg-1][1])
        pulsar[13]="SimPulse{}{}.gg".format(pulsar_arg_names[arg], str(pulsar_number))
        pulsar[arg]='{0:.2f}'.format(float(str(b1)))
        subprocess.run(pulsar)
        try:
            x = compare_pulsars_1d(pulsar_number, pulsar_arg_names[arg], exp)
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
            n+=1
    print("writing Results{} to file".format(pulsar_arg_names[arg]))
    np.savetxt('results{}.txt'.format(pulsar_arg_names[arg]), res, delimiter=',')

# def pulsar_worker_all(exp, N):
#     n = 1
#     res = []
#     E = [0.7, 0.72, 0.74, 0.76, 0.78, 0.80, 0.82, 0.84, 0.86, 0.88, 0.90]  # avoid weird floating point error
#     osm = [43, 44, 45, 46, 47]
#     while n<=N:
#         pulsar_number=str(n)
#         pulsar[1] = str(np.random.uniform(230, 250))  # 1
#         pulsar[2] = str(np.random.uniform(9, 12))  # 2
#         pulsar[3] = str(np.random.uniform(1, 2))  # 3
#         pulsar[5] = str(np.random.uniform(14,20))
#         pulsar[6] = str(round(np.random.choice(E), 2))
#         pulsar[7] = str(round(np.random.choice(osm), 0))
#         pulsar[8] = str(np.random.uniform(40, 80))  # 7
#         pulsar[9] = str(np.random.uniform(4, 8))  # 8
#         pulsar[10] = str(np.random.uniform(0.5, 1.5))  # 9
#         pulsar[12] = str(np.random.uniform(1,7))
#         pulsar[13] = "SimPulse{}N{}.gg".format(str(pulsar_number),str(N))
#         subprocess.run(pulsar)
#         try:
#             x = compare_pulsars_all(pulsar_number, N, exp)
#             result = []
#             for i in [x for x in range(1,13) if (x!=4 and x!=11)]:
#                 result.append(pulsar[i])
#             result.append(x)
#             res.append(result)
#         except Exception:
#             print("Skipping")
#             continue
#         finally:
#             print("cleaning up")
#             try:
#                 os.remove("SimPulse{}N{}.gg".format(str(pulsar_number),str(N)))
#                 os.remove("SimPulse{}N{}.gg.ASCII".format(str(pulsar_number),str(N)))
#             except FileNotFoundError as e:
#                 print("Pulsar number {} in N={} all variable run skipped.".format(str(pulsar_number), str(N)))
#             n += 1
#     print("writing Results to file")
#     np.savetxt('AllVarResults/results{}.txt'.format(N), res, delimiter=',')


df_exp = read_pulsar("weak.all37.p3fold.rebinned.ASCII")  # experimental p3fold here
intensities_exp = get_intensities(df_exp, 1)
intensities_RMS = np.array(df_exp.col4)
exp_croppedlist = ((intensities_RMS.reshape(50, 1123))[:, 600:700]).flatten()  # off pulse RMS noise
RMS_noise = np.var(exp_croppedlist)



def main():

    pool = mp.Pool(mp.cpu_count() + 2)

    #fire off workers
    start_time=time.time()
    for i in [x for x in range(1,13) if (x!=4 and x!=11)]:
        job = pool.apply_async(pulsar_worker_1d, (i, intensities_exp))

    # collect results from the workers through the pool result queue

    #now we are done, kill the listener
    pool.close()
    pool.join()
    print("----%s seconds----" % (time.time()-start_time))
if __name__ == "__main__":
    main()