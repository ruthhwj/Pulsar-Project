import multiprocessing as mp
import time
import subprocess
import glob, os
from astropy.io import ascii
import numpy as np
import csv
#should be working
pulsar=["./pulsar-getter.sh", "250", "10.5", "1", "15" , "17", "0.8", "45", "80", "7.7", "0.85", "15", "29", "refpulsar.gg"]
pulsar_arg_names = ["scriptname", "Cone1Intensity", "Cone1BeamAngle", "Cone1BeamletAngle","Cone1NumberOfSparks", "Cone1phi0", "Eccentricity", "Orientation", "Cone2Intensity",
                    "Cone2BeamAngle", "Cone2BeamletAngle","Cone2NumberOfSparks", "Cone2phi0", "Filename"]
pulsar_arg_ranges = [[230, 260], [8, 12], [1, 2], [15, 15], [10,20] , [0.5, 0.9], [40, 50], [60, 120], [6,10], [0.5,1.5], [15,15], [22,32]] #ranges over which to search for each variable
def read_pulsar(string): # Reads ASCII, returns dataframe  #"weak.all37.p3fold.ASCII"
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
    chi = 0
    for i in range(len(intensities_ref)):
        x1 = (intensities_img[i] - intensities_ref[i])
        chi += abs(x1 * x1)/ (RMS_noise * (50 * 300 - 2))
    return chi   #


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

# a
def compare_pulsars_all(pulsar_number, N, intensities_exp_flat):
    df_sim = read_pulsar("SimPulse{}N{}.gg.ASCII".format(str(pulsar_number), str(N)))
    intensities_sim = get_intensities(df_sim, 1)
    chi = fit_measure(intensities_exp_flat, intensities_sim)
    print("returning chi")
    return chi


def pulsar_worker_1d(arg, exp): # int argument,
    n = 1
    N=10
    res = []
    while n<=N:
        pulsar_number=str(n)
        b1=np.random.uniform(pulsar_arg_ranges[arg-1][0], pulsar_arg_ranges[arg-1][1])
        pulsar[13]="SimPulse{}{}.gg".format(pulsar_arg_names[arg], str(pulsar_number))
        pulsar[arg]='{0:.2f}'.format(float(str(b1)))
        for i in pulsar:
            print(i)
        subprocess.run(pulsar)
        try:
            x = compare_pulsars_1d(pulsar_number, pulsar_arg_names[arg], exp)
            result = [b1, x]
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

def pulsar_worker_all(exp, N):
    n = 1
    res = []

    while n<=N:
        pulsar_number=str(n)

        for arg in [x for x in range(1, 13) if (x != 4 and x != 11)]:
            b1 = np.random.uniform(pulsar_arg_ranges[arg - 1][0], pulsar_arg_ranges[arg - 1][1])
            pulsar[arg] = '{0:.2f}'.format(float(str(b1)))

        pulsar[13] = "SimPulse{}N{}.gg".format(str(pulsar_number),str(N))
        subprocess.run(pulsar)

        try:
            chi = compare_pulsars_all(pulsar_number, N, exp)
            print("Reduced chi squared =" + str(chi))
            result = []
            for i in [x for x in range(1,13) if (x!=4 and x!=11)]:
                result.append(pulsar[i])
            result.append(chi)
            res.append(result)
        except Exception:
            print("Skipping")
            continue
        finally:
            print("cleaning up")
            try:
                os.remove("SimPulse{}N{}.gg".format(str(pulsar_number),str(N)))
                os.remove("SimPulse{}N{}.gg.ASCII".format(str(pulsar_number),str(N)))
            except FileNotFoundError as e:
                print("Pulsar number {} in N={} all variable run skipped.".format(str(pulsar_number), str(N)))
            n += 1
    print("writing Results to file")
    with open("AllResults_N{}.csv".format(N), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(res)

df_exp = read_pulsar("weak.all37.p3fold.rebinned.ASCII")  # experimental p3fold here
intensities_exp = brighten(get_intensities(df_exp, 0)).flatten()
intensities_RMS = np.array(df_exp.col4)
exp_croppedlist = ((intensities_RMS.reshape(50, 1123))[:, 0:600]).flatten()  # off pulse RMS noise
RMS_noise = np.var(exp_croppedlist)



def main():
    pool = mp.Pool(mp.cpu_count() + 2)
    subprocess.check_output(pulsar)
    df_sim = read_pulsar("refpulsar.gg.ASCII")
    intensities_sim = get_intensities(df_sim, 1)
    chi = fit_measure(intensities_exp, intensities_sim)
    print("ref chi = {} ".format(chi))
    #fire off workers
    start_time=time.time()
    # N = [100,500,1000,5000,10000]
    # UNCOMMENT FOR 1D
    # for i in [x for x in range(1,13) if (x!=4 and x!=11)]:
    #     pool.apply_async(pulsar_worker_1d, (i, intensities_exp))
    #UNCOMMENT FOR 1D

    #COMMENT FOR 1D
    # for i in N:
    #  job = pool.apply_async(pulsar_worker_all, (intensities_exp,i))
    #COMMENT FOR 1D
    # collect results from the workers through the pool result queue

    #now we are done, kill the listener
    pool.close()
    pool.join()
    print("----%s seconds----" % (time.time()-start_time))
if __name__ == "__main__":
    main()