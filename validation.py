# Creek computing one step at time (like in real time).
# May 2022, Revisited March 2024, -> Criteria-Rainbo
# update July 2025

import matplotlib.pyplot as plt
import matplotlib.dates as md
import pandas as pd
import numpy as np
from scipy.ndimage import shift
from scipy.stats import pearsonr
from scipy.signal import find_peaks
import glob
from datetime import datetime, timedelta
import sys
sys.path.append("../")                      # cerca i moduli anche nella dir sopra
from Criteria_Rainbo_model import *


def nearest_date(items, pivot):
    if len(items) > 0:
        for t in range(len(items)):
            delta = (items[t] - pivot) / np.timedelta64(1, 'h')     # hours
            if (delta >= -0.5) and (delta <= 3.0):
                return items[t]
    return -1


# parameters for peaks recognitions
peak_hmin = 0.2  # [m] hmin for peak search
peak_prominence = 0.1  # [m] minimum peak prominence
peak_width = 2  # [timestep] minimal horizontal distance in samples between neighbouring peaks

basin = QUADERNA

if basin == QUADERNA:
    inputPath = "./INPUT/QUADERNA/"
    outputPath = "./OUTPUT/QUADERNA/"
    criteriaOutputFileName = inputPath + "CriteriaOutput/Quaderna.csv"
    shift_default = 1.5         # hours
    all_files = glob.glob(inputPath + "Quaderna_*.csv")
    precName = 'P30'
elif basin == RAVONE:
    inputPath = "./INPUT/RAVONE/"
    outputPath = "./OUTPUT/RAVONE/"
    criteriaOutputFileName = inputPath + "CriteriaOutput/Ravone.csv"
    shift_default = 0.5         # hours
    all_files = glob.glob(inputPath + "Test_*.csv")
    precName = 'P15'
    
# insert complete filename to read a single test case or wildcard for all cases
df_daily = pd.read_csv(criteriaOutputFileName)
df_daily.index = pd.to_datetime(df_daily['DATE'])

# loop on several cases
list_scores = []
for fileName in all_files:
    df_in = pd.read_csv(fileName)
    df_in.index = pd.to_datetime(df_in['Dataf'])
    del df_in['Dataf']

    # compute time step
    date0 = df_in.index[0]
    date1 = df_in.index[1]
    timeStep = (date1 - date0).seconds # [s]
    nrIntervals = int(3600 / timeStep)

    #print(fileName, date0)
    # [mm] water holding capacity from daily preprocessed data
    deficit35 = max(df_daily[df_daily.index == (date0 - timedelta(days=1)).strftime("%Y-%m-%d")].DEFICIT_35.values[0],0) #WHC = df.WHC[0]
    deficit90 = df_daily[df_daily.index == (date0 - timedelta(days=1)).strftime("%Y-%m-%d")].DEFICIT_90.values[0] 
 
    # Run Criteria-Rainbo model with df_in in input and wch35/whc90
    df = creek(basin, df_in, precName, deficit35, deficit90)
    
    positive_swc = df.index[df.swc > 0].strftime("%d-%m %H:%M").tolist()   
    r_start = positive_swc[0] if len(positive_swc) > 0 else 'No RunOff'
    raincum = df[precName].sum()

    # remove roows with empty data in obs
    df_clean = df[df['Livello'].notna()]

    # estimation array
    vest = df_clean.estLevel.values
    # observed array
    vobs = df_clean.Livello.values
    # index
    xo = df_clean.index

    # searching for multiple peaks and relative statistics
    peaks_obs, prop_peak_obs = find_peaks(vobs, height=peak_hmin, prominence=peak_prominence,width=peak_width)
    peaks_est, prop_peak_est = find_peaks(vest, height=peak_hmin, prominence=peak_prominence,width=peak_width)

    level_peaks_obs = df_clean.Livello[df_clean.index[peaks_obs]]
    level_peaks_est = df_clean.estLevel[df_clean.index[peaks_est]]

    # maximum values
    frame = {'maxOBS': level_peaks_obs, 'maxEST': level_peaks_est}
    df_max = pd.DataFrame(frame)
    df_max = df_max.dropna(axis=0, how='all')
    array_peak_obs = df_max.maxOBS.dropna()
    array_peak_est = df_max.maxEST.dropna()

    dt = []
    dm = []
    # loop sulle date dei massimi osservati per associarli con i previsti
    for timeIndex in array_peak_est.index:
        nearest = nearest_date(array_peak_obs.index, timeIndex)  # trova l'ora del picco più vicino a quello osservato
        if nearest != -1:
            dt.append((nearest - timeIndex) / np.timedelta64(1, 'h'))       # calcola il time shift
            dm.append(df_max.maxEST[timeIndex] - df_max.maxOBS[nearest])    # calcola l'errore

    # mean peaks characteristics
    mPeak_anti = round(np.mean(dt), 2)
    mPeak_err = round(np.mean(dm), 2)

    if (vest == vest[0]).all():
        # tutti dati uguali: no runoff event
        r = np.nan
        r_shift = np.nan
    else:
        # shift estimated data
        if mPeak_anti >= 0:
            shiftNr = round(mPeak_anti * nrIntervals)
        else:
            shiftNr = round(shift_default * nrIntervals)
        vest_shift = shift(vest, shiftNr)

        r, p_value = pearsonr(vobs, vest)
        r_shift, p_value = pearsonr(vobs, vest_shift)
        r = round(r, 3)
        r_shift = round(r_shift, 3)

        RMSE = np.sqrt(((vest - vobs) ** 2).mean())
        RMSE = round(RMSE, 3)

        string_ini = date0.strftime("%d-%m-%Y")
        val_evento = [string_ini, deficit35, r, r_shift, RMSE, mPeak_err, mPeak_anti]
        list_scores.append(val_evento)

    # print
    print("Evento: ", string_ini, "WHC35: ", deficit35, "\tWHC90: ", deficit90,
          "\tRaincum: ", round(raincum, 1), "\tRunoff start: ", r_start)
 
    # plot
    plt.figure(figsize=(10, 5))
    plt.subplots_adjust(bottom=0.2)
    plt.xticks(rotation=75)
    ax = plt.gca()
    xfmt = md.DateFormatter('%Y-%m-%d %H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    if basin == QUADERNA:
        ax.set_ylim([0, 2.5])
    else:
        # Ravone
        ax.set_ylim([-0.2, 4.0])
    ax.grid(linestyle=':')
    ax.plot(xo, vobs, 'r.', label='Observed')
    ax.plot(xo, vest, label='Estimated')
    ax.set_ylabel('water level [m]')
    plt.title('WHCini 35/90=' + str(round(deficit35, 0)) + '/' + str(round(deficit90, 0)) + '   R=' + str(r) + '   RMSE[m]=' + str(RMSE)
              + '   mPeak error[m]=' + str(mPeak_err) + '  mPeak shift[h]=' + str(mPeak_anti), size=12)
    plt.plot(df_max.maxOBS.index, df_max.maxOBS.values, "x")
    plt.plot(df_max.maxEST.index, df_max.maxEST.values, "x")
    plt.legend()
    plt.savefig(outputPath + "Prev_" + string_ini + ".png", bbox_inches='tight', dpi=100)
        
    # write csv out with level, whc, infiltration
    df_max.to_csv(outputPath + "Max_" + string_ini + ".csv")
    df.to_csv(outputPath + "Data_" + string_ini + ".csv", columns=[precName,'WHC90','swc','estLevel','Livello'])

df_out = pd.DataFrame(list_scores, columns=["date", "DEFICIT35", "R", "R_SHIFT", "RMSE", "mP_error", "mP_ant"])
df_out.to_csv(outputPath + "stat_tests.csv")  # salva su csv
print(df_out.describe())
