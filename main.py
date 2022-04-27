import matplotlib.pyplot as plt
import matplotlib.dates as md
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from scipy.signal import find_peaks
import glob


def assign_prec(w1, shift):
    prec = np.zeros(len(w1))
    for i in range(len(prec)):
        if w1[i] > 0:
            if w1[i - shift] <= 0:
                prec[i] = w1[i]
            else:
                prec[i] = w1[i] - w1[i - shift]
    return prec


def assign_runoff(surfaceWaterContent, currentPrec, maxInfiltration, start_index):
    runoff = np.zeros(len(currentPrec))
    for i in range(len(currentPrec)):
        if i >= start_index:
            runoff[i] = np.maximum(surfaceWaterContent[i] + currentPrec[i] - maxInfiltration, 0)
            for j in range(i+1, len(surfaceWaterContent)):
                surfaceWaterContent[j] = max(0, surfaceWaterContent[j] - runoff[i] * 0.04)
        else:
            runoff[i] = 0
    return runoff


def search_runoff_start(w1, p15):
    start = 0
    noPrecCount = 0
    for i in range(len(w1)):
        if start == 0 and w1[i] > 0:
            start = i
            noPrecCount = 0
        if p15[i] == 0:
            noPrecCount += 1
        else:
            if start > 0 and noPrecCount >= 2 and w1[i] < 4:
                start = i
            noPrecCount = 0
    return start


def nearest_date(items, pivot):
    if len(items) > 0:
        return min(items, key=lambda x: abs(x - pivot))
    else:
        return -1


# [mm/hour] infiltration
def getKSoil(currentWHC, currentDate):
    if currentDate.month > 6 and currentWHC > 30:
        # soil cracking
        return 5.0
    else:
        # ksat
        return 1.5


# [mm/hour] crop interception
def getKCrop(currentDate):
    if 3 < currentDate.month < 11:
        return 0.3
    else:
        return 0.0


def main():
    # parameters
    alpha = 0.1  # [-]  time factor
    peak_hmin = 0.1  # [m] hmin for maxlevel search
    peak_prominence = 0.1  # [m] prevalenza minima del picco
    peak_width = 0  # [timestep] minimal horizontal distance in samples between neighbouring peaks

    path = ".\\INPUT\\"
    # insert complete filename to read a single test case
    all_files = glob.glob(path + "/Test_*.csv")

    # loop on several cases
    list_scores = []
    for fileName in all_files:
        df = pd.read_csv(fileName)

        # index = date
        df['date'] = pd.to_datetime(df['Dataf'])
        df.index = df['date']
        date0 = df.date[0]

        # [mm] water holding capacity
        WHC = df.WHC[0]
        k_soil = getKSoil(WHC, date0)
        k_crop = getKCrop(date0)

        # clean dataset
        del df['WHC']
        del df['date']
        del df['Datai']
        del df['Dataf']

        # compute w1
        df['w1'] = df['P15'].cumsum() - WHC

        start_index = search_runoff_start(df.w1, df.P15)
        r_start = df.index[start_index]

        # time: nr hours after runoff start
        df['time'] = np.maximum(0, (df.index - r_start) / np.timedelta64(1, 'h'))  # [hours]
        # time factor
        df['factor'] = np.exp(-alpha * df.time)

        # compute runoff
        nrIntervals = 4
        waterOut = (k_soil + k_crop) * 0.25 * nrIntervals
        currentPrec = assign_prec(df.w1, nrIntervals)
        surface_wc = np.maximum(df.w1.shift(nrIntervals) - df.time.shift(nrIntervals) * (k_soil + k_crop), 0)
        runoff = np.maximum(currentPrec + surface_wc * df.factor.shift(nrIntervals) - waterOut, 0)
        # runoff = assign_runoff(surface_wc, currentPrec, waterOut, start_index)

        # forecast water level
        df['estLevel'] = 3.8 / (1 + 20 * np.exp(-0.15 * runoff)) - 0.15

        # clean dataset (only event)
        df_event = df[df.index >= r_start]
        df_clean = df_event[df_event['Livello'].notna()]

        # estimation vector
        df_est = df_clean[['estLevel']]
        # observed vector
        df_obs = df_clean[['Livello']]

        # searching for multiple peaks and relative statistics
        peaks_obs, prop_peak_obs = find_peaks(df_clean.Livello.values, height=peak_hmin, prominence=peak_prominence,
                                              width=peak_width)
        peaks_fc, prop_peak_fc = find_peaks(df_clean.estLevel.values, height=peak_hmin, prominence=peak_prominence,
                                            width=peak_width)

        maxOBS = df_clean.Livello[df_clean.index[peaks_obs]]
        maxFC = df_clean.estLevel[df_clean.index[peaks_fc]]

        # df dei valori massimi
        frame = {'maxOBS': maxOBS, 'maxFC': maxFC}
        df_max = pd.DataFrame(frame)
        df_max = df_max.dropna(axis=0, how='all')
        array_do = df_max.maxOBS.dropna()
        array_df = df_max.maxFC.dropna()

        print(df_max.head())

        dt = []
        dm = []
        listaoss = []
        listafo = []
        # loop sulle date dei massimi osservati per associarli con i previsti
        for i in array_do.index:
            nearest = nearest_date(array_df.index, i)  # trovo l'ora del picco più vicino a quello osservato
            if nearest != -1:
                listaoss.append(i)
                listafo.append(nearest)
                dt.append((nearest - i) / np.timedelta64(1, 'h'))  # calcolo il time shift
                dm.append(df_max.maxFC[nearest] - df_max.maxOBS[i])  # calcolo l'errore

        # mean peaks characteristics
        mPeak_anti = round(np.mean(dt), 2)
        mPeak_err = round(np.mean(dm), 2)

        df_visu = df_est.join(df_obs, how='outer').dropna()

        xo = df_visu.index
        yo = df_visu.Livello
        ye = df_visu.estLevel

        r, p_value = pearsonr(yo, ye)
        r = round(r, 2)
        RMSE = np.sqrt(((ye - yo) ** 2).mean())
        RMSE = round(RMSE, 2)

        string_ini = r_start.strftime("%Y_%m_%d")
        val_evento = [string_ini, WHC, k_soil, r, RMSE, mPeak_err, mPeak_anti]
        list_scores.append(val_evento)

        # print
        print("WHC: ", WHC, "  K: ", k_soil, "  Runoff start: ", r_start)
        print("Observed peaks: ", listaoss)
        print("Forecast peaks: ", listafo)

        # plot
        plt.figure(figsize=(10, 5))
        plt.subplots_adjust(bottom=0.2)
        plt.xticks(rotation=75)
        ax = plt.gca()
        xfmt = md.DateFormatter('%Y-%m-%d %H:%M')
        ax.xaxis.set_major_formatter(xfmt)
        ax.set_ylim([0.0, 2.0])
        ax.grid(linestyle=':')
        ax.plot(xo, yo, 'r.', label='Observed')
        ax.plot(xo, ye, label='Estimated')
        ax.set_ylabel('water level [m]')
        plt.title('Ksoil=' + str(k_soil) + '   R=' + str(r) + '   RMSE[m]=' + str(RMSE)
                  + '   mPeak error[m]=' + str(mPeak_err) + '  mPeak shift[h]=' + str(mPeak_anti), size=12)
        plt.plot(df_max.maxOBS.index, df_max.maxOBS.values, "x")
        plt.plot(df_max.maxFC.index, df_max.maxFC.values, "x")
        plt.legend()
        # plt.show()
        plt.savefig('Prev_' + string_ini + '.png', bbox_inches='tight', dpi=300)

    df_out = pd.DataFrame(list_scores, columns=["date", "WHC", "K_soil", "R", "RMSE", "mP_error", "mP_ant"])
    df_out.to_csv("Ravone_stat_tests.csv")  # salva su csv


main()
