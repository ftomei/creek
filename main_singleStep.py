import matplotlib.pyplot as plt
import matplotlib.dates as md
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from scipy.signal import find_peaks
import glob


def nearest_date(items, pivot):
    if len(items) > 0:
        return min(items, key=lambda x: abs(x - pivot))
    else:
        return -1


# [mm/hour]
def getSoilInfiltration(currentWHC, currentDate):
    if currentDate.month > 6 and currentWHC > 30:
        # soil cracking
        return 5.0
    else:
        # saturated soil conductivity
        return 1.5


# [mm/hour]
def getCropInterception(currentDate):
    if 3 < currentDate.month < 11:
        return 0.3
    else:
        return 0.0


def main():
    # parameters
    peak_hmin = 0.1  # [m] hmin for peak search
    peak_prominence = 0.1  # [m] minimum peak prominence
    peak_width = 0  # [timestep] minimal horizontal distance in samples between neighbouring peaks

    inputPath = ".\\INPUT\\"
    outputPath = ".\\OUTPUT\\"
    # insert complete filename to read a single test case
    all_files = glob.glob(inputPath + "Test_*.csv")

    # loop on several cases
    list_scores = []
    for fileName in all_files:
        df = pd.read_csv(fileName)
        df.index = pd.to_datetime(df['Dataf'])

        # compute time step
        date0 = df.index[0]
        date1 = df.index[1]
        timeStep = (date1 - date0).seconds          # [s]
        nrIntervals = int(3600 / timeStep)

        # [mm] water holding capacity
        WHC = df.WHC[0]

        # [mm/h] soil infiltration and crop interception
        hourlyInfiltration = getSoilInfiltration(WHC, date0)    # [mm/h]
        hourlyInterception = getCropInterception(date0)         # [mm/h]
        infiltration = hourlyInfiltration / nrIntervals
        interception = hourlyInterception / nrIntervals

        # [-] runoff decay factor
        alpha = 0.16 / nrIntervals

        # rainfall [mm]
        rainfall = df['P15']

        # estimation vector
        estLevel = np.zeros(len(rainfall))

        # surface water content [mm]
        swc = -WHC
        start_index = 0
        for i in range(len(rainfall)):
            if swc <= 0:
                # phase 1: soil saturation
                swc += rainfall[i]
                if swc > 0:
                    start_index = i
            else:
                # phase 2: runoff
                w0 = swc * (1 - alpha)
                prec = max(rainfall[i] - interception, 0)
                swc = max(w0 + prec - infiltration, 0)
                # water level [m]
                estLevel[i] = 3.8 / (1 + 20 * np.exp(-0.15 * swc)) - 0.15

        df['estLevel'] = estLevel
        r_start = df.index[start_index]

        # clean dataset (only event)
        # df_event = df[df.index >= r_start]
        df_clean = df[df['Livello'].notna()]

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

        # maximum values
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
        val_evento = [string_ini, WHC, infiltration, r, RMSE, mPeak_err, mPeak_anti]
        list_scores.append(val_evento)

        # print
        print("WHC: ", WHC, "  K: ", hourlyInfiltration, "  Runoff start: ", r_start)
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
        plt.title('Inf=' + str(hourlyInfiltration) + '   R=' + str(r) + '   RMSE[m]=' + str(RMSE)
                  + '   mPeak error[m]=' + str(mPeak_err) + '  mPeak shift[h]=' + str(mPeak_anti), size=12)
        plt.plot(df_max.maxOBS.index, df_max.maxOBS.values, "x")
        plt.plot(df_max.maxFC.index, df_max.maxFC.values, "x")
        plt.legend()
        # plt.show()
        plt.savefig(outputPath + "Prev_" + string_ini + ".png", bbox_inches='tight', dpi=300)

    df_out = pd.DataFrame(list_scores, columns=["date", "WHC", "Inf", "R", "RMSE", "mP_error", "mP_ant"])
    df_out.to_csv(outputPath + "Ravone_stat_tests.csv")  # salva su csv


main()