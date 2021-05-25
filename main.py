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
            if w1[i-shift] <= 0:
                prec[i] = w1[i]
            else:
                prec[i] = w1[i] - w1[i-shift]
    return prec


def search_runoff_start(w1, p15):
    start = 0
    count = 0
    for i in range(len(w1)):
        if w1[i] > 0:
            if start == 0:
                start = i
                count = 0
        if p15[i] == 0:
            count += 1
        else:
            if start > 0 and count >= 2 and w1[i] < 5:
                start = i
            count = 0
    return start

def nearest_date(items,pivot):
    nearest=min(items, key=lambda x: abs(x - pivot))
    return nearest


if __name__ == '__main__':
    # parameters
    delta = 0              # [hours] shift estimation - observed
    alpha = 0.1            # [-]  time factor
    peak_hmin = 0.1        # [m] hmin for maxlevel search
    peak_prominence = 0.1 # [m] prevalenza minima del picco
    peak_width = 0         # [timestep] minimal horizontal distance in samples between neighbouring peaks


    # read single test case
    #fileName = r'.\INPUT\Test_feb2015.csv'
    #df = pd.read_csv(fileName)

    path = ".\\INPUT\\"
    all_files = glob.glob(path + "/Test_*.csv")

    #loop on several cases
    list_scores = []
    for fileName in all_files:
            df = pd.read_csv(fileName)

            # WHC and k_soil
            WHC = df.WHC[0]     # [mm] water holding capacity
            kmin = 1.6 #[mm/hour]
            kmax = 2.6 #[mm/hour]
            k_soil = WHC/10
            k_soil = np.clip(k_soil, kmin, kmax)  # [mm/hour] infiltration coefficent function of initial WHC value (bounded between kmin and kmax mm/h)

            # index = date
            df['date'] = pd.to_datetime(df['Dataf'])
            df.index = df['date']

            # clean dataset
            del df['WHC']
            del df['date']
            del df['Datai']
            del df['Dataf']

            # compute w1
            df['w1'] = df['P15'].cumsum() - WHC

            start_index = search_runoff_start(df.w1, df.P15)
            r_start = df.index[start_index]
            print(r_start)


            # time: nr hours after runoff start
            df['time'] = np.maximum(0, (df.index - r_start) / np.timedelta64(1, 'h'))  # [hours]
            # time factor
            df['factor'] = np.exp(-alpha * df.time)

            # compute runoff
            nrIntervals = 4
            currentPrec = assign_prec(df.w1, nrIntervals)
            currentInfiltration = k_soil * 0.25 * nrIntervals
            surface_wc = np.maximum(df.w1.shift(nrIntervals) - df.time.shift(nrIntervals) * k_soil, 0)
            runoff = np.maximum(currentPrec + surface_wc * df.factor.shift(nrIntervals) - currentInfiltration, 0)

            # forecast water level
            df['estLevel'] = 3.8 / (1 + 20 * np.exp(-0.15 * runoff)) - 0.15

            # clean dataset (only event)
            df_event = df[df.index >= r_start]
            df_clean = df_event[df_event['Livello'].notna()]

            # estimation vector
            df_est = df_clean[['estLevel']]

            # observed vector
            df_shift = df_clean.shift(-delta, freq='h')
            df_obs = df_shift[['Livello']]

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

            print (df_max.head())

            dt = []
            dm = []
            listaoss = []
            listafo = []
            # loop sulle date dei massimi osservati per associarli con i previsti
            for i in array_do.index:
                nearest = nearest_date(array_df.index, i)  # trovo l'ora del picco più vicino a quello osservato
                listaoss.append(i)
                listafo.append(nearest)
                dt.append((nearest - i) / np.timedelta64(1, 'h'))  # calcolo il time shift
                dm.append(df_max.maxFC[nearest] - df_max.maxOBS[i])  # calcolo l'errore

            # mean peaks characteristics
            mPeak_anti = np.mean(dt)
            mPeak_err = np.mean(dm)

            df_visu = df_est.join(df_obs, how='outer').dropna()

            xo = df_visu.index
            yo = df_visu.Livello
            ye = df_visu.estLevel

            r, p_value = pearsonr(yo, ye)
            RMSE = np.sqrt(((ye - yo) ** 2).mean())

            string_ini = r_start.strftime("%d%m%Y")
            val_evento = [string_ini,WHC,k_soil,r, RMSE, mPeak_err, mPeak_anti]
            list_scores.append(val_evento)

            #print
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
            plt.title('R=' + str(round(r, 3)) + '   RMSE[m]=' + str(round(RMSE, 2))
                      + '   mPeak error[m]=' + str(round(mPeak_err, 2)) + '  mPeak anticipation[h]= '
                                               + str(round(mPeak_anti, 2)), size=12)
            plt.plot(df_max.maxOBS.index, df_max.maxOBS.values, "x")
            plt.plot(df_max.maxFC.index, df_max.maxFC.values, "x")
            plt.legend()
            #plt.show()
            plt.savefig('Prev_' + string_ini + '.png', bbox_inches='tight', dpi=300)

df_out = pd.DataFrame(list_scores, columns=["date","WHC","K_soil","R", "RMSE", "mP_error", "mP_ant"])
df_out.to_csv("Ravone_stat_tests.csv") # salva su csv
