# Criteria-Rainbo model
# calibration of logistic regression
# December 2025 F.Grazzini, F.Tomei

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import seaborn as sns
from scipy.optimize import curve_fit

from basin import *


# Compute water level [m] from surface water content (swc)
def estimateLevel(swc, swc0, k):
    if basin == rainbo.QUADERNA:
        alphaRunoff, zeroIdro, hMax, m, k_est, swc0_est = rainbo.getBasinParameters_Quaderna()
    elif basin == rainbo.ZENA:
        alphaRunoff, zeroIdro, hMax, m, k_est, swc0_est = rainbo.getBasinParameters_Zena()
    else:
        # default: ravone
        alphaRunoff, zeroIdro, hMax, m, k_est, swc0_est = rainbo.getBasinParameters_Ravone()

    d0 = hMax / (1.0 + m * np.exp(k * swc0)) - zeroIdro
    waterLevel = hMax / (1.0 + m * np.exp(-k * (swc - swc0))) - d0 * np.maximum(0.0, 1.0 - swc / swc0)
    return waterLevel


def main():
    # insert complete filename to read a single test case
    files = glob.glob(f"{outputPath}/Data_*.csv")

    # Create an empty dataframe to store the combined data
    combined_df = pd.DataFrame()

    # Loop through each CSV file and append its contents to the combined dataframe
    for csv_file in files:
        df = pd.read_csv(csv_file)
        combined_df = pd.concat([combined_df, df])

    # index
    combined_df.index = pd.to_datetime(combined_df['Dataf'])

    # compute time step
    date0 = combined_df.index[0]
    date1 = combined_df.index[1]
    timeStep = (date1 - date0).seconds      # [s]
    nrIntervals = int(3600 / timeStep)

    shiftNr = -round(shift_hours * nrIntervals)
    combined_df['level'] = combined_df.Livello.shift(shiftNr)
    df_clean = combined_df[['Dataf', 'swc', 'level']].copy()

    df = df_clean.dropna().copy()
    df['Ltrend'] = np.where(df['level'].diff() > 0, 'salita', 'discesa')
    # mantiene solo dati nella fase di crescita
    df = df[df['Ltrend'] == 'salita']

    xdata = df.swc.values
    ydata = df['level'].values

    df_plot = df[df.swc >= -10]

    popt, pcov = curve_fit(estimateLevel, xdata, ydata, None, method='dogbox', maxfev=5000)
    print("swc0 = ", popt[0])
    print("k = ", popt[1])

    plt.figure(figsize=(12, 7))

    sns.scatterplot(data=df_plot, x='swc', y='level', hue='Ltrend', style=df_plot.index.year)
    xteo = np.arange(0, 60, 1)
    y = estimateLevel(xteo, *popt)
    plt.plot(xteo, y, label='fit')
    plt.legend()

    estLevel = estimateLevel(xdata, *popt)
    r, p_value = pearsonr(estLevel, ydata)
    print("r = ", r)

    ax = plt.gca()
    # Ravone
    ax.set_ylim([-0.2, 4.0])

    if basin == rainbo.QUADERNA:
        ax.set_ylim([0, 2.5])
    elif basin == rainbo.ZENA:
        ax.set_ylim([0, 3.0])

    plt.show()

main()
