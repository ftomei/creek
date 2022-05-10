# creek model
import matplotlib.pyplot as plt
import matplotlib.dates as md
import pandas as pd
import numpy as np
import glob


# hourly soil infiltration [mm/hour]
def getSoilInfiltration(WHC, currentDate):
    if currentDate.month > 6 and WHC > 30:
        # soil cracking
        return 5.0
    else:
        # saturated soil conductivity
        return 1.5


# hourly crop interception [mm/hour]
def getCropInterception(currentDate):
    if 3 < currentDate.month < 11:
        return 0.3
    else:
        return 0.0


# estimate water level [m]
# swc [mm] = surface water content
# WHC [mm] = water holding capacity at the start of the day
# prec [mm] = precipitation during the timeStep [s]
def getWaterLevel(swc, prec, WHC, currentDate, timeStep):
    nrIntervals = int(3600 / timeStep)

    # [-] runoff decay factor
    alpha = 0.16 / nrIntervals

    interception = getCropInterception(currentDate) / nrIntervals
    infiltration = getSoilInfiltration(WHC, currentDate) / nrIntervals

    if swc < 0:
        # phase 1: soil saturation
        swc += prec
        if swc > 0:
            waterLevel = 3.8 / (1 + 20 * np.exp(-0.15 * swc)) - 0.15
        else:
            waterLevel = 0
    else:
        # phase 2: runoff
        w0 = swc * (1 - alpha)
        waterInput = max(prec - interception, 0)
        swc = max(w0 + waterInput - infiltration, 0)
        waterLevel = 3.8 / (1 + 20 * np.exp(-0.15 * swc)) - 0.15

    return swc, waterLevel


def main():
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
        timeStep = (date1 - date0).seconds  # [s]

        # [mm] water holding capacity
        WHC = df.WHC[0]
        dateStr = date0.strftime("%Y_%m_%d")
        print(dateStr, "   WHC: ", WHC)

        # [mm] precipitation
        precipitation = df['P15']

        # estimation vector
        estLevel = np.zeros(len(precipitation))

        # surface water content [mm]
        swc = -WHC
        for i in range(len(precipitation)):
            swc, waterLevel = getWaterLevel(swc, precipitation[i], WHC, date0, timeStep)
            estLevel[i] = waterLevel

        df['estLevel'] = estLevel

        # clean dataset
        df_clean = df[df['Livello'].notna()]
        # estimation vector
        df_est = df_clean[['estLevel']]
        # observed vector
        df_obs = df_clean[['Livello']]

        df_visu = df_est.join(df_obs, how='outer').dropna()

        xo = df_visu.index
        yo = df_visu.Livello
        ye = df_visu.estLevel

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
        plt.legend()
        plt.savefig(outputPath + "Prev_" + dateStr + ".png", bbox_inches='tight', dpi=300)


main()
