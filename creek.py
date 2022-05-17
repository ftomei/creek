# creek model
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md

NODATA = -9999


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


# update swc and estimate current water level [m]
# swc [mm] = surface water content
# WHC [mm] = water holding capacity at the start of the day
# prec [mm] = precipitation during the timeStep [s]
def computeWaterLevel(swc, prec, WHC, currentDate, timeStep):
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


def initializeArray(myArray):
    for i in range(len(myArray)):
        myArray[i] = NODATA


def main():
    inputPath = ".\\INPUT\\"
    outputPath = ".\\OUTPUT\\"
    fileName = inputPath + "Test_2022-04-21.csv"

    df = pd.read_csv(fileName)
    df.index = pd.to_datetime(df['Dataf'])

    # time step
    timeStep = (df.index[1] - df.index[0]).seconds  # [s]

    # [mm] precipitation
    precipitation = df['P15']

    # [mm] previous precipitation vector
    nrIntervals = int(3600 / timeStep)
    previousPrec = np.zeros(24 * nrIntervals + 1)
    initializeArray(previousPrec)
    indexPreviousPrec = 0

    # [m] estimated Level vector
    estLevel = np.zeros(len(precipitation))

    # initialize with first value
    currentDate = df.index[0]
    dateStr = currentDate.strftime("%Y_%m_%d")
    currentWHC = df.WHC[0]
    # [mm] current surface  water content
    swc = -currentWHC
    # [mm] surface water content at 00:00
    swc0 = swc

    # main cycle
    for i in range(len(precipitation)):
        currentDate = df.index[i]

        # 00:00 initialize previous precipitation and save current swc
        if currentDate.hour == 0 and currentDate.minute == 0:
            initializeArray(previousPrec)
            indexPreviousPrec = 0
            swc0 = swc

        # new value of WHC: update current swc using previous precipitation
        if not pd.isna(df.WHC[i]):
            currentWHC = df.WHC[i]
            dateStr = currentDate.strftime("%Y_%m_%d")
            print(dateStr, "   WHC: ", currentWHC, "   SWC0:", round(swc0, 2))
            if swc0 <= 0:
                swc = -currentWHC
                for j in range(indexPreviousPrec):
                    if previousPrec[j] != NODATA:
                        swc, waterLevel = computeWaterLevel(swc, previousPrec[j], currentWHC,
                                                            currentDate, timeStep)

        # compute current surface water content and water level
        swc, waterLevel = computeWaterLevel(swc, precipitation[i], currentWHC, currentDate, timeStep)
        estLevel[i] = waterLevel

        # save previous precipitation
        previousPrec[indexPreviousPrec] = precipitation[i]
        indexPreviousPrec += 1

    # clean dataset
    df['estLevel'] = estLevel
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
