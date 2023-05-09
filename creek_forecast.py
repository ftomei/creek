# creek model
# forecast (hourly) version
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
    alpha = 0.18 / nrIntervals

    interception = getCropInterception(currentDate) / nrIntervals
    infiltration = getSoilInfiltration(WHC, currentDate) / nrIntervals

    if swc < 0:
        # phase 1: soil saturation
        swc += max(prec - interception, 0)
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
    fileName = inputPath + "Forecast_2023_05_10_x2.csv"

    df = pd.read_csv(fileName)
    df.index = pd.to_datetime(df['Datetime'])

    # time step
    timeStep = 3600  # [s]
    nrIntervals = 1

    # [mm] precipitation
    precipitation = df['prec']

    # [m] estimated Level vector
    estLevel = np.zeros(len(precipitation))

    # initialize with first value
    currentDate = df.index[0]
    dateStr = currentDate.strftime("%Y_%m_%d")
    currentWHC = df.WHC[0]
    # [mm] current surface  water content
    swc = -currentWHC

    # main cycle
    for i in range(len(precipitation)):
        currentDate = df.index[i]

        # compute current surface water content and water level
        swc, waterLevel = computeWaterLevel(swc, precipitation[i], currentWHC, currentDate, timeStep)
        estLevel[i] = waterLevel

    # estimated dataset
    df['estLevel'] = estLevel
    xo = df.index
    ye = df.estLevel

    # plot
    plt.figure(figsize=(10, 5))
    plt.subplots_adjust(bottom=0.2)
    plt.xticks(rotation=75)
    ax = plt.gca()
    xfmt = md.DateFormatter('%Y-%m-%d %H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    ax.set_ylim([0.0, 2.0])
    ax.grid(linestyle=':')
    ax.plot(xo, ye, label='Estimated')
    ax.set_ylabel('water level [m]')
    plt.legend()
    plt.savefig(outputPath + "Forecast_" + dateStr + ".png", bbox_inches='tight', dpi=300)


main()
