# creek model
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as md
from Criteria_Rainbo_model import *


# main looping over precipitation a calling other functions
def creek(basin, df, precFieldName, swc35, deficit90):
    # [mm] precipitation
    precipitation = df[precFieldName]

    # [m] estimated Level vector
    nrData = len(precipitation)
    estLevel = np.zeros(nrData)

    # initialize with first value
    timeStep = (df.index[1] - df.index[0]).total_seconds()

    # [mm] current water storages (swc: surface and first soil layer)
    index = 0
    currentSWC = min(swc35[index], 0)
    currentDeficit90 = deficit90[index]
    LeafIntercepted = 0     # TODO transition

    # main cycle
    previousDate = df.index[0]
    for i in range(nrData):
        currentDate = df.index[i]
        if currentDate.date() != previousDate.date():
            index += 1
            if index < len(deficit90):
                currentDeficit90 = deficit90[index]
            if index < len(swc35):
                if currentSWC < 0:
                    currentSWC = min(swc35[index], 0)
            previousDate = currentDate

        # compute current surface water content and water level
        waterLevel, currentSWC, currentDeficit90, \
        LeafIntercepted = computeWaterLevel(basin, currentDate, timeStep, precipitation[i],
                                            currentSWC, currentDeficit90, LeafIntercepted)
        estLevel[i] = waterLevel

    # estimated dataset
    df['estLevel'] = estLevel
    return df


def main():
    basin = QUADERNA

    if basin == QUADERNA:
        inputPath = ".\\INPUT\\QUADERNA\\"
        outputPath = ".\\OUTPUT\\QUADERNA\\"

        precipitationField = "P30"
        fileName = "Quaderna_2023_05_01.csv"
        #fileName = "Quaderna_2023_05_16.csv"
        #fileName = "Quaderna_2023_05_16_radar_corrected.csv"

        # swc35: deficit35 con segno invertito
        swc35 = [-21.1, 11.4, 27.8]             # dati 01-03 maggio 23 Quaderna
        deficit90 = [47.4, 15.5, -4.1]

        #swc35 = [20.1, 31.31]                   # dati 16 17 maggio 23 Quaderna suolo MGG
        #deficit90 = [-34.21, -49.82]

    else:
        inputPath = ".\\INPUT\\RAVONE\\"
        outputPath = ".\\OUTPUT\\RAVONE\\"

        precipitationField = "P15"
        fileName = "Ravone_2023_05_01.csv"

        # swc35: deficit35 con segno invertito
        swc35 = [-34, -4, +14]          # dati 01-03 maggio 2023 Ravone
        deficit90 = [84, 51, -10]

    fullFileName = inputPath + fileName

    df = pd.read_csv(fullFileName)
    df.index = pd.to_datetime(df['Dataf'])

    currentDate = df.index[0]
    dateStr = currentDate.strftime("%Y_%m_%d")

    df = creek(basin, df, precipitationField, swc35, deficit90)

    # estimation vector
    df_est = df[['estLevel']]

    # observed vector
    df_obs = df[['Livello']]

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
    ax.set_ylim([0.0, 3.0])
    ax.grid(linestyle=':')
    ax.plot(xo, yo, 'r.', label='Observed')
    ax.plot(xo, ye, label='Estimated')
    ax.set_ylabel('water level [m]')
    plt.legend()

    outputFileName = outputPath + fileName[:6] + "_" + dateStr + ".png"
    plt.savefig(outputFileName, bbox_inches='tight', dpi=300)

    print("Output file: ", outputFileName)


main()
