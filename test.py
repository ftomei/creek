# creek model
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
from Criteria_Rainbo_model import computeWaterLevel, computeDischarge


# main looping over precipitation a calling other functions
def creek(df, precFieldName, swc35, deficit90):
    # [mm] precipitation
    precipitation = df[precFieldName]

    # [m] estimated Level vector
    nrData = len(precipitation)
    estLevel = np.zeros(nrData)

    # initialize with first value
    timeStep = (df.index[1] - df.index[0]).total_seconds()

    # [mm] current water storages (swc: surface and first soil layer)
    index_swc = 0
    currentSWC = max(swc35[index_swc], 0)
    currentDeficit90 = deficit90[index_swc]
    LeafIntercepted = 0                         # TODO transition

    # main cycle
    previousDate = df.index[0]
    for i in range(nrData):
        currentDate = df.index[i]
        if currentDate.date() != previousDate.date():
            index_swc += 1
            if index_swc <= len(swc35)-1:
                currentSWC = (swc35[index_swc] + currentSWC) * 0.5
                currentDeficit90 = deficit90[index_swc]
            previousDate = currentDate

        # compute current surface water content and water level
        waterLevel, currentSWC, currentDeficit90, LeafIntercepted = computeWaterLevel(currentDate, timeStep,
                                              precipitation[i], currentSWC, currentDeficit90, LeafIntercepted)
        estLevel[i] = waterLevel

    # estimated dataset
    df['estLevel'] = estLevel

    # [m3/s] discharge
    df['estQ'] = df.estLevel.apply(computeDischarge)

    return df


def main():
    inputPath = ".\\INPUT\\"
    outputPath = ".\\OUTPUT\\"

    fileName = inputPath + "Test_2023-05-01.csv"
    #fileName = inputPath + "Test_2015-03-25.csv"

    df = pd.read_csv(fileName)
    df.index = pd.to_datetime(df['Dataf'])

    currentDate = df.index[0]
    dateStr = currentDate.strftime("%Y_%m_%d")

    # swc: deficit con segno invertito
    swc35 = [-34, -4, +14]      # dati di maggio 23
    deficit90 = [84, 51, -10]

    #swc35 = [-1]                 # dati di marzo 2015
    #deficit90 = [-11]

    df = creek(df, 'P15', swc35, deficit90)

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
    ax.set_ylim([0.0, 2.0])
    ax.grid(linestyle=':')
    ax.plot(xo, yo, 'r.', label='Observed')
    ax.plot(xo, ye, label='Estimated')
    ax.set_ylabel('water level [m]')
    plt.legend()
    plt.savefig(outputPath + "Test_" + dateStr + ".png", bbox_inches='tight', dpi=300)


main()
