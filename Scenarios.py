# creek model - Tomei & Grazzini

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import timedelta
from matplotlib.ticker import MultipleLocator

from basin import *

NODATA = -9999


# select case
fileName = ""
if basin == rainbo.RAVONE:
    fileName = 'Test_2024-10-19.csv'
elif basin == rainbo.QUADERNA:
    fileName = 'Quaderna_2024_10_19.csv'
    #fileName = 'Quaderna_2015_03_25.csv'
elif basin == rainbo.ZENA:
    fileName = 'Test_2024-12-08.csv'
else:
    print("Wrong basin")
    exit()


# select case
df_in = pd.read_csv(inputPath + fileName)
df_in.index = pd.to_datetime(df_in['Dataf'])
del df_in['Dataf']

# compute time step
date0 = df_in.index[0]
date1 = df_in.index[1]
timeStep = (date1 - date0).seconds  # [s]
nrIntervals = int(3600 / timeStep)

# [mm] water deficit from CRITERIA-1D output
df_daily = pd.read_csv(criteriaOutputFileName)
df_daily.index = pd.to_datetime(df_daily['DATE'])
deficit35 = df_daily[df_daily.index == (date0 - timedelta(days=1)).strftime("%Y-%m-%d")].DEFICIT_35.values[0]
deficit90 = df_daily[df_daily.index == (date0 - timedelta(days=1)).strftime("%Y-%m-%d")].DEFICIT_90.values[0]

isBaseFlow = (deficit35 < 0 and deficit90 < 0)
# cut surface storage
deficit35 = max(deficit35, -2.0)  # [mm]

# compute
df_out = rainbo.creek(basin, df_in, precName, deficit35, deficit90, isBaseFlow)

# prepare the figure
plt.subplots(figsize=(15, 7))
plt.subplots_adjust(bottom=0.2)
plt.xticks(rotation=80)

ax = sns.lineplot(data=df_out, x=df_out.index.strftime("%d/%m %H:%M"),
                  y='Livello', linewidth=1.5, label='Obs level', color='black')

sns.lineplot(data=df_out, x=df_out.index.strftime("%d/%m %H:%M"),
             y='estLevel', linewidth=1.5, label='Real-time fcst', color='blue')

# y axes
ax.set_ylabel('Water level [m]')
ax.grid(linestyle='')
levelObsMax = int(df_out['Livello'].dropna().max()) + 1
levelEstMax = int(max(df_out['estLevel'])) + 1
levelMax = max(levelObsMax, levelEstMax)
levelMax = max(levelMax, 3.0)                   # [m]
ax.set_ylim([0.0, levelMax])

# secondary axes: prec
# duplicate x axes
axp = ax.twinx()
precMax = int(max(df_out[precName])) + 1
precMax = max(precMax, 5)                       # [mm]
axp.set_ylim([0.0, precMax])

sns.barplot(data=df_out, x=df_out.index.strftime("%d/%m %H:%M"),
            y=precName, alpha=0.5, color='steelblue', ax=axp)

axp.set_ylabel('Rainfall [mm]')

# generate sensitivy changing values of WHC
whc90 = [-50, 0, 50, 100, 150, 200]
colors = ['black', 'red', 'orange', 'lightgreen', 'green', 'pink']
for i in range(len(whc90)):
    # initial surface content
    whc35 = whc90[i] * 0.3
    whc35 = max(whc35, -2.0)  # [mm]
    isBaseFlow = (whc90[i] < 0)
    df_out = rainbo.creek(basin, df_in, precName, whc35, whc90[i], isBaseFlow)
    sns.lineplot(data=df_out, x=df_out.index.strftime("%d/%m %H:%M"), y='estLevel',
                 label='Deficit = ' + str(whc90[i]), color=colors[i], linewidth=2, linestyle="dotted", ax=ax)

ax.axhline(alarmLevels[0], linestyle='dashed', color='yellow', label='warning')
ax.axhline(alarmLevels[1], linestyle='dashed', color='orange', label='prealarm')
ax.axhline(alarmLevels[2], linestyle='dashed', color='red', label='alarm')

# x axis
ax.xaxis.set_major_locator(MultipleLocator(6))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.set(xlabel=None)

# title
firstDate = date0.strftime("%Y-%m-%d")
plt.title('Sensitivity soil state ' + " " + firstDate + " - Deficit (90cm) = " + str(deficit90))
sns.move_legend(ax, "upper left")

outputFileName = outputPath + "Scenarios_" + firstDate + ".png"
plt.savefig(outputFileName, bbox_inches='tight', dpi=100)

plt.show()

