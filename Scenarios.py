# creek model - Tomei & Grazzini

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import timedelta
from matplotlib.ticker import MultipleLocator
import Criteria_Rainbo_model as rainbo

NODATA = -9999

# select case
basin = rainbo.QUADERNA
#fileName = 'Test_2024-10-19.csv'   # Ravone
fileName = 'Quaderna_2023_05_01.csv'
timeName = 'Dataf'

if basin == rainbo.RAVONE:
    inputPath = "./INPUT/RAVONE/"
    outputPath = "./OUTPUT/RAVONE/"
    criteriaOutputFileName = inputPath + "CriteriaOutput/Ravone.csv"
    precName = 'P15'
    alarmLevels = [0.4, 1.4, 2.0]

if basin == rainbo.QUADERNA:
    inputPath = "./INPUT/QUADERNA/"
    outputPath = "./OUTPUT/QUADERNA/"
    criteriaOutputFileName = inputPath + "CriteriaOutput/Quaderna.csv"
    precName = 'P30'
    alarmLevels = [0.9, 1.3, 1.7]

# select case
df_in = pd.read_csv(inputPath + fileName)
df_in.index = pd.to_datetime(df_in[timeName])
del df_in[timeName]

# compute time step
date0 = df_in.index[0]
date1 = df_in.index[1]
timeStep = (date1 - date0).seconds  # [s]
nrIntervals = int(3600 / timeStep)

# [mm] water holding capacity from Criteria1D data
df_daily = pd.read_csv(criteriaOutputFileName)
df_daily.index = pd.to_datetime(df_daily['DATE'])
deficit35 = max(df_daily[df_daily.index == (date0 - timedelta(days=1)).strftime("%Y-%m-%d")].DEFICIT_35.values[0], 0)
deficit90 = df_daily[df_daily.index == (date0 - timedelta(days=1)).strftime("%Y-%m-%d")].DEFICIT_90.values[0]

# compute
df_out = rainbo.creek(basin, df_in, precName, deficit35, deficit90)

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

# secondary axes: prec
# duplicate x axes
axp = ax.twinx()
precMax = int(max(df_out[precName])) + 1
precMax = max(precMax, 5)
axp.set_ylim([0.0, precMax])

sns.barplot(data=df_out, x=df_out.index.strftime("%d/%m %H:%M"),
            y=precName, alpha=0.5, color='steelblue', ax=axp)

axp.set_ylabel('Rainfall [mm]')

# generate sensitivy changing values of WHC
whc90 = [0, 50, 100, 150, 200]
colors = ['red', 'orange', 'lightgreen', 'green', 'pink']
for i in range(len(whc90)):
    whc35 = whc90[i] * 0.4
    df_out = rainbo.creek(basin, df_in, precName, whc35, whc90[i])
    sns.lineplot(data=df_out, x=df_out.index.strftime("%d/%m %H:%M"), y='estLevel',
                 label='Deficit = ' + str(whc90[i]), color=colors[i], linewidth=2, linestyle="dotted", ax=ax)

ax.axhline(alarmLevels[0], linestyle='dashed', color='yellow', label='warning')
ax.axhline(alarmLevels[1], linestyle='dashed', color='orange', label='prealarm')
ax.axhline(alarmLevels[2], linestyle='dashed', color='red', label='alarm')

# x axis
ax.xaxis.set_major_locator(MultipleLocator(6))
ax.xaxis.set_minor_locator(MultipleLocator(1))
#ax.grid(linestyle=':')
ax.set(xlabel=None)

# title
firstDate = date0.strftime("%Y-%m-%d")
plt.title('Sensitivity soil state ' + " " + firstDate + " - Deficit (90cm) = " + str(deficit90))
sns.move_legend(ax, "upper left")

outputFileName = outputPath + "Scenarios_" + firstDate + ".png"
plt.savefig(outputFileName, bbox_inches='tight', dpi=100)

plt.show()

