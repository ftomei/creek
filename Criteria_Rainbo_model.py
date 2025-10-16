# functions of Criteria-Rainbo model
# June 2025 F.Grazzini, F.Tomei
# prova 

import numpy as np

RAVONE = 1
QUADERNA = 2


# Parameters of the sigmoid function found by fitting with observations
# Ravone basin
def getBasinParameters_Ravone():
    zeroIdro = -0.2         # [m] minimum water level
    hMax = 4.5              # [m] maximum water level
    k = 0.10                # factor controlling signal response (higher, increase level)
    referenceLevel = 1.25   # [m]
    swc0 = 22               # [mm] swc value to be associated with the reference Level
    m = (hMax - (referenceLevel - zeroIdro)) / (referenceLevel - zeroIdro)
    return zeroIdro, hMax, m, k, swc0

# Quaderna basin
def getBasinParameters_Quaderna():
    zeroIdro = 0.1          # [m] minimum water level
    hMax = 2.5              # [m] maximum water level
    k = 0.10                # factor controlling signal response (higher, increase level)
    referenceLevel = 1.25    # [m]
    swc0 = 16              # [mm] swc value to be associated with the reference level
    m = (hMax - (referenceLevel - zeroIdro)) / (referenceLevel - zeroIdro)
    return zeroIdro, hMax, m, k, swc0


# Water infiltration in deep soil layer [mm/hour]
# deficit90: current water deficit in 90 cm of soil
def getSoilInfiltration(basin, deficit90):
    if basin == QUADERNA:
        infMax = 2.0        # mm/hour representative of very dry soil
        infMin = 0.2        # mm/hour representative of saturated soil
    else:
        infMax = 6.0        # mm/hour representative of very dry soil
        infMin = 0.2        # mm/hour representative of saturated soil
    deficit90max = 100
    deficit90min = -40
    if deficit90 > deficit90max:
        currentInf = infMax
    elif deficit90 < deficit90min:
        currentInf = infMin
    else:
        ratio = (deficit90 - deficit90min) / (deficit90max - deficit90min)
        currentInf = infMin + ratio*ratio * infMax
    return currentInf


# Vegetation maximum water storage [mm]
def maxCropInterception(currentDate):
    # {month:val} maximum water storage vegetation [mm]
    vegetationStorage = {1: 2, 2: 2, 3: 3, 4: 5, 5: 6, 6: 7, 7: 8, 8: 7, 9: 7, 10: 5, 11: 4, 12: 3}
    return vegetationStorage[currentDate.month]


# Compute water level [m] from surface water content (swc) with basin specific parameters
def estimateLevel(swc, hMax, m, k, zeroIdro, swc0):
    waterLevel = hMax / (1 + m * np.exp(-k * (swc - swc0))) + zeroIdro
    return waterLevel


# Main function transforming inflows in outflows
def computeWaterLevel(basin, currentDate, timeStep, rainfall, currentSwc, currentDeficit90, currentLeafIntercepted):
    alpha = 0.18     # runoff decay factor, % of runoff that leaves the system in one hour
    nrIntervals = 3600 / timeStep

    # [mm] seasonal max crop interception 
    maxInterception = max(0, maxCropInterception(currentDate) - currentLeafIntercepted)
    currentLeafInterception = min(rainfall * 0.2, maxInterception)
    newLeafIntercepted = currentLeafIntercepted + currentLeafInterception

    # rain reaching the soil [mm]
    rainReachingSoil = rainfall - currentLeafInterception
    # maximum amount of water that can infiltrate into deep soil [mm]
    maxDeepInfiltration = getSoilInfiltration(basin, currentDeficit90) / nrIntervals
    # current deep infiltration [mm]
    currentDeepInfiltration = min(rainReachingSoil, maxDeepInfiltration)

    if currentSwc < 0:
        # phase 1: rain reaching the ground infiltrates completely until it fills the surface storage (first 35 cm)
        newSwc = currentSwc + rainReachingSoil - currentDeepInfiltration
        newDeficit90 = currentDeficit90 - rainReachingSoil
    else:
        # phase 2: rain only partially infiltrates and begins to produce runoff
        runoff = currentSwc * (alpha / nrIntervals)
        newSwc = max(currentSwc + rainReachingSoil - runoff - currentDeepInfiltration, 0)
        newDeficit90 = currentDeficit90 - currentDeepInfiltration

    # basin parameters
    if basin == RAVONE:
        zeroIdro, hMax, m, k, swc0 = getBasinParameters_Ravone()
    elif basin == QUADERNA:
        zeroIdro, hMax, m, k, swc0 = getBasinParameters_Quaderna()
    else:
        zeroIdro, hMax, m, k, swc0 = getBasinParameters_Ravone()

    if newSwc > 0:
        waterLevel = estimateLevel(newSwc, hMax, m, k, zeroIdro, swc0)      # [m]
    else:
        waterLevel = zeroIdro                                               # [m]

    return waterLevel, newSwc, newDeficit90, newLeafIntercepted


# main looping over precipitation a calling other functions
def creek(basin, df_in, precFieldName, deficit35, deficit90):
    # initialize
    df_out = df_in

    # [mm] precipitation
    precipitation = df_in[precFieldName]

    # [m] estimated Level vector
    nrData = len(precipitation)
    estLevel = np.zeros(nrData)
    swcout = np.zeros(nrData)
    whc90out = np.zeros(nrData)

    # initialize with first value
    currentDate = df_in.index[0]
    timeStep = (df_in.index[1] - df_in.index[0]).total_seconds()
    dateStr = currentDate.strftime("%Y_%m_%d")

    # [mm] current water storages (swc: surface and first soil layer)
    swc = min(-deficit35, 0)
    currentWHC90 = deficit90
    LeafIntercepted = 0

    # main cycle
    for j in range(nrData):
        currentDate = df_in.index[j]

        # compute current surface water content and water level
        waterLevel, swc, currentWHC90, LeafIntercepted = computeWaterLevel(basin, currentDate, timeStep, precipitation[j], swc,
                                                                           currentWHC90, LeafIntercepted)
        estLevel[j] = waterLevel
        swcout[j] = swc
        whc90out[j] = currentWHC90

    # estimated datasets
    df_out['estLevel'] = estLevel
    df_out['swc'] = swcout
    df_out['WHC90'] = whc90out

    return df_out
