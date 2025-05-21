# Package containing all functions and main of Criteria-Rainbo model
# March 2025 F.Grazzini, F. Tomei

import numpy as np


# Parameters of the sigmoid function found by fitting with observations
# Ravone basin fit on march 2025 with test_logistic_Rainbo.ipynb
def getBasinParameters():
    zeroIdro = -0.2     # [m] minimum water level
    hMax = 4.5          # [m] maximum water level
    m = 2.1             # factor controlling base level (decreasing, increase level)
    k = 0.11            # factor controlling signal response (higher, increase level)
    x0 = 20             # swc value to be associated with a particolar observed level (in this case L=1.25 m, obtained true confrontation with data)
    return zeroIdro, hMax, m, k, x0


# Water infiltration in deep soil layer [mm/hour]
# deficit90: current water deficit in 90 cm of soil
def getSoilInfiltrationNew(deficit90):
    infMax = 10.0       # mm/hour representative of very dry soil with cracks
    infMin = 0.5        # mm/hour representative of saturated soil
    deficit90max = 100
    deficit90min = -50
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
    vegetationStorage = {1: 2, 2: 2, 3: 3, 4: 5, 5: 6, 6: 8, 7: 10, 8: 9, 9: 7, 10: 6, 11: 5, 12: 3}
    return vegetationStorage[currentDate.month]


# Compute water level [m] from surface water content (swc) with basin specific parameters
def estimateLevel(swc, maxLevel, m, k, zeroIdro, x0):
    waterLevel = maxLevel / (1 + m * np.exp(-k * (swc - x0))) + zeroIdro
    return waterLevel


# Main function transforming inflows in outflows
def computeWaterLevel(currentDate, timeStep, rainfall, currentSwc, currentDeficit90, currentLeafIntercepted):
    zeroIdro, maxLevel, m, k, x0 = getBasinParameters()
    alpha = 0.18  # runoff decay factor, % of runoff that leaves the system in one hour
    nrIntervals = 3600 / timeStep

    # [mm] seasonal max crop interception 
    maxInterception = max(0, maxCropInterception(currentDate) - currentLeafIntercepted)
    currentLeafInterception = min(rainfall * 0.2, maxInterception)
    newLeafIntercepted = currentLeafIntercepted + currentLeafInterception

    # rain reaching the soil [mm]
    rainReachingSoil = rainfall - currentLeafInterception
    # maximum amount of water that can infiltrate into deep soil [mm]
    maxDeepInfiltration = getSoilInfiltrationNew(currentDeficit90) / nrIntervals
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

    if newSwc > 0:
        waterLevel = estimateLevel(newSwc, maxLevel, m, k, zeroIdro, x0)
    else:
        waterLevel = zeroIdro

    return waterLevel, newSwc, newDeficit90, newLeafIntercepted


# scala di deflusso Michele Di Lorenzo Giugno 2023, aggiustato fondo -0.20 dopo alluvione ottobre 2024
def computeDischarge(myLevel):
    return 3.5 * (myLevel + 0.20) ** 1.8

