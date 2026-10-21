# Criteria-Rainbo model
# Sept. 2026 F.Grazzini, F.Tomei


import numpy as np

RAVONE = 1
QUADERNA = 2
ZENA = 3

NODATA = -9999


# basin parameters
def getBasinParameters(basin):

    if basin == RAVONE:
      alpha_runoff = 0.22     # decay factor: % of runoff that leaves the system in one hour
      zeroIdro = -0.2         # [m] minimum water level
      hMax = 4.0              # [m] maximum water level
      k = 0.10                # factor controlling signal response (higher, increase level)
      referenceLevel = 1.25   # [m]
      swc0 = 18               # [mm] swc value to be associated with the reference Level
      baseFlow = 0.2          # [mm] hourly base flow when basin is over FC

    if basin == QUADERNA:
      alpha_runoff = 0.18     # decay factor: % of runoff that leaves the system in one hour
      zeroIdro = 0.05         # [m] minimum water level
      hMax = 2.6              # [m] maximum water level
      k = 0.10                # factor controlling signal response (higher, increase level)
      referenceLevel = 1.25   # [m]
      swc0 = 17               # [mm] swc value to be associated with the reference level
      baseFlow = 0.4          # [mm] hourly base flow when basin is over FC

    if basin == ZENA:
      alpha_runoff = 0.18     # decay factor: % of runoff that leaves the system in one hour
      zeroIdro = 0.0          # [m] minimum water level
      hMax = 4.0              # [m] maximum water level
      k = 0.11                # factor controlling signal response (higher, increase level)
      referenceLevel = 1.25   # [m]
      swc0 = 16               # [mm] swc value to be associated with the reference level
      baseFlow = 0.2          # [mm] hourly base flow when basin is over FC

    m = (hMax - referenceLevel) / referenceLevel
    
    return alpha_runoff, zeroIdro, hMax, m, k, swc0, baseFlow
    

# Water infiltration in deep soil layer [mm/hour]
# deficit90: current water deficit in 90 cm of soil
def getSoilInfiltration(basin, deficit90):
    if basin == QUADERNA:
        infMax = 2.0        # mm/hour representative of very dry soil
        infMin = 0.2        # mm/hour representative of saturated soil
    if basin == ZENA:
        infMax = 2.0        # mm/hour representative of very dry soil
        infMin = 0.2        # mm/hour representative of saturated soil
    if basin == RAVONE:
        infMax = 6.0        # mm/hour representative of very dry soil
        infMin = 0.2        # mm/hour representative of saturated soil
        
    deficit90max = 150
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
    d0 = hMax / (1.0 + m * np.exp(k * swc0)) - zeroIdro
    waterLevel = hMax / (1.0 + m * np.exp(-k * (swc - swc0))) - d0 * (1.0 - min(swc, swc0)/swc0)
    return waterLevel
    
    
# Estimate surface water content (swc) from water level [m] with basin specific parameters
def estimateSWC(basin, waterLevel):
    alpha_runoff, zeroIdro, hMax, m, k, swc0, baseFlow = getBasinParameters(basin)
    
    if waterLevel <= zeroIdro:
        return 0.0
    
    bestSWC = 0
    w0 = estimateLevel(bestSWC, hMax, m, k, zeroIdro, swc0)
    minDiff = abs(w0 - waterLevel)
    
    for swc in np.arange(0.5, 20.0, 0.5):
        w = estimateLevel(swc, hMax, m, k, zeroIdro, swc0)
        diff = abs(w - waterLevel)

        if diff < minDiff:
            bestSWC = swc
            minDiff = diff
        
    return bestSWC


# Main function transforming inflows in outflows
def computeWaterLevel(basin, currentDate, timeStep, rainfall, currentSwc, currentWHC, currentLeafIntercepted, isBaseFlow):
    nrIntervals = 3600 / timeStep

    alphaRunoff, zeroIdro, hMax, m, k, swc0, baseFlow = getBasinParameters(basin)

    if not isBaseFlow:
        baseFlow = 0.0
        
    baseFlow /= nrIntervals

    # [mm] seasonal max crop interception 
    maxInterception = max(0, maxCropInterception(currentDate) - currentLeafIntercepted)
    leafInterception = min(rainfall * 0.2, maxInterception)
    newLeafIntercepted = currentLeafIntercepted + leafInterception

    # rain reaching the soil [mm]
    rainReachingSoil = rainfall - leafInterception
    # maximum amount of water that can infiltrate into deep soil [mm]
    maxDeepInfiltration = getSoilInfiltration(basin, currentWHC) / nrIntervals
    # current deep infiltration [mm]
    currentDeepInfiltration = min(rainReachingSoil, maxDeepInfiltration)
    
    # Soglia per precipitazioni estreme, rendere specifica per bacino
    THRESHOLD_HOUR = 30.0                              # mm h-1
    threshold = THRESHOLD_HOUR / nrIntervals 

    runoff_fraction = 0.0
    if rainReachingSoil > threshold and currentSwc < -threshold:
        runoff_fraction = 0.1 * (rainReachingSoil - threshold)

    if currentSwc < 0:
        # phase 1: rain reaching the ground infiltrates completely until surface storage is under field capacity
        newSwc = currentSwc + rainReachingSoil - currentDeepInfiltration - runoff_fraction
        newWHC = currentWHC - rainReachingSoil - runoff_fraction
    else:
        # phase 2: rain only partially infiltrates and begins to produce runoff
        runoff = currentSwc * (alphaRunoff / nrIntervals)
        newSwc = currentSwc + rainReachingSoil + baseFlow - runoff - currentDeepInfiltration
        newWHC = currentWHC - currentDeepInfiltration - baseFlow

    if newSwc > 0.0:
        waterLevel = estimateLevel(newSwc + runoff_fraction, hMax, m, k, zeroIdro, swc0)      # [m]
    elif runoff_fraction > 0.:
        waterLevel = estimateLevel(runoff_fraction, hMax, m, k, zeroIdro, swc0)               # [m]
    else:
        waterLevel = zeroIdro                                                                 # [m]

    return waterLevel, newSwc, newWHC, newLeafIntercepted


# scala di deflusso Michele Di Lorenzo Giugno 2023, aggiustato fondo -0.20 dopo alluvione ottobre 2024
def computeDischargeRavone(myLevel):
    return 3.5 * (myLevel + 0.20) ** 1.8


# main function
def creek(basin, df_in, precFieldName, deficit35, deficit90, startLevel):
    # initialize
    df_out = df_in

    # [mm] precipitation
    precipitation = df_in[precFieldName]

    # [m] estimated Level vector
    nrData = len(precipitation)
    estLevel = np.zeros(nrData)
    swc_out = np.zeros(nrData)
    whc_out = np.zeros(nrData)

    # initialize with first value
    currentDate = df_in.index[0]
    timeStep = (df_in.index[1] - df_in.index[0]).total_seconds()
    dateStr = currentDate.strftime("%Y_%m_%d")
    
    # [mm] current water storages (swc: surface and first soil layer)
    currentWHC = deficit90
    LeafIntercepted = 0. 

    if deficit35 > 0:
        # under FC
        currentSWC = -deficit35           
    else:
        # over FC
        if startLevel != NODATA:
            currentSWC = estimateSWC(basin, startLevel)
        else:
            currentSWC = 0

    isBaseFlow = False 
    if deficit35 <= 0 and deficit90 <= 0:
        isBaseFlow = True
            
    # main cycle
    for j in range(nrData):
        currentDate = df_in.index[j]

        # compute current surface water content and water level
        waterLevel, currentSWC, currentWHC, LeafIntercepted = computeWaterLevel(basin, currentDate, timeStep, precipitation[j],
                                                              currentSWC, currentWHC, LeafIntercepted, isBaseFlow)
        estLevel[j] = waterLevel
        swc_out[j] = currentSWC
        whc_out[j] = currentWHC

    # estimated datasets
    df_out['estLevel'] = estLevel
    df_out['swc'] = swc_out
    df_out['WHC90'] = whc_out

    # [m3/s] discharge
    if basin == RAVONE:
        df_out['estQ'] = df_out.estLevel.apply(computeDischargeRavone)

    return df_out
