import glob
import Criteria_Rainbo_model as rainbo

basin = rainbo.ZENA

if basin == rainbo.RAVONE:
    inputPath = "./INPUT/RAVONE/"
    outputPath = "./OUTPUT/RAVONE/"
    criteriaOutputFileName = inputPath + "CriteriaOutput/Ravone.csv"
    all_files = glob.glob(inputPath + "Test_*.csv")
    precName = 'P15'
    alarmLevels = [0.4, 1.4, 2.0]
    shift_hours = 0.25


if basin == rainbo.QUADERNA:
    inputPath = "./INPUT/QUADERNA/"
    outputPath = "./OUTPUT/QUADERNA/"
    criteriaOutputFileName = inputPath + "CriteriaOutput/Quaderna.csv"
    all_files = glob.glob(inputPath + "Quaderna_*.csv")
    precName = 'P30'
    alarmLevels = [0.9, 1.3, 1.7]
    shift_hours = 1.0


if basin == rainbo.ZENA:
    inputPath = "./INPUT/ZENA/"
    outputPath = "./OUTPUT/ZENA/"
    criteriaOutputFileName = inputPath + "CriteriaOutput/Zena.csv"
    all_files = glob.glob(inputPath + "Test_*.csv")
    precName = 'P30'
    alarmLevels = [0.7, 1.3, 2.2]
    shift_hours = 1.0
