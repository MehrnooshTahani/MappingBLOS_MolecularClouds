"""
This is the fifth stage of the BLOSMapping method where the dependence on density is assessed.
    - All parameters except for density are held constant, and the magnetic field is calculated with electron abundances
    corresponding to changes of 0 and +/- 1, 2.5, 5, 10, 20, 20, 40, and 50 % the fiducial input density
"""
from LocalLibraries.CalculateB import CalculateB
import os
from MolecularClouds.LocalLibraries.RegionOfInterest import Region
import pandas as pd
import LocalLibraries.config as config

# -------- CHOOSE THE REGION OF INTEREST --------
cloudName = config.cloud
regionOfInterest = Region(cloudName)
# -------- CHOOSE THE REGION OF INTEREST. --------

# -------- DEFINE FILES AND PATHS --------
MatchedRMExtincPath = os.path.join(config.dir_root, config.dir_fileOutput, config.cloud, config.prefix_RMExtinctionMatch + config.cloud + '.txt')
RefPointPath = os.path.join(config.dir_root, config.dir_fileOutput, config.cloud, config.prefix_selRefPoints + config.cloud + '.txt')
saveFileDir = os.path.join(config.dir_root, config.dir_fileOutput, config.cloud, config.dir_densitySensitivity)
# -------- DEFINE FILES AND PATHS. --------

# -------- READ REFERENCE POINT TABLE --------
refPointTable = pd.read_csv(RefPointPath)
# -------- READ REFERENCE POINT TABLE. --------

# -------- CALCULATE BLOS AS A FUNCTION OF PERCENT OF THE INPUT DENSITY --------
p = [1, 2.5, 5, 10, 20, 30, 40, 50]  # Percents of the input density

# Calculate BLOS at +/- these percents:
# eg ['-50', '-40', '-30', '-20', '-10', '-5', '-2.5', '-1', '0', '+1', '+2.5', '+5', '+10', '+20', '+30', '+40', '+50']
percent = ['-{}'.format(i) for i in p[::-1]] + ['0'] + ['+{}'.format(i) for i in p]

for value in percent:
    AvAbundanceName = 'Av_T0_n' + value
    AvAbundancePath = regionOfInterest.AvFileDir + os.sep + AvAbundanceName + '.out'
    saveFilePath = saveFileDir + os.sep + 'B_' + AvAbundanceName + '.txt'
    B = CalculateB(AvAbundancePath, MatchedRMExtincPath, refPointTable, saveFilePath)

print('Saving calculated magnetic field values in the folder: '+saveFileDir)
# -------- CALCULATE BLOS AS A FUNCTION OF PERCENT OF THE INPUT DENSITY. --------
