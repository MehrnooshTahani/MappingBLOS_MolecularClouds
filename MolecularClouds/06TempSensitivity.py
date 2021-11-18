"""
This is the sixth stage the BLOSMapping method where the dependence on temperature is assessed
    - All parameters except for temperature are held constant, and the magnetic field is calculated with electron
    abundances corresponding to a change of +/-  5, 10, and 20 % the fiducial input temperature
"""
from CalculateB import CalculateB
import os
from MolecularClouds.Classes.RegionOfInterest import Region
import pandas as pd

# -------- CHOOSE THE REGION OF INTEREST --------
cloudName = input("Enter the name of the region of interest: ")
cloudName = cloudName.capitalize()  # Ensure only the first letter is capitalized
regionOfInterest = Region(cloudName)
# -------- CHOOSE THE REGION OF INTEREST. --------

# -------- DEFINE FILES AND PATHS --------
currentDir = os.path.abspath(os.getcwd())
MatchedRMExtincPath = os.path.join(currentDir, 'FileOutput/'+cloudName+'/MatchedRMExtinction'+cloudName+'.txt')
RefPointPath = os.path.join(currentDir, 'FileOutput/'+cloudName+'/RefPoints'+cloudName+'.txt')
saveFileDir = os.path.join(currentDir, 'FileOutput/'+cloudName+'/TemperatureSensitivity/')
# -------- DEFINE FILES AND PATHS. --------

# -------- READ REFERENCE POINT TABLE --------
refPointTable = pd.read_csv(RefPointPath)
# -------- READ REFERENCE POINT TABLE. --------

# -------- CALCULATE BLOS AS A FUNCTION OF PERCENT OF THE INPUT TEMPERATURE --------
p = [5, 10, 20]  # Percents of the input temperature

# Calculate BLOS at +/- these percents:
# eg ['-20', '-10', '-5', '0', '+5', '+10', '+20']
percent = ['-{}'.format(i) for i in p[::-1]] + ['0'] + ['+{}'.format(i) for i in p]

for value in percent:
    AvAbundanceName = 'Av_T' + value + '_n0'
    AvAbundancePath = regionOfInterest.AvFileDir + AvAbundanceName + '.out'
    saveFilePath = saveFileDir + 'B_' + AvAbundanceName + '.txt'
    B = CalculateB(AvAbundancePath, MatchedRMExtincPath, refPointTable, saveFilePath)

print('Saving calculated magnetic field values in the folder: '+saveFileDir)
# -------- CALCULATE BLOS AS A FUNCTION OF PERCENT OF THE INPUT TEMPERATURE. --------
