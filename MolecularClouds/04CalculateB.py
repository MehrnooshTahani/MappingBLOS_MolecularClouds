"""
This file is for calculating BLOS

This file was made for debugging / testing code - I do not think it is needed anymore
See instead 0XOptimalNumRefPoints
"""
import os
from Classes.CalculateB import CalculateB
from MolecularClouds.Classes.RegionOfInterest import Region

# -------- CHOOSE THE REGION OF INTEREST --------
cloudName = 'OrionB'
# cloudName = input("Enter the name of the region of interest: ")
cloudName = cloudName.capitalize()  # Ensure only the first letter is capitalized
regionOfInterest = Region(cloudName)
# -------- CHOOSE THE REGION OF INTEREST. --------

# -------- SET NUMBER OF REFERENCE POINTS --------
NRef = 5
# -------- SET NUMBER OF REFERENCE POINTS. --------

# -------- DEFINE FILES AND PATHS --------
currentDir = os.path.abspath(os.getcwd())
MatchedRMExtincPath = os.path.join(currentDir, 'FileOutput/'+cloudName+'/MatchedRMExtinction'+cloudName+'.txt')
RefPointPath = os.path.join(currentDir, 'FileOutput/'+cloudName+'/RefPoints/NRef'+str(NRef)+'.txt')

saveFilePath = os.path.join(currentDir, 'FileOutput/'+cloudName+'/BLOS/BLOSNRef'+str(NRef)+'.txt')
# -------- DEFINE FILES AND PATHS. --------

CalculateB(regionOfInterest.AvFilePath, MatchedRMExtincPath, RefPointPath, saveFilePath)
