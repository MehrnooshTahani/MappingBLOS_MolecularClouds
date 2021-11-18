"""
This is the third stage of the BLOSMapping method where magnetic field values are calculated.

We must also calculate the optimal number of off positions

THIS IS THE OLD METHOD
"""
import os
from Classes.BScalingXeAverageLayer0 import ZoomedScalingXeAverage
from MolecularClouds.Classes.RegionOfInterest import Region

# -------- CHOOSE THE REGION OF INTEREST --------
cloudName = 'California'
# cloudName = input("Enter the name of the region of interest: ")
cloudName = cloudName.capitalize()  # Ensure only the first letter is capitalized
regionOfInterest = Region(cloudName)
# -------- CHOOSE THE REGION OF INTEREST. --------

# -------- SET TOTAL NUMBER OF REFERENCE POINTS --------
# Set the total number of reference points:
TotalNRefPoints = 2
# -------- SET TOTAL NUMBER OF REFERENCE POINTS. --------

# -------- DEFINE FILES AND PATHS --------
currentDir = os.path.abspath(os.getcwd())
MatchedRMExtincPath = os.path.join(currentDir, 'FileOutput/'+cloudName+'/MatchedRMExtinction'+cloudName+'.txt')
AvAbundancePath = os.path.join(currentDir, 'Data/ChemicalAbundance/n1.0e3_T12.0_G1/Av_T0_n0.out')
saveFilePath = os.path.join(currentDir, 'FileOutput/'+cloudName+'/NRef'+str(TotalNRefPoints)+'.txt')
# -------- DEFINE FILES AND PATHS. --------

BData = ZoomedScalingXeAverage(AvAbundancePath, MatchedRMExtincPath, TotalNRefPoints, saveFilePath)

print(BData.finalFiducial)
print(BData.finalFiducialExt)
print(BData.finalRMRefErrAvg)
print(BData.FinalFidStandardDev)

print(BData.fiducialIndex)

'''
37.17999999999999
0.462934554
10.3
46.14166880380466
'''



