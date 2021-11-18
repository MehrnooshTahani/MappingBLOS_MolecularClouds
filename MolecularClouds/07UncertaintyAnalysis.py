"""
This is the seventh stage of the BLOSMapping method where the uncertainties in the BLOS values are calculated
"""
import os
import pandas as pd
from Classes.RegionOfInterest import Region


# -------- FUNCTION DEFINITION --------
def extinctionChemUncertainties(B, BHigher, BLower):
    if max(B, BHigher, BLower) == B:
        BLowerValue = min(BHigher, BLower)
        BHigherValue = B
    elif min(B, BHigher, BLower) == B:
        BHigherValue = max(BHigher, BLower)
        BLowerValue = B
    else:
        BHigherValue = max(BHigher, BLower)
        BLowerValue = min(BHigher, BLower)

    upperDelta = BHigherValue - B
    lowerDelta = B - BLowerValue

    return upperDelta, lowerDelta
# -------- FUNCTION DEFINITION. --------


# -------- CHOOSE THE REGION OF INTEREST --------
cloudName = input("Enter the name of the region of interest: ")
cloudName = cloudName.capitalize()  # Ensure only the first letter is capitalized
regionOfInterest = Region(cloudName)
# -------- CHOOSE THE REGION OF INTEREST. --------

# -------- DEFINE FILES AND PATHS --------
currentDir = os.path.abspath(os.getcwd())

BFilePath = os.path.join(currentDir, 'FileOutput/' + cloudName + '/BLOSPoints' + cloudName + '.txt')
BData_Density50IncreasePath = os.path.join(currentDir,
                                           'FileOutput/' + cloudName + '/DensitySensitivity/B_Av_T0_n+50.txt')
BData_Density50DecreasePath = os.path.join(currentDir,
                                           'FileOutput/' + cloudName + '/DensitySensitivity/B_Av_T0_n-50.txt')
BData_Temp20IncreasePath = os.path.join(currentDir,
                                        'FileOutput/' + cloudName + '/TemperatureSensitivity/B_Av_T+20_n0.txt')
BData_Temp20DecreasePath = os.path.join(currentDir,
                                        'FileOutput/' + cloudName + '/TemperatureSensitivity/B_Av_T-20_n0.txt')
saveFilePath = os.path.join(currentDir, 'FileOutput/' + cloudName + '/FinalBLOSResults' + cloudName + '.txt')
# -------- DEFINE FILES AND PATHS --------

# -------- READ BLOS DATA--------
BData = pd.read_csv(BFilePath)
# -------- READ BLOS DATA. --------

# -------- CREATE A TABLE FOR THE UNCERTAINTY DATA --------
cols = ['ID#', 'Ra(deg)', 'Dec(deg)', 'Extinction', 'Magnetic_Field(uG)',
        'TotalUpperBUncertainty', 'TotalLowerBUncertainty']
FinalBLOSResults = pd.DataFrame(columns=cols)
# -------- CREATE A TABLE FOR THE UNCERTAINTY DATA. --------

# -------- COPY OVER B DATA --------
# To round, use .round() eg:) FinalBLOSResults['Ra(deg)'] = BData['Ra(deg)'].round(2)
FinalBLOSResults['ID#'] = BData['ID#']
FinalBLOSResults['Ra(deg)'] = BData['Ra(deg)']
FinalBLOSResults['Dec(deg)'] = BData['Dec(deg)']
FinalBLOSResults['Extinction'] = BData['Extinction']
FinalBLOSResults['Magnetic_Field(uG)'] = BData['Magnetic_Field(uG)']
# -------- COPY OVER B DATA. --------

TotalRMErrStDevinB = (BData['Magnetic_Field(uG)']) * (BData['TotalRMScaledErrWithStDev'] /
                                                                         BData['Scaled_RM'])

# -------- CALCULATE UNCERTAINTIES --------
BTotalUpperUncertainty = []
BTotalLowerUncertainty = []

BChemDens50Increase = list(pd.read_csv(BData_Density50IncreasePath)['Magnetic_Field(uG)'])
BChemDens50Decrease = list(pd.read_csv(BData_Density50DecreasePath)['Magnetic_Field(uG)'])
BChemTemp20Increase = list(pd.read_csv(BData_Temp20IncreasePath)['Magnetic_Field(uG)'])
BChemTemp20Decrease = list(pd.read_csv(BData_Temp20DecreasePath)['Magnetic_Field(uG)'])

for index in range(len(BData)):
    upperDeltaBExt, lowerDeltaBExt = extinctionChemUncertainties(
        BData['Magnetic_Field(uG)'][index], BData['BField_of_Min_Extinction'][index],
        BData['BField_of_Max_Extinction'][index])

    upperDeltaBChemDens, lowerDeltaBChemDens = extinctionChemUncertainties(
        BData['Magnetic_Field(uG)'][index], BChemDens50Increase[index], BChemDens50Decrease[index])

    upperDeltaBChemTemp, lowerDeltaBChemTemp = extinctionChemUncertainties(
        BData['Magnetic_Field(uG)'][index], BChemTemp20Increase[index], BChemTemp20Decrease[index])

    # Calculate uncertainties
    BTotalUpperUncertainty.append("{0:.0f}".format(round(((TotalRMErrStDevinB[index]) ** 2 + upperDeltaBExt ** 2
                                                          + upperDeltaBChemDens ** 2 + upperDeltaBChemTemp ** 2) ** (
                                                                 1 / 2), 0)))

    BTotalLowerUncertainty.append("{0:.0f}".format(round(((TotalRMErrStDevinB[index]) ** 2 + lowerDeltaBExt ** 2
                                                          + lowerDeltaBChemDens ** 2 + lowerDeltaBChemTemp ** 2) ** (
                                                                 1 / 2), 0)))

FinalBLOSResults['TotalUpperBUncertainty'] = BTotalUpperUncertainty
FinalBLOSResults['TotalLowerBUncertainty'] = BTotalLowerUncertainty
# -------- CALCULATE UNCERTAINTIES. --------

# -------- SAVE FINAL BLOS RESULTS --------
FinalBLOSResults.to_csv(saveFilePath, index=False)
print('Saving calculated magnetic field values and associated uncertainties to '+saveFilePath)
# -------- SAVE FINAL BLOS RESULTS --------
