"""
This is the second part of the sixth stage of the BLOSMapping method where the dependence on temperature is assessed.
    - In this part, the differences in the original BLOS and the BLOS calculated with varying electron abundances
    are plotted
"""
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from MolecularClouds.Classes.RegionOfInterest import Region
import Classes.config as config

# -------- CHOOSE THE REGION OF INTEREST --------
cloudName = input("Enter the name of the region of interest: ")
cloudName = cloudName.capitalize()  # Ensure only the first letter is capitalized
#cloudName = config.cloud
regionOfInterest = Region(cloudName)
# -------- CHOOSE THE REGION OF INTEREST. --------

# -------- DEFINE FILES AND PATHS --------
BScaledFileDir = os.path.join(config.dir_root, config.dir_fileOutput, cloudName, config.dir_temperatureSensitivity)
InitialPath = os.path.join(config.dir_root, config.dir_fileOutput, cloudName, config.dir_temperatureSensitivity + os.sep + 'B_Av_T0_n0.txt')
saveFigurePath = os.path.join(config.dir_root, config.dir_fileOutput, cloudName, config.dir_plots + os.sep + 'BTemperatureSensitivity.png')
# -------- DEFINE FILES AND PATHS. --------

# -------- EXTRACT ORIGINAL BLOS VALUES --------
InitialBData = pd.read_csv(InitialPath, sep='\t')
B = list(InitialBData['Magnetic_Field(uG)'])
# -------- EXTRACT ORIGINAL BLOS VALUES. --------

# -------- EXTRACT BLOS FOR EACH PERCENT OF THE INPUT DENSITY -------
p = [5, 10, 20]  # Percents of the input temperature
percent = ['-{}'.format(i) for i in p[::-1]] + ['0'] + ['+{}'.format(i) for i in p]

errPercent = []
errPercentFiles = []
for i, value in enumerate(percent):
    AvAbundanceName = 'Av_T' + value + '_n0'
    BScaledFilePath = BScaledFileDir + os.sep + 'B_' + AvAbundanceName + '.txt'
    try:
        BScaledTemp = list(pd.read_csv(BScaledFilePath, sep='\t')['Magnetic_Field(uG)'])
    except:
        errPercent.append(value)
        errPercentFiles.append(BScaledFilePath)
percent = [item for item in percent if item not in errPercent]

if len(errPercentFiles) > 0:
    print('-------------------------------------------------------------------------------')
    print('Warning: The following data have not been loaded due to an error.')
    print('{}'.format(errPercentFiles))
    print('Please review the results.')
    print('-------------------------------------------------------------------------------')

# Each row is a BLOS point, each column is the BLOS value corresponding to each percent of the input density
AllBScaled = np.zeros([len(B), len(percent)])

for i, value in enumerate(percent):
    AvAbundanceName = 'Av_T' + value + '_n0'
    BScaledFilePath = BScaledFileDir + os.sep + 'B_' + AvAbundanceName + '.txt'
    BScaledTemp = list(pd.read_csv(BScaledFilePath, sep='\t')['Magnetic_Field(uG)'])
    AllBScaled[:, i] = BScaledTemp[:]
    errPercentFiles.append(BScaledFilePath)

# -------- EXTRACT BLOS FOR EACH PERCENT OF THE INPUT DENSITY. -------

# -------- CHOOSE INDICES OF BLOS POINTS TO PLOT -------
numToPlot = int(input("Choose the number of points to plot: "))
AvMinToPlot = float(input("Choose the minimum extinction to plot: "))
AvMaxToPlot = float(input("Choose the maximum extinction to plot: "))

indMin = list(np.where(np.array(InitialBData['Scaled_Extinction']) >= AvMinToPlot)[0])
indMax = list(np.where(np.array(InitialBData['Scaled_Extinction']) <= AvMaxToPlot)[0])
ind = list(set(indMin).intersection(indMax))
if len(ind) < numToPlot:
    indToPlot = ind
else:
    indToPlot = ind[:numToPlot]
# -------- CHOOSE INDICES OF BLOS POINTS TO PLOT. -------

# -------- CREATE A FIGURE -------
fig = plt.figure(figsize=(12, 12), dpi=120, facecolor='w', edgecolor='k')

plt.ylabel('Magnetic Field Difference (' + r'$ \mu G$)', fontsize=16)
plt.xlabel(r'$\frac{\Delta T}{T_0} (\%)$', fontsize=16, labelpad=20)
plt.title('B$_{LOS}$ Variation, ' + cloudName + ', n = n$_0$', fontsize=16)

x = np.arange(0, len(percent))
plt.xticks(x, percent)

for i, ind in enumerate(indToPlot):
    barloc = [item + 0.15 * int(i) for item in x]
    plt.bar(barloc, B[ind] - AllBScaled[ind], width=0.25, edgecolor='k', label=indToPlot[i])

# ---- Style the legend
plt.legend(loc='upper center', ncol=2)
# ---- Style the legend.
# -------- CREATE A FIGURE. -------

# ---- Display or save the figure
# plt.show()
plt.savefig(saveFigurePath)
# ---- Display or save the figure.
print('Saving figure to '+saveFigurePath)
# -------- CREATE A FIGURE. --------
