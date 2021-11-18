"""
This file contains the class to calculate BLOS values.

"""
import pandas as pd
import numpy as np
import config


# -------- CLASS DEFINITION --------
class CalculateB:
    def __init__(self, AvAbundancePath, ExtincRMPath, RefPointTable, saveFilePath):
        """
        Takes files containing extinction, rotation measure data, and reference point data for the region of interest
        and calculates BLOS

        :param AvAbundancePath:  Path to extinction data produced by chemical evolution code
        :param ExtincRMPath: Path to matched extinction and rotation measure data (produced in stage 02)
        :param RefPointTable: Table (pandas dataframe)  of potential reference points
        :param saveFilePath: Path to which output is saved
        """

        conversionFactor = config.VExtinct_2_Hcol  # to convert extinction to H column density
        pcTocm = config.pcTocm

        # -------- LOAD REFERENCE POINTS --------
        refData = RefPointTable
        # -------- LOAD REFERENCE POINTS. --------

        # -------- LOAD MATCHED RM AND EXTINCTION DATA
        AllMatchedRMExtinctionData = pd.read_csv(ExtincRMPath, sep='\t')
        # -------- LOAD MATCHED RM AND EXTINCTION DATA.

        # -------- LOAD ABUNDANCE DATA
        Av, eAbundance = np.loadtxt(AvAbundancePath, usecols=(1, 2), unpack=True, skiprows=2)
        # -------- LOAD ABUNDANCE DATA.

        # -------- FIND FIDUCIAL REFERENCE VALUES --------
        self.fiducialRM = np.mean(refData['Rotation_Measure(rad/m2)'])
        self.fiducialRMAvgErr = np.mean(refData['RM_Err(rad/m2)'])
        # Standard error of the sampled mean:
        self.fiducialRMStd = np.std(refData['Rotation_Measure(rad/m2)'], ddof=1) / np.sqrt(
            len(refData['Rotation_Measure(rad/m2)']))
        self.fiducialExtinction = np.mean(refData['Extinction_Value'])
        # -------- FIND FIDUCIAL REFERENCE VALUES. --------

        # -------- REMOVE REFERENCE POINTS FROM THE MATCHED RM AND EXTINCTION DATA --------
        # The rm points used as reference points should not be used to calculate BLOS
        ind = refData['ID#']  # Indices of the reference points
        RMExtinctionData = AllMatchedRMExtinctionData.copy().drop(ind).reset_index(drop=True)
        # -------- REMOVE REFERENCE POINTS FROM THE MATCHED RM AND EXTINCTION DATA. --------

        # -------- CREATE BLOS TABLE --------
        cols = ['ID#', 'Ra(deg)', 'Dec(deg)', 'RM_Raw_Value', 'RM_Raw_Err', 'Scaled_RM', 'TotalRMScaledErrWithStDev',
                'TotalRMScaledErrWithAvgErr', 'Extinction', 'Scaled_Extinction', 'eAbundance',
                'BScaled_RM_ERR_with_0.05Uncty&StDev', 'BScaled_RM_ERR_with_0.05Uncty',
                'Raw_Magnetic_FieldMagnetic_Field(uG)', 'Magnetic_Field(uG)',
                'Reference_BField_RMErr(\u00B1)', 'BField_of_Min_Extinction', 'BField_of_Max_Extinction']
        self.BLOSData = pd.DataFrame(columns=cols)
        # -------- CREATE BLOS TABLE. --------

        # -------- ADD TO BLOS TABLE --------
        self.BLOSData['Ra(deg)'] = RMExtinctionData['Ra(deg)']
        self.BLOSData['Dec(deg)'] = RMExtinctionData['Dec(deg)']
        self.BLOSData['ID#'] = RMExtinctionData[
            'ID#']  # Want to keep track of the 'original' index as well
        self.BLOSData['RM_Raw_Value'] = RMExtinctionData['Rotation_Measure(rad/m2)']
        self.BLOSData['RM_Raw_Err'] = RMExtinctionData['RM_Err(rad/m2)']
        self.BLOSData['Extinction'] = RMExtinctionData['Extinction_Value']
        # -------- ADD TO BLOS TABLE. --------

        # -------- SCALE THE RM AND EXTINCTION DATA --------
        self.BLOSData['Scaled_RM'] = RMExtinctionData['Rotation_Measure(rad/m2)'] - self.fiducialRM

        self.BLOSData['Scaled_Extinction'] = RMExtinctionData['Extinction_Value'] - self.fiducialExtinction
        Scaled_Min_Extinction_Value = RMExtinctionData['Min_Extinction_Value'] - self.fiducialExtinction
        Scaled_Max_Extinction_Value = RMExtinctionData['Max_Extinction_Value'] - self.fiducialExtinction
        # -------- SCALE THE RM AND EXTINCTION DATA. --------

        # -------- CALCULATE THE RM ERROR --------
        self.BLOSData['TotalRMScaledErrWithStDev'] = RMExtinctionData['RM_Err(rad/m2)'] + self.fiducialRMStd
        self.BLOSData['TotalRMScaledErrWithAvgErr'] = RMExtinctionData['RM_Err(rad/m2)'] + self.fiducialRMAvgErr
        # -------- CALCULATE THE RM ERROR. --------

        # -------- FIND THE LAYER OF INTEREST --------
        eAbundanceMatched = []
        indLayerOfInterest = []
        indLayerOfInterest_MinExt = []
        indLayerOfInterest_MaxExt = []

        # For each BLOS point:
        for i in range(0, len(self.BLOSData['Scaled_Extinction'])):
            '''We want to find the layer that is closest in value and greater than the half the scaled 
             extinction value. Since Av is a list ordered from least to greatest, this corresponds to the first location 
             where Av is greater than half the scaled extinction value 
             '''

            ind = np.where(Av >= self.BLOSData['Scaled_Extinction'][i] / 2)[0][0]
            indLayerOfInterest.append(ind)
            eAbundanceMatched.append(eAbundance[ind])

            indMin = np.where(Av >= Scaled_Min_Extinction_Value[i] / 2)[0][0]
            indLayerOfInterest_MinExt.append(indMin)

            indMax = np.where(Av >= Scaled_Max_Extinction_Value[i] / 2)[0][0]
            indLayerOfInterest_MaxExt.append(indMax)

        self.BLOSData['eAbundance'] = eAbundanceMatched
        # -------- FIND THE LAYER OF INTEREST. --------

        # -------- CALCULATE THE TOTAL ELECTRON COLUMN DENSITY --------
        LayerNe = []
        LayerNeMinExt = []
        LayerNeMaxExt = []

        # -------- Matched Extinction Value
        # For each BLOS point:
        for i, val in enumerate(indLayerOfInterest):
            # Temporary value
            tempSumAvSubXe = 0
            if val != 0:
                for index3 in range(1, val):
                    tempSumAvSubXe = tempSumAvSubXe + ((Av[index3] - Av[index3 - 1]) * eAbundance[index3])
                # Interpolate:
                xp = [Av[val], Av[val - 1]]
                fp = [eAbundance[val], eAbundance[val - 1]]
                interpAv = (self.BLOSData['Scaled_Extinction'][i] / 2)
                interpEAbund = np.interp(interpAv, xp, fp)
                tempSumAvSubXe = tempSumAvSubXe + ((Av[0]) * eAbundance[0]) + (interpAv - Av[val - 1]) * interpEAbund
            else:
                tempSumAvSubXe = tempSumAvSubXe + ((Av[0]) * eAbundance[0])
            LayerNe.append(tempSumAvSubXe * conversionFactor)
        # -------- Matched Extinction Value.

        # -------- Minimum Extinction Value
        # For each BLOS point:
        for i, val in enumerate(indLayerOfInterest_MinExt):
            # Temporary value
            tempSumAvSubXe = 0
            if val != 0:
                for index3 in range(1, val):
                    tempSumAvSubXe = tempSumAvSubXe + ((Av[index3] - Av[index3 - 1]) * eAbundance[index3])
                # Interpolate:
                xp = [Av[val], Av[val - 1]]
                fp = [eAbundance[val], eAbundance[val - 1]]
                interpAv = (Scaled_Min_Extinction_Value[i] / 2)
                interpEAbund = np.interp(interpAv, xp, fp)
                tempSumAvSubXe = tempSumAvSubXe + ((Av[0]) * eAbundance[0]) + (interpAv - Av[val - 1]) * interpEAbund
            else:
                tempSumAvSubXe = tempSumAvSubXe + ((Av[0]) * eAbundance[0])
            LayerNeMinExt.append(tempSumAvSubXe * conversionFactor)
        # -------- Minimum Extinction Value.

        # -------- Maximum Extinction Value
        # For each BLOS point:
        for i, val in enumerate(indLayerOfInterest_MaxExt):
            # Temporary value
            tempSumAvSubXe = 0
            if val != 0:
                for index3 in range(1, val):
                    tempSumAvSubXe = tempSumAvSubXe + ((Av[index3] - Av[index3 - 1]) * eAbundance[index3])
                # Interpolate:
                xp = [Av[val], Av[val - 1]]
                fp = [eAbundance[val], eAbundance[val - 1]]
                interpAv = (Scaled_Max_Extinction_Value[i] / 2)
                interpEAbund = np.interp(interpAv, xp, fp)
                tempSumAvSubXe = tempSumAvSubXe + ((Av[0]) * eAbundance[0]) + (interpAv - Av[val - 1]) * interpEAbund
            else:
                tempSumAvSubXe = tempSumAvSubXe + ((Av[0]) * eAbundance[0])
            LayerNeMaxExt.append(tempSumAvSubXe * conversionFactor)
        # -------- Maximum Extinction Value.
        # -------- CALCULATE THE TOTAL ELECTRON COLUMN DENSITY. -------

        # -------- CALCULATE THE MAGNETIC FIELD --------
        self.BLOSData['Raw_Magnetic_FieldMagnetic_Field(uG)'] = self.BLOSData['RM_Raw_Value'] / (
                0.812 * np.array(LayerNe) * pcTocm * 2)

        self.BLOSData['Magnetic_Field(uG)'] = self.BLOSData['Scaled_RM'] / (
                0.812 * np.array(LayerNe) * pcTocm * 2)

        self.BLOSData['Reference_BField_RMErr(\u00B1)'] = (self.BLOSData['RM_Raw_Err'] / self.BLOSData['RM_Raw_Value']) \
                                                     * self.BLOSData['Magnetic_Field(uG)']

        self.BLOSData['BField_of_Min_Extinction'] = self.BLOSData['Scaled_RM'] / (0.812 * np.array(LayerNeMinExt) * pcTocm * 2)
        self.BLOSData['BField_of_Max_Extinction'] = self.BLOSData['Scaled_RM'] / (0.812 * np.array(LayerNeMaxExt) * pcTocm * 2)

        self.BLOSData['BScaled_RM_ERR_with_0.05Uncty&StDev'] = ((0.05 + self.fiducialRMStd)
                                                           * self.BLOSData['Magnetic_Field(uG)']) / \
                                                          self.BLOSData['Scaled_RM']

        self.BLOSData['BScaled_RM_ERR_with_0.05Uncty'] = (0.1 * self.BLOSData['Magnetic_Field(uG)']) \
                                                    / self.BLOSData['Scaled_RM']
        # -------- CALCULATE THE MAGNETIC FIELD. --------

        # -------- SAVE BLOS DATA --------
        if saveFilePath != 'none':
            self.BLOSData.to_csv(saveFilePath, index=False)
        # -------- SAVE BLOS DATA. --------
