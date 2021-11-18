"""
This file contains the class which calculates BLOS

THIS IS OLD
"""
import csv
from itertools import zip_longest
import numpy as np
import pandas as pd


# -------- CLASS DEFINITION --------
class ZoomedScalingXeAverage:
    def __init__(self, AvAbundancePath, ExtincRMPath, NRefPoints, saveFilePath,):
        """
        Takes files containing extinction and rotation measure data for the region of interest and calculates BLOS

        :param AvAbundancePath: Path to extinction data produced by chemical evolution code
        :param ExtincRMPath:  Path to matched extinction and rotation measure data (produced in stage 02)
        :param NRefPoints:  Number of reference points
        :param saveFilePath:  Path to which output is saved
        """

        # -------- LOAD ABUNDANCE DATA
        Av, eAbundance = np.loadtxt(AvAbundancePath, usecols=(1, 2), unpack=True, skiprows=2)
        # -------- LOAD ABUNDANCE DATA.

        # -------- LOAD AND UNPACK MATCHED RM AND EXTINCTION DATA
        matchedRMExtinctionData = pd.read_csv(ExtincRMPath, sep='\t')
        '''THE FOLLOWING WILL HAVE INCORRECT COLUMNS - MUST CHANGE '''
        Ra = list(matchedRMExtinctionData[matchedRMExtinctionData.columns[2]])
        Dec = list(matchedRMExtinctionData[matchedRMExtinctionData.columns[3]])
        RM = list(matchedRMExtinctionData[matchedRMExtinctionData.columns[4]])
        RMErr = list(matchedRMExtinctionData[matchedRMExtinctionData.columns[5]])
        Extinction = list(matchedRMExtinctionData[matchedRMExtinctionData.columns[8]])
        MinExtinctionRanged = list(matchedRMExtinctionData[matchedRMExtinctionData.columns[10]])
        MaxExtinctionRanged = list(matchedRMExtinctionData[matchedRMExtinctionData.columns[13]])
        # -------- LOAD AND UNPACK MATCHED RM AND EXTINCTION DATA.

        MagneticField = []
        ScaledRM = []
        ScaledExt = []
        ScaledMinExt = []
        ScaledMaxExt = []
        MagneticFieldRMScaled = []

        LayerNe = []
        LayerNeMinExt = []
        LayerNeMaxExt = []
        BScaledRMERRW05UnctyStDev = []
        BScaledRMERRW05Uncty = []

        ConversionFactor = 2.21 * (10 ** 21)  # to convert extinction to H column density
        pcTocm = 3.24078e-19

        self.fiducialIndex = []
        self.fiducialExt = []
        self.fiducialRM = []
        self.fiducialRMErr = []
        self.finalFiducial = 0  # The reference rotation measure ( RM_{ref} )  for the region of interest
        self.finalFiducialExt = 0  # The reference extinction ( A_{v,ref} ) for the region of interest
        self.finalRMRefErrAvg = 0  # The reference rotation measure error  for the region of interest
        self.FinalFidStandardDev = 0

        TotalRMScaledErrWithStDev = []
        TotalRMScaledErrWithAvgErr = []

        # Errors in Magnetic field values with only taking the RM errors into account:
        MagneticFieldScaledRMErr = []
        # The value of Scaled Magnetic Field by using the min Extinction in the uncertainty location:
        MagneticFieldScaledMinExtinc = []
        # The value of Scaled Magnetic Field by using the max Extinction in the uncertainty location:
        MagneticFieldScaledMaxExtinc = []
        eAbundMatched = []

        indLayerOfInterest = []  # REMOVE LATERE
        test = [] # Remove later

        # -------- FUNCTION DEFINITION --------
        def nextMin(lastMinIndex, extinc):  # later define lastMin and lastMins as a list
            """
            Function to find the next minimum value in the list
            :param lastMinIndex: The ind of the previous minimum (or next greatest) extinction value
            :param extinc: The list of extinction values
            :return: Index of the next minimum value in the list
            """
            lastMinColDens = extinc[lastMinIndex]
            diff = 2
            nextMinIndex = 0
            for ind in range(len(extinc)):
                if abs(extinc[ind] - lastMinColDens) <= diff and extinc[ind] > lastMinColDens:
                    diff = abs(extinc[ind] - lastMinColDens)
                    nextMinIndex = ind

            return nextMinIndex
        # -------- FUNCTION DEFINITION. --------

        # -------- FIND THE FIRST REFERENCE POINT --------
        # Find the Reference point for RM based on lowest Column density (i.e extinction )
        minColDensitIndex = Extinction.index(min(Extinction))
        self.fiducialIndex.append(minColDensitIndex)
        self.fiducialRM.append(RM[minColDensitIndex])
        self.fiducialRMErr.append(RMErr[minColDensitIndex])
        self.fiducialExt.append(Extinction[minColDensitIndex])
        self.finalFiducial = self.fiducialRM[0]  # 0 for taking the first ref out
        self.finalFiducialExt = self.fiducialExt[0]
        self.finalRMRefErrAvg = self.fiducialRMErr[0]
        # -------- FIND THE FIRST REFERENCE POINT. --------

        # -------- FIND THE NEXT REFERENCE POINTS --------
        for index in range(1, NRefPoints):
            self.fiducialIndex.append(nextMin(self.fiducialIndex[index - 1], Extinction))
            self.fiducialRM.append(RM[self.fiducialIndex[index]])
            self.fiducialRMErr.append(RMErr[self.fiducialIndex[index]])
            self.fiducialExt.append(Extinction[self.fiducialIndex[index]])
            self.finalFiducial = self.finalFiducial + self.fiducialRM[index]
            self.finalFiducialExt = self.finalFiducialExt + self.fiducialExt[index]
            self.finalRMRefErrAvg = self.finalRMRefErrAvg + self.fiducialRMErr[index]
        # -------- FIND THE NEXT REFERENCE POINTS. --------

        # -------- FIND THE AVERAGE REFERENCE POINTS (RM_REF) --------
        self.finalFiducial = self.finalFiducial / NRefPoints
        self.finalFiducialExt = self.finalFiducialExt / NRefPoints
        self.finalRMRefErrAvg = self.finalRMRefErrAvg / NRefPoints
        # -------- FIND THE AVERAGE REFERENCE POINTS (RM_REF). --------

        # -------- FIND THE STANDARD DEVIATION IN REFERENCE POINTS -------
        if NRefPoints != 1:
            FidStandardDevTemp = 0
            for index in range(NRefPoints):
                FidStandardDevTemp = (self.fiducialRM[index] - self.finalFiducial) ** 2 + FidStandardDevTemp
            # Standard error of the sampled mean:
            self.FinalFidStandardDev = (FidStandardDevTemp / ((NRefPoints - 1) * NRefPoints)) ** 0.5
            # Standard deviation of the population:
            # self.FinalFidStandardDev = (FidStandardDevTemp / ( (NRefPoints ) )) ** 0.5
        else:
            self.FinalFidStandardDev = 0
        # -------- FIND THE STANDARD DEVIATION IN REFERENCE POINTS. --------

        # For each rotation measure (and its matched extinction value) within the region of interest:
        for index in range(len(Ra)):

            # -------- SCALE THE DATA --------
            # The rotation measure and extinction values must be corrected for to account for the galactic contribution

            # Subtract the reference rm from the rm value:
            ScaledRM.append(RM[index] - self.finalFiducial)

            # Subtract the reference extinction from the extinction value:
            ScaledExt.append(Extinction[index] - self.finalFiducialExt)
            ScaledMinExt.append(MinExtinctionRanged[index] - self.finalFiducialExt)
            ScaledMaxExt.append(MaxExtinctionRanged[index] - self.finalFiducialExt)

            # To find the average Error of the RM of the Off Points:
            TotalRMScaledErrWithStDev.append(RMErr[index] + self.FinalFidStandardDev)
            # To find the standard deviation of the Ref Points:
            TotalRMScaledErrWithAvgErr.append(RMErr[index] + self.finalRMRefErrAvg)
            # -------- SCALE THE DATA. --------

            '''
            For each rotation measure we must determine a corresponding total electron column density.
            
                - In file 02, rotation measure values were matched with visual extinction values.
                - The chemical evolution code provides electron abundance as a function of cloud layer
            
            To calculate the electron column density corresponding to a given rotation measure with a given visual 
            extinction, we must account for the contribution to the total electron column density from each layer of the 
            cloud.
            
            First we match half of the extinction value from the visual extinction map to the extinction value from the 
            chemical evolution code.  This allows us to find the index corresponding to the "layer of interest" 
            
            We then calculate the total electron column density by summing over all layers up to the layer of interest.
            To account for any discrepancy between the matched visual extinction and the visual extinction corresponding 
            to the layer of interest, the final value is extrapolated
            
            In order to estimate errors, this process is repeated with the minimum and maximum extinction values 
            determined in file 02
            '''

            # -------- FIND THE LAYER OF INTEREST --------
            eAbundMatched.append(1)
            # to find the match in chem e abundance
            minAvDifference = 1
            minAvDifferenceMinExt = 1  # for min extinction value in the uncertainty region
            minAvDifferenceMaxExt = 1  # for max extinction value in the uncertainty region

            minIndex = 0
            minIndexMinExt = 0
            minIndexMaxExt = 0

            # For each extinction value and its corresponding electron column density in the abundance file:
            for index2 in range(len(Av)):

                if abs(((ScaledExt[index]) / 2) - Av[index2]) <= minAvDifference and Av[index2] >= (
                        (ScaledExt[index]) / 2):
                    minAvDifference = abs(((ScaledExt[index]) / 2) - Av[index2])
                    minIndex = index2  # Index of the layer of interest

                if abs(((ScaledMinExt[index]) / 2) - Av[index2]) <= minAvDifferenceMinExt and Av[index2] >= (
                        (ScaledMinExt[index]) / 2):
                    minAvDifferenceMinExt = abs(((ScaledMinExt[index]) / 2) - Av[index2])
                    minIndexMinExt = index2

                if abs(((ScaledMaxExt[index]) / 2) - Av[index2]) <= minAvDifferenceMaxExt and Av[index2] >= (
                        (ScaledMaxExt[index]) / 2):
                    minAvDifferenceMaxExt = abs(((ScaledMaxExt[index]) / 2) - Av[index2])
                    minIndexMaxExt = index2
            # -------- FIND THE LAYER OF INTEREST. --------

            indLayerOfInterest.append(minIndex) #  REMOVE LATERE

            # -------- CALCULATE THE TOTAL ELECTRON COLUMN DENSITY --------
            # Temporary values
            tempSumAvSubXe = 0
            tempSumAvSubXeMinExt = 0
            tempSumAvSubXeMaxExt = 0

            # Matched Extinction Value
            if minIndex != 0:
                for index3 in range(1, minIndex):
                    tempSumAvSubXe = tempSumAvSubXe + ((Av[index3] - Av[index3 - 1]) * eAbundance[index3])

                # Interpolate:
                xp = [Av[minIndex], Av[minIndex - 1]]
                fp = [eAbundance[minIndex], eAbundance[minIndex - 1]]
                interpAv = (ScaledExt[index] / 2)
                interpEAbund = np.interp(interpAv, xp, fp)

                tempSumAvSubXe = tempSumAvSubXe + ((Av[0]) * eAbundance[0]) + (
                        interpAv - Av[minIndex - 1]) * interpEAbund
            else:
                tempSumAvSubXe = tempSumAvSubXe + ((Av[0]) * eAbundance[0])

            # Minimum Extinction Value
            if minIndexMinExt != 0:
                for index5 in range(1, minIndexMinExt):
                    tempSumAvSubXeMinExt = tempSumAvSubXeMinExt + ((Av[index5] - Av[index5 - 1]) * eAbundance[index5])

                # Interpolate:
                xpMinExt = [Av[minIndexMinExt], Av[minIndexMinExt - 1]]
                fpMinExt = [eAbundance[minIndexMinExt], eAbundance[minIndexMinExt - 1]]
                interpAvMinExt = (ScaledMinExt[index] / 2)
                interpEAbundMinExt = np.interp(interpAvMinExt, xpMinExt, fpMinExt)

                tempSumAvSubXeMinExt = tempSumAvSubXeMinExt + ((Av[0]) * eAbundance[0]) + (
                        interpAvMinExt - Av[minIndexMinExt - 1]) * interpEAbundMinExt
            else:
                tempSumAvSubXeMinExt = tempSumAvSubXeMinExt + ((Av[0]) * eAbundance[0])

            # Maximum Extinction Value
            if minIndexMaxExt != 0:
                for index6 in range(1, minIndexMaxExt):
                    tempSumAvSubXeMaxExt = tempSumAvSubXeMaxExt + ((Av[index6] - Av[index6 - 1]) * eAbundance[index6])

                # Interpolate:
                xpMaxExt = [Av[minIndexMaxExt], Av[minIndexMaxExt - 1]]
                fpMaxExt = [eAbundance[minIndexMaxExt], eAbundance[minIndexMaxExt - 1]]
                interpAvMaxExt = ScaledMaxExt[index] / 2
                interpEAbundMaxExt = np.interp(interpAvMaxExt, xpMaxExt, fpMaxExt)

                tempSumAvSubXeMaxExt = tempSumAvSubXeMaxExt + ((Av[0]) * eAbundance[0]) + (
                        interpAvMaxExt - Av[minIndexMaxExt - 1]) * interpEAbundMaxExt
            else:
                tempSumAvSubXeMaxExt = tempSumAvSubXeMaxExt + ((Av[0]) * eAbundance[0])

            LayerNe.append(tempSumAvSubXe * ConversionFactor)
            test.append((LayerNe[index]))  # REMOVE LATER
            LayerNeMinExt.append(tempSumAvSubXeMinExt * ConversionFactor)
            LayerNeMaxExt.append(tempSumAvSubXeMaxExt * ConversionFactor)
            # -------- CALCULATE THE TOTAL ELECTRON COLUMN DENSITY. --------

            # -------- CALCULATE THE MAGNETIC FIELD --------
            MagneticField.append(RM[index] / (0.812 * LayerNe[index] * pcTocm * 2))
            MagneticFieldRMScaled.append(ScaledRM[index] / (0.812 * LayerNe[index] * pcTocm * 2))
            MagneticFieldScaledRMErr.append((RMErr[index] / RM[index]) * MagneticFieldRMScaled[index])
            MagneticFieldScaledMinExtinc.append(ScaledRM[index] / (0.812 * LayerNeMinExt[index] * pcTocm * 2))
            MagneticFieldScaledMaxExtinc.append(ScaledRM[index] / (0.812 * LayerNeMaxExt[index] * pcTocm * 2))

            BScaledRMERRW05UnctyStDev.append(
                ((0.05 + self.FinalFidStandardDev) * MagneticFieldRMScaled[index]) / ScaledRM[index])

            BScaledRMERRW05Uncty.append((0.1 * MagneticFieldRMScaled[index]) / ScaledRM[index])
            # -------- CALCULATE THE MAGNETIC FIELD. --------

            # -------- REMOVING REFERENCE POINTS --------
            for index6 in range(NRefPoints):
                if index == self.fiducialIndex[index6]:
                    MagneticFieldRMScaled[index] = 0
                    MagneticFieldScaledMaxExtinc[index] = 0
                    MagneticFieldScaledMinExtinc[index] = 0
                    MagneticFieldScaledRMErr[index] = 0
                    BScaledRMERRW05UnctyStDev[index] = 0
                    BScaledRMERRW05Uncty[index] = 0
            # -------- REMOVING REFERENCE POINTS. --------


        # -------- WRITE TO A FILE --------
        with open(saveFilePath, 'w') as f:
            f.write('Ra' + '\t' + 'Dec' + '\t' + 'RM_Raw_Value' + '\t' + 'RM_Raw_Err' + '\t' + 'Scaled_RM' + '\t' +
                    'TotalRMScaledErrWithStDev' + '\t' + 'TotalRMScaledErrWithAvgErr' + '\t' + 'Extinction' + '\t' +
                    'Scaled_Extinction' + '\t' + 'eAbundance' + '\t' + 'BScaled_RM_ERR_with_0.05Uncty&StDev' + '\t' +
                    'BScaled_RM_ERR_with_0.05Uncty' + '\t' + 'Raw_Magnetic_FieldMagnetic_Field(uG)' + '\t' +
                    'Magnetic_Field(uG)' + '\t' + 'Reference_BField_RMErr(\u00B1)' + '\t' +
                    'BField_of_Min_Extinction' + '\t' + 'BField_of_Max_Extinction' + '\t' + 'indLayerOfInterest' +'\t' + 'test' + '\n')
            writer = csv.writer(f, delimiter='\t')
            writer.writerows(zip_longest(Ra, Dec, RM, RMErr, ScaledRM,
                                         TotalRMScaledErrWithStDev, TotalRMScaledErrWithAvgErr, Extinction,
                                         ScaledExt, eAbundMatched, BScaledRMERRW05UnctyStDev, BScaledRMERRW05Uncty,
                                         MagneticField, MagneticFieldRMScaled, MagneticFieldScaledRMErr,
                                         MagneticFieldScaledMinExtinc, MagneticFieldScaledMaxExtinc,indLayerOfInterest,test, fillvalue=''))
            f.close()
        # -------- WRITE TO A FILE. --------

# -------- CLASS DEFINITION. --------
