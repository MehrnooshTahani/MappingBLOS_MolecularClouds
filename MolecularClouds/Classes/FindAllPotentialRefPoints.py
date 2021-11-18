"""
This file contains the class to find "all potential reference points" for a given region of interest.
It is used in the third stage of the BLOS Mapping Method.
    - "All potential reference points" are defined to be all points within the region of interest with extinction values
    less than the specified threshold
"""
import pandas as pd
import numpy as np
import os


# -------- CLASS DEFINITION --------
class FindAllPotentialReferencePoints:
    def __init__(self, cloudName, threshold, saveFilePath='none'):
        """
          Takes matched extinction and rotation measure data for the region of interest
        and finds all of the potential reference points

        - All potential reference points have Av < threshold.  The threshold value defaults to 1.

        :param cloudName: Name of the region of interest
        :param threshold: Threshold extinction value.
        :param saveFilePath: Path to which output is saved. If 'none', no file will be saved
        """

        # -------- LOAD AND UNPACK MATCHED RM AND EXTINCTION DATA --------
        currentDir = os.path.abspath(os.getcwd())
        MatchedRMExtincPath = os.path.join(currentDir, 'FileOutput/' + cloudName + '/MatchedRMExtinction'
                                           + cloudName + '.txt')
        matchedRMExtinctionData = pd.read_csv(MatchedRMExtincPath, sep='\t')
        # -------- LOAD AND UNPACK MATCHED RM AND EXTINCTION DATA. --------

        # -------- FIND ALL POTENTIAL REFERENCE POINTS --------

        # ---- Criterion: Av < threshold
        '''
        We will only consider points with visual extinction less than the specified threshold value as potential 
        reference points
            - Here we extract these points and sort the resulting dataframe from smallest to greatest extinction 
        '''
        # Indices where the extinction is less than the threshold:
        ind_extinction = np.where(matchedRMExtinctionData['Extinction_Value'] <= threshold)[0]

        # All potential reference points are all reference points with extinction less than one:
        self.AllRefPoints = matchedRMExtinctionData.loc[ind_extinction].sort_values('Extinction_Value',
                                                                                    ignore_index=True)
        # ---- Criterion: Av < threshold.

        self.numAllRefPoints = len(self.AllRefPoints)
        # -------- FIND ALL POTENTIAL REFERENCE POINTS.  --------

        # -------- SAVE REFERENCE POINT DATA AS A TABLE --------
        if saveFilePath != 'none':
            self.AllRefPoints.to_csv(saveFilePath, index=False)
        # -------- SAVE REFERENCE POINT DATA AS A TABLE. --------

# -------- CLASS DEFINITION. --------
