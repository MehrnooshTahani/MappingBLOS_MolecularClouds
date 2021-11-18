__author__ = 'Mehrnoosh'

import csv
from itertools import zip_longest
import numpy as np


class RefPointsInfo():
    def __init__(self, ExtincRMPath, NRefPoints):

        Ra, RaErr, Dec, DecErr, RM, RMErr, Extinction, MinExtinctionRanged, MaxExtinctionRanged = \
            np.loadtxt(ExtincRMPath, usecols=( 2, 3, 4, 5, 6, 7, 12, 13, 16), unpack = True, skiprows=1)

        RaZoomed=[]
        DecZoomed=[]
        RMZoomed=[]
        ExtinctionZoomed=[]
        RMErrZoomed=[]
        MinExtinctionRangedZoomed=[]
        MaxExtinctionRangedZoomed=[]

        fiducialIndex=[]
        self.fiducialExt=[]
        self.fiducialRM=[]
        fiducialRMErr=[]
        self.finalFiducial=0

        #~~~~Function to find the next minimum value in the list~~~~
        def nextMin(lastMinIndex, extinc): #later define lastMin and lastMins as a list
            lastMinColDens = extinc[lastMinIndex]
            diff = 1
            nextMinIndex=0
            for index in range(len(extinc)):
                if abs(extinc[index] - lastMinColDens) <= diff and extinc[index] > lastMinColDens:
                    diff = abs(extinc[index] - lastMinColDens)
                    nextMinIndex = index
                    # print(index)
            return nextMinIndex

        #~~~~~~~~~~~~Finding the zoomed region~~~~~~~~~
        for indexZoomed in range(len(Ra)):
            #~~Perseus~~
            if Ra[indexZoomed] > 50 and Ra[indexZoomed] < 58 and Dec[indexZoomed] > 28 and Dec[indexZoomed]<34:
            #~~Perseus -side Filament~~
            # if Ra[indexZoomed] > 53 and Ra[indexZoomed] < 58 and Dec[indexZoomed] > 28 and Dec[indexZoomed]<34:
            #~~California~~
            # if Ra[indexZoomed] > 58 and Ra[indexZoomed] < 70 and Dec[indexZoomed] < 42 and Dec[indexZoomed]>34:
            #~~Real Taurus~~
            # if Ra[indexZoomed] > 62 and Ra[indexZoomed] < 72 and Dec[indexZoomed] > 21.9 and Dec[indexZoomed]<30:
            #~~~Bottom Bottom Filament~~~~~
            # if Ra[indexZoomed] > 62 and Ra[indexZoomed] < 68 and Dec[indexZoomed] > 25.5 and Dec[indexZoomed]<30:
            #~~~Bottom Top Filament~~~
            # if Ra[indexZoomed] > 62 and Ra[indexZoomed] < 70 and Dec[indexZoomed] > 21.9 and Dec[indexZoomed]<26:
            #~~~Bottom Left Filament~~~
            # if Ra[indexZoomed] > 68.5 and Ra[indexZoomed] < 74 and Dec[indexZoomed] > 24 and Dec[indexZoomed]<28:
            #~~~Jouni Filament~~
            # if Ra[indexZoomed] > 64 and Ra[indexZoomed] < 76 and Dec[indexZoomed] > 25 and Dec[indexZoomed]<30:
            #~~~Left Part~~~
            # if Ra[indexZoomed] > 67.5 and Ra[indexZoomed] < 72 and Dec[indexZoomed] < 37 and Dec[indexZoomed]>30:
                RaZoomed.append(Ra[indexZoomed])
                DecZoomed.append(Dec[indexZoomed])
                RMZoomed.append(RM[indexZoomed])
                RMErrZoomed.append(RMErr[indexZoomed])
                ExtinctionZoomed.append(Extinction[indexZoomed])
                MinExtinctionRangedZoomed.append(MinExtinctionRanged[indexZoomed])
                MaxExtinctionRangedZoomed.append(MaxExtinctionRanged[indexZoomed])

        #~~~~~~~~~~Finding the Reference point for RM based on lowest Column density~~~~~~~~~~~
        minColDensitIndex = ExtinctionZoomed.index(min(ExtinctionZoomed))
        fiducialIndex.append(minColDensitIndex)
        self.fiducialRM.append(RMZoomed[minColDensitIndex])
        self.fiducialExt.append(ExtinctionZoomed[minColDensitIndex])
        self.finalFiducial = self.fiducialRM[0]#0 for taking the first ref out
        self.finalFiducialExt = self.fiducialExt[0]

        #~~~Second Reference point~~~~~
        for index in range(1, NRefPoints):
            fiducialIndex.append(nextMin(fiducialIndex[index -1], ExtinctionZoomed))
            self.fiducialRM.append(RMZoomed[fiducialIndex[index]])
            self.fiducialExt.append(ExtinctionZoomed[fiducialIndex[index]])
            self.finalFiducial = self.finalFiducial + self.fiducialRM[index]
            self.finalFiducialExt= self.finalFiducialExt + self.fiducialExt[index]

        #~~~~~~for taking the first ref out~~~~~~~
        # if NRefPoints!= 1:
        self.finalFiducial = self.finalFiducial/(NRefPoints) # NRefPoint for all, and NRefPoint-1  for excluding first ref point.
        self.finalFiducialExt = self.finalFiducialExt/(NRefPoints)
        # else:
        #     finalFiducial = finalFiducial

        # for index in range(len(fiducialIndex)):
        #     print(index, "new point",fiducialIndex[index], "coordinate", RaZoomed[fiducialIndex[index]],
        #           DecZoomed[fiducialIndex[index]], HColDensityZoomed[fiducialIndex[index]])
        #~~Finding Standard Deviation of Reference Points~~
        if NRefPoints != 1:
            FidStandardDevTemp=0
            for index in range(1, NRefPoints): # Standard dev for 10 ref points
                FidStandardDevTemp = (self.fiducialRM[index] - self.finalFiducial)**2 + FidStandardDevTemp
            self.FinalFidStandardDev = (FidStandardDevTemp/((NRefPoints-1)*NRefPoints))**(0.5)

