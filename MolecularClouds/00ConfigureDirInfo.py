import os
import math
from configparser import ConfigParser

# -------- DEFINE STARTING VARIABLES. --------
configStartSettings = ConfigParser()
configStartSettings['Cloud'] = {
    'Cloud': 'Oriona',
    }
configStartSettings['Judgement'] = {
    'Off Disk Latitude': 15.,
    'On Disk Extinction Threshold': 2.,
    'Off Disk Extinction Threshold': 1.,
    'Near High Extinction Multiplier': 2,
    'Far High Extinction Multiplier': 14,
    'High Extinction Threshold Multiplier': 5,
    'Anomalous Values Standard Deviation': 3.
    }
with open('configStartSettings.ini', 'w') as output_file:
    configStartSettings.write(output_file)
# -------- DEFINE STARTING VARIABLES. --------

# -------- DEFINE DIRECTORIES AND NAMES. --------
configDirectoryAndNames = ConfigParser()

configDirectoryAndNames['Output File Locations'] = {
    'Root': os.path.abspath(os.getcwd()),
    'File Output': 'FileOutput',
    'Plots': 'Plots',
    'Density Sensitivity': 'DensitySensitivity',
    'Temperature Sensitivity': 'TemperatureSensitivity'
    }

configDirectoryAndNames['Output File Prefixes'] = {
    'RM Mapping': 'RMMapping',
    'RM-Extinction Matching': 'MatchedRMExtinction',
    'All Potential Reference Points': 'AllPotentialRefPoints',
    'Selected Reference Points': 'RefPoints',
    'Reference Data': 'ReferenceData',
    'BLOS Point Data': 'BLOSPoints',
    'BLOS Point Figure': 'BLOSPointMap',
    'BLOS Uncertainties': 'FinalBLOSResults',
    'Optimal Reference Points': 'DataNoRef'
    }

configDirectoryAndNames['Input File Locations'] = {
    'Input Data': 'Data',
    'Cloud Parameter Data': 'CloudParameters',
    'Chemical Abundance Data': 'ChemicalAbundance'
    }

configDirectoryAndNames['Input File Names'] = {
    'RM Catalogue': 'RMCatalogue.txt'
    }

with open('configDirectoryAndNames.ini', 'w') as output_file:
    configDirectoryAndNames.write(output_file)
# -------- DEFINE DIRECTORIES AND NAMES. --------

# -------- DEFINE CONSTANTS. --------
configConstants = ConfigParser()

configConstants['Conversion Factors'] = {
    'Visual Extinction to Hydrogen Column Density': 2.21e21,
    'Parsecs to Centimeters': 3.24078e-19
    }

with open('configConstants.ini', 'w') as output_file:
    configConstants.write(output_file)
# -------- DEFINE CONSTANTS. --------

cloudInfoExport = ConfigParser()
cloudInfoExport['Cloud Info'] = {
    'distance': 0,
    'cloudJeansLength': 1,

    'fitsFileName': '',
    'fitsDataType': 'HydrogenColumnDensity',

    'xmin': math.nan,
    'xmax': math.nan,
    'ymin': math.nan,
    'ymax': math.nan,

    'n0': '',
    'T0': '',
    'G0': '',

    'cloudLatitude': 0,

    'raHoursMax': 0,
    'raMinsMax': 0,
    'raSecMax': 0,
    'raHoursMin': 0,
    'raMinsMin': 0,
    'raSecMin': 0,
    'decDegMax': 0,
    'decDegMin': 0,
}
with open('template.ini', 'w') as output_file:
    cloudInfoExport.write(output_file)