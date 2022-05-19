"""
This is the zeroth stage of the BLOSMapping method where the necessary directories are made

    - When run, the file will ask for the name of a region of interest.  It will then check to see of this region of
    interest has a folder yet, if not it will make the needed folders and sub-folders.
"""
import os
import LocalLibraries.config as config

# -------- CHOOSE THE REGION OF INTEREST --------
cloudName = config.cloud
# -------- CHOOSE THE REGION OF INTEREST. --------

# -------- MAKE DIRECTORIES FOR THE REGION OF INTEREST --------
# Check if there is a FileOutput directory already; if not make one
#   - this will house the results for all of regions of interest
os.chdir(config.dir_root)
InHome = os.listdir(config.dir_root)
if config.dir_fileOutput not in InHome:
    os.mkdir(config.dir_fileOutput)

# Move into the FileOutput directory
os.chdir(config.dir_fileOutput)

# Check if there is a directory for the specified region of interest already; if not make one
InFileOutput = os.listdir()
InFileOutput = [item.lower() for item in InFileOutput]
if cloudName.lower() not in InFileOutput:
    os.mkdir(cloudName)
    # Move into the directory for the region of interest and make subsequent directories to house its results
    os.chdir(cloudName)
    os.mkdir(config.dir_plots)
    os.mkdir(config.dir_densitySensitivity)
    os.mkdir(config.dir_temperatureSensitivity)

print('Folder \'' + cloudName + '\' with sub-folders: \'' + config.dir_plots + '\', \'' + config.dir_densitySensitivity + '\', and \'' + config.dir_temperatureSensitivity + '\''
                            ' created in {}'.format(os.path.join(config.dir_root, config.dir_fileOutput)))
# -------- MAKE DIRECTORIES FOR THE REGION OF INTEREST. --------
