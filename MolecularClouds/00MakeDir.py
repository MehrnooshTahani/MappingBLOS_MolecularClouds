"""
This is the zeroth stage of the BLOSMapping method where the necessary directories are made

    - When run, the file will ask for the name of a region of interest.  It will then check to see of this region of
    interest has a folder yet, if not it will make the needed folders and sub-folders.
"""
import os

currentDir = os.path.abspath(os.getcwd())

# -------- CHOOSE THE REGION OF INTEREST --------
cloudName = input("Enter the name of the region of interest: ")
cloudName = cloudName.capitalize()  # Ensure only the first letter is capitalized
# -------- CHOOSE THE REGION OF INTEREST. --------

# -------- MAKE DIRECTORIES FOR THE REGION OF INTEREST --------
# Check if there is a FileOutput directory already; if not make one
#   - this will house the results for all of regions of interest
InHome = os.listdir()
if 'FileOutput' not in InHome:
    os.mkdir('FileOutput')

# Move into the FileOutput directory
os.chdir('FileOutput')

# Check if there is a directory for the specified region of interest already; if not make one
InFileOutput = os.listdir()
InFileOutput = [item.lower() for item in InFileOutput]
if cloudName.lower() not in InFileOutput:
    os.mkdir(cloudName)
    # Move into the directory for the region of interest and make subsequent directories to house its results
    os.chdir(cloudName)
    os.mkdir('Plots')
    os.mkdir('DensitySensitivity')
    os.mkdir('TemperatureSensitivity')

print('Folder \''+cloudName+'\' with sub-folders: \'Plots\', \'DensitySensitivity\', and \'TemperatureSensitvity\''
                            ' created in {}'.format(currentDir+'/FileOutput/'))
# -------- MAKE DIRECTORIES FOR THE REGION OF INTEREST. --------
