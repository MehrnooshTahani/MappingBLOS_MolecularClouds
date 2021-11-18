"""
The zeroth stage of the BLOSMapping method is to define regions of interest.
 - The boundaries and parameters defined here will be used throughout the analysis

TEMPLATE (with defaults)
        elif regionName.lower() == CLOUDNAME:
            # Distance to the region of interest:
            self.distance =  # [pc]
            # Path to the fits file containing to the region of interest:
            self.fitsFilePath = os.path.join(currentDir, '')
            self.fitsDataType = # 'HydrogenColumnDensity' or 'VisualExtinciton'
            # Pixel limits of the region of interest in the fits file:
            self.xmin = 'none'
            self.xmax = 'none'
            self.ymin = 'none'
            self.ymax = 'none'
            # Path to the fiducial extinction and electron abundance for the region of interest:
            self.n0 = ''
            self.T0 = ''
            self.G0 = ''
            Parameters = 'n' + self.n0 + '_T' + self.T0 + '_G' + self.G0
            self.AvFileDir = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/')
            self.AvFilePath = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/Av_T0_n0.out')
            # Galactic Latitude of the region of interest:
            self.cloudLatitude =   # [deg]
            # Boundaries of the region of interest:
            self.raHoursMax =
            self.raMinsMax =
            self.raSecMax =
            self.raHoursMin =
            self.raMinsMin =
            self.raSecMin =
            self.decDegMax =
            self.decDegMin =
"""
import os
from sys import exit
currentDir = os.path.abspath(os.getcwd())


class Region:
    def __init__(self, regionName):

        if regionName.lower() == 'aquila':
            """Parameters corresponding to the region containing the Aquila molecular cloud """
            # Distance to the region of interest:
            self.distance = 400  # CHECK THIS [pc]
            # Path to the fits file containing to the region of interest:
            self.fitsFilePath = os.path.join(currentDir, 'Data/HGBS_aquilaM2_column_density_map.fits')
            self.fitsDataType = 'HydrogenColumnDensity'
            # Pixel limits of the region of interest in the fits file:
            self.xmin = 'none'
            self.xmax = 'none'
            self.ymin = 'none'
            self.ymax = 'none'
            # Path to the fiducial extinction and electron abundance for the region of interest:
            self.n0 = '-1'
            self.T0 = '-1'
            self.G0 = '-1' # can be 1 for most clouds unless clouds with many type o and b stars
            Parameters = 'n' + self.n0 + '_T' + self.T0 + '_G' + self.G0
            self.AvFileDir = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/')
            self.AvFilePath = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/Av_T0_n0.out')
            # Galactic Latitude of the region of interest:
            self.cloudLatitude = 4  # [deg]
            # Boundaries of the region of interest:
            self.raHoursMax = 18
            self.raMinsMax = 40
            self.raSecMax = 0
            self.raHoursMin = 18
            self.raMinsMin = 20
            self.raSecMin = 0
            self.decDegMax = 0
            self.decDegMin = -5

        elif regionName.lower() == 'california':
            """Parameters corresponding to the region containing the California molecular cloud """
            # Distance to the region of interest:
            self.distance = 450.  # CHECK THIS [pc] 470 +/- 2 pc (Zucker et al. 2019)
            # Path to the fits file containing to the region of interest:
            self.fitsFilePath = os.path.join(currentDir, 'Data/2015_02_CalTauPer_toBernsteinCooper.fits')
            self.fitsDataType = 'VisualExtinction'
            # Pixel limits of the region of interest in the file:
            self.xmin = 330
            self.xmax = 760
            self.ymin = 470
            self.ymax = 840
            # Path to the fiducial extinction and electron abundance for the region of interest:
            self.n0 = '4.5e2'
            self.T0 = '10.0'
            self.G0 = '1'
            Parameters = 'n' + self.n0 + '_T' + self.T0 + '_G' + self.G0
            self.AvFileDir = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/')
            self.AvFilePath = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/Av_T0_n0.out')
            # Galactic Latitude of the region of interest:
            self.cloudLatitude = -8.5  # [deg]
            # Boundaries of the region of interest:
            self.raHoursMax = 4
            self.raMinsMax = 36
            self.raSecMax = 0
            self.raHoursMin = 3
            self.raMinsMin = 54
            self.raSecMin = 0
            self.decDegMax = 42
            self.decDegMin = 34.5

        elif regionName.lower() == 'cepheus':
            """Parameters corresponding to the region containing the Cepheus molecular cloud """
            # Distance to the region of interest:
            self.distance = 350  # CHECK THIS [pc]
            # Path to the fits file containing to the region of interest:
            self.fitsFilePath = os.path.join(currentDir, 'Data/HGBS_cep1251_column_density_map.fits')
            self.fitsDataType = 'HydrogenColumnDensity'
            # Pixel limits of the region of interest in the fits file:
            self.xmin = 'none'
            self.xmax = 'none'
            self.ymin = 'none'
            self.ymax = 'none'
            # Path to the fiducial extinction and electron abundance for the region of interest:
            self.n0 = '-1'
            self.T0 = '-1'
            self.G0 = '-1'
            Parameters = 'n' + self.n0 + '_T' + self.T0 + '_G' + self.G0
            self.AvFileDir = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/')
            self.AvFilePath = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/Av_T0_n0.out')
            # Galactic Latitude of the region of interest:
            self.cloudLatitude = 14.5  # [deg]
            # Boundaries of the region of interest:
            self.raHoursMax = 22
            self.raMinsMax = 50
            self.raSecMax = 0
            self.raHoursMin = 22
            self.raMinsMin = 15
            self.raSecMin = 0
            self.decDegMax = 76
            self.decDegMin = 74

        elif regionName.lower() == 'coronaaustralis':
            """Parameters corresponding to the region containing the Corona Australis molecular cloud """
            # Distance to the region of interest:
            self.distance = 131  # CHECK THIS [pc]
            # Path to the fits file containing to the region of interest:
            self.fitsFilePath = os.path.join(currentDir, 'Data/HGBS_craNS_column_density_map.fits')
            self.fitsDataType = 'HydrogenColumnDensity'
            # Pixel limits of the region of interest in the fits file:
            self.xmin = 'none'
            self.xmax = 'none'
            self.ymin = 'none'
            self.ymax = 'none'
            # Path to the fiducial extinction and electron abundance for the region of interest:
            self.n0 = '-1'
            self.T0 = '-1'
            self.G0 = '-1'
            Parameters = 'n' + self.n0 + '_T' + self.T0 + '_G' + self.G0
            self.AvFileDir = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/')
            self.AvFilePath = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/Av_T0_n0.out')
            # Galactic Latitude of the region of interest:
            self.cloudLatitude = -18.5  # [deg]
            # Boundaries of the region of interest:
            self.raHoursMax = 19
            self.raMinsMax = 24
            self.raSecMax = 0
            self.raHoursMin = 18
            self.raMinsMin = 55
            self.raSecMin = 0
            self.decDegMax = -35
            self.decDegMin = -40

        elif regionName.lower() == 'ic5146':
            """Parameters corresponding to the region containing the ic5146 molecular cloud """
            # Distance to the region of interest:
            self.distance = 760  # CHECK THIS [pc]
            # Path to the fits file containing to the region of interest:
            self.fitsFilePath = os.path.join(currentDir, 'Data/HGBS_ic5146_column_density_map.fits')
            self.fitsDataType = 'HydrogenColumnDensity'
            # Pixel limits of the region of interest in the fits file:
            self.xmin = 'none'
            self.xmax = 'none'
            self.ymin = 'none'
            self.ymax = 'none'
            # Path to the fiducial extinction and electron abundance for the region of interest:
            self.n0 = '-1'
            self.T0 = '-1'
            self.G0 = '-1'
            Parameters = 'n' + self.n0 + '_T' + self.T0 + '_G' + self.G0
            self.AvFileDir = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/')
            self.AvFilePath = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/Av_T0_n0.out')
            # Galactic Latitude of the region of interest:
            self.cloudLatitude = -4.5  # [deg]
            # Boundaries of the region of interest:
            self.raHoursMax = 21
            self.raMinsMax = 57
            self.raSecMax = 0
            self.raHoursMin = 21
            self.raMinsMin = 41
            self.raSecMin = 0
            self.decDegMax = 48
            self.decDegMin = 46

        elif regionName.lower() == 'musca':
            """Parameters corresponding to the region containing the Musca molecular cloud """
            # Distance to the region of interest:
            self.distance = 150  # CHECK THIS [pc]
            # Path to the fits file containing to the region of interest:
            self.fitsFilePath = os.path.join(currentDir, 'Data/HGBS_musca_column_density_map.fits')
            self.fitsDataType = 'HydrogenColumnDensity'
            # Pixel limits of the region of interest in the fits file:
            self.xmin = 'none'
            self.xmax = 'none'
            self.ymin = 'none'
            self.ymax = 'none'
            # Path to the fiducial extinction and electron abundance for the region of interest:
            self.n0 = '-1'
            self.T0 = '-1'
            self.G0 = '-1'
            Parameters = 'n' + self.n0 + '_T' + self.T0 + '_G' + self.G0
            self.AvFileDir = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/')
            self.AvFilePath = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/Av_T0_n0.out')
            # Galactic Latitude of the region of interest:
            self.cloudLatitude = -9  # [deg]
            # Boundaries of the region of interest:
            self.raHoursMax = 12
            self.raMinsMax = 45
            self.raSecMax = 0
            self.raHoursMin = 12
            self.raMinsMin = 8
            self.raSecMin = 0
            self.decDegMax = -70
            self.decDegMin = -73

        elif regionName.lower() == 'ophiuchus':
            """Parameters corresponding to the region containing the Ophiuchus molecular cloud """
            # Distance to the region of interest:
            self.distance = 131  # CHECK THIS [pc]
            # Path to the fits file containing to the region of interest:
            self.fitsFilePath = os.path.join(currentDir, 'Data/HGBS_oph_l1688_column_density_map.fits')
            self.fitsDataType = 'HydrogenColumnDensity'
            # Pixel limits of the region of interest in the fits file:
            self.xmin = 'none'
            self.xmax = 'none'
            self.ymin = 'none'
            self.ymax = 'none'
            # Path to the fiducial extinction and electron abundance for the region of interest:
            self.n0 = '-1'
            self.T0 = '-1'
            self.G0 = '-1'
            Parameters = 'n' + self.n0 + '_T' + self.T0 + '_G' + self.G0
            self.AvFileDir = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/')
            self.AvFilePath = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/Av_T0_n0.out')
            # Galactic Latitude of the region of interest:
            self.cloudLatitude = 16.5  # [deg]
            # Boundaries of the region of interest:
            self.raHoursMax = 16
            self.raMinsMax = 36
            self.raSecMax = 0
            self.raHoursMin = 16
            self.raMinsMin = 17
            self.raSecMin = 0
            self.decDegMax = -21
            self.decDegMin = -26

        elif regionName.lower() == 'oriona':
            """Parameters corresponding to the region containing the Orion A molecular cloud """
            # Distance to the region of interest:
            self.distance = 400.  # CHECK THIS [pc] 432 +/- 2 pc (Zucker et al. 2019).
            # Path to the fits file containing to the region of interest:
            self.fitsFilePath = os.path.join(currentDir, 'Data/2016_07_Orion_Plume.fits')
            self.fitsDataType = 'VisualExtinction'
            # Pixel limits of the region of interest in the file:
            self.xmin = 'none'
            self.xmax = 'none'
            self.ymin = 0
            self.ymax = 350
            # Path to the fiducial extinction and electron abundance for the region of interest:
            self.n0 = '1.0e4'
            self.T0 = '25.0'
            self.G0 = '10000'
            Parameters = 'n' + self.n0 + '_T' + self.T0 + '_G' + self.G0
            self.AvFileDir = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/')
            self.AvFilePath = os.path.join(currentDir, 'Data/ChemicalAbundance/'+Parameters+'/Av_T0_n0.out')
            # Galactic Latitude of the region of interest:
            self.cloudLatitude = -19.5  # [deg]
            # Boundaries of the region of interest:
            self.raHoursMax = 5
            self.raMinsMax = 52
            self.raSecMax = 0
            self.raHoursMin = 5
            self.raMinsMin = 22
            self.raSecMin = 0
            self.decDegMax = -4
            self.decDegMin = -11.5

        elif regionName.lower() == 'orionb':
            """Parameters corresponding to the region containing the Orion B molecular cloud """
            # Distance to the region of interest:
            self.distance = 400.  # CHECK THIS [pc]
            # Path to the fits file containing to the region of interest:
            self.fitsFilePath = os.path.join(currentDir, 'Data/2016_07_Orion_Plume.fits')
            self.fitsDataType = 'VisualExtinction'
            # Pixel limits of the region of interest in the file:
            self.xmin = 'none'
            self.xmax = 'none'
            self.ymin = 295
            self.ymax = 700
            # Path to the fiducial extinction and electron abundance for the region of interest:
            self.n0 = '1.0e4'
            self.T0 = '25.0'
            self.G0 = '10000'
            Parameters = 'n' + self.n0 + '_T' + self.T0 + '_G' + self.G0
            self.AvFileDir = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/')
            self.AvFilePath = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/Av_T0_n0.out')
            # Galactic Latitude of the region of interest:
            self.cloudLatitude = -15.5  # [deg]
            # Boundaries of the region of interest:
            self.raHoursMax = 6
            self.raMinsMax = 3
            self.raSecMax = 0
            self.raHoursMin = 5
            self.raMinsMin = 30
            self.raSecMin = 0
            self.decDegMax = 4
            self.decDegMin = -4

        elif regionName.lower() == 'perseus':
            """Parameters corresponding to the region containing the Perseus molecular cloud """
            # Distance to the region of interest:
            self.distance = 250.  # CHECK THIS [pc] 294+/- 10 pc (Zucker et al. 2019).
            # Path to the fits file containing to the region of interest:
            self.fitsFilePath = os.path.join(currentDir, 'Data/2015_02_CalTauPer_toBernsteinCooper.fits')
            self.fitsDataType = 'VisualExtinction'
            # Pixel limits of the region of interest in the file:
            self.xmin = 720
            self.xmax = 1110
            self.ymin = 215
            self.ymax = 560
            # Path to the fiducial extinction and electron abundance for the region of interest:
            self.n0 = '1.0e3'
            self.T0 = '12.0'
            self.G0 = '1'
            Parameters = 'n' + self.n0 + '_T' + self.T0 + '_G' + self.G0
            self.AvFileDir = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/')
            self.AvFilePath = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/Av_T0_n0.out')
            # Galactic Latitude of the region of interest:
            self.cloudLatitude = -19.5  # [deg]
            # Boundaries of the region of interest:
            self.raHoursMax = 3
            self.raMinsMax = 52
            self.raSecMax = 0
            self.raHoursMin = 3
            self.raMinsMin = 20
            self.raSecMin = 0
            self.decDegMax = 34
            self.decDegMin = 28

        elif regionName.lower() == 'pipe':
            """Parameters corresponding to the region containing the Pipe Nebula """
            # Distance to the region of interest:
            self.distance = 200  # CHECK THIS [pc]
            # Path to the fits file containing to the region of interest:
            self.fitsFilePath = os.path.join(currentDir, 'Data/HGBS_pipe_column_density_map.fits')
            self.fitsDataType = 'HydrogenColumnDensity'
            # Pixel limits of the region of interest in the fits file:
            self.xmin = 'none'
            self.xmax = 'none'
            self.ymin = 'none'
            self.ymax = 'none'
            # Path to the fiducial extinction and electron abundance for the region of interest:
            self.n0 = '-1'
            self.T0 = '-1'
            self.G0 = '-1'
            Parameters = 'n' + self.n0 + '_T' + self.T0 + '_G' + self.G0
            self.AvFileDir = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/')
            self.AvFilePath = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/Av_T0_n0.out')
            # Galactic Latitude of the region of interest:
            self.cloudLatitude = 5  # [deg]
            # Boundaries of the region of interest:
            self.raHoursMax = 17
            self.raMinsMax = 36
            self.raSecMax = 0
            self.raHoursMin = 17
            self.raMinsMin = 9
            self.raSecMin = 0
            self.decDegMax = -24
            self.decDegMin = -29

        elif regionName.lower() == 'polaris':
            """Parameters corresponding to the region containing the Polaris Molecular Cloud """
            # Distance to the region of interest:
            self.distance = 500  # CHECK THIS [pc]
            # Path to the fits file containing to the region of interest:
            self.fitsFilePath = os.path.join(currentDir, 'Data/HGBS_polaris_column_density_map.fits')
            self.fitsDataType = 'HydrogenColumnDensity'
            # Pixel limits of the region of interest in the fits file:
            self.xmin = 0
            self.xmax = 6500
            self.ymin = 790
            self.ymax = 6000
            # Path to the fiducial extinction and electron abundance for the region of interest:
            self.n0 = '-1'
            self.T0 = '-1'
            self.G0 = '-1'
            Parameters = 'n' + self.n0 + '_T' + self.T0 + '_G' + self.G0
            self.AvFileDir = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/')
            self.AvFilePath = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/Av_T0_n0.out')
            # Galactic Latitude of the region of interest:
            self.cloudLatitude = 25  # [deg]
            # Boundaries of the region of interest:
            self.raHoursMax = 24
            self.raMinsMax = 0
            self.raSecMax = 0
            self.raHoursMin = 0
            self.raMinsMin = 0
            self.raSecMin = 0
            self.decDegMax = 90
            self.decDegMin = 85

        elif regionName.lower() == 'serpens':
            """Parameters corresponding to the region containing the Serpens Molecular Cloud """
            # Distance to the region of interest:
            self.distance = 400  # CHECK THIS [pc]
            # Path to the fits file containing to the region of interest:
            self.fitsFilePath = os.path.join(currentDir, 'Data/HGBS_serpens_column_density_map.fits')
            self.fitsDataType = 'HydrogenColumnDensity'
            # Pixel limits of the region of interest in the fits file:
            self.xmin = 'none'
            self.xmax = 'none'
            self.ymin = 'none'
            self.ymax = 'none'
            # Path to the fiducial extinction and electron abundance for the region of interest:
            self.n0 = '-1'
            self.T0 = '-1'
            self.G0 = '-1'
            Parameters = 'n' + self.n0 + '_T' + self.T0 + '_G' + self.G0
            self.AvFileDir = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/')
            self.AvFilePath = os.path.join(currentDir, 'Data/ChemicalAbundance/' + Parameters + '/Av_T0_n0.out')
            # Galactic Latitude of the region of interest:
            self.cloudLatitude = 4  # [deg]
            # Boundaries of the region of interest:
            self.raHoursMax = 18
            self.raMinsMax = 45
            self.raSecMax = 0
            self.raHoursMin = 18
            self.raMinsMin = 24
            self.raSecMin = 0
            self.decDegMax = 2
            self.decDegMin = -2

        else:
            print('There is no region of interest with this name.  Check spelling or create a new region of interest'
                  ' then try again.')
            exit()
