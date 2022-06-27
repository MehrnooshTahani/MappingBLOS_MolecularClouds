import math


def getBoxBounds(data, boxXMin, boxXMax, boxYMin, boxYMax):
    xmin = 0
    xmax = data.shape[1]
    ymin = 0
    ymax = data.shape[0]
    if not math.isnan(boxXMin):
        xmin = int(boxXMin)
    if not math.isnan(boxXMax):
        xmax = int(boxXMax)
    if not math.isnan(boxYMin):
        ymin = int(boxYMin)
    if not math.isnan(boxYMax):
        ymax = int(boxYMax)
    return xmin, xmax, ymin, ymax
