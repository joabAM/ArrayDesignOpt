
"""

@author: Joab Apaza
@email: japaza@igp.gob.pe

"""


from numpy import nanargmax,  nanmin,log10, sqrt, ogrid


def circularMask(h, w, center=None, radius=None, radius2=None):
    flag_ellipse = False
    if center is None: # use the middle of the array
        center = (int(w/2), int(h/2))
    if radius is None: # use the smallest distance between the center and array edges
        radius = min(center[0], center[1], w-center[0], h-center[1])
    if radius2 is not None:
        flag_ellipse = True
    yo, xo = ogrid[:h, :w]

    if flag_ellipse:
        dist_from_center = ((xo - center[0])**2)/radius**2 + ((yo-center[1])**2)/radius2**2
        mask = dist_from_center <= 1

    else:
        dist_from_center = sqrt((xo - center[0])**2 + (yo-center[1])**2)
        mask = dist_from_center <= radius
    return mask

from scipy.signal import argrelmin,argrelmax

#global findBeamRadius
def findBeamRadius(data, nx, ny, nOr=1, z=1):
    #NUMPY ARR[Y,X], ROW,COLUMS
    #normalized data, 0dB maximum
    logData = 10*log10(data)
    #plt.imshow(logData)
    MAXVAL = nanargmax(logData)
    imax_x, imax_y = int(MAXVAL%ny), (MAXVAL//nx)
    #[imax_y,imax_x] = argwhere(logData == 0.0)[0]
    lefty=False
    leftx=False
    if len(logData[imax_y,0:imax_x]) > len(logData[imax_y,imax_x:]):
        xdata = logData[imax_y,0:imax_x]
        leftx = True
    else:
        xdata = logData[imax_y,imax_x:]

    if len(logData[0:imax_y,imax_x]) > len(logData[imax_y:,imax_x]):
        ydata = logData[0:imax_y,imax_x]
        lefty = True
    else:
        ydata = logData[imax_y:,imax_x]

    minx = nanmin(xdata)
    miny = nanmin(ydata)
    xdata -= minx
    ydata -= miny
    xmin = argrelmin(xdata, order=nOr, mode='wrap')[0]
    ymin = argrelmin(ydata, order=nOr, mode='wrap')[0]

    ixmin = 0
    iymin = 0
    if leftx:
        ixmin=len(xdata) - xmin[-1]
    else:
        ixmin=xmin[0]
    if lefty:
        iymin=len(ydata) - ymin[-1]
    else:
        iymin=ymin[0]

    if len(xmin)>1:
        #r1 = (nx-xmin[int(len(xmin)/2)])*(2*z)/nx -z
        r1 = (ixmin)*2*z/nx
    else:
        r1=1
    if len(ymin)>1:
        #r2 = (ny-ymin[int(len(ymin)/2)])*(2*z)/ny -z
        r2 = (iymin)*2*z/ny
    else:
        r2=1

    #print("radius: ", r1,r2,imax_x, imax_y)
    return r1, r2, imax_x, imax_y