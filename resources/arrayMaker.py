


"""
@author: Joab Apaza
@email: japaza@igp.gob.pe

Some usage examples:
antCoordinates = circularArray(radius=50,spacing=6*0.9,nlayers=10, fit=256)
antCoordinates = hexArray(spacing=, nlevels=10, fit=256)
antCoordinates = circleArray(radius=5, nelements=60)


"""
###########################################################################################################################
###########################################################################################################################
from numpy import pi as PI

from numpy import sqrt, floor, sin , cos, arcsin, asarray #, round, abs



###########################################################################################################################
###########################################################################################################################

def hexArray(spacing=10.0, nlevels=7, fit=None):
    # The number of antennas
    # nlevels = 1,2,3...
    SCALE = spacing
    NCEN = nlevels
    ARRAY = []
    XARRAY = []
    YARRAY = []
    XDEL = []
    YDEL = []
    # for I in range(20):
    #     NTOT = 3*(I*I + I) + 1
    #     if NTOT >= nelements:
    #         NCEN = I
    #         break
    #number of rows at the upper part not counting the central
    NROW = NCEN
    DELTAX = 1.0 / (2.0*NCEN)
    DELTAX =  SCALE
    DELTAY = DELTAX * sqrt(3.0) / 2
    NXX = 2*NCEN + 2
    YYY = -DELTAY
    X0 = -DELTAX * NCEN - DELTAX/2
    IND = 0
    for K in range(NROW):
        NXX = NXX -1
        X0 = X0 + DELTAX/2
        YYY = YYY + DELTAY
        XXX = X0 - DELTAX
        for I  in range(NXX):
            XXX = XXX + DELTAX
            IND = IND + 1

            XARRAY.append(XXX)
            YARRAY.append(YYY)
            if K !=0 :
                IND = IND + 1
                XARRAY.append(XXX)
                YARRAY.append(-YYY)

    if fit!=None:
        N = len(XARRAY)

        while (fit < N ):
            ixmax =  XARRAY.index(max(XARRAY))
            XDEL.append(XARRAY[ixmax])
            YDEL.append(YARRAY[ixmax])
            XARRAY.pop(ixmax)
            YARRAY.pop(ixmax)

            iymax =  YARRAY.index(max(YARRAY))
            XDEL.append(XARRAY[iymax])
            YDEL.append(YARRAY[iymax])
            XARRAY.pop(iymax)
            YARRAY.pop(iymax)

            ixmin =  XARRAY.index(min(XARRAY))
            XDEL.append(XARRAY[ixmin])
            YDEL.append(YARRAY[ixmin])
            XARRAY.pop(ixmin)
            YARRAY.pop(ixmin)

            iymin =  YARRAY.index(min(YARRAY))
            XDEL.append(XARRAY[iymin])
            YDEL.append(YARRAY[iymin])
            XARRAY.pop(iymin)
            YARRAY.pop(iymin)

            N = len(XARRAY)

        if fit > N:
            NEELEM =  fit - N
            NDELELEM = len(XDEL)
            if NEELEM > NDELELEM: NEELEM=NDELELEM
            for K in range(NEELEM):
                XARRAY.append(XDEL[ -1*(K+1)])
                YARRAY.append(YDEL[ -1*(K+1)])

    ARRAY = [ [X,Y] for X,Y in zip(XARRAY,YARRAY)]
    ARRAY = asarray(ARRAY).T

    return ARRAY

def circleArray(radius=10.0, nelements=37):
    NARRAY = nelements
    SCALE = radius
    STEPAN = (2*PI) / NARRAY
    ARRAY = []
    for IARRAY in range(NARRAY):
        ANGLE = STEPAN * (IARRAY - 1)
        X = 0.5 * cos(ANGLE)
        Y = 0.5 * sin(ANGLE)
        ARRAY.append([X,Y])
    ARRAY = asarray(ARRAY).T

    return ARRAY*SCALE

def circularArray(radius=10.0, spacing=5.0, nlayers=10, fit=None):
    # The number of antennas at the array should be equal 5,13,39.....= -3(1-2^)
    # NC = 1,2,3...
    SCALE = radius
    ARRAY = []
    XARRAY = [0]
    YARRAY = [0]
    XDEL = []
    YDEL = []
    #number of rows at the upper part not counting the central
    if nlayers > 10:
        NROW = 10
    else:
        NROW = nlayers
    DELTAX = 1.0 / NROW
    DELTAX =  spacing
    RADIUS = DELTAX
    for K in range(NROW):
        if K == 0:
            ARRAY.append([0,0])
            continue
        STEPAN = 2 * arcsin(DELTAX/(2*RADIUS))
        NELEMROW =floor((2*PI)/STEPAN)
        #redefine STEPAN
        STEPAN = 2 * PI /NELEMROW
        for IE in range(int(NELEMROW)):
            ANGLE = STEPAN * IE
            X = RADIUS * cos(ANGLE)
            Y = RADIUS * sin(ANGLE)
            XARRAY.append(X)
            YARRAY.append(Y)
        RADIUS +=DELTAX

    if fit!=None:
        N = len(XARRAY)

        while (fit < N ):
            ixmax =  XARRAY.index(max(XARRAY))
            XDEL.append(XARRAY[ixmax])
            YDEL.append(YARRAY[ixmax])
            XARRAY.pop(ixmax)
            YARRAY.pop(ixmax)

            iymax =  YARRAY.index(max(YARRAY))
            XDEL.append(XARRAY[iymax])
            YDEL.append(YARRAY[iymax])
            XARRAY.pop(iymax)
            YARRAY.pop(iymax)

            ixmin =  XARRAY.index(min(XARRAY))
            XDEL.append(XARRAY[ixmin])
            YDEL.append(YARRAY[ixmin])
            XARRAY.pop(ixmin)
            YARRAY.pop(ixmin)

            iymin =  YARRAY.index(min(YARRAY))
            XDEL.append(XARRAY[iymin])
            YDEL.append(YARRAY[iymin])
            XARRAY.pop(iymin)
            YARRAY.pop(iymin)

            N = len(XARRAY)

        if fit > N:
            NEELEM =  fit - N
            NDELELEM = len(XDEL)
            if NEELEM > NDELELEM: NEELEM=NDELELEM
            for K in range(NEELEM):
                XARRAY.append(XDEL[ -1*(K+1)])
                YARRAY.append(YDEL[ -1*(K+1)])

    ARRAY = [ [X,Y] for X,Y in zip(XARRAY,YARRAY)]
    ARRAY = asarray(ARRAY).T

    return ARRAY

def squaredArray(spacing=5.0, nelements=256):
    SIDE = int(sqrt(nelements))
    DXY = spacing
    ARRAY = []
    for IARRAY in range(SIDE):
        for KARRAY in range(SIDE):
            X = DXY/2 + DXY * KARRAY - (SIDE * DXY)/2
            Y = DXY/2 + DXY * IARRAY - (SIDE * DXY)/2
            ARRAY.append([X,Y])
    ARRAY = asarray(ARRAY).T

    return ARRAY

def rectArray(DX=5.0,DY=5.0, NX=32, NY=64, DX_OFF=None, DY_OFF=None, TGRID=0.0):
    OFF_X = DX/2
    OFF_Y = DY/2
    if DX_OFF != None:
        OFF_X = DX_OFF
    if DY_OFF != None:
        OFF_Y = DY_OFF
    ARRAY = []
    for IARRAY in range(NY):
        for KARRAY in range(NX):
            X = OFF_X + DX * KARRAY - (NX * DX)/2 + (IARRAY%2)*TGRID
            Y = OFF_Y + DY * IARRAY - (NY * DY)/2
            ARRAY.append([X,Y])
    ARRAY = asarray(ARRAY).T

    return ARRAY

from numpy.random import uniform

def generateRandomArray(RADIUS=50,NXY=100):


    # Generate a random angle in radians.

    ARRAY = []
    for N in range(NXY):
        ANGLE = uniform(0, 2 * PI)
        DISTANCE = uniform(0, RADIUS)
        X = DISTANCE * cos(ANGLE)
        Y = DISTANCE * sin(ANGLE)
        ARRAY.append([X,Y])
    ARRAY = asarray(ARRAY).T

    return ARRAY