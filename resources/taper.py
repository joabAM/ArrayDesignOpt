#######################___ TAPERING___ #######################################
#window functions from:
# https://www.osti.gov/servlets/purl/1365510
# Ref: Catalog of Window Taper Functions for Sidelobe Control
# Armin W. Doerry , Sandia National Laboratories

#writen by Joab Apaza adapted to conformal planar arrays

"""

@author: Joab Apaza
@email: japaza@igp.gob.pe

"""
from numpy import pi as PI

from numpy import sqrt, floor, sin , cos, arcsin, asarray , radians


def taperTri(XPos,YPos,DIAMETER, X0=None, Y0=None):
    "Shift X and Y coordinates to eliminate negative values"
    X = XPos + DIAMETER/2
    Y = YPos + DIAMETER/2
    K = DIAMETER
    if X0==None and Y0==None:
        #adjusting the center
        X0 = K/2
        Y0 = K/2
    k = sqrt( (X-X0)**2 + (Y-Y0)**2)
    FT =  abs(k)/K
    fm , fn = FT.max(), FT.min()
    return  (FT-fn)/(fm-fn)

def taperHam(XPos,YPos,DIAMETER, a=0.54):
    "Shift X and Y coordinates to eliminate negative values"
    X = XPos + DIAMETER/2
    Y = YPos + DIAMETER/2
    K = DIAMETER
    #adjusting the new center
    X0 = K/2
    Y0 = K/2
    k = sqrt( (X-X0)**2 + (Y-Y0)**2)
    FH = 1 + ((1-a)/a)**cos(2*PI*k/K)
    fm , fn = FH.max(), FH.min()
    return  (FH-fn)/(fm-fn)

def taperCos2(XPos,YPos,DIAMETER):
    #The raised-cosine or von Hann window
    "Shift X and Y coordinates to eliminate negative values"
    X = XPos + DIAMETER/2
    Y = YPos + DIAMETER/2
    K = DIAMETER
    #adjusting the new center
    X0 = K/2
    Y0 = K/2
    k = sqrt( (X-X0)**2 + (Y-Y0)**2)
    FC2 = 0.5+0.5*cos(2*PI*k/K)
    return FC2
    #fm , fn = FC2.max(), FC2.min()
    #return  (FC2-fn)/(fm-fn)

def taperBlackM(XPos,YPos,DIAMETER):
    "Shift X and Y coordinates to eliminate negative values"
    X = XPos + DIAMETER/2
    Y = YPos + DIAMETER/2
    K = DIAMETER
    #adjusting the new center
    X0 = K/2
    Y0 = K/2
    k = sqrt( (X-X0)**2 + (Y-Y0)**2)
    FBM = 1+ (0.5/0.42)*cos(2*PI*k/K) + (0.08/0.42)*cos(2*PI*k/K)
    fm , fn = FBM.max(), FBM.min()
    return  (FBM-fn)/(fm-fn)

def taperParabolic(XPos,YPos,DIAMETER):
    "Shift X and Y coordinates to eliminate negative values"
    X = XPos + DIAMETER/2
    Y = YPos + DIAMETER/2
    K = DIAMETER
    #adjusting the new center
    X0 = K/2
    Y0 = K/2
    k = sqrt( (X-X0)**2 + (Y-Y0)**2)
    FP = 1 - (abs(k)/K)**2
    fm , fn = FP.max(), FP.min()
    return  (FP-fn)/(fm-fn)

def taperWelch(XPos,YPos,DIAMETER):
    #or Riesz, Bochner, Parzen
    "Shift X and Y coordinates to eliminate negative values"
    X = XPos + DIAMETER/2
    Y = YPos + DIAMETER/2
    K = DIAMETER
    #adjusting the new center
    X0 = K/2
    Y0 = K/2
    k = sqrt( (X-X0)**2 + (Y-Y0)**2)
    FW = 1.5*(1 - 4*k**2)/K
    fm , fn = FW.max(), FW.min()
    return  (FW-fn)/(fm-fn)

################################################################################
################################################################################
################################################################################

def taperTri_UV(XPos, YPos, DIAMETER, AZI=None, ELEV=None):
    X,Y, X0, Y0 = normTaper(XPos, YPos,DIAMETER, AZI=AZI, ELEV=ELEV)
    k = sqrt( (X-X0)**2 + (Y-Y0)**2)
    FT = 1 - abs(k)
    fm , fn = FT.max(), FT.min()
    return  (FT-fn)/(fm-fn)

def taperHam_UV(XPos,YPos,DIAMETER, AZI=None, ELEV=None):
    X,Y, X0, Y0 = normTaper(XPos, YPos,DIAMETER, AZI=AZI, ELEV=ELEV)
    k = sqrt( (X-X0)**2 + (Y-Y0)**2)
    a = 0.54
    FH = 1 + ((1-a)/a)**cos(2*PI*k)
    fm , fn = FH.max(), FH.min()
    return  1 - (FH-fn)/(fm-fn)

def taperCos2_UV(XPos,YPos,DIAMETER, AZI=None, ELEV=None):
    #The raised-cosine or von Hann window
    X,Y, X0, Y0 = normTaper(XPos, YPos,DIAMETER, AZI=AZI, ELEV=ELEV)
    k = sqrt( (X-X0)**2 + (Y-Y0)**2)
    a = 0.5
    FC2 = 1 + ((1-a)/a)*cos(2*PI*k)
    #return  FC2
    fm , fn = FC2.max(), FC2.min()
    return  (FC2-fn)/(fm-fn)

def taperBlackM_UV(XPos,YPos,DIAMETER, AZI=None, ELEV=None):
    X,Y, X0, Y0 = normTaper(XPos, YPos,DIAMETER, AZI=AZI, ELEV=ELEV)

    k = sqrt( (X-X0)**2 + (Y-Y0)**2)
    FBM = 1+ (0.5/0.42)*cos(2*PI*k) + (0.08/0.42)*cos(2*PI*k)
    #return  FBM
    fm , fn = FBM.max(), FBM.min()
    return  (FBM-fn)/(fm-fn)

def taperParabolic_UV(XPos,YPos,DIAMETER, AZI=None, ELEV=None):
    X,Y, X0, Y0 = normTaper(XPos, YPos,DIAMETER, AZI=AZI, ELEV=ELEV)
    k = sqrt( (X-X0)**2 + (Y-Y0)**2)
    FP = 1 - abs(k)**2
    #return  FP
    fm , fn = FP.max(), FP.min()
    return  (FP-fn)/(fm-fn)

def taperWelch_UV(XPos,YPos,DIAMETER, AZI=None, ELEV=None):
    #or Riesz, Bochner, Parzen
    X,Y, X0, Y0 = normTaper(XPos, YPos,DIAMETER, AZI=AZI, ELEV=ELEV)
    k = sqrt( (X-X0)**2 + (Y-Y0)**2)
    FW = 1.5*(1 - 4*k**2)
    return  FW
    #fm , fn = FW.max(), FW.min()
    #return  (FW-fn)/(fm-fn)


################################################################################
def normTaper(XPos, YPos,DIAMETER, AZI=None, ELEV=None):
    "Shift X and Y coordinates to eliminate negative values"
    X = XPos + DIAMETER/2
    Y = YPos + DIAMETER/2
    #normalize
    X /=DIAMETER
    Y /=DIAMETER

    if AZI==None and ELEV==None:
        cx0 = None
        cy0 = None
    else:
        cx0 = cos(radians(ELEV))*sin(radians(AZI))
        cy0 = cos(radians(ELEV))*cos(radians(AZI))

    if cx0==None and cy0==None:
        #adjusting the center
        X0 = 0.5
        Y0 = 0.5
    else:
        X0 = (cx0 +1)/2
        Y0 = (cy0 +1)/2

    return X, Y, X0, Y0