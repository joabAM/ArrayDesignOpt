
"""

@author: Joab Apaza
@email: japaza@igp.gob.pe

"""
###########################################################################################
###########################################################################################
from numpy import sin,cos,arccos, sqrt
from numpy import NaN
from numpy import pi as PI

def dipoleAntennaFactor(Cy, Cz, ground=True):
    "M.Milla dipole implementation, Cx,Cy grid of cosine values"
    Cy2 = 1 - Cy**2
    Cy2[Cy2<=0] = NaN
    sinth = sqrt(Cy2)
    FDIP = cos(PI/2*Cy)/sinth
    FGND= 2j*sin(PI/2*Cz)
    if ground:
        return FDIP*FGND
    else:
        return FDIP

# def cosineAntennaFactor(Cx, Cy, EF=1.35):
#     "cosine antenna based on Arik.D. Brown book converted to sinespace by Joab Apaza"
#     Ef = EF
#     theta = arccos(sqrt(1-(Cx**2 + Cy**2)))
#     return cos(theta)**(Ef/2)

def cosineAntennaFactor(Cx, Cy, EF=1.35):
    "cosine antenna matlab page I use the same parameter for both"
    #https://www.mathworks.com/help/phased/ug/cosine-antenna-element.html
    Ef = EF
    theta = arccos(sqrt(1-(Cx**2 + Cy**2)))
    return cos(theta)**(Ef)
