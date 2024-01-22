

""" Short description of this Python module.
Longer description of this module.


This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.
This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>.

This library was developed and tested with matplotlib v3.6.0

"""

__author__ = "Joab Apaza"
__contact__ = "japaza@igp.gob.pe"
__copyright__ = "Geophysical Institute of Peru"
#__credits__ = ["Joab Apaza"]
__deprecated__ = False
__email__ =  "japaza@igp.gob.pe"
#__license__ = "GPLv3"
__maintainer__ = "developer"
__status__ = "dev"
__version__ = "0.0.1"



import  shapely
import traceback
from matplotlib import pyplot as plt
import warnings
import matplotlib.cbook

from IPython.display import display, clear_output

from time import sleep
from os import getcwd
from os import path
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
from matplotlib.patches import Ellipse

from numpy import pi as PI
from numpy import NaN
from numpy import isfinite, nanargmax, nanmax, nansum, nanmin, fmax

from numpy import log10, power, sqrt, exp, floor #, round, abs

from numpy import sin , cos, tan,arccos, arcsin, radians, degrees

from numpy import linspace, arange, zeros, ones, asarray, meshgrid, ogrid, loadtxt,savetxt, where, flip, sign, roll

from numpy import mean, maximum, argwhere, arctan2, uint8

from numpy import column_stack, unique, argsort,isnan



################################################################################
################################################################################

from resources.arrayMaker import *
from utils.elementPatterns import *
from utils.fbeam import *

################################################################################
################################################################################
class AntennaArray(object):

    antCoordinates = []
    Xdata = []
    Ydata = []

    def  __init__(self, type="squared", N=256,  file=None, **kwargs,):
        if file != None:
            try:
                self.antCoordinates = loadtxt(file).T
                self.Xdata = self.antCoordinates[1]
                self.Ydata = self.antCoordinates[2]
            except Exception:
                traceback.print_exc()
        else:
            if type=="hex":
                self.antCoordinates = hexArray(**kwargs)
            elif type=="circle":
                self.antCoordinates = circleArray(**kwargs)
            elif type=="circular":
                self.antCoordinates = circularArray(**kwargs)
            elif type=="squared":
                self.antCoordinates = squaredArray(**kwargs)
            elif type=="rect":
                self.antCoordinates = rectArray(**kwargs)
            elif type=="rand":
                self.antCoordinates = generateRandomArray(**kwargs)
            else:
                print("Invalid array type selected")
                return 0
            self.Xdata = self.antCoordinates[0]
            self.Ydata = self.antCoordinates[1]

    def info(self):
        pass

    def setElement(self, Cy, Cz, element="cosine", EF=1.1, ground=True):
        pass

         
    def saveLocs2File(self,Xdata=[], Ydata=[], filename=""):
        if len(Xdata)==0 and len(Ydata)==0:
            XPOS = self.Xdata
            YPOS = self.Ydata
        else:
            XPOS = Xdata
            YPOS = Ydata
        if len(XPOS) != len(YPOS):
            print("mismatch data length")
            return
        N = len(XPOS)
        IDX = linspace(1,N,num=N,endpoint=True)
        LOCS = zeros( (3, N) )
        LOCS[0] = IDX
        LOCS[1] = XPOS
        LOCS[2] = YPOS
        if filename=="":
            aux = getcwd()
            filename = path.join(aux,"arrayDesignOpt/txt/antLocation.txt")
        savetxt(filename, LOCS.T, fmt='%3.5e')

    def rotateCoordinates(self, X, Y, TH=-6.0, SAVE=None, header=""):
        N = len(self.Xdata)
        X = self.Xdata
        Y = self.Ydata
        LOCS = zeros( (3, self.Xdata) )
        IDX = linspace(1,N,num=N,endpoint=True)
        LOCS[0] = IDX
        LOCS[1] = X*cos(radians(TH)) - Y*sin(radians(TH))  # X rot
        LOCS[2] = Y*cos(radians(TH)) + X*sin(radians(TH))  # Y rot
        if SAVE!=None:
            savetxt(SAVE, LOCS.T, fmt='%3.5e')
        return LOCS[1], LOCS[2]

    def showArrayCoordinates(self, xmin=None, xmax=None, ymin=None, ymax=None, radius=0,  width=0, height=0):
        fig = plt.figure("Antenna Element Locations",figsize=(6, 6))
        #plt.axes().set_aspect('equal')
        plt.scatter(self.Xdata, self.Ydata, marker='x')
        plt.minorticks_on()
        plt.grid(b=True, which='major', color='0.9', linestyle='--')
        plt.grid(b=True, which='minor', color='0.9', linestyle='--')

        plt.ylabel("meters in Y (v) direction")
        plt.xlabel("meters in X (u) direction")
        
        ax = plt.gca()

        if radius!=0:
            circle=plt.Circle((0,0),radius,  color='k', fill=False)
            ax.add_patch(circle)
        
        if width!=0 and height!=0:
            ellipse = Ellipse(xy=(0, 0), width=100, height=120, edgecolor='k', fc='None', lw=1)
            ax.add_patch(ellipse)

        #rect = plt.Rectangle((-RAD,-RAD),2*RAD, 2*RAD,  color='k', fill=False)
        #ax.add_patch(rect)
            
        if xmin!=None and xmax!=None:
            plt.xlim([xmin,xmax])
        if ymin!=None and ymax!=None:
            plt.xlim([ymin,ymax])

        plt.title('Antenna Array %d elements'%len(self.Xdata))
                   
        plt.show()

    def getPattern(self, XPOS=[], YPOS=[], WGT=[],FC=50, ELEV=90, AZI=0, MIN_DB=-50, MAX_DB=0, XMAX=1, XMIN=-1, YMAX=1, YMIN=-1,ANTENNA="COSINE", PLOT_RADMIN=False,PLOT_LEVELS=False,
               POINTSXY=1025, EF=1.1, GND=True, MOD_Z=1.0):
        if len(XPOS)==0 and len(YPOS)==0:
            XPOS = self.Xdata
            YPOS = self.Ydata

        c = 2.997962*10**8
        f = FC              #frequency in MHz
        LAMBDA =  c/(f*10**6)  #wavelenght (m)
        XF ,YF = XPOS/LAMBDA, YPOS/LAMBDA
        # Elevation = 32
        # Azimuth = 273.27
        # ELEV = 62.8
        # AZI = 346.33
        CX0 = cos(radians(ELEV))*sin(radians(AZI))
        CY0 = cos(radians(ELEV))*cos(radians(AZI))

        NARRAY = len(XPOS)
        NX = POINTSXY
        NY = POINTSXY
        Z = (MOD_Z)
        #print(Z)
        dcx = (2.0*Z)/NX
        dcy = (2.0*Z)/NY
        # cosx = dcx*(arange(NX)-floor(NX/2))
        # cosy = dcy*(arange(NY)-floor(NY/2))
        cosx = linspace(-Z, Z, num=NX, endpoint=True)
        cosy = linspace(-Z, Z, num=NY, endpoint=True)
        if len(WGT)<1:
            WGT = ones(NARRAY)
        #M.Milla conversion
        CX, CY = meshgrid(cosx, cosy)
        Cz2 = 1-(CX**2+CY**2)
        Cz2[Cz2<=0] = NaN
        Cz = sqrt(Cz2)
        domg = 1./Cz
        domg[~isfinite(domg)] = 0
        if Z!=1:
            ANTENNA=""
            XMAX = Z
            XMIN = -Z
            YMAX = Z
            YMIN = -Z
        FANT = 0
        if ANTENNA=="DIPOLE":
            FANT = dipoleAntennaFactor(CY,Cz,ground=GND)
        elif ANTENNA=="COSINE":
            FANT = cosineAntennaFactor(CX,CY,EF=EF)
        else:
            print("Missing antenna element")
            FANT = 1

        FARR = 0
        for IARRAY in range(NARRAY):

            FARR += WGT[IARRAY]*exp( 1j*2*PI*(XF[IARRAY]*(CX-CX0)  + YF[IARRAY]*(CY-CY0)))

        #antenna pattern
        FARR = FARR*FANT
        PATT = power(abs(FARR),2)
        maskUV = circularMask(NX, NY, radius=int(NX*0.5))
        PATT[~maskUV] = NaN
        PATT_NORM = PATT/nanmax(PATT)

        #Integrated gain (M.Milla)
        G1W = abs(FARR)**2
        G1W = (4*PI/((dcx*dcy)*nansum(G1W*domg)))*G1W
        G2W = G1W #* G1W #¿Why?check  getBeamINT



        RADMIN1, RADMIN2, xcenter, ycenter  = findBeamRadius(PATT_NORM, NX, NY, nOr=1, z=Z)
        WZENX, WZENY = 0,0
        mask = circularMask(NX, NY, radius=int(NX*RADMIN1/2), radius2=int(NX*RADMIN2/2))
        maskedG2W = G2W.copy()
        maskedG2W[mask] = NaN
        try:
            IMAX = (nanargmax(maskedG2W))
            _x, _y = (IMAX%NX), (IMAX/NY)
            WZENX, WZENY = (_x/NX)*(2*Z) -Z , (_y/NX)*(2*Z) -Z
        except:
            PLOT_RADMIN = False

        #print(IMAX, _x, _y, WZENX, WZENY)

        fig = plt.figure("Power Pattern", figsize=(6, 6))
        plt.axes().set_aspect('equal')
        im = plt.pcolormesh(cosy,cosx,10*log10(PATT_NORM),shading='auto',cmap='jet')

        plt.colorbar(im, fraction=0.046, pad=0.04)
        plt.clim([MIN_DB,MAX_DB])


        print("RADMIN1, RADMIN2 " , RADMIN1, RADMIN2)
        if PLOT_RADMIN:
            #circle1=plt.Circle(( (ycenter/NX)*(2*Z) -Z , (xcenter/NY)*(2*Z) -Z),abs(RADMIN1),  color='k', fill=False)
            # circle2=plt.Circle((CX0 ,CY0),abs(RADMIN1),  color='k', fill=False)
            ellipse=Ellipse((CX0,CY0),2*RADMIN1, 2*RADMIN2, color='k', fill=False)

            ax = plt.gca()
            #ax.add_patch(circle1)
            # ax.add_patch(circle2)
            ax.add_patch(ellipse)
            plt.scatter(WZENX, WZENY, marker='x', color='k')

        if PLOT_LEVELS:
            plt.contour(cosx,cosy,10*log10(PATT_NORM),levels=[ -30, -20, -10, -3],colors=['k', 'w','g','w'],linestyles=['dashed'])
        plt.xlabel('cosx')
        plt.ylabel('cosy')
        plt.title('Radiation Pattern {} elements @ {} points'.format(NARRAY,POINTSXY))
        plt.ylim((YMIN, YMAX))
        plt.xlim((XMIN, XMAX))
        plt.show()
        return PATT_NORM, domg, cosy, RADMIN1, G2W

    def getFirstNull(self,data, reverse=False):
        F = data
        N = len(F)
        if reverse:
            F = flip(data)
        slopes = [ (F[k+1]-F[k]) for k in arange(0,N-1)]
        d_sign = sign(slopes)
        d_sign[0] = d_sign[1] #first element slope depends on the second
        signchange = ((roll(d_sign, 1) - d_sign) != 0).astype(int)
        signchange[0]=0  #first element sign doesnt have change
        firstNulls = where(signchange==1)[0]
        iFN = where(d_sign[firstNulls]>0)[0][0]
        if reverse:
            return (N-1) - firstNulls[iFN]
        return firstNulls[iFN]

    def getSLL_HPBW_FNBW(self,IN_PATTERN, RADMIN, HYF=True,NAME="", NORM=True, DEGREES=False, XMIN=None, XMAX=None, YMIN=-60, YMAX=0.0, SHOW=True, DECIMALS=2):
        
        print("\n------------------------------------------------------------------------------------\n")

        PATTERN = IN_PATTERN.copy()
        NX, NY = PATTERN.shape
        OFFSET = 0
        LVL = -0.0 #level offset, to find....
        if not NORM: #get directivity
            MAX_REF_LEVEL = nanmax(PATTERN)
            OFFSET = 10*log10(MAX_REF_LEVEL)
            if SHOW:
                print("Maximum Integrated Gain (Directivity):",round(OFFSET,2),"dBi")
                print("")
            PATTERN /= nanmax(PATTERN)#normalize the pattern
            if YMAX == 0:
                YMAX = 40
        RADMIN1, RADMIN2, IMAX_X, IMAX_Y = findBeamRadius(PATTERN, NX, NY)

        #print(IMAX_X, IMAX_Y)

        figpC, axs = plt.subplots(2)
        figpC.suptitle(NAME+" Pattern")
        figpC.canvas.manager.set_window_title('Pattern Cuts')
        figpC.set_size_inches(10,6)
        axs[0].grid()
        axs[1].grid()
        axs[0].minorticks_on()
        axs[1].minorticks_on()
        PATTERNXDB = 10*log10(PATTERN[IMAX_Y,:]) + OFFSET
        PATTERNYDB = 10*log10(PATTERN[:,IMAX_X]) + OFFSET

        COSX = linspace(-1, 1, num=NX)
        COSY = linspace(-1, 1, num=NY)


        ANGLEX = degrees(arcsin(COSX))
        ANGLEY = degrees(arcsin(COSX))
        if SHOW:
            if DEGREES:
                axs[0].plot(ANGLEX,PATTERNXDB) # v = cte
                axs[1].plot(ANGLEY,PATTERNYDB) # u = cte
                axs[0].set_xlabel('° deg x')
                axs[1].set_xlabel('° deg y')
                if XMIN==None and XMAX==None:
                    axs[0].set_xlim(-90,90)
                    axs[1].set_xlim(-90,90)
                else:
                    axs[0].set_xlim(XMIN,XMAX)
                    axs[1].set_xlim(XMIN,XMAX)
            else:
                axs[0].plot(COSX,PATTERNXDB) # v = cte
                axs[1].plot(COSY,PATTERNYDB) # u = cte
                axs[0].set_xlabel('cosx')
                axs[1].set_xlabel('cosy')
                if XMIN==None and XMAX==None:
                    axs[0].set_xlim(-1.0,1.0)
                    axs[1].set_xlim(-1.0,1.0)
                else:
                    axs[0].set_xlim(XMIN,XMAX)
                    axs[1].set_xlim(XMIN,XMAX)
        axs[0].set_ylim(YMIN, YMAX)
        axs[1].set_ylim(YMIN, YMAX)



        SLLX1 = PATTERN[IMAX_Y, 0:int(IMAX_X - RADMIN2*NX/2),   ]
        SLLX2 = PATTERN[IMAX_Y, int(IMAX_X + RADMIN2*NX/2): ]

        SLLY1 = PATTERN[0:int(IMAX_Y - RADMIN1*NY/2), IMAX_X]
        SLLY2 = PATTERN[int(IMAX_Y + RADMIN1*NY/2):, IMAX_X]
        ##########################################################
        try:
            xSLL1 = nanargmax(SLLX1)
        except:
            xSLL1 = None
        try:
            xSLL2 = nanargmax(SLLX2)
        except:
            xSLL2 = None
        ###############
        if xSLL1!=None:
            SLLXDB1 = round(10*log10(SLLX1[xSLL1])+ LVL,2)
        else:
            SLLXDB1 = -10000
        if xSLL2!=None:
            SLLXDB2 = round(10*log10(SLLX2[xSLL2])+ LVL,2)
        else:
            SLLXDB2 = -10000


        ##########################################################
        try:
            ySLL1 = nanargmax(SLLY1)
        except:
            ySLL1 = None
        try:
            ySLL2 = nanargmax(SLLY2)
        except:
            ySLL2 = None
        ###############
        if ySLL1!=None:
            SLLYDB1 = round(10*log10(SLLY1[ySLL1])+ LVL,2)
        else:
            SLLYDB1 = -10000
        if ySLL2!=None:
            SLLYDB2 = round(10*log10(SLLY2[ySLL2])+ LVL,2)
        else:
            SLLYDB2 = -10000
        ##########################################################


        SLLXDB = max(SLLXDB1, SLLXDB2)
        SLLYDB = max(SLLYDB1, SLLYDB2)
        # print(repr(10*log10(SLLY1 + LVL)))
        # print(repr(10*log10(SLLY2 + LVL)))
        if SHOW:
            print("Side Lobe Levels")
            print("SLL X = " , SLLXDB, "dB")
            print("SLL Y = " , SLLYDB, "dB")

        ##########################################################
        ##########################################################
        #SLL
        ##########################################################
        ##########################################################

        #masked the main beam
        mask = circularMask(NX, NY,center=(IMAX_X,IMAX_Y), radius=int( NX * RADMIN1/2),  radius2=int( NY * RADMIN2/2))
        maskedPATTERN = PATTERN.copy()
        maskedPATTERN[mask] = NaN
        # plt.figure()
        # cosxy = linspace(-1, 1, num=NX, endpoint=True)
        # plt.pcolormesh(cosxy, cosxy,maskedPATTERN)
        # plt.show()
        IMAX = (nanargmax(maskedPATTERN))
        _x, _y = (IMAX%NX), (IMAX//NY)
        #SLL location in the uv plane
        iSSLx, iSLLy = (_x/NX)*2 -1 , (_y/NX)*2 -1  #to -1 to 1

        SSLOUT = nanmax(maskedPATTERN/nanmax(PATTERN))
        #SSLOUT = PATTERN[ _y, _x]
        G_SLL = round(10*log10(SSLOUT)+LVL,DECIMALS)

        if not HYF:
            return G_SLL

        if SHOW:
            print("SLL: ", G_SLL ,"dB")# #, "\t \t u, v:", iSSLx, iSLLy)
            print("")
            #print(iSSLx, iSLLy)
        ##########################################################
        ##########################################################
        #HPBW
        ##########################################################
        ##########################################################
        PATTERNXDB -= OFFSET
        PATTERNYDB -= OFFSET
        locs3dBx1 = where(PATTERNXDB >= -3)[0][0]
        locs3dBx2 = NX - where(flip(PATTERNXDB) >= -3)[0][0]

        locs3dBy1 = where(PATTERNYDB >= -3)[0][0]
        locs3dBy2 = NY - where(flip(PATTERNYDB) >= -3)[0][0]

        HPBWX = round(abs(ANGLEX[locs3dBx1] - ANGLEX[locs3dBx2]),2)
        HPBWY = round(abs(ANGLEY[locs3dBy1] - ANGLEY[locs3dBy2]),2)
        if SHOW:
            print("Half Power Beam Width")
            print("HPBW X: ",HPBWX)
            print("HPBW Y: ",HPBWY)
            print("")

        ##########################################################
        ##########################################################
        #FNBW
        ##########################################################
        ##########################################################
        #
        # # R1, R2, xcenter, ycenter  = findBeamRadius(PATTERN, NX, NY, nOr=1)
        # # xmb = ((xcenter/NX)*2 - 1)
        # # ymb = ((ycenter/NY)*2 - 1)
        # # # print(R1, R2, xmb, ymb)
        # # FNBWX_l = degrees(arcsin(xmb - R1))
        # # FNBWX_h = degrees(arcsin(xmb + R1))
        # # FNBWY_l = degrees(arcsin(ymb - R2))
        # # FNBWY_h = degrees(arcsin(ymb + R2))
        # # # print(FNBWX_l, FNBWX_h, FNBWY_l, FNBWY_h)
        # # FNBWX = round(abs(FNBWX_h-FNBWX_l),2)
        # # FNBWY = round(abs(FNBWY_h-FNBWY_l),2)
        # # print("First Null Beam Width")
        # # print("FNBW X: ",FNBWX)
        # # print("FNBW Y: ",FNBWY)

        xDB1 = PATTERNXDB[0:IMAX_X]
        xDB2 = PATTERNXDB[IMAX_X:]
        yDB1 = PATTERNYDB[0:IMAX_Y]
        yDB2 = PATTERNYDB[IMAX_Y:]
        try:
            i_FNBWX_l = self.getFirstNull(xDB1, reverse=True)
        except:
            i_FNBWX_l = 0
        try:
            i_FNBWX_h = self.getFirstNull(xDB2, reverse=False)
        except:
            i_FNBWX_h = NX - 1
            IMAX_X = 0
        FNBWX_l = ANGLEX[i_FNBWX_l]
        FNBWX_h = ANGLEX[i_FNBWX_h + IMAX_X]
        try:
            i_FNBWY_l = self.getFirstNull(yDB1, reverse=True)
        except:
            i_FNBWY_l = 0
        try:
            i_FNBWY_h = self.getFirstNull(yDB2, reverse=False)
        except:
            i_FNBWY_h = NY
        FNBWY_l = ANGLEY[NY-1 - i_FNBWY_l]
        FNBWY_h = ANGLEY[NY-1 - (i_FNBWY_h + IMAX_Y)]
        FNBWX = round(abs(FNBWX_h-FNBWX_l),2)
        FNBWY = round(abs(FNBWY_h-FNBWY_l),2)
        # print(FNBWX_l, FNBWX_h, FNBWY_l, FNBWY_h)
        if SHOW:
            print("First Null Beam Width")
            print("FNBW X: ",FNBWX)
            print("FNBW Y: ",FNBWY)


        if SHOW:
            print("\nBeam Simmetry")
            ANG_X_MAX = ANGLEX[IMAX_X]
            ANG_Y_MAX = ANGLEY[NY-1 - IMAX_Y]
            print("X: ",round(abs(ANG_X_MAX-FNBWX_l),2),"|", round(abs(ANG_X_MAX-FNBWX_h),2) )
            print("Y: ",round(abs(ANG_Y_MAX-FNBWY_l),2),"|", round(abs(ANG_Y_MAX-FNBWY_h),2) )


        if SHOW:
            print("\nAxis Ratio")
            print("HPBW X/Y: ",round(HPBWX/HPBWY,2))
            print("FNBW X/Y: ",round(FNBWX/FNBWY, 2))
            print("")
        plt.show()
        return G_SLL, SLLXDB, SLLYDB, FNBWX, FNBWY, _x, _y


    def getBeamINT (self, GAIN_INT, DOMG, COSXY, STEP=0.01, MIN=0.1, MAX=10.0, XMIN=-7,XMAX=7, TWOWAYS=False):
        GAIN_INT_NORM = GAIN_INT/nanmax(GAIN_INT)
        if TWOWAYS:
            GAIN_INT2 = GAIN_INT* GAIN_INT
            GAIN_INT2_NORM = GAIN_INT2/nanmax(GAIN_INT2)
            GAIN_INT_NORM = GAIN_INT2_NORM

        G2W_INT1 =nansum(GAIN_INT_NORM*DOMG,axis=0) #theta x
        G2W_INT2 = nansum(GAIN_INT_NORM*DOMG,axis=1) #theta y
        PEAK1 = G2W_INT1.max()
        PEAK2 = G2W_INT2.max()
        iPEAK1 = where(G2W_INT1==G2W_INT1.max())[0][0]
        iPEAK2 = where(G2W_INT2==G2W_INT2.max())[0][0]
        #print(iPEAK1,iPEAK2)
        err_min=1e6
        beam1 = 0
        beam2 = 0
        for beam in arange(MIN, MAX, STEP):
            #f = PEAK1*exp(-(90* (COSXY- COSXY[iPEAK1]))**2/(2*beam**2))
            f = PEAK1*exp(-(90* COSXY)**2/(2*beam**2))
            err = abs(G2W_INT1 - f).sum()
            if err < err_min:
                err_min = err
                beam1 = beam
        err_min=1e6
        for beam in arange(MIN, MAX, STEP):
            #f = PEAK2*exp(-(90*(COSXY- COSXY[iPEAK2]))**2/(2*beam**2))
            f = PEAK2*exp(-(90* COSXY)**2/(2*beam**2))
            err = abs(G2W_INT2 - f).sum()
            if err < err_min:
                err_min = err
                beam2 = beam

        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(6, 3))
        #M.Milla model
        #Obs: 90*Cosxy, probably needs a correction
        # axis 1
        axs[0].plot(90*COSXY,G2W_INT1 ,label='int pattern 1')
        axs[0].plot(90*COSXY, PEAK1*exp(-(90*COSXY)**2/(2*beam1**2)) ,label='fitted pattern 1')
        axs[0].set_xlabel('$\\theta_x$ [deg]')
        axs[0].set_title('Integrated pattern 1')
        axs[0].legend()
        axs[0].set_xlim([XMIN+90*COSXY[iPEAK1],XMAX+90*COSXY[iPEAK1]])

        # second axis
        axs[1].plot(90*COSXY,G2W_INT2 ,label='int pattern 2')
        axs[1].plot(90*COSXY, PEAK2*exp(-(90*COSXY)**2/(2*beam2**2)), label='fitted pattern 2')
        axs[1].set_xlabel('$\\theta_y$ [deg]')
        axs[1].set_title('Integrated pattern 2')
        axs[1].legend()
        axs[1].set_xlim([XMIN+90*COSXY[iPEAK2],XMAX+90*COSXY[iPEAK2]])

        B1, B2 = round(beam1, 2), round(beam2, 2)
        print("Beam widths:")
        print("\u03B8 x:", B1)
        print("\u03B8 y:", B2)
        return B1, B2
    


    def getSLL_HPBW_RADIAL(self, IN_PATTERN, RADMIN, HYF=True,NAME="", NORM=True, DEGREES=False, XMIN=None, XMAX=None, YMIN=-60, YMAX=0.0, SHOW=True, DECIMALS=2):
        PATTERN = IN_PATTERN.copy()
        NX, NY = PATTERN.shape
        OFFSET = 0
        LVL = -0.0 #level offset, to find....
        print("\n------------------------------------------------------------------------------------\n")
        if not NORM: #get directivity
            MAX_REF_LEVEL = nanmax(PATTERN)
            OFFSET = 10*log10(MAX_REF_LEVEL)
            if SHOW:
                print("Maximum Integrated Gain (Directivity):",round(OFFSET,2),"dBi")
                print("")
            PATTERN /= nanmax(PATTERN)#normalize the pattern
            if YMAX == 0:
                YMAX = 40
        RADMIN1, RADMIN2, IMAX_X, IMAX_Y = findBeamRadius(PATTERN, NX, NY)
        IX0 = int(NX/2)
        IY0 = int(NY/2)
        if IMAX_X==IX0 and IMAX_Y==IY0:
            print("Steering angle must be different from 90° in elevation")
            return
        # print(IMAX_X, IMAX_Y)

        figpRC, axrc = plt.subplots(1)
        figpRC.suptitle("" + NAME+" Pattern")
        figpRC.set_size_inches(10,3)
        figpRC.canvas.manager.set_window_title('Pattern Radial Cuts')

        axrc.grid()
        axrc.minorticks_on()


        ######*****************************************
        # print(IX0, IY0)
        lsy = []
        ly=[]
        if IMAX_Y!=IY0 and IMAX_X!=IX0:
            m = (IMAX_Y - IY0)/(IMAX_X - IX0)
            b = IMAX_Y- m*IMAX_X
            dIY = 3
            # print("Y= {} X + {}".format(m,b))
            allIY = [i for i in range(NY)]
            allIY = asarray(allIY)

            for i in range(NX):
                y = int(m*i + b)
                if y >= NY: y = NY-1
                if y < 0 : y = 0
                if y not in ly:
                    ly.append(y)
                    lsy.append([y, i])
        elif IMAX_Y==IY0 and IMAX_X!=IX0:
            lsy = [[int(NY/2),i] for i in range(NY)]
        elif IMAX_X==IX0 and IMAX_Y!=IY0:
            lsy = [[i,int(NX/2)] for i in range(NY)]
        iii = asarray(lsy)

        # print(iii.shape)
        PATTERN_RCUT= [10*log10(PATTERN[i[0], i[1]]) for i in iii  ]
        PATTERN_RCUT = asarray(PATTERN_RCUT) + OFFSET
        PATTERN_RCUT = PATTERN_RCUT[~isnan(PATTERN_RCUT)]
        N = len(PATTERN_RCUT)
        COSXY = linspace(-1, 1, num=N)
        ANGLES = degrees(arcsin(COSXY))
        ######*****************************************
        IMAX_R = nanargmax(PATTERN_RCUT)

        if SHOW:
            if DEGREES:
                axrc.plot(ANGLES,PATTERN_RCUT) #
                axrc.set_xlabel('° deg ')
                if XMIN==None and XMAX==None:
                    axrc.set_xlim(-90,90)

                else:
                    axrc.set_xlim(XMIN,XMAX)

            else:
                axrc.plot(COSXY,PATTERN_RCUT) # v = cte

                axrc.set_xlabel('cos')

                if XMIN==None and XMAX==None:
                    axrc.set_xlim(-1.0,1.0)

                else:
                    axrc.set_xlim(XMIN,XMAX)

        axrc.set_ylim(YMIN, YMAX)


        if (IMAX_R - RADMIN2*N/2) > 0:
            SLLR1 = PATTERN_RCUT[0:int(IMAX_R - RADMIN2*N/2)  ]
        else:
            SLLR1 = PATTERN_RCUT[0]

        if (IMAX_R + RADMIN2*N/2) > N:
            SLLR2 = PATTERN_RCUT[N]
        else:
            SLLR2 = PATTERN_RCUT[int(IMAX_R + RADMIN2*N/2): ]

        ##########################################################
        try:
            rSLL1 = nanargmax(SLLR1)
        except:
            rSLL1 = None
        try:
            rSLL2 = nanargmax(SLLR2)
        except:
            rSLL2 = None
        ###############
        print()
        if rSLL1!=None:
            if not isinstance(SLLR1,float):
                SLLRDB1 = round(SLLR1[rSLL1]+ LVL,2)
            else:
                SLLRDB1 = SLLR1
        else:
            SLLRDB1 = -10000
        if rSLL2!=None:
            if not isinstance(SLLR2,float):
                SLLRDB2 = round(SLLR2[rSLL2]+ LVL,2)
            else:
                SLLRDB2 = SLLR2
        else:
            SLLRDB2 = -10000


        ##########################################################

        SLLRDB = max(SLLRDB1, SLLRDB2)

        if SHOW:
            print("Side Lobe Level")
            print("SLL Radial = " , SLLRDB, "dB")

        ##########################################################
        ##########################################################
        #SLL General
        ##########################################################
        ##########################################################

        #masked the main beam
        mask = circularMask(NX, NY,center=(IMAX_X,IMAX_Y), radius=int( NX * RADMIN1/2),  radius2=int( NY * RADMIN2/2))
        maskedPATTERN = PATTERN.copy()
        maskedPATTERN[mask] = NaN
        # plt.figure()
        # cosxy = linspace(-1, 1, num=NX, endpoint=True)
        # plt.pcolormesh(cosxy, cosxy,maskedPATTERN)
        # plt.show()
        IMAX = (nanargmax(maskedPATTERN))
        _x, _y = (IMAX%NX), (IMAX//NY)
        #SLL location in the uv plane
        iSSLx, iSLLy = (_x/NX)*2 -1 , (_y/NX)*2 -1  #to -1 to 1

        SSLOUT = nanmax(maskedPATTERN/nanmax(PATTERN))
        #SSLOUT = PATTERN[ _y, _x]
        G_SLL = round(10*log10(SSLOUT)+LVL,DECIMALS)

        if not HYF:
            return G_SLL

        if SHOW:
            print("SLL: ", G_SLL ,"dB")# #, "\t \t u, v:", iSSLx, iSLLy)
            print("")
            #print(iSSLx, iSLLy)
        ##########################################################
        ##########################################################
        #HPBW
        ##########################################################
        ##########################################################
        PATTERNRDB = PATTERN_RCUT - OFFSET


        locs3dBr1 = where(PATTERNRDB >= -3)[0][0]
        locs3dBr2 = N - where(flip(PATTERNRDB) >= -3)[0][0]

        steerAngle = ANGLES[IMAX_R]
        #print(steerAngle)
        m1 = (0-PATTERNRDB[locs3dBr1])/(steerAngle-ANGLES[locs3dBr1])
        m2 = (0-PATTERNRDB[locs3dBr2])/(steerAngle-ANGLES[locs3dBr2])
        x3db1 = -3/m1 + steerAngle
        x3db2 = -3/m2 + steerAngle
        #print(x3db1, x3db2)
        #HPBWR = round(abs(ANGLES[locs3dBr1] - ANGLES[locs3dBr2]),2)
        HPBWR = round(abs(x3db1 - x3db2),2)
        if SHOW:
            print("Half Power Beam Width")
            print("HPBW R: ",HPBWR )
            print("")

        ##########################################################
        ##########################################################
        #FNBW
        ##########################################################
        ##########################################################


        rDB1 = PATTERN_RCUT[0:IMAX_R]
        rDB2 = PATTERN_RCUT[IMAX_R:]

        try:
            i_FNBWR_l = self.getFirstNull(rDB1, reverse=True)
        except:
            i_FNBWR_l = 0
        try:
            i_FNBWR_h = self.getFirstNull(rDB2, reverse=False)
        except:
            i_FNBWR_h = N - 1
            IMAX_R = 0
        # print(i_FNBWR_l, i_FNBWR_h)
        FNBWR_l = ANGLES[i_FNBWR_l]
        FNBWR_h = ANGLES[i_FNBWR_h + IMAX_R]
        # print(FNBWR_l, FNBWR_h)
        FNBWR = round(abs(FNBWR_h-FNBWR_l),2)

        if SHOW:
            print("First Null Beam Width")
            print("FNBW R: ",FNBWR)
        plt.show()
        return G_SLL, SLLRDB, FNBWR, _x, _y, [ANGLES,PATTERN_RCUT]

        ################################################################################
        ################################################################################
        ################################################################################

    def minimizeSLL(self, radmin, Xlocs=[], Ylocs=[],  gain=0.01, freq=50, maxIters=100, spacing=5,
                    maxRadius=50,maxSide=None, maxRadius2=None, plotPattern=False, save=None,Elev=90.0, Azi=0.0,
                    resolutionPatternPoints=512, minDB=-50, antenna="cosine", EF=1.2, GND=True, MOD_Z=1.0):
        """
        Xlocs, Ylocs    : antenna locations in meters from the center
        radmin          : main beam radius
        """
        if len(Xlocs)==0 and len(Ylocs)==0:
            Xlocs = self.Xdata
            Ylocs = self.Ydata

        c = 2.997962*10**8

        NITER = maxIters
        Z = MOD_Z
        LAMBDA =  c/(freq*10**6)  #wavelenght (m)
        print("LAMBDA ",LAMBDA )

        BYRADIUS = False
        BYELLIPSE = False
        MAXSIDE = 0
        MAXRADIUS = 0
        MAXRADIUS2 = 0
        if maxRadius2!=None:
            MAXRADIUS = maxRadius
            MAXRADIUS2 = maxRadius2
            BYELLIPSE = True
        elif maxRadius!=None:
            MAXRADIUS = maxRadius
            BYRADIUS = True
        else:
            if maxSide!=None:
                MAXSIDE = maxSide
                BYRADIUS = False
            else:
                MAXRADIUS = max(Xlocs.max(),Ylocs.max())

        RADMIN = 1.01*abs(radmin)/2
        DLIMIT = spacing
        X = Xlocs.copy()/LAMBDA
        Y = Ylocs.copy()/LAMBDA
        NARRAY = len(Xlocs)
        GAIN = gain
        print("NORM GAIN ",gain/NARRAY**2)
        Cx0 = cos(radians(Elev))* sin(radians(Azi))
        Cy0 =  cos(radians(Elev))* cos(radians(Azi))

        if len(Xlocs) != len(Ylocs):
            print("X and Y coordinates must be the same length!")
            return None, None
        ##--------------------------------------------------------------------------
        Nx = resolutionPatternPoints
        Ny = resolutionPatternPoints
        dcx = (2.0*Z)/Nx
        dcy = (2.0*Z)/Ny
        # cosx = dcx*(arange(Nx)-floor(Nx/2))
        # cosy = dcy*(arange(Ny)-floor(Ny/2))
        cosx = linspace(-Z, Z, num=Nx, endpoint=True)
        cosy = linspace(-Z, Z, num=Ny, endpoint=True)
        Cx, Cy = meshgrid(cosx, cosy)
        Cz2 = 1-(Cx**2+Cy**2)
        Cz2[Cz2<=0] = NaN
        Cz = sqrt(Cz2)
        domg = 1./Cz
        domg[~isfinite(domg)] = 0
        ##--------------------------------------------------------------------------
        if Z!=1:
            antenna=""
        FANT = 0
        if antenna=="dipole":
            FANT = dipoleAntennaFactor(Cy,Cz,ground=GND)
        elif antenna=="cosine":
            FANT = cosineAntennaFactor(Cx,Cy,EF=EF)
        else:
            print("Missing antenna element")
            FANT = 1
        ##--------------------------------------------------------------------------
        Xbest = zeros( NARRAY )
        Ybest = zeros( NARRAY )
        DX = zeros(NARRAY)
        DY = zeros(NARRAY)
        itBest=0
        SLLBEST = 0
        it = 0
        WXY_TRY = 0
        WZENX, WZENY = 0,0
        WZENX_LAST, WZENY_LAST = 1000,1000
        minSLL = 10000

        fig, axs = plt.subplots(1, 2, figsize=(14, 6))
        plt.ion()
        fig.suptitle('Optimization process', fontsize=16)
        #plt.axes().set_aspect('equal')
        while(it <= NITER):
            if WXY_TRY > 100:
                print("Break, same WZENX and WZENY it:{}".format(it))
                break
            FARR = 0
            #array factor
            for IARRAY in range(NARRAY):
                FARR += exp(1j*2*PI*X[IARRAY]*(Cx-Cx0) + 1j*2*PI*Y[IARRAY]*(Cy-Cy0) )
            #gain array
            #antenna pattern
            FARR = FARR*FANT
            PATT = (abs(FARR))**2
            # # mask array from -1 to 1 or limited from -z to +z
            maskUV = circularMask(Nx, Ny, radius=int(Nx*0.5))
            PATT[~maskUV] = NaN
            BEAZEN = PATT/nanmax(PATT)
            #mask the main beam
            #MAXVAL = nanargmax(BEAZEN)
            #IMAX_X, IMAX_Y = int(MAXVAL%Nx), (MAXVAL//Ny)
            IMAX_X, IMAX_Y =int(((Cx0+Z)/(2*Z))*Nx ) , int( ((Cy0+Z)/(2*Z))*Ny )
            mask = circularMask(Nx, Ny,center=(IMAX_X,IMAX_Y), radius=int(Nx * RADMIN))
            maskedPATT = PATT.copy()
            maskedPATT[mask] = NaN
            #maximum sidelobe location
            IMAX = (nanargmax(maskedPATT))
            _x, _y = (IMAX%Nx), (IMAX//Ny)
            #worst direction in the uv plane
            WZENX, WZENY = (_x/Nx)*(2*Z) -Z , (_y/Nx)*(2*Z) -Z  #to -1 to 1  or -z to +z
            if round(WZENX,2) == WZENX_LAST and round(WZENY,2) == WZENY_LAST:
                WXY_TRY += 1
            else:
                WZENX_LAST = round(WZENX,2)
                WZENY_LAST = round(WZENY,2)
                WXY_TRY = 0

            #sidelobe level
            #store the best distribution
            newSLL = nanmax(maskedPATT/nanmax(PATT))
            #newSLL = BEAZEN[_x, _y]
            #print(nanmax(BEAZEN))
            if it == 0:
                minSLL = newSLL
            SLL = round(10*log10(newSLL),3)
            if newSLL <= minSLL :

                minSLL = newSLL
                Xbest = X.copy()
                Ybest = Y.copy()
                itBest = it
                SLLBEST = SLL

            if plotPattern:

                axs[0].clear()
                axs[1].clear()

                graph1 = axs[0].pcolormesh(cosx,cosy,10*log10(BEAZEN),shading='auto',cmap='jet', vmax=0, vmin=minDB)
                #axs[0].clim([-60,0])
                circle=plt.Circle( ( round(Cx0,4),round(Cy0,4) ),RADMIN*2,  color='w', fill=False)
                axs[0].add_patch(circle)
                axs[0].scatter(WZENX, WZENY, marker='x', color='k')
                title = "Correcting positions, norm gain: {:3.4f}, it: {}, SLL: {:2.2f}".format(GAIN/NARRAY**2,it,SLL )
                axs[0].set_title (title)
                #axs[0].set_aspect('equal')

                graph2 = axs[1].scatter(X*LAMBDA, Y*LAMBDA,marker='x')
                axs[1].minorticks_on()
                axs[1].grid(b=True, which='major', color='0.8', linestyle='-')
                axs[1].grid(b=True, which='minor', color='0.8', linestyle='--')
                #axs[1].set_aspect('equal')
                if BYELLIPSE:
                    axs[1].set_title ("Antenna locations restricted to {} and {} m of radius".format(MAXRADIUS,MAXRADIUS2 ))
                elif BYRADIUS:
                    axs[1].set_title ("Antenna locations restricted to {} m of radius".format(MAXRADIUS))
                else:
                    axs[1].set_title ("Antenna locations restricted to {} m of side".format(MAXSIDE))
                axs[1].set_ylabel("meters in Y (v) direction")
                axs[1].set_xlabel("meters in X (u) direction")

                if BYELLIPSE:
                    ellipse=Ellipse((0,0),2*MAXRADIUS, 2*MAXRADIUS2, color='k', fill=False)
                    axs[1].add_patch(ellipse)
                elif BYRADIUS:
                    circle=plt.Circle((0,0),MAXRADIUS,  color='k', fill=False)
                    axs[1].add_patch(circle)
                else:
                    rect = plt.Rectangle((-MAXSIDE,-MAXSIDE),2*MAXSIDE, 2*MAXSIDE,  color='k', fill=False)
                    axs[1].add_patch(rect)

               
                fig.canvas.draw()
                
                plt.pause(0.001)
            

            for IARRAY in range(NARRAY):
                ANUMER = 0
                for KARRAY in range(NARRAY):
                    ANUMER += sin(2 * PI *( WZENX*(X[KARRAY]- X[IARRAY]) + WZENY*(Y[KARRAY]- Y[IARRAY]) ) )
                DPPX = 2 * 2*PI * WZENX * ANUMER / (NARRAY * NARRAY)
                DPPY = 2 * 2*PI * WZENY * ANUMER / (NARRAY * NARRAY)
                #
                DX[IARRAY] = - GAIN * DPPX
                DY[IARRAY] = - GAIN * DPPY

            if DX.mean()> 0.1 or DY.mean()>0.1:
                print("To much gain value!!")
                print("Displacements dX:{}, dY:{} ".format(DX.mean(), DY.mean()))

            #check boundaries
            for IARRAY in range(NARRAY):
                XTEM = X[IARRAY] + DX[IARRAY]
                YTEM = Y[IARRAY] + DY[IARRAY]

                RT = sqrt(  (XTEM*LAMBDA)**2 + (YTEM*LAMBDA)**2) #in meters
                PE = 0
                if BYELLIPSE:
                    PE = ((XTEM*LAMBDA)**2 /MAXRADIUS**2) + ((YTEM*LAMBDA)**2/MAXRADIUS2**2)
                if  PE > 1 and BYELLIPSE:
                    continue
                if RT > MAXRADIUS and BYRADIUS:
                    continue
                if not BYRADIUS and not BYELLIPSE and (abs(XTEM*LAMBDA) > MAXSIDE or abs(YTEM*LAMBDA) > MAXSIDE) :
                    continue

                flag_toclose = False
                for KARRAY in range(NARRAY):
                    if (KARRAY == IARRAY):
                        continue
                    NEWDIS = sqrt( ((XTEM - X[KARRAY])*LAMBDA)**2 + ((YTEM - Y[KARRAY])*LAMBDA)**2) #in meters
                    if  NEWDIS <= (DLIMIT):
                        #print("x, y, r, d ", XTEM,YTEM,NEWDIS,DLIMIT)
                        flag_toclose = True
                        break
                if not flag_toclose:
                    X[IARRAY] = XTEM
                    Y[IARRAY] = YTEM

            it += 1


            if save != None:
                figname =save+"/{}".format(it)
                plt.savefig(figname+'.png')

        #final pattern

        FARR = 0
        for IARRAY in range(NARRAY):
                FARR += exp(1j*2*PI*Xbest[IARRAY]*(Cx-Cx0) + 1j*2*PI*Ybest[IARRAY]*(Cy-Cy0))

        FARR = FARR*FANT
        #gain array
        PATT = abs(FARR)**2
        maskUV = circularMask(Nx, Ny, radius=int(Nx*0.5))
        PATT[~maskUV] = NaN
        BEAZEN = PATT/nanmax(PATT)
        fig = plt.figure("Optimized Pattern", figsize=(6, 6))
        #plt.axes().set_aspect('equal')
        im = plt.pcolormesh(cosx,cosy,10*log10(BEAZEN),shading='auto',cmap='jet')
        plt.title ("Best in {} iterations with gain={}, best iter ={}, SLL={}".format(NITER,GAIN/NARRAY**2, itBest, SLLBEST))
        plt.colorbar(im,fraction=0.046, pad=0.04)

        plt.clim([minDB,0])
        plt.show() 
        return minSLL, Xbest*LAMBDA, Ybest*LAMBDA, BEAZEN, X*LAMBDA, Y*LAMBDA
