#!/usr/bin/env python

import sys
import os
import getopt
import shutil

import numpy
from matplotlib import pyplot,rc,rcdefaults,cm
from math import exp,log, log10, cos, sin, acos, asin
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
import scipy.ndimage as ndimage
import matplotlib.colors as colors

# custom-made functions
#----------------------------------------------------------------------
def open_data(inpath,filelist):

    Dict = {}

    for i,namefile in enumerate(filelist):

        # open file, read all lines
        inputpath = os.path.join(inpath,namefile)
        f=open(inputpath,'rU') 
        lines=f.readlines()
        f.close()

        # storing headers in list headerow
        headerow=lines[0].strip().split()

        # deleting rows that are not data (first and last rows of the file)
        del lines[0]

        # transforming data from string to float type, storing it in array 'data'
        datafloat=[]
        for line in lines:
            datafloat.append(line.strip().split())
        data=numpy.array(datafloat,float)

        # creating one dictionnary and storing the float data in it
        dictnamelist= {}
        for j,varname in enumerate(headerow):
            dictnamelist[varname]=data[:,j]
        Dict[namefile] = dictnamelist

    return Dict

#----------------------------------------------------------------------
def make_pcolormesh_box_boundaries(arrayX, arrayY):

    stepx = (arrayX[1]-arrayX[0])/2.
    stepy = (arrayY[1]-arrayY[0])/2.

    boundX = arrayX-stepx
    boundX = numpy.append(boundX, boundX[len(boundX)-1]+2.*stepx)
    boundY = arrayY-stepy
    boundY = numpy.append(boundY, boundY[len(boundY)-1]+2.*stepy)
   
    return boundX, boundY

#----------------------------------------------------------------------
def remake_axis(var, rangevar, newvar, resultsfile, Dict):

#    print var, rangevar, newvar, resultsfile, Dict.keys()
    W = [0.]*len(rangevar)

    for i,item in enumerate(rangevar):
        V = []
        for j,newitem in enumerate(Dict[resultsfile][newvar]):
            if (abs(Dict[resultsfile][var][j]-item) < 1e-6):
                V = V + [newitem]
        print item, V  # uncomment this line if you want to check the result of the function
        W[i] = numpy.mean(numpy.array(V))

    Dict[resultsfile]['grouped_%s'%newvar] = [0.]*len(Dict[resultsfile][var])

    for i,item in enumerate(rangevar):
        for j,newitem in enumerate(Dict[resultsfile][newvar]):
            if (abs(Dict[resultsfile][var][j]-item) < 1e-6):
                Dict[resultsfile]['grouped_%s'%newvar][j] = W[i]

    return Dict[resultsfile]['grouped_%s'%newvar], W
    
#----------------------------------------------------------------------
def open_output(path,inputfoldernamelist):
    """
Opening results file
    """
    Dict = {}
    
    for i,fil in enumerate(inputfoldernamelist):
        inputpath = os.path.join(path,fil)
        f=open(inputpath,'rU') 
        lines=f.readlines()
        f.close()
        ################ header row ################
        if (fil.endswith('res.dat')):
            headerow=lines[6].strip().split()
            del lines[0:8]
        else: 
            headerow=lines[0].strip().split()
            del lines[0]
        #### last lines with text or empty lines ####
        if (fil.startswith('output_dyn') or fil.endswith('output_dyn') or fil.startswith('output_sca') or fil.endswith('output_sca')): 
            if lines[-1].startswith('Saturation level'): del lines[-1:]
        if (fil.endswith('res.dat')):
            del lines[-4:]
        datafloat=[]
        for line in lines:
            datafloat.append(line.strip().split())
        data=numpy.array(datafloat, dtype=float)
        dictnamelist = {}
        for j,varname in enumerate(headerow):
            dictnamelist[varname]=data[:,j]
        Dict[fil] = dictnamelist

    return Dict

#----------------------------------------------------------------------
def open_csv(inpath,namefilelist):

    import csv

    Dict = {}

    for i,namefile in enumerate(namefilelist):
        inputpath=os.path.join(inpath,namefile)
        f=open(inputpath,'rU')
        reader=csv.reader(f, delimiter=',', skipinitialspace=True)
        all=[]
        for row in reader:
            all.append(row)
        headerow=all[0]
        del all[0]
        datafloat=[]
        for row in all:
            datafloat.append(map(float,row))
        data=numpy.array(datafloat)
        dictnamelist = {}
        for j,varname in enumerate(headerow):
            dictnamelist[varname]=data[:,j]
        Dict[namefile] = dictnamelist

    return Dict

#----------------------------------------------------------------------
def sigmoid(x, slope):
    y = 3.8 / (1+ numpy.exp(-slope*(x-185.)))
    return y
    
#----------------------------------------------------------------------
def getth(daynr,lat,lon,thour):

    pi   = acos(-1.)
    piby = pi/ 180.
    lonr = lon*piby
    latr = lat*piby
    obliq = 23.45 * piby
    deday = 4.88 + 2*pi/365* daynr
    delta = asin(sin(obliq)*sin(deday))
    houra = lonr - pi + thour* (2.*pi/24.)
    angle = acos(sin(delta)*sin(latr) +  cos(delta)*cos(latr)*cos(houra))

    return angle

#----------------------------------------------------------------------
def oneD_post_process(inpath,folderlist):
    
    Dict = {}
    Var   = 'wg' 
    Range = numpy.arange(0.06,0.1501,0.0009)
    Range2 = numpy.arange(0.,1.005,0.01)

    for i,folderr in enumerate(folderlist):

        print "\nOpening %s\n"%folderr
        analysisdir      = os.path.join(inpath, folderr)
	listofnamoptions = [f for f in os.listdir( analysisdir ) if ( f.startswith('nam') and (Var[0:2] in f) ) ]

        smi = []
        gs  = []
        LE  = []
        NPP = []
        ef  = []
        wue = []
	
        for i,tag in enumerate(listofnamoptions):
	
            sensdir = 'SENS_'+tag[11:16]
	    print "Reading %s"%sensdir
	    smi     = smi + [Range2[i]]
            Dummy   = open_output(os.path.join(analysisdir,sensdir),['output_sca','output_dyn'])
	    
	    ID_12   = numpy.where(Dummy['output_dyn']['UTC(hours)']==12.)[0][0]

	    gs      = gs  + [1000./Dummy['output_sca']['rs(s.m-1)'][ID_12]]
	    
	    LE      = LE  + [Dummy['output_dyn']['LE(W.m-2)'][ID_12]]	    
	    NPP     = NPP + [-Dummy['output_sca']['An(mgCO2.m-2.s-1)'][ID_12]]
	    
	    evafra = Dummy['output_dyn']['LE(W.m-2)']/(Dummy['output_dyn']['LE(W.m-2)']+Dummy['output_dyn']['SH(W.m-2)'])
	    ef      = ef  + [evafra[ID_12]]

            NEE     =  Dummy['output_sca']['An(mgCO2.m-2.s-1)'] + Dummy['output_sca']['Resp(mgCO2.m-2.s-1)']
            WUEeco  = -0.001* 12./44.* NEE / (Dummy['output_dyn']['LE(W.m-2)'] / 2.5e6)
            wue     = wue + [numpy.mean(WUEeco)]

        Dict[folderr] = {}
	
        Dict[folderr]['SMI'] = smi
        Dict[folderr]['gs']  = gs
        Dict[folderr]['LE']  = LE
        Dict[folderr]['NPP'] = NPP
        Dict[folderr]['EF']  = ef
        Dict[folderr]['WUE'] = wue
	
    return Dict


#######################################################################
########################### START OF SCRIPT ###########################

if __name__=="__main__":


    try:                                
        opts, args = getopt.getopt(sys.argv[1:], "-h")
    except getopt.GetoptError:           
        print "Error"
        sys.exit(2)

    for options in opts:
        options=options[0].lower()
        if options == '-h':
            helptext = """
This script plots the figures of the methods section of Article 2

You do not need to specify any arguments to run this script:
./Article2_methods_figures.py
                
                """
            
            print helptext
            
            sys.exit(2)


    ######### TO USER: define which figures you want to plot:#################

    # Article no. 2 figures:
    
    figure_2  = False
    figure_3  = False
    figure_4a = False
    figure_4b = False
    figure_5a = True
    figure_5b = True
    figure_6a = False
    figure_6b = False
    figure_7a = True
    figure_7b = True
    figure_8  = False
    figure_9a = False
    figure_9b = False
    figure_9c = False
    
    # EGU poster figures:
    
    figure_10  = False
    figure_11a = False
    figure_11b = False
    figure_12  = False
    figure_13a = False
    figure_13b = False
    figure_13c = False
    figure_14a = False
    figure_14b = False
    
    ##########################################################################
    
    
    if (figure_2==False and figure_3==False and figure_4a==False and figure_4b==False and figure_5a==False 
        and figure_5b==False and figure_8==False and figure_9a==False and figure_9b==False and figure_9c==False and figure_10==False and figure_11a==False and figure_11b==False and figure_12==False and figure_13a==False and figure_13b==False and figure_13c==False and figure_14a==False and figure_14b==False): 
        print "\nYou didn't ask to plot any figure! Open the script and change that.\n"

    # Working directories
    currentdir = os.getcwd()
    csvdatadir = '/Users/mariecombe/Modeling/MXL/MLmodel/Data_files/'
    storagedir = '/Storage/MXL/Marie/MXL-Ags_Sensitivity'

    # First close all opened figures in python
    pyplot.close('all')

    # general plot properties
    rcdefaults()

    rc('xtick', labelsize=16)
    rc('xtick.major', size=6)
    rc('ytick', labelsize=16)
    rc('ytick.major', size=6)
    rc('contour', negative_linestyle = 'solid')
    rc('lines', linewidth=1)
    rc('axes', labelsize=16, linewidth=0.5)
    rc('legend', frameon=False, handlelength=3., fontsize=16)
    rc('grid', linewidth=0.)

    # Ranges of variation of the drivers 
    
    RangeDict = {'thetam0': numpy.arange(284.,290.01,0.06), 
                 'wg': numpy.arange(0.060,0.363,0.003),
                 'cc': numpy.arange(0.000,0.303,0.003), 
                 'ga': numpy.arange(0.002,0.00806,0.00006),
                 'smi': numpy.arange(0.0,1.001,0.01), #numpy.arange(0.4,0.601,0.00667), #numpy.arange(0.4,0.601,0.002) 
                 'smi2': numpy.arange(0.4,0.601,0.002), 
		 'curvature': numpy.arange(0.,100.1,1.),
		 'curvature2': numpy.linspace(0.,100.,21.),
		 'wsls': numpy.arange(0.,0.0000401,0.0000004),
		 'deltatheta'  : numpy.arange(0.5,5.01,0.045),
	         'gammatheta'  : numpy.arange(0.002,0.00806,0.00006), 
		 'day': numpy.arange(130.,230.5,1.), #numpy.arange(200.,230.5,1.)
		 'day2': numpy.linspace(216.,236.,21.)
		 #'gammatheta': numpy.arange(0.002,0.00806,0.00006)
		 } ;


    if (figure_2 == True): #if you want to plot figure 1

        wwp = 0.0
        wfc = 1.0
        wsat = 3.3

        # open a figure
        fig = pyplot.figure('stress')
        fig.subplots_adjust(0.12,0.12,0.95,0.97,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
    
        # calculate the water stress response for each P4
        x = numpy.arange(-0.2,2.001,0.008)
    
        # defining the levels of P4
        X = numpy.arange(2.17,4.01,0.0366)
        P3 = [0.]*len(X)
        for i,val in enumerate(X):
	    P3[i] = (2.**val)-0.999999
        P0 = numpy.arange(0.000001,1.6,0.064)
        P1 = numpy.arange(1.6,3.5,0.076) 
        P4 = numpy.concatenate((P0,P1,P3), axis=0)   
    
        for param in P4[::5]:
            y = [0.]*len(x)
            for i,val in enumerate(x):
                if (val <= wwp):
                    y[i] = 0.
                elif (val >= wfc):
                    y[i] = 1.
                else:
                    y[i] = ( 1-exp(-(param)*(val-wwp)/(wfc-wwp)) )/( 1-exp(-(param)) )
    
        # plot each water stress function
            pyplot.plot(x,y,lw=1.5,c='k',label='P4 = %5.2f'%(param))
	
        # finish the figure
        #pyplot.legend(loc='lower right', labelspacing=0.2)
        pyplot.ylim(0.,1.1)
        pyplot.xlim(-0.2,1.2)
        pyplot.xticks([0.,0.2,0.4,0.6,0.8,1.,1.2])
        pyplot.fill_between(x, 1.0, 1.2, where=(x<=0.), color='red', alpha=0.3)
        pyplot.annotate('Dry', xy=(0.07, 0.95), xycoords='axes fraction',fontsize=20,horizontalalignment='center', verticalalignment='center')
        pyplot.fill_between(x, 1.0, 1.2, where=((x>0.)&(x<1.)), color='yellow', alpha=0.5)
        pyplot.annotate('Water stressed', xy=(0.5, 0.95), xycoords='axes fraction',fontsize=20,horizontalalignment='center', verticalalignment='center')
        pyplot.fill_between(x, 1.0, 1.2, where=(x>=1.), color='blue', alpha=0.2)
        pyplot.annotate('Wet', xy=(0.93, 0.95), xycoords='axes fraction',fontsize=20,horizontalalignment='center', verticalalignment='center')
        pyplot.axvline(x=0., lw=1.2, color='black', linestyle='--')
        pyplot.axvline(x=1., lw=1.2, color='black', linestyle='--')
        pyplot.ylabel(r'$\beta$', fontsize=22)
        pyplot.xlabel('SMI (-)', fontsize=20)
        #pyplot.title('Water stress responses\n', fontsize=18)
    
        savepath = os.path.join(currentdir,'f02.png')
        fig.savefig(savepath,dpi=300)


    if (figure_3 == True): #if you want to plot figure 3

        # opening LAI data files 
    
        CropObsDict          = open_csv(csvdatadir,['LAI_obs_6-aug-2007.csv','MeanT150cm_obs_clim_1951-2000.csv',
                                                     'MinT150cm_obs_clim_1951-2000.csv','MaxT150cm_obs_clim_1951-2000.csv'])
        

        # fitting a sigmoid curve through the observed LAI data 
       
        xA=CropObsDict['LAI_obs_6-aug-2007.csv']['DOY']
        yA=CropObsDict['LAI_obs_6-aug-2007.csv']['LAI']
        popt, pcov = curve_fit(sigmoid, xA, yA)    

        # interpolating the mean monthly temperature 
       
        xB=CropObsDict['MeanT150cm_obs_clim_1951-2000.csv']['DOY']
        yB=CropObsDict['MeanT150cm_obs_clim_1951-2000.csv']['average_T']-2.
        thetam_interp = UnivariateSpline(xB, yB, s=1, k=3)

        x3=numpy.linspace(130.,230.,100.)

        LAI = sigmoid(x3, *popt)
        CVEG = 1.- numpy.exp(-LAI)    
        THETAM0 = thetam_interp(x3)+273.15
        ALPHA = CVEG * 0.198 + (1-CVEG) * 0.15

        day = 130
        lat = 51.59
        lon = 5.38
        thour = 12.
        cc = 0.225 
    
        SWIN = []
        SWNET = []
        for i,day in enumerate(x3):
            angle = max(0.0,cos(getth(1.0*day,lat,lon,thour)))
            Tr    = (0.6 + 0.2 * angle) * (1 - 0.4 * cc)
	    insw  = 1368. * Tr * angle
	    netsw = insw * (1-ALPHA[i])
            SWIN  = SWIN + [insw]
	    SWNET = SWNET + [netsw]
    

        # open a figure
        fig = pyplot.figure('cycle', figsize=(9,7))
        fig.subplots_adjust(0.1,0.1,0.97,0.97,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
    
        ax = fig.add_subplot(2,2,1)
        pyplot.scatter(xA,yA,edgecolor='k',facecolor='k')
        pyplot.plot(x3,LAI,lw=2,ls='-',c='k',label='LAI')
        pyplot.plot(x3,CVEG,lw=2,ls='--',c='k',label=r'$f_{veg}$')
        ax.legend(loc='upper left', ncol=1)
        pyplot.xlim(130.,230.)
        pyplot.ylim(0.,4.)
        pyplot.yticks([0.,1.,2.,3.,4.])
        pyplot.annotate('(a)', xy=(0.92, 0.08), xycoords='axes fraction', fontsize=16, horizontalalignment='center', verticalalignment='center')
        ax.set_ylabel(r'Index (-)', labelpad=7)
    
        ax = fig.add_subplot(2,2,2)
        pyplot.plot(x3,ALPHA,lw=2,ls='-',c='k')
        pyplot.xlim(130.,230.)
        pyplot.ylim(0.14,0.21)
        pyplot.yticks([0.15,0.16,0.17,0.18,0.19,0.20])
        pyplot.annotate('(b)', xy=(0.92, 0.08), xycoords='axes fraction', fontsize=16, horizontalalignment='center', verticalalignment='center')
        ax.set_ylabel(r'$\alpha$ (-)', fontsize='18', labelpad=5)

        ax = fig.add_subplot(2,2,3)
        pyplot.plot(x3,THETAM0,lw=2,ls='-',c='k')
        pyplot.scatter(xB,yB+273.15,edgecolor='k',facecolor='k')
        pyplot.xlim(130.,230.)
        pyplot.ylim(283.,289.)
        pyplot.yticks([283.,284.,285.,286.,287.,288.])
        pyplot.annotate('(c)', xy=(0.92, 0.08), xycoords='axes fraction', fontsize=16, horizontalalignment='center', verticalalignment='center')
        ax.set_xlabel('time (DOY)', labelpad=7)
        ax.set_ylabel(r'$\theta_0$ (K)', labelpad=3)

        ax = fig.add_subplot(2,2,4)
        pyplot.plot(x3,SWIN,lw=2,ls='-',c='k',label=r'SW$_{in}$')
        pyplot.plot(x3,SWNET,lw=2,ls='--',c='k',label=r'SW$_{net}$')
        ax.legend(loc='lower left', ncol=1)
        pyplot.xlim(130.,230.)
        pyplot.ylim(550.,900.)
        pyplot.yticks([600.,650.,700.,750.,800.,850.])
        pyplot.annotate('(d)', xy=(0.92, 0.08), xycoords='axes fraction', fontsize=16, horizontalalignment='center', verticalalignment='center')
        ax.set_xlabel('time (DOY)', labelpad=7)
        ax.set_ylabel(r'max. SW (W m$^{-2}$)', labelpad=5)


        savepath = os.path.join(currentdir,'f03.png')
        fig.savefig(savepath,dpi=300)


    # Reboot general plot properties
    rcdefaults()

    seq_dash = [3,4]
    seq_dash1 = [6,3]
    seq_dash2 = [4,2,1,1,1,2]
    seq_dash4 = [1,1]#[3,2,1,1,1,1,1,1,1,2]
    seq_dash3 = [3,2]

    rc('xtick', labelsize=8)
    rc('xtick.major', size=3)
    rc('ytick', labelsize=8)
    rc('ytick.major', size=3)
    rc('contour', negative_linestyle = 'solid')
    rc('lines', linewidth=1)
    rc('axes', labelsize=8, linewidth=0.5)
    rc('legend', frameon=True, handlelength=3., fontsize=8)
    rc('grid', linewidth=0.)


    
    if (figure_4a == True): #if you want to plot figure 4a

        #results_filename = '2DIM_SA_thetam0_284-290_cc_0-30_SMI=0.8/Results_thetam0_cc.dat'
        #ResDict = open_data(storagedir,[results_filename])
        results_filename = '2DIM_SA_thetam0_284-290_cc_0-30_SMI=0.55/Results_thetam0_cc.dat'
        ResDict = open_data(storagedir,[results_filename])

        # defining X, Y, Z variable names, and the results file where these will be plotted from
        xvar = 'thetam0'
        yvar = 'cc'
        zvar1 = 'NEE_int'
        zvar2 = 'wce_int'

        # Rebuilding Y axis to produce the new Y range
        #newaxis = remake_axis('cc', numpy.arange(0.000,0.303,0.003), 'SWnet_12h', results_filename, ResDict)        
        #RangeDict[yvar] = newaxis[1]

        # values of the X, Y, Z variables for plotting
        X = list(RangeDict[xvar])
        Y = list(RangeDict[yvar])
        X, Y = numpy.meshgrid(X,Y)

        Z1 = numpy.reshape(ResDict[results_filename][zvar1], (-1,101))
        Z1 = Z1.T
        Z1 = numpy.ma.masked_invalid(Z1)

        Z2 = numpy.reshape(ResDict[results_filename][zvar2], (-1,101))
        Z2 = Z2.T
        Z2 = numpy.ma.masked_invalid(Z2)
        Z2 = ndimage.gaussian_filter(Z2, sigma=1.0, order=0)

        bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar])
 
        # plotting

        barlim=numpy.linspace(-71.,-59.,200) # NEE colorbars for high SMI
        barnorm = colors.BoundaryNorm(barlim,256)
        barlim=[-70.,-68.,-66.,-64.,-62.,-60.]  # NEE colorbar ticks

        fig = pyplot.figure('CO2_budgeta', figsize=(3.2,2.41))
        fig.subplots_adjust(0.15,0.16,0.90,0.96,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
        pyplot.pcolormesh(bound[0], bound[1], Z1) #norm=barnorm, cmap=cm.jet_r)
        cbar = pyplot.colorbar()#ticks=barlim)
        cbar.set_label(r'NEE (g$_{\mathrm{CO2}}$ m$^{-2}$ d$^{-1}$)') 
        pyplot.annotate('(a)', xy=(1.2, -0.15), xycoords='axes fraction', annotation_clip=False, size=10, color='black', va="center", ha="center")

        pyplot.xlim(min(bound[0]),max(bound[0]))
        pyplot.ylim(min(bound[1]),max(bound[1]))
        pyplot.xticks([284.,285.5,287.,288.5,290.])

        levels = [-21.,-20.,-19.,-18.,-17.,-16.,-15.,-14.,-13.,-12.,-11.,-10.]
        cs= pyplot.contour(X, Y, Z2, colors='k', linewidths=0.5) 
        pyplot.clabel(cs,  fontsize=8, fmt='%.0f') 

        pyplot.xlabel(r'$\theta$$_0$ (K)')
        pyplot.ylabel(r'max. SW$_{\mathrm{net}}$ (W m$^{-2}$)')

        savepath = os.path.join(currentdir,'f04a.jpeg')
        fig.savefig(savepath,dpi=1000)


    if (figure_4b == True): #if you want to plot figure 4b

        #results_filename = '2DIM_SA_thetam0_284-290_cc_0-30_SMI=0.8/Results_thetam0_cc.dat'
        #ResDict = open_data(storagedir,[results_filename])
        results_filename = '2DIM_SA_thetam0_284-290_cc_0-30_SMI=0.55/Results_thetam0_cc.dat'
        ResDict = open_data(storagedir,[results_filename])

        # defining X, Y, Z variable names, and the results file where these will be plotted from
        xvar = 'thetam0'
        yvar = 'cc'
        zvar1 = 'co2_min'
        zvar2 = 'h_max'
        
        # Rebuilding Y axis to produce the new Y range
        #newaxis = remake_axis('cc', numpy.arange(0.000,0.303,0.003), 'SWnet_12h', results_filename, ResDict)        
        #RangeDict[yvar] = newaxis[1]

        # values of the X, Y, Z variables for plotting
        X = list(RangeDict[xvar])
        Y = list(RangeDict[yvar])
        X, Y = numpy.meshgrid(X,Y)

        Z1 = numpy.reshape(ResDict[results_filename][zvar1], (-1,101))
        Z1 = Z1.T
        Z1 = numpy.ma.masked_invalid(Z1)
	Z1 = ndimage.gaussian_filter(Z1, sigma=1.0, order=0)

        Z2 = numpy.reshape(ResDict[results_filename][zvar2], (-1,101))
        Z2 = Z2.T
        Z2 = numpy.ma.masked_invalid(Z2)
        Z2 = ndimage.gaussian_filter(Z2, sigma=1.0, order=0)

        bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar])
 
        # plotting

        barlim=numpy.linspace(344.,349.,200) # cmin colorbars for high SMI
        barnorm = colors.BoundaryNorm(barlim,256)
        barlim=[344.,345.,346.,347.,348.,349.]  # cmin colorbar ticks

        fig = pyplot.figure('CO2_budgetb', figsize=(3.2,2.41))
        fig.subplots_adjust(0.15,0.16,0.90,0.96,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
        pyplot.pcolormesh(bound[0], bound[1], Z1)#, norm=barnorm, cmap=cm.jet_r)
        cbar = pyplot.colorbar()#ticks=barlim)
        cbar.set_label(r'Diurnal CO$_2$ minimum (ppm)') 
        pyplot.annotate('(b)', xy=(1.2, -0.15), xycoords='axes fraction', annotation_clip=False, size=10, color='black', va="center", ha="center")

        pyplot.xlim(min(bound[0]),max(bound[0]))
        pyplot.ylim(min(bound[1]),max(bound[1]))
        pyplot.xticks([284.,285.5,287.,288.5,290.])

        levels = [1050.,1100.,1150.,1200.,1250.,1300.,1350.]
        cs= pyplot.contour(X, Y, Z2, colors='k', linewidths=0.5) 
        pyplot.clabel(cs,  fontsize=8, fmt='%.0f') 

        pyplot.xlabel(r'$\theta$$_0$ (K)')
        pyplot.ylabel(r'max. SW$_{\mathrm{net}}$ (W m$^{-2}$)')

        savepath = os.path.join(currentdir,'f04b.jpeg')
        fig.savefig(savepath,dpi=1000)


    if (figure_5a == True): #if you want to plot figure 5a

        results_filename = '2DIM_SA_cc_0-30_gammatheta_2-8_SMI=0.2/Results_cc_gammatheta.dat'
        ResDict = open_data(storagedir,[results_filename])

        # defining X, Y, Z variable names, and the results file where these will be plotted from
        xvar = 'ga'
        yvar = 'cc'
        zvar1 = 'WUEeco_ave'
        zvar2 = 'VPD_ave'

        # Rebuilding Y axis to produce the new Y range
        #newaxis = remake_axis('cc', RangeDict['cc'], 'SWnet_12h', results_filename, ResDict)        
        #RangeDict[yvar] = newaxis[1]
	#sys.exit(2)

        # values of the X, Y, Z variables for plotting
        X = list(RangeDict[xvar])
        Y = list(RangeDict[yvar])

        #Z1 = build_Z_grid(xvar, X, yvar, Y, zvar1, results_filename)
        Z1 = numpy.reshape(ResDict[results_filename][zvar1], (-1,101))
        Z1 = Z1
        Z1 = numpy.ma.masked_invalid(Z1)
        #import scipy.ndimage as ndimage
        #Z1 = ndimage.gaussian_filter(Z1, sigma=1.0, order=0)

        #Z2 = build_Z_grid(xvar, X, yvar, Y, zvar2, results_filename)
        Z2 = numpy.reshape(ResDict[results_filename][zvar2], (-1,101))
        Z2 = Z2
        Z2 = numpy.ma.masked_invalid(Z2)
        #import scipy.ndimage as ndimage
        #Z2 = ndimage.gaussian_filter(Z2, sigma=1.0, order=0)

        # box boundaries should be used with pcolormesh(), and box center with contour()
        bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar])
        X, Y = numpy.meshgrid(X,Y)

        # plotting
	
        barlim=numpy.linspace(-4.95,-4.50,200)  # max WUE colorbars: 
        barnorm = colors.BoundaryNorm(barlim,256)
	barlim = [-4.95,-4.90,-4.85,-4.80,-4.75,-4.70,-4.65,-4.60,-4.55,-4.50]
 
        fig = pyplot.figure('WUE1a', figsize=(3.2,2.41))
        fig.subplots_adjust(0.14,0.16,0.89,0.96,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
        pyplot.ticklabel_format(style='sci',axis='x', scilimits=(0,0))
        pyplot.pcolormesh(bound[0], bound[1], Z1)#,norm=barnorm)#, cmap=cm.jet_r)

        cbar = pyplot.colorbar()#ticks=barlim)
        cbar.set_label(r'WUE$_{\mathrm{eco}}$ (g$_{\mathrm{C}}$ kg$_{\mathrm{H2O}}$$^{-1}$')  

        pyplot.xlim(min(bound[0]),max(bound[0]))
        pyplot.ylim(min(bound[1]),max(bound[1]))
    
        levels = [11.5,11.6,11.7,11.8,11.9,12.0,12.1,12.2,12.3,12.4]
        cs= pyplot.contour(X, Y, Z2, colors='k', linewidths=0.5, linestyles='-')
        pyplot.clabel(cs, fontsize=8, fmt='%.1f')
    
        pyplot.xlabel(r'$\gamma_{\theta}$ (K m$^{-1}$)', labelpad=1)
        pyplot.ylabel(r'$cc$ (-)', labelpad=3)#r'max. SW$_{\mathrm{net}}$ (W m$^{-2}$)', labelpad=3)
        pyplot.annotate('(a)', xy=(1.2, -0.15), xycoords='axes fraction', annotation_clip=False, size=10, color='black', va="center", ha="center")

        savepath = os.path.join(currentdir,'f05a.jpeg')
        fig.savefig(savepath,dpi=1000)


    if (figure_5b == True): #if you want to plot figure 5a

        results_filename = '2DIM_SA_cc_0-30_gammatheta_2-8_SMI=0.2/Results_cc_gammatheta.dat'
        ResDict = open_data(storagedir,[results_filename])

        # defining X, Y, Z variable names, and the results file where these will be plotted from
        xvar = 'ga'
        yvar = 'cc'
        zvar1 = 'tm_max'
        zvar2 = 'Rs_ave'

        # Rebuilding Y axis to produce the new Y range
        #newaxis = remake_axis('cc', RangeDict['cc'], 'SWnet_12h', results_filename, ResDict)        
        #RangeDict[yvar] = newaxis[1]

        # values of the X, Y, Z variables for plotting
        X = list(RangeDict[xvar])
        Y = list(RangeDict[yvar])

        #Z1 = build_Z_grid(xvar, X, yvar, Y, zvar1, results_filename)
        Z1 = numpy.reshape(ResDict[results_filename][zvar1], (-1,101))
        Z1 = Z1
        Z1 = numpy.ma.masked_invalid(Z1)
        #import scipy.ndimage as ndimage
        #Z1 = ndimage.gaussian_filter(Z1, sigma=1.0, order=0)

        #Z2 = build_Z_grid(xvar, X, yvar, Y, zvar2, results_filename)
        Z2 = numpy.reshape(ResDict[results_filename][zvar2], (-1,101))
        Z2 = Z2
        Z2 = numpy.ma.masked_invalid(Z2)
        Z2 = 1000./Z2
        #import scipy.ndimage as ndimage
        #Z2 = ndimage.gaussian_filter(Z2, sigma=1.0, order=0)

        # box boundaries should be used with pcolormesh(), and box center with contour()
        bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar])
        X, Y = numpy.meshgrid(X,Y)

        # plotting
	
        barlim=numpy.linspace(295.0,303.0,200)  # max temp colorbars: 
        barnorm = colors.BoundaryNorm(barlim,256)
	barlim=[295.,296.,297.,298.,299.,300.,301.,302.,303.]
 
        fig = pyplot.figure('WUE1b', figsize=(3.2,2.41))
        fig.subplots_adjust(0.14,0.16,0.89,0.96,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
        pyplot.ticklabel_format(style='sci',axis='x', scilimits=(0,0))
        pyplot.pcolormesh(bound[0], bound[1], Z1)#,norm=barnorm, cmap=cm.jet_r)

        cbar = pyplot.colorbar()#ticks=barlim)
        cbar.set_label(r'Maximum ABL temperature (K)')   

        pyplot.xlim(min(bound[0]),max(bound[0]))
        pyplot.ylim(min(bound[1]),max(bound[1]))
    
        levels = [8.3,8.7,9.1]
        cs= pyplot.contour(X, Y, Z2, colors='k', linewidths=0.5, linestyles='-')
        pyplot.clabel(cs, fontsize=8, fmt='%.1f')
    
        pyplot.xlabel(r'$\gamma_{\theta}$ (K m$^{-1}$)', labelpad=1)
        pyplot.ylabel(r'$cc$ (-)', labelpad=3)#r'max. SW$_{\mathrm{net}}$ (W m$^{-2}$)', labelpad=3)
        pyplot.annotate('(b)', xy=(1.2, -0.15), xycoords='axes fraction', annotation_clip=False, size=10, color='black', va="center", ha="center")

        savepath = os.path.join(currentdir,'f05b.jpeg')
        fig.savefig(savepath,dpi=1000)
	
    if (figure_6a == True):

        results_filename = '2DIM_SA_curved_thetam0_284-290_cc_0-30_SMI=0.2/Results_thetam0_cc.dat'
        ResDict = open_data(storagedir,[results_filename])

        # defining X, Y, Z variable names, and the results file where these will be plotted from
        xvar = 'thetam0'
        yvar = 'cc'
        zvar1 = 'NEE_int'
        zvar2 = 'wce_int'

        # Rebuilding Y axis to produce the new Y range
        #newaxis = remake_axis('cc', numpy.arange(0.000,0.303,0.003), 'SWnet_12h', results_filename, ResDict)        
        #RangeDict[yvar] = newaxis[1]

        # values of the X, Y, Z variables for plotting
        X = list(RangeDict[xvar])
        Y = list(RangeDict[yvar])
        X, Y = numpy.meshgrid(X,Y)

        Z1 = numpy.reshape(ResDict[results_filename][zvar1], (-1,101))
        Z1 = Z1.T
        Z1 = numpy.ma.masked_invalid(Z1)
#        import scipy.ndimage as ndimage
#        Z1 = ndimage.gaussian_filter(Z1, sigma=1.0, order=0)

        Z2 = numpy.reshape(ResDict[results_filename][zvar2], (-1,101))
        Z2 = Z2.T
        Z2 = numpy.ma.masked_invalid(Z2)
        Z2 = ndimage.gaussian_filter(Z2, sigma=1.0, order=0)

        bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar])
 
        # plotting

        barlim=numpy.linspace(-21.,-9.,200) # NEE colorbars for low SMI
        barnorm = colors.BoundaryNorm(barlim,256)
        barlim=[-20.,-18.,-16.,-14.,-12.,-10.]  # NEE colorbar ticks

        fig = pyplot.figure('CO2_budget2a', figsize=(3.2,2.41))
        fig.subplots_adjust(0.15,0.16,0.90,0.96,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
        pyplot.pcolormesh(bound[0], bound[1], Z1)#, norm=barnorm, cmap=cm.jet_r)
        cbar = pyplot.colorbar()#ticks=barlim)
        cbar.set_label(r'NEE (g$_{\mathrm{CO2}}$ m$^{-2}$ d$^{-1}$)') 
        pyplot.annotate('(a)', xy=(1.2, -0.15), xycoords='axes fraction', annotation_clip=False, size=10, color='black', va="center", ha="center")

        pyplot.xlim(min(bound[0]),max(bound[0]))
        pyplot.ylim(min(bound[1]),max(bound[1]))
        pyplot.xticks([284.,285.5,287.,288.5,290.])

        levels = [-45.,-46.,-47.,-48.,-49.,-50.,-51.,-52.]
        cs= pyplot.contour(X, Y, Z2, colors='k', linewidths=0.5) 
        pyplot.clabel(cs,  fontsize=8, fmt='%.0f') 

        pyplot.xlabel(r'$\theta$$_0$ (K)')
        pyplot.ylabel(r'max. SW$_{\mathrm{net}}$ (W m$^{-2}$)')

        savepath = os.path.join(currentdir,'f06a.jpeg')
        fig.savefig(savepath,dpi=1000)

    
    if (figure_6b == True):

        results_filename = '2DIM_SA_curved_thetam0_284-290_cc_0-30_SMI=0.2/Results_thetam0_cc.dat'
        ResDict = open_data(storagedir,[results_filename])

        # defining X, Y, Z variable names, and the results file where these will be plotted from
        xvar = 'thetam0'
        yvar = 'cc'
        zvar1 = 'co2_min'
        zvar2 = 'h_max'
        
        # Rebuilding Y axis to produce the new Y range
        #newaxis = remake_axis('cc', numpy.arange(0.000,0.303,0.003), 'SWnet_12h', results_filename, ResDict)        
        #RangeDict[yvar] = newaxis[1]

        # values of the X, Y, Z variables for plotting
        X = list(RangeDict[xvar])
        Y = list(RangeDict[yvar])
        X, Y = numpy.meshgrid(X,Y)

        Z1 = numpy.reshape(ResDict[results_filename][zvar1], (-1,101))
        Z1 = Z1.T
        Z1 = numpy.ma.masked_invalid(Z1)
	Z1 = ndimage.gaussian_filter(Z1, sigma=1.0, order=0)

        Z2 = numpy.reshape(ResDict[results_filename][zvar2], (-1,101))
        Z2 = Z2.T
        Z2 = numpy.ma.masked_invalid(Z2)
        Z2 = ndimage.gaussian_filter(Z2, sigma=1.0, order=0)

        bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar])
 
        # plotting

        barlim=numpy.linspace(344.,349.,200) # cmin colorbars for high SMI
        barnorm = colors.BoundaryNorm(barlim,256)
        barlim=[344.,345.,346.,347.,348.,349.]  # cmin colorbar ticks

        fig = pyplot.figure('CO2_budget2b', figsize=(3.2,2.41))
        fig.subplots_adjust(0.15,0.16,0.90,0.96,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
        pyplot.pcolormesh(bound[0], bound[1], Z1)#, norm=barnorm, cmap=cm.jet_r)
        cbar = pyplot.colorbar()#ticks=barlim)
        cbar.set_label(r'Diurnal CO$_2$ minimum (ppm)') 
        pyplot.annotate('(b)', xy=(1.2, -0.15), xycoords='axes fraction', annotation_clip=False, size=10, color='black', va="center", ha="center")

        pyplot.xlim(min(bound[0]),max(bound[0]))
        pyplot.ylim(min(bound[1]),max(bound[1]))
        pyplot.xticks([284.,285.5,287.,288.5,290.])

        levels = [1050.,1100.,1150.,1200.,1250.,1300.,1350.]
        cs= pyplot.contour(X, Y, Z2, colors='k', linewidths=0.5) 
        pyplot.clabel(cs,  fontsize=8, fmt='%.0f') 

        pyplot.xlabel(r'$\theta$$_0$ (K)')
        pyplot.ylabel(r'max. SW$_{\mathrm{net}}$ (W m$^{-2}$)')

        savepath = os.path.join(currentdir,'f06b.jpeg')
        fig.savefig(savepath,dpi=1000)


    if (figure_8 == True):
    
        print "Printing figure 8!"
        results_filename = '2DIM_SA_cbeta_0-100_SMI_0-1/Results_wg_P4.dat'
        ResDict = open_data(storagedir,[results_filename])

        # defining X, Y, Z variable names, and the results file where these will be plotted from
        xvar = 'smi'
        yvar = 'curvature'
        zvar1 = 'Rs_ave'
        zvar2 = 'wce_h_int'
        zvar3 = 'wcs_h_int'

        # values of the X, Y, Z variables for plotting
        X = list(RangeDict[xvar])
        Y = list(RangeDict[yvar])

        Z1 = numpy.reshape(ResDict[results_filename][zvar1], (-1,101))
        Z1 = 1000./Z1
        Z1 = numpy.ma.masked_invalid(Z1)
        import scipy.ndimage as ndimage
        Z1 = ndimage.gaussian_filter(Z1, sigma=1.0, order=0)

        Z2 = numpy.reshape(ResDict[results_filename][zvar2], (-1,101))
        Z2 = numpy.ma.masked_invalid(Z2)

        Z3 = numpy.reshape(ResDict[results_filename][zvar3], (-1,101))
        Z3 = numpy.ma.masked_invalid(Z3)
	
	Z2 = Z2/Z3
        import scipy.ndimage as ndimage
        Z2 = ndimage.gaussian_filter(Z2, sigma=1.0, order=0)

        # box boundaries should be used with pcolormesh(), and box center with contour()
        bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar])
        X, Y = numpy.meshgrid(X,Y)

        # plotting
	
        barlim=numpy.linspace(2.,17.,200)  # max temp colorbars: 
        barnorm = colors.BoundaryNorm(barlim,256)
	barlim=[2.,5.,8.,11.,14.,17.]
 
        fig = pyplot.figure('Gs', figsize=(3.2,2.41))
        fig.subplots_adjust(0.14,0.16,0.89,0.96,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
        pyplot.ticklabel_format(style='sci',axis='x', scilimits=(0,0))
        pyplot.pcolormesh(bound[0], bound[1], Z1, vmin=2., norm=barnorm)

        cbar = pyplot.colorbar(extend='min',ticks=barlim)
        cbar.set_label(r'$g_s$ (mm s$^{-1}$)')  

        pyplot.xlim(min(bound[0]),max(bound[0]))
        pyplot.ylim(min(bound[1]),max(bound[1]))
    
        levels = [0.5,1.,2.,4.,8.,16.]
        cs= pyplot.contour(X, Y, Z2, levels, colors=['k','w','w','w','w','w','w'], linewidths=[0.5], linestyles=['-'])
        pyplot.clabel(cs, fontsize=8, fmt='%.1f', colors=['k','w','w','w','w','w','w'])
       
        levels = [2.]
        cs= pyplot.contour(X, Y, Z2, levels, colors=['r'], linewidths=[1.5], linestyles=['--'])
        pyplot.clabel(cs, fontsize=0, fmt='%.0f', colors=['r'], manual=[(-100.,-100.)])

        pyplot.scatter(0.2,0.,marker='o',s=8,lw='0.5',edgecolor='black',facecolor='white', zorder=3)
        pyplot.annotate("DL",xy=(0.2,0.), xytext=(0.2,5.), textcoords='data',size=10, weight='bold', color='w', va="center", ha="center")
        pyplot.scatter(0.8,0.,marker='o',s=8,lw='0.5',edgecolor='black',facecolor='black', zorder=3)
        pyplot.annotate("WL",xy=(0.8,0.), xytext=(0.8,5.), textcoords='data',size=10, weight='bold', color='k', va="center", ha="center")
        pyplot.scatter(0.2,100.,marker='o',s=8,lw='0.5',edgecolor='black',facecolor='black', zorder=3)
        pyplot.annotate("DC",xy=(0.2,100.), xytext=(0.2,95.), textcoords='data',size=10, weight='bold', color='k', va="center", ha="center")
    
        pyplot.xlabel(r'SMI (-)', labelpad=1)
        pyplot.ylabel(r'$C_{\beta}$ (%)', labelpad=3)
        pyplot.annotate('(b)', xy=(1.2, -0.15), xycoords='axes fraction', annotation_clip=False, size=10, color='black', va="center", ha="center")

        savepath = os.path.join(currentdir,'f08.jpeg')
        fig.savefig(savepath,dpi=1000)

    if (figure_9a == True):
    
        print "Printing figure 9a!"
        results_filename = '2DIM_SA_day_130-230_SMI_0-1/Results_day_wg.dat'
        ResDict = open_data(storagedir,[results_filename])

        # defining X, Y, Z variable names, and the results file where these will be plotted from
        xvar = 'smi'
        yvar = 'day'
        zvar1 = 'Anet_SUM'
        zvar2 = 'Rs_ave'

        # values of the X, Y, Z variables for plotting
        X = list(RangeDict[xvar])
        Y = list(RangeDict[yvar])

        Z1 = numpy.reshape(ResDict[results_filename][zvar1], (-1,101))
        Z1 = numpy.ma.masked_invalid(Z1)
	Z1 = Z1*12./44.*2.33

        Z2 = numpy.reshape(ResDict[results_filename][zvar2], (-1,101))
        Z2 = numpy.ma.masked_invalid(Z2)
        Z2 = 1000./Z2
        import scipy.ndimage as ndimage
        Z2 = ndimage.gaussian_filter(Z2, sigma=1.0, order=0)

        # box boundaries should be used with pcolormesh(), and box center with contour()
        bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar])
        X, Y = numpy.meshgrid(X,Y)
	
        # plotting
	
        fig = pyplot.figure('DM_yield', figsize=(3.2,2.41))
        fig.subplots_adjust(0.14,0.16,0.89,0.96,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
        pyplot.pcolormesh(bound[0], bound[1], Z1)

        cbar = pyplot.colorbar()
        cbar.set_label(r'DCY (g$_{\mathrm{DM}}$ m$^{-2}$ d$^{-1}$)') 

        pyplot.xlim(min(bound[0]),max(bound[0]))
        pyplot.ylim(min(bound[1]),max(bound[1]))
    
        levels = [0.5,2.0,5.0,10.,15.0]
        cs= pyplot.contour(X, Y, Z2, levels, colors=['w','k','k','k','k','k'], linewidths=[0.5,1.5,0.5,0.5,0.5,0.5], linestyles=['-','--','-','-','-','-'])
        pyplot.clabel(cs, fontsize=8, fmt='%.1f', colors=['w','k','k','k','k','k'])
       
        pyplot.scatter(0.2,216.,marker='o',s=8,lw='0.5',edgecolor='k',facecolor='k', zorder=3)
        pyplot.annotate("DL",xy=(0.2,216.), xytext=(0.2,210.), textcoords='data',size=10, weight='bold', color='k', va="center", ha="center")
        pyplot.scatter(0.8,216.,marker='o',s=8,lw='0.5',edgecolor='k',facecolor='k', zorder=3)
        pyplot.annotate("WL",xy=(0.8,216.), xytext=(0.8,210.), textcoords='data',size=10, weight='bold', color='k', va="center", ha="center")
    
        pyplot.xlabel(r'SMI (-)', labelpad=1)
        pyplot.ylabel(r'time (DOY)', labelpad=3)
        pyplot.annotate('(a)', xy=(1.2, -0.15), xycoords='axes fraction', annotation_clip=False, size=10, color='black', va="center", ha="center")

        savepath = os.path.join(currentdir,'f09a.jpeg')
        fig.savefig(savepath,dpi=1000)

    if (figure_9b == True):
    
        print "Printing figure 9b!"
        results_filename = '2DIM_SA_day_130-230_SMI_0-1/Results_day_wg.dat'
        ResDict = open_data(storagedir,[results_filename])

        # defining X, Y, Z variable names, and the results file where these will be plotted from
        xvar = 'smi'
        yvar = 'day'
        zvar1 = 'DTU_SUM'
        zvar2 = 'tm_max'

        # values of the X, Y, Z variables for plotting
        X = list(RangeDict[xvar])
        Y = list(RangeDict[yvar])

        Z1 = numpy.reshape(ResDict[results_filename][zvar1], (-1,101))
        Z1 = numpy.ma.masked_invalid(Z1)
#        import scipy.ndimage as ndimage
#        Z1 = ndimage.gaussian_filter(Z1, sigma=1.0, order=0)

        Z2 = numpy.reshape(ResDict[results_filename][zvar2], (-1,101))
        Z2 = numpy.ma.masked_invalid(Z2)
        import scipy.ndimage as ndimage
        Z2 = ndimage.gaussian_filter(Z2, sigma=1.0, order=0)

        # box boundaries should be used with pcolormesh(), and box center with contour()
        bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar])
        X, Y = numpy.meshgrid(X,Y)
	
        # plotting
	
        fig = pyplot.figure('RDVR', figsize=(3.2,2.41))
        fig.subplots_adjust(0.14,0.16,0.89,0.96,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
        pyplot.pcolormesh(bound[0], bound[1], Z1)

        cbar = pyplot.colorbar(format='%.2f')
        cbar.set_label('RDVR (-)') #r'max. $h$ (m)') 

        pyplot.xlim(min(bound[0]),max(bound[0]))
        pyplot.ylim(min(bound[1]),max(bound[1]))
    
        levels = [298.,300.,303.,302.,304.,306.,308.]
        cs= pyplot.contour(X, Y, Z2, colors=['k'], linewidths=[0.5], linestyles=['-'])
        pyplot.clabel(cs, fontsize=8, fmt='%.0f', colors=['k'])
       
        pyplot.scatter(0.2,216.,marker='o',s=8,lw='0.5',edgecolor='k',facecolor='k', zorder=3)
        pyplot.annotate("DL",xy=(0.2,216.), xytext=(0.2,210.), textcoords='data',size=10, weight='bold', color='k', va="center", ha="center")
        pyplot.scatter(0.8,216.,marker='o',s=8,lw='0.5',edgecolor='k',facecolor='k', zorder=3)
        pyplot.annotate("WL",xy=(0.8,216.), xytext=(0.8,210.), textcoords='data',size=10, weight='bold', color='k', va="center", ha="center")
    
        pyplot.xlabel(r'SMI (-)', labelpad=1)
        pyplot.ylabel(r'time (DOY)', labelpad=3)
        pyplot.annotate('(b)', xy=(1.2, -0.15), xycoords='axes fraction', annotation_clip=False, size=10, color='black', va="center", ha="center")

        savepath = os.path.join(currentdir,'f09b.jpeg')
        fig.savefig(savepath,dpi=1000)


    if (figure_9c == True):
    
        print "Printing figure 9c!"
        results_filename = '2DIM_SA_day_130-230_SMI_0-1/Results_day_wg.dat'
        ResDict = open_data(storagedir,[results_filename])

        # defining X, Y, Z variable names, and the results file where these will be plotted from
        xvar = 'smi'
        yvar = 'day'
        zvar1 = 'h_max'
        zvar2 = 'SH_int'

        # values of the X, Y, Z variables for plotting
        X = list(RangeDict[xvar])
        Y = list(RangeDict[yvar])

        Z1 = numpy.reshape(ResDict[results_filename][zvar1], (-1,101))
        Z1 = numpy.ma.masked_invalid(Z1)

        Z2 = numpy.reshape(ResDict[results_filename][zvar2], (-1,101))
        Z2 = numpy.ma.masked_invalid(Z2)
        import scipy.ndimage as ndimage
        Z2 = ndimage.gaussian_filter(Z2, sigma=1.0, order=0)

        # box boundaries should be used with pcolormesh(), and box center with contour()
        bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar])
        X, Y = numpy.meshgrid(X,Y)
	
        # plotting
	
        fig = pyplot.figure('max_h', figsize=(3.2,2.41))
        fig.subplots_adjust(0.14,0.16,0.89,0.96,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
        pyplot.pcolormesh(bound[0], bound[1], Z1)

        cbar = pyplot.colorbar()
        cbar.set_label(r'max. $h$ (m)') 

        pyplot.xlim(min(bound[0]),max(bound[0]))
        pyplot.ylim(min(bound[1]),max(bound[1]))
    
        cs= pyplot.contour(X, Y, Z2, colors=['w','k','k','k','k','k'], linewidths=[0.5], linestyles=['-'])
        pyplot.clabel(cs, fontsize=8, fmt='%.2f', colors=['w','k','k','k','k','k'])
       
        pyplot.scatter(0.2,216.,marker='o',s=8,lw='0.5',edgecolor='k',facecolor='k', zorder=3)
        pyplot.annotate("DL",xy=(0.2,216.), xytext=(0.2,210.), textcoords='data',size=10, weight='bold', color='k', va="center", ha="center")
        pyplot.scatter(0.8,216.,marker='o',s=8,lw='0.5',edgecolor='k',facecolor='w', zorder=3)
        pyplot.annotate("WL",xy=(0.8,216.), xytext=(0.8,210.), textcoords='data',size=10, weight='bold', color='w', va="center", ha="center")
    
        pyplot.xlabel(r'SMI (-)', labelpad=1)
        pyplot.ylabel(r'time (DOY)', labelpad=3)
        pyplot.annotate('(c)', xy=(1.2, -0.15), xycoords='axes fraction', annotation_clip=False, size=10, color='black', va="center", ha="center")

        savepath = os.path.join(currentdir,'f09c.jpeg')
        fig.savefig(savepath,dpi=1000)
	

    if (figure_10 == True): #if you want to plot figure 10

        print "Printing figure 10! (EGU poster)"
        results_filename = '2DIM_SA_cbeta_0-100_day_216-236/Results_P4_day.dat'
        ResDict = open_data(storagedir,[results_filename])

        # defining X, Y, Z variable names, and the results file where these will be plotted from
        xvar = 'day2'
        yvar = 'curvature2'
        zvar1 = 'co2_min'
        zvar2 = 'wce_h_int'
	zvar3 = 'wcs_h_int'
	zvar4 = 'w2_0'

        # Rebuilding Y axis to produce the new Y range
        #newaxis = remake_axis('cc', RangeDict['cc'], 'SWnet_12h', results_filename, ResDict)        
        #RangeDict[yvar] = newaxis[1]

        # values of the X, Y, Z variables for plotting
        X = list(RangeDict[xvar])
        Y = list(RangeDict[yvar])

        Z1 = numpy.reshape(ResDict[results_filename][zvar1], (-1,21))
        Z1 = numpy.ma.masked_invalid(Z1)
        import scipy.ndimage as ndimage
        Z1 = ndimage.gaussian_filter(Z1, sigma=1.0, order=0)

        Z2 = numpy.reshape(ResDict[results_filename][zvar2], (-1,21))
        Z2 = numpy.ma.masked_invalid(Z2)

        Z3 = numpy.reshape(ResDict[results_filename][zvar3], (-1,21))
        Z3 = numpy.ma.masked_invalid(Z3)

        Z2 = Z2/Z3
        import scipy.ndimage as ndimage
        Z2 = ndimage.gaussian_filter(Z2, sigma=1.0, order=0)

        Z4 = numpy.reshape(ResDict[results_filename][zvar4], (-1,21))
        Z4 = numpy.ma.masked_invalid(Z4)
	Z4 = (Z4-0.06)/(0.15-0.06)
        import scipy.ndimage as ndimage
        Z4 = ndimage.gaussian_filter(Z4, sigma=1.0, order=0)
	
        # box boundaries should be used with pcolormesh(), and box center with contour()
        bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar])
        X, Y = numpy.meshgrid(X,Y)

        # plotting

        barlim=numpy.linspace(297.,303.,200) # cmin colorbars for high SMI
        barnorm = colors.BoundaryNorm(barlim,256)
        barlim=[297.,298.,299.,300.,301.,302.,303.]  # cmin colorbar ticks

        fig = pyplot.figure('poster_5a', figsize=(3.2,2.41))
        fig.subplots_adjust(0.15,0.16,0.90,0.96,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
        pyplot.pcolormesh(bound[0], bound[1], Z1, cmap=cm.jet_r)#, norm=barnorm)
        cbar = pyplot.colorbar()#ticks=barlim)
        cbar.set_label(r'Diurnal CO$_2$ minimum (ppm)') 
        pyplot.annotate('(a)', xy=(1.2, -0.15), xycoords='axes fraction', annotation_clip=False, size=10, color='black', va="center", ha="center")

        pyplot.xlim(min(bound[0]),max(bound[0]))
        pyplot.ylim(min(bound[1]),max(bound[1]))
        #pyplot.xticks([284.,285.5,287.,288.5,290.])

        levels = [0.5,0.4,0.3,0.2,0.1,0.0] #[0.06,0.07,0.08,0.09,0.1]
        cs= pyplot.contour(X, Y, Z4, levels, colors='k', linewidths=0.5) 
        pyplot.clabel(cs,  fontsize=8, fmt='%.2f') 

        levels = [4.]
        cs= pyplot.contour(X, Y, Z2, levels, colors='k', linewidths=2.5, linestyles='--') 
        pyplot.clabel(cs,  fontsize=0, fmt='%.0f', manual=[(-100.,-100.)]) 

        pyplot.xlabel(r'time (DOY)')
        pyplot.ylabel(r'C$_{\beta}$ (%)')

        savepath = os.path.join(currentdir,'f10_posterEGU.jpeg')
        fig.savefig(savepath,dpi=1000)

    if (figure_11a == True): #if you want to plot figure 10

        print "Printing figure 11a! (EGU poster)"
	results_filename = '2DIM_SA_cbeta_0-100_day_216-236/Results_P4_day.dat'
        ResDict = open_data(storagedir,[results_filename])

        # defining X, Y, Z variable names, and the results file where these will be plotted from
        xvar = 'day2'
        yvar = 'curvature2'
        zvar1 = 'co2_min'
        zvar2 = 'wce_h_int'
	zvar3 = 'wcs_h_int'
	zvar4 = 'w2_0'

       # values of the X, Y, Z variables for plotting
        X = list(RangeDict[xvar])
        Y = list(RangeDict[yvar])

        Z1 = numpy.reshape(ResDict[results_filename][zvar1], (-1,21))
        Z1 = numpy.ma.masked_invalid(Z1)
        import scipy.ndimage as ndimage
        Z1 = ndimage.gaussian_filter(Z1, sigma=1.0, order=0)

        Z2 = numpy.reshape(ResDict[results_filename][zvar2], (-1,21))
        Z2 = numpy.ma.masked_invalid(Z2)

        Z3 = numpy.reshape(ResDict[results_filename][zvar3], (-1,21))
        Z3 = numpy.ma.masked_invalid(Z3)

        Z2 = Z2/Z3
        import scipy.ndimage as ndimage
        Z2 = ndimage.gaussian_filter(Z2, sigma=1.0, order=0)

        Z4 = numpy.reshape(ResDict[results_filename][zvar4], (-1,21))
        Z4 = numpy.ma.masked_invalid(Z4)
	Z4 = (Z4-0.06)/(0.15-0.06)
        import scipy.ndimage as ndimage
        Z4 = ndimage.gaussian_filter(Z4, sigma=1.0, order=0)
	
        # box boundaries should be used with pcolormesh(), and box center with contour()
        bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar])
        X, Y = numpy.meshgrid(X,Y)

        # plotting

        fig = pyplot.figure('poster_11a', figsize=(3.2,2.41))
        fig.subplots_adjust(0.15,0.16,0.90,0.96,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
        pyplot.pcolormesh(bound[0], bound[1], Z1, cmap=cm.jet_r)
        cbar = pyplot.colorbar()
        cbar.set_label(r'Diurnal CO$_2$ minimum (ppm)') 

        pyplot.xlim(min(bound[0]),max(bound[0]))
        pyplot.ylim(min(bound[1]),max(bound[1]))

        levels = [0.5,0.4,0.3,0.2,0.1,0.0] #[0.06,0.07,0.08,0.09,0.1]
        cs= pyplot.contour(X, Y, Z4, levels, colors='k', linewidths=0.5) 
        pyplot.clabel(cs,  fontsize=8, fmt='%.1f') 

        levels = [4.]
        cs= pyplot.contour(X, Y, Z2, levels, colors='k', linewidths=2.5, linestyles='--') 
        pyplot.clabel(cs,  fontsize=0, fmt='%.0f', manual=[(-100.,-100.)]) 

        pyplot.xlabel(r'time (DOY)')
        pyplot.ylabel(r'C$_{\beta}$ (%)')

        savepath = os.path.join(currentdir,'f11a_posterEGU.jpeg')
        fig.savefig(savepath,dpi=1000)

    if (figure_11b == True): #if you want to plot figure 10

        print "Printing figure 11b! (EGU poster)"
	results_filename = '2DIM_SA_cbeta_0-100_day_216-236/Results_P4_day.dat'
        ResDict = open_data(storagedir,[results_filename])

        # defining X, Y, Z variable names, and the results file where these will be plotted from
        xvar = 'day2'
        yvar = 'curvature2'
        zvar1 = 'tm_max'
        zvar2 = 'wce_h_int'
	zvar3 = 'wcs_h_int'
	zvar4 = 'w2_0'

        # Rebuilding Y axis to produce the new Y range
        #newaxis = remake_axis('cc', RangeDict['cc'], 'SWnet_12h', results_filename, ResDict)        
        #RangeDict[yvar] = newaxis[1]

        # values of the X, Y, Z variables for plotting
        X = list(RangeDict[xvar])
        Y = list(RangeDict[yvar])

        Z1 = numpy.reshape(ResDict[results_filename][zvar1], (-1,21))
        Z1 = numpy.ma.masked_invalid(Z1)
        import scipy.ndimage as ndimage
        Z1 = ndimage.gaussian_filter(Z1, sigma=1.0, order=0)

        Z2 = numpy.reshape(ResDict[results_filename][zvar2], (-1,21))
        Z2 = numpy.ma.masked_invalid(Z2)

        Z3 = numpy.reshape(ResDict[results_filename][zvar3], (-1,21))
        Z3 = numpy.ma.masked_invalid(Z3)

        Z2 = Z2/Z3
        import scipy.ndimage as ndimage
        Z2 = ndimage.gaussian_filter(Z2, sigma=1.0, order=0)

        Z4 = numpy.reshape(ResDict[results_filename][zvar4], (-1,21))
        Z4 = numpy.ma.masked_invalid(Z4)
	Z4 = (Z4-0.06)/(0.15-0.06)
        import scipy.ndimage as ndimage
        Z4 = ndimage.gaussian_filter(Z4, sigma=1.0, order=0)
	
        # box boundaries should be used with pcolormesh(), and box center with contour()
        bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar])
        X, Y = numpy.meshgrid(X,Y)

        # plotting

        barlim=numpy.linspace(297.,303.,200) # cmin colorbars for high SMI
        barnorm = colors.BoundaryNorm(barlim,256)
        barlim=[297.,298.,299.,300.,301.,302.,303.]  # cmin colorbar ticks

        fig = pyplot.figure('poster_11b', figsize=(3.2,2.41))
        fig.subplots_adjust(0.15,0.16,0.90,0.96,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
        pyplot.pcolormesh(bound[0], bound[1], Z1)#, norm=barnorm, cmap=cm.jet_r)
        cbar = pyplot.colorbar(ticks=barlim)
        cbar.set_label(r'Diurnal maximum $\theta$ (K)') #r'Diurnal CO$_2$ minimum (ppm)') 

        pyplot.xlim(min(bound[0]),max(bound[0]))
        pyplot.ylim(min(bound[1]),max(bound[1]))

        levels = [0.5,0.4,0.3,0.2,0.1,0.0] #[0.06,0.07,0.08,0.09,0.1]
        cs= pyplot.contour(X, Y, Z4, levels, colors='k', linewidths=0.5) 
        pyplot.clabel(cs,  fontsize=8, fmt='%.1f') 

        levels = [4.]
        cs= pyplot.contour(X, Y, Z2, levels, colors='k', linewidths=2.5, linestyles='--') 
        pyplot.clabel(cs,  fontsize=0, fmt='%.0f', manual=[(-100.,-100.)]) 

        pyplot.xlabel(r'time (DOY)')
        pyplot.ylabel(r'C$_{\beta}$ (%)')

        savepath = os.path.join(currentdir,'f11b_posterEGU.jpeg')
        fig.savefig(savepath,dpi=1000)
		
    if (figure_12 == True):
    
        print "Printing figure 12 (EGU poster)!"
        results_filename = '2DIM_SA_cbeta_0-100_SMI_0-1/Results_wg_P4.dat'
        ResDict = open_data(storagedir,[results_filename])

        # defining X, Y, Z variable names, and the results file where these will be plotted from
        xvar = 'smi'
        yvar = 'curvature'
        zvar1 = 'Rs_ave'
        zvar2 = 'wce_h_int'
        zvar3 = 'wcs_h_int'

        # values of the X, Y, Z variables for plotting
        X = list(RangeDict[xvar])
        Y = list(RangeDict[yvar])

        Z1 = numpy.reshape(ResDict[results_filename][zvar1], (-1,101))
        Z1 = 1000./Z1
        Z1 = numpy.ma.masked_invalid(Z1)
        import scipy.ndimage as ndimage
        Z1 = ndimage.gaussian_filter(Z1, sigma=1.0, order=0)

        Z2 = numpy.reshape(ResDict[results_filename][zvar2], (-1,101))
        Z2 = numpy.ma.masked_invalid(Z2)

        Z3 = numpy.reshape(ResDict[results_filename][zvar3], (-1,101))
        Z3 = numpy.ma.masked_invalid(Z3)
	
	Z2 = Z2/(Z2+Z3)
        import scipy.ndimage as ndimage
        Z2 = ndimage.gaussian_filter(Z2, sigma=1.0, order=0)

        # box boundaries should be used with pcolormesh(), and box center with contour()
        bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar])
        X, Y = numpy.meshgrid(X,Y)

        # plotting
	
        barlim=numpy.linspace(2.,17.,200)  # max temp colorbars: 
        barnorm = colors.BoundaryNorm(barlim,256)
	barlim=[2.,5.,8.,11.,14.,17.]
 
        fig = pyplot.figure('Gs', figsize=(3.2,2.41))
        fig.subplots_adjust(0.14,0.16,0.89,0.96,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
        pyplot.ticklabel_format(style='sci',axis='x', scilimits=(0,0))
        pyplot.pcolormesh(bound[0], bound[1], Z1, vmin=2., norm=barnorm)

        cbar = pyplot.colorbar(extend='min',ticks=barlim)
        cbar.set_label(r'$g_s$ (mm s$^{-1}$)')  

        pyplot.xlim(min(bound[0]),max(bound[0]))
        pyplot.ylim(min(bound[1]),max(bound[1]))
    
        levels = [0.5,1.,2.,4.,8.,16.]
        cs= pyplot.contour(X, Y, Z2, colors=['k','k','w','w','w','w','w'], linewidths=[0.5], linestyles=['-'])
        pyplot.clabel(cs, fontsize=8, fmt='%.1f', colors=['k','k','w','w','w','w','w'])
       
        levels = [0.8]
        cs= pyplot.contour(X, Y, Z2, levels, colors=['r'], linewidths=[1.5], linestyles=['--'])
        pyplot.clabel(cs, fontsize=0, fmt='%.0f', colors=['r'], manual=[(-100.,-100.)])

        pyplot.scatter(0.8,0.,marker='o',s=8,lw='0.5',edgecolor='black',facecolor='black', zorder=3)
        pyplot.annotate("Wet S.",xy=(0.8,0.), xytext=(0.8,5.), textcoords='data',size=10, weight='bold', color='k', va="center", ha="center")
        pyplot.scatter(0.2,100.,marker='o',s=8,lw='0.5',edgecolor='black',facecolor='black', zorder=3)
        pyplot.annotate("Dry I.",xy=(0.23,100.), xytext=(0.23,95.), textcoords='data',size=10, weight='bold', color='k', va="center", ha="center")
    
        pyplot.xlabel(r'SMI (-)', labelpad=1)
        pyplot.ylabel(r'$C_{\beta}$ (%)', labelpad=3)
        pyplot.annotate('(b)', xy=(1.2, -0.15), xycoords='axes fraction', annotation_clip=False, size=10, color='black', va="center", ha="center")

        savepath = os.path.join(currentdir,'f12_posterEGU.jpeg')
        fig.savefig(savepath,dpi=1000)    

    if (figure_13a == True):
    
        print "Printing figure 13a (EGU poster)!"
        results_filename = '2DIM_SA_day_130-230_SMI_0-1/Results_day_wg.dat'
        ResDict = open_data(storagedir,[results_filename])

        # defining X, Y, Z variable names, and the results file where these will be plotted from
        xvar = 'day'
        yvar = 'smi'
        zvar1 = 'Anet_SUM'
        zvar2 = 'Rs_ave'

        # values of the X, Y, Z variables for plotting
        X = list(RangeDict[xvar])
        Y = list(RangeDict[yvar])

        Z1 = numpy.reshape(ResDict[results_filename][zvar1], (-1,101))
        Z1 = numpy.ma.masked_invalid(Z1)
	Z1 = Z1.T
	Z1 = Z1*12./44.*2.33

        Z2 = numpy.reshape(ResDict[results_filename][zvar2], (-1,101))
        Z2 = numpy.ma.masked_invalid(Z2)
	Z2 = Z2.T
        Z2 = 1000./Z2
        import scipy.ndimage as ndimage
        Z2 = ndimage.gaussian_filter(Z2, sigma=1.0, order=0)

        # box boundaries should be used with pcolormesh(), and box center with contour()
        bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar])
        X, Y = numpy.meshgrid(X,Y)
	
        # plotting
	
        fig = pyplot.figure('DM_yield', figsize=(3.2,2.41))
        fig.subplots_adjust(0.14,0.16,0.89,0.96,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
        pyplot.pcolormesh(bound[0], bound[1], Z1)

        cbar = pyplot.colorbar()
        cbar.set_label(r'DCY (g$_{\mathrm{DM}}$ m$^{-2}$ d$^{-1}$)') 

        pyplot.xlim(min(bound[0]),max(bound[0]))
        pyplot.ylim(min(bound[1]),max(bound[1]))
    
        levels = [0.5,2.0,5.0,10.,15.0]
        cs= pyplot.contour(X, Y, Z2, levels, colors=['w','k','k','k','k','k'], linewidths=[0.5,1.5,0.5,0.5,0.5,0.5], linestyles=['-','--','-','-','-','-'])
        pyplot.clabel(cs, fontsize=8, fmt='%.1f', colors=['w','k','k','k','k','k'])
       
        pyplot.ylabel(r'SMI (-)', labelpad=1)
        pyplot.xlabel(r'time (DOY)', labelpad=3)

        savepath = os.path.join(currentdir,'f13a_posterEGU.jpeg')
        fig.savefig(savepath,dpi=1000)

    if (figure_13b == True):
    
        print "Printing figure 13b (EGU poster)!"
        results_filename = '2DIM_SA_day_130-230_SMI_0-1/Results_day_wg.dat'
        ResDict = open_data(storagedir,[results_filename])

        # defining X, Y, Z variable names, and the results file where these will be plotted from
        xvar = 'day'
        yvar = 'smi'
        zvar1 = 'DTU_SUM'
        zvar2 = 'tm_max'

        # values of the X, Y, Z variables for plotting
        X = list(RangeDict[xvar])
        Y = list(RangeDict[yvar])

        Z1 = numpy.reshape(ResDict[results_filename][zvar1], (-1,101))
        Z1 = numpy.ma.masked_invalid(Z1)
	Z1 = Z1.T
#        import scipy.ndimage as ndimage
#        Z1 = ndimage.gaussian_filter(Z1, sigma=1.0, order=0)

        Z2 = numpy.reshape(ResDict[results_filename][zvar2], (-1,101))
        Z2 = numpy.ma.masked_invalid(Z2)
	Z2 = Z2.T
        import scipy.ndimage as ndimage
        Z2 = ndimage.gaussian_filter(Z2, sigma=1.0, order=0)

        # box boundaries should be used with pcolormesh(), and box center with contour()
        bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar])
        X, Y = numpy.meshgrid(X,Y)
	
        # plotting
	
        fig = pyplot.figure('RDVR', figsize=(3.2,2.41))
        fig.subplots_adjust(0.14,0.16,0.89,0.96,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
        pyplot.pcolormesh(bound[0], bound[1], Z1)

        cbar = pyplot.colorbar(format='%.2f')
        cbar.set_label('RDVR (-)') #r'max. $h$ (m)') 

        pyplot.xlim(min(bound[0]),max(bound[0]))
        pyplot.ylim(min(bound[1]),max(bound[1]))
    
        levels = [298.5,300.,301.5,303.,304.5,306.,307.5]
        cs= pyplot.contour(X, Y, Z2, levels, colors=['k'], linewidths=[0.5], linestyles=['-'])
        pyplot.clabel(cs, fontsize=8, fmt='%.1f', colors=['k'])
       
        pyplot.ylabel(r'SMI (-)', labelpad=1)
        pyplot.xlabel(r'time (DOY)', labelpad=3)

        savepath = os.path.join(currentdir,'f13b_posterEGU.jpeg')
        fig.savefig(savepath,dpi=1000)


    if (figure_13c == True):
    
        print "Printing figure 13c (EGU poster)!"
        results_filename = '2DIM_SA_day_130-230_SMI_0-1/Results_day_wg.dat'
        ResDict = open_data(storagedir,[results_filename])

        # defining X, Y, Z variable names, and the results file where these will be plotted from
        xvar = 'day'
        yvar = 'smi'
        zvar1 = 'h_max'
        zvar2 = 'SH_int'

        # values of the X, Y, Z variables for plotting
        X = list(RangeDict[xvar])
        Y = list(RangeDict[yvar])

        Z1 = numpy.reshape(ResDict[results_filename][zvar1], (-1,101))
        Z1 = numpy.ma.masked_invalid(Z1)
	Z1 = Z1.T

        Z2 = numpy.reshape(ResDict[results_filename][zvar2], (-1,101))
        Z2 = numpy.ma.masked_invalid(Z2)
	Z2 = Z2.T
        import scipy.ndimage as ndimage
        Z2 = ndimage.gaussian_filter(Z2, sigma=1.0, order=0)

        # box boundaries should be used with pcolormesh(), and box center with contour()
        bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar])
        X, Y = numpy.meshgrid(X,Y)
	
        # plotting
	
        fig = pyplot.figure('max_h', figsize=(3.2,2.41))
        fig.subplots_adjust(0.14,0.16,0.89,0.96,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
        pyplot.pcolormesh(bound[0], bound[1], Z1)

        cbar = pyplot.colorbar()
        cbar.set_label(r'max. $h$ (m)') 

        pyplot.xlim(min(bound[0]),max(bound[0]))
        pyplot.ylim(min(bound[1]),max(bound[1]))
    
        cs= pyplot.contour(X, Y, Z2, colors=['w','k','k','k','k','k'], linewidths=[0.5], linestyles=['-'])
        pyplot.clabel(cs, fontsize=8, fmt='%.0f', colors=['w','k','k','k','k','k'])
       
        pyplot.ylabel(r'SMI (-)', labelpad=1)
        pyplot.xlabel(r'time (DOY)', labelpad=3)

        savepath = os.path.join(currentdir,'f13c_posterEGU.jpeg')
        fig.savefig(savepath,dpi=1000)
    
    if (figure_14a == True):
    
        print "Printing figure 14a (EGU poster)!"
        results_filename = '2DIM_SA_div_0-4x10-5_SMI_0.4-0.6/Results_wsls_wg2.dat'
        ResDict = open_data(storagedir,[results_filename])

        # defining X, Y, Z variable names, and the results file where these will be plotted from
        xvar = 'smi2'
        yvar = 'wsls'
        zvar1 = 'co2_min'

        # values of the X, Y, Z variables for plotting
        X = list(RangeDict[xvar])
        Y = list(RangeDict[yvar])

        Z1 = numpy.reshape(ResDict[results_filename][zvar1], (-1,101))
        Z1 = numpy.ma.masked_invalid(Z1)
	Z1 = Z1.T

        # box boundaries should be used with pcolormesh(), and box center with contour()
        bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar])
        X, Y = numpy.meshgrid(X,Y)

        # plotting

        barlim=numpy.linspace(341.,361.,200) # cmin colorbars for high SMI
        barnorm = colors.BoundaryNorm(barlim,256)
        barlim=[342.,345.,348.,351.,354.,357.,360.]  # cmin colorbar ticks

        fig = pyplot.figure('poster_14a', figsize=(3.2,2.41))
        fig.subplots_adjust(0.15,0.14,0.90,0.94,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
        pyplot.pcolormesh(bound[0], bound[1], Z1, cmap=cm.jet_r, norm=barnorm)
        pyplot.ticklabel_format(style='sci',axis='y', scilimits=(0,0))
        cbar = pyplot.colorbar(ticks=barlim)
        cbar.set_label(r'Diurnal CO$_2$ minimum (ppm)') 

        pyplot.xlim(min(bound[0]),max(bound[0]))
        pyplot.ylim(min(bound[1]),max(bound[1]))

        pyplot.xlabel(r'SMI (-)', labelpad=1)
        pyplot.ylabel(r'$D$ (s$^{-1}$)', labelpad=1)

        savepath = os.path.join(currentdir,'f14a_posterEGU.jpeg')
        fig.savefig(savepath,dpi=1000)

    if (figure_14b == True):
    
        print "Printing figure 14b (EGU poster)!"
        results_filename = '../2DIM_SA_gammatheta_2-8_deltatheta_0.5-5/Results_gammatheta_deltatheta.dat'
        ResDict = open_data(currentdir,[results_filename])

        # defining X, Y, Z variable names, and the results file where these will be plotted from
        xvar = 'gammatheta'
        yvar = 'deltatheta'
        zvar1 = 'co2_min'

        # values of the X, Y, Z variables for plotting
        X = list(RangeDict[xvar])
        Y = list(RangeDict[yvar])

        Z1 = numpy.reshape(ResDict[results_filename][zvar1], (-1,101))
        Z1 = numpy.ma.masked_invalid(Z1)
	Z1 = Z1

        # box boundaries should be used with pcolormesh(), and box center with contour()
        bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar])
        X, Y = numpy.meshgrid(X,Y)

        # plotting

        barlim=numpy.linspace(347.,355.,200) # cmin colorbars for high SMI
        barnorm = colors.BoundaryNorm(barlim,256)
        barlim=[347.,348.,349.,350.,351.,352.,353.,354.,355.]  # cmin colorbar ticks

        fig = pyplot.figure('poster_14b', figsize=(3.2,2.41))
        fig.subplots_adjust(0.15,0.16,0.90,0.96,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
        pyplot.pcolormesh(bound[0], bound[1], Z1, cmap=cm.jet_r, norm=barnorm)
        pyplot.ticklabel_format(style='sci',axis='x', scilimits=(0,0))
        cbar = pyplot.colorbar(ticks=barlim)
        cbar.set_label(r'Diurnal CO$_2$ minimum (ppm)') 

        pyplot.xlim(min(bound[0]),max(bound[0]))
        pyplot.ylim(min(bound[1]),max(bound[1]))

        pyplot.xlabel(r'$\gamma_{\theta}$ (K m$^{-1}$)', labelpad=1)
        pyplot.ylabel(r'$\Delta{\theta}_0$ (K)', labelpad=1)

        savepath = os.path.join(currentdir,'f14b_posterEGU.jpeg')
        fig.savefig(savepath,dpi=1000)


