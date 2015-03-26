#!/usr/bin/env python
# sensitivity.py


"""
Started 10 oct 2013
Author: Marie Combe

Plotting script for the sensitivity analysis results - Article 1

"""

# importing the necessary modules
import os
import sys
import getopt
import numpy		# for arrays handling
from matplotlib import gridspec
from scipy.stats.mstats import mquantiles # for calculating the 1st, 2nd (median), and 3rd quartiles of a data set
from matplotlib import pyplot
from matplotlib import rc
from matplotlib import rcdefaults
from matplotlib import colors
from matplotlib import cm
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA


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
def remake_axis(var, rangevar, newvar, resultsfile):

    W = [0.]*len(rangevar)

    for i,item in enumerate(rangevar):
        V = []
        for j,newitem in enumerate(ResDict[resultsfile][newvar]):
            if (abs(ResDict[resultsfile][var][j]-item) < 1e-6):
                V = V + [newitem]
        # print item, V  # uncomment this line if you want to check the result of the function
        W[i] = numpy.mean(numpy.array(V))

    ResDict[resultsfile]['grouped_%s'%newvar] = [0.]*len(ResDict[resultsfile][var])

    for i,item in enumerate(rangevar):
        for j,newitem in enumerate(ResDict[resultsfile][newvar]):
            if (abs(ResDict[resultsfile][var][j]-item) < 1e-6):
                ResDict[resultsfile]['grouped_%s'%newvar][j] = W[i]

    return ResDict[resultsfile]['grouped_%s'%newvar], W

#----------------------------------------------------------------------
def build_Z_grid(xvar, Xrange, yvar, Yrange, zvar, resultsfile):

    Z=[[0.]*len(Xrange)]*len(Yrange)
    Z=numpy.array(Z,float)
    for i,valx in enumerate(Xrange):  # for each x of Range1
        for j, valy in enumerate(Yrange): # for each y of Range2
            index=-1
            for h in range(0,len(ResDict[resultsfile][zvar])): # we match the z value
                if (((abs(ResDict[resultsfile][xvar][h]-valx))<1e-6) and ((abs(ResDict[resultsfile][yvar][h]-valy))<1e-6)): 
                    Z[j,i]= ResDict[resultsfile][zvar][h] # we assign a z value in the 2D array Z

    return Z

 


#######################################################################
########################### START OF SCRIPT ###########################

if __name__ == "__main__":

    try:                                
        opts, args = getopt.getopt(sys.argv[1:], "-h")
    except getopt.GetoptError:           
        print "Error"
        sys.exit(2)

    for options in opts:
        options=options[0].lower()
        if options == '-h':
            helptext = """
This script produces a sensitivity analysis plot of the MXL-AGs model.
This is a 3-variables plot, presented in 2D (2 variables on x- and y-axis, and one 
result variable in colors).

You should specify not specify any arguments to run this script:
./Article1figures_part2.py
                
                """
            
            print helptext
            
            sys.exit(2)

    # Working directories
    currentdir = os.getcwd()
    storagedir = '/Storage/MXL/Marie/MXL-Ags_Sensitivity'

    print 'Opening results files...\n'

    # Open results
    ResDict = open_data(storagedir,[
                         '2DIM_SA_thetam0_284-290_cc_0-30_SMI=0.8/Results_thetam0_cc.dat',
                         '2DIM_SA_thetam0_284-290_cc_0-30_SMI=0.2/Results_thetam0_cc.dat',
			 '2DIM_SA_cc_0-30_gammatheta_2-8_SMI=0.8/Results_cc_ga.dat',
     			 '2DIM_SA_cc_0-30_gammatheta_2-8_SMI=0.2/Results_cc_ga.dat',
			 '2DIM_SA_cbeta_0-100_SMI_0-1/Results_wg_P4.dat',
			 '2DIM_SA_day_130-230_SMI_0-1/Results_day_wg.dat'
			])
    
    RangeDict = {'thetam0': numpy.arange(284.,290.01,0.06), 
                 'wg': numpy.arange(0.060,0.363,0.003),
                 'cc': numpy.arange(0.000,0.303,0.003), 
                 'ga': numpy.arange(0.002,0.00806,0.00006),
                 'smi': numpy.arange(0.0,1.001,0.01), #numpy.arange(0.4,0.601,0.00667), #numpy.arange(0.4,0.601,0.002) 
		 'curvature': numpy.arange(0.,100.1,1.),
		 'wsls': numpy.arange(0.,0.0000401,0.0000004),
		 'day': numpy.arange(130.,230.5,1.) #numpy.arange(200.,230.5,1.)
		 #'gammatheta': numpy.arange(0.002,0.00806,0.00006)
		 } ;


    # First close all opened figures in python
    pyplot.close('all')

    # general plot properties
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
##
    rc('contour', negative_linestyle = 'solid')
    rc('lines', linewidth=1)
    rc('axes', labelsize=8, linewidth=0.5)
    rc('legend', frameon=True, handlelength=3., fontsize=8)
    rc('grid', linewidth=0.)


    ################################################################
    # figure 1: CMIN with x=thetam0, y=cc(SWin), and contour plot is BL_h

    # defining X, Y, Z variable names, and the results file where these will be plotted from
    xvar = 'thetam0'
    yvar = 'cc'
    zvar1 = 'co2_min'
    zvar2 = 'h_max'
    zvar3 = 'tm_max'

    results_filename = '2DIM_SA_thetam0_284-290_cc_0-30_SMI=0.2/Results_thetam0_cc.dat'

    # Rebuilding Y axis to produce the new Y range
    newaxis = remake_axis('cc', RangeDict['cc'], 'SWnet_12h', results_filename)        
    RangeDict[yvar] = newaxis[1]


    print 'Plotting %s as colors, and %s as contour lines...' %(zvar1, zvar2)


    # values of the X, Y, Z variables for plotting
    X = list(RangeDict[xvar])
    Y = list(RangeDict[yvar])
    X, Y = numpy.meshgrid(X,Y)

    #Z1 = build_Z_grid(xvar, X, yvar, Y, zvar1, results_filename)
    Z1 = numpy.reshape(ResDict[results_filename][zvar1], (-1,101))
    Z1 = Z1.T
    Z1 = numpy.ma.masked_invalid(Z1)
    import scipy.ndimage as ndimage
    Z1 = ndimage.gaussian_filter(Z1, sigma=1.0, order=0)
    #Z1 = 1000./Z1

    #Z2 = build_Z_grid(xvar, X, yvar, Y, zvar2, results_filename)
    Z2 = numpy.reshape(ResDict[results_filename][zvar2], (-1,101))
    Z2 = Z2.T
    Z2 = numpy.ma.masked_invalid(Z2)
    
    # use a Gaussian filter to smooth jagged contours: (entrainment flux only)
    import scipy.ndimage as ndimage
    Z2 = ndimage.gaussian_filter(Z2, sigma=1.0, order=0)
    #Z2 = 1000./Z2


#    Z3 = numpy.reshape(ResDict[results_filename][zvar3], (-1,101))
#    Z3 = Z3.T
#    Z3 = numpy.ma.masked_invalid(Z3)
    #Z2 = 100.*Z3/(Z3+Z2)
    

    # limitation of the parameter space
#    limth1 = 33 # theta_0 = 283.96
#    limth2 = 85 # theta_0 = 290.08
#    X=X[:,limth1:limth2]
#    Y=Y[:,limth1:limth2]
#    Z1 = Z1[:,limth1:limth2]
#    Z2 = Z2[:,limth1:limth2]
#    Z3 = Z3[:,limth1:limth2]
    
    # box boundaries should be used with pcolormesh(), and box center with contour()
    #bound = make_pcolormesh_box_boundaries(RangeDict[xvar][limth1:limth2],RangeDict[yvar])
	
    # box boundaries should be used with pcolormesh(), and box center with contour()
    bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar])
 
    # plotting

# CO2 colorbars:
#    barlim=numpy.linspace(344.,349.,200) # high SMI
    barlim=numpy.linspace(363.,368.,200) # low SMI
# NEE colorbars:
#    barlim=numpy.linspace(-71.,-59.,200)
#    barlim=numpy.linspace(-21.,-9.,200)
    import matplotlib.colors as colors
    barnorm = colors.BoundaryNorm(barlim,256)

    fig = pyplot.figure('CMIN1', figsize=(3.2,2.41))
    fig.subplots_adjust(0.15,0.16,0.90,0.96,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
    pyplot.pcolormesh(bound[0], bound[1], Z1, norm=barnorm)#, cmap=cm.jet_r)

# CO2 colorbar ticks:
#    barlim=[344.,345.,346.,347.,348.,349.]
    barlim=[363.,364.,365.,366.,367.,368.]
# NEE colorbar ticks:
#    barlim=[-70.,-68.,-66.,-64.,-62.,-60.]#[-66.,-64.,-62.,-60.,-58.,-56.]
#    barlim=[-20.,-18.,-16.,-14.,-12.,-10.]

    cbar = pyplot.colorbar(ticks=barlim)
    cbar.set_label(r'Diurnal CO$_2$ minimum (ppm)') #r'NEE (g$_{\mathrm{CO2}}$ m$^{-2}$ d$^{-1}$)') 
    pyplot.annotate('(b)', xy=(1.2, -0.15), xycoords='axes fraction', annotation_clip=False, size=10, color='black', va="center", ha="center")

    pyplot.xlim(min(bound[0]),max(bound[0]))
    pyplot.ylim(min(bound[1]),max(bound[1]))
    #pyplot.xticks([280.,283.,286.,289.,292.])
    pyplot.xticks([284.,285.5,287.,288.5,290.])

#    levels = [-11.0,-10.5,-10.0,-9.5,-9.0]  
#    levels = [-45.,-46.,-47.,-48.,-49.,-50.,-51.,-52.]
#    levels = [-21.,-20.,-19.,-18.,-17.,-16.,-15.,-14.,-13.,-12.,-11.,-10.]
    levels = [1500.,1550.,1600.,1650.,1700.]
#    levels = [1050.,1100.,1150.,1200.,1250.,1300.,1350.]
    cs= pyplot.contour(X, Y, Z2, levels, colors='k', linewidths=0.5) #levels,
    pyplot.clabel(cs,  fontsize=8, fmt='%.0f') #levels[1::2],

#    #levels2 = [9.,9.5,10.,10.5,11.]    
#    cs2= pyplot.contour(X, Y, Z3, colors='k', linestyles='--', linewidths=0.5)
#    pyplot.clabel(cs2,  fontsize=8, fmt='%.0f') #levels[1::2],

    #levels3 = [900.,1000.,1100.,1200.,1300.,1400.]    
    #cs3= pyplot.contour(X, Y, Z4, colors='k', linestyles=':', linewidths=0.5)
    #pyplot.clabel(cs3,  fontsize=8, fmt='%.2f') #levels[1::2],
    
#    pyplot.scatter(286.,783.87,marker='s', s=8,lw='0.5',edgecolor='black',facecolor='white',zorder=4)
    pyplot.xlabel(r'$\theta$$_0$ (K)')
    pyplot.ylabel(r'max. SW$_{\mathrm{net}}$ (W m$^{-2}$)')


    savepath = os.path.join('/Users/mariecombe/Modeling/MXL/MLmodel/CMIN1.jpeg')
    fig.savefig(savepath,dpi=1000)





    ################################################################
    # figure 2: WUE with x=wg0, y=gammaq, and contour plot is BL_h

#    # defining X, Y, Z variable names, and the results file where these will be plotted from
#    xvar = 'ga'
#    yvar = 'cc'
#    zvar1 = 'tm_max'
#    zvar2 = 'Rs_ave'
#    zvar3 = 'NEE_int'
#    results_filename = '2DIM_SA_cc_0-30_gammatheta_2-8_SMI=0.8/Results_cc_ga.dat'#'Results_new_smi0.8_cc_ga.dat' #'Results_beta=0.20_cc_gammatheta.dat'
#
#    # Rebuilding Y axis to produce the new Y range
#    newaxis = remake_axis('cc', RangeDict['cc'], 'SWnet_12h', results_filename)        
#    RangeDict[yvar] = newaxis[1]
#
#
#    print 'Plotting %s as colors, and %s as contour lines...' %(zvar1, zvar2)
#
#
#    # values of the X, Y, Z variables for plotting
#    X = list(RangeDict[xvar])
#    Y = list(RangeDict[yvar])
#
#    #Z1 = build_Z_grid(xvar, X, yvar, Y, zvar1, results_filename)
#    Z1 = numpy.reshape(ResDict[results_filename][zvar1], (-1,101))
#    Z1 = Z1
#    Z1 = numpy.ma.masked_invalid(Z1)
#    #import scipy.ndimage as ndimage
#    #Z1 = ndimage.gaussian_filter(Z1, sigma=1.0, order=0)
#
#    #Z2 = build_Z_grid(xvar, X, yvar, Y, zvar2, results_filename)
#    Z2 = numpy.reshape(ResDict[results_filename][zvar2], (-1,101))
#    Z2 = Z2
#    Z2 = numpy.ma.masked_invalid(Z2)
#    Z2 = 1000./Z2
#    #import scipy.ndimage as ndimage
#    #Z2 = ndimage.gaussian_filter(Z2, sigma=1.0, order=0)
#
#    #Z1 = (1000./(Z1+1.6*Z2))/(1000./(Z1+Z2))
#    #Z1 = (Z1/numpy.min(Z1)-1.)*100.
#    #Z2 = (Z2/numpy.min(Z2)-1.)*100.
#    #Z1 = ((Z5 / 12.691) - 1. )*100.
#    #baseH = 3.83364/2.5
#    #baseH = 10.37169/2.5
#    #Z2 = Z2 / 2.5 - baseH # to go from MJ.m-2.day-1 to kgH2O.m-2.day-1 (or mm.day-1), divide LE_int by 2.5
#    
#    #Z3 = build_Z_grid(xvar, X, yvar, Y, zvar3, results_filename)
#    #Z3 = numpy.reshape(ResDict[results_filename][zvar3], (-1,101))
#    #Z3 = Z3
#    #Z3 = numpy.ma.masked_invalid(Z3)
#    #baseC = 14.61138*12./44.
#    #baseC = 63.26387*12./44.
#    #Z3 = -12./44.*Z3 - baseC
#    
#    #Z4 = Z3/Z2 # WUEeco in gC.kgH2O-1 per day
#    #Z4 = numpy.ma.masked_greater(Z4, 10.)
#    #Z4 = numpy.ma.masked_less(Z4, -10.)
#        
#    # box boundaries should be used with pcolormesh(), and box center with contour()
#    bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar])
#    X, Y = numpy.meshgrid(X,Y)
#
## WUE colorbars:
##    barlim=numpy.linspace(-3.20,-2.75,200)
##    barlim=numpy.linspace(-4.95,-4.50,200)
## max temp colorbars:    
#    barlim=numpy.linspace(295.0,303.0,200)   
## gs colorbars:
##    barlim=numpy.linspace(11.5,12.4,200)
##    barlim=numpy.linspace(1.8,2.7,200)
#    import matplotlib.colors as colors
#    barnorm = colors.BoundaryNorm(barlim,256)
#    # plotting
#    fig = pyplot.figure('WUE1', figsize=(3.2,2.41))
#    fig.subplots_adjust(0.14,0.16,0.89,0.96,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
#    pyplot.ticklabel_format(style='sci',axis='x', scilimits=(0,0))
#    pyplot.pcolormesh(bound[0], bound[1], Z1,norm=barnorm)#, cmap=cm.jet_r)
#
## WUE colorbars:
##    barlim=[2.75,2.80,2.85,2.90,2.95,3.00,3.05,3.10,3.15,3.20]
##    barlim = [-3.20,-3.15,-3.10,-3.05,-3.00,-2.95,-2.90,-2.85,-2.80,-2.75]
##    barlim = [-4.95,-4.90,-4.85,-4.80,-4.75,-4.70,-4.65,-4.60,-4.55,-4.50]
## max temp colorbars:
#    barlim=[295.,296.,297.,298.,299.,300.,301.,302.,303.]
## gs colorbars:
##    barlim=[11.5,11.6,11.7,11.8,11.9,12.0,12.1,12.2,12.3,12.4]
##    barlim=[1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7]
#    cbar = pyplot.colorbar(ticks=barlim)
#    cbar.set_label(r'Maximum ABL temperature (K)') #r'WUE$_{\mathrm{eco}}$ (g$_{\mathrm{C}}$ kg$_{\mathrm{H2O}}$$^{-1}$')  
#
#    pyplot.xlim(min(bound[0]),max(bound[0]))
#    pyplot.ylim(min(bound[1]),max(bound[1]))
#    
#    #levels = [71.5,72.0,72.5,73.0,73.5]
#    #levels= [900.,1200.,1500.,1800.,2100.,2400.,2700.,3000.,3300.,3600.]
#    #levels= [3.50,3.75,4.00,4.25]
#    levels = [11.5,11.6,11.7,11.8,11.9,12.0,12.1,12.2,12.3,12.4]
##    levels = [8.3,8.7,9.1]
##    levels = [13.0,13.4,13.8,14.2,14.6,15.0,15.4,15.8,16.2,16.6,17.0]
##    levels = [1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7]
##    levels = [0.,1.,2.,3.]
#    cs= pyplot.contour(X, Y, Z2, levels, colors='k', linewidths=0.5, linestyles='-')
#    pyplot.clabel(cs, fontsize=8, fmt='%.1f')
#    
##    #levels = [-2.00,-1.5,-1.0,-0.5]
##    cs2= pyplot.contour(X, Y, Z3, linestyles='--',colors='k', linewidths=0.5)
##    pyplot.clabel(cs2, fontsize=8, fmt='%.2f')
#    
#    #pyplot.scatter(0.008,783.87,marker='s',s=8,lw='0.5',edgecolor='black',facecolor='white', zorder=3)
##    pyplot.xlabel(r'SMI (-)', labelpad=1)
##    pyplot.ylabel(r'curvature (%)', labelpad=3)
#
#    pyplot.xlabel(r'$\gamma_{\theta}$ (K m$^{-1}$)', labelpad=1)
#    pyplot.ylabel(r'max. SW$_{\mathrm{net}}$ (W m$^{-2}$)', labelpad=3)
#    pyplot.annotate('(a)', xy=(1.2, -0.15), xycoords='axes fraction', annotation_clip=False, size=10, color='black', va="center", ha="center")
#
##    host = host_subplot(111, axes_class=AA.Axes)
##    sndy = host.twinx()
##    new_fixed_axis = sndy.get_grid_helper().new_fixed_axis
##    sndy.axis["right"] = new_fixed_axis(loc="left", axes=sndy, offset=(-45,0))
##    sndy.axis["right"].toggle(all=True)
##    sndy.set_ylim(30,0)
##    sndy.set_ylabel("cloud cover (%)")
#
#
#    savepath = '/Users/mariecombe/Modeling/MXL/MLmodel/WUE1.jpeg'
#    fig.savefig(savepath,dpi=1000)


    ################################################################
    # figure 3: water stress effect on Rs and wce/wch

#    # defining X, Y, Z variable names, and the results file where these will be plotted from
#    xvar = 'smi'
#    yvar = 'curvature'
#    zvar1 = 'Rs_ave'
#    zvar2 = 'wce_h_int'
#    zvar3 = 'wcs_h_int'
#    results_filename = '2DIM_SA_cbeta_0-100_SMI_0-1/Results_wg_P4.dat'#'Results_wg_P4.dat'
#
#    # Rebuilding Y axis to produce the new Y range
#    #newaxis = remake_axis('cc', RangeDict['cc'], 'Swin_max', results_filename)        
#    #RangeDict[yvar] = newaxis[1]
#
#
#    print 'Plotting %s as colors, and %s as contour lines...' %(zvar1, zvar2)
#
#
#    # values of the X, Y, Z variables for plotting
#    X = list(RangeDict[xvar])
#    Y = list(RangeDict[yvar])
#
#    #Z1 = build_Z_grid(xvar, X, yvar, Y, zvar1, results_filename)
#    Z1 = numpy.reshape(ResDict[results_filename][zvar1], (-1,101))
#    #Z1 = Z1.T
#    Z1 = 1000./Z1
#    #Z1 = 1000./(1.6*Z1)
#    #Z1 = (Z1/10.289971393879524 - 1.)*100.
#    #Z1 = (Z1/16.46395423 - 1.)*100.
#    #Z1 = (Z1/7.3677999999999999 - 1.)*100.
#    #Z1 = (Z1/120.12719 - 1.)*100.
#    Z1 = numpy.ma.masked_invalid(Z1)
#
#    #Z2 = build_Z_grid(xvar, X, yvar, Y, zvar2, results_filename)
#    Z2 = numpy.reshape(ResDict[results_filename][zvar2], (-1,101))
#    #Z2 = Z2.T
##    Z2 = 1./Z2
##    Z1 = Z1/Z2
#    Z2 = numpy.ma.masked_invalid(Z2)
#    
#    #Z3 = build_Z_grid(xvar, X, yvar, Y, zvar3, results_filename)
#    Z3 = numpy.reshape(ResDict[results_filename][zvar3], (-1,101))
#    #Z3 = Z3.T
#    Z3 = numpy.ma.masked_invalid(Z3)
#    Z2 = Z2/Z3
#    import scipy.ndimage as ndimage
#    Z2 = ndimage.gaussian_filter(Z2, sigma=1.0, order=0)
#    
#    # box boundaries should be used with pcolormesh(), and box center with contour()
#    bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar])
#    X, Y = numpy.meshgrid(X,Y)
#
## Gs colorbars:
#    barlim=numpy.linspace(2.,17.,200)
#    import matplotlib.colors as colors
#    barnorm = colors.BoundaryNorm(barlim,256)
#    # plotting
#    fig = pyplot.figure('WUE1', figsize=(3.2,2.41))
#    fig.subplots_adjust(0.15,0.15,0.93,0.97,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
#    #pyplot.ticklabel_format(style='sci',axis='x', scilimits=(0,0))
#    pyplot.pcolormesh(bound[0], bound[1], Z1,vmin=2.,norm=barnorm)
#
## Gs colorbars:
#    barlim=[2.,5.,8.,11.,14.,17.]
#    cbar = pyplot.colorbar(extend='min',ticks=barlim)
#    cbar.set_label(r'$g_s$ (mm s$^{-1}$)') 
#    pyplot.xlim(min(bound[0]),max(bound[0]))
#    pyplot.ylim(min(bound[1]),max(bound[1]))
#    
#    levels= [0.5,1.,2.,4.,8.,16.]
#    cs= pyplot.contour(X, Y, Z2, levels, colors=['k','w','w','w','w','w','w'], linewidths=[0.5], linestyles=['-'])
#    pyplot.clabel(cs, fontsize=8, fmt='%.1f', colors=['k','w','w','w','w','w','w'])
#
#    # adding a decoupling threshold dotted line:
#    levelsz= [2.]
#    csz= pyplot.contour(X, Y, Z1, levelsz, colors=['r'], linewidths=[1.5], linestyles=['--'])
#    pyplot.clabel(csz, fontsize=0, fmt='%.0f', colors=['r'], manual=[(-100.,-100.)])
#
#    
##    #levels = [-2.00,-1.5,-1.0,-0.5]
##    cs2= pyplot.contour(X, Y, Z3, linestyles='--',colors='k', linewidths=0.5)
##    pyplot.clabel(cs2, fontsize=8, fmt='%.2f')
#
#    pyplot.scatter(0.2,0.,marker='o',s=8,lw='0.5',edgecolor='black',facecolor='white', zorder=3)
#    pyplot.annotate("DL",xy=(0.2,0.), xytext=(0.2,5.), textcoords='data',size=10, weight='bold', color='w', va="center", ha="center")
#    pyplot.scatter(0.8,0.,marker='o',s=8,lw='0.5',edgecolor='black',facecolor='black', zorder=3)
#    pyplot.annotate("WL",xy=(0.8,0.), xytext=(0.8,5.), textcoords='data',size=10, weight='bold', color='k', va="center", ha="center")
#    pyplot.scatter(0.2,100.,marker='o',s=8,lw='0.5',edgecolor='black',facecolor='black', zorder=3)
#    pyplot.annotate("DC",xy=(0.2,100.), xytext=(0.2,95.), textcoords='data',size=10, weight='bold', color='k', va="center", ha="center")
#    
#    #pyplot.scatter(0.008,783.87,marker='s',s=8,lw='0.5',edgecolor='black',facecolor='white', zorder=3)
#
#    pyplot.xlabel(r'SMI (-)', labelpad=1)
#    pyplot.ylabel(r'$C_{\beta}$ (%)', labelpad=3)
#
##    host = host_subplot(111, axes_class=AA.Axes)
##    sndy = host.twinx()
##    new_fixed_axis = sndy.get_grid_helper().new_fixed_axis
##    sndy.axis["right"] = new_fixed_axis(loc="left", axes=sndy, offset=(-45,0))
##    sndy.axis["right"].toggle(all=True)
##    sndy.set_ylim(30,0)
##    sndy.set_ylabel("cloud cover (%)")
#
#
#    savepath = '/Users/mariecombe/Modeling/MXL/MLmodel/RS1.jpeg'
#    fig.savefig(savepath,dpi=1000)

    ################################################################
    # figure X: crop yield in gC.m-2

#    # defining X, Y, Z variable names, and the results file where these will be plotted from
#    xvar = 'smi'
#    yvar = 'day'
#    zvar1 = 'DTU_SUM'
#    zvar2 = 'tm_max'
#    zvar3 = 'Rs_ave'
#    zvar4 = 'wce_h_int'
#    zvar5 = 'wcs_h_int'
#
#    results_filename = '2DIM_SA_day_130-230_SMI_0-1/Results_day_wg.dat'#'Saved.SENS_doy_smi/Results_day_wg.dat' #'Results_DOY_SMI.dat'
#
#    # Rebuilding Y axis to produce the new Y range
##    newaxis = remake_axis('cc', RangeDict['cc'], 'Swin_max', results_filename)        
##    RangeDict[yvar] = newaxis[1]
#
#
#    print 'Plotting %s as colors, and %s as contour lines...' %(zvar1, zvar2)
#
#
#    # values of the X, Y, Z variables for plotting
#    X = list(RangeDict[xvar])
#    Y = list(RangeDict[yvar])
#    X, Y = numpy.meshgrid(X,Y)
#
#    Z1 = numpy.reshape(ResDict[results_filename][zvar1], (-1,len(RangeDict[xvar])))
#    #Z1 = 1000./Z1
#    #Z1 = Z1/2.45
#    #Z1 = -Z1*12./44.*2.33
#    #Z1 = 1.1 + Z1
#    Z1 = numpy.ma.masked_invalid(Z1)
#
#    Z2 = numpy.reshape(ResDict[results_filename][zvar2], (-1,len(RangeDict[yvar])))
#    Z2 = numpy.ma.masked_invalid(Z2)
#    #Z2 = 1000./Z2
#    #Z2 = Z2/2.45
#    import scipy.ndimage as ndimage
#    Z2 = ndimage.gaussian_filter(Z2, sigma=1.0, order=0)
#
##    A = numpy.ma.masked_less(Z1,2.)
##    M = numpy.ma.getmask(A)
##    B = numpy.ma.array(Z2,mask=M)
#
#
#    Z3 = numpy.reshape(ResDict[results_filename][zvar3], (-1,len(RangeDict[yvar])))
#    Z3 = numpy.ma.masked_invalid(Z3)
#    Z3 = 1000./Z3
#    Z3 = ndimage.gaussian_filter(Z3, sigma=1.0, order=0)
##
##    Z4 = numpy.reshape(ResDict[results_filename][zvar4], (-1,len(RangeDict[yvar])))
##    Z4 = numpy.ma.masked_invalid(Z4)
##    Z5 = numpy.reshape(ResDict[results_filename][zvar5], (-1,len(RangeDict[yvar])))
##    Z5 = numpy.ma.masked_invalid(Z5)
##    Z6 = Z4/Z5
##    Z6 = numpy.ma.masked_outside(Z6,-50.,50.)
##    #Z6 = ndimage.gaussian_filter(Z6, sigma=1.0, order=0)
#
#
##    Z4 = numpy.abs(Z2)/(numpy.abs(Z2)+numpy.abs(Z3)) * 100.
#    
#    
#    # limitation of the parameter space
##    limday = 50 # DOY 180
##    limstyle = 'down'
##    if limstyle == 'down':
##        X=X[limday:,:]
##        Y=Y[limday:,:]
##        Z1 = Z1[limday:,:]
##        Z2 = Z2[limday:,:]
##        Z3 = Z3[limday:,:]
##        Z4 = Z4[limday:,:]
##        # box boundaries should be used with pcolormesh(), and box center with contour()
##        bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar][limday:])
##	
##    elif limstyle == 'up':
##        X=X[:limday,:]
##        Y=Y[:limday,:]
##        Z1 = Z1[:limday,:]
##        Z2 = Z2[:limday,:]
##        Z3 = Z3[:limday,:]
##        Z4 = Z4[:limday,:]
##        # box boundaries should be used with pcolormesh(), and box center with contour()
##        bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar][:limday])
#	
#    bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar])
#    
## Gs colorbars:
##    barlim=numpy.linspace(2.,17.,200)
##    import matplotlib.colors as colors
##    barnorm = colors.BoundaryNorm(barlim,256)
#
#    # plotting
#    fig = pyplot.figure('ANET', figsize=(3.2,2.41))
#    fig.subplots_adjust(0.15,0.15,0.93,0.97,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
##    pyplot.ticklabel_format(style='sci',axis='x', scilimits=(0,0))
#    pyplot.pcolormesh(bound[0], bound[1], Z1)#,vmin=0.,vmax=1.)
## Gs colorbars:
##    barlim=[2.,5.,8.,11.,14.,17.]
#    
#    cbar = pyplot.colorbar(format='%.2f')#extend='min',ticks=barlim)
#    cbar.set_label('RDVR (-)') #r'DCY (g$_{\mathrm{DM}}$ m$^{-2}$ d$^{-1}$)')  #r'max. $h$ (m)') #'EF (%)') #r'LEveg (mm day$^{-1}$)') #r'Crop transpiration (MJ m$^{-2}$ d$^{-1}$)') #Minimum CO$_2$ mole fraction')
#    pyplot.xlim(min(bound[0]),max(bound[0]))
#    pyplot.ylim(min(bound[1]),max(bound[1]))
#    
##    levels= [1.,2.,4.,8.,16.,32.]
#    #levels = [0.5,2.0,5.0,10.,15.0]
#    #levels = [710.]
#    #levels = [0.10,0.5,1.5,3.0]
#    colorzz = ['w','k','k','k','k','k'] #['k','k','k','k','w','w']
#    cs= pyplot.contour(X, Y, Z2, colors='k', linewidths=0.5, linestyles='-')
#    pyplot.clabel(cs, fontsize=8, fmt='%.0f',colors='k')
#    # add a decoupling line:
#    levelz = [2.]
#    cs2= pyplot.contour(X, Y, Z3, levelz, colors='k', linewidths=1.5, linestyles='--')
#    pyplot.clabel(cs2, fontsize=0, fmt='%.0f',colors='k', manual=[(-10.,-10.)])
#
#    
##    levels = [298.,300.,302.,304.,306.,308.]
#
#    #levels = [0.25,0.5,0.75,1.]
##    cs2= pyplot.contour(X, Y, Z4, linestyles='--',colors='k', linewidths=0.5)
##    pyplot.clabel(cs2, fontsize=8, fmt='%.2f')
#    
#    #pyplot.scatter(0.008,783.87,marker='s',s=8,lw='0.5',edgecolor='black',facecolor='white', zorder=3)
##    pyplot.xlabel(r'SMI (-)', labelpad=1)
##    pyplot.ylabel(r'curvature (%)', labelpad=3)
#
#
#    pyplot.scatter(0.2,216.,marker='o',s=8,lw='0.5',edgecolor='k',facecolor='k', zorder=3)
#    pyplot.annotate("DL",xy=(0.2,216.), xytext=(0.2,210.), textcoords='data',size=10, weight='bold', color='k', va="center", ha="center")
#    pyplot.scatter(0.8,216.,marker='o',s=8,lw='0.5',edgecolor='k',facecolor='k', zorder=3)
#    pyplot.annotate("WL",xy=(0.8,216.), xytext=(0.8,210.), textcoords='data',size=10, weight='bold', color='k', va="center", ha="center")
#
#    pyplot.annotate('(b)', xy=(1.2, -0.15), xycoords='axes fraction', annotation_clip=False, size=10, color='black', va="center", ha="center")
#
#    pyplot.xlabel(r'SMI (-)', labelpad=1)
#    pyplot.ylabel(r'time (DOY)', labelpad=3)
#
##    host = host_subplot(111, axes_class=AA.Axes)
##    sndy = host.twinx()
##    new_fixed_axis = sndy.get_grid_helper().new_fixed_axis
##    sndy.axis["right"] = new_fixed_axis(loc="left", axes=sndy, offset=(-45,0))
##    sndy.axis["right"].toggle(all=True)
##    sndy.set_ylim(30,0)
##    sndy.set_ylabel("cloud cover (%)")
#
#
#    savepath = '/Users/mariecombe/Modeling/MXL/MLmodel/ANET.png'
#    fig.savefig(savepath,dpi=300)

    ################################################################
    # figure 3: ?? with x=wg0, y=gammaq, and contour plot is BL_h

#    # defining X, Y, Z variable names, and the results file where these will be plotted from
#    xvar = 'smi'
#    yvar = 'wsls'
#    zvar1 = 'WUEint_14h'
#    zvar2 = 'evafra_14h'
#    results_filename = 'Results_wg_wsls.dat'
#
#
#    print 'Plotting %s as colors, and %s as contour lines...' %(zvar1, zvar2)
#
#
#    # values of the X, Y, Z variables for plotting
#    X = list(RangeDict[xvar])
#    Y = list(RangeDict[yvar])
#
#    #Z1 = build_Z_grid(xvar, X, yvar, Y, zvar1, results_filename)
#    Z1 = numpy.reshape(ResDict[results_filename][zvar1], (-1,101))
#    Z1 = Z1.T
#    Z1 = numpy.ma.masked_invalid(Z1)
#
#    #Z2 = build_Z_grid(xvar, X, yvar, Y, zvar2, results_filename)
#    Z2 = numpy.reshape(ResDict[results_filename][zvar2], (-1,101))
#    Z2 = Z2.T/100.
#    #Z2 = (Z2-0.06)/(0.15-0.06)
#    Z2 = numpy.ma.masked_invalid(Z2)
#
#        
#    # box boundaries should be used with pcolormesh(), and box center with contour()
#    bound = make_pcolormesh_box_boundaries(RangeDict[xvar],RangeDict[yvar])
#    X, Y = numpy.meshgrid(X,Y)
#
## WUE colorbars:
##    barlim=numpy.linspace(2.75,3.20,200)
##    barlim=numpy.linspace(4.50,4.95,200)
## max temp colorbars:    
##    barlim=numpy.linspace(295.0,303.0,200)   
## gs colorbars:
##    barlim=numpy.linspace(11.5,12.4,200)
##    barlim=numpy.linspace(1.8,2.7,200)
##    import matplotlib.colors as colors
##    barnorm = colors.BoundaryNorm(barlim,256)
#
#    # plotting
#    fig = pyplot.figure('wsls', figsize=(3.1,2.41))
#    fig.subplots_adjust(0.15,0.14,0.95,0.95,0.4,0.) # left,bottom,right,top, vertical spaces, horizontal spaces
#    pyplot.ticklabel_format(style='sci',axis='y', scilimits=(0,0))
#    pyplot.pcolormesh(bound[0], bound[1], Z1)#,norm=barnorm)#, cmap=cm.jet_r)#)
## WUE colorbars:
##    barlim=[2.75,2.80,2.85,2.90,2.95,3.00,3.05,3.10,3.15,3.20]
##    barlim=[4.50,4.55,4.60,4.65,4.70,4.75,4.80,4.85,4.90,4.95]
## max temp colorbars:
##    barlim=[295.,296.,297.,298.,299.,300.,301.,302.,303.]
## gs colorbars:
##    barlim=[11.5,11.6,11.7,11.8,11.9,12.0,12.1,12.2,12.3,12.4]
##    barlim=[1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7]
#    pyplot.colorbar()#ticks=barlim)
#    pyplot.xlim(min(bound[0]),max(bound[0]))
#    pyplot.ylim(min(bound[1]),max(bound[1]))
#    
#    cs= pyplot.contour(X, Y, Z2, colors='k', linewidths=0.5, linestyles='-')
#    pyplot.clabel(cs, fontsize=8, fmt='%.1f')
#    
#    pyplot.scatter(0.556,0.000007,marker='o',s=8,lw='0.5',edgecolor='black',facecolor='black', zorder=3)
#    pyplot.annotate("C",xy=(0.556,0.000007), xytext=(0.556, 0.000009), textcoords='data',size=10, weight='bold', color='k', va="center", ha="center")
#    pyplot.scatter(0.506,0.000007,marker='o',s=8,lw='0.5',edgecolor='black',facecolor='RoyalBlue', zorder=3)
#    pyplot.annotate("D",xy=(0.556,0.000007), xytext=(0.506, 0.000009), textcoords='data',size=10, weight='bold', color='RoyalBlue', va="center", ha="center")
#    pyplot.scatter(0.556,0.00004,marker='o',s=8,lw='0.5',edgecolor='black',facecolor='MediumVioletRed', zorder=3)
#    pyplot.annotate("H",xy=(0.556,0.000007), xytext=(0.556, 0.0000375), textcoords='data',size=10, weight='bold', color='MediumVioletRed', va="center", ha="center")
#    
#    
#    #pyplot.scatter(0.008,783.87,marker='s',s=8,lw='0.5',edgecolor='black',facecolor='white', zorder=3)
#    pyplot.xlabel(r'SMI (-)', labelpad=1)
#    pyplot.ylabel(r'$D$ (s$^{-1}$)', labelpad=1)
#
#
#    savepath = os.path.join('/Users/mariecombe/Modeling/MXL/MLmodel','wsls.png')
#    fig.savefig(savepath,dpi=300)
    
    print "done!"
