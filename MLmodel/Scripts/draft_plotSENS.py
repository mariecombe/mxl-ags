#!/usr/bin/env python
# sensitivity.py


"""
Started 25 june 2013
Author: Marie Combe

Plotting script for the sensitivity analysis results

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
from matplotlib import colors
from matplotlib.backends.backend_pdf import PdfPages



# custom-made functions
#----------------------------------------------------------------------
def open_data(inpath,namefile):

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

    return dictnamelist


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

You should specify the following arguments to run this script:
    - Results file path

As an example:
./plotSENS.py Results_wg_thetam0
                
                """
            
            print helptext
            
            sys.exit(2)

    if len(args)<1:
        print "FAILURE\nYou forgot to write (some) script arguments!! See help text (call it with -h)\n"
        sys.exit(2)

    # Get arguments
    ResultsFile = args[0]

    if not os.path.exists(ResultsFile): # if path to result file does not exist
        print "Does not exist! (%s)" % ResultsFile
        sys.exit(2)

    # Define variables of sensitivity analysis
    #Var1 = 'wg'
    #Range1 = numpy.arange(0.096,0.1231,0.00027) # it will skip the last value
    #Var2 = 'P4'
    #Range2 = numpy.arange(0.,15.15,0.15)
    #Var2 = 'ga'
    #Range2 = numpy.arange(-0.005,0.00005,0.00005) # it will skip the last value
    #Range2 = numpy.arange(0.002,0.00806,0.00006) #gammatheta range
    #Var1 = 'thetam0'    #'thetam0'                                                            
    #Range1 = numpy.arange(280.,295.15,0.15)   #(280.,291.,1.) # it will skip the last value
    #Var1 = 'cc'
    #Range1 = numpy.arange(0.000,0.303,0.003)
    #Var2 = 'la'
    #Range2 = numpy.arange(1.00,4.01,0.03)

    # Working directories
    currentdir = os.getcwd()
    #sensdir = 'MXLmodel/COUPLED_RUNS/SENS_'+Var1+'_'+Var2

    # Open results
    ResDict = open_data(currentdir,ResultsFile)


    # Plot results

    # First close all opened figures in python
    pyplot.close('all')

    xvar = Var1
    yvar = Var2
    if (yvar == 'day'): yvar = 'Qnet_max' # if we have varied the DOY, we replace 'day' by 'Qnet'

    # We build a grouped column in the case we define the x or y axis of our plot different from Var1 or Var2 
    # this new column is made of groups of values of the changed x or y axis. We have as many classes/groups of values as Var1/Var2 values.
    # this way we end up having identical dimensions for the X and Y arrays, and the meshgrid function will work properly. 

    for axis in [xvar,yvar]:
        if (axis==xvar and axis != Var1): # if we have defined the x-axis as something else than Var1 or Var2
            Range_ = Range1
            Var_ = Var1 
        elif (axis==yvar and axis != Var2): # if we have defined the y-axis as something else than Var1 or Var2
            Range_ = Range2
            Var_ = Var2
        else:
            continue
        print "Rebuilding the axis '%s' to replace the original '%s' axis"%(axis, Var_)


        W = [0.]*len(Range_)
        for i,day in enumerate(Range_):
            V = []
            for j,qnet in enumerate(ResDict[axis]):
                if (ResDict[Var_][j]==day):
                    V = V + [qnet]
            W[i] = numpy.mean(numpy.array(V))

        ResDict['grouped_%s'%axis] = [0.]*len(ResDict[axis])
        for i,day in enumerate(Range_):
            for j,qnet in enumerate(ResDict[axis]):
                if (ResDict[Var_][j]==day):
                    ResDict['grouped_%s'%axis][j] = W[i]
        
    if (yvar != Var2): 
        yvar = 'grouped_%s'%yvar
	Range2 = W
    
    # open pdf file
    pdf_pages = PdfPages(os.path.join(currentdir,'SA_report_'+Var1+'_'+Var2+'.pdf'))

    headers = ['Qnet_max','tm_range', 'tm_max', 'qm_range', 'qm_ave', 'co2_range', 'co2_min', 'h_max', 'lcl_h_min', 'TDIF_range', 'TDIF_max', 'VPD_range', 
               'VPD_max', 'CDIF_range', 'CDIF_max', 'WUEplt_12h', 'WUEeco_12h', 'Anet_min', 'Resp_max', 'wce_max', 'NEE_min', 'Ra_12h', 'Rs_12h', 'SH_max', 'LE_max', 
               'GR_max', 'beta_12h', 'evafra_12h','wg_range']
    units =   ['[W.m-2]','[K]', '[K]', '[g.kg-1]', '[g.kg-1]', '[ppm]', '[ppm]', '[m]', '[m]', '[K]', '[K]', '[g.kg-1]', 
               '[g.kg-1]', '[ppm]', '[ppm]', '[???]', '[???]', '[mgCO2.m-2.s-1]', '[mgCO2.m-2.s-1]', '[mgCO2.m-2.s-1]', '[mgCO2.m-2.s-1]', '[s.m-1]', '[s.m-1]', 
               '[W.m-2]', '[W.m-2]', '[W.m-2]', '[-]', '[-]', '[cm3.cm-3]']
    for l,var in enumerate(headers):
        if (var!=Var1 and var!=Var2 and not var.startswith('grouped') and not var.startswith('Swin') and not var.startswith('Qnet')):
            zvar = var
            print 'Plotting ',zvar
            pyplot.close('all')

            X = list(Range1)
            Y = list(Range2)

            # Now that we have proper X and Y arrays we build the 2D array containing the values of zvar 

            Z= numpy.reshape(ResDict[zvar], (-1,101))
	    Z = Z.T


            Z = numpy.ma.masked_invalid(Z) # we want to avoid plotting NaN, Inf, -Inf...
            realminZ = round(numpy.ma.min(Z), 2)
            realmaxZ = round(numpy.ma.max(Z), 2)
            realmeanZ = round(numpy.mean(Z), 2)

            QT = mquantiles(Z) # quartiles
            IQR = QT[2]-QT[0]  #  interquartile range
            #Z = numpy.ma.masked_less(Z, QT[0] - IQR*1.5)   # masking the outliers on the plot
            #Z = numpy.ma.masked_greater(Z, QT[2] + IQR*1.5)

            X, Y = numpy.meshgrid(X,Y)

            # Plot the 3D plot in 2D

            rc('xtick', labelsize=16)
            rc('ytick', labelsize=16)


            fig = pyplot.figure(l, figsize=(9,6))
            fig.subplots_adjust(0.12,0.10,0.98,0.78,0.3,0.2) # left,bottom,right,top, vertical spaces, horizontal spaces
            fig.suptitle('MXL-A-Gs: '+zvar+' '+units[l] + '\n\nmin: '+str(realminZ)+
                         ', max: '+str(realmaxZ)+', mean = '+ str(realmeanZ) +'\n', fontsize=16)

            gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 
            sub1 = pyplot.subplot(gs[0])
            pyplot.pcolormesh(X, Y, Z)

            pyplot.colorbar()
            pyplot.xlim(min(ResDict[xvar]),max(ResDict[xvar]))
            pyplot.ylim(min(ResDict[yvar]),max(ResDict[yvar]))
            pyplot.xlabel(xvar, fontsize=16)
            pyplot.ylabel(yvar, fontsize=16)
            pyplot.title('Space\n', fontsize=16)

            sub2 = pyplot.subplot(gs[1])
            pyplot.boxplot(numpy.ravel(Z))
            pyplot.scatter(1,numpy.mean(Z), c='red')
            pyplot.xlabel(zvar, fontsize=16)
            pyplot.title('Distribution\n', fontsize=16)

            pdf_pages.savefig()
            #pyplot.show()

    pdf_pages.close()
    print 'SUCCESS!! Pdf saved in %s' %currentdir



