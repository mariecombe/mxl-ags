#!/usr/bin/env python
# Create_Sensitivity_File.py


"""
Started 15 octobre 2014
Author: Marie Combe

Script to create the sensitivity analysis range for DOY

"""

# importing the necessary modules

import sys
import os
import getopt
import shutil
#from csv import *		# reading .csv files
import numpy		# for arrays handling
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# custom-made functions
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
def half_sigmoid(x, slope):
    y = 1.9 / (1+ numpy.exp(-slope*(x-185.)))
    return y

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
    Script that generates namoptions file for the sensitivity analysis of DOY
                
                """
            
            print helptext
            
            sys.exit(2)      


    # Defining working directories

    currentdir = os.getcwd()
    gecrosdir  = '/Users/mariecombe/Modeling/gecros-mxl/branches/marie-gecrosmxl-branch/'
    
    # Define variables of sensitivity analysis
    Var = 'day'       
    start = 130.
    end = 230.
    step = 1.

    # list of namoptions files for sensitivity analysis
    
    listofnamoptions = [f for f in os.listdir( currentdir ) if ( f.startswith('namoptions_') ) ]

    # variables to vary with consistent changes

    Range = numpy.arange(start,end+step/2.,step)   # it will skip the last value

    # opening LAI data files 
    
    Mxl_gecros_resdat    = open_output(currentdir,['uncoupled_res.dat'])
    CropObsDict          = open_csv(currentdir,['LAI_obs_6-aug-2007.csv','MeanT150cm_obs_clim_1951-2000.csv',
                                                'MaxT150cm_obs_clim_1951-2000.csv','MinT150cm_obs_clim_1951-2000.csv'])
        
    # fitting a sigmoid curve / interpolating: GECROS LAI data 
       
    xA=Mxl_gecros_resdat['uncoupled_res.dat']['TIME']
    yA=Mxl_gecros_resdat['uncoupled_res.dat']['LAI']
    lai_interpA = UnivariateSpline(xA, yA, s=1, k=2)
    poptA, pcovA = curve_fit(sigmoid, xA, yA)    

    # fitting a sigmoid curve / interpolating: observed LAI data 
       
    xB=CropObsDict['LAI_obs_6-aug-2007.csv']['DOY']
    yB=CropObsDict['LAI_obs_6-aug-2007.csv']['LAI']
    lai_interpB = UnivariateSpline(xB, yB, s=1, k=2)
    poptB, pcovB = curve_fit(sigmoid, xB, yB)
    
    # interpolating GECROS developmental stage number
    
    xF=Mxl_gecros_resdat['uncoupled_res.dat']['TIME']
    yF=Mxl_gecros_resdat['uncoupled_res.dat']['DVS']
    dvs_interp = UnivariateSpline(xF, yF, s=1, k=1)

    # interpolating: observed mean temperature 
       
    xC=CropObsDict['MeanT150cm_obs_clim_1951-2000.csv']['DOY']
    yC=CropObsDict['MeanT150cm_obs_clim_1951-2000.csv']['average_T']
    lai_interpC1 = UnivariateSpline(xC, yC, s=1, k=2)
    lai_interpC2 = UnivariateSpline(xC, yC, s=1, k=3)

    # interpolating: observed max temperature 
       
    xD=CropObsDict['MaxT150cm_obs_clim_1951-2000.csv']['DOY']
    yD=CropObsDict['MaxT150cm_obs_clim_1951-2000.csv']['average_T']
    lai_interpD1 = UnivariateSpline(xD, yD, s=1, k=2)
    lai_interpD2 = UnivariateSpline(xD, yD, s=1, k=3)
    
    # interpolating: observed min temperature 
       
    xE=CropObsDict['MinT150cm_obs_clim_1951-2000.csv']['DOY']
    yE=CropObsDict['MinT150cm_obs_clim_1951-2000.csv']['average_T']
    lai_interpE1 = UnivariateSpline(xE, yE, s=1, k=2)
    lai_interpE2 = UnivariateSpline(xE, yE, s=1, k=3)


    x2=numpy.linspace(1., 365., 365.)
    x3=numpy.linspace(130.,230.,100.)
    
    # plotting the LAI curves

    plt.close('all')

#    plt.figure(1)
#    plt.plot(xA,yA,'o')
#    plt.plot(x2,sigmoid(x2, *poptA),'-')
#    plt.plot(x2,lai_interpA(x2),'--')
#    plt.title('GECROS LAI')
#    plt.ylim(-1.,6.)
#    plt.legend(['data', 'sigmoid fit', 'smoothed 2nd order spline interpolation'], loc='best')

#    plt.figure(2)
#    plt.plot(xB,yB, 'o')
#    plt.plot(x2,sigmoid(x2, *poptB),lw=2,ls='-',c='r')
##    plt.plot(x2,lai_interpB(x2),'--')
#    plt.ylim(0.,4.)
#    plt.yticks([0.,1.,2.,3.,4.])
#    plt.xlim(120.,240.)
##    plt.title('Observed LAI')
#    plt.legend(['observed LAI','modeled LAI'], loc='best', fontsize=14)

#    plt.figure(3)
##    plt.plot(xA, yA,'o',c='blue')
##    plt.plot(xB, yB,'^',c='red')
#    plt.plot(x2,sigmoid(x2, *poptA),'-', lw=2, c='blue')
#    plt.plot(x2,sigmoid(x2, *poptB),'-', lw=2, c='red')
#    plt.plot(x2,lai_interpA(x2),'--', lw=2, c='blue')
#    plt.plot(x2,lai_interpB(x2),'--', lw=2, c='red')
#    plt.ylim(-1.,7.)
#    plt.title('LAI')
#    plt.legend([
##                'LAI data GECROS', 'LAI data obs', 
#		'sigmoid fit GECROS', 'sigmoid fit obs', 
#                'smoothed 2nd order spline interpolation GECROS', 'smoothed 2nd order spline interpolation obs'], loc='best')
#    
    plt.figure(4)
    plt.plot(xD,yD,'o', c='r', label='obs max T')
#    plt.plot(x2,lai_interpD1(x2),'-')
    plt.plot(x2,lai_interpD2(x2),'--', c='r')

    plt.plot(xC,yC,'o',c='green', label = 'obs mean T')
#    plt.plot(x2,lai_interpC1(x2),'-')
    plt.plot(x2,lai_interpC2(x2),'--',c='green')

#    plt.plot(xC,yC-2.,'^',c='k')
    plt.plot(x2,lai_interpC2(x2)-2.,'-', c='k', label='modeled 6:00UTC T')
        
    plt.plot(xE,yE,'o', c='b', label='obs min T')
#    plt.plot(x2,lai_interpE1(x2),'-')
    plt.plot(x2,lai_interpE2(x2),'--', c='b')
    
    plt.ylim(5.,30.)
    plt.xlim(120.,240.)
    plt.xlabel('DOY')
    plt.legend(loc='best')
#    plt.title('Climatic averages (1951-2000) of\nmax, min, mean monthly temperatures')

#    plt.figure(5)
#    plt.plot(xF,yF,'o',label='data')
#    plt.plot(x3,dvs_interp(x3),'-', label='interpolation')
#    plt.legend(loc='best')
#    plt.title('Crop Developmental stage')
#    
    plt.show()

    # NB: from the exercise done above, I chose to use: (16 Oct 2014)
    # LAI as a sigmoid curve: y = 3.8 / (1 + numpy.exp(-slope*(DOY-185.))), with a fitted slope of 0.12523755
    # thetam0 as a smoothed spline interpolation of the 3rd order of the mean monthly T minus 2K
