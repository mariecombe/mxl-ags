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
from scipy.interpolate import InterpolatedUnivariateSpline, interp1d, UnivariateSpline, splrep, splev
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
    CropObsDict          = open_csv(currentdir,['LAI_obs_6-aug-2007.csv'])
        

    # fitting a sigmoid curve / interpolating: observed LAI data 
       
    xB=CropObsDict['LAI_obs_6-aug-2007.csv']['DOY']
    yB=CropObsDict['LAI_obs_6-aug-2007.csv']['LAI']
    popt, pcov = curve_fit(sigmoid, xB, yB)    



    if (Var=='day'):
        LAI = sigmoid(Range, *popt)
	CVEG = 1.- numpy.exp(-LAI)
    
    for x,variant in enumerate(Range):
        for i, filename in enumerate(listofnamoptions):
            print filename+'%s%03d'%(Var[0:2],x+1)
	    
            # read the content of the old file - and close it
            f=open(os.path.join(currentdir,filename),'rU') 
            lines=f.readlines()
            f.close()
	    
            # copy the old file to a new name with the variable tag
            shutil.copyfile(os.path.join(currentdir,filename),os.path.join(currentdir,filename+'%s%03d'%(Var[0:2],x+1) ))
	    
            # write in the new file all old lines except line with variable with new value
            f=open(os.path.join(currentdir,filename+'%s%03d'%(Var[0:2],x+1) ), 'w')
            datafloat=[]
            for j,line in enumerate(lines):
                if line.startswith('outdir  '): 
		    if "OBS04AUG2007" in line: 
		        f.write("outdir     = 'SENS_%s%03d'\n"%(Var[0:2],x+1))
		    elif 'SENS' in line: 
		        a= line.strip().split()
			f.write("outdir     = 'SENS_%s%03d%s'\n"%(Var[0:2],x+1,a[2][6:11]))
                elif line.startswith(Var): f.write(Var+'    =  %d\n'%(variant) )
                elif (Var=='day' and line.startswith('LAI  ')): f.write('LAI      =  %10.5f\n'%(LAI[x]) )
		elif (Var=='day' and line.startswith('cveg  ')): f.write('cveg      =  %10.5f\n'%(CVEG[x]) )
                else: f.write(line)
            f.close()


