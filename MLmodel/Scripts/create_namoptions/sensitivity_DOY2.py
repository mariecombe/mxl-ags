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
    obsdir = '/Users/mariecombe/Modeling/MXL/MLmodel/Files_csv'
    
    # Define variables of sensitivity analysis
    Var = 'day'       
    start = 130.
    end = 230.
    step = 1.
    temp_separate = True # if this switch is true, we make thetam0 vary in the x-axis, separate from LAI and DOY

    # list of namoptions files for sensitivity analysis
    
    listofnamoptions = [f for f in os.listdir( currentdir ) if ( f.startswith('namoptions_') ) ]

    # variables to vary with consistent changes

    Range = numpy.arange(start,end+step/2.,step)   # it will skip the last value

    # opening LAI data files 
    
#    Mxl_gecros_resdat    = open_output(obsdir,['uncoupled_res.dat'])
    CropObsDict          = open_csv(obsdir,['LAI_obs_6-aug-2007.csv','MeanT150cm_obs_clim_1951-2000.csv'])
        

    # fitting a sigmoid curve through the observed LAI data 
       
    xA=CropObsDict['LAI_obs_6-aug-2007.csv']['DOY']
    yA=CropObsDict['LAI_obs_6-aug-2007.csv']['LAI']
    popt, pcov = curve_fit(sigmoid, xA, yA)    

    # interpolating the mean monthly temperature 
       
    xB=CropObsDict['MeanT150cm_obs_clim_1951-2000.csv']['DOY']
    yB=CropObsDict['MeanT150cm_obs_clim_1951-2000.csv']['average_T']-2.
    thetam_interp = UnivariateSpline(xB, yB, s=1, k=3)

    #------------------------------------------------
    # Create consistent changes in DOY, LAI, thetam0:
    #------------------------------------------------
    
    LAI = sigmoid(Range, *popt)
    CVEG = 1.- numpy.exp(-LAI)
    
    THETAM0 = thetam_interp(Range)+273.15
    TSOIL = THETAM0 + 2.
    T2 = THETAM0 + 3.
    DTU0 = numpy.array([0.]*len(Range))
    ALPHA = CVEG * 0.198 + (1-CVEG) * 0.15
    
    # if we assume that sowing happened at DOY 135, flowering at DOY 195, and maturity at DOY 282:
#    for i,doy in enumerate(Range):
#        if (doy < 135.): DTU0[i] = 0.
#	elif (135. <= doy <= 195.): DTU0[i] = (doy-135.)/60.
#        elif (195. < doy <= 282.): DTU0[i] = (doy-195.)/87.+1.
#        else: DTU0[i] = 2.
    
    # because we want a constant relative humidity over the runs, we modify the following:
    QM0 = 1000.* 0.942 * ( 0.622 * 0.611 * numpy.exp( (2.45e6/461.)*( (1./273.15)-(1./ THETAM0) ) ) ) /102.2

    
    # Create namoption files:
    
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
                elif line.startswith(Var):                          f.write(Var+'    =  %d\n'%(variant) )
                elif (Var=='day' and line.startswith('LAI  ')):     f.write('LAI     =  %10.5f\n'%(LAI[x]) )
		elif (Var=='day' and line.startswith('cveg  ')):    f.write('cveg    =  %10.5f\n'%(CVEG[x]) )
		elif (Var=='day' and line.startswith('thetam0  ')): f.write('thetam0 =  %10.5f\n'%(THETAM0[x]) )
		elif (Var=='day' and line.startswith('Tsoil  ')):   f.write('Tsoil   =  %10.5f\n'%(TSOIL[x]) )
		elif (Var=='day' and line.startswith('T2  ')):      f.write('T2      =  %10.5f\n'%(T2[x]) )
		elif (Var=='day' and line.startswith('qm0  ')):     f.write('qm0     =  %10.5f\n'%(QM0[x]) )
		elif (Var=='day' and line.startswith('albedo  ')):  f.write('albedo  =  %10.5f\n'%(ALPHA[x]) )
#		elif (Var=='day' and line.startswith('DS0  ')):     f.write('DS0     =  %10.5f\n'%(DTU0[x]) )
                else: f.write(line)
            f.close()


