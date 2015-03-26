#!/usr/bin/env python
# Create_Sensitivity_File.py


"""
Started 22 octobre 2013
Author: Marie Combe

Script to create the sensitivity analysis files

"""

# importing the necessary modules

import sys
import os
import getopt
import shutil
#from csv import *		# reading .csv files
import numpy		# for arrays handling


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
    Script that generates namoptions file for sensitivity analysis
                
                """
            
            print helptext
            
            sys.exit(2)      


    # Defining working directories

    currentdir = os.getcwd()

    # Define variables of sensitivity analysis
    Var = 'thetam0'       
    start = 284.
    end = 290.
    step = 0.06 #0.12

    # list of namoptions files for sensitivity analysis
    listofnamoptions = [f for f in os.listdir( currentdir ) if ( f.startswith('namoptions_') ) ]

    # variables to vary with consistent changes

    Range = numpy.arange(start,end+step/2.,step)   # it will skip the last value

    if (Var == 'thetam0'):
        Tsoil = numpy.arange(start+2., end+2.+step/2.,step)
        T2 = numpy.arange(start+3., end+3.+step/2.,step)
        #qm0 = numpy.array([0.]*len(Range))
        qm0 = 1000.* 0.942 * ( 0.622 * 0.611 * numpy.exp( (2.45e6/461.)*( (1./273.15)-(1./ Range) ) ) ) /102.2

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
                if line.startswith('outdir  '): f.write("outdir     = 'SENS_%s%03d'\n"%(Var[0:2],x+1))
                elif line.startswith(Var): f.write(Var+'    =  %10.5f\n'%(variant) )
                elif (Var=='thetam0' and line.startswith('Tsoil  ')): f.write('Tsoil      =  %10.2f\n'%(Tsoil[x]) )
                elif (Var=='thetam0' and line.startswith('T2  ')): f.write('T2         =  %10.2f\n'%(T2[x]) )
                elif (Var=='thetam0' and line.startswith('qm0  ')): f.write('qm0        =  %5.1f\n'%(qm0[x]) )
                else: f.write(line)
            f.close()

