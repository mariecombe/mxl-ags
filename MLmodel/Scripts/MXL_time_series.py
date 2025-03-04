#!/usr/bin/env python
# ShowResGECROS.py


"""
Started 30 April 2015
Author: Marie Combe

Time series script for the MXL res.dat file

"""

import sys
import os
import getopt
import numpy		# for arrays handling
from matplotlib import pyplot as plt

#----------------------------------------------------------------------
def open_MXLoutput(curdir,outputfileslist,folder,dictnamelist):
    """
Routine to open a series of MXL result files (output_dyn files)
You need to specify :
        curdir              :   the name of the current directory
        outputfileslist :   a list of folder names in which to find the output
                                e.g. ['DOY_110','DOY_111']
        dictnamelist        :   a list of empty dictionnaries (1 per 
                                result file)
The open_output_dyn function will return one dictionnary with the data
from all result files
    """
    dictnamelist = {}
    a=len(outputfileslist)

    # for each result file:
    for i in range (0,a):
        # open file, read all lines
        inputpath = os.path.join(currentdir,folder,outputfileslist[i])
        f=open(inputpath,'rU') 
        lines=f.readlines()
        f.close()
        # storing headers in list headerow
        headerow=lines[0].strip().split()
        # deleting rows that are not data (first and last rows of the file)
        del lines[0:3]
        del lines[len(lines)-1]
        # transforming data from string to float type, storing it in array 'data'
        datafloat=[]
        for line in lines:
            a= line.strip().split()
            datafloat.append(line.strip().split())
        for j in range(0,len(datafloat)):
            for k in range(0,len(datafloat[j])):
                if (('-' in datafloat[j][k] or '+' in datafloat[j][k]) and not 'E' in datafloat[j][k]): datafloat[j][k] = 'NaN'
        data=numpy.array(datafloat,float)
        # creating one dictionnary and storing the float data in it
        dictnamelist[outputfileslist[i]] = {}
        for j,varname in enumerate(headerow):
            dictnamelist[outputfileslist[i]][varname]=data[:,j]
    
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
This script plots the figures of the methods section of Article 2

You do not need to specify any arguments to run this script:
./Article2_methods_figures.py
                
                """
            
            print helptext
            
            sys.exit(2)

    # Current directory
    
    currentdir = "/Users/mariecombe/Modeling/MXL/MLmodel/DROUGHT"
    folderlist = [f for f in os.listdir( currentdir ) if ( f.startswith('SENS_') ) ]
    outputfileslist = ['output_dyn','output_sca']

    
    cases = ['P4001','P4101']
    subtitles = [r'Drought-sensitive vegetation (linear $\beta$)',r'Drought-insensitive vegetation (curved $\beta$)']

    # Reading the output
    
    Results = {}
    for case in cases:
        for folder in folderlist:
            if (case in folder):
                # Read all results files and store all dictionnaries into one dictionnary
                Results[folder] = open_MXLoutput(currentdir,outputfileslist,folder,[])
                Results[folder]['output_dyn']['new_time(h)'] = Results[folder]['output_dyn']['RT(hours)'] + (float(folder[12:15])-1.)*24.
    print "Read all output! \n"
    
    plt.close('fig1')
    plt.close('fig2')

    
    # Figure 1: LE and wg
    
    fig = plt.figure('fig1', figsize=(16,8)) 
    a=210
    for i,case in enumerate(cases):
        a = a + 1
	ax1 = plt.subplot(a)
        plt.ylabel(r'heat sflux (W m$^{-2}$)')
        plt.ylim([0.,400.])
        ax2 = ax1.twinx()
    
        for folder in folderlist:
            if ((case in folder) and ("wg001" in folder)): 
	        ax1.plot(Results[folder]['output_dyn']['new_time(h)'], Results[folder]['output_dyn']['LE(W.m-2)'], c='k', label='LE')
	        ax1.plot(Results[folder]['output_dyn']['new_time(h)'], Results[folder]['output_dyn']['LEveg(W.m-2)'], ls='--', c='r', label='LEveg')
    	        ax2.plot(Results[folder]['output_dyn']['new_time(h)'], Results[folder]['output_sca']['w2(cm3.cm-3)'], c='b', label='w2')
    	        #ax2.plot(Results[folder]['output_dyn']['new_time(h)'], Results[folder]['output_sca']['w2(cm3.cm-3)'], ls='--', c='b', label='w2')
	    elif (case in folder):
	        ax1.plot(Results[folder]['output_dyn']['new_time(h)'], Results[folder]['output_dyn']['LE(W.m-2)'], c='k')
	        ax1.plot(Results[folder]['output_dyn']['new_time(h)'], Results[folder]['output_dyn']['LEveg(W.m-2)'], ls='--', c='r')
	        ax2.plot(Results[folder]['output_dyn']['new_time(h)'], Results[folder]['output_sca']['w2(cm3.cm-3)'], c='b')
    	        #ax2.plot(Results[folder]['output_dyn']['new_time(h)'], Results[folder]['output_sca']['w2(cm3.cm-3)'], ls='--', c='b')
        plt.xlabel('time (h)')
        plt.ylabel(r'soil moisture (cm$^{3}$ cm$^{-3}$)')
        ax1.legend(loc='upper right')
        ax2.legend(loc='center right')
        plt.xlim([0.,530.])
        plt.ylim([0.06,0.11])
        plt.title(subtitles[i])

    plt.xlabel('time (h)')

    # Figure 2: h and theta
    
    fig = plt.figure('fig2', figsize=(16,8)) 
    a=210
    for i,case in enumerate(cases):
        a = a + 1
	ax1 = plt.subplot(a)
        plt.ylabel('ABL height (m)')
        plt.ylim([0.,1800.])
        ax2 = ax1.twinx()
    
        for folder in folderlist:
            if ((case in folder) and ("wg001" in folder)): 
	        ax1.plot(Results[folder]['output_dyn']['new_time(h)'], Results[folder]['output_dyn']['zi(m)'], c='k', label='h')
    	        ax2.plot(Results[folder]['output_dyn']['new_time(h)'], Results[folder]['output_dyn']['thetam(K)'], c='b', label=r'$\theta$')
	    elif (case in folder):
	        ax1.plot(Results[folder]['output_dyn']['new_time(h)'], Results[folder]['output_dyn']['zi(m)'], c='k')
	        ax2.plot(Results[folder]['output_dyn']['new_time(h)'], Results[folder]['output_dyn']['thetam(K)'], c='b')
        plt.xlabel('time (h)')
        plt.ylabel(r'$\theta$ (K)')
        ax1.legend(loc='upper right')
        ax2.legend(loc='center right')
        plt.xlim([0.,530.])
        #plt.ylim([0.06,0.11])
        plt.title(subtitles[i])

    plt.xlabel('time (h)')

	
    plt.show()
	
	

