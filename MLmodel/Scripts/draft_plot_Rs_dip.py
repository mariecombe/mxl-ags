#!/usr/bin/env python
# ShowResGECROS.py


"""
Started 14 april 2014
Author: Marie Combe

Plotting script for MXL-A-gs multiple runs

"""

import getopt
import sys
import os
import numpy
from matplotlib.pyplot import figure,plot,show,close,title,xlabel,ylabel,xlim,ylim
from matplotlib import rc

#----------------------------------------------------------------------
def open_output(path,inputfoldernamelist):
    """
Opening results of MXL-AGS runs
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
        for j in range(0,len(datafloat)):
            for k in range(0,len(datafloat[j])):
                if (('-' in datafloat[j][k] or '+' in datafloat[j][k]) and not 'E' in datafloat[j][k]): datafloat[j][k] = 'NaN'
        data=numpy.array(datafloat, dtype=float)
        dictnamelist = {}
        for j,varname in enumerate(headerow):
            dictnamelist[varname]=data[:,j]
        Dict[fil] = dictnamelist

    return Dict
    
    
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
                                
                """
            
            print helptext
            
            sys.exit(2)      
    
    for items in args:
        items = items.lower()
    if len(args)>=1:
        print "FAILURE\nNo script arguments required. Define the concerned runs within the script!\n"
        sys.exit(2)

    # Current directory
    currentdir = os.getcwd()

    folderz = ['Saved.SENS_beta=0.20_thetam0_cc/']
    var1 = 'cc'
    var2 = 'thetam0'

    # close all opened figures in python
    close('all')
    
    # set the default parameters of the plotting routines
    rc('xtick', labelsize=18)
    rc('ytick', labelsize=18)
    rc('axes', labelsize=18)
    rc('legend', fontsize=14)
    rc('grid', linewidth=0.)
 
    
    # open a new figure
    fig = figure()
    fig.subplots_adjust(0.15,0.12,0.97,0.85,0.3,0.16)
    
    subfolderlist = []
    for i in range(0,101): # cloud cover indices
        for j in range (100,101): # temperature indices
	
            subfoldername = 'SENS_%s%03d%s%03d'%(var1[0:2],(i+1),var2[0:2],(j+1))
	    print subfoldername
            f=open_output(os.path.join(currentdir,folderz[0],subfoldername),['output_sca'])
            plot(f['output_sca']['UTC(hours)'],1000./f['output_sca']['rs(s.m-1)'])
    title('Parameter space %s x %s\n' r'$\theta _{m,0}$=295K'%(var1,r'$\theta _{m,0}$'), fontsize=18)
#    title('Parameter space %s x %s\ncc=30%%'%(var1,r'$\theta _{m,0}$'), fontsize=18)
    xlabel('Time (UTC)')
    ylim(0.,5.)
    ylabel(r'$g_s$ (mm s$^{-1}$)')
    show()
