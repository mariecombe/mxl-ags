#!/usr/bin/env python

import sys
import os
import getopt
import shutil

import numpy
from matplotlib import pyplot,rc, rcdefaults
from math import exp,log, log10, cos, sin, acos, asin
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline


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

    # Working directories
    currentdir = os.getcwd()
    storagedir = '/Storage/MXL/Marie/MXL-Ags_Sensitivity'

    # First close all opened figures in python
    pyplot.close('all')

    # general plot properties
    rcdefaults()

    rc('xtick', labelsize=16)
    rc('xtick.major', size=6)
    rc('ytick', labelsize=16)
    rc('ytick.major', size=6)
##
    rc('contour', negative_linestyle = 'solid')
    rc('lines', linewidth=1)
    rc('axes', labelsize=16, linewidth=0.5)
    rc('legend', frameon=False, handlelength=3., fontsize=16)
    rc('grid', linewidth=0.)


    figure = True
    
    if (figure == True): #if you want to plot figure 2

        x1    = numpy.linspace(0.,1.,101.)
        beta1 = ( 1-numpy.exp(-(0.00001)*x1) )/( 1-exp(-(0.00001)) )
        beta2 = ( 1-numpy.exp(-(15.)*x1) )/( 1-exp(-(15.)) )

        folderlist = ["1DIM_SA_SMI_0.-1.","1DIM_SA_SMI_0.-1._P4=15"]
        Results = oneD_post_process(storagedir,folderlist)
    
        fig = pyplot.figure('beta', figsize=(12,7))
        fig.subplots_adjust(0.10,0.10,0.97,0.95,0.5,0.3) # left,bottom,right,top, vertical spaces, horizontal spaces

        ax = fig.add_subplot(2,3,1)
        pyplot.axhline(y=1., lw=1., color='black', linestyle='--')
        pyplot.plot(x1,beta1,lw=2.,c='g',label='linear')
        pyplot.plot(x1,beta2,lw=2.,c='r',label='curved')
        ax.set_ylim(0.,1.1)
        ax.set_xlim(0.,1.)
        ax.set_xticks([0.,0.2,0.4,0.6,0.8,1.])
        ax.set_ylabel(r'$\beta$  (applied on $A_g$)')
        ax.set_xlabel('SMI (-)')
        ax.legend(loc='best', ncol=1, fontsize=14)

        ax = fig.add_subplot(2,3,4)
        pyplot.plot(Results["1DIM_SA_SMI_0.-1."]['SMI'],Results["1DIM_SA_SMI_0.-1."]['gs'],lw=2.,c='g')
        pyplot.plot(Results["1DIM_SA_SMI_0.-1._P4=15"]['SMI'],Results["1DIM_SA_SMI_0.-1._P4=15"]['gs'],lw=2.,c='r')
        ax.set_xlim(0.,1)
        ax.set_ylim(0.,22.)
        ax.set_xticks([0.,0.2,0.4,0.6,0.8,1.])
        ax.set_ylabel(r'$g_{l,s}$ @ $12 UTC$'+'\n'+r'(mm s$^{-1}$)')
        ax.set_xlabel('SMI (-)')
    
        ax = fig.add_subplot(2,3,2)
        pyplot.plot(Results["1DIM_SA_SMI_0.-1."]['SMI'],Results["1DIM_SA_SMI_0.-1."]['LE'],lw=2.,c='g')
        pyplot.plot(Results["1DIM_SA_SMI_0.-1._P4=15"]['SMI'],Results["1DIM_SA_SMI_0.-1._P4=15"]['LE'],lw=2.,c='r')
        ax.set_xlim(0.,1)
        ax.set_ylim(0.,430.)
        ax.set_xticks([0.,0.2,0.4,0.6,0.8,1.])
        ax.set_yticks([0.,100.,200.,300.,400.])
        ax.set_ylabel(r'LE @ $12 UTC$'+'\n'+r'(W m$^{-2}$)')
        ax.set_xlabel('SMI (-)')

        ax = fig.add_subplot(2,3,5)
        pyplot.plot(Results["1DIM_SA_SMI_0.-1."]['SMI'],Results["1DIM_SA_SMI_0.-1."]['NPP'],lw=2.,c='g')
        pyplot.plot(Results["1DIM_SA_SMI_0.-1._P4=15"]['SMI'],Results["1DIM_SA_SMI_0.-1._P4=15"]['NPP'],lw=2.,c='r')
        ax.set_xlim(0.,1)
        ax.set_ylim(0.,2.7)
        ax.set_xticks([0.,0.2,0.4,0.6,0.8,1.])
        ax.set_ylabel(r'NPP @ $12 UTC$'+'\n'+r'(mg$_{CO2}$ m$^{-2}$ s$^{-1}$)')
        ax.set_xlabel('SMI (-)')

        ax = fig.add_subplot(2,3,3)
        pyplot.plot(Results["1DIM_SA_SMI_0.-1."]['SMI'],Results["1DIM_SA_SMI_0.-1."]['EF'],lw=2.,c='g')
        pyplot.plot(Results["1DIM_SA_SMI_0.-1._P4=15"]['SMI'],Results["1DIM_SA_SMI_0.-1._P4=15"]['EF'],lw=2.,c='r')
        ax.set_xlim(0.,1)
        ax.set_ylim(0.,0.9)
        ax.set_xticks([0.,0.2,0.4,0.6,0.8,1.])
        ax.set_yticks([0.,0.2,0.4,0.6,0.8])
        ax.set_ylabel(r'EF @ $12 UTC$ (-)')
        ax.set_xlabel('SMI (-)')

        ax = fig.add_subplot(2,3,6)
        pyplot.plot(Results["1DIM_SA_SMI_0.-1."]['SMI'],Results["1DIM_SA_SMI_0.-1."]['WUE'],lw=2.,c='g')
        pyplot.plot(Results["1DIM_SA_SMI_0.-1._P4=15"]['SMI'],Results["1DIM_SA_SMI_0.-1._P4=15"]['WUE'],lw=2.,c='r')
        ax.set_xlim(0.,1)
        ax.set_ylim(0.,6.)
        ax.set_xticks([0.,0.2,0.4,0.6,0.8,1.])
        ax.set_ylabel('daytime average WUEeco\n'+r'(g$_C$ Kg$_{H2O}$$^{-1}$)')
        ax.set_xlabel('SMI (-)')


        pyplot.show()
