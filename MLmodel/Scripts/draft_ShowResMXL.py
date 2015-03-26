#!/usr/bin/env python
# ShowResGECROS.py


"""
Started 10 may 2012
Author: Marie Combe

Plotting script for the MXL res.dat file

"""

# importing the necessary modules

import sys
import os
import getopt
#from csv import *		# reading .csv files
import numpy		# for arrays handling
import matplotlib
from matplotlib import *
from matplotlib.pyplot import *
from math import exp

# defining custom functions
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

#----------------------------------------------------------------------
def plot_MXL(no,x,Vars_list,color,xaxislabel,yaxislabel,titleplot,ann):
# fig = plot_GECROS(j,times,vars,colorr,tlabel,out_vars[item],titl,'False')
    """
Routine to create a plot from GECROS output data
Arguments to input are:
        no:             id number of the simulation
        x:              variable on x axis
        Vars_list:      variables on y axis
        color:          one color assigned to all the y variables
        xaxislabel:     label for x axis
        yaxislabel:     label for y axis
        titleplot:      title of the plot
        ann:            True or False statement, makes an annotation 
                        to point at anthesis time on the plot
        
if annotation is requested, the anthesis time must be known
To know anthesis time, run the search_anth routine
        
        """
    a=len(Vars_list) # a is the number of variables to plot
    plot(x[0], Vars_list[0], c=color, linewidth=2, linestyle='solid',
         label=no)
    xlabel(xaxislabel)
    ylabel(yaxislabel)
    title(titleplot, fontsize=12)
    # putting a legend only if there are more than one variable to plot
    if a>1:
        legend(loc='upper left')	
    grid()
    # Annotation if desired: one vertical dashed line and an arrow with text
    if ann=='True': 
        axvline(x=search_anth()[1], linewidth=1.5, linestyle='dashed', color='black')
        annotate('anthesis starts', xy=(anthesis[1],900), xytext=(200,1100), 
                 arrowprops=dict(facecolor='k', width=1.5, headwidth=6, shrink=0.05))
    # need to annotate at absolute coordinates (left-center?)
    else:
        pass
    return fig

#----------------------------------------------------------------------
def open_ObsHaarweg(inpath,namefile):

    import csv

    inputpath=inpath+namefile
    f=open(inputpath,'rU')
    reader=csv.reader(f, delimiter=',', skipinitialspace=True)
    all=[]
    for row in reader:
        all.append(row)
    headerow=all[0]
    unitsrow=all[1]
    del all[0:2]
    datafloat=[]
    for row in all:
        datafloat.append(map(float,row))
    data=numpy.array(datafloat)
    dictnamelist = {}
    for j,varname in enumerate(headerow):
        dictnamelist[varname]=data[:,j]
    return dictnamelist

#----------------------------------------------------------------------
def CreateDirs(dirname,forceclean=False):
    """
        Create a directory and report success, only create if non-existent or forceclean is used
        
        """
    import shutil
    
    if forceclean:
        try:
            shutil.rmtree(dirname)
        except:
            pass
    
    if not os.path.exists(dirname):
        os.makedirs(dirname)
        msg='Creating new directory %s' % dirname
        #        logging.info(msg)
        print msg
    
    else:
        msg='Using existing directory %s' % dirname
        #        logging.info(msg)
        print msg
    
    return None

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
                Script for ploting the MXL model output (output_dyn, output_sca files). 
                Execute it with the command:
                
                []>./ShowResMXL.py folder titleplot
                folder    :  name of the folder where the output files to use are located.
                titleplot :  title to be writen on all plots (e.g. name of the case?)
                
                The folder will then get a PLOT folder containing the results.
                
                """
            
            print helptext
            
            sys.exit(2)      
    
    for items in args:
        items = items.lower()
    if not os.path.exists(args[0]): # if path to result file does not exist
        print "Does not exist! (%s)" % items
        sys.exit(2)
    if len(args)<1:
        print "FAILURE\nYou forgot to write script arguments!! See help text (call it with -h)\n"
        sys.exit(2)

    # Current directory
    currentdir = os.getcwd()

    folder = args[0]
    if len(args)==2:
      folder2= args[1]
    if len(args)==1:
      folder2 = 'OBS04AUG2007'
   

    # output variables and theirs labels for plots

    out_vars = {'output_dyn':{'zi(m)':"BL height [m]",
                              'thetam(K)':"Mixed layer potential temperature [K]",
                              'we(m.s-1)':"Entrainment velocity [m/s]",
                              'wte(Km.s-1)':"Entrainment heat flux [K.m/s]",
                              'dtheta(K)':"Potential temperature jump at top of BL [K]",
                              'beta':"Heat entrainment ratio (= entrainment flux/subsidence flux)",
                              'dv(m.s-1)':"v-wind jump at BL top [m/s]",
                              'RT(hours)':"Time since start of run [h]",
                              'UTC(hours)':"Time [h UTC]",
                              'vm(m.s-1)':"Mixed-layer v wind speed [m/s]",
                              'um(m.s-1)':"Mixed-layer u wind speed [m/s]",
                              'ws(m.s-1)':"Large scale vertical velocity [m/s]",
                              'du(m.s-1)':"u-wind jump at BL top [m/s]",
                              'wts(Km.s-1)':"Surface heat flux [K.m/s]",
                              'SH(W.m-2)':"Sensible heat flux [W.m-2]",
                              'LE(W.m-2)':"Latent heat flux [W.m-2]",
                              'GR(W.m-2)':"Ground heat flux [W.m-2]",
                              'Qnet(W.m-2)':"Net radiation at surface [W.m-2]"},
                'output_sca':{'zi(m)':"BL height [m]",
                              'betac':"CO2 entrainment ratio (= entrainment flux/subsidence flux)",
                              'wqs':"Surface moisture flux [g.kg-1.m.s-1]",
                              'cm(ppm)':"Mixed-layer CO2 mixing ratio [ppm]",
                              'qm(g.kg-1)':"Mixed-layer moisture [g/kg]",
                              'dc(ppm)':"CO2 mixing ratio jump at BL top [ppm]",
                              'UTC(hours)':"Time [h UTC]",
                              'wce':"Entrainment CO2 flux [ppm.m/s]",
                              'betaq':"Moisture entrainment ratio (= entrainment flux/subsidence flux)",
                              'RT(hours)':"Time since start of run [h]",
                              'wqe':"Entrainment moisture flux [g.kg-1.m.s-1]",
                              'dq(g.kg-1)':"Moisture jump at BL top [g/kg]",
                              'wcs':"Surface CO2 flux [ppm.m/s]",
                              'rs':"Surface resistance [s.m-1]",
                              'Resp':"CO2 release by soil respiration [mg.m-2.s-1]",
                              'An':"Net plant assimilation of CO2 (A-Rd) [mg.m-2.s-1]"}}

    # keep only variables to plot in out_vars
    del out_vars['output_dyn']['UTC(hours)']
    del out_vars['output_dyn']['RT(hours)']
    del out_vars['output_sca']['UTC(hours)']
    del out_vars['output_sca']['RT(hours)']    
    
    # OPENING RESULT FILES :

    csvdir = '/Users/mariecombe/Modeling/MXL/MLmodel/Files_csv/'
    
    # opening observations file only if needed
    if folder=='OBS19JUN2007': 
        OBS = 'Dijkgraaf_19-june-2007'
    elif folder=='OBS04AUG2007': 
        OBS = 'Dijkgraaf_04-august-2007'
    else:
        OBS = ''

    OBS = 'Dijkgraaf_04-august-2007'
    
    if OBS == 'Haarweg_2008': 
        ObsDict = open_ObsHaarweg(csvdir,'nanana.csv')
        ObsDict['Hour(UTC)']= numpy.zeros(len(ObsDict['Hr']))
        for i in range (0, len(ObsDict['Hr'])):
            if ObsDict['Mn'][i]==30:
                ObsDict['Hour(UTC)'][i]=ObsDict['Hr'][i]+0.5
            else:
                ObsDict['Hour(UTC)'][i]=ObsDict['Hr'][i] 
    elif OBS=='Dijkgraaf_19-june-2007':
        ObsDict = open_ObsHaarweg(csvdir,'Dijkgraaf_19-june-2007.csv')
    elif OBS=='Dijkgraaf_04-august-2007':
        ObsDict = open_ObsHaarweg(csvdir,'Dijkgraaf_04-august-2007.csv')
        HDict = open_ObsHaarweg(csvdir,'DeBilt_zi_4-aug-2007.csv')
        

    # Retrieve result foldername (passed as script argument) and create
    # 1 empty dictionnary to store the data:
    folderdictlist = []
    outputfileslist = ['output_dyn','output_sca']
    

    # Read all results files and store all dictionnaries into one dictionnary
    ResDict = open_MXLoutput(currentdir,outputfileslist,folder,folderdictlist)
    
    ResDict2 = open_MXLoutput(currentdir,outputfileslist,folder2,folderdictlist)
    #gecdir = '/Users/mariecombe/Modeling/gecros-mxl/branches/marie-gecrosmxl-branch/MXLmodel/COUPLED_RUNS'
    #ResDict2 = open_output(gecdir,['Mxl_Ags_output_dyn','Mxl_Ags_output_sca'])
    #ResDict2['output_dyn'] = ResDict2.pop('Mxl_Ags_output_dyn')
    #ResDict2['output_sca'] = ResDict2.pop('Mxl_Ags_output_sca')
    
    Mxl_gecros_output = open_output('/Users/mariecombe/Modeling/gecros-mxl/branches/marie-gecrosmxl-branch/MXLmodel/COUPLED_RUNS/DOY_216',['output_dyn','output_sca'])

    # PLOTTING :
   
    # close all opened figures in python
    matplotlib.pyplot.close('all')

    # make figures directory:
    plotfolder = os.path.join(currentdir,folder,'PLOTS')
    CreateDirs(plotfolder,forceclean=False)

#        for i,item in enumerate(out_vars.keys()): # for each output file
#            
#            #subplots_adjust(wspace=0.3, hspace=0.4)
#            for var in out_vars[item].keys(): # for each variable
#                # then plot variables
#                fig = figure()
#                print "Plotting %s" % var
#                #times = []
#                #vars = []
#            
#                #for i in range (0,a-1):
#                #    times.append(ResDict[item]['UTC(hours)'])
#                #    vars.append(ResDict[item][var])
#                colorr = 'black'
#                titl = out_vars[item][var]
#
#                #subplot(nb_row,nb_col,j+1)
#                plot(ResDict[item]['UTC(hours)'], ResDict[item][var], c='black', linewidth=2, linestyle='solid')
#                xlabel('UTC(hours)')
#                ylabel(var)
#                title(out_vars[item][var], fontsize=12)
#            
#                # save the figure
#                #fig.suptitle(out_group[h+1], fontsize=20)
#                fig_titl = var + '.png'
#                savepath = os.path.join(currentdir,folder,'PLOTS',fig_titl)
#                fig.savefig(savepath)

################ Figure 1 #################

    fig1=figure(1, figsize=(8,8)) 
    fig1.subplots_adjust(0.09,0.06,0.91,0.92,0.1,0.3)

    fig1.suptitle('SMI = 0.2, cc= 30%\n'+r'$\gamma_{\theta}$ = 2 K km$^{-1}$ (red line) vs. $\gamma_{\theta}$ = 8 K km$^{-1}$ (black line)', fontsize=12)

    sub1=subplot(421)
    plot(ResDict2['output_dyn']['UTC(hours)'], ResDict2['output_dyn']['Qnet(W.m-2)'],'b-', ls='-', lw=2,color='k')
    plot(ResDict['output_dyn']['UTC(hours)'], ResDict['output_dyn']['Qnet(W.m-2)'],'b-', ls='--',lw=2,color='r')
    axhline(y=0., lw=1, color='black')
#    scatter(ObsDict['UTC(hours)'], ObsDict['Qnet(W.m-2)'], color='k')
    ylabel('Qnet [W.m-2]')
    #xlim(0.,23.)
    ylim(0.,600.)
    xlim(6.,18.)
    grid()
    #matplotlib.pyplot.setp(sub1.get_xticklabels(), visible=False)

    a=max(ResDict['output_dyn']['GR(W.m-2)'])
    b=max(ResDict['output_dyn']['SH(W.m-2)'])
    c=max(ResDict['output_dyn']['LE(W.m-2)'])
    if OBS: e=max(ObsDict['LE(W.m-2)'])
    d=max([a,b,c])
    if OBS: d=max([a,b,c,e])
    Rad=[0.,100.,200.,300.,400.,500.,600.,700.,800.,900.]
    diff = Rad-d
    diff = np.where(diff<0.,10000.,diff)
    idx=diff.argmin()
    Max_GR_SH_LE = Rad[idx]
    
    sub2 = subplot(422)
    sub2.yaxis.tick_right()
    sub2.yaxis.set_label_position('right')
    plot(ResDict2['output_dyn']['UTC(hours)'], ResDict2['output_dyn']['GR(W.m-2)'],'b-',lw=2,ls='-',color='k')
    plot(ResDict['output_dyn']['UTC(hours)'], ResDict['output_dyn']['GR(W.m-2)'],'b-',lw=2,ls='-',color='r')
#    scatter(ObsDict['UTC(hours)'], ObsDict['GR1(W.m-2)'], color='k')
    ylabel('GR [W.m-2]')
    axhline(y=0., lw=1, color='black')
    #yticks(Rad)
    xlim(6.,18.)
    #ylim(0.,Max_GR_SH_LE)
    ylim(0.,100.)
    grid()
    #matplotlib.pyplot.setp(sub2.get_xticklabels(), visible=False)

    sub3=subplot(423)
    plot(ResDict2['output_dyn']['UTC(hours)'], ResDict2['output_dyn']['SH(W.m-2)'],'b-', ls='-',lw=2,color='k')
    plot(ResDict['output_dyn']['UTC(hours)'], ResDict['output_dyn']['SH(W.m-2)'],'b-', ls='-',lw=2,color='r')
#    plot(Mxl_gecros_output['output_dyn']['UTC(hours)'], Mxl_gecros_output['output_dyn']['SH(W.m-2)'],
#         label='MXL-GECROS',c='k',lw=2,linestyle='-')
    #scatter(ObsDict['UTC(hours)'], ObsDict['SH(W.m-2)'], color='red')
#    if (OBS=='Dijkgraaf_04-august-2007'): scatter(ObsDict['UTC(hours)'], ObsDict['correctedSH(W.m-2)'], color='k')
    axhline(y=0., lw=1, color='black')
#    yticks(Rad)
    xlim(6.,18.)
    #ylim(0.,Max_GR_SH_LE)
    ylim(0.,400.)
    ylabel('SH [W.m-2]')
    grid()
    #matplotlib.pyplot.setp(sub3.get_xticklabels(), visible=False)

    sub4 = subplot(424)
    sub4.yaxis.tick_right()
    sub4.yaxis.set_label_position('right')
    plot(ResDict2['output_dyn']['UTC(hours)'], ResDict2['output_dyn']['LE(W.m-2)'],'b-',ls='-',lw=2,color='k')
    plot(ResDict['output_dyn']['UTC(hours)'], ResDict['output_dyn']['LE(W.m-2)'],'b-',ls='-',lw=2,color='r')
#    plot(Mxl_gecros_output['output_dyn']['UTC(hours)'], Mxl_gecros_output['output_dyn']['LE(W.m-2)'],
#         label='MXL-GECROS',c='k',lw=2,linestyle='-')
#    scatter(ObsDict['UTC(hours)'], ObsDict['LE(W.m-2)'], color='red')
#    if (OBS=='Dijkgraaf_04-august-2007'): scatter(ObsDict['UTC(hours)'], ObsDict['correctedLE(W.m-2)'], color='k')
    axhline(y=0., lw=1, color='black')
    yticks(Rad)
    xlim(6.,18.)
    #ylim(0.,Max_GR_SH_LE)
    ylim(0.,400.)
    ylabel('LE [W.m-2]')
    grid()
    #matplotlib.pyplot.setp(sub4.get_xticklabels(), visible=False)

    sub5=subplot(425)
#    plot(ResDict['output_dyn']['UTC(hours)'], ResDict['output_dyn']['thetam(K)'],linestyle='-', lw=2,color='black')
    plot(ResDict2['output_dyn']['UTC(hours)'], ResDict2['output_dyn']['T2m(K)'],ls='-',lw=2,color='k')
    plot(ResDict['output_dyn']['UTC(hours)'], ResDict['output_dyn']['T2m(K)'],ls='-',lw=2,color='r')
    #plot(ResDict['output_dyn']['UTC(hours)'], ResDict['output_dyn']['Tskin(K)'],'b-', linestyle=':',lw=2,color='black')
#    scatter(ObsDict['UTC(hours)'], ObsDict['thetam(K)'], color='k')
    #xlabel('time UTC [h]')
    ylabel('theta 2m (K)')
    ylim(285., 305.)
    xlim(6.,18.)
    grid()
    #matplotlib.pyplot.setp(sub5.get_xticklabels(), visible=False)

    sub6= subplot(426)
    sub6.yaxis.tick_right()
    sub6.yaxis.set_label_position('right')
#    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['qm(g.kg-1)'],ls='-',lw=2,color='k')
    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['q2m(g.kg-1)'],ls='-',lw=2,color='k')
#    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['qm(g.kg-1)'],ls='--',lw=2,color='r',label='qm')
    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['q2m(g.kg-1)'],ls='-',lw=2,color='r',label='q2m')
    #plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['qsatTs(g.kg-1)'],'b-', linestyle='--',lw=2,color='black')
#    scatter(ObsDict['UTC(hours)'], ObsDict['qm(g.kg-1)'], color='k')
    ylim([6.0,11.0])
#    legend(loc='lower right', prop={'size':12})
#    yticks([8.0,9.0,10.0,11.0,12.0])
    ylabel('q2m (g.kg-1)')
    xlim(6.,18.)
    grid()
    #matplotlib.pyplot.setp(sub6.get_xticklabels(), visible=False)

    sub7=subplot(427)
#    plot(ResDict['output_dyn']['UTC(hours)'], ResDict['output_dyn']['lcl(m)'],linestyle='--', label='lcl', lw=2,color='black')
    plot(ResDict2['output_dyn']['UTC(hours)'], ResDict2['output_dyn']['zi(m)'],'b-', ls='-',lw=2,color='k')
    plot(ResDict['output_dyn']['UTC(hours)'], ResDict['output_dyn']['zi(m)'],'b-', label='h', ls='-',lw=2,color='r')
#    plot(ResDict2['output_dyn']['UTC(hours)'], ResDict2['output_dyn']['lcl(m)'],linestyle='--', lw=2,color='blue')
#    scatter(HDict['UTC(hours)'], HDict['zi(m)'], color='k')
    ylabel('z[m]')
    ylim(0.,3500.)
    xlim(6.,18.)
    grid()
#    legend(loc='lower right', prop={'size':12})
    xlabel('time UTC [h]')

    sub8=subplot(428)
    sub8.yaxis.tick_right()
    sub8.yaxis.set_label_position('right')
    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['cm(ppm)'], 'b-', ls='-',lw=2,color='k')
    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['cm(ppm)'], 'b-', ls='-',lw=2,color='r')
#    scatter(ObsDict['UTC(hours)'], ObsDict['cm(ppm)'], color='k')
    xlim(6.,18.)
    ylim(330.,430.)
    ylabel('[CO2] [ppm]')
    xlabel('time UTC [h]')
    grid()

#    fig1bis= figure(100, figsize=(8,4)) 
#    fig1bis.subplots_adjust(0.1,0.30,0.91,0.97,0.1,0.2)
#    seq_dash1 = [6,3]
#        
#    sub1=fig1bis.add_subplot(121)
#    line1 = sub1.plot(ResDict['output_dyn']['UTC(hours)'], ResDict['output_dyn']['SH(W.m-2)'],linestyle='-',lw=2,color='k',label=r'MXL-Ag$_{\mathrm{s}}$')
#    line  = sub1.plot(Mxl_gecros_output['output_dyn']['UTC(hours)'], Mxl_gecros_output['output_dyn']['SH(W.m-2)'],
#            label='MXL-GECROS',c='k',lw=2,linestyle='-')
#    line[0].set_dashes(seq_dash1)
#    scatter1 = sub1.scatter(ObsDict['UTC(hours)'], ObsDict['SH(W.m-2)'], color='r',label='OBS')
#    scatter2 = sub1.scatter(ObsDict['UTC(hours)'], ObsDict['correctedSH(W.m-2)'], color='green', label=r'$\beta$ correction')
#    scatter3 = sub1.scatter(ObsDict['UTC(hours)'], ObsDict['corrected2SH(W.m-2)'], color='b', label='EF correction')
#    axhline(y=0., lw=1, color='black')
#    yticks([-50.,0.,50.,100.,150.,200.,250.])
#    xlim(6.,18.)
#    ylim(-50.,250.)
#    ylabel(r'SH (W m$^{-2}$)')
#    xlabel('time UTC (h)')
#    legend = sub1.legend([scatter1,line[0], scatter2,line1[0],scatter3], 
#             ['OBS','MXL-GECROS',r'$\beta$ correction',r'MXL-A-gs','EF correction'],
#             fontsize=10, loc='upper right',ncol=3,bbox_to_anchor = (1.75, -0.22))
#
#
#    sub2=fig1bis.add_subplot(122)
#    sub2.yaxis.tick_right()
#    sub2.yaxis.set_label_position('right')
#    plot(ResDict['output_dyn']['UTC(hours)'], ResDict['output_dyn']['LE(W.m-2)'], lw=2,linestyle='-',color='k',label='LE MXL-AGS')
#    line=sub2.plot(Mxl_gecros_output['output_dyn']['UTC(hours)'], Mxl_gecros_output['output_dyn']['LE(W.m-2)'],
#         label='LE MXL-GECROS',c='k',lw=2,linestyle='-')
#    line[0].set_dashes(seq_dash1)
#    scatter(ObsDict['UTC(hours)'], ObsDict['LE(W.m-2)'], color='r')
#    scatter(ObsDict['UTC(hours)'], ObsDict['correctedLEbis(W.m-2)'], color='green')
#    scatter(ObsDict['UTC(hours)'], ObsDict['corrected2LE(W.m-2)'], color='b')
#    axhline(y=0., lw=1, color='black')
#    yticks([0.,100.,200.,300.,400.])
#    xlim(6.,18.)
#    #ylim(0.,Max_GR_SH_LE)
#    ylim(0.,450.)
#    ylabel(r'LE (W m$^{-2}$)')
#    xlabel('time UTC (h)')
#    
#    fig1ter = figure (101, figsize=(5,4))
#    fig1ter.subplots_adjust(0.13,0.25,0.95,0.97,0.1,0.2)
#    seq_dash1 = [6,3]
#        
#    sub1=fig1ter.add_subplot(111)
#    line1 = sub1.plot(ResDict2['output_dyn']['UTC(hours)'], ResDict2['output_dyn']['GR(W.m-2)'],linestyle='-',lw=2,color='k',label=r'MXL-Ag$_{\mathrm{s}}$')
#    line  = sub1.plot(Mxl_gecros_output['output_dyn']['UTC(hours)'], Mxl_gecros_output['output_dyn']['GR(W.m-2)'],
#            label='MXL-GECROS',c='k',lw=2,linestyle='-')
#    line[0].set_dashes(seq_dash1)
#    scatter1 = sub1.scatter(ObsDict['UTC(hours)'], ObsDict['GR1(W.m-2)'], color='r',label='OBS')
#    #scatter2 = sub1.scatter(ObsDict['UTC(hours)'], ObsDict['correctedSH(W.m-2)'], color='green', label=r'$\beta$ correction')
#    #scatter3 = sub1.scatter(ObsDict['UTC(hours)'], ObsDict['corrected2SH(W.m-2)'], color='b', label='EF correction')
#    axhline(y=0., lw=1, color='black')
#    yticks([0.,10.,20.,30.,40.,50.,60.,70.])
#    xlim(6.,18.)
#    ylim(0.,75.)
#    ylabel(r'G (W m$^{-2}$)')
#    xlabel('time UTC (h)')
#    legend = sub1.legend([scatter1,line[0],line1[0]], 
#             ['OBS','MXL-GECROS',r'MXL-A-gs'],
#             fontsize=10, loc='upper right',ncol=3,bbox_to_anchor = (0.95, -0.22))

    
################ Figure 2 #################

#    fig2= figure(2, figsize=(8,6)) 
#    fig2.subplots_adjust(0.1,0.08,0.91,0.92,0.1,0.2)
#
#    fig2.suptitle(folder, fontsize=12)
#
#    sub1=subplot(221)
#    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['rs(s.m-1)'],'b-', lw=2,color='black')
#    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['rs(s.m-1)'],'b-', lw=2,color='blue')
#    scatter(ObsDict['UTC(hours)'], ObsDict['rs(s.m-1)'], color='red')
#    ylabel('rs [s.m-1]')
#    ylim(0.,1500.)
#    xlim(6.,18.)
#    grid()
#    #matplotlib.pyplot.setp(sub1.get_xticklabels(), visible=False)
#
#    sub2=subplot(222)
#    sub2.yaxis.tick_right()
#    sub2.yaxis.set_label_position('right')
#    axhline(y=0., lw=1, color='black')
#    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['An(mgCO2.m-2.s-1)'], linestyle='--', c='black', lw=2,label='An')
#    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['Resp(mgCO2.m-2.s-1)'], linestyle=':', c='black', lw=2,label='Resp')
#    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['An(mgCO2.m-2.s-1)'], linestyle='--', c='blue', lw=2)
#    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['Resp(mgCO2.m-2.s-1)'], linestyle=':', c='blue', lw=2)
#    sub2.legend(loc='center', prop={'size':12})
#    ylabel('Fc [mgCO2.m-2.s-1]')
#    xlim(6.,18.)
#    ylim(-2.5,0.5)
#    grid()
#   # matplotlib.pyplot.setp(sub2.get_xticklabels(), visible=False)
#	
#    NEE = numpy.zeros(len(ResDict['output_sca']['An(mgCO2.m-2.s-1)']))
#    for i in range(0, len(ResDict['output_sca']['An(mgCO2.m-2.s-1)'])):
#        NEE[i] = float(ResDict['output_sca']['An(mgCO2.m-2.s-1)'][i]) + float(ResDict['output_sca']['Resp(mgCO2.m-2.s-1)'][i])
#
#    if OBS:
#        NEE2 = numpy.zeros(len(ResDict2['output_sca']['An(mgCO2.m-2.s-1)']))
#        for i in range(0, len(ResDict2['output_sca']['An(mgCO2.m-2.s-1)'])):
#            NEE2[i] = float(ResDict2['output_sca']['An(mgCO2.m-2.s-1)'][i]) + float(ResDict2['output_sca']['Resp(mgCO2.m-2.s-1)'][i])
#
#    sub3=subplot(223)
#    axhline(y=0., lw=1, color='black')
#    plot(ResDict['output_sca']['UTC(hours)'], NEE, 'b-', lw=2,label='NEE',color='black')
#    plot(ResDict2['output_sca']['UTC(hours)'], NEE2, 'b-', lw=2,label='NEE',color='blue')
#    scatter(ObsDict['UTC(hours)'], ObsDict['NEE(mgCO2.m-2.s-1)'], color='red')
#    ylabel('NEE [mgCO2.m-2.s-1]')
#    xlabel('time UTC [h]')
#    xlim(6.,18.)
#    ylim(-3.,1.)
#    grid()
#
#    sub4=subplot(224)
#    sub4.yaxis.tick_right()
#    sub4.yaxis.set_label_position('right')
#    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['cm(ppm)'], 'b-', lw=2,color='black')
#    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['cm(ppm)'], 'b-', lw=2,color='blue')
#    scatter(ObsDict['UTC(hours)'], ObsDict['cm(ppm)'], color='red')
#    xlim(6.,18.)
#    ylim(330.,450.)
#    ylabel('[CO2] [ppm]')
#    xlabel('time UTC [h]')
#    grid()

################ Figure 3 #################

    fig3 = figure()
    #plot(ResDict['output_dyn']['UTC(hours)'], (-0.001* 12./44.* NEE/(ResDict['output_dyn']['LE(W.m-2)']/2.5e6)), 'b-', lw=2,color='black')
    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['wg(cm3.cm-3)'], 'b-', ls='--',lw=2,color='r')
    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['wg(cm3.cm-3)'], 'b-',ls='-',lw=2,color='k')
    ylabel(r'Soil water content $\theta$ [cm$^{3}$ cm$^{-3}$]')
    xlabel('time UTC [h]')
    title(folder, fontsize=12)

################# Figure 3 #################
#
#    fig5 = figure(5, figsize=(6,5))
#    fig5.subplots_adjust(0.2,0.15,0.97,0.9,0.3,0.16)
#    fig5.suptitle(folder, fontsize=12)
#    
#    sub1 = subplot(211)
#    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['thetasurf(K)'], lw=2,color='black')
#    ylabel(r'$\theta _{surf}$ [K]')
#    xlabel('time UTC [h]')
#    xlim(6.,18.)
#    
#    sub1 = subplot(212)
#    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['qsurf(g.kg-1)'], lw=2,color='black')
#    ylabel(r'$q_{surf}$ [g.kg-1]')
#    xlabel('time UTC [h]')
#    xlim(6.,18.)
#    xlabel('time UTC [h]')
#    
    #fig3 = figure()
    #!y=
    #!plot(ResDict['output_dyn']['UTC(hours)'], (-0.001* 12./44.* NEE/(ResDict['output_dyn']['LE(W.m-2)']/2.5e6)), 'b-', lw=2,color='black')

################ Figure 4 #################

#    fig4bis = figure(40, figsize=(7,8))
#    fig4bis.subplots_adjust(0.17,0.08,0.97,0.92,0.5,0.2)
#    fig4bis.suptitle('high subsidence (blue) vs. control case (red)', fontsize=12)
#
#    sub0 = subplot(321)
#    axhline(y=0., lw=1, color='black')
#    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['Ds'], ls='-', lw=2, color='b')
#    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['Ds'], ls='--', lw=2, color='red')
##    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['Dstar'], ls='-', lw=2, color='b',label=r'$D_*$')
##    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['Dstar'], ls='--', lw=2, color='red')
#    ylabel(r'$D_s$ [kPa]')
##    sub0.legend(loc='upper left', fontsize=12.)
#    xlim(6.,18.)
#
#    sub2 = subplot(322)
#    plot(-1,-1,ls='-',lw=1,color='k',label=r'$c_{a}$')
#    plot(-1,-1,ls='-',lw=2,color='k',label=r'$c_{i}$')
#    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['cm(ppm)'], ls='-', lw=1,color='b')
#    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['ci(ppm)'], ls='-', lw=2,color='b')
#    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['cm(ppm)'], ls='--', lw=1,color='r')
#    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['ci(ppm)'], ls='--', lw=2,color='r')
#    ylabel(r'CO$_2$ concentrations [ppm]')
#    ylim(0.,422.)
#    sub2.legend(loc='best')
#    xlim(6.,18.)
#
#    sub2 = subplot(323)
#    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['ci(ppm)']-ResDict['output_sca']['cm(ppm)'], ls='-', lw=2,color='b')
#    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['ci(ppm)']-ResDict2['output_sca']['cm(ppm)'], ls='--', lw=2,color='r')
#    ylabel(r'$c_i$-$c_a$ [ppm]')
#    #ylim(0.,422.)
##    sub2.legend(loc='best')
#    xlim(6.,18.)
#
#    sub10 = subplot(3,2,4)
#    axhline(y=0., lw=1, color='black')
##    plot(ResDict['output_sca']['UTC(hours)'], 1000.*ResDict['output_sca']['gcco2(s.m-1)'], ls='-', lw=2,color='b',label=r'$g_{c,co2}$')
#    plot(ResDict['output_sca']['UTC(hours)'], 1000.*1.6*ResDict['output_sca']['gcco2(s.m-1)'], ls='-', lw=2,color='b')
##    plot(ResDict2['output_sca']['UTC(hours)'], 1000.*ResDict2['output_sca']['gcco2(s.m-1)'], ls='--', lw=2,color='r')
#    plot(ResDict2['output_sca']['UTC(hours)'], 1000.*1.6*ResDict2['output_sca']['gcco2(s.m-1)'], ls='--', lw=2,color='r')
##    sub10.legend(loc='best')
#    ylabel('Canopy conductance for\n'+ r'CO$_2$ transfer [mm s$^{-1}$]')
#    xlabel('time UTC [h]')
#    xlim(6.,18.)
#
#    sub11 = subplot(3,2,5)
##    axhline(y=0., lw=1, color='black')
#    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['An(mgCO2.m-2.s-1)'], lw=2,color='b', label='$A_{net}$')
#    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['An(mgCO2.m-2.s-1)'], ls='--', lw=2,color='r')
#    ylabel(r'Net canopy assimilation'+'\n[mg m$^{-2}$ s$^{-1}$]')
#    xlabel('time UTC [h]')
#    xlim(6.,18.)

#################

#    fig4 = figure(4, figsize=(17,12))
#    fig4.subplots_adjust(0.05,0.08,0.97,0.92,0.3,0.1)
#    fig4.suptitle(folder, fontsize=12)
#
#    sub0 = subplot(341)
#    axhline(y=0., lw=1, color='black')
#    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['Ds'], ls='-', lw=2, color='b',label=r'$D_s$')
#    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['Ds'], ls='--', lw=2, color='red')
#    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['Dstar'], ls='-', lw=2, color='b',label=r'$D_*$')
#    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['Dstar'], ls='--', lw=2, color='red')
#    ylabel(r'VPD [kPa]')
#    sub0.legend(loc='upper left', fontsize=12.)
#    xlim(6.,18.)
#
#    sub1 = subplot(342)
#    axhline(y=0., lw=1, color='black')
#    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['Cfrac'], 'b-', lw=2,color='black')
#    ylabel(r'$f$')
#    xlim(6.,18.)
#
#    sub2 = subplot(343)
#    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['cm(ppm)'], ls='-', lw=2,color='b',label=r'$C_{ext}$')
#    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['ci(ppm)'], ls='-', lw=2,color='b', label=r'$C_i$')
#    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['CO2comp(ppm)'], ls='-', lw=2,color='b',label=r'$\Gamma$')
#    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['cm(ppm)'], ls='--', lw=2,color='r')
#    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['ci(ppm)'], ls='--', lw=2,color='r')
#    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['CO2comp(ppm)'], ls='--', lw=2,color='r')
#    ylabel(r'CO$_2$ concentrations [ppm]')
#    ylim(0.,422.)
#    sub2.legend(loc='best')
#    xlim(6.,18.)
#
#    sub3 = subplot(344)
#    axhline(y=0., lw=1, color='black')
#    plot(ResDict['output_sca']['UTC(hours)'], 1000.*ResDict['output_sca']['gmesophyl'], 'b-', lw=2,color='black', label=r'$g_{m}$')
#    ylabel('Mesophyll conductance [mm s$^{-1}$]')
#    sub3.legend(loc='best')
#    xlim(6.,18.)
#
#    sub4 = subplot(345)
#    axhline(y=0., lw=1, color='black')
#    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['Ammax'], 'b-', lw=2,color='black',label=r'$A_{m,max}$')
#    ylabel('Maximum leaf photosynthetic\nrate ' r'[mg m$^{-2}$ s$^{-1}$]')
#    sub4.legend(loc='best')
#    xlim(6.,18.)
#
#    sub5 = subplot(346)
#    axhline(y=0., lw=1, color='black')
#    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['Am(mgCO2.m-2.s-1)'], lw=2,color='black',label=r'$A_{m}$')
#    plot(ResDict['output_sca']['UTC(hours)'], (1./9.)*ResDict['output_sca']['Am(mgCO2.m-2.s-1)'], linestyle='--', lw=2,color='black',label=r'$R_{dark}$')
#    ylabel('Leaf photosynthetic rate\nand respiration ' r'[mg m$^{-2}$ s$^{-1}$]')
#    sub5.legend(loc='best')
#    xlim(6.,18.)
#
#    sub6 = subplot(347)
#    axhline(y=0., lw=1, color='black')
#    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['scaleAg_star(-)'], 'b-', lw=2,color='black')
#    ylabel(r'Upscaling factor for $A_m$ + $R_{dark}$')
#    xlim(6.,18.)
#
#
#    sub7 = subplot(348)
#    axhline(y=0., lw=1, color='black')
#    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['Ag_star(mgCO2.m-2.s-1)'], lw=2,color='black', label='$A_g$*')
#    ylabel('Unstressed gross canopy\nassimilation ' r'[mg m$^{-2}$ s$^{-1}$]')
#    sub7.legend(loc='best')
#    xlim(6.,18.)
#
#    smi = (ResDict['output_sca']['wg(cm3.cm-3)']-0.06)/(0.15-0.06)
#    beta = [0.]*len(smi)
#    for i,val in enumerate(smi):
#        if (val <= 0.):
#            beta[i] = 1.0e-3
#        elif (val >= 1.0):
#            beta[i] = 1.0
#        else:
#            beta[i] = ( 1.- exp(-0.000001*val) )/( 1.- exp(-0.000001) )
#
#    
#    sub8 = subplot(349)
#    plot(ResDict['output_sca']['UTC(hours)'], beta, 'b-', lw=2,color='black')
#    ylabel(r'Water stress $\beta$ [-]')
#    xlabel('time UTC [h]')
#    xlim(6.,18.)
#
#    sub9 = subplot(3,4,10)
#    axhline(y=0., lw=1, color='black')
#    plot(ResDict['output_sca']['UTC(hours)'], beta*ResDict['output_sca']['Ag_star(mgCO2.m-2.s-1)'], lw=2,color='black', label='$A_g$')
#    ylabel('Water-stressed gross canopy\nassimilation ' r'[mg m$^{-2}$ s$^{-1}$]')
#    xlim(6.,18.)
#    sub9.legend(loc='best')
#    xlabel('time UTC [h]')    
#
#    sub10 = subplot(3,4,11)
#    axhline(y=0., lw=1, color='black')
#    plot(ResDict['output_sca']['UTC(hours)'], 1000.*ResDict['output_sca']['gcco2(s.m-1)'], ls='-', lw=2,color='b',label=r'$g_{c,co2}$')
#    plot(ResDict['output_sca']['UTC(hours)'], 1000.*1.6*ResDict['output_sca']['gcco2(s.m-1)'], ls='-', lw=2,color='b',label=r'$g_{c,w}$')
#    plot(ResDict2['output_sca']['UTC(hours)'], 1000.*ResDict2['output_sca']['gcco2(s.m-1)'], ls='--', lw=2,color='r')
#    plot(ResDict2['output_sca']['UTC(hours)'], 1000.*1.6*ResDict2['output_sca']['gcco2(s.m-1)'], ls='--', lw=2,color='r')
#    sub10.legend(loc='best')
#    ylabel(r'Canopy conductances [mm s$^{-1}$]')
#    xlabel('time UTC [h]')
#    xlim(6.,18.)
#
#    sub11 = subplot(3,4,12)
#    axhline(y=0., lw=1, color='black')
#    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['An(mgCO2.m-2.s-1)'], lw=2,color='b', label='$A_{net}$')
#    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['An(mgCO2.m-2.s-1)'], ls='--', lw=2,color='r')
#    ylabel(r'Net canopy assimilation [mg m$^{-2}$ s$^{-1}$]')
#    xlabel('time UTC [h]')
#    xlim(6.,18.)

################ Figure 6 #################

    fig6 = figure(6, figsize=(9,7))
    fig6.subplots_adjust(0.15,0.08,0.85,0.95,0.0,0.)
    #fig6.suptitle(r'$\beta$ = 0.2, CC = 0%, and thetam0 = 295K', fontsize=12)
    
    
    wte_h = -ResDict['output_dyn']['wte(Km.s-1)']/ResDict['output_dyn']['zi(m)']
    wts_h = ResDict['output_dyn']['wts(Km.s-1)']/ResDict['output_dyn']['zi(m)']
    advt = numpy.array([0.0003]*238+[0.]*479)

    wte_h2 = -ResDict2['output_dyn']['wte(Km.s-1)']/ResDict2['output_dyn']['zi(m)']
    wts_h2 = ResDict2['output_dyn']['wts(Km.s-1)']/ResDict2['output_dyn']['zi(m)']
    advt2 = numpy.array([0.0003]*238+[0.]*479)

    wte_int = sum(wte_h) * (ResDict['output_dyn']['UTC(hours)'][1]-ResDict['output_dyn']['UTC(hours)'][0])*3600.
    wts_int = sum(wts_h) * (ResDict['output_dyn']['UTC(hours)'][1]-ResDict['output_dyn']['UTC(hours)'][0])*3600.
    advt_int = sum(advt) * (ResDict['output_dyn']['UTC(hours)'][1]-ResDict['output_dyn']['UTC(hours)'][0])*3600.
    #print '  t: entr: %6.1f      K, surf: %6.1f      K, adv: %6.1f      K, total: %6.1f      K // %6.1f      K'%( wte_int, wts_int,advt_int, wte_int + wts_int + advt_int, ResDict['output_dyn']['thetam(K)'][716]-ResDict['output_dyn']['thetam(K)'][0])

    ax = fig6.add_subplot(2,2,1)
    axhline(0.,lw=1,color='k')
    plot(ResDict2['output_dyn']['UTC(hours)'], -10.*ResDict2['output_dyn']['wte(Km.s-1)'], lw=2,color='black',label='entrainment')
    plot(ResDict2['output_dyn']['UTC(hours)'], 10.*ResDict2['output_dyn']['wts(Km.s-1)'], lw=2,color='black', linestyle='--',label='surface')
    plot(ResDict2['output_dyn']['UTC(hours)'], 10.*advt2*ResDict2['output_dyn']['zi(m)'], lw=2,color='black', linestyle=':',label='advection')
    plot(ResDict['output_dyn']['UTC(hours)'], -10.*ResDict['output_dyn']['wte(Km.s-1)'], lw=2,color='red')
    plot(ResDict['output_dyn']['UTC(hours)'], 10.*ResDict['output_dyn']['wts(Km.s-1)'], lw=2,color='red', linestyle='--')
    plot(ResDict['output_dyn']['UTC(hours)'], 10.*advt*ResDict['output_dyn']['zi(m)'], lw=2,color='red', linestyle=':')
    #ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    setp(ax.get_xticklabels(), visible=False)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    ax.xaxis.set_ticks_position('both')
    for tick in ax.get_xaxis().get_major_ticks():
        tick.set_pad(2.)
        tick.label1 = tick._get_text1()
    ylabel('heat fluxes\n' r'[1 e$^{-1}$ K m s$^{-1}$]')
    #ylim([-0.0001,0.0006])
    ax.legend(loc='best',fontsize=10)
    xlim(6.,18.)
    ylim(-1.,8.)
    xticks([6.,8.,10.,12.,14.,16.])
    yticks([0.,1.,2.,3.,4.,5.,6.,7.])
     
    ax = fig6.add_subplot(2,2,2)
    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['zi(m)'], lw=2,color='red')
    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['zi(m)'], lw=2,color='black')
    #ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    setp(ax.get_xticklabels(), visible=False)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    ax.xaxis.set_ticks_position('both')
    for tick in ax.get_xaxis().get_major_ticks():
        tick.set_pad(2.)
        tick.label1 = tick._get_text1()
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.yaxis.set_ticks_position('both')
    ylabel('ABL height [m]')
    xlim(6.,18.)
    ylim(0.,3500.)
    xticks([8.,10.,12.,14.,16.,18.])
   
    ax = fig6.add_subplot(2,2,3)
    axhline(0.,lw=1,color='k')
    plot(ResDict['output_dyn']['UTC(hours)'], 10000.*wte_h, lw=2,color='red')
    plot(ResDict['output_dyn']['UTC(hours)'], 10000.*wts_h, lw=2,color='red', linestyle='--')
    plot(ResDict['output_dyn']['UTC(hours)'], 10000.*advt, lw=2,color='red', linestyle=':')
    plot(ResDict2['output_dyn']['UTC(hours)'], 10000.*wte_h2, lw=2,color='black',label='entrainment')
    plot(ResDict2['output_dyn']['UTC(hours)'], 10000.*wts_h2, lw=2,color='black', linestyle='--',label='surface')
    plot(ResDict2['output_dyn']['UTC(hours)'], 10000.*advt2, lw=2,color='black', linestyle=':',label='advection')
    #ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ylabel('tendencies\n' r'[1 e$^{-4}$ K s$^{-1}$]')
    ylim([-1.,6.])
    ax.legend(loc='best',fontsize=10)
    xlim(6.,18.)
    ylim(-1.,6.)
    xticks([6.,8.,10.,12.,14.,16.])
    yticks([0.,1.,2.,3.,4.,5.])
    xlabel('time UTC [h]')
    
    ax = fig6.add_subplot(2,2,4)
    plot(ResDict['output_dyn']['UTC(hours)'], ResDict['output_dyn']['thetam(K)'], lw=2,color='red')
    plot(ResDict2['output_dyn']['UTC(hours)'], ResDict2['output_dyn']['thetam(K)'], lw=2,color='black')
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.yaxis.set_ticks_position('both')
    ylabel(r'$\theta _m$ [K]')
    xlabel('time UTC [h]')
    ylim([286.,302.])
    xlim(6.,18.)
    xticks([8.,10.,12.,14.,16.,18.])
    yticks([286.,288.,290.,292.,294.,296.,298.,300.])
    xlabel('time UTC [h]')


################# Figure 7 #################

    fig7 = figure(7, figsize=(9,7))
    fig7.subplots_adjust(0.15,0.08,0.85,0.95,0.0,0.)
    #fig7.suptitle(r'$\beta$ = 0.2, CC = 0%, and thetam0 = 295K', fontsize=12)

    wqe_h = -ResDict['output_sca']['wqe']/ResDict['output_sca']['zi(m)']
    wqs_h = ResDict['output_sca']['wqs']/ResDict['output_sca']['zi(m)']
    advq = numpy.array([0.00035]*88+[0.]*629)

    wqe_h2 = -ResDict2['output_sca']['wqe']/ResDict2['output_sca']['zi(m)']
    wqs_h2 = ResDict2['output_sca']['wqs']/ResDict2['output_sca']['zi(m)']
    advq2 = numpy.array([0.00035]*88+[0.]*629)

    wqe_int = sum(wqe_h) * (ResDict['output_sca']['UTC(hours)'][1]-ResDict['output_sca']['UTC(hours)'][0])*3600.
    wqs_int = sum(wqs_h) * (ResDict['output_sca']['UTC(hours)'][1]-ResDict['output_sca']['UTC(hours)'][0])*3600.
    advq_int = sum(advq) * (ResDict['output_sca']['UTC(hours)'][1]-ResDict['output_sca']['UTC(hours)'][0])*3600.
    #print '  q: entr: %6.1f g.kg-1, surf: %6.1f g.kg-1, adv: %6.1f g.kg-1, total: %6.1f g.kg-1 // %6.1f g.kg-1'%( wqe_int, wqs_int,advq_int, wqe_int + wqs_int + advq_int, ResDict['output_sca']['qm(g.kg-1)'][716]-ResDict['output_sca']['qm(g.kg-1)'][0])

    ax = fig7.add_subplot(2,2,1)
    axhline(0.,lw=1,color='k')
    plot(ResDict['output_sca']['UTC(hours)'], -10.*ResDict['output_sca']['wqe'], lw=2,color='red')
    plot(ResDict['output_sca']['UTC(hours)'], 10.*ResDict['output_sca']['wqs'], lw=2,color='red', linestyle='--')
    plot(ResDict['output_sca']['UTC(hours)'], 10.*advq*ResDict['output_sca']['zi(m)'], lw=2,color='red', linestyle=':')
    plot(ResDict2['output_sca']['UTC(hours)'], -10.*ResDict2['output_sca']['wqe'], lw=2,color='black',label='entrainment')
    plot(ResDict2['output_sca']['UTC(hours)'], 10.*ResDict2['output_sca']['wqs'], lw=2,color='black', linestyle='--',label='surface')
    plot(ResDict2['output_sca']['UTC(hours)'], 10.*advq2*ResDict2['output_sca']['zi(m)'], lw=2,color='black', linestyle=':',label='advection')
    #ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    setp(ax.get_xticklabels(), visible=False)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    ax.xaxis.set_ticks_position('both')
    for tick in ax.get_xaxis().get_major_ticks():
        tick.set_pad(2.)
        tick.label1 = tick._get_text1()
    ylabel('moisture fluxes\n' r'[1 e$^{-1}$ g kg$^{-1}$ m s$^{-1}$]')
    #ylim([-0.0001,0.0006])
    ax.legend(loc='best',fontsize=10)
    xlim(6.,18.)
    ylim(-8.,2.)
    xticks([6.,8.,10.,12.,14.,16.])
    yticks([-6.,-4.,-2.,0.,2.])

    ax = fig7.add_subplot(2,2,2)
    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['zi(m)'], lw=2,color='red')
    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['zi(m)'], lw=2,color='black')
    #ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    setp(ax.get_xticklabels(), visible=False)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    ax.xaxis.set_ticks_position('both')
    for tick in ax.get_xaxis().get_major_ticks():
        tick.set_pad(2.)
        tick.label1 = tick._get_text1()
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.yaxis.set_ticks_position('both')
    ylabel('ABL height [m]')
    xlim(6.,18.)
    ylim(0.,3500.)
    xticks([8.,10.,12.,14.,16.,18.])
    
    ax = fig7.add_subplot(2,2,3)
    axhline(0.,lw=1,color='k')
    plot(ResDict['output_sca']['UTC(hours)'], 1000.*wqe_h, lw=2,color='red')
    plot(ResDict['output_sca']['UTC(hours)'], 1000.*wqs_h, lw=2,color='red', linestyle='--')
    plot(ResDict['output_sca']['UTC(hours)'], 1000.*advq, lw=2,color='red', linestyle=':')
    plot(ResDict2['output_sca']['UTC(hours)'], 1000.*wqe_h2, lw=2,color='black',label='entrainment')
    plot(ResDict2['output_sca']['UTC(hours)'], 1000.*wqs_h2, lw=2,color='black', linestyle='--',label='surface')
    plot(ResDict2['output_sca']['UTC(hours)'], 1000.*advq2, lw=2,color='black', linestyle=':',label='advection')
    #ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ylabel('tendencies\n' r'[1 e$^{-3}$ g kg$^{-1}$ s$^{-1}$]')
    ax.legend(loc='best',fontsize=10)
    xlim(6.,18.)
    ylim(-2.0,0.5)
    xticks([6.,8.,10.,12.,14.,16.])
    yticks([-2.0,-1.5,-1.0,-0.5,0.0])
    xlabel('time UTC [h]')
    
    ax = fig7.add_subplot(2,2,4)
    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['qm(g.kg-1)'], lw=2,color='red')
    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['qm(g.kg-1)'], lw=2,color='black')
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.yaxis.set_ticks_position('both')
    ylabel(r'$q_m$ [g kg$^{-1}$]')
    xlabel('time UTC [h]')
    ylim([6.5,10.5])
    xlim(6.,18.)
    xticks([8.,10.,12.,14.,16.,18.])
    yticks([6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0])
    xlabel('time UTC [h]')

################# Figure 8 #################

    fig8 = figure(8, figsize=(9,7))
    fig8.subplots_adjust(0.15,0.08,0.85,0.87,0.0,0.)
    fig8.suptitle(r'CO$_2$'+' budget at SMI = 0.2, cc= 30%\n'+r'$\gamma_{\theta}$ = 2 K km$^{-1}$ (blue line) vs. $\gamma_{\theta}$ = 8 K km$^{-1}$ (red line)', fontsize=12)
    
    
    wce_h = -ResDict['output_sca']['wce(ppm.m.s-1)']/ResDict['output_sca']['zi(m)']
    wcs_h = ResDict['output_sca']['wcs(ppm.m.s-1)']/ResDict['output_sca']['zi(m)']
    advc = numpy.array([0.]*len(ResDict['output_sca']['wce(ppm.m.s-1)']))

    wce_h2 = -ResDict2['output_sca']['wce(ppm.m.s-1)']/ResDict2['output_sca']['zi(m)']
    wcs_h2 = ResDict2['output_sca']['wcs(ppm.m.s-1)']/ResDict2['output_sca']['zi(m)']
    advc2 = numpy.array([0.]*len(ResDict['output_sca']['wce(ppm.m.s-1)']))

    wce_int = sum(wce_h) * (ResDict['output_sca']['UTC(hours)'][1]-ResDict['output_sca']['UTC(hours)'][0])*3600.
    wcs_int = sum(wcs_h) * (ResDict['output_sca']['UTC(hours)'][1]-ResDict['output_sca']['UTC(hours)'][0])*3600.
    advc_int = sum(advc) * (ResDict['output_sca']['UTC(hours)'][1]-ResDict['output_sca']['UTC(hours)'][0])*3600.
    #print 'CO2: entr: %6.1f    ppm, surf: %6.1f    ppm, adv: %6.1f    ppm, total: %6.1f    ppm // %6.1f    ppm'%(wce_int, wcs_int,advc_int, wce_int + wcs_int + advc_int, ResDict['output_sca']['cm(ppm)'][716]-ResDict['output_sca']['cm(ppm)'][0])
    
    ax = fig8.add_subplot(2,2,1)
    axhline(0.,lw=1,color='k')
    plot(-1.,-1.,lw=2,color='black', label='entrainment')
    plot(-1.,-1.,lw=2,linestyle=':',color='black',label='surface')
    plot(ResDict['output_sca']['UTC(hours)'], -ResDict['output_sca']['wce(ppm.m.s-1)'], lw=2,color='blue')
    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['wcs(ppm.m.s-1)'], lw=2,color='blue', linestyle='--')
    #plot(ResDict['output_sca']['UTC(hours)'], advc*ResDict['output_sca']['zi(m)'], lw=2,color='blue', linestyle=':')
    plot(ResDict2['output_sca']['UTC(hours)'], -ResDict2['output_sca']['wce(ppm.m.s-1)'], lw=2,color='red')
    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['wcs(ppm.m.s-1)'], lw=2,color='red', linestyle='--')
    #plot(ResDict2['output_sca']['UTC(hours)'], advc2*ResDict2['output_sca']['zi(m)'], lw=2,color='red', linestyle=':',label='advection')
    setp(ax.get_xticklabels(), visible=False)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    ax.xaxis.set_ticks_position('both')
    for tick in ax.get_xaxis().get_major_ticks():
        tick.set_pad(2.)
        tick.label1 = tick._get_text1()
    ylabel(r'CO$_2$' ' fluxes\n' r'[ppm m s$^{-1}$]')
    #ylim([-0.0001,0.0006])
    ax.legend(loc='lower right',fontsize=8)
    xlim(6.,18.)
    ylim(-6.,0.5)
    xticks([6.,8.,10.,12.,14.,16.])
    yticks([-4.,-3.,-2.,-1.,0.])
    
    ax = fig8.add_subplot(2,2,2)
    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['zi(m)'], lw=2,color='blue', label=r'low cc and $\theta_0$')
    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['zi(m)'], lw=2,color='red', label=r'high cc and $\theta_0$')
    setp(ax.get_xticklabels(), visible=False)
    ax.legend(loc='best',fontsize=12)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    ax.xaxis.set_ticks_position('both')
    for tick in ax.get_xaxis().get_major_ticks():
        tick.set_pad(2.)
        tick.label1 = tick._get_text1()
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.yaxis.set_ticks_position('both')
    ylabel('ABL height [m]')
    xlim(6.,18.)
    ylim(0.,3500.)
    xticks([8.,10.,12.,14.,16.,18.])
 
    ax = fig8.add_subplot(2,2,3)
    axhline(0.,lw=1,color='k')
    plot(-1.,-1.,lw=2,color='black', label='entrainment')
    plot(-1.,-1.,lw=2,linestyle=':',color='black',label='surface')
    plot(ResDict['output_sca']['UTC(hours)'], 1000.*wce_h, lw=2,color='blue')
    plot(ResDict['output_sca']['UTC(hours)'], 1000.*wcs_h, lw=2,color='blue', linestyle='--')
    #plot(ResDict['output_sca']['UTC(hours)'], 1000.*advc, lw=2,color='blue', linestyle=':')
    plot(ResDict2['output_sca']['UTC(hours)'], 1000.*wce_h2, lw=2,color='red')
    plot(ResDict2['output_sca']['UTC(hours)'], 1000.*wcs_h2, lw=2,color='red', linestyle='--')
    #plot(ResDict2['output_sca']['UTC(hours)'], 1000.*advc2, lw=2,color='red', linestyle=':',label='advection')
    ylabel('tendencies\n' r'[ppb s$^{-1}$]')
    ax.legend(loc='lower right',fontsize=8)
    xlim(6.,18.)
    ylim(-15.,0.)
    xticks([6.,8.,10.,12.,14.,16.])
    yticks([-15.,-12.5,-10.,-7.5,-5.,-2.5,0.])
    xlabel('time UTC [h]')
    
    ax = fig8.add_subplot(2,2,4)
    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['cm(ppm)'], lw=2,color='blue')
    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['cm(ppm)'], lw=2,color='red')
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.yaxis.set_ticks_position('both')
    ylabel(r'$c_m$ [ppm]')
    xlim(6.,18.)
    ylim(340.,420.)
    yticks([340.,350.,360.,370.,380.,390.,400.,410.])
    xticks([8.,10.,12.,14.,16.,18.])
    xlabel('time UTC [h]')


################# Figure 8 #################

    fig9 = figure(9, figsize=(9,7))
    fig9.subplots_adjust(0.15,0.08,0.85,0.95,0.0,0.)
    #fig9.suptitle(r'$\beta$ = 0.2, CC = 0%, and thetam0 = 295K', fontsize=12)
        
    ax = fig9.add_subplot(2,2,1)
    axhline(0.,lw=1,color='k')
    plot(ResDict['output_sca']['UTC(hours)'], -ResDict['output_sca']['wce(ppm.m.s-1)'], lw=2,color='blue')
    #plot(ResDict['output_sca']['UTC(hours)'], advc*ResDict['output_sca']['zi(m)'], lw=2,color='blue', linestyle=':')
    plot(ResDict2['output_sca']['UTC(hours)'], -ResDict2['output_sca']['wce(ppm.m.s-1)'], lw=2,color='red')
    #plot(ResDict2['output_sca']['UTC(hours)'], advc2*ResDict2['output_sca']['zi(m)'], lw=2,color='red', linestyle=':',label='advection')
    setp(ax.get_xticklabels(), visible=False)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    ax.xaxis.set_ticks_position('both')
    for tick in ax.get_xaxis().get_major_ticks():
        tick.set_pad(2.)
        tick.label1 = tick._get_text1()
    ylabel(r'CO$_2$' ' entrainment\n' r'[ppm m s$^{-1}$]')
    #ylim([-0.0001,0.0006])
    xlim(6.,18.)
    ylim(-4.,0.5)
    xticks([6.,8.,10.,12.,14.,16.])
    yticks([-3.,-2.,-1.,0.])
    
    ax = fig9.add_subplot(2,2,2)
    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['zi(m)'], lw=2,color='blue', label=r'low cc and $\theta_0$')
    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['zi(m)'], lw=2,color='red', label=r'high cc and $\theta_0$')
    setp(ax.get_xticklabels(), visible=False)
    ax.legend(loc='best',fontsize=12)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    ax.xaxis.set_ticks_position('both')
    for tick in ax.get_xaxis().get_major_ticks():
        tick.set_pad(2.)
        tick.label1 = tick._get_text1()
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.yaxis.set_ticks_position('both')
    ylabel('ABL height [m]')
    xlim(6.,18.)
    ylim(0.,2000.)
    xticks([8.,10.,12.,14.,16.,18.])
 
    ax = fig9.add_subplot(2,2,3)
    axhline(0.,lw=1,color='k')
#    plot(-1.,-1.,lw=2,color='black', label='entrainment')
#    plot(-1.,-1.,lw=2,linestyle=':',color='black',label='surface')   
    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['dc(ppm)'], lw=2,color='blue')
    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['dc(ppm)'], lw=2,color='red')    
    ylabel(r'CO$_2$ jump'+'\n[ppm]')
    xlim(6.,18.)
    #ylim(-10.,0.)
    xticks([6.,8.,10.,12.,14.,16.])
    #yticks([-10.,-7.5,-5.,-2.5,0.])
    xlabel('time UTC [h]')
    
#    ax = fig9.add_subplot(2,2,4)
#    plot(ResDict['output_sca']['UTC(hours)'], ResDict['output_sca']['cm(ppm)'], lw=2,color='blue')
#    plot(ResDict2['output_sca']['UTC(hours)'], ResDict2['output_sca']['cm(ppm)'], lw=2,color='red')
#    ax.yaxis.tick_right()
#    ax.yaxis.set_label_position("right")
#    ax.yaxis.set_ticks_position('both')
#    ylabel(r'$c_m$ [ppm]')
#    xlim(6.,18.)
#    ylim(340.,420.)
#    yticks([340.,350.,360.,370.,380.,390.,400.,410.])
#    xticks([8.,10.,12.,14.,16.,18.])
#    xlabel('time UTC [h]')



    savepath = os.path.join(currentdir,folder,'PLOTS')
    fig1.savefig(os.path.join(savepath,'figure1.png'),dpi=300)
#    fig1bis.savefig(os.path.join(savepath,'figure1bis.png'),dpi=300)
#    fig1ter.savefig(os.path.join(savepath,'figure1ter.png'),dpi=300)
#    fig2.savefig(os.path.join(savepath,'figure2.png'),dpi=300)
    fig3.savefig(os.path.join(savepath,'figure3.png'),dpi=300)
#    fig4.savefig(os.path.join(savepath,'figure4.png'),dpi=300)
#    fig4bis.savefig(os.path.join(savepath,'fig4bis.png'),dpi=300)
#    fig5.savefig(os.path.join(savepath,'figure5.png'),dpi=300)
    fig6.savefig(os.path.join(savepath,'figure6.png'),dpi=300)    
    fig7.savefig(os.path.join(savepath,'figure7.png'),dpi=300)
    fig8.savefig(os.path.join(savepath,'figure8.png'),dpi=300)
    fig9.savefig(os.path.join(savepath,'figure9.png'),dpi=300)
    # Report success of plotting routine
    print "\nMission accomplished! :D\n"

    #show()
