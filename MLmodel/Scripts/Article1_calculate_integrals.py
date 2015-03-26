#!/usr/bin/env python
# ShowResGECROS.py


"""
Started 29 November 2013
Author: Marie Combe

Plotting script for the MXL res.dat file

"""

# importing the necessary modules

import sys
import os
import getopt
import numpy		# for arrays handling

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
        if (fil.startswith('Mxl_Ags_output_dyn') or fil.endswith('Mxl_Ags_output_dyn') or fil.startswith('Mxl_Ags_output_sca') or fil.endswith('Mxl_Ags_output_sca')): 
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

    # Current directory

    gecrosdir = '/Users/mariecombe/Modeling/gecros-mxl/branches/marie-gecrosmxl-branch/'
    agsdir = '/Users/mariecombe/Modeling/MXL/MLmodel/'

    #AgsDict = open_output(os.path.join(agsdir,'OBS04AUG2007'),['output_dyn','Mxl_Ags_output_sca'])
    AgsDict = open_output(os.path.join(gecrosdir,'MXLmodel/COUPLED_RUNS/'),['Mxl_Ags_output_dyn','Mxl_Ags_output_sca'])
    GecDict = open_output(os.path.join(gecrosdir,'MXLmodel'),['MXL_weather'])
    GecDict2 = open_output(os.path.join(gecrosdir,'MXLmodel/COUPLED_RUNS/'),['Haarweg_weather'])
    ObsDict = open_ObsHaarweg(agsdir,'Dijkgraaf_04-august-2007.csv')
    cases = ['Observations','MXL-AGs','MXL-GECROS']
    
    
    for i,case in enumerate(cases):
    
        SumSwin = 0.
        SumSwnet = 0.
        SumLwnet = 0.
	SumQnet = 0.
	SumLE = 0.
	SumSH = 0.
	SumGR = 0.
	SumNEE = 0.
	n = 0.

	if (i==0):
	    dt = ObsDict['UTC(hours)'][1] - ObsDict['UTC(hours)'][0]
	    for j,val in enumerate(ObsDict['Swin(W.m-2)']):
                if (val > 0. and 8. <= ObsDict['UTC(hours)'][j] <= 18.):
	            SumSwin  = SumSwin   + val
		    SumSwnet = SumSwnet  + ObsDict['SWnet(W.m-2)'][j]
		    SumLwnet = SumLwnet  + ObsDict['LWnet(W.m-2)'][j]
		    SumQnet  = SumQnet   + ObsDict['Qnet(W.m-2)'][j]
		    if (ObsDict['correctedLE(W.m-2)'][j] >0.): 
		    	SumLE    = SumLE     + ObsDict['correctedLE(W.m-2)'][j]
		    SumSH    = SumSH     + ObsDict['correctedSH(W.m-2)'][j]
		    SumGR    = SumGR     + ObsDict['GR1(W.m-2)'][j]
		    SumNEE   = SumNEE    + ObsDict['NEE(mgCO2.m-2.s-1)'][j]/1000.
		    n        = n         + 1.
		elif (ObsDict['UTC(hours)'][j]>18.):
	            SumSwin  = SumSwin   - ObsDict['Swin(W.m-2)'][j-1]
		    SumSwnet = SumSwnet  - ObsDict['SWnet(W.m-2)'][j-1]
		    SumLwnet = SumLwnet  - ObsDict['LWnet(W.m-2)'][j-1]
		    SumQnet  = SumQnet   - ObsDict['Qnet(W.m-2)'][j-1]
		    SumLE    = SumLE     - ObsDict['correctedLE(W.m-2)'][j-1]
		    SumSH    = SumSH     - ObsDict['correctedSH(W.m-2)'][j-1]
		    SumGR    = SumGR     - ObsDict['GR1(W.m-2)'][j-1]
		    SumNEE   = SumNEE    - ObsDict['NEE(mgCO2.m-2.s-1)'][j-1]/1000.
		    break
            SumSwin   = SumSwin * dt * 3600.
            SumSwnet  = SumSwnet * dt * 3600.
            SumLwnet  = SumLwnet * dt * 3600.
	    SumQnet   = SumQnet * dt * 3600.
	    SumLE     = SumLE * dt * 3600.
	    SumSH     = SumSH * dt * 3600.
	    SumGR     = SumGR * dt * 3600.
	    SumNEE    = SumNEE * dt * 3600.#*12./44.
	
	    print '\nObservations from 6am to 6pm:'
	    print '%12s%12s%12s%12s%12s%12s%12s%12s'%('Sum SWin','Sum SWnet','Sum LWnet','Sum Qnet','Sum cor.LE','Sum cor.SH','Sum GR','Sum NEE')
	    print '%12.2f%12.2f%12.2f%12.2f%12.2f%12.2f%12.2f%12.2f'%(SumSwin/1000000., SumSwnet/1000000.,SumLwnet/1000000.,SumQnet/1000000.,SumLE/1000000.,SumSH/1000000.,SumGR/1000000.,SumNEE)

            SumSwin = 0.
	    SumQnet = 0.
	    SumLE = 0.
	    SumSH = 0.
	    SumNEE = 0.
	    n = 0.

	    for j,val in enumerate(ObsDict['Swin(W.m-2)']):
                if (val > 0.):
	            SumSwin  = SumSwin   + val
		    SumQnet  = SumQnet   + ObsDict['Qnet(W.m-2)'][j]
		    if (ObsDict['correctedLE(W.m-2)'][j] >0.): 
		    	SumLE    = SumLE     + ObsDict['correctedLE(W.m-2)'][j]
		    SumSH    = SumSH     + ObsDict['correctedSH(W.m-2)'][j]
		    SumNEE   = SumNEE    + ObsDict['NEE(mgCO2.m-2.s-1)'][j]/1000.
		    n        = n         + 1.
		elif (val <= 0. and ObsDict['UTC(hours)'][j]>12.):
	            SumSwin  = SumSwin   - ObsDict['Swin(W.m-2)'][j-1]
		    SumQnet  = SumQnet   - ObsDict['Qnet(W.m-2)'][j-1]
		    SumLE    = SumLE     - ObsDict['correctedLE(W.m-2)'][j-1]
		    SumSH    = SumSH     - ObsDict['correctedSH(W.m-2)'][j-1]
		    SumNEE   = SumNEE    - ObsDict['NEE(mgCO2.m-2.s-1)'][j-1]/1000.
		    break
            SumSwin   = SumSwin * dt * 3600.
	    SumQnet   = SumQnet * dt * 3600.
	    SumLE     = SumLE * dt * 3600.
	    SumSH     = SumSH * dt * 3600.
	    SumNEE    = SumNEE * dt * 3600.

	    print '\nObservations from sunrise to sunset:     Sum SWin: calcul. = %6.2f, Jacobs = %6.2f'%(SumSwin/1000000.,24.88)
	    print '                                         Sum Qnet: calcul. = %6.2f, Jacobs = %6.2f'%(SumQnet/1000000.,15.45)
	    print '                                         Sum LE:   calcul. = %6.2f, Jacobs = %6.2f'%(SumLE/1000000.,9.55)
	    print '                                         Sum SH:   calcul. = %6.2f, Jacobs = %6.2f'%(SumSH/1000000.,2.63)
	    print '                                         Sum NEE:  calcul. = %6.2f, Jacobs = %6.2f'%(SumNEE,-12.62*44./12.)

            print 'dt =', dt
############
# MXL-AGS
############

        if (i==1):
	    dt = AgsDict['Mxl_Ags_output_dyn']['UTC(hours)'][1] - AgsDict['Mxl_Ags_output_dyn']['UTC(hours)'][0]
	    for j,val in enumerate(AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)']):
                if (val > 0. and 8. <= AgsDict['Mxl_Ags_output_dyn']['UTC(hours)'][j] <= 18.):
	            SumSwin  = SumSwin   + val
		    SumSwnet = SumSwnet  + (AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'][j]-AgsDict['Mxl_Ags_output_dyn']['Swout(W.m-2)'][j])
		    SumLwnet = SumLwnet  + (AgsDict['Mxl_Ags_output_dyn']['Lwin(W.m-2)'][j]-AgsDict['Mxl_Ags_output_dyn']['Lwout(W.m-2)'][j])
		    SumQnet  = SumQnet   + AgsDict['Mxl_Ags_output_dyn']['Qnet(W.m-2)'][j]
		    if (AgsDict['Mxl_Ags_output_dyn']['LE(W.m-2)'][j] >0.): 
		    	SumLE    = SumLE     + AgsDict['Mxl_Ags_output_dyn']['LE(W.m-2)'][j]
		    SumSH    = SumSH     + AgsDict['Mxl_Ags_output_dyn']['SH(W.m-2)'][j]
		    SumGR    = SumGR     + AgsDict['Mxl_Ags_output_dyn']['GR(W.m-2)'][j]
		    SumNEE   = SumNEE    + (AgsDict['Mxl_Ags_output_sca']['An(mgCO2.m-2.s-1)'][j]+AgsDict['Mxl_Ags_output_sca']['Resp(mgCO2.m-2.s-1)'][j])/1000.
		    n        = n         + 1.
		elif (AgsDict['Mxl_Ags_output_dyn']['UTC(hours)'][j] > 18.):
	            SumSwin  = SumSwin   - AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'][j-1]
		    SumSwnet = SumSwnet  - (AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'][j-1]-AgsDict['Mxl_Ags_output_dyn']['Swout(W.m-2)'][j-1])
		    SumLwnet = SumLwnet  - (AgsDict['Mxl_Ags_output_dyn']['Lwin(W.m-2)'][j-1]-AgsDict['Mxl_Ags_output_dyn']['Lwout(W.m-2)'][j-1])
		    SumQnet  = SumQnet   - AgsDict['Mxl_Ags_output_dyn']['Qnet(W.m-2)'][j-1]
		    SumLE    = SumLE     - AgsDict['Mxl_Ags_output_dyn']['LE(W.m-2)'][j-1]
		    SumSH    = SumSH     - AgsDict['Mxl_Ags_output_dyn']['SH(W.m-2)'][j-1]
		    SumGR    = SumGR     - AgsDict['Mxl_Ags_output_dyn']['GR(W.m-2)'][j-1]
		    SumNEE   = SumNEE    - (AgsDict['Mxl_Ags_output_sca']['An(mgCO2.m-2.s-1)'][j-1]+AgsDict['Mxl_Ags_output_sca']['Resp(mgCO2.m-2.s-1)'][j-1])/1000.
		    break
            SumSwin   = SumSwin * dt * 3600.
            SumSwnet  = SumSwnet * dt * 3600.
            SumLwnet  = SumLwnet * dt * 3600.
	    SumQnet   = SumQnet * dt * 3600.
	    SumLE     = SumLE * dt * 3600.
	    SumSH     = SumSH * dt * 3600.
	    SumGR     = SumGR * dt * 3600.
	    SumNEE    = SumNEE * dt * 3600.#*12./44.
	
	    print '##################'
	    print '\nMXL-A-Gs model from 6am to 6pm:'
	    print '%12s%12s%12s%12s%12s%12s%12s%12s'%('Sum SWin','Sum SWnet','Sum LWnet','Sum Qnet','Sum LE','Sum SH','Sum GR','Sum NEE')
	    print '%12.2f%12.2f%12.2f%12.2f%12.2f%12.2f%12.2f%12.2f'%(SumSwin/1000000., SumSwnet/1000000.,SumLwnet/1000000.,SumQnet/1000000.,SumLE/1000000.,SumSH/1000000.,SumGR/1000000.,SumNEE)

            #linear interpolation between sunrise and 6hUTC for Qnet, LE, SH, NEE:

            # On DOY 216, lon=5.38 lat=51.59, cc=22.5%
            # hour start SWin: 4.133333206177
            # hour end SWin:    19.149999618530
            # SWin:
            X2_X1 = AgsDict['Mxl_Ags_output_dyn']['UTC(hours)'][0] - 0.41333E+01
	    Y2_Y1 = (AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'][0] - 0.65782E+00 )
	    Y1    = 0.65782E+00
	    SWinsum = Y1
	    for i in range(X2_X1/dt): 
	        SWinsum = SWinsum + ( Y1 + ((i+1)*dt) * Y2_Y1/X2_X1 )
            SWinsum = SWinsum*dt*3600.

	    X2_X1 = 0.19150E+02 - AgsDict['Mxl_Ags_output_dyn']['UTC(hours)'][len(AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'])-1]
	    Y2_Y1 = (0.58628E+00 - AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'][len(AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'])-1] )
	    Y1    = AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'][len(AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'])-1]
	    SWinsum2 = Y1
	    for i in range(X2_X1/dt): 
	        SWinsum2 = SWinsum2 + ( Y1 + ((i+1)*dt) * Y2_Y1/X2_X1 )
            SWinsum2 = SWinsum2*dt*3600.
	    
	    print '\nMXL-A-gs model from sunrise to sunset:   Sum SWin: %6.2f + %6.2f + %6.2f = %6.2f'%(SumSwin/1000000.,SWinsum/1000000.,SWinsum2/1000000.,(SumSwin+SWinsum+SWinsum2)/1000000.)


            #Qnet: 
            X2_X1 = AgsDict['Mxl_Ags_output_dyn']['UTC(hours)'][0] - 0.41333E+01
	    Y2_Y1 = (AgsDict['Mxl_Ags_output_dyn']['Qnet(W.m-2)'][0] - (-50.) )
	    Y1    = -50.
	    Qnetsum = Y1
	    for i in range(X2_X1/dt): 
	        Qnetsum = Qnetsum + ( Y1 + ((i+1)*dt) * Y2_Y1/X2_X1 )
            Qnetsum = Qnetsum*dt*3600.

	    X2_X1 = 0.19150E+02 - AgsDict['Mxl_Ags_output_dyn']['UTC(hours)'][len(AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'])-1]
	    Y2_Y1 = (-50. - AgsDict['Mxl_Ags_output_dyn']['Qnet(W.m-2)'][len(AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'])-1] )
	    Y1    = AgsDict['Mxl_Ags_output_dyn']['Qnet(W.m-2)'][len(AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'])-1]
	    Qnetsum2 = Y1
	    for i in range(X2_X1/dt): 
	        Qnetsum2 = Qnetsum2 + ( Y1 + ((i+1)*dt) * Y2_Y1/X2_X1 )
            Qnetsum2 = Qnetsum2*dt*3600.
	    
	    print '                                         Sum Qnet: %6.2f + %6.2f + %6.2f = %6.2f'%(SumQnet/1000000.,Qnetsum/1000000.,Qnetsum2/1000000.,(SumQnet+Qnetsum+Qnetsum2)/1000000.)

	    # LE:
	    X2_X1 = AgsDict['Mxl_Ags_output_dyn']['UTC(hours)'][0] - 0.41333E+01
	    Y2_Y1 = AgsDict['Mxl_Ags_output_dyn']['LE(W.m-2)'][0] - 0
	    Y1    = 0
	    LEsum = Y1
	    for i in range(X2_X1/dt): 
	        LEsum = LEsum + ( Y1 + ((i+1)*dt) * Y2_Y1/X2_X1 )
            LEsum = LEsum*dt*3600.

	    X2_X1 = 0.19150E+02 - AgsDict['Mxl_Ags_output_dyn']['UTC(hours)'][len(AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'])-1]
	    Y2_Y1 = 0 - AgsDict['Mxl_Ags_output_dyn']['LE(W.m-2)'][len(AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'])-1]
	    Y1    = AgsDict['Mxl_Ags_output_dyn']['LE(W.m-2)'][len(AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'])-1]
	    LEsum2 = Y1
	    for i in range(X2_X1/dt): 
	        LEsum2 = LEsum2 + ( Y1 + ((i+1)*dt) * Y2_Y1/X2_X1 )
            LEsum2 = LEsum2*dt*3600.
	    
	    print '                                         Sum LE:   %6.2f + %6.2f + %6.2f = %6.2f'%(SumLE/1000000.,LEsum/1000000.,LEsum2/1000000.,(SumLE+LEsum+LEsum2)/1000000.)


	    # SH:
	    X2_X1 = AgsDict['Mxl_Ags_output_dyn']['UTC(hours)'][0] - 0.41333E+01
	    Y2_Y1 = AgsDict['Mxl_Ags_output_dyn']['SH(W.m-2)'][0] - (-50.)
	    Y1    = -50.
	    SHsum = Y1
	    for i in range(X2_X1/dt): 
	        SHsum = SHsum + ( Y1 + ((i+1)*dt) * Y2_Y1/X2_X1 )
            SHsum = SHsum*dt*3600.

	    X2_X1 = 0.19150E+02 - AgsDict['Mxl_Ags_output_dyn']['UTC(hours)'][len(AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'])-1]
	    Y2_Y1 = -50. - AgsDict['Mxl_Ags_output_dyn']['SH(W.m-2)'][len(AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'])-1]
	    Y1    = AgsDict['Mxl_Ags_output_dyn']['SH(W.m-2)'][len(AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'])-1]
	    SHsum2 = Y1
	    for i in range(X2_X1/dt): 
	        SHsum2 = SHsum2 + ( Y1 + ((i+1)*dt) * Y2_Y1/X2_X1 )
            SHsum2 = SHsum2*dt*3600.
	    
	    print '                                         Sum SH:   %6.2f + %6.2f + %6.2f = %6.2f'%(SumSH/1000000.,SHsum/1000000.,SHsum2/1000000.,(SumSH+SHsum+SHsum2)/1000000.)
	
            #NEE:
            X2_X1 = AgsDict['Mxl_Ags_output_sca']['UTC(hours)'][0] - 0.41333E+01
	    Y2_Y1 = (AgsDict['Mxl_Ags_output_sca']['An(mgCO2.m-2.s-1)'][0] + AgsDict['Mxl_Ags_output_sca']['Resp(mgCO2.m-2.s-1)'][0])/1000. - 0
	    Y1    = 0.
	    NEEsum = Y1
	    for i in range(X2_X1/dt): 
	        NEEsum = NEEsum + ( Y1 + ((i+1)*dt) * Y2_Y1/X2_X1 )
            NEEsum = NEEsum*dt*3600.

	    X2_X1 = 0.19150E+02 - AgsDict['Mxl_Ags_output_sca']['UTC(hours)'][len(AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'])-1]
	    Y2_Y1 = 0 - (AgsDict['Mxl_Ags_output_sca']['An(mgCO2.m-2.s-1)'][len(AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'])-1] + AgsDict['Mxl_Ags_output_sca']['Resp(mgCO2.m-2.s-1)'][len(AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'])-1])/1000. 
	    Y1    = (AgsDict['Mxl_Ags_output_sca']['An(mgCO2.m-2.s-1)'][len(AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'])-1] + AgsDict['Mxl_Ags_output_sca']['Resp(mgCO2.m-2.s-1)'][len(AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'])-1])/1000.
	    NEEsum2 = Y1
	    for i in range(X2_X1/dt): 
	        NEEsum2 = NEEsum2 + ( Y1 + ((i+1)*dt) * Y2_Y1/X2_X1 )
            NEEsum2 = NEEsum2*dt*3600.
	    
	    print '                                         Sum NEE:  %6.2f + %6.2f + %6.2f = %6.2f'%(SumNEE,NEEsum,NEEsum2,(SumNEE+NEEsum+NEEsum2))

            print 'dt =', dt	    
############
# MXL-GECROS
############

        if (i==2):
	    dt = GecDict['MXL_weather']['UTC(hours)'][1] - GecDict['MXL_weather']['UTC(hours)'][0]
	    for j,val in enumerate(GecDict['MXL_weather']['SWin(W.m-2)']):
                if (val > 0. and 8. <= GecDict['MXL_weather']['UTC(hours)'][j] <= 18.):
	            SumSwin  = SumSwin   + val
		    SumSwnet = SumSwnet  + (GecDict['MXL_weather']['ATRJSU(W.m-2)'][j]+GecDict['MXL_weather']['ATRJSH(W.m-2)'][j]+GecDict['MXL_weather']['ATRJsoil(W.m-2)'][j])
		    SumLwnet = SumLwnet  + (GecDict['MXL_weather']['LW_SU(W.m-2)'][j]+GecDict['MXL_weather']['LW_SH(W.m-2)'][j]+GecDict['MXL_weather']['LW_soil(W.m-2)'][j])
		    SumQnet  = SumQnet   + (GecDict['MXL_weather']['Qnet_SU(W.m-2)'][j]+GecDict['MXL_weather']['Qnet_SH(W.m-2)'][j]+GecDict['MXL_weather']['Qnet_soil(W.m-2)'][j])
		    if (GecDict['MXL_weather']['LE_Turb(W.m-2)'][j] >0.): 
		    	SumLE    = SumLE     + GecDict['MXL_weather']['LE_Turb(W.m-2)'][j]
		    SumSH    = SumSH     + GecDict['MXL_weather']['SH_Turb(W.m-2)'][j]
		    SumGR    = SumGR     + GecDict['MXL_weather']['GR_Turb(W.m-2)'][j]
		    SumNEE   = SumNEE    + (-GecDict['MXL_weather']['Anet_Turb(gCO2.m-2.s-1)'][j]/1000.+GecDict['MXL_weather']['FromDailyTot.To.Instant'][j]*(5.98+27.40))
		    n        = n         + 1.
		elif (GecDict['MXL_weather']['UTC(hours)'][j] > 18.):
	            SumSwin  = SumSwin   - GecDict['MXL_weather']['SWin(W.m-2)'][j-1]
		    SumSwnet = SumSwnet  - (GecDict['MXL_weather']['ATRJSU(W.m-2)'][j-1]+GecDict['MXL_weather']['ATRJSH(W.m-2)'][j-1]+GecDict['MXL_weather']['ATRJsoil(W.m-2)'][j-1])
		    SumLwnet = SumLwnet  - (GecDict['MXL_weather']['LW_SU(W.m-2)'][j-1]+GecDict['MXL_weather']['LW_SH(W.m-2)'][j-1]+GecDict['MXL_weather']['LW_soil(W.m-2)'][j-1])
		    SumQnet  = SumQnet   - (GecDict['MXL_weather']['Qnet_SU(W.m-2)'][j-1]+GecDict['MXL_weather']['Qnet_SH(W.m-2)'][j-1]+GecDict['MXL_weather']['Qnet_soil(W.m-2)'][j-1])
		    SumLE    = SumLE     - GecDict['MXL_weather']['LE_Turb(W.m-2)'][j-1]
		    SumSH    = SumSH     - GecDict['MXL_weather']['SH_Turb(W.m-2)'][j-1]
		    SumGR    = SumGR     - GecDict['MXL_weather']['GR_Turb(W.m-2)'][j-1]
		    SumNEE   = SumNEE    - (-GecDict['MXL_weather']['Anet_Turb(gCO2.m-2.s-1)'][j-1]/1000.+GecDict['MXL_weather']['FromDailyTot.To.Instant'][j-1]*(5.98+27.40))
		    break
            SumSwin   = SumSwin * dt * 3600.
            SumSwnet  = SumSwnet * dt * 3600.
            SumLwnet  = SumLwnet * dt * 3600.
	    SumQnet   = SumQnet * dt * 3600.
	    SumLE     = SumLE * dt * 3600.
	    SumSH     = SumSH * dt * 3600.
	    SumGR     = SumGR * dt * 3600.
	    SumNEE    = SumNEE * dt * 3600.#*12./44.
	
	    print '##################'
	    print '\nMXL-GECROS model from 6am to 6pm:'
	    print '%12s%12s%12s%12s%12s%12s%12s%12s'%('Sum SWin','Sum SWnet','Sum LWnet','Sum Qnet','Sum LE','Sum SH','Sum GR','Sum NEE')
	    print '%12.2f%12.2f%12.2f%12.2f%12.2f%12.2f%12.2f%12.2f'%(SumSwin/1000000., SumSwnet/1000000.,SumLwnet/1000000.,SumQnet/1000000.,SumLE/1000000.,SumSH/1000000.,SumGR/1000000.,SumNEE)

            #linear interpolation between sunrise and 6hUTC for Qnet, LE, SH, NEE:

            #SWin:
	    dt = AgsDict['Mxl_Ags_output_dyn']['UTC(hours)'][1] - AgsDict['Mxl_Ags_output_dyn']['UTC(hours)'][0]
            X2_X1 = AgsDict['Mxl_Ags_output_dyn']['UTC(hours)'][0] - 0.41333E+01
	    Y2_Y1 = (AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'][0] - 0.65782E+00 )
	    Y1    = 0.65782E+00
	    SWinsum = Y1
	    for i in range(X2_X1/dt): 
	        SWinsum = SWinsum + ( Y1 + ((i+1)*dt) * Y2_Y1/X2_X1 )
            SWinsum = SWinsum*dt*3600.

	    X2_X1 = 0.19150E+02 - AgsDict['Mxl_Ags_output_dyn']['UTC(hours)'][len(AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'])-1]
	    Y2_Y1 = (0.58628E+00 - AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'][len(AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'])-1] )
	    Y1    = AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'][len(AgsDict['Mxl_Ags_output_dyn']['Swin(W.m-2)'])-1]
	    SWinsum2 = Y1
	    for i in range(X2_X1/dt): 
	        SWinsum2 = SWinsum2 + ( Y1 + ((i+1)*dt) * Y2_Y1/X2_X1 )
            SWinsum2 = SWinsum2*dt*3600.
	    
	    print '\nMXL-GECROS model from sunrise to sunset: Sum SWin: %6.2f + %6.2f + %6.2f = %6.2f'%(SumSwin/1000000.,SWinsum/1000000.,SWinsum2/1000000.,(SumSwin+SWinsum+SWinsum2)/1000000.)


            #Qnet:
	    dt = GecDict['MXL_weather']['UTC(hours)'][1] - GecDict['MXL_weather']['UTC(hours)'][0]
            X2_X1 = GecDict['MXL_weather']['UTC(hours)'][0] - GecDict2['Haarweg_weather']['UTC(hours)'][355]
	    Y2_Y1 = (GecDict['MXL_weather']['Qnet_SU(W.m-2)'][0]+GecDict['MXL_weather']['Qnet_SH(W.m-2)'][0]+GecDict['MXL_weather']['Qnet_soil(W.m-2)'][0] - 
	             (GecDict2['Haarweg_weather']['Qnet_SU(W.m-2)'][355]+GecDict2['Haarweg_weather']['Qnet_SH(W.m-2)'][355]+GecDict2['Haarweg_weather']['Qnet_soil(W.m-2)'][355]) )
	    Y1    = (GecDict2['Haarweg_weather']['Qnet_SU(W.m-2)'][355]+GecDict2['Haarweg_weather']['Qnet_SH(W.m-2)'][355]+GecDict2['Haarweg_weather']['Qnet_soil(W.m-2)'][355])
	    Qnetsum = Y1
	    for i in range(X2_X1/dt): 
	        Qnetsum = Qnetsum + ( Y1 + ((i+1)*dt) * Y2_Y1/X2_X1 )
            Qnetsum = Qnetsum*dt*3600.

	    X2_X1 = GecDict2['Haarweg_weather']['UTC(hours)'][359] - GecDict['MXL_weather']['UTC(hours)'][143]
	    Y2_Y1 = (GecDict2['Haarweg_weather']['Qnet_SU(W.m-2)'][359]+GecDict2['Haarweg_weather']['Qnet_SH(W.m-2)'][359]+GecDict2['Haarweg_weather']['Qnet_soil(W.m-2)'][359] - 
	            (GecDict['MXL_weather']['Qnet_SU(W.m-2)'][143]+GecDict['MXL_weather']['Qnet_SH(W.m-2)'][143]+GecDict['MXL_weather']['Qnet_soil(W.m-2)'][143]) )
	    Y1    = GecDict['MXL_weather']['Qnet_SU(W.m-2)'][143]+GecDict['MXL_weather']['Qnet_SH(W.m-2)'][143]+GecDict['MXL_weather']['Qnet_soil(W.m-2)'][143]
	    Qnetsum2 = Y1
	    for i in range(X2_X1/dt): 
	        Qnetsum2 = Qnetsum2 + ( Y1 + ((i+1)*dt) * Y2_Y1/X2_X1 )
            Qnetsum2 = Qnetsum2*dt*3600.
	    
	    print '                                         Sum Qnet: %6.2f + %6.2f + %6.2f = %6.2f'%(SumQnet/1000000.,Qnetsum/1000000.,Qnetsum2/1000000.,(SumQnet+Qnetsum+Qnetsum2)/1000000.)

	    # LE:
	    X2_X1 = GecDict['MXL_weather']['UTC(hours)'][0] - GecDict2['Haarweg_weather']['UTC(hours)'][355]
	    Y2_Y1 = GecDict['MXL_weather']['LE_Turb(W.m-2)'][0] - GecDict2['Haarweg_weather']['LE_Turb(W.m-2)'][355]
	    Y1    = GecDict2['Haarweg_weather']['LE_Turb(W.m-2)'][355]
	    LEsum = Y1
	    for i in range(X2_X1/dt): 
	        LEsum = LEsum + ( Y1 + ((i+1)*dt) * Y2_Y1/X2_X1 )
            LEsum = LEsum*dt*3600.

	    X2_X1 = GecDict2['Haarweg_weather']['UTC(hours)'][359] - GecDict['MXL_weather']['UTC(hours)'][143]
	    Y2_Y1 = GecDict2['Haarweg_weather']['LE_Turb(W.m-2)'][359] - GecDict['MXL_weather']['LE_Turb(W.m-2)'][143]
	    Y1    = GecDict['MXL_weather']['LE_Turb(W.m-2)'][143]
	    LEsum2 = Y1
	    for i in range(X2_X1/dt): 
	        LEsum2 = LEsum2 + ( Y1 + ((i+1)*dt) * Y2_Y1/X2_X1 )
            LEsum2 = LEsum2*dt*3600.
	    
	    print '                                         Sum LE:   %6.2f + %6.2f + %6.2f = %6.2f'%(SumLE/1000000.,LEsum/1000000.,LEsum2/1000000.,(SumLE+LEsum+LEsum2)/1000000.)


	    # SH:
	    X2_X1 = GecDict['MXL_weather']['UTC(hours)'][0] - GecDict2['Haarweg_weather']['UTC(hours)'][355]
	    Y2_Y1 = GecDict['MXL_weather']['SH_Turb(W.m-2)'][0] - GecDict2['Haarweg_weather']['SH_Turb(W.m-2)'][355]
	    Y1    = GecDict2['Haarweg_weather']['SH_Turb(W.m-2)'][355]
	    SHsum = Y1
	    for i in range(X2_X1/dt): 
	        SHsum = SHsum + ( Y1 + ((i+1)*dt) * Y2_Y1/X2_X1 )
            SHsum = SHsum*dt*3600.

	    X2_X1 = GecDict2['Haarweg_weather']['UTC(hours)'][359] - GecDict['MXL_weather']['UTC(hours)'][143]
	    Y2_Y1 = GecDict2['Haarweg_weather']['SH_Turb(W.m-2)'][359] - GecDict['MXL_weather']['SH_Turb(W.m-2)'][143]
	    Y1    = GecDict['MXL_weather']['SH_Turb(W.m-2)'][143]
	    SHsum2 = Y1
	    for i in range(X2_X1/dt): 
	        SHsum2 = SHsum2 + ( Y1 + ((i+1)*dt) * Y2_Y1/X2_X1 )
            SHsum2 = SHsum2*dt*3600.
	    
	    print '                                         Sum SH:   %6.2f + %6.2f + %6.2f = %6.2f'%(SumSH/1000000.,SHsum/1000000.,SHsum2/1000000.,(SumSH+SHsum+SHsum2)/1000000.)
	
            #NEE:
            X2_X1 = GecDict['MXL_weather']['UTC(hours)'][0] - GecDict2['Haarweg_weather']['UTC(hours)'][355]
	    Y2_Y1 = ( (-GecDict['MXL_weather']['Anet_Turb(gCO2.m-2.s-1)'][0]/1000.+GecDict['MXL_weather']['FromDailyTot.To.Instant'][0]*(5.98+27.40)) - 
	             (-GecDict2['Haarweg_weather']['Anet_Turb(gCO2.m-2.s-1)'][355]/1000.+GecDict2['Haarweg_weather']['FromDailyTot.To.Instant'][355]*(5.98+27.40)) )
	    Y1    = (-GecDict2['Haarweg_weather']['Anet_Turb(gCO2.m-2.s-1)'][355]/1000.+GecDict2['Haarweg_weather']['FromDailyTot.To.Instant'][355]*(5.98+27.40))
	    NEEsum = Y1
	    for i in range(X2_X1/dt): 
	        NEEsum = NEEsum + ( Y1 + ((i+1)*dt) * Y2_Y1/X2_X1 )
            NEEsum = NEEsum*dt*3600.

	    X2_X1 = GecDict2['Haarweg_weather']['UTC(hours)'][359] - GecDict['MXL_weather']['UTC(hours)'][143]
	    Y2_Y1 = ( (-GecDict2['Haarweg_weather']['Anet_Turb(gCO2.m-2.s-1)'][359]/1000.+GecDict2['Haarweg_weather']['FromDailyTot.To.Instant'][359]*(5.98+27.40)) - 
	            (-GecDict['MXL_weather']['Anet_Turb(gCO2.m-2.s-1)'][143]/1000.+GecDict['MXL_weather']['FromDailyTot.To.Instant'][143]*(5.98+27.40)) )
	    Y1    = (-GecDict['MXL_weather']['Anet_Turb(gCO2.m-2.s-1)'][143]/1000.+GecDict['MXL_weather']['FromDailyTot.To.Instant'][143]*(5.98+27.40))
	    NEEsum2 = Y1
	    for i in range(X2_X1/dt): 
	        NEEsum2 = NEEsum2 + ( Y1 + ((i+1)*dt) * Y2_Y1/X2_X1 )
            NEEsum2 = NEEsum2*dt*3600.
	    
	    print '                                         Sum NEE:  %6.2f + %6.2f + %6.2f = %6.2f'%(SumNEE,NEEsum,NEEsum2,(SumNEE+NEEsum+NEEsum2))

            print 'dt =', dt,'\n'
