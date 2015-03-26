#!/usr/bin/env python
# sensitivity.py


"""
Started 5 nov 2014
Author: Marie Combe

Plotting script for the sensitivity analysis results - Article 1

"""

# importing the necessary modules
import os
import sys
import getopt
import numpy as np		# for arrays handling


# custom-made functions
#----------------------------------------------------------------------
def open_data(inpath,filelist):

    Dict = {}

    for i,namefile in enumerate(filelist):

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
        data=np.array(datafloat,float)

        # creating one dictionnary and storing the float data in it
        dictnamelist= {}
        for j,varname in enumerate(headerow):
            dictnamelist[varname]=data[:,j]
        Dict[namefile] = dictnamelist

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
This script produces a sensitivity analysis plot of the MXL-AGs model.
This is a 3-variables plot, presented in 2D (2 variables on x- and y-axis, and one 
result variable in colors).

You should specify not specify any arguments to run this script:
./calculate_article2_text_nbs.py
                
                """
            
            print helptext
            
            sys.exit(2)

    # Working directories
    currentdir = os.getcwd()



    # Create a summary output file - erase old file if exists
    resfile = 'article2_co2_nbs.csv'
    if (os.path.isfile(os.path.join(currentdir,resfile))): os.remove(os.path.join(currentdir,resfile))
    Results = open(os.path.join(currentdir,resfile), 'w')

    # writing the header row
    Results.write('Section 4.3 Modifications under water stress\n')
    Results.write('Section 4.3.1 Atmospheric CO2 budget\n')
    Results.write(',CO2_min(ppm),,gs_12h(mm.s-1),,Anet_int(gCO2.m-2.d-1),,LE_int(MJ.m-2.d-1),,SH_int(MJ.m-2.d-1),,h_max(m),,wce_int(gCO2.m-2.d-1),,tm_max(K),,Resp_int(gCO2.m-2.d-1),,wce_h_int(ppm),,wcs_h_int,\n')
    Results.write(',mean,std,mean,std,mean,std,mean,std,mean,std,mean,std,mean,std,mean,std,mean,std,mean,std,mean,std\n')

    # Open results
    filelist = ['2DIM_SA_thetam0_284-290_cc_0-30_SMI=0.8/Results_thetam0_cc.dat','2DIM_SA_thetam0_284-290_cc_0-30_SMI=0.2/Results_thetam0_cc.dat']    
    ResDict = open_data(currentdir,filelist)
    
    for i,fil in enumerate(filelist):
        z1 = np.ma.masked_invalid(ResDict[fil]['co2_min'])
        z2 = 1000./np.ma.masked_invalid(ResDict[fil]['Rs_14h'])
        z3 = np.ma.masked_invalid(ResDict[fil]['Anet_int'])
        z4 = np.ma.masked_invalid(ResDict[fil]['LE_int'])
        z5 = np.ma.masked_invalid(ResDict[fil]['SH_int'])
        z6 = np.ma.masked_invalid(ResDict[fil]['h_max'])
        z7 = np.ma.masked_invalid(ResDict[fil]['wce_int'])
        z8 = np.ma.masked_invalid(ResDict[fil]['tm_max'])
	z9 = np.ma.masked_invalid(ResDict[fil]['Resp_int'])
	z10 = np.ma.masked_invalid(ResDict[fil]['wce_h_int'])
	z11 = np.ma.masked_invalid(ResDict[fil]['wcs_h_int'])
        Results.write('%s,%6.1f,%6.1f,%6.1f,%6.1f,%6.1f,%6.1f,%6.1f,%6.1f,%6.1f,%6.1f,%6.0f,%6.0f,%6.1f,%6.1f,%6.1f,%6.1f,%6.1f,%6.1f,%6.1f,%6.1f,%6.1f,%6.1f\n'%(fil,
	                                                          np.mean(z1),np.std(z1),np.mean(z2),np.std(z2),np.mean(z3),np.std(z3),
	                                                          np.mean(z4),np.std(z4),np.mean(z5),np.std(z5),np.mean(z6),np.std(z6),
								  np.mean(z7),np.std(z7),np.mean(z8),np.std(z8),np.mean(z9),np.std(z9),
								  np.mean(z10),np.std(z10),np.mean(z11),np.std(z11) )) 

    Results.write('\n\n\n')

    # writing the header row
    Results.write('Section 4.3.2 Ecosystem WUE\n')
    Results.write(',WUEeco_ave(g{C}.kg{H2O}-1),,AF(%),,TF(%),,SH_int(MJ.m-2.d-1),,h_max(m),,tm_max(K),,WUEplt_ave(?),,NEE(gCO2.m-2.d-1),,LE_int(MJ.m-2.d-1),,Anet(gCO2.m-2.d-1),,LEveg_int(MJ.m-2.d-1),,Resp(gCO2.m-2.d-1),,Evap(MJ.m-2.d-1),,CO2_diff,,VPD,,gcco2_gc,,qsat,,qmxl(g.kg-1),\n')
    Results.write(',mean,std,mean,std,mean,std,mean,std,mean,std,mean,std,mean,std,mean,std,mean,std,mean,std,mean,std,mean,std,mean,std,mean,std,mean,std,mean,std,mean,std,mean,std\n')

    # Open results
    filelist = ['2DIM_SA_cc_0-30_gammatheta_2-8_SMI=0.8/Results_cc_ga.dat','2DIM_SA_cc_0-30_gammatheta_2-8_SMI=0.2/Results_cc_ga.dat']    
    ResDict = open_data(currentdir,filelist)

    for i,fil in enumerate(filelist):
        z1 = np.ma.masked_invalid(ResDict[fil]['WUEeco_ave'])
        z2 = 100.*np.ma.masked_invalid(   np.abs(ResDict[fil]['Anet_int'])  /  ( np.abs(ResDict[fil]['Anet_int']) + np.abs(ResDict[fil]['Resp_int']) )   )
        z3 = 100.*np.ma.masked_invalid(   ResDict[fil]['LEveg_int']  /  ResDict[fil]['LE_int']   )
        z4 = np.ma.masked_invalid(ResDict[fil]['SH_int'])
        z5 = np.ma.masked_invalid(ResDict[fil]['h_max'])
        z6 = np.ma.masked_invalid(ResDict[fil]['tm_max'])
        z7 = np.ma.masked_invalid( np.abs(ResDict[fil]['Anet_int']) / ResDict[fil]['LEveg_int'])
        z8 = np.ma.masked_invalid(ResDict[fil]['NEE_int'])
	z9 = np.ma.masked_invalid(ResDict[fil]['LE_int'])
	z10 = np.ma.masked_invalid(ResDict[fil]['Anet_int'])
	z11 = np.ma.masked_invalid(ResDict[fil]['LEveg_int'])
        z12 = np.ma.masked_invalid(ResDict[fil]['Resp_int'])
        z13 = np.ma.masked_invalid(ResDict[fil]['LEsoil_int'])
	
        z14 = np.ma.masked_invalid(ResDict[fil]['CD_ave'])
        z15 = np.ma.masked_invalid(ResDict[fil]['VPD_ave'])
        z16 = np.ma.masked_invalid(   (ResDict[fil]['Ra_14h']+1.6*ResDict[fil]['Rs_14h'])  /  ( ResDict[fil]['Ra_14h'] + ResDict[fil]['Rs_14h'] )   )
	z17 = np.ma.masked_invalid(ResDict[fil]['qm_ave'])
	z18 = np.ma.masked_invalid(ResDict[fil]['qm_ave'])
        Results.write('%s,%6.2f,%6.2f,%6.1f,%6.1f,%6.1f,%6.1f,%6.2f,%6.2f,%6.0f,%6.0f,%6.1f,%6.1f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f\n'%(fil,
	                 np.mean(z1),np.std(z1),np.mean(z2),np.std(z2),np.mean(z3),np.std(z3),np.mean(z4),np.std(z4),np.mean(z5),np.std(z5),np.mean(z6),np.std(z6),
			 np.mean(z7),np.std(z7),np.mean(z8),np.std(z8),np.mean(z9),np.std(z9),np.mean(z10),np.std(z10),np.mean(z11),np.std(z11),np.mean(z12),np.std(z12),
			 np.mean(z13),np.std(z13),np.mean(z14),np.std(z14),np.mean(z15),np.std(z15),np.mean(z16),np.std(z16),np.mean(z17),np.std(z17),np.mean(z18),np.std(z18) )) 

    Results.write('\n\n\n')


    Results.close()
