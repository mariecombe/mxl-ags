# test model
# run mlmodel

from pylab import *
from model_ags import *

run1input = modelinput()

run1input.dt         = 60.       # time step [s]
run1input.runtime    = 57600.    # total run time [s]

# mixed-layer input
run1input.sw_ml      = True      # mixed-layer model switch
run1input.h          = 200.      # initial ABL height [m]
run1input.Ps         = 101300.   # surface pressure [Pa]
run1input.ws         = 1.5e-5        # large scale vertical velocity [m s-1]
run1input.fc         = 1.e-4     # Coriolis parameter [m s-1]                        # not in namoption?

run1input.theta      = 283.0     # initial mixed-layer potential temperature [K]     
run1input.dtheta     = 1.        # initial temperature jump at h [K]
run1input.gammatheta = 0.006     # free atmosphere potential temperature lapse rate [K m-1]
run1input.advtheta   = 0.        # advection of heat [K s-1]
run1input.beta       = 0.20      # entrainment ratio for virtual heat [-]
run1input.wtheta     = 0.        # surface kinematic heat flux [K m s-1]              # in case flux is prescribed?

run1input.q          = 0.008     # initial mixed-layer specific humidity [kg kg-1]
run1input.dq         = -1.0e-3   # initial specific humidity jump at h [kg kg-1]
run1input.gammaq     = 0.        # free atmosphere specific humidity lapse rate [kg kg-1 m-1]
run1input.advq       = 0.        # advection of moisture [kg kg-1 s-1]
run1input.wq         = 0.        # surface kinematic moisture flux [kg kg-1 m s-1]   # in case flux is prescribed?

run1input.co2        = 422.      # initial mixed-layer carbon dioxide [ppm]
run1input.dco2       = -44.      # initial carbon dioxide jump at h [ppm]
run1input.gammaco2   = 0.        # free atmosphere carbon dioxide lapse rate [ppm m-1]
run1input.advco2     = 0.        # advection of carbon dioxide [ppm s-1]             # not in namoption?
run1input.wco2       = 0.        # surface kinematic carbon dioxide flux [ppm m s-1]

run1input.sw_wind    = True      # prognostic wind switch
run1input.u          = 5.        # initial mixed-layer u-wind speed [m s-1]
run1input.du         = 5.        # initial u-wind jump at h [m s-1]                  # == ug - um0
run1input.gammau     = 0.        # free atmosphere u-wind speed lapse rate [s-1]
run1input.advu       = 0.        # advection of u-wind [m s-2]                       # not in namoption?

run1input.v          = -5.       # initial mixed-layer u-wind speed [m s-1]
run1input.dv         = 5.        # initial v-wind jump at h [m s-1]
run1input.gammav     = 0.        # free atmosphere v-wind speed lapse rate [s-1]
run1input.advv       = 0.        # advection of v-wind [m s-2]                       # not in namoption?

# surface layer input
run1input.sw_sl      = True      # surface layer switch
run1input.ustar      = 0.3       # surface friction velocity [m s-1]                 # not in namoption?
run1input.z0m        = 0.02      # roughness length for momentum [m]
run1input.z0h        = 0.002     # roughness length for scalars [m]

# radiation parameters
run1input.sw_rad     = True      # radiation switch
run1input.lat        = 51.97     # latitude [deg]
run1input.lon        = 5.67     # longitude [deg]
run1input.doy        = 153.      # day of the year [-]
run1input.tstart     = 4.0       # time of the day [h UTC]
run1input.cc         = 0.0       # cloud cover fraction [-]
run1input.Q          = 400.      # net radiation [W m-2] 

# land surface parameters
run1input.sw_ls      = True      # land surface switch
run1input.wg         = 0.25      # volumetric water content top soil layer [m3 m-3]
run1input.w2         = 0.25      # volumetric water content deeper soil layer [m3 m-3]
run1input.cveg       = 1.        # vegetation fraction [-]
run1input.Tsoil      = 285.      # temperature top soil layer [K]
run1input.T2         = 285.      # temperature deeper soil layer [K]
run1input.a          = 0.219     # Clapp and Hornberger retention curve parameter a
run1input.b          = 4.9       # Clapp and Hornberger retention curve parameter b
run1input.p          = 4.        # Clapp and Hornberger retention curve parameter c
run1input.CGsat      = 3.56e-6   # saturated soil conductivity for heat

run1input.wsat       = 0.472     # saturated volumetric water content ECMWF config [-]
run1input.wfc        = 0.323     # volumetric water content field capacity [-]
run1input.wwilt      = 0.171     # volumetric water content wilting point [-]

run1input.C1sat      = 0.132
run1input.C2ref      = 1.8

run1input.LAI        = 3.        # leaf area index [-]
run1input.gD         = 0.0       # correction factor transpiration for VPD [-]
run1input.rsmin      = 110.      # minimum resistance transpiration [s m-1]
run1input.rssoilmin  = 0.        # minimun resistance soil evaporation [s m-1]
run1input.alpha      = 0.25      # surface albedo [-]

run1input.Ts         = 283.      # initial surface temperature [K]

run1input.Wmax       = 2.e-4     # thickness of water layer on wet vegetation [m]      # parameter in bulk_chemistry.f90
run1input.Wl         = 0.0       # equivalent water layer depth for wet vegetation [m]

run1input.Lambda     = 5.9       # thermal diffusivity skin layer [-]

# Plant physiology model (A-gs)
run1input.CO2comp298 = 68.5                      # CO2 compensation concentration [mg m-3]
run1input.Q10CO2     = 1.5                       # function parameter to calculate CO2 compensation concentration [-]
run1input.gm298      = 7.0                       # mesophyill conductance at 298 K [mm s-1]
run1input.Ammax298   = 2.2                       # CO2 maximal primary productivity [mg m-2 s-1]
run1input.Q10gm      = 2.0                       # function parameter to calculate mesophyll conductance [-] 
run1input.T1gm       = 278.                      # reference temperature to calculate mesophyll conductance gm [K]
run1input.T2gm       = 301.                      # reference temperature to calculate mesophyll conductance gm [K]
run1input.Q10Am      = 2.0                       # function parameter to calculate maximal primary profuctivity Ammax
run1input.T1Am       = 281.                      # reference temperature to calculate maximal primary profuctivity Ammax [K] 
run1input.T2Am       = 311.                      # reference temperature to calculate maximal primary profuctivity Ammax [K]
run1input.f0         = 0.89                      # maximum value Cfrac [-]
run1input.ad         = 0.07                      # regression coefficient to calculate Cfrac [kPa-1]
run1input.alpha0     = 0.017                     # initial low light conditions [mg J-1]
run1input.Kx         = 0.7                       # extinction coefficient PAR [-]
run1input.gmin       = 0.25 / 1000.              # cuticular (minimum) conductance [m s-1]

# Soil respiration model (coupled to A-gs)
run1input.Cw         = 0.0016                    # constant water stress correction (eq. 13 Jacobs et al. 2007) [-]
run1input.wmax       = 0.55                      # upper reference value soil water [-] 
run1input.wmin       = 0.005                     # lower reference value soil water [-] 
run1input.R10        = 0.23                      # respiration at 10 C [mg CO2 m-2 s-1] 
run1input.E0         = 53.3e3                    # activation energy [53.3 kJ kmol-1] 
 
run1ags = model(run1input)
run1ags.runmodel()
#
run2ags = model(run1input)
run2ags.input.dtheta= 1.       # initial temperature jump at h [K]
run2ags.runmodel()

# plotting the results

matplotlib.pyplot.close('all')

#hdata0        = loadtxt('data/BLheight.txt', usecols=[1,2])
#fluxdata0     = loadtxt('data/cabsurf_surface_flux_200309-24-25-26.lot', usecols=[2,3,4],skiprows=4)
#fluxco20      = loadtxt('data/cabsurf_surface_flux_200309-24-25-26.lot', usecols=[8],skiprows=4)
#tempdata0     = loadtxt('data/caboper_air_temperature_200309-24-25-26.lot', usecols=[2,3,4,5,6,7,8],skiprows=4)
#tddata0       = loadtxt('data/caboper_dew_point_200309-24-25-26.lot', usecols=[2,3,4,5,6,7,8],skiprows=4)
#CO2data0_060  = loadtxt('data/CO2-September2003.txt', usecols=[3],skiprows=4)
#CO2data0_120  = loadtxt('data/CO2-September2003.txt', usecols=[4],skiprows=4)
#CO2data0_200  = loadtxt('data/CO2-September2003.txt', usecols=[5],skiprows=4)

#hdata         = hdata0[:,1]
#hdatatime     = hdata0[:,0]
#fluxdata      = fluxdata0[180:252,1::]
#fluxco2       = fluxco20[162:252]
#tempdata080   = tempdata0[180:252,3]
#tempdata200   = tempdata0[180:252,1]
#tddata080     = tddata0[180:252,3]
#tddata200     = tddata0[180:252,1]
cabtime       = arange(6. + 1./12., 18.+ 1./12. , 1./6.)
#cabtimeco2    = arange(3. + 1./12., 18.+ 1./12. , 1./6.)
#CO2datatime   = arange(6., 19., 1.)
#CO2data060    = CO2data0_060[29:42]
#CO2data120    = CO2data0_120[29:42]
#CO2data200    = CO2data0_200[29:42]
print(cabtime,run1ags.out.h)
print(cabtime,run1ags.out.rs)
print(cabtime,run1ags.out.theta)

#print(cabtimeco2,tempdata200)

figure()
subplot(221)
plot(run1ags.out.t, run1ags.out.h, 'b-')
plot(run2ags.out.t, run2ags.out.h, 'b--')
#plot(hdatatime, hdata, 'b^')
xlabel('time UTC [h]')
ylabel('h [m]')
#axis([7.5, 15.7, 0, 1400])

subplot(222)
plot(run1ags.out.t, run1ags.out.H, 'r-')
plot(run2ags.out.t, run2ags.out.H, 'r--')
#plot(cabtime,fluxdata[:,0], 'ro')

#plot(run1ags.out.t, run1ags.out.LE, 'b-')
#plot(run2ags.out.t, run2ags.out.LE, 'b--')
#plot(cabtime,fluxdata[:,1], 'b^')
xlabel('time UTC [h]')
ylabel(u'heat flux [W m\u207B\u00B2]')
#axis([7.5, 15.7, 0, 300])

#for n in range(len(cabtime)):
#  if(tempdata080[n] == -9999.0):
#    tempdata080[n] = NaN
#  if(tddata080[n] == -9999.0):
#    tddata080[n] = NaN
#  if(tempdata200[n] == -9999.0):
#    tempdata200[n] = NaN
#  if(tddata200[n] == -9999.0):
#    tddata200[n] = NaN
#
#tempdata080 = tempdata080 + 273.16
#tddata080   = tddata080   + 273.16
#tempdata200 = tempdata200 + 273.16
#tddata200   = tddata200   + 273.16
#
#pottempdata080 = tempdata080 + 9.81 / 1005. * 080.
#pottempdata200 = tempdata200 + 9.81 / 1005. * 200.
#
subplot(223)
plot(run1ags.out.t, run1ags.out.theta, 'b-')
plot(run2ags.out.t, run2ags.out.theta, 'b--')
##plot(cabtime, pottempdata080, 'b^', label='080 m')
##plot(cabtime, pottempdata200, 'r^', label='200 m')
##leg1 = legend(loc=0)
##leg1.draw_frame(False)
xlabel('time UTC [h]')
ylabel(u'theta [K]')
#axis([7.5, 15.7, 284., 292.])
#
#edata080 = 610.7 * exp ( (17.2694*(tddata080 - 273.16)) / (tddata080 - 35.86) )
#qdata080 = 0.622 * edata080 / (102500. - 1.2 * 9.81 * 080.) * 1000.
#edata200 = 610.7 * exp ( (17.2694*(tddata200 - 273.16)) / (tddata200 - 35.86) )
#qdata200 = 0.622 * edata200 / (102500. - 1.2 * 9.81 * 200.) * 1000.
#
subplot(224)
plot(run1ags.out.t, 1000.*run1ags.out.q, 'b-')
plot(run2ags.out.t, 1000.*run2ags.out.q, 'b--')
##plot(cabtime, qdata080, 'b^', label='080 m')
##plot(cabtime, qdata200, 'r^', label='200 m')
##leg1 = legend(loc=0)
##leg1.draw_frame(False)
xlabel('time UTC [h]')
ylabel(u'q [g kg\u207B\u00B9]')
#axis([7.5, 15.7, 0.0, 6.])
#
figure()
subplot(211)
plot(run1ags.out.t, run1ags.out.rsAgs, 'g-')
plot(run2ags.out.t, run2ags.out.rsAgs, 'g--')
xlabel('time UTC [h]')
ylabel(u'rsags [mm s-1]')
#axis([7.5, 15.7,0.0,200.])
#
subplot(212)
plot(run1ags.out.t, run1ags.out.An,   'g-', label='An')
plot(run2ags.out.t, run2ags.out.An,   'g--')
plot(run1ags.out.t, run1ags.out.Resp, 'r-', label='Rs')
plot(run2ags.out.t, run2ags.out.Resp, 'r--')
##plot(cabtimeco2,fluxco2, 'r^')
leg1 = legend(loc=0)
leg1.draw_frame(False)
xlabel('time UTC [h]')
ylabel(u'wCO2 [mgC m-2 s-1]')
##axis([7.5, 15.7,-1.4,0.6])
#axis([3.0, 15.7,-1.4,0.6])
#
#figure()
#subplot(211)
#plot(run1ags.out.t, run1ags.out.awco2, 'g-', label='Control')
#plot(run2ags.out.t, run2ags.out.awco2, 'g--',label='dth = 2.1 K')
#leg1 = legend(loc=0)
#leg1.draw_frame(False)
##plot(cabtimeco2,fluxco2, 'g^')
##plot(CO2datatime, CO2data, 'go')
#xlabel('time UTC [h]')
#ylabel(u'wCO2 [mgC m-2 s-1]')
#axis([3.0, 15.7,-1.0,0.1])
#
#subplot(212)
#plot(run1ags.out.t, run1ags.out.co2, 'g-')
#plot(run2ags.out.t, run2ags.out.co2, 'g--')
##leg1 = legend(loc=0)
##leg1.draw_frame(False)
##legend()
##plot(CO2datatime, CO2data060, 'ro', label='060 m')
##plot(CO2datatime, CO2data120, 'go', label='120 m')
##plot(CO2datatime, CO2data200, 'bo', label='200 m')
##leg1 = legend(loc=0)
##leg1.draw_frame(False)
#xlabel('time UTC [h]')
#ylabel(u'CO2 [ppm]')
#axis([7.5, 15.7,370,420])
###savefig('sensdtheta.eps',dpi=300)

show()
