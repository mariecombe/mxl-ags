# test model
# run mlmodel

from pylab import *
from model_ags import *

run1input = modelinput()

run1input.dt         = 60.       # time step [s]
run1input.runtime    = 36000.    # total run time [s]

# mixed-layer input
run1input.sw_ml      = True      # mixed-layer model switch
run1input.h          = 175.      # initial ABL height [m]
run1input.Ps         = 102900.   # surface pressure [Pa]
run1input.ws         = 0.        # large scale vertical velocity [m s-1]
run1input.fc         = 1.e-4     # Coriolis parameter [m s-1]

run1input.theta      = 284.0     # initial mixed-layer potential temperature [K]
run1input.dtheta     = 4.2       # initial temperature jump at h [K]
run1input.gammatheta = 0.0036    # free atmosphere potential temperature lapse rate [K m-1]
run1input.advtheta   = 0.        # advection of heat [K s-1]
run1input.beta       = 0.20      # entrainment ratio for virtual heat [-]
run1input.wtheta     = 0.1       # surface kinematic heat flux [K m s-1]     ## what does that corresponds to?

run1input.q          = 0.0049    # initial mixed-layer specific humidity [kg kg-1]
run1input.dq         = -8.0e-4   # initial specific humidity jump at h [kg kg-1]
run1input.gammaq     = -1.2e-6   # free atmosphere specific humidity lapse rate [kg kg-1 m-1]
run1input.advq       = 0.        # advection of moisture [kg kg-1 s-1]
run1input.wq         = 0.1       # surface kinematic moisture flux [kg kg-1 m s-1]  ## what does that corresponds to?

run1input.co2        = 422.      # initial mixed-layer carbon dioxide [ppm]
run1input.dco2       = -48.      # initial carbon dioxide jump at h [ppm]
run1input.gammaco2   = -0.0e-03  # free atmosphere carbon dioxide lapse rate [ppm m-1]
run1input.advco2     = 0.        # advection of carbon dioxide [ppm s-1]
run1input.wco2       = -0.0e-04  # surface kinematic carbon dioxide flux [ppm m s-1]

run1input.sw_wind    = True      # prognostic wind switch
run1input.u          = 5.        # initial mixed-layer u-wind speed [m s-1]
run1input.du         = 3.        # initial u-wind jump at h [m s-1]
run1input.gammau     = 0.002     # free atmosphere u-wind speed lapse rate [s-1]
run1input.advu       = 0.        # advection of u-wind [m s-2]

run1input.v          = 0.0       # initial mixed-layer u-wind speed [m s-1]
run1input.dv         = 0.0       # initial v-wind jump at h [m s-1]
run1input.gammav     = 0.        # free atmosphere v-wind speed lapse rate [s-1]
run1input.advv       = 0.        # advection of v-wind [m s-2]

# surface layer input
run1input.sw_sl      = False     # surface layer switch
run1input.ustar      = 0.3       # surface friction velocity [m s-1]  #### where do I input this??
run1input.z0m        = 0.05      # roughness length for momentum [m]
run1input.z0h        = 0.01      # roughness length for scalars [m]

# radiation parameters
run1input.sw_rad     = True      # radiation switch        
run1input.lat        = 51.97     # latitude [deg]
run1input.lon        = -4.93     # longitude [deg]
run1input.doy        = 268.      # day of the year [-]
run1input.tstart     = 6.0       # time of the day [h UTC]
run1input.cc         = 0.0       # cloud cover fraction [-]
run1input.Q          = 400.      # net radiation [W m-2] 

# land surface parameters
run1input.sw_ls      = True      # land surface switch
run1input.wg         = 0.42      # volumetric water content top soil layer [m3 m-3]
run1input.w2         = 0.42      # volumetric water content deeper soil layer [m3 m-3]
run1input.cveg       = 0.9       # vegetation fraction [-]
run1input.Tsoil      = 282.      # temperature top soil layer [K]
run1input.T2         = 285.      # temperature deeper soil layer [K]
run1input.a          = 0.083     # Clapp and Hornberger retention curve parameter a
run1input.b          =11.4       # Clapp and Hornberger retention curve parameter b
run1input.p          =12.        # Clapp and Hornberger retention curve parameter c
run1input.CGsat      = 3.6e-6    # saturated soil conductivity for heat

run1input.wsat       = 0.6       # saturated volumetric water content ECMWF config [-]
run1input.wfc        = 0.491     # volumetric water content field capacity [-]
run1input.wwilt      = 0.314     # volumetric water content wilting point [-]

run1input.C1sat      = 0.342
run1input.C2ref      = 0.3

run1input.LAI        = 2.        # leaf area index [-]
run1input.gD         = 0.0       # correction factor transpiration for VPD [-]
run1input.rsmin      =110.       # minimum resistance transpiration [s m-1]
run1input.rssoilmin  = 50.       # minimun resistance soil evaporation [s m-1]
run1input.alpha      = 0.25      # surface albedo [-]

run1input.Ts         = 284.      # initial surface temperature [K]

run1input.Wmax       = 0.00020   # thickness of water layer on wet vegetation [m]     #### to set in the fortran source code
run1input.Wl         = 0.00014   # equivalent water layer depth for wet vegetation [m]

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
#run2ags             = model(run1input)
#run3ags             = model(run1input)
#run4ags             = model(run1input)
#
# run 2
#run2ags.input.co2   = 750.      # initial mixed-layer carbon dioxide [ppm]
# run 3
#run3ags.input.theta = 286.0     # initial mixed-layer potential temperature [K]
#run3ags.input.q     = 0.0056    # initial mixed-layer specific humidity [kg kg-1]
#
# run 4
#run4ags.input.co2   = 750.      # initial mixed-layer carbon dioxide [ppm]
#run4ags.input.theta = 286.0     # initial mixed-layer potential temperature [K]
#run4ags.input.q     = 0.0056    # initial mixed-layer specific humidity [kg kg-1]
#
#run2ags.runmodel()
#run3ags.runmodel()
#run4ags.runmodel()
# plotting the results

#print(cabtime,fluxdata[:,1])
#print(cabtime)
#print(cabtimeco2,tempdata200)

fig1= figure(1, figsize=(8,8)) 
fig1.subplots_adjust(0.09,0.08,0.91,0.94,0.02,0.02)
sub1=subplot(321)
#sub1.xaxis.tick_top()
#sub1.set_xticks_top((7,8,9,10,11,12,13,14,15))
#subplot(321)
#
plot(run1ags.out.t, run1ags.out.H,'b-', lw=2, label='YEAR2003' ,color='b')
#plot(run1ags.out.t, run2ags.out.H,'g-', lw=2, label='YEAR2100-C')
#plot(run1ags.out.t, run3ags.out.H,'r-', lw=2, label='YEAR2100-T')
#plot(run1ags.out.t, run4ags.out.H,'k-', lw=2 ,label='YEAR2100')
#plot(run1ags.out.t, run1ags.out.h)
#plot(run1ags.out.t, run2ags.out.h,'b--',  lw=2)
#plot(run1ags.out.t, run3ags.out.h,'b:', lw=2 )
#plot(hdatatime, hdata, 'b^')
for tick in gca().xaxis.iter_ticks():
 tick[0].label2On=True
 tick[0].label1On=False
leg1 = legend(loc=1)
leg1.draw_frame(False)
xlabel('time UTC [h]')
#plt.setp(sub1.get_xticklabels(), visible=False)
ylabel(u'SH [W m\u207B\u00B2]')
#ylabel('h [m]')
#axis([7.5, 15.7, 0, 1400])
axis([6.0, 16., -50, 200])
text(7.5,170,'(a)',fontsize=12)

sub2 = subplot(322)
sub2.yaxis.tick_right()
sub2.yaxis.set_label_position('right')
#sub2.xaxis.tick_top()
plot(run1ags.out.t, run1ags.out.LE,'b-',lw=2)
#plot(run1ags.out.t, run2ags.out.LE,'g-', lw=2)
#plot(run1ags.out.t, run3ags.out.LE,'r-', lw=2)
#plot(run1ags.out.t, run4ags.out.LE,'k-',lw=2)
for tick in gca().xaxis.iter_ticks():
 tick[0].label2On=True
 tick[0].label1On=False
ylabel(u'LE [W m\u207B\u00B2]')
#leg1 = legend(loc=1)
#leg1.draw_frame(False)
axis([6.0, 16., 0, 300])
text(7.5,260,'(b)',fontsize=12)

sub3=subplot(323)
plot(run1ags.out.t, run1ags.out.An,   'b-', lw=2)
#plot(run1ags.out.t, run4ags.out.An,   'k-', lw=2)
plot(run1ags.out.t, run1ags.out.Resp, 'b-', lw=2)
#plot(run1ags.out.t, run4ags.out.Resp, 'k-', lw=2)
#plot(run1ags.out.t, run2ags.out.An,   'g-', lw=2)
#plot(run1ags.out.t, run3ags.out.An,   'r-', lw=2)
#plot(run1ags.out.t, run2ags.out.Resp, 'g-', lw=2)
#plot(run1ags.out.t, run3ags.out.Resp, 'r-', lw=2)
plot(run1ags.out.t, run1ags.out.awco2,'b-', lw=2)
#plot(run1ags.out.t, run4ags.out.awco2,'k-', lw=2)
#plot(run1ags.out.t, run2ags.out.awco2,'g-', lw=2)
#plot(run1ags.out.t, run3ags.out.awco2,'r-', lw=2)
#plot(CO2datatime, CO2data, 'go')
#xlabel('time UTC [h]')
ylabel(u'Fc [mgC m\u207B\u00B2 s\u207B\u00B9]')
plt.setp(sub3.get_xticklabels(), visible=False)
#sub3.xaxis.tick_bottom()
axis([6.0, 16.,-1.5,0.99])
text(7.5,0.72,'(c)',fontsize=12)
text(11.5,-0.80,'Phot',fontsize=15)
text(11.5, 0.1,'NEE',fontsize=15)
text(11.5, 0.5,'Resp',fontsize=15)

sub4=subplot(324)
sub4.yaxis.tick_right()
sub4.yaxis.set_label_position('right')
#sub2.xaxis.tick_top()
plot(run1ags.out.t, run1ags.out.rsAgs,'b-', lw=2)
#plot(run1ags.out.t, run2ags.out.rsAgs,'g-', lw=2)
#plot(run1ags.out.t, run3ags.out.rsAgs,'r-', lw=2) 
#plot(run1ags.out.t, run4ags.out.rsAgs,'k-', lw=2) 
leg1 = legend(loc=0)
#leg1.draw_frame(False)
xlabel('time UTC [h]')
ylabel(u'r$_s$ [mm s\u207B\u00B9]')
plt.setp(sub4.get_xticklabels(), visible=False)
axis([6.0, 16.,100,299])
text(7.7,280,'(d)',fontsize=12)

#fig2= figure(2, figsize=(8,6)) 
#fig2.subplots_adjust(0.09,0.08,0.91,0.94,0.02,0.02)
sub5=subplot(325)
#sub1.xaxis.tick_top()
plot(run1ags.out.t, run1ags.out.h,'b-', lw=2)
#plot(run1ags.out.t, run2ags.out.h,'g-', lw=2)
#plot(run1ags.out.t, run3ags.out.h,'r-', lw=2 )
#plot(run1ags.out.t, run4ags.out.h,'k-', lw=2)
#plot(hdatatime, hdata, 'b^')
#leg1 = legend(loc=3)
#leg1.draw_frame(False)
xlabel('time UTC [h]')
ylabel('h [m]')
axis([6.0, 16., 0, 1399])
text(7.5,1200,'(e)',fontsize=12)
#axis([7.5, 15.7,-1.4,0.6])

sub6=subplot(326)
sub6.yaxis.tick_right()
sub6.yaxis.set_label_position('right')
#sub2.xaxis.tick_top()
plot(run1ags.out.t, run1ags.out.LCL,'b-', label='YEAR2003', lw=2)
#plot(run1ags.out.t, run2ags.out.LCL,'g-', label='YEAR2100-C', lw=2,)
#plot(run1ags.out.t, run3ags.out.LCL,'r-', label='YEAR2100-T', lw=2) 
#plot(run1ags.out.t, run4ags.out.LCL,'k-', label='YEAR2100', lw=2)
leg1 = legend(loc=4)
leg1.draw_frame(False)
xlabel('time UTC [h]')
ylabel(u'r$_s$ [mm s\u207B\u00B9]')
ylabel('LCL [m]')
axis([6.0, 16.,500,1999])
text(7.5,1780,'(f)',fontsize=12)
savefig('figure02.eps',dpi=300)
savefig('figure02.png',dpi=300)
