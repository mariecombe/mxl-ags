# Coupled land-atmosphere-cloud-chemistry model
# Chiel van Heerwaarden & Jordi Vila-Guerau de Arellano, 2009
"""
Open issues:

sensible heat flux always peaks too early in the morning, maybe other Jarvis functions?
"""

import copy
import numpy
#import ribtol

class model:
  def __init__(self, modelinput):
    # initialize the different components of the model
    self.input = copy.deepcopy(modelinput)

  def runmodel(self):
    # initialize model variables
    self.initmodel()

    # time integrate model 
    for self.t in range(self.tsteps):
      
      # time integrate components
      self.timestep()
      
      # store output for current time step
      self.store()

    # delete unnecessary variables from memory
    self.exitmodel()
    
  def initmodel(self):
    # assign variables from input data
    # initialize constants
    self.Lv         =  2.45e6                # heat of vaporization [J kg-1]
    self.cp         =  1005.                 # specific heat of dry air [J kg-1 K-1]
    self.rho        =  1.2                   # density of air [kg m-3]
    self.k          =  0.4                   # Von Karman constant [-]
    self.g          =  9.81                  # gravity acceleration [m s-2]
    self.Rd         =  287.                  # gas constant for dry air [J kg-1 K-1]
    self.Rv         =  461.5                 # gas constant for moist air [J kg-1 K-1]
    self.bolz       =  5.67e-8               # Bolzman constant [-]
    self.rhow       =  1000.                 # density of water [kg m-3]
    self.S0         =  1368.                 # solar constant [W m-2]
    #self.S0         =   980.                 # solar constant [W m-2]
    self.mco2       =  46.                   # molecular weight CO2 [g mol -1]
    self.mair       =  28.9                  # molecular weight air [g mol -1]
    self.nuco2q     =  1.6                   # ratio molecular viscosity water to carbon dioxide
    self.Cw         =  0.0016                # constant water stress correction (eq. 13 Jacobs et al. 2007) [-]
    self.wmax       =  0.55                  # upper reference value soil water [-]
    self.wmin       =  0.005                 # lower reference value soil water [-]
    self.R10        =  0.23                  # respiration at 10 C [mg CO2 m-2 s-1]
    self.E0         =  53.3e3                # activation energy [53.3 kJ kmol-1]

    # initialize mixed-layer
    self.sw_ml      =  self.input.sw_ml      # mixed-layer model switch
    self.h          =  self.input.h          # initial ABL height [m]
    self.Ps         =  self.input.Ps         # surface pressure [Pa]
    self.ws         =  self.input.ws         # large scale vertical velocity [m s-1]
    self.fc         =  self.input.fc         # coriolis parameter [s-1]
    self.we         =  -1.                   # entrainment velocity [m s-1]
    
    self.theta      =  self.input.theta      # initial mixed-layer potential temperature [K]
    self.dtheta     =  self.input.dtheta     # initial temperature jump at h [K]
    self.gammatheta =  self.input.gammatheta # free atmosphere potential temperature lapse rate [K m-1]
    self.advtheta   =  self.input.advtheta   # advection of heat [K s-1]
    self.beta       =  self.input.beta       # entrainment ratio for virtual heat [-]
    self.wtheta     =  self.input.wtheta     # surface kinematic heat flux [K m s-1]

    self.T2m        =  -1.                   # 2m temperature [K]
    self.q2m        =  -1.                   # 2m specific humidity [kg kg-1]
    self.e2m        =  -1.                   # 2m vapor pressure [Pa]
    self.esat2m     =  -1.                   # 2m saturated vapor pressure [Pa]
    self.u2m        =  -1.                   # 2m u-wind [m s-1]
    self.v2m        =  -1.                   # 2m v-wind [m s-1]

    self.thetasurf  =  self.input.theta      # surface potential temperature [K]
    self.thetav     =  -1.                   # initial mixed-layer potential temperature [K]
    self.dthetav    =  -1.                   # initial virtual temperature jump at h [K]
    self.thetavsurf =  -1.                   # surface virtual potential temperature [K]
    self.qsurf      =  -1.                   # surface specific humidity [g kg-1]
    self.wthetav    =  -1.                   # surface kinematic virtual heat flux [K m s-1]
    
    self.q          =  self.input.q          # initial mixed-layer specific humidity [kg kg-1]
    self.dq         =  self.input.dq         # initial specific humidity jump at h [kg kg-1]
    self.gammaq     =  self.input.gammaq     # free atmosphere specific humidity lapse rate [kg kg-1 m-1]
    self.advq       =  self.input.advq       # advection of moisture [kg kg-1 s-1]
    self.wq         =  self.input.wq         # surface kinematic moisture flux [kg kg-1 m s-1]

    self.qsat       =  -1.                   # mixed-layer saturated specific humidity [kg kg-1]
    self.esat       =  -1.                   # mixed-layer saturated vapor pressure [Pa]
    self.e          =  -1.                   # mixed-layer vapor pressure [Pa]
    self.qsatsurf   =  -1.                   # surface saturated specific humidity [g kg-1]
    self.dqsatdT    =  -1.                   # slope saturated specific humidity curve [g kg-1 K-1]
    self.rh2m       =  -1.                   # relative humidity at 2 m 
    self.rhtop      =  -1.                   # relative humidity top boundary layer 
    self.lcl        =  -1.                   # lifting condensation level [m]
    
    self.co2        =  self.input.co2        # initial mixed-layer carbon dioxide [ppm]
    self.dco2       =  self.input.dco2       # initial carbon dioxide jump at h [ppm]
    self.gammaco2   =  self.input.gammaco2   # free atmosphere carbon dioxide lapse rate [ppm m-1]
    self.advco2     =  self.input.advco2     # advection of carbon dioxide [ppm s-1]
    self.wco2       =  self.input.wco2       # surface kinematic carbon dioxide flux [ppm m s-1]
    self.awco2      =  self.input.awco2      # surface carbon dioxide flux [mg C m-2 s-1]
    
    self.sw_wind    =  self.input.sw_wind    # prognostic wind switch
    self.u          =  self.input.u          # initial mixed-layer u-wind speed [m s-1]
    self.du         =  self.input.du         # initial u-wind jump at h [m s-1]
    self.gammau     =  self.input.gammau     # free atmosphere u-wind speed lapse rate [s-1]
    self.advu       =  self.input.advu       # advection of u-wind [m s-2]
    
    self.v          =  self.input.v          # initial mixed-layer u-wind speed [m s-1]
    self.dv         =  self.input.dv         # initial u-wind jump at h [m s-1]
    self.gammav     =  self.input.gammav     # free atmosphere v-wind speed lapse rate [s-1]
    self.advv       =  self.input.advv       # advection of v-wind [m s-2]

    self.htend      =  -1.                   # tendency of CBL [m s-1]
    self.thetatend  =  -1.                   # tendency of mixed-layer potential temperature [K s-1]
    self.dthetatend =  -1.                   # tendency of potential temperature jump at h [K s-1]
    self.qtend      =  -1.                   # tendency of mixed-layer specific humidity [kg kg-1 s-1]
    self.dqtend     =  -1.                   # tendency of specific humidity jump at h [kg kg-1 s-1]
    self.co2tend    =  -1.                   # tendency of mixed-layer carbon dioxide [ppm s-1]
    self.dco2tend   =  -1.                   # tendency of carbon dioxide jump at h [ppm s-1]
    self.utend      =  -1.                   # tendency of u-wind [m s-1 s-1]
    self.dutend     =  -1.                   # tendency of u-wind jump at h [m s-1 s-1]
    self.vtend      =  -1.                   # tendency of v-wind [m s-1 s-1]
    self.dvtend     =  -1.                   # tendency of v-wind jump at h [m s-1 s-1]

    # initialize surface layer
    self.sw_sl      =  self.input.sw_sl      # surface layer switch
    self.ustar      =  self.input.ustar      # surface friction velocity [m s-1]
    self.uw         =  -1.                   # surface momentum flux in u-direction [m2 s-2]
    self.vw         =  -1.                   # surface momentum flux in v-direction [m2 s-2]
    self.z0m        =  self.input.z0m        # roughness length for momentum [m]
    self.z0h        =  self.input.z0h        # roughness length for scalars [m]
    self.Cm         =  -1.                   # drag coefficient for momentum [-]
    self.Cs         =  -1.                   # drag coefficient for scalars [-]
    self.L          =  -1.                   # Obukhov length [m]
    self.Rib        =  -1.                   # bulk Richardson number [-]
    self.ra         =  -1.                   # aerodynamic resistance [s m-1]

    # initialize radiation
    self.sw_rad     =  self.input.sw_rad     # radiation switch
    self.lat        =  self.input.lat        # latitude [deg]
    self.lon        =  self.input.lon        # longitude [deg]
    self.doy        =  self.input.doy        # day of the year [-]
    self.tstart     =  self.input.tstart     # time of the day [-]
    self.cc         =  self.input.cc         # cloud cover fraction [-]
    self.Swin       =  -1.                   # incoming short wave radiation [W m-2]
    self.Swout      =  -1.                   # outgoing short wave radiation [W m-2]
    self.Lwin       =  -1.                   # incoming long wave radiation [W m-2]
    self.Lwout      =  -1.                   # outgoing long wave radiation [W m-2]
    self.Q          =  self.input.Q          # net radiation [W m-2]

    # initialize land surface
    self.sw_ls      =  self.input.sw_ls      # land surface switch
    self.wg         =  self.input.wg         # volumetric water content top soil layer [m3 m-3]
    self.w2         =  self.input.w2         # volumetric water content deeper soil layer [m3 m-3]
    self.Tsoil      =  self.input.Tsoil      # temperature top soil layer [K]
    self.T2         =  self.input.T2         # temperature deeper soil layer [K]
                       
    self.a          =  self.input.a          # Clapp and Hornberger retention curve parameter a [-]
    self.b          =  self.input.b          # Clapp and Hornberger retention curve parameter b [-]
    self.p          =  self.input.p          # Clapp and Hornberger retention curve parameter p [-]
    self.CGsat      =  self.input.CGsat      # saturated soil conductivity for heat
                       
    self.wsat       =  self.input.wsat       # saturated volumetric water content ECMWF config [-]
    self.wfc        =  self.input.wfc        # volumetric water content field capacity [-]
    self.wwilt      =  self.input.wwilt      # volumetric water content wilting point [-]
                       
    self.C1sat      =  self.input.C1sat      
    self.C2ref      =  self.input.C2ref      
    
    self.LAI        =  self.input.LAI        # leaf area index [-]
    self.gD         =  self.input.gD         # correction factor transpiration for VPD [-]
    self.rsmin      =  self.input.rsmin      # minimum resistance transpiration [s m-1]
    self.rssoilmin  =  self.input.rssoilmin  # minimum resistance soil evaporation [s m-1]
    self.alpha      =  self.input.alpha      # surface albedo [-]

    self.rs         =  600                  # resistance transpiration [s m-1]
    self.rssoil     =  1.e6                  # resistance soil [s m-1]
                       
    self.Ts         =  self.input.Ts         # surface temperature [K]
                       
    self.cveg       =  self.input.cveg       # vegetation fraction [-]
    self.Wmax       =  self.input.Wmax       # thickness of water layer on wet vegetation [m]
    self.Wl         =  self.input.Wl         # equivalent water layer depth for wet vegetation [m]
    self.cliq       =  -1.                   # wet fraction [-]
                       
    self.Lambda     =  self.input.Lambda     # thermal diffusivity skin layer [-]

    self.Tsoiltend  =  -1.                   # soil temperature tendency [K s-1]
    self.wgtend     =  -1.                   # soil moisture tendency [m3 m-3 s-1]
    self.Wltend     =  -1.                   # equivalent liquid water tendency [m s-1]


    self.H          =  -1.                   # sensible heat flux [W m-2]
    self.LE         =  -1.                   # evapotranspiration [W m-2]
    self.LEliq      =  -1.                   # open water evaporation [W m-2]
    self.LEveg      =  -1.                   # transpiration [W m-2]
    self.LEsoil     =  -1.                   # soil evaporation [W m-2]
    self.LEpot      =  -1.                   # potential evaporation [W m-2]
    self.LEref      =  -1.                   # reference evaporation using rs = rsmin / LAI [W m-2]
    self.G          =  -1.                   # ground heat flux [W m-2]
    
    # initialize plant physilogical model (A-gs)
    self.ci         =  -1.                   # CO2 concentration inside the leaf [mg m-3] 
    self.cfrac      =  -1.                   # CO2 concentration fraction [-] 
    self.Ds         =  -1.                   # vapor pressure deficit [kPa]
    self.D0         =  -1.                   # vapor pressure deficit stomata closes [kPa]
    self.gm         =  -1.                   # mesophyll conducatnce [mm s-1] 
    self.fmin       =  -1.                   # minimum value cfrac [-]
    self.fmin0      =  -1.                   # function to calculate fmin [-]
    self.Ammax      =  -1.                   # CO2 maximal primary productivity [mg m-2 s-1] 
    self.Am         =  -1.                   # CO2 primray productivity [mg m-2 s-1] 
    self.An         =  -1.                   # net CO2 flow into the plant [mg m-2 s-1] 
    self.Rdark      =  -1.                   # CO2 dark respiration [mg m-2 s-1]
    self.PAR        =  -1.                   # Photosyntetically Active Radiation [W m-2] 
    self.gcCo2      =  -1.                   # CO2 conductance at canopy level [mm s-1] 
    self.rsAgs      =  -1.                   # surface resistance moisture [s mm-1] 
    self.rsCO2      =  -1.                   # surface resistance carbon dioxide [s mm-1] 
    self.betaw      =  -1.                   # function depending on soil moisture content to calculate stress function [-] 
    self.fstr       =  -1.                   # stress function included in canopy conductance [-] 
    
    # initialize constants depending C3 or C4 plants
    self.CO2comp298 =  self.input.CO2comp298 # CO2 compensation concentration [mg m-3]
    self.Q10CO2     =  self.input.Q10CO2     # function parameter to calculate CO2 compensation concentration [-]
    self.gm298      =  self.input.gm298      # mesophyill conductance at 298 K [mm s-1]
    self.Ammax298   =  self.input.Ammax298   # CO2 maximal primary productivity [mg m-2 s-1]
    self.Q10gm      =  self.input.Q10gm      # function parameter to calculate mesophyll conductance [-]
    self.T1gm       =  self.input.T1gm       # reference temperature to calculate mesophyll conductance gm [K]
    self.T2gm       =  self.input.T2gm       # reference temperature to calculate mesophyll conductance gm [K]
    self.Q10Am      =  self.input.Q10Am      # function parameter to calculate maximal primary profuctivity Ammax
    self.T1Am       =  self.input.T1Am       # reference temperature to calculate maximal primary profuctivity Ammax [K]
    self.T2Am       =  self.input.T2Am       # reference temperature to calculate maximal primary profuctivity Ammax [K]
    self.f0         =  self.input.f0         # maximum value Cfrac [-]
    self.ad         =  self.input.ad         # regression coefficient to calculate Cfrac [kPa-1]
    self.alpha0     =  self.input.alpha0     # initial low light conditions [mg J-1]
    self.frveg      =  self.input.frveg      # fraction of the shortwve radiation contributing to PAR [-]
    self.Kx         =  self.input.Kx         # extinction coefficient PAR [-]
    self.gmin       =  self.input.gmin       # cuticular (minimum) conductance [mm s-1]

    # initialize soil  -1. ration model (coupled to A-gs) 
    self.fw         =  -1.                   # water stress correction function [-]
    self.Resp       =  -1.                   # soil surface carbon dioxide flux [mg m-2 s-1] 

    # initialize constants soil

    
    # initialize time variables
    self.tsteps = int(numpy.floor(self.input.runtime / self.input.dt))
    self.dt     = self.input.dt
    self.t      = 0

    # initialize output
    self.out = modeloutput(self.tsteps)

    # calculate initial diagnostic variables
    if(self.sw_rad):
      self.runradmodel()

    if(self.sw_sl):
      self.runslmodel()

    if(self.sw_ls):
      self.runlsmodel()
    
    if(self.sw_ml):
      self.runmlmodel()

    # store initial values in output
    self.store()

  def timestep(self):
    # run radiation model
    if(self.sw_rad):
      self.runradmodel()

    # run surface layer model
    if(self.sw_sl):
      self.runslmodel()
    
    # run land surface model
    if(self.sw_ls):
      self.runlsmodel()

    # run mixed-layer model
    if(self.sw_ml):
      self.runmlmodel()

    # time integrate land surface model
    if(self.sw_ls):
      self.intlsmodel()

    # time integrate mixed-layer model
    if(self.sw_ml):
      self.intmlmodel()
   
  def runmlmodel(self):
    if(not self.sw_sl):
      # decompose ustar along the wind components
      self.uw       = - (self.ustar ** 4. / (self.v ** 2. / self.u ** 2. + 1.)) ** (0.5)
      self.vw       = - (self.ustar ** 4. / (self.u ** 2. / self.v ** 2. + 1.)) ** (0.5)
    
    # compute mixed-layer tendencies
    # first compute necessary virtual temperature units
    self.thetav   = self.theta  + 0.61 * self.theta * self.q
    self.wthetav  = self.wtheta + 0.61 * self.theta * self.wq
    self.dthetav  = (self.theta + self.dtheta) * (1. + 0.61 * (self.q + self.dq)) - self.theta * (1. + 0.61 * self.q)
    
    # compute tendencies
    #self.we    = (self.beta * self.wthetav) / self.dthetav
    self.we     = (self.beta * self.wthetav + 5. * self.ustar ** 3. * self.thetav / (self.g * self.h)) / self.dthetav

    self.htend       = self.we + self.ws
    
#    tutc = 0.
#    tutc          = self.t * self.dt / 3600. + self.tstart
#    if(self.t < 252.):
#      self.advtheta = -1.00e-04 
#      self.advq     = -4.17e-08*max(1. - 0.18 * (tutc - 6.), 0.)
#    else:
#      tutc          = self.t * self.dt / 3600. + self.tstart
#      self.advtheta = -1.00e-04*max(1. - 0.36 * (tutc - 6.), -1.)
#      self.advq     =  0. 

#    print(self.t,tutc,self.advtheta)
    
    self.thetatend   = (self.wtheta + self.we * self.dtheta) / self.h + self.advtheta 
    self.qtend       = (self.wq     + self.we * self.dq)     / self.h + self.advq
    self.co2tend     = (self.wco2   + self.we * self.dco2)   / self.h + self.advco2
    
    self.dthetatend  = self.gammatheta * self.we - self.thetatend
    self.dqtend      = self.gammaq     * self.we - self.qtend
    self.dco2tend    = self.gammaco2   * self.we - self.co2tend
   
    # assume u + du = ug, so ug - u = du
    if(self.sw_wind):
      self.utend       = -self.fc * self.dv + (self.uw + self.we * self.du)  / self.h + self.advu
      #self.utend       = (self.uw + self.we * self.du)  / self.h + self.advu
      self.vtend       =  self.fc * self.du + (self.vw + self.we * self.dv)  / self.h + self.advv
      #self.vtend       =  0. 

      self.dutend      = self.gammau * self.we - self.utend
      self.dvtend      = self.gammav * self.we - self.vtend
 
  def intmlmodel(self):
    # set values previous time step
    h0      = self.h
    
    theta0  = self.theta
    dtheta0 = self.dtheta
    q0      = self.q
    dq0     = self.dq
    co20    = self.co2
    dco20   = self.dco2
    
    u0      = self.u
    du0     = self.du
    v0      = self.v
    dv0     = self.dv

    # integrate mixed-layer equations
    self.h        = h0      + self.dt * self.htend

    self.theta    = theta0  + self.dt * self.thetatend
    self.dtheta   = dtheta0 + self.dt * self.dthetatend
    self.q        = q0      + self.dt * self.qtend
    self.dq       = dq0     + self.dt * self.dqtend
    self.co2      = co20    + self.dt * self.co2tend
    self.dco2     = dco20   + self.dt * self.dco2tend

    if(self.sw_wind):
      self.u        = u0      + self.dt * self.utend
      self.du       = du0     + self.dt * self.dutend
      self.v        = v0      + self.dt * self.vtend
      self.dv       = dv0     + self.dt * self.dvtend

    
    if(self.h > 900. and self.h < 950.):
      self.gammatheta = 0.001 + (self.h - 900.) / 50. * 0.01

    if(self.h > 950.):
      self.gammatheta = 0.015
  
#    print(self.t,self.h,self.gammatheta)
 
  def runradmodel(self):
    sda    = 0.409 * numpy.cos(2. * numpy.pi * (self.doy - 173.) / 365.)
    sinlea = numpy.sin(2. * numpy.pi * self.lat / 360.) * numpy.sin(sda) - numpy.cos(2. * numpy.pi * self.lat / 360.) * numpy.cos(sda) * numpy.cos(2. * numpy.pi * (self.t * self.dt + self.tstart * 3600.) / 86400. - 2. * numpy.pi * self.lon / 360.)
    sinlea = max(sinlea, 0.0001)
    
    Ta  = self.theta * ((self.Ps - 0.1 * self.h * self.rho * self.g) / self.Ps ) ** (self.Rd / self.cp)

    Tr  = (0.6 + 0.2 * sinlea) * (1. - 0.4 * self.cc)

    #self.Swin  = self.S0 * numpy.sin(sinlea)
    #self.Swout = self.alpha * self.S0 * numpy.sin(sinlea)
    
    self.Swin  = self.S0 * Tr * sinlea
    self.Swout = self.alpha * self.S0 * Tr * sinlea
    self.Lwin  = 0.8 * self.bolz * Ta ** 4.
    self.Lwout = self.bolz * self.Ts ** 4.
      
    self.Q     = self.Swin - self.Swout + self.Lwin - self.Lwout

  def runslmodel(self):

    ueff           = numpy.sqrt(self.u ** 2. + self.v ** 2.)
    self.thetasurf = self.theta + self.wtheta / (self.Cs * ueff)
    esatsurf       = 0.611e3 * numpy.exp(17.2694 * (self.thetasurf - 273.16) / (self.thetasurf - 35.86))
    qsatsurf       = 0.622 * esatsurf / self.Ps
    cq             = (1. + self.Cs * ueff * self.rs) ** -1.
    self.qsurf     = (1. - cq) * self.q + cq * qsatsurf

    self.thetavsurf = self.thetasurf * (1. + 0.61 * self.qsurf)
    self.thetav     = self.theta * (1. + 0.61 * self.q)
    self.wthetav    = self.wtheta + 0.61 * self.theta * self.wq

    zsl  = 0.1 * self.h
          
    self.Rib  = self.g / self.thetav * 0.1 * self.h * (self.thetav - self.thetavsurf) / (self.u ** 2. + self.v ** 2.)
    self.Rib  = min(self.Rib, 0.2)

    if(self.Rib > 0.):
      L    = 1.
      L0   = 2.
    else:
      L  = -1.
      L0 = -2.
    
    while (abs(L - L0) > 0.001):
      L0      = L
      fx      = self.Rib - zsl / L * (numpy.log(zsl / self.z0h) - self.psih(zsl / L) + self.psih(self.z0h / L)) / (numpy.log(zsl / self.z0m) - self.psim(zsl / L) + self.psim(self.z0m / L)) ** 2.
      Lstart  = L - 0.001*L
      Lend    = L + 0.001*L
      fxdif   = ( (- zsl / Lstart * (numpy.log(zsl / self.z0h) - self.psih(zsl / Lstart) + self.psih(self.z0h / Lstart)) / (numpy.log(zsl / self.z0m) - self.psim(zsl / Lstart) + self.psim(self.z0m / Lstart)) ** 2.) \
          - (-zsl / Lend * (numpy.log(zsl / self.z0h) - self.psih(zsl / Lend) + self.psih(self.z0h / Lend)) / (numpy.log(zsl / self.z0m) - self.psim(zsl / Lend) + self.psim(self.z0m / Lend)) ** 2.) ) / (Lstart - Lend)
      L       = L - fx / fxdif

    self.L = L
    #self.L    = ribtol.ribtol(self.Rib, zsl, self.z0m, self.z0h)

    self.Cm   = self.k ** 2. / (numpy.log(zsl / self.z0m) - self.psim(zsl / self.L) + self.psim(self.z0m / self.L)) ** 2.
    self.Cs   = self.k ** 2. / (numpy.log(zsl / self.z0m) - self.psim(zsl / self.L) + self.psim(self.z0m / self.L)) / (numpy.log(zsl / self.z0h) - self.psih(zsl / self.L) + self.psih(self.z0h / self.L))

      
    #if(wthetav > 0.):
    #  wstar     = (g / thetav * h * wthetav) ** (1./3.)
    #else:
    #  wstar     = 0.0001
    #ueff      = numpy.sqrt(u ** 2. + wstar ** 2.)
    
    self.ustar = numpy.sqrt(self.Cm) * ueff
    self.uw    = - self.Cm * ueff * self.u
    self.vw    = - self.Cm * ueff * self.v

    # diagnostic meteorological variables
    self.T2m    = self.thetasurf - self.wtheta / self.ustar / self.k * (numpy.log(2. / self.z0h) - self.psih(2. / self.L) + self.psih(self.z0h / self.L))
    self.q2m    = self.qsurf     - self.wq     / self.ustar / self.k * (numpy.log(2. / self.z0h) - self.psih(2. / self.L) + self.psih(self.z0h / self.L))
    self.u2m    =                - self.uw     / self.ustar / self.k * (numpy.log(2. / self.z0m) - self.psim(2. / self.L) + self.psim(self.z0m / self.L))
    self.v2m    =                - self.vw     / self.ustar / self.k * (numpy.log(2. / self.z0m) - self.psim(2. / self.L) + self.psim(self.z0m / self.L))
    self.esat2m = 0.611e3 * numpy.exp(17.2694 * (self.T2m - 273.16) / (self.T2m - 35.86))
    self.e2m    = self.q2m * self.Ps / 0.622
    self.rh2m   = self.e2m / self.esat2m
    
    # Calculating RH at the top of the boundary layer
    Ptop        = self.Ps / numpy.exp((self.g* self.h) / (self.Rd * self.theta))
    Ttop        = self.theta * (Ptop / self.Ps ) ** (self.Rd / self.cp)
    esattop     = 0.611e3 * numpy.exp(17.2694 * (Ttop - 273.16) / (Ttop - 35.86)) 
    etop        = self.q * Ptop/ 0.622 
    self.rhtop  = etop / esattop

    # Calculating the lifting condenstaion level (Bolton 2008)
    Tdew        = 1. / ((1./273.16) - (self.Rv/self.Lv) * numpy.log(self.e2m/611))
    Tlcl        = 1. / ((1./(Tdew - 56.)) + (numpy.log(self.theta/Tdew)/800.)) + 56
    self.lcl    = 0. - (self.cp * (Tlcl-self.theta) / self.g) 

  def psim(self, zeta):
    if(zeta <= 0):
      #x     = (1. - 16. * zeta) ** (0.25)
      #psim  = 3.14159265 / 2. - 2. * numpy.arctan(x) + numpy.log( (1.+x) ** 2. * (1. + x ** 2.) / 8.)
      x     = (1. + 3.6 * abs(zeta) ** (2./3.)) ** (-0.5)
      psim = 3. * numpy.log( (1. + 1. / x) / 2.)
    else:
      psim  = -2./3. * (zeta - 5./0.35) * numpy.exp(-0.35 * zeta) - zeta - (10./3.) / 0.35
    return psim
    
  def psih(self, zeta):
    if(zeta <= 0):
      #x     = (1. - 16. * zeta) ** (0.25)
      #psih  = 2. * numpy.log( (1. + x ** 2.) / 2. )
      x     = (1. + 7.9 * abs(zeta) ** (2./3.)) ** (-0.5)
      psih  = 3. * numpy.log( (1. + 1. / x) / 2.)
    else:
      psih  = -2./3. * (zeta - 5./0.35) * numpy.exp(-0.35 * zeta) - (1. + (2./3.) * zeta) ** (1.5) - (10./3.) / 0.35 + 1.
    return psih

  def factorial(self, k):
    factorial = 1
    for n in range(2, k + 1):
      factorial = factorial * n
    return factorial

  def E1(self, x):
    E1sum = 0.

    for k in range(1, 100):
      E1sum = E1sum + (-1.) ** (k + 0.0) * x ** (k + 0.0) / ( (k + 0.0) * self.factorial(k) )

    E1 = -0.57721566490153286060 - numpy.log(x) - E1sum
    #Eisum = 0.

    #for k in range(0, 100):
      #Eisum + factorial(k) / (-x) ** k
      #print(factorial(k) / (-x) ** k)

    #Ei = numpy.exp(-x) / x * Eisum

    return E1

  def runlsmodel(self):
    # compute ra
    ueff       = numpy.sqrt(self.u ** 2. + self.v ** 2.)
    if(self.sw_sl):
      self.ra    = (self.Cm * ueff) ** (-1.)
    else:
      self.ra    = ueff / self.ustar ** 2.

    # first calculate essential thermodynamic variables
    self.esat    = 0.611e3 * numpy.exp(17.2694 * (self.theta - 273.16) / (self.theta - 35.86))
    self.qsat    = 0.622 * self.esat / self.Ps
    desatdT      = self.esat * (17.2694 / (self.theta - 35.86) - 17.2694 * (self.theta - 273.16) / (self.theta - 35.86)**2.)
    self.dqsatdT = 0.622 * desatdT / self.Ps
    self.e = self.q * self.Ps / 0.622

    # calculate surface resistances using Jarvis-Stewart model
    if(self.sw_rad):
      f1          = 1. / ((0.004 * self.Swin + 0.05) / (0.81 * (0.004 * self.Swin + 1.)))
    else:
      f1          = 1.

    #fpar        = 0.55 * self.Swin / 100. * 2. / self.LAI
    #f1new       = (1. + fpar) / (fpar + self.rsmin / 10000.)
    if(self.w2 > self.wwilt):# and self.w2 <= self.wfc):
      f2        = (self.wfc - self.wwilt) / (self.w2 - self.wwilt)
    else:
      f2        = 1.e8

    f3          = 1. / numpy.exp(- self.gD * (self.esat2m - self.e2m) / 100.)

    f4          = 1./ (1. - 0.0016 * (298.0 - self.T2m) ** 2.)

    self.rs     = self.rsmin / self.LAI * f1 * f2 * f3
  
    sw_rsAgs    = True 
 
    if(sw_rsAgs): 
      # calculate surface resistances using plant physiological (A-gs) model
        # calculate CO2 compensation concentration
      self.CO2comp = self.CO2comp298 * self.Q10CO2 ** ( 0.1 * (self.thetasurf - 298.) )    # CO2 compensation concentration
      self.CO2comp = self.CO2comp * self.rho
      #self.CO2comp = self.CO2comp*(self.mco2/self.mair)*self.rho                           # conversion mumol mol-1 (ppm) to mgCO2 m3 

        # calculate mesophyll conductance 
      self.gm    = self.gm298 *  self.Q10gm ** ( 0.1 * (self.thetasurf - 298.) ) / ( (1. + numpy.exp(0.3 * (self.T1gm - self.thetasurf))) * (1. + numpy.exp(0.3 * (self.thetasurf - self.T2gm)))) # mesophyill conductance
      self.gm    = self.gm / 1000.   # conversion from mm s-1 to m s-1 

      #print(self.t,self.thetasurf,self.CO2comp,self.gm)
        # calculate CO2 concentration inside the leaf (ci)
      #self.fmin0 = self.gmin/self.nuco2q - (1. / 9.) * self.gm * 1000.    
      self.fmin0 = self.gmin/self.nuco2q - 1. / 9. * self.gm
      #self.fmin  = (-self.fmin0 + ( self.fmin0 ** 2. + 4 * self.gmin/self.nuco2q * self.gm * 1000. ) ** (0.5)) / (2. * self.gm * 1000.) 
      self.fmin  = (-self.fmin0 + ( self.fmin0 ** 2. + 4 * self.gmin/self.nuco2q * self.gm ) ** (0.5)) / (2. * self.gm) 
      #print(self.t,self.gmin/self.nuco2q,1. / 9. * self.gm,self.fmin0,self.fmin)
      #print(self.t,self.gm,self.fmin0,self.fmin)
      
      esatsurf   = 0.611e3 * numpy.exp(17.2694 * (self.Ts - 273.16) / (self.Ts - 35.86))
      self.Ds    = (esatsurf - self.e)   / 1000.     # in kPa
      self.D0    = (self.f0 - self.fmin) / self.ad
      #print(self.t,esatsurf, self.e, self.Ds,self.D0)

      self.cfrac = self.f0 * (1. - self.Ds / self.D0) + self.fmin * (self.Ds / self.D0)
      co2abs     = self.co2*(self.mco2/self.mair)*self.rho                                                   # conversion mumol mol-1 (ppm) to mgCO2 m3
      self.ci    = self.cfrac * (co2abs - self.CO2comp) + self.CO2comp          
      #print(self.t,self.cfrac, co2abs, self.ci)

        # calculate maximal gross primary production in high light conditions (Ag)
      self.Ammax = self.Ammax298 *  self.Q10Am ** ( 0.1 * (self.thetasurf - 298.) ) / ( (1. + numpy.exp(0.3 * (self.T1Am - self.thetasurf))) * (1. + numpy.exp(0.3 * (self.thetasurf - self.T2Am))))
        
      #print(self.t,self.thetasurf,self.CO2comp,self.gm,self.Ammax)
        # calculate effect of soil moisture stress on gross assimilation rate
      self.betaw = max(1e-3, min(1,(self.wg - self.wwilt)/(self.wfc - self.wwilt)))
      #print(self.t,self.betaw)
      #self.betaw = 0.26 

        # calculate stress function
      #self.fstr  = 2.*self.betaw - self.betaw**2.
      self.fstr  = self.betaw

        # calculate gross assimilation rate (Am) 
      self.Am    = self.Ammax * ( 1. - numpy.exp( -(self.gm * (self.ci - self.CO2comp) / self.Ammax) ) )
      #print(self.t,self.Am)

      self.Rdark = (1./ 9.) * self.Am
 
      #self.PAR   = max(1e-1,self.Swin*self.cveg)
      self.PAR   = 0.40*max(1e-1,self.Swin*self.cveg)

        # calculate  light use efficiency
      self.alphac= self.alpha0 * (co2abs - self.CO2comp) / (co2abs + 2. * self.CO2comp)
        
        # calculate gross primary productivity 
      self.Ag    = (self.Am + self.Rdark)*(1 - numpy.exp(self.alphac*self.PAR/(self.Am + self.Rdark)))

        # 1.- calculate upscaling from leaf to canopy: net flow CO2 into the plant (An)

      y       =  self.alphac * self.Kx * self.PAR / (self.Am + self.Rdark)
      self.An = (self.Am + self.Rdark) * (1. - 1. / (self.Kx * self.LAI) * (self.E1( y * numpy.exp(-self.Kx * self.LAI)) - self.E1(y)))
      #print(self.t,y, self.An)

        # 2.- calculate upscaling from leaf to canopy: CO2 conductance at canopy level 

      a1    = 1. / (1. - self.f0)
      Dstar = self.D0 / (a1 - 1.)

      self.gcco2 = self.LAI * (self.gmin/self.nuco2q + a1 * self.fstr * self.An / ((co2abs - self.CO2comp) * (1. + self.Ds / Dstar)))
      #print(self.t,a1 * self.fstr * self.An, self.gcco2)

        # calculate surface resistance for moisture and carbon dioxide

      self.rsAgs = 1. / (1.6 * self.gcco2)
      #print(self.t, self.rsAgs)
      self.rsCO2 = 1. /self.gcco2

        # calculate net flux of CO2 into the plant (An)

      self.An    = -(co2abs - self.ci) / (self.ra + self.rsAgs)
      #print(self.t, self.An)

      # CO2 soil respiration surface flux

      self.fw     = self.Cw * self.wmax / (self.wg + self.wmin)
      self.Resp   = self.R10 * (1. - self.fw) * numpy.exp( self.E0 / (283.15 * 8.314) * (1. - 283.15 / (self.thetasurf)))  
      #self.Resp   = self.R10 * (1. - self.fw) * numpy.exp( self.E0 / (283.15 * 8.314) * (1. - 283.15 / (self.Tsoil)))  

      # CO2 flux

#      print(self.An,self.Resp)
      self.awco2  = (self.An + self.Resp)                                      # conversion mgCO2 m3 to mumol mol-1 (ppm)
      self.wco2   = (self.An + self.Resp)*(self.mair/self.mco2)*(1./self.rho)  # conversion mgCO2 m3 to mumol mol-1 (ppm)
      #print(self.t*self.dt / 3600. + self.tstart, self.An, self.Resp, self.awco2)

      # replace Jarvis rs with Ags rs
      self.rs     = self.rsAgs

    # recompute f2 using wg instead of w2
    if(self.wg > self.wwilt):# and self.w2 <= self.wfc):
      f2          = (self.wfc - self.wwilt) / (self.wg - self.wwilt)
    else:
      f2        = 1.e8
    self.rssoil = self.rssoilmin * f2 

    Wlmx = self.LAI * self.Wmax
    self.cliq = min(1., self.Wl / Wlmx) 
   
    # calculate skin temperature implictly
    self.Ts   = (self.Q  + self.rho * self.cp / self.ra * self.theta \
        + self.cveg * (1. - self.cliq) * self.rho * self.Lv / (self.ra + self.rs) * (self.dqsatdT * self.theta - self.qsat + self.q) \
        + (1. - self.cveg) * self.rho * self.Lv / (self.ra + self.rssoil) * (self.dqsatdT * self.theta - self.qsat + self.q) \
        + self.cveg * self.cliq * self.rho * self.Lv / self.ra * (self.dqsatdT * self.theta - self.qsat + self.q) + self.Lambda * self.Tsoil) \
      * (self.rho * self.cp / self.ra + self.cveg * (1. - self.cliq) * self.rho * self.Lv / (self.ra + self.rs) * self.dqsatdT + (1. - self.cveg) * self.rho * self.Lv / (self.ra + self.rssoil) * self.dqsatdT + self.cveg * self.cliq * self.rho * self.Lv / self.ra * self.dqsatdT + self.Lambda) ** (-1.)

    esatsurf  = 0.611e3 * numpy.exp(17.2694 * (self.Ts - 273.16) / (self.Ts - 35.86))
    self.qsatsurf  = 0.622 * esatsurf / self.Ps
    
    self.LEveg  = (1. - self.cliq) * self.cveg * self.rho * self.Lv / (self.ra + self.rs) * (self.dqsatdT * (self.Ts - self.theta) + self.qsat - self.q)
    self.LEliq  = self.cliq * self.cveg * self.rho * self.Lv / self.ra * (self.dqsatdT * (self.Ts - self.theta) + self.qsat - self.q)
    self.LEsoil = (1. - self.cveg) * self.rho * self.Lv / (self.ra + self.rssoil) * (self.dqsatdT * (self.Ts - self.theta) + self.qsat - self.q)

    self.Wltend      = - self.LEliq / (self.rhow * self.Lv)

    self.LE     = self.LEsoil + self.LEveg + self.LEliq
    self.H      = self.rho * self.cp / self.ra * (self.Ts - self.theta)
    self.G      = self.Lambda * (self.Ts - self.Tsoil)
    self.LEpot  = (self.dqsatdT * (self.Q - self.G) + self.rho * self.cp / self.ra * (self.qsat - self.q)) / (self.dqsatdT + self.cp / self.Lv)
    self.LEref  = (self.dqsatdT * (self.Q - self.G) + self.rho * self.cp / self.ra * (self.qsat - self.q)) / (self.dqsatdT + self.cp / self.Lv * (1. + self.rsmin / self.LAI / self.ra))
    
    CG          = self.CGsat * (self.wsat / self.w2) ** (self.b / (2. * numpy.log(10.)))

    self.Tsoiltend   = CG * self.G - 2. * numpy.pi / 86400. * (self.Tsoil - self.T2)
 
    d1          = 0.1
    C1          = self.C1sat * (self.wsat / self.wg) ** (self.b / 2. + 1.)
    C2          = self.C2ref * (self.w2 / (self.wsat - self.w2) )
    wgeq        = self.w2 - self.wsat * self.a * ( (self.w2 / self.wsat) ** self.p * (1. - (self.w2 / self.wsat) ** (8. * self.p)) )
    self.wgtend = - C1 / (self.rhow * d1) * self.LEsoil / self.Lv - C2 / 86400. * (self.wg - wgeq)

    # calculate kinematic heat fluxes
    self.wtheta   = self.H  / (self.rho * self.cp)
    self.wq       = self.LE / (self.rho * self.Lv)

  def intlsmodel(self):
    # integrate soil equations
    Tsoil0        = self.Tsoil
    wg0           = self.wg
    Wl0           = self.Wl

    self.Tsoil    = Tsoil0  + self.dt * self.Tsoiltend
    self.wg       = wg0     + self.dt * self.wgtend
    self.Wl       = Wl0     + self.dt * self.Wltend


  # store model output
  def store(self):
    t                      = self.t
    self.out.t[t]          = t * self.dt / 3600. + self.tstart

    self.out.h[t]          = self.h
    self.out.ws[t]         = self.ws
    
    self.out.theta[t]      = self.theta
    self.out.thetav[t]     = self.thetav
    self.out.dtheta[t]     = self.dtheta
    self.out.dthetav[t]    = self.dthetav
    self.out.gammatheta[t] = self.gammatheta
    self.out.advtheta[t]   = self.advtheta
    self.out.beta[t]       = self.beta
    self.out.wtheta[t]     = self.wtheta
    self.out.wthetav[t]    = self.wthetav
    
    self.out.q[t]          = self.q
    self.out.qsat[t]       = self.qsat
    self.out.e[t]          = self.e
    self.out.esat[t]       = self.esat
    self.out.dq[t]         = self.dq
    self.out.gammaq[t]     = self.gammaq
    self.out.advq[t]       = self.advq
    self.out.wq[t]         = self.wq
    
    self.out.co2[t]        = self.co2
    self.out.dco2[t]       = self.dco2
    self.out.gammaco2[t]   = self.gammaco2
    self.out.advco2[t]     = self.advco2
    self.out.wco2[t]       = self.wco2
    self.out.awco2[t]      = self.awco2
    
    self.out.u[t]          = self.u
    self.out.du[t]         = self.du
    self.out.gammau[t]     = self.gammau
    self.out.advu[t]       = self.advu
    
    self.out.v[t]          = self.v
    self.out.dv[t]         = self.dv
    self.out.gammav[t]     = self.gammav
    self.out.advv[t]       = self.advv
    
    self.out.T2m[t]        = self.T2m
    self.out.q2m[t]        = self.q2m
    self.out.u2m[t]        = self.u2m
    self.out.v2m[t]        = self.v2m
    self.out.e2m[t]        = self.e2m
    self.out.esat2m[t]     = self.esat2m
    self.out.rh2m[t]       = self.rh2m
    self.out.rhtop[t]      = self.rhtop
    self.out.lcl[t]        = self.lcl
    
    self.out.thetasurf[t]  = self.thetasurf
    self.out.thetavsurf[t] = self.thetavsurf
    self.out.qsurf[t]      = self.qsurf
    self.out.ustar[t]      = self.ustar
    self.out.Cm[t]         = self.Cm
    self.out.Cs[t]         = self.Cs
    self.out.L[t]          = self.L
    self.out.Rib[t]        = self.Rib

    self.out.Swin[t]       = self.Swin
    self.out.Swout[t]      = self.Swout
    self.out.Lwin[t]       = self.Lwin
    self.out.Lwout[t]      = self.Lwout
    self.out.Q[t]          = self.Q

    self.out.ra[t]         = self.ra
    self.out.rs[t]         = self.rs
    self.out.H[t]          = self.H
    self.out.LE[t]         = self.LE
    self.out.LEliq[t]      = self.LEliq
    self.out.LEveg[t]      = self.LEveg
    self.out.LEsoil[t]     = self.LEsoil
    self.out.LEpot[t]      = self.LEpot
    self.out.LEref[t]      = self.LEref
    self.out.G[t]          = self.G

    #self.out.ci[t]         = self.ci 
    #self.out.Am[t]         = self.Am 
    self.out.An[t]         = self.An 
    self.out.Rdark[t]      = self.Rdark 
    self.out.Resp[t]       = self.Resp 
    #self.out.gcco2[t]      = self.gcco2
    self.out.rsAgs[t]      = self.rsAgs
    self.out.rsCO2[t]      = self.rsCO2

  # delete class variables to facilitate analysis in ipython
  def exitmodel(self):
    del(self.Lv)
    del(self.cp)
    del(self.rho)
    del(self.k)
    del(self.g)
    del(self.Rd)
    del(self.Rv)
    del(self.bolz)
    del(self.S0)
    del(self.rhow)

    del(self.t)
    del(self.dt)
    del(self.tsteps)
     
    del(self.h)          
    del(self.Ps)        
    del(self.fc)        
    del(self.ws)
    del(self.we)
    
    del(self.theta)
    del(self.dtheta)
    del(self.gammatheta)
    del(self.advtheta)
    del(self.beta)
    del(self.wtheta)
  
    del(self.T2m)
    del(self.q2m)
    del(self.e2m)
    del(self.esat2m)
    del(self.u2m)
    del(self.v2m)
    del(self.rh2m)
    del(self.rhtop)
    del(self.lcl)
    
    del(self.thetasurf)
    del(self.qsatsurf)
    del(self.thetav)
    del(self.dthetav)
    del(self.thetavsurf)
    del(self.qsurf)
    del(self.wthetav)
    
    del(self.q)
    del(self.qsat)
    del(self.dqsatdT)
    del(self.e)
    del(self.esat)
    del(self.dq)
    del(self.gammaq)
    del(self.advq)
    del(self.wq)
    
    del(self.co2)
    del(self.dco2)
    del(self.gammaco2)
    del(self.advco2)
    del(self.wco2)
    del(self.awco2)
    
    del(self.u)
    del(self.du)
    del(self.gammau)
    del(self.advu)
    
    del(self.v)
    del(self.dv)
    del(self.gammav)
    del(self.advv)

    del(self.htend)
    del(self.thetatend)
    del(self.dthetatend)
    del(self.qtend)
    del(self.dqtend)
    del(self.utend)
    del(self.dutend)
    del(self.vtend)
    del(self.dvtend)
   
    del(self.Tsoiltend) 
    del(self.wgtend)  
    del(self.Wltend) 

    del(self.ustar)
    del(self.uw)
    del(self.vw)
    del(self.z0m)
    del(self.z0h)        
    del(self.Cm)         
    del(self.Cs)
    del(self.L)
    del(self.Rib)
    del(self.ra)

    del(self.lat)
    del(self.lon)
    del(self.doy)
    del(self.tstart)
 
    del(self.Swin)
    del(self.Swout)
    del(self.Lwin)
    del(self.Lwout)
    del(self.cc)

    del(self.wg)
    del(self.w2)
    del(self.cveg)
    del(self.cliq)
    del(self.Tsoil)
    del(self.T2)
    del(self.a)
    del(self.b)
    del(self.p)
    del(self.CGsat)

    del(self.wsat)
    del(self.wfc)
    del(self.wwilt)

    del(self.C1sat)
    del(self.C2ref)

    del(self.LAI)
    del(self.rs)
    del(self.rssoil)
    del(self.rsmin)
    del(self.rssoilmin)
    del(self.alpha)
    del(self.gD)

    del(self.Ts)

    del(self.Wmax)
    del(self.Wl)

    del(self.Lambda)
    
    del(self.Q)
    del(self.H)
    del(self.LE)
    del(self.LEliq)
    del(self.LEveg)
    del(self.LEsoil)
    del(self.LEpot)
    del(self.LEref)
    del(self.G)

    del(self.sw_ls)
    del(self.sw_rad)
    del(self.sw_sl)
    del(self.sw_wind)


# class for storing mixed-layer model output data
class modeloutput:
  def __init__(self, tsteps):
    self.t          = numpy.zeros(tsteps)    # time [s]

    # mixed-layer variables
    self.h          = numpy.zeros(tsteps)    # initial ABL height [m]
    self.Ps         = numpy.zeros(tsteps)    # surface pressure [Pa]
    self.ws         = numpy.zeros(tsteps)    # large scale vertical velocity [m s-1]
    
    self.theta      = numpy.zeros(tsteps)    # initial mixed-layer potential temperature [K]
    self.thetav     = numpy.zeros(tsteps)    # initial mixed-layer virtual potential temperature [K]
    self.dtheta     = numpy.zeros(tsteps)    # initial potential temperature jump at h [K]
    self.dthetav    = numpy.zeros(tsteps)    # initial virtual potential temperature jump at h [K]
    self.gammatheta = numpy.zeros(tsteps)    # free atmosphere potential temperature lapse rate [K m-1]
    self.advtheta   = numpy.zeros(tsteps)    # advection of heat [K s-1]
    self.beta       = numpy.zeros(tsteps)    # entrainment ratio for virtual heat [-]
    self.wtheta     = numpy.zeros(tsteps)    # surface kinematic heat flux [K m s-1]
    self.wthetav    = numpy.zeros(tsteps)    # surface kinematic virtual heat flux [K m s-1]
    
    self.q          = numpy.zeros(tsteps)    # mixed-layer specific humidity [kg kg-1]
    self.qsat       = numpy.zeros(tsteps)    # mixed-layer saturated specific humidity [kg kg-1]
    self.e          = numpy.zeros(tsteps)    # mixed-layer vapor pressure [Pa]
    self.esat       = numpy.zeros(tsteps)    # mixed-layer saturated vapor pressure [Pa]
    self.dq         = numpy.zeros(tsteps)    # initial specific humidity jump at h [kg kg-1]
    self.gammaq     = numpy.zeros(tsteps)    # free atmosphere specific humidity lapse rate [kg kg-1 m-1]
    self.advq       = numpy.zeros(tsteps)    # advection of moisture [kg kg-1 s-1]
    self.wq         = numpy.zeros(tsteps)    # surface kinematic moisture flux [kg kg-1 m s-1]
    
    self.co2        = numpy.zeros(tsteps)    # initial mixed-layer carbon dioxide [ppm]
    self.dco2       = numpy.zeros(tsteps)    # initial carbon dioxide jump at h [ppm]
    self.gammaco2   = numpy.zeros(tsteps)    # free atmosphere carbon dioxide lapse rate [ppm m-1]
    self.advco2     = numpy.zeros(tsteps)    # advection of carbon dioxide [ppm s-1]
    self.wco2       = numpy.zeros(tsteps)    # surface kinematic carbon dioxide flux [ppm m s-1]
    self.awco2      = numpy.zeros(tsteps)    # surface carbon dioxide flux [mgC m-2 s-1]
    
    self.u          = numpy.zeros(tsteps)    # initial mixed-layer u-wind speed [m s-1]
    self.du         = numpy.zeros(tsteps)    # initial u-wind jump at h [m s-1]
    self.gammau     = numpy.zeros(tsteps)    # free atmosphere u-wind speed lapse rate [s-1]
    self.advu       = numpy.zeros(tsteps)    # advection of u-wind [m s-2]
    
    self.v          = numpy.zeros(tsteps)    # initial mixed-layer u-wind speed [m s-1]
    self.dv         = numpy.zeros(tsteps)    # initial u-wind jump at h [m s-1]
    self.gammav     = numpy.zeros(tsteps)    # free atmosphere v-wind speed lapse rate [s-1]
    self.advv       = numpy.zeros(tsteps)    # advection of v-wind [m s-2]

    # diagnostic meteorological variables
    self.T2m        = numpy.zeros(tsteps)    # 2m temperature [K]   
    self.q2m        = numpy.zeros(tsteps)    # 2m specific humidity [kg kg-1]
    self.u2m        = numpy.zeros(tsteps)    # 2m u-wind [m s-1]    
    self.v2m        = numpy.zeros(tsteps)    # 2m v-wind [m s-1]    
    self.e2m        = numpy.zeros(tsteps)    # 2m vapor pressure [Pa]
    self.esat2m     = numpy.zeros(tsteps)    # 2m saturated vapor pressure [Pa]
    self.rh2m       = numpy.zeros(tsteps)    # relative humidity 2m [-] 
    self.rhtop      = numpy.zeros(tsteps)    # relative humidity top of the boundary layer [-] 
    self.lcl        = numpy.zeros(tsteps)    # lifting condensation level [m] 

    # surface-layer variables
    self.thetasurf  = numpy.zeros(tsteps)    # surface potential temperature [K]
    self.thetavsurf = numpy.zeros(tsteps)    # surface virtual potential temperature [K]
    self.qsurf      = numpy.zeros(tsteps)    # surface specific humidity [kg kg-1]
    self.ustar      = numpy.zeros(tsteps)    # surface friction velocity [m s-1]
    self.z0m        = numpy.zeros(tsteps)    # roughness length for momentum [m]
    self.z0h        = numpy.zeros(tsteps)    # roughness length for scalars [m]
    self.Cm         = numpy.zeros(tsteps)    # drag coefficient for momentum []
    self.Cs         = numpy.zeros(tsteps)    # drag coefficient for scalars []
    self.L          = numpy.zeros(tsteps)    # Obukhov length [m]
    self.Rib        = numpy.zeros(tsteps)    # bulk Richardson number [-]

    # radiation variables
    self.Swin       = numpy.zeros(tsteps)    # incoming short wave radiation [W m-2]
    self.Swout      = numpy.zeros(tsteps)    # outgoing short wave radiation [W m-2]
    self.Lwin       = numpy.zeros(tsteps)    # incoming long wave radiation [W m-2]
    self.Lwout      = numpy.zeros(tsteps)    # outgoing long wave radiation [W m-2]
    self.Q          = numpy.zeros(tsteps)    # net radiation [W m-2]

    # land surface variables
    self.ra         = numpy.zeros(tsteps)    # aerodynamic resistance [s m-1]
    self.rs         = numpy.zeros(tsteps)    # surface resistance [s m-1]
    self.H          = numpy.zeros(tsteps)    # sensible heat flux [W m-2]
    self.LE         = numpy.zeros(tsteps)    # evapotranspiration [W m-2]
    self.LEliq      = numpy.zeros(tsteps)    # open water evaporation [W m-2]
    self.LEveg      = numpy.zeros(tsteps)    # transpiration [W m-2]
    self.LEsoil     = numpy.zeros(tsteps)    # soil evaporation [W m-2]
    self.LEpot      = numpy.zeros(tsteps)    # potential evaporation [W m-2]
    self.LEref      = numpy.zeros(tsteps)    # reference evaporation at rs = rsmin / LAI [W m-2]
    self.G          = numpy.zeros(tsteps)    # ground heat flux [W m-2]

    self.rsAgs      = numpy.zeros(tsteps)    # surface resistance Ags [s m-1]
    self.rsCO2      = numpy.zeros(tsteps)    # surface resistance CO2 [s m-1]
    self.An         = numpy.zeros(tsteps)    # 
    self.Rdark      = numpy.zeros(tsteps)    # 
    self.Resp       = numpy.zeros(tsteps)    # 

# class for storing mixed-layer model input data
class modelinput:
  def __init__(self):
    # general model variables
    self.runtime    = -1. # duration of model run [s]
    self.dt         = -1. # time step [s]

    # mixed-layer variables
    self.sw_ml      = True # mixed-layer model switch
    self.h          = -1. # initial ABL height [m]
    self.Ps         = -1. # surface pressure [Pa]
    self.ws         = -1. # large scale vertical velocity [m s-1]
    self.fc         = -1. # Coriolis parameter [s-1]
    
    self.theta      = -1. # initial mixed-layer potential temperature [K]
    self.dtheta     = -1. # initial temperature jump at h [K]
    self.gammatheta = -1. # free atmosphere potential temperature lapse rate [K m-1]
    self.advtheta   = -1. # advection of heat [K s-1]
    self.beta       = -1. # entrainment ratio for virtual heat [-]
    self.wtheta     = -1. # surface kinematic heat flux [K m s-1]
    
    self.q          = -1. # initial mixed-layer specific humidity [kg kg-1]
    self.dq         = -1. # initial specific humidity jump at h [kg kg-1]
    self.gammaq     = -1. # free atmosphere specific humidity lapse rate [kg kg-1 m-1]
    self.advq       = -1. # advection of moisture [kg kg-1 s-1]
    self.wq         = -1. # surface kinematic moisture flux [kg kg-1 m s-1]
    
    self.co2        = -1. # initial mixed-layer carbon dioxide [ppm]
    self.dco2       = -1. # initial carbon dioxide jump at h [ppm]
    self.gammaco2   = -1. # free atmosphere carbon dioxide lapse rate [ppm m-1]
    self.advco2     = -1. # advection of carbon dioxide [ppm s-1]
    self.wco2       = -1. # surface kinematic carbon dioxide flux [ppm m s-1]
    self.awco2      = -1. # surface carbon dioxide flux [mg C m-2 s-1]
    
    self.sw_wind    = False# prognostic wind switch
    self.u          = -1. # initial mixed-layer u-wind speed [m s-1]
    self.du         = -1. # initial u-wind jump at h [m s-1]
    self.gammau     = -1. # free atmosphere u-wind speed lapse rate [s-1]
    self.advu       = -1. # advection of u-wind [m s-2]

    self.v          = -1. # initial mixed-layer u-wind speed [m s-1]
    self.dv         = -1. # initial u-wind jump at h [m s-1]
    self.gammav     = -1. # free atmosphere v-wind speed lapse rate [s-1]
    self.advv       = -1. # advection of v-wind [m s-2]

    # surface layer variables
    self.sw_sl      = False # surface layer switch
    self.ustar      = -1. # surface friction velocity [m s-1]
    self.z0m        = -1. # roughness length for momentum [m]
    self.z0h        = -1. # roughness length for scalars [m]
    self.Cm         = -1. # drag coefficient for momentum [-]
    self.Cs         = -1. # drag coefficient for scalars [-]
    self.L          = -1. # Obukhov length [-]
    self.Rib        = -1. # bulk Richardson number [-]

    # radiation parameters
    self.sw_rad     = False # radiation switch
    self.lat        = -1. # latitude [deg]
    self.lon        = -1. # longitude [deg]
    self.doy        = -1. # day of the year [-]
    self.tstart     = -1  # time of the day [h UTC]
    self.cc         = -1  # cloud cover fraction [-]

    # land surface parameters
    self.sw_ls      = False # land surface switch
    self.wg         = -1. # volumetric water content top soil layer [m3 m-3]
    self.w2         = -1. # volumetric water content deeper soil layer [m3 m-3]
    self.Tsoil      = -1. # temperature top soil layer [K]
    self.T2         = -1. # temperature deeper soil layer [K]
    
    self.a          = -1. # Clapp and Hornberger retention curve parameter a
    self.b          = -1. # Clapp and Hornberger retention curve parameter b
    self.p          = -1. # Clapp and Hornberger retention curve parameter p 
    self.CGsat      = -1. # saturated soil conductivity for heat
    
    self.wsat       = -1. # saturated volumetric water content ECMWF config [-]
    self.wfc        = -1. # volumetric water content field capacity [-]
    self.wwilt      = -1. # volumetric water content wilting point [-]
    
    self.C1sat      = -1. 
    self.C2ref      = -1.
    
    self.LAI        = -1. # leaf area index [-]
    self.gD         = -1. # correction factor transpiration for VPD [-]
    self.rsmin      = -1. # minimum resistance transpiration [s m-1]
    self.rssoilmin  = -1. # minimum resistance soil evaporation [s m-1]
    self.alpha      = -1. # surface albedo [-]
    
    self.Ts         = -1. # initial surface temperature [K]
    
    self.cveg       = -1. # vegetation fraction [-]
    self.Wmax       = -1. # thickness of water layer on wet vegetation [m]
    self.Wl         = -1. # equivalent water layer depth for wet vegetation [m]
    
    self.Lambda     = -1. # thermal diffusivity skin layer [-]

    # plant physilogical parameters (A-gs)
#    self.ci         = -1. # CO2 concentration inside the leaf [mg m-3] 
#    self.cfrac      = -1. # CO2 concentration fraction [-] 
#    self.Ds         = -1. # vapor pressure deficit [kPa]
#    self.D0         = -1. # vapor pressure deficit stomata closes [kPa]
#    self.gm         = -1. # mesophyll conducatnce [mm s-1] 
#    self.fmin       = -1. # minimum value cfrac [-]
#    self.fmin0      = -1. # function to calculate fmin [-]
#    self.Ammax      = -1. # CO2 maximal primary productivity [mg m-2 s-1] 
#    self.Am         = -1. # CO2 primary productivity [mg m-2 s-1] 
#    self.An         = -1. # net CO2 flow into the plant [mg m-2 s-1] 
#    self.Rdark      = -1. # CO2 dark respiration [mg m-2 s-1]
#    self.PAR        = -1. # Photosyntetically Active Radiation [W m-2] 
#    self.gcCo2      = -1. # CO2 conductance at canopy level [mm s-1] 
#    self.rsAgs      = -1. # surface resistance moisture [s mm-1] 
#    self.rsCO2      = -1. # surface resistance carbon dioxide [s mm-1] 
#    self.betaw      = -1. # function depending on soil moisture content to calculate stress function [-] 
#    self.fstr       = -1. # stress function included in canopy conductance [-] 
    
    # initialize constants depending C3 or C4 plants
    self.CO2comp298 =  -1. # CO2 compensation concentration [mg m-3]
    self.Q10CO2     =  -1. # function parameter to calculate CO2 compensation concentration [-]
    self.gm298      =  -1. # mesophyill conductance at 298 K [mm s-1]
    self.Ammax298   =  -1. # CO2 maximal primary productivity [mg m-2 s-1]
    self.Q10gm      =  -1. # function parameter to calculate mesophyll conductance [-]
    self.T1gm       =  -1. # reference temperature to calculate mesophyll conductance gm [K]
    self.T2gm       =  -1. # reference temperature to calculate mesophyll conductance gm [K]
    self.Q10Am      =  -1. # function parameter to calculate maximal primary profuctivity Ammax
    self.T1Am       =  -1. # reference temperature to calculate maximal primary profuctivity Ammax [K]
    self.T2Am       =  -1. # reference temperature to calculate maximal primary profuctivity Ammax [K]
    self.f0         =  -1. # maximum value Cfrac [-]
    self.ad         =  -1. # regression coefficient to calculate Cfrac [kPa-1]
    self.alpha0     =  -1. # initial low light conditions [mg J-1]
    self.frveg      =  -1. # fraction of the shortwve radiation contributing to PAR [-]
    self.Kx         =  -1. # extinction coefficient PAR [-]
    self.gmin       =  -1. # cuticular (minimum) conductance [mm s-1]

    # soil respiration parameters (coupled to A-gs) 
#    self.fw         = -1. # water stress correction function [-]
#    self.Resp       = -1. # soil surface carbon dioxide flux [mg m-2 s-1] 
