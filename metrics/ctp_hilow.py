'''
---------------------------------------------------------------------------------
  Purpose:
 
  !DESCRIPTION:  
  This subroutine calculates the convective triggering potential (CTP) given early
  morning sounding profiles and also includes a low-level humidity index (Hi_low).
  The CTP returns useful information regarding which profile is primed for convection
  under high surface senisble heat flux and which profiles are likely to trigger
  convection under enhanced latent heat flux. 
  Note that this does not apply to established boundary layer profiles because 
  assumptions are made regarding the presences of inversion.
 
    CTP     =  integral of curve between the moist adiabat and environmental lapse
               rate from 100 mb above the ground to 300mb above the ground
                _       _                _       _
               |         |              |         |
    Hi_low  =  | T - T_d |          -   | T - T_d |
               |_       _|50mb abg      |_       _|150mb abg
 
               where T_d = dew point temperature [K] ;  T = air temperature [K]
                     abg = above ground
 
  
   References:  ** Original Method **
                Findell, K.L., Eltahir, E.A.B. 2003: Atmospheric Controls on 
                Soil Moisture-Boundary Layer Interactions. Part I: Framework
                Development. Journal of Hydrometeorology
 
                Findell, K.L., Eltahir, E.A.B. 2003: Atmospheric Controls on
                Soil Moisture-Boundary Layer Interactions. Part II: Feedbacks
                within the Continental United States. Journal of Hydrometeorology
 
                ** Application and Evaluation **
                Craig R. Ferguson and Eric F. Wood, 2011: Observed Land–Atmosphere
                Coupling from Satellite Remote Sensing and Reanalysis. J. Hydrometeor,
 
  Author and Revision History: 
  Original C code       -- Craig Ferguson 
  Converted to F90 code -- Joshua Roundy on Apr 2015
  Modified by           -- A.B. Tawfik on June 2015
 
 
 ---------------------------------------------------------------------------------
'''
subroutine ctp_hi_low ( nlev_in, tlev_in, qlev_in, plev_in,  &
                          t2m_in , q2m_in , psfc_in, CTP, HILOW , missing )

   implicit none
!
! Input/Output Variables
!
   integer, intent(in )                      ::  nlev_in   ! *** # of pressure levels
   real(4), intent(in )                      ::  missing   ! *** missing value - useful for obs
   real(4), intent(in ), dimension(nlev_in)  ::  tlev_in   ! *** air temperature at pressure levels [K]
   real(4), intent(in ), dimension(nlev_in)  ::  qlev_in   ! *** specific humidity levels [kg/kg]
   real(4), intent(in ), dimension(nlev_in)  ::  plev_in   ! *** pressure levels [Pa]
   real(4), intent(in )                      ::  psfc_in   ! *** surface pressure [Pa]
   real(4), intent(in )                      ::  t2m_in    ! *** 2-meter temperature [K]
   real(4), intent(in )                      ::  q2m_in    ! *** 2-meter specific humidity [kg/kg]
   real(4), intent(out)                      ::  CTP       ! *** Convective Triggering Potential [K]
   real(4), intent(out)                      ::  HILOW     ! *** Low-level humidity [K]
!
! Local variables
!
   real(4), parameter   ::  Rd=287.04, cp=1005.7, R_cp=Rd/cp, C2K = 273.15
   real(4), parameter   ::  Lv=2.5e6 , Lv_cp=Lv/cp, Rv=461.5, grav= 9.81
   real(4), parameter   ::  ep=0.622
   integer, parameter   ::  nsegments = 20
 
   integer                        ::  nn, nlev
   real(4), dimension(nlev_in+1)  ::  ilev 
   real(4), dimension(nlev_in+1)  ::  tlev 
   real(4), dimension(nlev_in+1)  ::  qlev 
   real(4), dimension(nlev_in+1)  ::  plev 
   real(4), dimension(nlev_in+1)  ::  allmissing
   logical, dimension(nlev_in+1)  ::  notmissing 
   real(4)                        ::  psfc
   real(4)                        ::  t2m
   real(4)                        ::  q2m
   real(4)                        ::  p50   , p100   , p150   , p300
   real(4)                        ::  temp50, temp100, temp150, temp300
   real(4)                        ::  qhum50, qhum100, qhum150, qhum300
   real(4)                        ::  tdew50, tdew150
   real(4)                        ::  p_old
   real(4)                        ::  tseg_old
   real(4)                        ::  tpar_old
   real(4)                        ::  tseg_mid
   real(4)                        ::  tpar_mid
   real(4)                        ::  tpar
   real(4)                        ::  p_segment
   real(4)                        ::  t_segment
   real(4)                        ::  q_segment
   real(4)                        ::  p_increment
   integer                        ::  ilo  , iup
   integer                        ::  lo50 , up50
   integer                        ::  lo100, up100
   integer                        ::  lo150, up150
   integer                        ::  lo300, up300
   real(4)                        ::  x_up, x_lo, y_up, y_lo
   real(4)                        ::  moist_lapse, dz
   real(4)                        ::  pmid, tmid, qmid
   real(4)                        ::  qseg_old
   real(4)                        ::  qsat
!-----------------------------------------------------------------------------


      !--------------------------------------------
      !--- Initialization and preliminary calculations
      !--------------------------------------------
      !--------------------------------
      !-- initialize output arrays
      !--------------------------------
      CTP    =  missing
      HILOW  =  missing
    
      !--------------------------------
      !-- initialize working arrays
      !--------------------------------
      nlev       =  nlev_in + 1
      plev(2:)   =  plev_in
      tlev(2:)   =  tlev_in
      qlev(2:)   =  qlev_in
      psfc       =  psfc_in
      t2m        =  t2m_in
      q2m        =  q2m_in
      plev(1)    =  psfc
      tlev(1)    =  t2m
      qlev(1)    =  q2m

      allmissing   =  missing
      notmissing   =  .false.

      !--------------------------------
      !-- define the level index
      !-- this is used for locating 
      !-- indices below and above 
      !-- desired pressure level
      !--------------------------------
      level_index: do nn=1,nlev
         ilev(nn)  =  real(nn)
      end do level_index

      !-----------------------------------------------------------------------------------
      !-- Make sure there are no missing critical inputs  --> if so return CTP HILOW missing
      !-----------------------------------------------------------------------------------
      if(     psfc.eq.missing   .or.  all(plev.eq.missing)  .or.  &
          all(tlev.eq.missing)  .or.  all(qlev.eq.missing)        )   return

#-----------------------------------------------------------------------------------
#-- Check pressure units --> need pressure to be in Pascals for plev and psfc
#-----------------------------------------------------------------------------------
gf.check_pressure_units(psfc)
gf.check_pressure_units(plev)


#-----------------------------------------------------------------------------
#-- Make sure there are no higher pressure levels than that given by surface pressure
#-- Important for models that have atmo levels below the surface
#-----------------------------------------------------------------------------
plev = plev.where( plev.le.psfc, np.nan )



      !-----------------------------------------------------------------------------------
      !-- Get 50, 100, 150, and 300mb levels above the ground
      !-----------------------------------------------------------------------------------
      p50   =   psfc - 5000.
      p100  =   psfc - 10000.
      p150  =   psfc - 15000.
      p300  =   psfc - 30000.


      !-----------------------------------------------------------------------------------
      !-- Find the index above and below each of the desired pressure levels (50,100,150,300)
      !-----------------------------------------------------------------------------------
      lo50   =   maxloc( ilev, DIM = 1, MASK = plev.ne.missing .and. plev-p50.ge.0 )
      up50   =   minloc( ilev, DIM = 1, MASK = plev.ne.missing .and. plev-p50.lt.0 )

      lo100  =   maxloc( ilev, DIM = 1, MASK = plev.ne.missing .and. plev-p100.ge.0 )
      up100  =   minloc( ilev, DIM = 1, MASK = plev.ne.missing .and. plev-p100.lt.0 )

      lo150  =   maxloc( ilev, DIM = 1, MASK = plev.ne.missing .and. plev-p150.ge.0 )
      up150  =   minloc( ilev, DIM = 1, MASK = plev.ne.missing .and. plev-p150.lt.0 )

      lo300  =   maxloc( ilev, DIM = 1, MASK = plev.ne.missing .and. plev-p300.ge.0 )
      up300  =   minloc( ilev, DIM = 1, MASK = plev.ne.missing .and. plev-p300.lt.0 )



      !--------------------------------------------------------------------------------
      !
      !              !!!  Low-Level Humidity (Hi_low) Section  !!!
      !
      !--------------------------------------------------------------------------------
      !--------------------------------------------
      !--- Perform linear interpolation to extract
      !--- each value at the desired level
      !--- This is done by finding the y-intercept
      !--- equal to the pressure level minus the
      !--- desired level (basically find the temp
      !--- and dew point corresponding to the level
      !--------------------------------------------
      x_up    =   plev(up50)-p50
      x_lo    =   plev(lo50)-p50
      y_up    =   tlev(up50)
      y_lo    =   tlev(lo50)
      temp50  =   y_up - ((y_up-y_lo)/(x_up-x_lo))*x_up

      y_up    =   qlev(up50)
      y_lo    =   qlev(lo50)
      qhum50  =   y_up - ((y_up-y_lo)/(x_up-x_lo))*x_up

      x_up    =   plev(up150)-p150
      x_lo    =   plev(lo150)-p150
      y_up    =   tlev(up150)
      y_lo    =   tlev(lo150)
      temp150 =   y_up - ((y_up-y_lo)/(x_up-x_lo))*x_up

      y_up    =   qlev(up150)
      y_lo    =   qlev(lo150)
      qhum150 =   y_up - ((y_up-y_lo)/(x_up-x_lo))*x_up


      !-----------------------------------------------------------------------------------
      !-- Calculate dew point temperature (K)
      !-----------------------------------------------------------------------------------
      call dew_point( qhum50 , p50 , tdew50 , missing )
      call dew_point( qhum150, p150, tdew150, missing )


      !-----------------------------------------------------------------------------------
      !-- Return Hi_low variable
      !-----------------------------------------------------------------------------------
      HILOW   =   (temp50 - tdew50)  +  (temp150 - tdew150)


      !--------------------------------------------------------------------------------
      !
      !            !!!  Convective Triggering Potential (CTP) Section  !!!
      !
      !--------------------------------------------------------------------------------
      !--------------------------------------------
      !--- Perform linear interpolation to extract
      !--- each value at the desired level
      !--- This is done by finding the y-intercept
      !--- equal to the pressure level minus the
      !--- desired level (basically find the temp
      !--- and dew point corresponding to the level
      !--------------------------------------------
      x_up    =   plev(up100)-p100
      x_lo    =   plev(lo100)-p100
      y_up    =   tlev(up100)
      y_lo    =   tlev(lo100)
      temp100 =   y_up - ((y_up-y_lo)/(x_up-x_lo))*x_up

      y_up    =   qlev(up100)
      y_lo    =   qlev(lo100)
      qhum100 =   y_up - ((y_up-y_lo)/(x_up-x_lo))*x_up

      x_up    =   plev(up300)-p300
      x_lo    =   plev(lo300)-p300
      y_up    =   tlev(up300)
      y_lo    =   tlev(lo300)
      temp300 =   y_up - ((y_up-y_lo)/(x_up-x_lo))*x_up

      y_up    =   qlev(up300)
      y_lo    =   qlev(lo300)
      qhum300 =   y_up - ((y_up-y_lo)/(x_up-x_lo))*x_up



      !----------------------------------------------------
      !--- Chop up the integration into 20 segments from
      !--- 100mb to 300mb above the ground
      !----------------------------------------------------
      p_increment  =  (p100 - p300) / real(nsegments)
      p_old        =  p100
      tseg_old     =  temp100
      tpar_old     =  temp100
      qseg_old     =  qhum100
      CTP          =  0.
      ctp_depth: do nn=1,nsegments

         !----------------------------------------------------
         !-- Pressure increment inbetween defined levels (Pa)
         !----------------------------------------------------
         p_segment  =  p100 - (p_increment * real(nn))
         iup        =  minloc( ilev, DIM = 1, MASK = plev.ne.missing .and. plev.lt.p_segment )
         ilo        =  maxloc( ilev, DIM = 1, MASK = plev.ne.missing .and. plev.gt.p_segment )

         !----------------------------------------------------
         !--- Perform another linear interpolation to get 
         !--- temperature and spc. humidity at the increment
         !----------------------------------------------------
         x_up    =   plev(iup) - p_segment
         x_lo    =   plev(ilo) - p_segment
         y_up    =   tlev(iup)
         y_lo    =   tlev(ilo)
         t_segment =   y_up - ((y_up-y_lo)/(x_up-x_lo))*x_up

         y_up    =   qlev(iup)
         y_lo    =   qlev(ilo)
         q_segment =   y_up - ((y_up-y_lo)/(x_up-x_lo))*x_up


         !----------------------------------------------------
         !--- Get moist adiabatic lapse rate [K/m] and 
         !--- the depth of the layer from lower to upper level
         !----------------------------------------------------
         pmid  =  ( (p_segment*log(p_segment)  +  p_old   *log(p_old))  /  log(p_segment*p_old) )
         tmid  =  ( (t_segment*log(p_segment)  +  tseg_old*log(p_old))  /  log(p_segment*p_old) )
         qmid  =  ( (q_segment*log(p_segment)  +  qseg_old*log(p_old))  /  log(p_segment*p_old) )
         dz    =  (p_old-p_segment) / (grav * pmid /(Rd*tmid*((1.+(qmid/ep)) / (1. + qmid))))
         qsat  =  saturation_specific_humidity( pmid, tmid, missing)
 
         moist_lapse  =  (grav/cp) * ( (1. + (Lv    * qsat)/(   Rd*tmid   )) /   &
                                       (1. + (Lv**2 * qsat)/(cp*Rv*tmid**2)) )


         !----------------------------------------------------
         !--- Get the parcel temperature
         !----------------------------------------------------
         tpar     =  tpar_old - moist_lapse*dz

         !----------------------------------------------------
         !--- Get mid-point temps from environment and parcel
         !----------------------------------------------------
         tpar_mid =  0.5 * (tpar      + tpar_old)
         tseg_mid =  0.5 * (t_segment + tseg_old)

         !----------------------------------------------------
         !--- Integrate from old increment to increment level
         !----------------------------------------------------
         CTP  =  CTP  +  (Rd * (tpar_mid-tseg_mid) * log(p_old/p_segment))

         !----------------------------------------------------
         !--- Update the last increment values 
         !----------------------------------------------------
         tpar_old  =  tpar
         tseg_old  =  t_segment
         qseg_old  =  q_segment
         p_old     =  p_segment

      end do ctp_depth


      return

 end subroutine ctp_hi_low
!---------------------------------------------------------------------------------










!---------------------------------------------------------------------------------
!
! subroutines:  Calculates the dew point temperature following the Arden Buck Eq.
!
!---------------------------------------------------------------------------------
 subroutine dew_point( qlev, plev, tdew, missing )

   implicit none

!
! Input/Output Variables
!
   real(4), intent(in ) :: missing  !** missing value
   real(4), intent(in ) :: qlev     !** specific humidity [kg/kg]
   real(4), intent(in ) :: plev     !** pressure [Pa]
   real(4), intent(out) :: tdew     !** dew point temp. [K]
!
! Local variables
!
   real(4)              ::  e
   real(4), parameter   ::  A=610.8, B=237.3, C=17.2693882

!---------------------------------------------------------------------------------

      !--------------------------------------------
      !--- Initialization and preliminary calculations
      !--------------------------------------------
      tdew  =  missing

      !--------------------------------------------
      !--- Vapor pressure and convert to Pa
      !--------------------------------------------
      e     =  (qlev*(plev/1e2))/(0.622+0.378*qlev)
      e     =  e*1e2

      !--------------------------------------------
      !--- Vapor pressure and convert to Pa
      !--------------------------------------------
      tdew  =  ( (log(e/A)*B) / (C-log(e/A)) ) + 273.15

 end subroutine dew_point
!---------------------------------------------------------------------------------









!---------------------------------------------------------------------------------
!
! subroutines:  Calculates the saturation specific humidity [kg/kg]
!               following the AMS glossary definition
!
!---------------------------------------------------------------------------------
 real(4) function saturation_specific_humidity(  p, t, missing )

   implicit none

!
! Input/Output Variables
!
   real(4), intent(in ) :: missing  !** missing value
   real(4), intent(in ) :: p        !** pressure [Pa]
   real(4), intent(in ) :: t        !** dry bulb temp. [K]
!
! Local variables
!
   real(4)              ::  press, esat
   real(4)              ::  numerator, denomenator
   real(4), parameter   ::  t0    = 273.15, ep=0.622, es0=6.11, a=17.269, b=35.86
   real(4), parameter   ::  onemep= 1.0 - ep

!---------------------------------------------------------------------------------

      !--------------------------------------------
      !--- Initialization and preliminary calculations
      !--------------------------------------------
      saturation_specific_humidity  =  missing

      !--------------------------------------------
      !--- Perform a quick check for missing values
      !--------------------------------------------
      if( t.eq.missing .or. p.eq.missing ) return

      !--------------------------------------------
      !--- Convert pressure to hPa
      !--------------------------------------------
      press       =  p/1e2

      !--------------------------------------------
      !--- Split out the numerator and denomenator 
      !--------------------------------------------
      numerator   =  ep* (es0*exp((a*( t-t0))/( t-b))) 
      denomenator =  press-onemep*(es0*exp((a*( t-t0))/( t-b)))

      !--------------------------------------------
      !--- Intermediate calculation
      !--------------------------------------------
      esat  =  numerator/denomenator

      !--------------------------------------------
      !--- Vapor pressure and co
      !--------------------------------------------
      saturation_specific_humidity  =  esat / (1 + esat)

 end function saturation_specific_humidity
!---------------------------------------------------------------------------------




end module Conv_Trig_Pot_Mod
