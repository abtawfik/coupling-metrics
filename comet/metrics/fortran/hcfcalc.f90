!---------------------------------------------------------------------------------
!
! Description:
!
! Provides quantities for 1) the potential temperature required to initiate moist convection
! as defined by the buoyant condensation level, 2) the height and pressure the PBL needs to
! reach to trigger moist convection, and 3) the potential temperature deficit between the
! current 2-m potential temperature and potential temp threshold.
! This method can be used as a condition for convective triggering (i.e. trigger when TDEF=0)
! Note the current formulation does not distinguish between shallow and deep convection
! but simply points to convective initiation.  This technique can also be used to identify
! locally triggered convection versus transient events using a single profile, where 
! locally triggered == where TDEF=0 and otherwise observed cloud cover is non-local
!
! Variables returned can be categorized into two groups:
! 1) THRESHOLD VARIABLES -- The first set of variales are those associated with 
!                           triggering convection and refer to TBM, BCLH, BCLP, and TDEF.
! 2) EVALUATION VARIABLE -- These variables return detailed information about the 
!                           current convective regime.  Namely, the amount of sensible
!                           and latent heat necessary for initiation (SHDEF_M and LHDEF_M),
!                           the most energy efficient pathway among the two (EADV_M), and
!                           the height, pressure, and potential temperature at which the
!                           transition from one energy advantange to another occurs
!                           (TRAN_H, TRAN_P, TRAN_T). These energy transition variables
!                           may not be always present because some vertical profiles
!                           may lie entirely within one regime. These variables will be
!                           returned as missing if TDEF=0 because there is no longer a
!                           physical meaning behind "energy advantage" and "transition".
!
! Note: Programatically evaluation variables are not required so the user can decided whether
!       they wish to include their calculation.  This is because Evaluation variables add
!       a noticable (but not excessive amount computation time.  Also, calculation of 
!       Evaluation variables is all or nothing, meaning you have to include all in the subroutine
!       call.  If you only include one Evaluation when calling the subroutine then the output will 
!       be returned as missing.
!
!
! Reference: --Methods--
!            Tawfik and Dirmeyer 2014 GRL: A processed based framework for 
!            quantifying the atmospheric background state of surface triggered
!            convection
!
!            Tawfik, A., P.A. Dirmeyer, J.A. Santanello Jr. (2015), The Heated
!            Condensation Framework. Part I: Description and Southern Great Plains 
!            Case Study, Journal of Hydromet. doi:10.1175/JHM-D-14-0117.1
!  
!            --Climatological Evaluation--
!            Tawfik, A., P.A. Dirmeyer, J.A. Santanello Jr. (2015), The Heated 
!            Condensation Framework. Part II:  Climatological behavior of convective
!            initiation and land-atmosphere coupling over the Continental United States,
!            Journal of Hydromet. doi:10.1175/JHM-D-14-0118.1
!
! Author and Revision History: 
! A.B. Tawfik on Aug 2014
!
!---------------------------------------------------------------------------------
module HCF_vars_calc

     !
     ! subroutines names (base name) 
     !
     public hcfcalc          ! calculates the full suite of HCF variables

!---------------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------------



!---------------------------------------------------------------------------------
!
! subroutines:  calculates buoyant condensation level and basic variables (THETA_BM; TDEF)
!
!---------------------------------------------------------------------------------
subroutine hcfcalc ( nlev1  , missing, tmp_in, press_in, qhum_in, hgt_in,  &
                     t2m    , psfc   , q2m   , h2m     ,                   &
                     TBM    , BCLH   , BCLP  , TDEF                        )

   implicit none
!
! Input/Output Variables
!
   integer, intent(in   )                       ::  nlev1       ! *** # of atmospheric levels
   real(4), intent(in   )                       ::  missing     ! *** Missing values
   real(4), intent(in   ), dimension(nlev1)     ::  tmp_in      ! *** Temperature (level), [K]
   real(4), intent(in   ), dimension(nlev1)     ::  hgt_in      ! *** Geometric Height above ground (level) [m]
   real(4), intent(in   ), dimension(nlev1)     ::  qhum_in     ! *** Specific Humidity (level) [kg/kg]
   real(4), intent(in   ), dimension(nlev1)     ::  press_in    ! *** Pressure (level) [Pa]
   real(4), intent(in   )                       ::  t2m         ! *** lowest level temperature  (typically 2-m) [K]
   real(4), intent(in   )                       ::  h2m         ! *** lowest level height       (typically 2-m) [m]
   real(4), intent(in   )                       ::  q2m         ! *** lowest level sp. humidity (typically 2-m) [kg/kg]
   real(4), intent(in   )                       ::  psfc        ! *** lowest level pressure     (typically sfc) [Pa]
   real(4), intent(out  )                       ::  TBM         ! *** buoyant mixing pot. temp (convective threshold) [K]
   real(4), intent(out  )                       ::  BCLH        ! *** height above ground of convective threshold [m]
   real(4), intent(out  )                       ::  BCLP        ! *** pressure of convective threshold level [Pa]
   real(4), intent(out  )                       ::  TDEF        ! *** potential temperature deficit need to initiate [K]
!
! Local variables
!
   real(4), parameter   ::  p_ref = 1e5 , Lv=2.5e6 , cp=1005.7, R_cp=287.04/1005.7
   real(4), parameter   ::  grav  = 9.81, Rd=287.04
   real(4), parameter   ::  pi    = 4.*atan(1.),  cp_g=cp/grav,  Lv_g=Lv/grav
   real(4), parameter   ::  r2d   = 180./pi
   real(4), parameter   ::  by100 = 1e2
   real(4), parameter   ::  t0    = 273.15, ep=0.622, es0=6.11, a=17.269, b=35.86
   real(4), parameter   ::  onemep= 1.0 - ep

   integer                       ::  zz, cc
   integer                       ::  nlev
   logical, dimension(nlev1+1)   ::  notmissing
   real(4), dimension(nlev1+1)   ::  rhoh
   real(4), dimension(nlev1+1)   ::  allmissing
   real(4), dimension(nlev1+1)   ::  pbar, qdef, qmix, qsat, dpress, qbar, logp, hbar, tbar
   real(4), dimension(nlev1+1)   ::  tmp_k, press, pot_k
   real(4), dimension(nlev1+1)   ::  hgt, qhum

   real(4)                       ::  p_up, t_up, h_up, q_up, p_lo, t_lo, h_lo, q_lo, m_up, m_lo
   real(4)                       ::  pot2m, qbcl
   integer                       ::  i_unsat, i_sat, num_unsat, num_sat



!-----------------------------------------------------------------------------

      !-----------------------------------------------------------------------------
      !-- Initialize output variables
      !-----------------------------------------------------------------------------
      TBM      =   missing
      BCLH     =   missing 
      BCLP     =   missing
      TDEF     =   missing


      !-----------------------------------------------------------------------------
      !-- Check input arrays to make sure all are NOT missing
      !-----------------------------------------------------------------------------
      if( all(tmp_in  .eq.missing)  .or.  all(qhum_in.eq.missing)   .or.   &
          all(press_in.eq.missing)  .or.  all(hgt_in .eq.missing)   .or.   &
                   t2m.eq.missing   .or.          q2m.eq.missing    .or.   &
                  psfc.eq.missing   .or.          h2m.eq.missing           )  then
          return
      end if 

      !-----------------------------------------------------------------------------
      !-- Store temporary working arrays and initialize
      !-----------------------------------------------------------------------------
      nlev       =  nlev1 + 1
      tmp_k(2:)  =  tmp_in
      hgt  (2:)  =  hgt_in
      qhum (2:)  =  qhum_in
      press(2:)  =  press_in

      tmp_k(1)   =  t2m
      hgt  (1)   =  h2m
      qhum (1)   =  q2m
      press(1)   =  psfc


      !-----------------------------------------------------------------------------
      !-- Initialize middle level variables
      !-----------------------------------------------------------------------------
      qdef         =  missing
      allmissing   =  missing
      notmissing   =  .false.

      tbar    =  missing
      qbar    =  missing
      dpress  =  missing
      pbar    =  missing
      hbar    =  missing
      logp    =  missing
      qmix    =  missing
      qsat    =  missing


      !-----------------------------------------------------------------------------
      !-- Make sure there are no higher pressure levels than that given by surface pressure
      !-- Important for models that have atmo levels below the surface
      !-----------------------------------------------------------------------------
      where( press(2:).ge.psfc ) press(2:) = missing


      !-----------------------------------------------------------------------------
      !-- Collapse the missing values along the level dimension so that mid-points
      !-- can be calculated using levels on either side of the missing level.
      !-- This would avoid just ignoring the level
      !-- Use the PACK function
      !-----------------------------------------------------------------------------
      notmissing =  .not.( press.eq.missing .or. qhum.eq.missing .or. tmp_k.eq.missing .or. hgt.eq.missing )
      tmp_k      =  pack (tmp_k, notmissing, allmissing)
      qhum       =  pack (qhum , notmissing, allmissing)
      press      =  pack (press, notmissing, allmissing)
      hgt        =  pack (hgt  , notmissing, allmissing)


      !-----------------------------------------------------------------------------
      !--------
      !--------  Calculate column potential temperature (K)
      !--------
      !-----------------------------------------------------------------------------
      where( tmp_k.ne.missing .and. press.ne.missing )   pot_k  =  tmp_k * ((p_ref/press))**(R_cp)


      !-----------------------------------------------------------------------------
      !--------
      !--------  Ignore missing data levels when calculating midpoints
      !--------  Important when using observational data that can have missing levels
      !--------
      !-----------------------------------------------------------------------------
      hbar  =  hgt
      pbar  =  press
      tbar  =  tmp_k
      !-----------------------------------------------------------------------------
      !--------
      !--------  Calculate middle layer specific humidity average [kg/kg]
      !--------  1st layer = the 2m specific humidity above then layer averages above
      !--------
      !-----------------------------------------------------------------------------
      qbar  =  qhum
      where( qhum(1:nlev-1).ne.missing  .and.  press(1:nlev-1).ne.missing  .and. &
             qhum(2:nlev  ).ne.missing  .and.  press(2:nlev  ).ne.missing        )
             qbar(2:nlev)    =  ((qhum(2:nlev  )*log(press(2:nlev  ))  + &
                                  qhum(1:nlev-1)*log(press(1:nlev-1))) / &
                                  log(press(2:nlev)* press(1:nlev-1)))
      endwhere


      !-----------------------------------------------------------------------------
      !--------
      !--------  Calculate pressure difference of each layer
      !--------
      !-----------------------------------------------------------------------------
      if( dpress(1).le.0 ) then
          dpress(1)  =  1.    !*** set to 1 Pa because the h2m is likley zero
      else
          dpress(1)  =  (psfc / (Rd * t2m * ((1. + (q2m/ep)) / (1. + q2m)) )) * grav * h2m
      end if
      where( pbar(1:nlev-1).ne.missing  .and.  pbar(2:nlev).ne.missing )
             dpress(2:nlev)  =  press(1:nlev-1) - press(2:nlev)
      endwhere


      !-----------------------------------------------------------------------------
      !--------
      !--------  Calculate log pressure to linearize it for slope calculation
      !--------
      !-----------------------------------------------------------------------------
      where( pbar.ne.missing ) logp = log(pbar)


      !-----------------------------------------------------------------------------
      !--------
      !--------  Calculate mixed layer specific humidity  and column density [kg/kg]
      !--------
      !-----------------------------------------------------------------------------
      where( dpress.ne.missing .and. qbar.ne.missing )
             qmix  =  qbar  * dpress/grav
             rhoh  =  dpress/grav
      endwhere
      do zz = 2,nlev
         if( qmix  (zz).ne.missing .and. qmix  (zz-1).ne.missing ) qmix  (zz)  =  qmix  (zz-1) + qmix  (zz)
         if( rhoh  (zz).ne.missing .and. rhoh  (zz-1).ne.missing ) rhoh  (zz)  =  rhoh  (zz-1) + rhoh  (zz)
      end do


      !-----------------------------------------------------------------------------
      !--------
      !--------  Calculate saturation specific humidity at each level [kg/kg]
      !--------
      !-----------------------------------------------------------------------------
      pbar  =  pbar/1e2
      where( tbar.ne.missing .and. pbar.ne.missing )
         qsat  =  by100*0.01 *(ep* (es0*exp((a*( tbar-t0))/( tbar-b))) ) /  &
                      (pbar-onemep*(es0*exp((a*( tbar-t0))/( tbar-b))))
         qsat  =  qsat/(1.+qsat)
      endwhere
      pbar  =  pbar*1e2


      !-----------------------------------------------------------------------------
      !--------
      !--------  Calculate specific humidity deficit [kg/kg]
      !--------
      !-----------------------------------------------------------------------------
      where( qmix.ne.missing .and. rhoh.ne.missing )
         qmix  =  qmix / rhoh
         qdef  =  qsat - qmix
      endwhere

      !-----------------------------------------------------------------------------
      !--- Check to make sure qdef is always negative when outside of the tropo
      !--- Assume a tropopause height of 10 km; so BCL cannot be higher
      !-----------------------------------------------------------------------------
      where( hbar.ge.10000. )
         qdef  =  -1.
      endwhere


      !***********************************************************
      !***   Calculate slope of each variable to find the      ***
      !***   y-intercept for each variable;                    ***
      !***   Meaning locate the two data points surrounding    ***
      !***   the sign change in qdef and linearly interpolate  ***
      !***   to find the "zero point"                          ***
      !***********************************************************
      !----- Find the point where the sign first turns negative from the ground up
      !*** Highest unsaturated level ***
      num_unsat    =   count(qdef.ne.missing .and. qdef.gt.0, DIM = 1)
      if( num_unsat.gt.0 ) then
         i_unsat   =   maxloc( hbar, DIM = 1, MASK = qdef.ne.missing .and. qdef.gt.0 )
      else
         i_unsat   =   0
      end if
      !*** Lowest saturated level ***
      num_sat      =   count(qdef.ne.missing .and. qdef.le.0, DIM = 1)
      if( num_sat.gt.0 ) then
         i_sat     =   minloc( hbar, DIM = 1, MASK = qdef.ne.missing .and. qdef.le.0 )
      else
         i_sat     =   0
      end if

      !-----------------------------------------------------------------------------
      !--- If all levels are saturated then put the deficit to zero
      !-----------------------------------------------------------------------------
      if( num_unsat.eq.0 ) then
         pot2m  =    (t2m) * ((p_ref/psfc)**(R_cp))  *  (1. + 0.61*qmix(1))
         BCLP   =    psfc
         BCLH   =    h2m
         TBM    =    pot2m
         TDEF   =    0.
         return
      end if


      !-----------------------------------------------------------------------------
      !--- Check to see if first level is satured; (Foggy scenario)
      !--- If so then check the second and third layers to see if fog will dissipate
      !--- If the 2nd and/or 3rd are not saturated then recalc CONVECTIVE saturation
      !--- transition level
      !-----------------------------------------------------------------------------
      if( i_sat.gt.1 .and. i_unsat.gt.i_sat )   i_unsat  =  i_sat - 1
      if( i_sat.eq.1 ) then
          cc = 0
          do zz=2,nlev-1
             !**** make sure initiation level is below 100 hPa above the ground to ensure it is actually fog
             if( (psfc-pbar(zz))/1e2.gt.100 ) exit

             !**** If within the 100 hPa layer above the ground then try to erode the fog layer first
             !**** to determine the convective initiation layer
             i_sat  =  minloc( hbar(zz:), DIM = 1, MASK = qdef(zz:).ne.missing .and. qdef(zz:).le.0 )
             cc     =  cc + 1

             !**** If still saturated then cycle
             if( i_sat.eq.1 ) cycle
             i_sat  =  i_sat + cc
             exit
          end do
          i_unsat   =  i_sat - 1
      end if


      !-----------------------------------------------------------------------------
      !--- If all layers below 100 hPa above the ground are still saturated then call it all saturated
      !--- And use the 1st level stats and call it "convective" because fog is unlikely to be
      !--- deeper than 100 hPa above the ground
      !-----------------------------------------------------------------------------
      if( i_unsat.eq.0 ) then
         pot2m  =    (t2m) * ((p_ref/psfc)**(R_cp))  *  (1. + 0.61*qmix(1))
         BCLP   =    psfc
         BCLH   =    h2m
         TBM    =    pot2m
         TDEF   =    0.
         !do zz=1,nlev
         !   write(*,*)  zz,  pbar(zz)/1e2,  qmix(zz)*1e3,  qsat(zz)*1e3,  qdef(zz)*1e3
         !end do
         return
      end if


      !-----------------------------------------------------------------------------
      !--- Check to make sure these are adjacent layers
      !--- If not then there is a problem with this implementation
      !-----------------------------------------------------------------------------
      if( i_unsat.eq.0  .or.  i_sat.eq.0 ) then
          write(*,*) " =========== ERROR  in locating saturation profiles ============ "
          write(*,*) i_sat, i_unsat
          do zz=1,nlev
             write(*,*)  zz,  pbar(zz)/1e2,  qmix(zz)*1e3,  qsat(zz)*1e3,  qdef(zz)*1e3
          end do
          write(*,*)
          write(*,*)
          write(*,*)  psfc, h2m, t2m, q2m
          do zz=1,nlev
             write(*,*)  press_in(zz), hgt_in(zz), tmp_in(zz), qhum_in(zz)
          end do
          write(*,*)
          write(*,*)
          stop
      end if



      !-----------------------------------------------------------------------------
      !--- Get the upper and lower bounds for each variable to be calc'd at the BCL
      !-----------------------------------------------------------------------------
      p_up        =  logp(i_sat)
      t_up        =  tbar(i_sat)
      h_up        =  hbar(i_sat)
      q_up        =  qdef(i_sat)
      m_up        =  qmix(i_sat)

      p_lo        =  logp(i_unsat)
      t_lo        =  tbar(i_unsat)
      h_lo        =  hbar(i_unsat)
      q_lo        =  qdef(i_unsat)
      m_lo        =  qmix(i_unsat)

      !-----------------------------------------------------------------------------
      !--- Calculate output variables; BCL height, BCL pressure,
      !--- Buoyant Mixing Potential Temp, and Potential Temperature Deficit
      !-----------------------------------------------------------------------------
      BCLP  =  exp( p_up - ((p_up-p_lo)/(q_up-q_lo))*q_up )
      BCLH  =     ( h_up - ((h_up-h_lo)/(q_up-q_lo))*q_up )
      qbcl  =     ( m_up - ((m_up-m_lo)/(q_up-q_lo))*q_up )
      TBM   =     ( t_up - ((t_up-t_lo)/(q_up-q_lo))*q_up ) * ((p_ref/BCLP)**(R_cp))

      !-----------------------------------------------------------------------------
      !--- Virtual Potential Temperature (K) calculation using the mixed humidty,
      !--- *** THIS is an assumption!!  only influences TDEF but an important
      !--- effect because if pot2m is close to TBM then a slight change in qbcl
      !--- can mean the difference between initiation (TDEF=0) or not
      !--- This should only be an issue over very shallow pbls ...
      !-----------------------------------------------------------------------------
      pot2m  =  (t2m) * ((p_ref/psfc)**(R_cp))  *  (1. + 0.61*qbcl)
      TDEF   =  TBM  - pot2m
      if(TDEF.lt.0) TDEF=0.

      return


end subroutine hcfcalc





end module HCF_vars_calc
