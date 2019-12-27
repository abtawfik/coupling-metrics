!---------------------------------------------------------------------------------
! Purpose:
!
! !DESCRIPTION:  
! Separates the surface flux contribution (evaporative fraction) from
! the non-evaporative contribution to changes in relative humidity at the top of the 
! boundary layer.  Terms quantifying dry air entrainment, boundary layer growth, and 
! boundary layer heating contribution to relative humidity changes at the top of the
! boundary layer are also returned.
!  
! !NOTES:
! The module works ONLY if you have the heat entrainment ratio.  This can be estimated
! by using the Mixing Diagram incremental method (See Santanello et al. 2009 
! Section 2 b.2) 
! 
! References: Ek and Mahrt 1994,  Daytime Evolution of Relative humidity at the
!             Boundary Layer Top
!
!             Ek and Holtslag 2004,  Influence of Soil Moisture on Boundary Layer Cloud
!             Development 
!
! Author and Revision History: 
! Original Code -- A.B. Tawfik on Apr 2015
!
!---------------------------------------------------------------------------------
module RH_Tend_Mod

     !
     ! subroutine name 
     !
     public rh_tend_calc      ! the main RH-tendency interface

!---------------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
!
! subroutines:  main relative humidity tendency calculation routine
!
!---------------------------------------------------------------------------------
  subroutine rh_tend_calc ( nlev     ,  ntim                       ,         &  
                            tmp_in   ,  qhum_in,  hgt_in,  press_in,         &
                            pbl_h    ,  shf    ,  lhf   ,  dt      ,         &
                            ef       ,  ne     ,                             &
                            pbl_heating ,  pbl_growth, dry_entrain,          &
                            dRH_estimate,  missing                           )
   implicit none
!
! Input/Output Variables
!
   integer, intent(in   )                              ::  ntim, nlev
   real(4), intent(in   )                              ::  missing        ! *** missing value - useful for obs
   real(4), intent(in   )                              ::  dt             ! *** time increment of the hours dimension [seconds]
   real(4), intent(in   ), dimension(ntim,nlev)        ::  press_in       ! *** pressure - (time,nlev) [Pa]
   real(4), intent(in   ), dimension(ntim,nlev)        ::  hgt_in         ! *** height   - (time,nlev) [m]
   real(4), intent(in   ), dimension(ntim,nlev)        ::  tmp_in         ! *** temperature - (time,nlev) [K]
   real(4), intent(in   ), dimension(ntim,nlev)        ::  qhum_in        ! *** specific humidity - (time,nlev) [kg/kg]
   real(4), intent(in   ), dimension(ntim)             ::  pbl_h          ! *** boundary layer height - (time) [m]
   real(4), intent(in   ), dimension(ntim)             ::  shf, lhf       ! *** surface energy fluxes - (time) [W/m2]
   real(4), intent(out  ), dimension(ntim)             ::  ef             ! *** evaporative fraction - (time) [unitless]
   real(4), intent(out  ), dimension(ntim)             ::  ne             ! *** non-evaporative term - (time) [unitless]
   real(4), intent(out  ), dimension(ntim)             ::  pbl_growth     ! *** PBL growth contribution to d(RH)/dt [%/s]
   real(4), intent(out  ), dimension(ntim)             ::  pbl_heating    ! *** PBL heating contribution to d(RH)/dt [%/s]
   real(4), intent(out  ), dimension(ntim)             ::  dry_entrain    ! *** Dry air entrainment to d(RH)/dt [%/s]
   real(4), intent(out  ), dimension(ntim)             ::  dRH_estimate   ! *** d(RH)/dt estimated by RH-tend method [%/s]
!
! Local variables
!
   real(4), parameter   ::  p_ref  = 1e5  , Lv    = 2.5e6 , cp   = 1005.7, Rd  = 287.04 
   real(4), parameter   ::  grav   = 9.81 , Lv_cp = Lv/cp , R_cp = Rd/cp
   real(4), parameter   ::  Rv     = 461.5
   real(4), parameter   ::  pi     = 4.*atan(1.),  cp_g=cp/grav,  Lv_g=Lv/grav
   real(4), parameter   ::  r2d    = 180./pi
   real(4), parameter   ::  by100  = 1e2
   real(4), parameter   ::  t0     = 273.15, ep=0.622, es0=6.11, a=17.269, b=35.86
   real(4), parameter   ::  onemep = 1 - ep

   integer                       ::  tt, zz 
   real(4), dimension(ntim,nlev) ::  psfc_2d      
   real(4), dimension(ntim,nlev) ::  pblh_2d      
   real(4), dimension(ntim,nlev) ::  pot_k      
   real(4), dimension(ntim,nlev) ::  tmp_k
   real(4), dimension(ntim,nlev) ::  press
   real(4), dimension(ntim,nlev) ::  qhum
   real(4), dimension(ntim,nlev) ::  hgt   
   real(4), dimension(ntim,nlev) ::  qsat 
   real(4), dimension(ntim,nlev) ::  logp  

   integer, dimension(ntim)   ::  iabove_pbl, ibelow_pbl
   real(4), dimension(ntim)   ::  psfc
   real(4), dimension(ntim)   ::  C_theta      
   real(4), dimension(ntim)   ::  delta_q    
   real(4), dimension(ntim)   ::  g_theta    
   real(4), dimension(ntim)   ::  pbl_theta    
   real(4), dimension(ntim)   ::  pbl_t      
   real(4), dimension(ntim)   ::  pbl_q        
   real(4), dimension(ntim)   ::  pbl_qs      
   real(4), dimension(ntim)   ::  pbl_p        
   real(4), dimension(ntim)   ::  pblh_diff
   real(4), dimension(ntim)   ::  relh        
   real(4), dimension(ntim)   ::  c_1         
   real(4), dimension(ntim)   ::  c_2        
   real(4), dimension(ntim)   ::  pbl_h0        
   real(4), dimension(ntim)   ::  pblh_mid
   real(4), dimension(ntim)   ::  lhf0       
   real(4), dimension(ntim)   ::  shf0        
   real(4), dimension(ntim)   ::  rho

   real(4), dimension(ntim)   ::  rho2m
   real(4), dimension(ntim)   ::  q2m
   real(4), dimension(ntim)   ::  t2m

   real(4), dimension(ntim)   ::  bowen_s
   real(4), dimension(ntim)   ::  cp_theta_final
   real(4), dimension(ntim)   ::  cp_theta_initial
   real(4), dimension(ntim)   ::  cp_theta
   real(4), dimension(ntim)   ::  cp_deltaT
   real(4), dimension(ntim)   ::  Lv_qhum_initial
   real(4), dimension(ntim)   ::  Lv_qhum_final
   real(4), dimension(ntim)   ::  Lv_qhum
   real(4), dimension(ntim)   ::  Lv_q_0
   real(4), dimension(ntim)   ::  cp_T_0
   real(4), dimension(ntim)   ::  shf_ent

   real(4), dimension(ntim,nlev) ::  pbl_diff
   real(4), dimension(ntim)      ::  p_up, t_up, q_up, s_up, h_up
   real(4), dimension(ntim)      ::  p_lo, t_lo, q_lo, s_lo, h_lo
   integer, dimension(ntim)      ::  num_above, num_below, i_below, i_above
  
   character(len=24), parameter :: FMT1 = "(9(F12.4,2x))"   
   character(len=23), parameter :: FMT2 = "(4(I4,2x), 6(F12.4,2x))"   

!-----------------------------------------------------------------------------

      !--------------------------------------------
      !--- Initialization and preliminary calculations
      !--------------------------------------------
      !--
      !-- initialize output arrays
      !--
      ef           =  missing
      ne           =  missing
      pbl_growth   =  missing
      pbl_heating  =  missing
      dry_entrain  =  missing
      dRH_estimate =  missing

      !-----------------------------------------------------------------------------
      !-- Check input arrays to make sure all are NOT missing
      !-----------------------------------------------------------------------------
      if( all(tmp_in  .eq.missing)  .or.  all(qhum_in.eq.missing)   .or.   &
          all(press_in.eq.missing)  .or.  all(hgt_in .eq.missing)   .or.   &
          all(   pbl_h.eq.missing)  .or.  all(    lhf.eq.missing)   .or.   &
          all(     shf.eq.missing)                                         )  then
          return
      end if


      !--
      !-- initialize working arrays
      !--
      rho          =  missing
      C_theta      =  missing
      delta_q      =  missing
      g_theta      =  missing
      pbl_t        =  missing
      pbl_q        =  missing
      pbl_qs       =  missing
      pbl_p        =  missing 
      pbl_theta    =  missing
      pblh_2d      =  missing
      relh         =  missing
      c_1          =  missing
      c_2          =  missing
      psfc_2d      =  missing
      psfc         =  missing
      pot_k        =  missing
      pbl_h0       =  missing
      pblh_mid     =  missing
      pbl_diff     =  missing
      pblh_diff    =  missing
      shf0         =  missing
      lhf0         =  missing
      ibelow_pbl   =  -1
      iabove_pbl   =  -1

      bowen_s      =  missing
      shf_ent      =  missing
      rho2m        =  missing
      q2m          =  missing
      t2m          =  missing
      Lv_qhum      =  missing
      cp_theta     =  missing
      cp_T_0       =  missing
      Lv_q_0       =  missing
      Lv_qhum_initial   =  missing
      Lv_qhum_final     =  missing
      cp_theta_initial  =  missing
      cp_theta_final    =  missing

      p_up         =  missing
      t_up         =  missing 
      s_up         =  missing 
      q_up         =  missing 
      h_up         =  missing 

      p_lo         =  missing 
      t_lo         =  missing 
      s_lo         =  missing 
      q_lo         =  missing 
      h_lo         =  missing 

      pbl_p        =  missing         
      pbl_q        =  missing
      pbl_qs       =  missing
      pbl_t        =  missing

      tmp_k        =  missing        
      press        =  missing
      qhum         =  missing
      hgt          =  missing
      tmp_k        =  tmp_in        
      press        =  press_in
      qhum         =  qhum_in
      hgt          =  hgt_in


      !------------------------------------------------------------------
      !--                                                              --
      !-- Set surface pressure and remove strange zeros as a check     --
      !--                                                              --
      !------------------------------------------------------------------
      where( press.eq.0 )  press = missing
      psfc   =   press(:,1)
      t2m    =   tmp_k(:,1)
      q2m    =   qhum (:,1)

      !----------------------------------------------------------------------------------
      !-- Approximate air density (kg/m3)
      !----------------------------------------------------------------------------------
      where( t2m.ne.missing .and. psfc.ne.missing .and. q2m.ne.missing )
         rho2m  =  psfc / (Rd * t2m * ((1. + (q2m/ep)) / (1. + q2m)))
      endwhere

      !------------------------------------------------------------------
      !--                                                              --
      !-- Map to 2D and remove pressures greater than surface          --
      !-- Useful for models on constant pressure fields but with a     --
      !-- surface pressure output field                                --
      !--                                                              --
      !------------------------------------------------------------------
      do zz=1,nlev
         psfc_2d(:,zz)  =  psfc
         pblh_2d(:,zz)  =  pbl_h
      end do
      where( press.gt.psfc_2d ) 
          tmp_k  =  missing
          press  =  missing
          hgt    =  missing
          qhum   =  missing
      endwhere



      !-----------------------------------------------------------------------------
      !-- Check input arrays to make sure all are NOT missing
      !-- NEEDS to be done again just in case observational data are incorrect
      !-- and there are pressure measurement errors (as can be seen in ARM BE)
      !-----------------------------------------------------------------------------
      if( all(tmp_k.eq.missing)  .or.  all(qhum.eq.missing)   .or.   &
          all(press.eq.missing)  .or.  all( hgt.eq.missing)   .or.   &
          all(pbl_h.eq.missing)  .or.  all( lhf.eq.missing)   .or.   &
          all(  shf.eq.missing)                                         )  then
          return
      end if



      !---------------------------------------------------------------
      !--                                                           --
      !-- Calculate temperature, humidity, and pressure at the PBL  --
      !--                                                           --
      !---------------------------------------------------------------
      !-- Calculate log of pressure
      where( press.ne.missing ) logp = log(press)

      !-- Calucalte the saturation specific humidity (kg/kg)
      where( press.ne.missing )  press  =  press/1e2
      where( tmp_k.ne.missing .and. press.ne.missing )
         qsat  =  by100*0.01 *(ep* (es0*exp((a*( tmp_k-t0))/( tmp_k-b))) ) /  &
                     (press-onemep*(es0*exp((a*( tmp_k-t0))/( tmp_k-b))))
         qsat  =  qsat/(1.+qsat)
      endwhere
      where( press.ne.missing )  press  =  press*1e2


      !-- Height difference between boundary layer and height
      where( pblh_2d.ne.missing  .and.  hgt.ne.missing )  pbl_diff  =  pblh_2d  -  hgt


      !----- Find the point where the sign first turns negative from the ground up
      !*** Lowest height above the PBL ***
      num_above  =   0
      num_above  =   count(pbl_diff.ne.missing .and. pbl_diff.le.0, DIM = 2)
      i_above    =   minloc( hgt, DIM = 2, MASK = pbl_diff.ne.missing .and. pbl_diff.le.0 )

      !*** Highest height below the PBL ***
      num_below   =   0
      num_below   =   count(pbl_diff.ne.missing .and. pbl_diff.gt.0, DIM = 2)
      i_below     =   maxloc( hgt, DIM = 2, MASK = pbl_diff.ne.missing .and. pbl_diff.gt.0 )

      !--- Get the upper and lower bounds for each variable to be calc'd at the BCL
      do tt=1,ntim
        if( i_above(tt).eq.0 .or. i_below(tt).eq.0 ) then
           cycle
        end if
        p_up(tt)    =  logp    (tt,i_above(tt))
        t_up(tt)    =  tmp_k   (tt,i_above(tt))
        s_up(tt)    =  qsat    (tt,i_above(tt))
        q_up(tt)    =  qhum    (tt,i_above(tt))
        h_up(tt)    =  pbl_diff(tt,i_above(tt))

        p_lo(tt)    =  logp    (tt,i_below(tt))
        t_lo(tt)    =  tmp_k   (tt,i_below(tt))
        s_lo(tt)    =  qsat    (tt,i_below(tt))
        q_lo(tt)    =  qhum    (tt,i_below(tt))
        h_lo(tt)    =  pbl_diff(tt,i_below(tt))
      end do
 
      !--- Calculate pressure, sp. humidity, sat. sp. humid, and temp at PBL
      where( p_up .ne.missing  .and.  h_up  .ne.missing )  pbl_p    =  exp( p_up - ((p_up-p_lo)/(h_up-h_lo))*h_up )
      where( q_up .ne.missing  .and.  h_up  .ne.missing )  pbl_q    =     ( q_up - ((q_up-q_lo)/(h_up-h_lo))*h_up )
      where( s_up .ne.missing  .and.  h_up  .ne.missing )  pbl_qs   =     ( s_up - ((s_up-s_lo)/(h_up-h_lo))*h_up )
      where( t_up .ne.missing  .and.  h_up  .ne.missing )  pbl_t    =     ( t_up - ((t_up-t_lo)/(h_up-h_lo))*h_up )
      where( pbl_q.ne.missing  .and.  pbl_qs.ne.missing )  relh     =     (pbl_q/pbl_qs)


      !--
      !-- Potential temperature for the profile and at the PBL [K]
      !-- Reference pressure is at the surface and NOT 1000 hPa
      !--
      where( tmp_k.ne.missing .and. press.ne.missing )   pot_k      =  tmp_k * ((psfc_2d/press))**(R_cp)
      where( pbl_t.ne.missing .and. pbl_p.ne.missing )   pbl_theta  =  pbl_t * ((psfc   /pbl_p))**(R_cp)





      !-------------------------------------------------------
      !--                                                   --
      !--                  BEGIN SECTION                    --
      !-- Estimate heat entrainment ratio for each timestep --
      !--     using the mixing diagram stepwise method      --
      !-- calc'd using mixing_diagram.f90 reduce routine    --
      !--       Based off of Santanello et al. 2009         --
      !--                                                   --
      !-------------------------------------------------------
      !******* SHIFT FLUXES to get C_theta using mixing space
      pbl_h0(:ntim-1)   =   pbl_h(2:)
      shf0  (:ntim-1)   =   shf  (2:)
      lhf0  (:ntim-1)   =   lhf  (2:)

      !----------------------------------------------------------------------------------
      !-- Calculate 2-m potential temperature  * specific heat capacity (J/kg)
      !----------------------------------------------------------------------------------
      !  where( t2m.ne.missing .and. psfc.ne.missing )  cp_theta  =  cp * (t2m * ((p_ref/psfc))**(R_cp))
      where( t2m.ne.missing .and. psfc.ne.missing )  cp_theta  =  cp * t2m

      !----------------------------------------------------------------------------------
      !-- Calculate 2-m specific humidity  * latent heat of vaporization (J/kg)
      !----------------------------------------------------------------------------------
      where( q2m.ne.missing                       )  Lv_qhum   =  Lv * q2m

      !----------------------------------------------------------------------------------
      !-- Calculate surface bowen ratio (unitless)
      !----------------------------------------------------------------------------------
      where( lhf0.ne.0 .and. shf0.ne.missing .and. lhf0.ne.missing )  bowen_s  =  shf0/lhf0


      !----------------------------------------------------------------------------------
      !-- Surface heat and moisture vector components (see Eq 1 in Santanello et al. 2009)
      !----------------------------------------------------------------------------------
      where( shf0.ne.missing .and. rho2m.ne.missing .and. pbl_h0.ne.missing .and. rho2m*pbl_h0.ne.0 )
         cp_deltaT  =  (shf0*dt) / (rho2m*pbl_h0)
      elsewhere
         cp_deltaT  =  missing
      endwhere

      !----------------------------------------------------------------------------------
      !-- Define the start and end times of the increment for sensible and latent heat
      !----------------------------------------------------------------------------------
      cp_theta_initial         =   cp_theta        !*** start of time increment
      cp_theta_final(:ntim-1)  =   cp_theta(2:)    !*** end of time increment
      Lv_qhum_initial          =   Lv_qhum         !*** start of time increment
      Lv_qhum_final (:ntim-1)  =   Lv_qhum (2:)    !*** end of time increment

      where( Lv_qhum_initial.ne.missing .and. cp_deltaT.ne.missing .and. bowen_s.ne.missing )
             Lv_q_0  =  Lv_qhum_initial   +  (cp_deltaT/bowen_s)
      end where
      where( cp_theta_initial.ne.missing .and. cp_deltaT.ne.missing )
             cp_T_0  =  cp_theta_initial  +  cp_deltaT
      end where

      !--------------------------------------------------
      !    Calculate Entrainment Heat Flux in Wm-2
      ! !!! THIS IS WITHOUT THE ADVECTION VECTOR !!!
      ! !!! Entrainment ratio could incorrectly  !!!
      ! !!!    estimated in windy conditions     !!!   
      !--------------------------------------------------
      where( cp_theta_final.ne.missing .and. rho2m.ne.missing .and. pbl_h0.ne.missing .and. cp_T_0.ne.missing )
          shf_ent  =  (((cp * ((cp_theta_final/cp)  -  &
                        (cp_T_0          /cp)))) * (rho2m * pbl_h0))/(dt)
      endwhere
      where( shf_ent(:ntim-1).ne.missing  .and.  shf0(:ntim-1).ne.missing  .and.  shf0(:ntim-1).ne.0 )
          C_theta(2:)  =  shf_ent(:ntim-1)/shf0(:ntim-1)
      endwhere
      !-------------------------------------------------------
      !--                                                   --
      !--            END OF Entrainment SECTION             --
      !--   used to calculate C_theta (entrainment ratio)   --
      !--                                                   --
      !-------------------------------------------------------





      !-------------------------------------------------------
      !--
      !--           Get Midpoint average in time
      !--
      !-------------------------------------------------------
      shf0           =   missing
      lhf0           =   missing
      where( pbl_t(2:).ne.missing .and. pbl_t(:ntim-1).ne.missing )
         pbl_t    (2:)  =   0.5 * (pbl_t    (1:ntim-1) + pbl_t    (2:))
      endwhere

      where( pot_k(2:,:).ne.missing .and. pot_k(:ntim-1,:).ne.missing )
         pot_k    (2:,:)  =   0.5 * (pot_k    (1:ntim-1,:) + pot_k    (2:,:))
      endwhere

      where( pblh_2d(2:,:).ne.missing .and. pblh_2d(:ntim-1,:).ne.missing )
         pblh_2d  (2:,:)  =   0.5 * (pblh_2d(1:ntim-1,:) + pblh_2d(2:,:))
      endwhere

      where( hgt(2:,:).ne.missing .and. hgt(:ntim-1,:).ne.missing )
         hgt     (2:,:)  =   0.5 * (hgt    (1:ntim-1,:) + hgt    (2:,:))
      endwhere

      where( qhum(2:,:).ne.missing .and. qhum(:ntim-1,:).ne.missing )
         qhum    (2:,:)  =   0.5 * (qhum    (1:ntim-1,:) + qhum    (2:,:))
      endwhere

      where( pbl_h(2:).ne.missing .and. pbl_h(:ntim-1).ne.missing )
         pblh_mid (2:)  =   0.5 * (pbl_h    (1:ntim-1) + pbl_h    (2:))
      endwhere

      where( pbl_p(2:).ne.missing .and. pbl_p(:ntim-1).ne.missing )
         pbl_p    (2:)  =   0.5 * (pbl_p    (1:ntim-1) + pbl_p    (2:))
      endwhere

      where( pbl_q(2:).ne.missing .and. pbl_q(:ntim-1).ne.missing )
         pbl_q    (2:)  =   0.5 * (pbl_q    (1:ntim-1) + pbl_q    (2:))
      endwhere

      where( pbl_qs(2:).ne.missing .and. pbl_qs(:ntim-1).ne.missing )
         pbl_qs   (2:)  =   0.5 * (pbl_qs   (1:ntim-1) + pbl_qs   (2:))
      endwhere

      where( pbl_theta(2:).ne.missing .and. pbl_theta(:ntim-1).ne.missing )
         pbl_theta(2:)  =   0.5 * (pbl_theta(1:ntim-1) + pbl_theta(2:))
      endwhere

      where( relh(2:).ne.missing .and. relh(:ntim-1).ne.missing )
         relh     (2:)  =   0.5 * (relh     (1:ntim-1) + relh     (2:))
      endwhere

      where( shf(2:).ne.missing .and. shf(:ntim-1).ne.missing )
         shf0     (2:)  =   0.5 * (shf      (1:ntim-1) + shf      (2:))
      endwhere

      where( lhf(2:).ne.missing .and. lhf(:ntim-1).ne.missing )
         lhf0     (2:)  =   0.5 * (lhf      (1:ntim-1) + lhf      (2:))
      endwhere

      where( psfc(2:).ne.missing .and. psfc(:ntim-1).ne.missing )
         psfc     (2:)  =   0.5 * (psfc      (1:ntim-1) + psfc      (2:))
      endwhere

      psfc     (1)   =   missing
      pbl_t    (1)   =   missing
      pbl_p    (1)   =   missing
      pbl_q    (1)   =   missing
      pbl_qs   (1)   =   missing
      pbl_theta(1)   =   missing
      relh     (1)   =   missing
      shf0     (1)   =   missing
      lhf0     (1)   =   missing
      C_theta  (1)   =   missing
      pot_k    (1,:) =   missing
      hgt      (1,:) =   missing
      qhum     (1,:) =   missing
      pblh_2d  (1,:) =   missing

      !--                                                
      !-- Approximate air density at the PBL (kg/m3)
      !--                                                  
      where( pbl_t.ne.missing .and. pbl_p.ne.missing .and. pbl_q.ne.missing )
         rho  =  pbl_p / (Rd * pbl_t * ((1. + (pbl_q/ep)) / (1. + pbl_q)))
      endwhere
      


      !----------------------------------
      !--                              --
      !-- Evaporative Fraction Section --
      !--              ef              --
      !----------------------------------
      !--
      !-- Equation 2 in Ek and Holtslag 2004
      !--
      where( lhf0.ne.missing .and. shf0.ne.missing .and. lhf0+shf0.ne.0 )   ef   =   lhf0 / (shf0+lhf0)      


      !-----------------------------
      !--                         --
      !-- Non-Evaporative Section --
      !--            ne           --
      !-----------------------------
      !--
      !-- Calculate the C_1 and C_2 term from Equation A8 from Ek and Holtslag 2004
      !--
      where( pbl_qs.ne.missing  .and.  pbl_t.ne.missing  .and.  psfc.ne.missing  .and.  pbl_p.ne.missing )
         c_1  =   (Lv/Rv)  *  (pbl_qs/(pbl_t**2))  *  (pbl_p/psfc)**R_cp
      endwhere

      where( pbl_qs.ne.missing  .and.  pbl_t.ne.missing                                                  )
         c_2  =  ((Lv/Rv   *  (pbl_qs/(pbl_t**2))) -  (cp/Rd * pbl_qs/pbl_t))  *  (grav/cp)
      endwhere


      !--
      !-- Calculate vertical gradient of potential temperature above PBL - g_theta [K]
      !--                         ***  AND  ***
      !-- Calculate change in specific humidity across the PBL - delta_q [kg/kg]
      !--
      iabove_pbl   =   minloc( hgt, DIM = 2, MASK = hgt.gt.pblh_2d .and. hgt.ne.missing )
      ibelow_pbl   =   maxloc( hgt, DIM = 2, MASK = hgt.lt.pblh_2d .and. hgt.ne.missing )
      do tt=1,ntim
         if( iabove_pbl(tt).eq.0 .or. ibelow_pbl(tt).eq.0 ) then
            cycle
         end if
         if( pot_k(tt,iabove_pbl(tt)).ne.missing  .and.  pbl_theta(tt).ne.missing ) then
             g_theta(tt)  =  pot_k(tt,iabove_pbl(tt)) - pbl_theta(tt)
         end if
         if( qhum(tt,iabove_pbl(tt)).ne.missing  .and.  qhum(tt,ibelow_pbl(tt)).ne.missing ) then
             delta_q(tt)  =   qhum(tt,iabove_pbl(tt)) - qhum(tt,ibelow_pbl(tt))
         end if
      end do


      !--
      !-- Equation 3 in Ek and Holtslag 2004
      !--
      !*** Dry air entrainment term
      where( delta_q.ne.missing  .and.  pblh_mid.ne.missing  .and.  &
             g_theta.ne.missing  .and.  pblh_mid*g_theta.ne.0       )
          dry_entrain  =  delta_q/(pblh_mid*g_theta)
      endwhere

      !*** PBL growth term
      where( relh   .ne.missing  .and.  c_2  .ne.missing  .and.  &
             g_theta.ne.missing  .and.  g_theta.ne.0             )
          pbl_growth   =  relh * c_2 / g_theta
      endwhere

      !*** PBL heating term
      where( relh   .ne.missing  .and.  c_1  .ne.missing         )  
          pbl_heating  =  relh * c_1
      endwhere

      !*** Non-evaporative term (sum of the above terms multiplied by an entrainment ratio)
      where( C_theta   .ne.missing  .and.  dry_entrain.ne.missing  .and.  &
             pbl_growth.ne.missing  .and.  pbl_heating.ne.missing         )
          ne   =   Lv_cp * (1. + C_theta) * ( dry_entrain  +  pbl_growth  -  pbl_heating )
      endwhere

      !-- PBL heating contribution to Change in RH (unitless)
      where( C_theta.ne.missing  .and.  pbl_heating.ne.missing )
          pbl_heating   =   Lv_cp * (1. + C_theta) * ( -1. *  relh * c_1 )
      endwhere

      !-- PBL growth contribution to Change in RH (unitless)
      where( C_theta.ne.missing  .and.  pbl_growth.ne.missing )
          pbl_growth   =   Lv_cp * (1. + C_theta) * ( relh * c_2 / g_theta )
      endwhere

      !-- Dry air entrainment contribution to Change in RH (unitless)
      where( C_theta.ne.missing  .and.  dry_entrain.ne.missing )
          dry_entrain   =   Lv_cp * (1. + C_theta) * ( delta_q/(pblh_mid*g_theta) )
      endwhere


      !*** Estimate Relative Humidity Change (%/hr)
      where( ef .ne.missing  .and.  ne  .ne.missing  .and.  pblh_mid.ne.missing  .and.  &
             rho.ne.missing  .and.  lhf0.ne.missing  .and.  pbl_qs  .ne.missing         )
         dRH_estimate  =  ((shf0+lhf0) / (rho*Lv*pblh_mid*pbl_qs))  *  (ef + ne*(1.-ef))  *  3600. * 1e2
      endwhere


      return

end subroutine rh_tend_calc



end module RH_Tend_Mod
