!---------------------------------------------------------------------------------
! Purpose:
!
! !DESCRIPTION:  
! This subroutine calculates the full suite of mixing diagram variables when
! given a time series of of surface fluxes, surface 2m state variables, and 
! boundary layer height. Mixing diagram variables returned are surface bowen
! ratio, entrainment bowen ratio, advection ratios, and the various latent and 
! sensible heat fluxes associated with those ratios. The utility is to characterize 
! the coupling between surface fluxes and top of boundary layer fluxes in tandem
! with knowledge regarding soil moisture state.  More details regarding motivation
! can be found in the references below.
! 
!  References: Santanello et al. 2009,  A Modeling and Observational Framework 
!              for Diagnosing Local Land-Atmosphere Coupling on Diurnal Time Scales
!
!              Santanello et al. 2011,  Diagnosing the Sensitivity of Local 
!              Land-Atmosphere Coupling via the Soil Moisture-Boundary Layer Interaction
!
! Author and Revision History: 
! Original Author of NCL scripts -- Jatin Kala  on Nov 2013
! Converted to F90 module        -- A.B. Tawfik on Apr 2015
!
!---------------------------------------------------------------------------------
module Mixing_Diag_Mod

     !
     ! subroutine name 
     !
     public mixing_diag_daily

!---------------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------------



!---------------------------------------------------------------------------------
!
! subroutines:  calculates mixing diagram variables using the INCREMENTAL more precise approach
!               Assumes left most dimension is time
!               Output is dimensioned:  
!               (days, spatial dimension like lat/lon)
!---------------------------------------------------------------------------------
  subroutine mixing_diag_daily( dim2      , ntim    , nday   ,  nhr                               ,  &
                                t2m       , psfc    , q2m    ,  pbl_h  , shf    , lhf    , dt     ,  &
                                shf_ent   , lhf_ent , shf_sfc,  lhf_sfc, shf_tot, lhf_tot, missing   )
   implicit none

!
! Input/Output Variables
!
   integer, intent(in   )                         ::  dim2           ! *** missing value - useful for obs
   integer, intent(in   )                         ::  ntim           ! *** total number of time slices
   integer, intent(in   )                         ::  nday           ! *** # of days
   integer, intent(in   )                         ::  nhr            ! *** hours of the time increment per day [hours]
   real(4), intent(in   )                         ::  dt             ! *** time increment of the hours dimension [seconds]
 
   real(4), intent(in   )                         ::  missing        ! *** missing value - useful for obs
   real(4), intent(in   ), dimension(ntim,dim2)   ::  t2m, q2m, psfc ! *** 2m quantities - (days,hours) [K], [kg/kg], [Pa]
   real(4), intent(in   ), dimension(ntim,dim2)   ::  pbl_h          ! *** boundary layer height - (days,hours) [m]
   real(4), intent(in   ), dimension(ntim,dim2)   ::  shf, lhf       ! *** surface fluxes- dimensioned (days,hours) [W/m2]
   real(4), intent(out  ), dimension(nday,dim2)   ::  shf_ent        ! *** entrainment flux of sensible heat [W/m2]
   real(4), intent(out  ), dimension(nday,dim2)   ::  lhf_ent        ! *** entrainment flux of latent heat [W/m2]
   real(4), intent(out  ), dimension(nday,dim2)   ::  shf_sfc        ! *** surface flux of sensible heat [W/m2]
   real(4), intent(out  ), dimension(nday,dim2)   ::  lhf_sfc        ! *** surface flux of latent heat [W/m2]
   real(4), intent(out  ), dimension(nday,dim2)   ::  shf_tot        ! *** Total flux of sensible heat [W/m2]
   real(4), intent(out  ), dimension(nday,dim2)   ::  lhf_tot        ! *** Total flux of latent heat [W/m2]
!
! Local variables
!
   real(4), parameter   ::  p_ref=1e5 , Lv=2.5e6, cp=1005.7, R_cp=287.04/1005.7
   real(4), parameter   ::  grav = 9.81, Rd=287.04, ep = 0.622
 
   integer                           ::  tt, yy
   integer, dimension(dim2*nday)     ::  my_good
   real(4), dimension(ntim,dim2)     ::  rho                  
   real(4), dimension(ntim,dim2)     ::  bowen_s             
   real(4), dimension(ntim,dim2)     ::  cp_theta_final      
   real(4), dimension(ntim,dim2)     ::  cp_theta_initial  
   real(4), dimension(ntim,dim2)     ::  cp_theta           
   real(4), dimension(ntim,dim2)     ::  cp_deltaT          
   real(4), dimension(ntim,dim2)     ::  Lv_qhum_initial      
   real(4), dimension(ntim,dim2)     ::  Lv_qhum_final     
   real(4), dimension(ntim,dim2)     ::  Lv_qhum         
   real(4), dimension(ntim,dim2)     ::  Lv_q_0          
   real(4), dimension(ntim,dim2)     ::  cp_T_0             

   real(4), dimension(nhr,dim2*nday) ::  lhf0             
   real(4), dimension(nhr,dim2*nday) ::  shf0           
   real(4), dimension(nhr,dim2*nday) ::  lhfi0           
   real(4), dimension(nhr,dim2*nday) ::  shfi0          
   real(4), dimension(nhr,dim2*nday) ::  lhfs0      
   real(4), dimension(nhr,dim2*nday) ::  shfs0     
 
   real(4), dimension(ntim,dim2)     ::  lhf00             
   real(4), dimension(ntim,dim2)     ::  shf00           
   real(4), dimension(ntim,dim2)     ::  lhfi00           
   real(4), dimension(ntim,dim2)     ::  shfi00          
   real(4), dimension(ntim,dim2)     ::  lhfs00      
   real(4), dimension(ntim,dim2)     ::  shfs00     
 
   real(4), dimension(dim2*nday)     ::  lhf_tots
   real(4), dimension(dim2*nday)     ::  shf_tots
   real(4), dimension(dim2*nday)     ::  lhf_ents
   real(4), dimension(dim2*nday)     ::  shf_ents
   real(4), dimension(dim2*nday)     ::  lhf_sfcs
   real(4), dimension(dim2*nday)     ::  shf_sfcs

!-----------------------------------------------------------------------------


      !--------------------------------------------
      !--- Initialization and preliminary calculations
      !--------------------------------------------
      !----------------------------------------------------------------------------------
      !-- initialize output arrays
      !----------------------------------------------------------------------------------
      shf_ent    =  missing
      lhf_ent    =  missing
      shf_sfc    =  missing
      lhf_sfc    =  missing
      shf_tot    =  missing
      lhf_tot    =  missing


      !-----------------------------------------------------------------------------
      !-- Check input arrays to make sure all are NOT missing
      !-----------------------------------------------------------------------------
      if( all(t2m  .eq.missing)   .or.  all(q2m  .eq.missing)    .or.   &
          all(psfc .eq.missing)   .or.  all(pbl_h.eq.missing)    .or.   &
          all(lhf  .eq.missing)   .or.  all(shf  .eq.missing)           )  then
          return
      end if

    
      !----------------------------------------------------------------------------------
      !-- initialize working arrays
      !----------------------------------------------------------------------------------
      bowen_s    =  missing
      rho        =  missing
      Lv_qhum    =  missing
      cp_theta   =  missing
      shf0       =  missing
      lhf0       =  missing
      shfi0      =  missing
      lhfi0      =  missing
      shfs0      =  missing
      lhfs0      =  missing

      shf_ents   =  missing
      lhf_ents   =  missing
      shf_sfcs   =  missing
      lhf_sfcs   =  missing
      shf_tots   =  missing
      lhf_tots   =  missing
    
      cp_T_0     =  missing
      Lv_q_0     =  missing

      Lv_qhum_initial   =  missing
      Lv_qhum_final     =  missing
      cp_theta_initial  =  missing
      cp_theta_final    =  missing
      my_good           =  0

      !----------------------------------------------------------------------------------
      !-- Approximate air density (kg/m3)
      !----------------------------------------------------------------------------------
      where( t2m.ne.missing .and. psfc.ne.missing .and. q2m.ne.missing ) 
         rho  =  psfc / (Rd * t2m * ((1. + (q2m/ep)) / (1. + q2m)))
      endwhere

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
      where( lhf.ne.0 .and. shf.ne.missing .and. lhf.ne.missing )  bowen_s  =  shf/lhf

      
      !--------------------------------------------
      !--- Proceed to calculating output variables
      !--------------------------------------------
      !----------------------------------------------------------------------------------
      !--
      !-- Surface heat and moisture vector components (see Eq 1 in Santanello et al. 2009)
      !--
      !----------------------------------------------------------------------------------
      where( shf.ne.missing .and. rho.ne.missing .and. pbl_h.ne.missing .and. rho*pbl_h.ne.0 ) 
         cp_deltaT  =  (shf*dt) / (rho*pbl_h)
      elsewhere 
         cp_deltaT  =  missing
      endwhere

      !----------------------------------------------------------------------------------
      !-- Define the start and end times of the increment for sensible and latent heat
      !----------------------------------------------------------------------------------
      cp_theta_initial           =   cp_theta        !*** start of time increment
      cp_theta_final(:ntim-1,:)  =   cp_theta(2:,:)  !*** end of time increment
      Lv_qhum_initial            =   Lv_qhum         !*** start of time increment
      Lv_qhum_final (:ntim-1,:)  =   Lv_qhum (2:,:)  !*** end of time increment

      where( Lv_qhum_initial.ne.missing .and. cp_deltaT.ne.missing .and. bowen_s.ne.missing )
             Lv_q_0  =  Lv_qhum_initial   +  (cp_deltaT/bowen_s)
      end where
      where( cp_theta_initial.ne.missing .and. cp_deltaT.ne.missing )
             cp_T_0  =  cp_theta_initial  +  cp_deltaT
      end where


      !--------------------------------------------------
      ! output arrays
      !--------------------------------------------------
      !**************************************************
      !*****
      !*****  Output daily average quantities 
      !*****
      !**************************************************

          !--------------------------------------------------
          ! Calculate Heat budget in Wm-2
          !--------------------------------------------------
          where( cp_theta_final.ne.missing .and. rho.ne.missing .and. pbl_h.ne.missing .and. cp_T_0.ne.missing ) 
             shf00   =  (((cp * ((cp_theta_final/cp)  -  &
                          (cp_theta_initial/cp)))) * (rho * pbl_h))/(dt)
             shfi00  =  (((cp * ((cp_theta_final/cp)  -  &
                          (cp_T_0          /cp)))) * (rho * pbl_h))/(dt)
          endwhere
          shf0   =  reshape( shf00 , (/nhr,dim2*nday/) )
          shfi0  =  reshape( shfi00, (/nhr,dim2*nday/) )

          !--------------------------------------------------
          !****  Total sensible heat Wm-2
          !--------------------------------------------------
          my_good    =   count( shf0.ne.missing, DIM = 1 )
          where( my_good.gt.0 )
             shf_tots =   sum(shf0, DIM = 1, MASK = shf0.ne.missing)    /  my_good
          endwhere
          shf_tot  =  reshape( shf_tots, (/nday,dim2/) )      

          !--------------------------------------------------
          !****  Entrainment sensible heat Wm-2
          !--------------------------------------------------
          my_good     =   count( shfi0.ne.missing, DIM = 1 )
          where( my_good.gt.0 )
             shf_ents =   sum(shfi0, DIM = 1, MASK = shfi0.ne.missing)  /  my_good
          endwhere
          shf_ent  =  reshape( shf_ents, (/nday,dim2/) )      

          !--------------------------------------------------
          !****  Surface sensible heat Wm-2
          !--------------------------------------------------
          where( cp_theta_initial.ne.missing .and. rho.ne.missing .and. pbl_h.ne.missing .and. cp_T_0.ne.missing ) 
             shfs00   =  (((cp * ((cp_T_0        /cp)  -  &
                           (cp_theta_initial/cp)))) * (rho * pbl_h))/(dt)
          endwhere
          shfs0    =   reshape( shfs00, (/nhr,dim2*nday/) )
          my_good  =   count  ( shfs0.ne.missing, DIM = 1 )
          where( my_good.gt.0 )
             shf_sfcs =  sum(shfs0, DIM = 1, MASK = shfs0.ne.missing)   /  my_good
          endwhere
          shf_sfc  =  reshape( shf_sfcs, (/nday,dim2/) )      


          !--------------------------------------------------
          !Calculate Moisture Budget in Wm-2
          !--------------------------------------------------
          where( Lv_qhum_final.ne.missing .and. rho.ne.missing .and. pbl_h.ne.missing .and. Lv_q_0.ne.missing ) 
             lhf00   =  (((Lv * ((Lv_qhum_final/Lv)  -  &
                          (Lv_qhum_initial/Lv)))) * (rho * pbl_h))/(dt)
             lhfi00  =  (((Lv * ((Lv_qhum_final/Lv)  -  &
                          (Lv_q_0         /Lv)))) * (rho * pbl_h))/(dt)
          endwhere
          lhf0    =  reshape( lhf00 , (/nhr,dim2*nday/) )
          lhfi0   =  reshape( lhfi00, (/nhr,dim2*nday/) )

          !--------------------------------------------------
          !****  Total latent heat Wm-2
          !--------------------------------------------------
          my_good     =   count( lhf0.ne.missing, DIM = 1 )
          where( my_good.gt.0 )
             lhf_tots =   sum(lhf0, DIM = 1, MASK = lhf0.ne.missing)   /  my_good
          endwhere
          lhf_tot  =  reshape( lhf_tots, (/nday,dim2/) )      

          !--------------------------------------------------
          !****  Entrainment latent heat Wm-2
          !--------------------------------------------------
          my_good     =   count( lhfi0.ne.missing, DIM = 1 )
          where( my_good.gt.0 )
             lhf_ents =   sum(lhfi0, DIM = 1, MASK = lhfi0.ne.missing) /  my_good
          endwhere
          lhf_ent  =  reshape( lhf_ents, (/nday,dim2/) )      

          !--------------------------------------------------
          !****  Surface latent heat Wm-2
          !--------------------------------------------------
          where( Lv_qhum_initial.ne.missing .and. rho.ne.missing .and. pbl_h.ne.missing .and. Lv_q_0.ne.missing ) 
             lhfs00   =  (((Lv * ((Lv_q_0        /Lv)  -  &
                           (Lv_qhum_initial/Lv)))) * (rho * pbl_h))/(dt)
          endwhere
          lhfs0   =  reshape( lhfs00, (/nhr,dim2*nday/) )

          my_good    =   count( lhfs0.ne.missing, DIM = 1 )
          where( my_good.gt.0 )
             lhf_sfcs  =   sum(lhfs0, DIM = 1, MASK = lhfs0.ne.missing)  /  my_good
          endwhere
          lhf_sfc  =  reshape( lhf_sfcs, (/nday,dim2/) )      



      return

end subroutine mixing_diag_daily



end module Mixing_Diag_Mod
