
!---------------------------------------------------------------------------------------------
!
!  Purpose:  This program runs a quick sample data test against several observed profiles from 
!            the Intensive Observation Period on June 6th, 2002 from the Southern Great Plains
!            Atmosphere Radiation Radiation Central Facility. 
!            **** NOTE --  This is NOT a generic program but specific to this input file structure
!
!

Program Coupling_metrics

       use RH_Tend_Mod
       use Mixing_Diag_Mod
       use HCF_vars_calc
       implicit none

   !   
   ! Local variables 

   ! dimension sizes
   integer, parameter                        ::  nlev  =  23
   integer, parameter                        ::  ntim  =  8

   ! Loop indices
   integer                                   ::  tt, zz

   ! file name prefix for each file that is looped over
   character(len=256)                        ::  input_profile, input_fluxes
   integer                                   ::  unit_profile , unit_fluxes

   ! Input variables read in from file 
   real(4), dimension(ntim,nlev)             ::  tlev, plev, hlev, qlev
   real(4), dimension(ntim,1)                ::  t2m , p2m , h2m , q2m
   real(4), dimension(ntim,1)                ::  shf , lhf , pblh


   real(4), dimension(ntim)                  :: tbm    ,  bclh   , bclp ,  tdef,    &
                                                hadv   ,  tranp  , tadv ,           &
                                                shdef  ,  lhdef  , eadv
   real(4), dimension(1,1)                   :: sh_ent, lh_ent, sh_sfc, lh_sfc, sh_tot,  lh_tot
   real(4), dimension(ntim)                  :: ef, ne, heating, growth, dry, drh_dt

   real(4)                                   :: yr1, mn1, hr1
   real(4)                                   :: yr2, mn2, hr2
   real(4)                                   :: yr3, mn3, hr3
   real(4)                                   :: yr4, mn4, hr4
   real(4)                                   :: yr5, mn5, hr5
   real(4)                                   :: yr6, mn6, hr6
   real(4)                                   :: yr7, mn7, hr7

   ! missing value read in from file
   real(4), parameter                        ::  missing = -9999.

   ! Format statements
   character(len=24), parameter :: FMT1 = "(4(F12.4,2x))"   
   character(len=24), parameter :: FMT2 = "(4(A,2x))"   


!---------------------------------------------------------------------------------

         !********************************  
         !***                               
         !***    Initialize inputs to missing values
         !***                               
         !********************************  
         plev  =  missing
         tlev  =  missing
         qlev  =  missing
         hlev  =  missing
         t2m   =  missing
         q2m   =  missing
         h2m   =  missing
         p2m   =  missing

         shf   =  missing
         lhf   =  missing
         pblh  =  missing

         tbm   =  missing
         tdef  =  missing
         shdef =  missing
         bclp  =  missing
         bclh  =  missing
         lhdef =  missing
         eadv  =  missing
         hadv  =  missing
         tranp =  missing
         tadv  =  missing
         sh_ent  =  missing
         lh_ent  =  missing
         sh_sfc  =  missing
         lh_sfc  =  missing
         sh_tot  =  missing
         lh_tot  =  missing

         ef      =  missing
         ne      =  missing
         heating =  missing
         growth  =  missing
         dry     =  missing
         drh_dt  =  missing


         !********************************  
         !***                               
         !***    Read sample profile data
         !***                               
         !********************************  
         unit_profile  =  10
         open( unit=unit_profile, file='Sample_profile.txt' )
         do zz=1,nlev
            read(unit_profile,*) yr1,mn1,hr1,plev(1,zz),hlev(1,zz),tlev(1,zz),qlev(1,zz),  &
                                 yr2,mn2,hr2,plev(2,zz),hlev(2,zz),tlev(2,zz),qlev(2,zz),  &
                                 yr3,mn3,hr3,plev(3,zz),hlev(3,zz),tlev(3,zz),qlev(3,zz),  &
                                 yr4,mn4,hr4,plev(4,zz),hlev(4,zz),tlev(4,zz),qlev(4,zz),  &
                                 yr5,mn5,hr5,plev(5,zz),hlev(5,zz),tlev(5,zz),qlev(5,zz),  &
                                 yr6,mn6,hr6,plev(6,zz),hlev(6,zz),tlev(6,zz),qlev(6,zz),  &
                                 yr7,mn7,hr7,plev(7,zz),hlev(7,zz),tlev(7,zz),qlev(7,zz)

         end do
         close(unit_fluxes)

         !********************************  
         !***                               
         !***    Read sample flux and surface data
         !***                               
         !********************************  
         unit_fluxes  =  11
         open( unit=unit_fluxes, file='Sample_fluxes.txt' )
         do tt=1,ntim-1
            read(unit_fluxes,*) p2m(tt,1),h2m(tt,1),t2m(tt,1),q2m(tt,1),shf(tt,1),lhf(tt,1),pblh(tt,1)
         end do
         close(unit_fluxes)
         where( p2m .ne.missing )  p2m  =  p2m  * 1e2
         where( pblh.ne.missing )  pblh =  pblh * 1e3

write(*,*)
write(*,*)
write(*,*)
write(*,*) "Sensible Heat Fluxes"
write(*,*) shf
write(*,*)
write(*,*)
write(*,*)
write(*,*) "Latent Heat Fluxes"
write(*,*) lhf
write(*,*)
write(*,*)
write(*,*)
write(*,*) "Boundary Layer Height"
write(*,*) pblh/1e3
write(*,*)
write(*,*)
write(*,*)
write(*,*) "2-m Temperature"
write(*,*) t2m
write(*,*)
write(*,*)
write(*,*)
write(*,*) "2-m Specific Humidity"
write(*,*) q2m
write(*,*)
write(*,*)
write(*,*)

         !---------------------------------------------
         !---
         !--- Heated Condensation Section
         !---
         !---------------------------------------------
         !**********************************    
         !*** Loop over time
         !**********************************    
         do tt = 1,ntim

                    !---------------------------------------------
                    !--- Call HCF variables calculate routines
                    !---------------------------------------------
                    call hcfcalc( nlev       ,  missing                            ,    &
                                  tlev (tt,:),  plev (tt,:), qlev(tt,:),  hlev(tt,:),    &
                                  t2m  (tt,1),  p2m  (tt,1), q2m (tt,1),  h2m (tt,1),    &
                                  tbm  (tt)  ,  bclh (tt)  , bclp(tt)  ,  tdef(tt)  ,    &
                                  hadv (tt)  ,  tranp(tt)  , tadv(tt)  ,                 &
                                  shdef(tt)  ,  lhdef(tt)  , eadv(tt)                    )


         end do  !*** end of time loop



         !---------------------------------------------
         !---
         !--- Mixing Diagrams Section
         !---
         !---------------------------------------------
         call mixing_diag_daily ( 1          , ntim       , 1          ,  8       ,          &  
                                  t2m (:,1:1), p2m(:,1:1) , q2m(:,1:1) ,                     &
                                  pblh(:,1:1), shf(:,1:1) , lhf(:,1:1) ,  8.*3600.,          &
                                  sh_ent     , lh_ent     ,                                  &
                                  sh_sfc     , lh_sfc     , sh_tot     ,  lh_tot,  missing )

         !---------------------------------------------
         !---
         !--- Relative Humidity Tendency Section
         !---
         !---------------------------------------------
         call rh_tend_calc ( nlev     ,  ntim   ,                          &
                             tlev     ,  qlev   ,  hlev   ,  plev    ,     &
                             pblh     ,  shf    ,  lhf    ,  8.*3600.,     &
                             ef       ,  ne     ,                          &
                             heating  ,  growth ,  dry    ,  drh_dt , missing  )



         write(*,*)
         write(*,*)
         write(*,*)  "  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  "
         write(*,*)  "  !!!!!!!!!!!      Heated Condensation Output      !!!!!!!!!!!!  "
         write(*,*)  "  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  "
         write(*,FMT2) "       TDEF",    "        TBM","   BCL Height",    "  BCL Pressure"
         do tt=1,ntim-1
            write(*,FMT1)  tdef(tt), tbm(tt), bclh(tt)/1e3, bclp(tt)/1e2
         end do
         write(*,*)
         write(*,FMT2)  "       SHDEF","         LHDEF","    Energy Advantage","   Transition Pressure"
         do tt=1,ntim-1
            write(*,FMT1)  shdef(tt)/1e6, lhdef(tt)/1e6, eadv(tt), tranp(tt)/1e2
         end do
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*)

         write(*,*)  "  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  "
         write(*,*)  "  !!!!!!!!!!!      Mixing Diagrams Output      !!!!!!!!!!!!  "
         write(*,*)  "  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  "
         write(*,FMT2) " SH Entrainment",    " SH Surface", "   SH Total",    "Entrainment Ratio"
         write(*,FMT1)  sh_ent, sh_sfc, sh_tot, sh_ent/sh_sfc
         write(*,*) 
         write(*,FMT2) " LH Entrainment",    " LH Surface", "   LH Total",    "Entrainment Ratio"
         write(*,FMT1)  lh_ent, lh_sfc, lh_tot, lh_ent/lh_sfc
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*)  "  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
         write(*,*)  "  !!!!!!!!!!!      RH Tendency Output      !!!!!!!!!!!!"
         write(*,*)  "  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
         write(*,FMT2) "       EF", "           NE", "      Estimate dRH/dt"
         do tt=1,ntim-1
            write(*,FMT1)  ef(tt), ne(tt), drh_dt(tt)
         end do
         write(*,*)
         write(*,FMT2)  "  PBL Drying","  PBL Heating","  PBL Growth"
         do tt=1,ntim-1
            write(*,FMT1)  dry(tt), heating(tt), growth(tt)
         end do
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*)


end Program Coupling_metrics
