
!---------------------------------------------------------------------------------------------
!
!  Purpose:  Interface that calculates several coupling metrics from the Atmosphere Radiation 
!            Measurement Southern Great Plains Best Estimate sites (ARM-BE SGP) hourly data.
!            Specifically the numerical weather prediction output is used for the profiles 
!            from the Rapid update Cycle (RUC) 1-hourly output as well as surface flux and 
!            state variable observed measurements at 2-meters. 
!            This is just used as a sample dataset and currently has been tested with the
!            intel fortran compiler (ifort). This requires netcdf libraries to function.
!
!
!

Program Coupling_metrics

       use RH_Tend_Mod
       use Mixing_Diag_Mod
       use HCF_vars_calc
       use netcdf
       use IFPORT
       implicit none
       include 'netcdf.inc'

   !   
   ! Local variables 
   !   
   ! Output variables after calculation 
   real(4), dimension(:)  , allocatable  ::  BCL_P, BCL_H, TBM, TDEF
   real(4), dimension(:)  , allocatable  ::  T_ADV45, H_ADV45, SHDEF_MIX, LHDEF_MIX, EADV_MIX
   real(4), dimension(:)  , allocatable  ::  SH_ENT, LH_ENT, SH_SFC, LH_SFC, SH_TOT, LH_TOT
   real(4), dimension(:)  , allocatable  ::  GROWTH, HEATING, DRYING, EF, NE, dRH_dt

   ! Input variables for calculation
   real(4), dimension(:,:), allocatable    ::  TK , QK , HGT
   real(4), dimension(:)  , allocatable    ::  T2M, Q2M, PSFC
   real(4)                                 ::  H2M
   real(4), dimension(:)  , allocatable    ::  T2M0, Q2M0, PSFC0, SH0, LH0, PBLH0

   ! Temporary variables used in Relative Humidity Tendency Calculations
   real(4), dimension(:,:), allocatable    ::  TK0, QK0, HGT0
   real(4), dimension(:,:), allocatable    ::  hgt2, tmp2, qhum2, press2
   real(4), dimension(:)  , allocatable    ::  shf2, lhf2, pbl2
   real(4), dimension(:)  , allocatable    ::  ef0, ne0, heating0, growth0, dry0, drh_dt0

   ! number of dimensions; known due to the interface selection
   integer, parameter                        ::  ndims = 1

   ! dimension sizes
   integer, parameter                        ::  nlev      =  37
   integer, parameter                        ::  hr_p_day  =  24
   real(4), parameter                        ::  dt        =  1. * 3600.  !*** [seconds] of hourly step
   integer                                   ::  ntim
   integer, parameter                        ::  nlev1     =  nlev + 1

   ! pressure levels (static)
   real(4), dimension(nlev)                  ::  plev

   ! Loop indices
   integer                                   ::  ii, tt, ff, cc, zz, dd, xx, yy, kk

   ! Bowen ratios coming out of mixing diagrams (output from mixing diag function)
   real(4), dimension(:,:)    , allocatable    ::  a_sh0, a_lh0
   real(4), dimension(:,:)    , allocatable    ::  sh_tot0, lh_tot0
   real(4), dimension(:,:)    , allocatable    ::  sh_ent0, lh_ent0
   real(4), dimension(:,:)    , allocatable    ::  lh_sfc0, sh_sfc0

   ! file name prefix for each file that is looped over
   character(len=256)                        ::  input0

   ! number of files to loop over
   integer                                   ::  nfiles

   ! missing value read in from file
   real(4), parameter                        ::  missing = -9999.

   ! dummy variables for bottom and top pressure levels used to check correct bottom-up order
   real(4), dimension(nlev)                  ::  qhum, tmpk, pres, ht

   ! I/O variables
   character(len=256)                        ::  NAME1, lname, units
   integer                                   ::  ncid_in, ncido
   integer                                   ::  dimid(ndims), dimid1, dimid2, dimid3, id_in, alldims(ndims)
   integer                                   ::  varid0, varid1, varid2, varid3, varid4, varid_time
   integer                                   ::  varid5, varid6, varid7, varid8

   ! Name of files that are to operated on
   character(len=256), allocatable           ::  in_files(:)

   ! Input variables names
   character(len=256), dimension(10)         ::  vnames


   ! Specific humidity deficit at the PBL [kg/kg]
   real(4)                                   ::  pbl_qdef, tran_p

   ! Output directory
   character(len=256)                        ::  out_dir, out_name, out_dir1, out_dir2

   ! Output choices
   logical                                   ::  MIXING_OUT, HCF_OUT, RH_OUT

   ! Input directory
   character(len=256)                        ::  inpath

   ! Command line execution commands and output stats
   character(len=10)                         ::  init_date
   integer                                   ::  istat

   ! Date integers
   real(8), allocatable                      ::  time_stamp(:)

   ! Temporary subroutine output variables   
   real(4), dimension(nlev+1)                ::  pmlt_0, pmlq_0, pmlqs_0, shdef_0, lhdef_0, eadv_0
   real(4), dimension(:,:), allocatable      ::  t2m_in, q2m_in, psf_in, shf_in, lhf_in, pbl_in

   ! Number of days per month 
   integer                                   ::  ndays

!---------------------------------------------------------------------------------


         !********************************  
         !***                               
         !***    Define dimensions and output
         !***                               
         !********************************  
         HCF_OUT     =  .true.
         MIXING_OUT  =  .true.
         RH_OUT      =  .true.

         out_name  =  "ARMBE"
         out_dir   =  "/glade/scratch/abtawfik/ARMBE/HCF_output/"
         out_dir1  =  "/glade/scratch/abtawfik/ARMBE/MIXING_output/"
         out_dir2  =  "/glade/scratch/abtawfik/ARMBE/RH_TEND_output/"

         vnames(1)  =  "T_nwp_p"
         vnames(2)  =  "qhum_nwp"
         vnames(3)  =  "hgt_nwp"
         vnames(4)  =  "p"
         vnames(5)  =  "T_sfc"
         vnames(6)  =  "q_sfc"
         vnames(7)  =  "p_sfc"
         vnames(8)  =  "SH"
         vnames(9)  =  "LH"
         vnames(10) =  "pblh_nwp"


         !********************************  
         !***                               
         !***    Read Command line arguments
         !***                               
         !********************************  
         nfiles      =  command_argument_count()
         write(*,*)  "  Total number of files   =   ",nfiles
         if(.not.allocated(in_files)) allocate(in_files(nfiles))
         in_files  =  'FAKE'
         do cc = 1, nfiles
            call get_command_argument(cc, in_files(cc) )
            write(*,*)  cc,"   =   ",trim(in_files(cc))
         end do





         !-----------------------------------------------------
         !--- 
         !--- Loop over all files in command line
         !--- 
         !-----------------------------------------------------         
         !*****************************
         !*** Sequence below is:
         !*** 1) Open netcdf input files
         !*** 2) Perform dimension inquiry
         !*** 3) Get the varids for netCDF variables.
         !*** 4) Read the data from the file.
         !*** 5) Close netcdf files
         !*** 6) Perform the COUPLING METRIC calculation
         !*** 7) Open netcdf output file
         !*** 8) Write Output variables to netcdf file
         !*****************************
         do ff = 1,nfiles

            !-- Open netcdf file
            call check( nf_open(trim(in_files(ff)), nf_nowrite, ncid_in) )

            !-- Input file name prefix (e.g. everything before the file extension)
            write(*,*) "-----    Current File   -----  ", trim(in_files(ff))

            !-- Get the month and year
            istat =  system('echo `(basename '//trim(in_files(ff))//') | cut -d "." -f 3` > '//trim(in_files(ff))//'.txt')
            open ( unit=10, position="rewind", file=trim(in_files(ff))//".txt")
            read ( 10, *) init_date
            close( 10)

            !-- Get the file path for input files
            !istat =  system('echo `(dirname '//trim(in_files(ff))//')` > '//trim(in_files(ff))//'1.txt')
            !open ( unit=11, position="rewind", file=trim(in_files(ff))//"1.txt")
            !read ( 11, '(A)') inpath
            !close( 11)


            !**********************************    
            !*** Allocate output arrays
            !**********************************    
            !--- Get number of time stamps
            call check( nf_inq_dimid (ncid_in, 'time', id_in) ) 
            call check( nf_inq_dimlen(ncid_in, id_in , ntim ) )
            if(.not.allocated(time_stamp)) allocate(time_stamp(ntim))
            call check( nf_inq_varid     (ncid_in, 'time', id_in      ) ) 
            call check( nf_get_var_double(ncid_in,  id_in, time_stamp ) )
            write(*,*) " ----- THIS IS THE DATE ----- ",trim(init_date),ntim
            call check( nf_close(ncid_in) )
            ndays    =  ntim/hr_p_day


            !**************************************
            !**************************************
            !**************************************
            !****                              ****
            !****         HCF Section          ****
            !****                              ****
            !**************************************
            !**************************************
            !**************************************
            if( HCF_OUT ) then 


                write(*,*) 
                write(*,*) 
                write(*,*) 
                write(*,*)  "  !!!!!!!!!!!      Heated Condensation Begins      !!!!!!!!!!!!  "
                !**********************************    
                !*** Allocate input arrays
                !**********************************    
                if(.not.allocated(TK  ))  allocate(TK  (nlev,ntim) )
                if(.not.allocated(QK  ))  allocate(QK  (nlev,ntim) ) 
                if(.not.allocated(HGT ))  allocate(HGT (nlev,ntim) )
                if(.not.allocated(T2M ))  allocate(T2M (ntim)      )
                if(.not.allocated(Q2M ))  allocate(Q2M (ntim)      )
                if(.not.allocated(PSFC))  allocate(PSFC(ntim)      )

                plev  =  (/ 1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, &
                             750, 725, 700, 675, 650, 625, 600, 575, 550, 525, &
                             500, 475, 450, 425, 400, 375, 350, 325, 300, 275, &
                             250, 225, 200, 175, 150, 125, 100 /)  * 1e2

                !-- Assume reference height is 2-meters
                H2M  =  2.

                if(.not.allocated(TBM  )) allocate(TBM  (ntim))
                if(.not.allocated(TDEF )) allocate(TDEF (ntim))
                if(.not.allocated(BCL_P)) allocate(BCL_P(ntim))
                if(.not.allocated(BCL_H)) allocate(BCL_H(ntim))

                if(.not.allocated(T_ADV45   )) allocate( T_ADV45    (ntim))
                if(.not.allocated(H_ADV45   )) allocate( H_ADV45    (ntim))
                if(.not.allocated(SHDEF_MIX )) allocate( SHDEF_MIX  (ntim))
                if(.not.allocated(LHDEF_MIX )) allocate( LHDEF_MIX  (ntim))
                if(.not.allocated(EADV_MIX  )) allocate( EADV_MIX   (ntim))

                !**********************************    
                !*** Initialize output variables
                !**********************************    
                BCL_P   =  missing   
                TBM     =  missing    
                BCL_H   =  missing  
                TDEF    =  missing
   
                T_ADV45   =  missing
                H_ADV45   =  missing
                SHDEF_MIX =  missing
                LHDEF_MIX =  missing  
                EADV_MIX  =  missing  


                !-- Open netcdf file
                call check( nf_open(trim(in_files(ff)), nf_nowrite, ncid_in) )

                !--------------------------
                !--  3 dimensional read  
                !-------------------------- 
                !-- Temperature profile (K)
                call check( nf_inq_varid    ( ncid_in, trim(vnames(1)), id_in) )
                call check( nf_get_var_real ( ncid_in, id_in, TK  ) )
 
                !-- Specific Humidity profile (kg/kg)
                call check( nf_inq_varid   (ncid_in, trim(vnames(2)), id_in) ) 
                call check( nf_get_var_real(ncid_in, id_in, QK  ) )
   
                !-- Geopotential Height (m)  ;  must calculate from hydrostatic assumption
                call check( nf_inq_varid   (ncid_in, trim(vnames(3)), id_in) ) 
                call check( nf_get_var_real(ncid_in, id_in, HGT  ) )
                

                !--------------------------
                !--  2 dimensional read  
                !-------------------------- 
                !-- 2-meter Temperature (K)
                call check( nf_inq_varid   (ncid_in, trim(vnames(5)), id_in) )   
                call check( nf_get_var_real(ncid_in, id_in, T2M  ) )
   
                !-- 2-meter Specific Humidity (kg/kg)
                call check( nf_inq_varid   (ncid_in, trim(vnames(6)) , id_in) ) 
                call check( nf_get_var_real(ncid_in, id_in, Q2M  ) )
   
                !-- Surface Pressure (Pa)
                call check( nf_inq_varid   (ncid_in, trim(vnames(7)), id_in) ) 
                call check( nf_get_var_real(ncid_in, id_in, PSFC  ) )
                where( PSFC.ne.missing )  PSFC  =  PSFC * 1e2
                call check( nf_close(ncid_in) )



                !**********************************    
                !*** Loop over time
                !**********************************    
                do tt = 1,ntim

                    !---------------------------------------------
                    !---
                    !--- Call HCF variables calculate routines
                    !---
                    !---------------------------------------------
                    tmpk  =  TK  (:,tt)
                    qhum  =  QK  (:,tt)
                    ht    =  HGT (:,tt)
                    pres  =  plev
                    call hcfcalc( nlev            ,  missing                               ,    &
                                  tmpk            ,  pres         , qhum        , ht       ,    &
                                  T2M      (tt)   ,  PSFC     (tt), Q2M     (tt), H2M      ,    &
                                  TBM      (tt)   ,  BCL_H    (tt), BCL_P   (tt), TDEF (tt),    &
                                  H_ADV45  (tt)   ,  tran_p       , T_ADV45 (tt)           ,    &
                                  SHDEF_MIX(tt)   ,  LHDEF_MIX(tt), EADV_MIX(tt)                )

                end do  !*** end of time loop



                !---------------------------------------------
                !---                                       ---
                !--- Begin the Output to NetCDF Section    ---
                !---                                       ---
                !---------------------------------------------   
                !*** Name the output file and open for writing
                write(*,*)" Beginning Output for file number "
                NAME1  = trim(out_dir)//"/"//trim(out_name)//"."//trim(init_date)//".HCF_VARS.nc"
                write(*,*) trim(NAME1)
                call check( nf_create(trim(NAME1), nf_64bit_offset, ncido) )
                
                !*** Define dimensions
                call check( nf_def_dim(ncido, "time", ntim, dimid1) )
                alldims = (/ dimid1 /)
      
                !*** Define variable
                call check( nf_def_var(ncido, "time"    , nf_double, 1, alldims(1), varid_time))

                call check( nf_def_var(ncido, "BCL_P"   , nf_float, ndims, alldims, varid0))
                call check( nf_def_var(ncido, "BCL_H"   , nf_float, ndims, alldims, varid2))
                call check( nf_def_var(ncido, "THETA_BM", nf_float, ndims, alldims, varid1))
                call check( nf_def_var(ncido, "TDEF"    , nf_float, ndims, alldims, varid3))
      
                call check( nf_def_var(ncido, "EADV"    , nf_float, ndims, alldims, varid4))
                call check( nf_def_var(ncido, "LHDEF"   , nf_float, ndims, alldims, varid6))
                call check( nf_def_var(ncido, "SHDEF"   , nf_float, ndims, alldims, varid5))
                call check( nf_def_var(ncido, "Tadv_45" , nf_float, ndims, alldims, varid7))
                call check( nf_def_var(ncido, "Hadv_45" , nf_float, ndims, alldims, varid8))

                !*** Assign long_name, missing values, and units
                lname  =  "Time offset from midnight"
                units  =  "seconds since 1998-1-1 0:00:00 0:00"
                call check( nf_put_att_text( ncido, varid_time, "long_name" , len_trim(lname), trim(lname) ))
                call check( nf_put_att_text( ncido, varid_time, "units"     , len_trim(units), trim(units) ))

                lname  =  "Buoyant Condensation Level Pressure"
                units  =  "Pa"
                call check( nf_put_att_text( ncido, varid0, "long_name" , len_trim(lname), trim(lname) ))
                call check( nf_put_att_text( ncido, varid0, "units"     , len_trim(units), trim(units) ))
                call check( nf_put_att_real( ncido, varid0, "_FillValue", nf_float, 1, missing ))
                lname  =  "Height at top of Buoyant Condensation Level"
                units  =  "m"
                call check( nf_put_att_text( ncido, varid2, "long_name", len_trim(lname), trim(lname) ))
                call check( nf_put_att_text( ncido, varid2, "units"    , len_trim(units), trim(units) ))
                call check( nf_put_att_real( ncido, varid2, "_FillValue", nf_float, 1, missing ))
                lname  =  "Potential temperature at the BCL"
                units  =  "K"
                call check( nf_put_att_text( ncido, varid1, "long_name", len_trim(lname), trim(lname) ))
                call check( nf_put_att_text( ncido, varid1, "units"    , len_trim(units), trim(units) ))
                call check( nf_put_att_real( ncido, varid1, "_FillValue", nf_float, 1, missing ))
                lname  =  "Potential temperature deficit needed to reach saturation"
                units  =  "K"
                call check( nf_put_att_text( ncido, varid3, "long_name", len_trim(lname), trim(lname) ))
                call check( nf_put_att_text( ncido, varid3, "units"    , len_trim(units), trim(units) ))
                call check( nf_put_att_real( ncido, varid3, "_FillValue", nf_float, 1, missing ))
      
                lname  =  "Latent Heat deficit for current Mixed Layer"
                units  =  "J/m2"
                call check( nf_put_att_text( ncido, varid6, "long_name", len_trim(lname), trim(lname) ))
                call check( nf_put_att_text( ncido, varid6, "units"    , len_trim(units), trim(units) ))
                call check( nf_put_att_real( ncido, varid6, "_FillValue", nf_float, 1, missing ))
                lname  =  "Sensible Heat deficit for current Mixed Layer"
                units  =  "J/m2"
                call check( nf_put_att_text( ncido, varid5, "long_name", len_trim(lname), trim(lname) ))
                call check( nf_put_att_text( ncido, varid5, "units"    , len_trim(units), trim(units) ))
                call check( nf_put_att_real( ncido, varid5, "_FillValue", nf_float, 1, missing ))
                lname  =  "Potential temperature at EADV=45 transition"
                units  =  "K"
                call check( nf_put_att_text( ncido, varid7, "long_name", len_trim(lname), trim(lname) ))
                call check( nf_put_att_text( ncido, varid7, "units"    , len_trim(units), trim(units) ))
                call check( nf_put_att_real( ncido, varid7, "_FillValue", nf_float, 1, missing ))
                lname  =  "Height at EADV=45 transition"
                units  =  "m"
                call check( nf_put_att_text( ncido, varid8, "long_name", len_trim(lname), trim(lname) ))
                call check( nf_put_att_text( ncido, varid8, "units"    , len_trim(units), trim(units) ))
                call check( nf_put_att_real( ncido, varid8, "_FillValue", nf_float, 1, missing ))
                lname  =  "Energy Advantage of the current mixed layer"
                units  =  "degrees"
                call check( nf_put_att_text( ncido, varid4, "long_name", len_trim(lname), trim(lname) ))
                call check( nf_put_att_text( ncido, varid4, "units"    , len_trim(units), trim(units) ))
                call check( nf_put_att_real( ncido, varid4, "_FillValue", nf_float, 1, missing ))
                call check( nf_enddef (ncido) )
           

                !*** Write data out
                call check( nf_put_var_double( ncido,  varid_time, time_stamp) )
                call check( nf_put_var_real( ncido,      varid0, BCL_P  ) )  
                call check( nf_put_var_real( ncido,      varid1, TBM    ) )
                call check( nf_put_var_real( ncido,      varid2, BCL_H  ) )
                call check( nf_put_var_real( ncido,      varid3, TDEF   ) )
                call check( nf_put_var_real( ncido,      varid4, EADV_MIX   ) )  
                call check( nf_put_var_real( ncido,      varid5, SHDEF_MIX  ) )
                call check( nf_put_var_real( ncido,      varid6, LHDEF_MIX  ) )
                call check( nf_put_var_real( ncido,      varid7, T_ADV45    ) )
                call check( nf_put_var_real( ncido,      varid8, H_ADV45    ) )
                call check( nf_close       ( ncido ) )


                !**** Deallocate variables that have a time dimension
                if(allocated(TBM       )) deallocate(TBM       )
                if(allocated(TDEF      )) deallocate(TDEF      )
                if(allocated(BCL_P     )) deallocate(BCL_P     )
                if(allocated(BCL_H     )) deallocate(BCL_H     )
                if(allocated(T_ADV45   )) deallocate( T_ADV45    )
                if(allocated(H_ADV45   )) deallocate( H_ADV45    )
                if(allocated(SHDEF_MIX )) deallocate( SHDEF_MIX  )
                if(allocated(LHDEF_MIX )) deallocate( LHDEF_MIX  )
                if(allocated(EADV_MIX  )) deallocate( EADV_MIX   )

                if(allocated(TK       ))  deallocate(TK  ) 
                if(allocated(QK       ))  deallocate(QK  ) 
                if(allocated(HGT      ))  deallocate(HGT ) 
                if(allocated(T2M      ))  deallocate(T2M ) 
                if(allocated(Q2M      ))  deallocate(Q2M ) 
                if(allocated(PSFC     ))  deallocate(PSFC)


            end if !*** 
                   !*** HCF SECTION
                   !***



            !**************************************
            !**************************************
            !**************************************
            !****                              ****
            !****        Mixing Diagram        ****
            !****                              ****
            !**************************************
            !**************************************
            !**************************************

            if( MIXING_OUT ) then

                !**********************************    
                !*** Allocate output arrays
                !**********************************    
                if(.not.allocated(SH_ENT)) allocate(SH_ENT(ndays))
                if(.not.allocated(LH_ENT)) allocate(LH_ENT(ndays))
                if(.not.allocated(SH_SFC)) allocate(SH_SFC(ndays))
                if(.not.allocated(LH_SFC)) allocate(LH_SFC(ndays))
                if(.not.allocated(SH_TOT)) allocate(SH_TOT(ndays))
                if(.not.allocated(LH_TOT)) allocate(LH_TOT(ndays))

                !**********************************    
                !*** Initialize output variables
                !**********************************    
                SH_ENT    =  missing
                LH_ENT    =  missing
                SH_SFC    =  missing
                LH_SFC    =  missing
                SH_TOT    =  missing
                LH_TOT    =  missing

                write(*,*) 
                write(*,*) 
                write(*,*) 
                write(*,*)  "  !!!!!!!!!!!      Mixing Diagram Begins      !!!!!!!!!!!!  "
                !-- allocate temporary arrays
                if(.not.allocated(T2M0 ))  allocate(T2M0 (ntim))
                if(.not.allocated(Q2M0 ))  allocate(Q2M0 (ntim))
                if(.not.allocated(PSFC0))  allocate(PSFC0(ntim))
                if(.not.allocated(SH0  ))  allocate(SH0  (ntim))
                if(.not.allocated(LH0  ))  allocate(LH0  (ntim))
                if(.not.allocated(PBLH0))  allocate(PBLH0(ntim))

                if(.not.allocated( t2m_in    ))  allocate( t2m_in    (ntim,1))
                if(.not.allocated( q2m_in    ))  allocate( q2m_in    (ntim,1))
                if(.not.allocated( psf_in    ))  allocate( psf_in    (ntim,1))
                if(.not.allocated( pbl_in    ))  allocate( pbl_in    (ntim,1))
                if(.not.allocated( shf_in    ))  allocate( shf_in    (ntim,1))
                if(.not.allocated( lhf_in    ))  allocate( lhf_in    (ntim,1))

                if(.not.allocated( sh_ent0   ))  allocate( sh_ent0   (ndays,1))
                if(.not.allocated( lh_ent0   ))  allocate( lh_ent0   (ndays,1))
                if(.not.allocated( sh_tot0   ))  allocate( sh_tot0   (ndays,1))
                if(.not.allocated( lh_tot0   ))  allocate( lh_tot0   (ndays,1))
                if(.not.allocated( sh_sfc0   ))  allocate( sh_sfc0   (ndays,1))
                if(.not.allocated( lh_sfc0   ))  allocate( lh_sfc0   (ndays,1))

                !--------------------------
                !--  2 dimensional read  
                !-------------------------- 
                !-- Open netcdf file
                call check( nf_open(trim(in_files(ff)), nf_nowrite, ncid_in) )

                !-- 2-meter Temperature (K)
                call check( nf_inq_varid(ncid_in, trim(vnames(5)), id_in) )   
                call check( nf_get_var_real(ncid_in, id_in, T2M0  ) )
   
                !-- 2-meter Specific Humidity (kg/kg)
                call check( nf_inq_varid(ncid_in, trim(vnames(6)) , id_in) ) 
                call check( nf_get_var_real(ncid_in, id_in, Q2M0  ) )
   
                !-- Surface Pressure (Pa)
                call check( nf_inq_varid(ncid_in, trim(vnames(7)), id_in) ) 
                call check( nf_get_var_real(ncid_in, id_in, PSFC0  ) )

                !-- Surface Sensible Heat Flux (W/m2)
                call check( nf_inq_varid(ncid_in, trim(vnames(8)) , id_in) ) 
                call check( nf_get_var_real(ncid_in, id_in, SH0   ) )
                !--- shift by 1 timestep because it is the 3-hour average over the past hour
                SH0(:ntim-1)   =   SH0(2:)  
                SH0(ntim)      =   missing

                !-- Surface Latent Heat Flux (W/m2)
                call check( nf_inq_varid(ncid_in, trim(vnames(9)) , id_in) ) 
                call check( nf_get_var_real(ncid_in, id_in, LH0   ) )
                !--- shift by 1 timestep because it is the 3-hour average over the past hour
                LH0(:ntim-1)   =   LH0(2:)  
                LH0(ntim)      =   missing

                !-- PBL height (m)
                call check( nf_inq_varid(ncid_in, trim(vnames(10)) , id_in) ) 
                call check( nf_get_var_real(ncid_in, id_in, PBLH0   ) )
                call check( nf_close(ncid_in) )
                !--- shift by 1 timestep because it is the 3-hour average over the past hour
                PBLH0(:ntim-1)   =   PBLH0(2:)  
                PBLH0(ntim)      =   missing


                !---------------------------------------------
                !---
                !--- Call HCF variables calculate routines
                !---
                !---------------------------------------------
                t2m_in(:,1)   =  T2M0  
                q2m_in(:,1)   =  Q2M0  
                psf_in(:,1)   =  PSFC0 
                pbl_in(:,1)   =  PBLH0 
                shf_in(:,1)   =  SH0    
                lhf_in(:,1)   =  LH0   
                call mixing_diag_daily ( 1       , ntim    , ndays  ,  hr_p_day,          &  
                                         t2m_in  , psf_in  , q2m_in ,                     &
                                         pbl_in  , shf_in  , lhf_in ,  dt,                &
                                         sh_ent0 , lh_ent0 ,                              &
                                         sh_sfc0 , lh_sfc0 , sh_tot0,  lh_tot0 ,  missing )
                SH_ENT   =    sh_ent0(:,1)
                LH_ENT   =    lh_ent0(:,1)
                SH_TOT   =    sh_tot0(:,1)
                LH_TOT   =    lh_tot0(:,1)
                SH_SFC   =    sh_sfc0(:,1)
                LH_SFC   =    lh_sfc0(:,1)


                !-- deallocate temporary arrays
                if(allocated(T2M0 ))  deallocate(T2M0 )
                if(allocated(Q2M0 ))  deallocate(Q2M0 )
                if(allocated(PSFC0))  deallocate(PSFC0)
                if(allocated(SH0  ))  deallocate(SH0  )
                if(allocated(LH0  ))  deallocate(LH0  )
                if(allocated(PBLH0))  deallocate(PBLH0)

                if(allocated( t2m_in    ))  deallocate( t2m_in )
                if(allocated( q2m_in    ))  deallocate( q2m_in )
                if(allocated( psf_in    ))  deallocate( psf_in )
                if(allocated( pbl_in    ))  deallocate( pbl_in )
                if(allocated( shf_in    ))  deallocate( shf_in )
                if(allocated( lhf_in    ))  deallocate( lhf_in )

                if(allocated( sh_ent0   ))  deallocate( sh_ent0   )
                if(allocated( lh_ent0   ))  deallocate( lh_ent0   )
                if(allocated( sh_tot0   ))  deallocate( sh_tot0   )
                if(allocated( lh_tot0   ))  deallocate( lh_tot0   )
                if(allocated( sh_sfc0   ))  deallocate( sh_sfc0   )
                if(allocated( lh_sfc0   ))  deallocate( lh_sfc0   )




                !---------------------------------------------
                !---                                       ---
                !--- Begin the Output to NetCDF Section    ---
                !---                                       ---
                !---------------------------------------------   
                !*** Name the output file and open for writing
                write(*,*)" Beginning Output for file number "
                NAME1  = trim(out_dir1)//"/"//trim(out_name)//"."//trim(init_date)//".Mixing_Diag.nc"
                write(*,*) trim(NAME1)
                call check( nf_create(trim(NAME1), nf_64bit_offset, ncido) )
                
                !*** Define dimensions
                call check( nf_def_dim(ncido, "time", ndays, dimid1) )
                alldims = (/ dimid1 /)
      
                !*** Define variable
                call check( nf_def_var(ncido, "time"   , nf_double, 1, alldims(1), varid_time))
                call check( nf_def_var(ncido, "SH_ent" , nf_float, ndims, alldims, varid2))
                call check( nf_def_var(ncido, "LH_ent" , nf_float, ndims, alldims, varid3))
                call check( nf_def_var(ncido, "SH_tot" , nf_float, ndims, alldims, varid4))
                call check( nf_def_var(ncido, "LH_tot" , nf_float, ndims, alldims, varid5))
                call check( nf_def_var(ncido, "SH_sfc" , nf_float, ndims, alldims, varid6))
                call check( nf_def_var(ncido, "LH_sfc" , nf_float, ndims, alldims, varid7))
      
                !*** Assign long_name, missing values, and units
                lname  =  "Entrainment Sensible Heat Flux"
                units  =  "W m-2"
                call check( nf_put_att_text( ncido, varid2, "long_name", len_trim(lname), trim(lname) ))
                call check( nf_put_att_text( ncido, varid2, "units"    , len_trim(units), trim(units) ))
                call check( nf_put_att_real( ncido, varid2, "_FillValue", nf_float, 1, missing ))
                lname  =  "Entrainment Latent Heat Flux"
                units  =  "W m-2"
                call check( nf_put_att_text( ncido, varid3, "long_name", len_trim(lname), trim(lname) ))
                call check( nf_put_att_text( ncido, varid3, "units"    , len_trim(units), trim(units) ))
                call check( nf_put_att_real( ncido, varid3, "_FillValue", nf_float, 1, missing ))
                lname  =  "Total Sensible Heat Flux"
                units  =  "W m-2"
                call check( nf_put_att_text( ncido, varid4, "long_name", len_trim(lname), trim(lname) ))
                call check( nf_put_att_text( ncido, varid4, "units"    , len_trim(units), trim(units) ))
                call check( nf_put_att_real( ncido, varid4, "_FillValue", nf_float, 1, missing ))
                lname  =  "Total Latent Heat Flux"
                units  =  "W m-2"
                call check( nf_put_att_text( ncido, varid5, "long_name", len_trim(lname), trim(lname) ))
                call check( nf_put_att_text( ncido, varid5, "units"    , len_trim(units), trim(units) ))
                call check( nf_put_att_real( ncido, varid5, "_FillValue", nf_float, 1, missing ))
                lname  =  "Surface Sensible Heat Flux"
                units  =  "W m-2"
                call check( nf_put_att_text( ncido, varid6, "long_name", len_trim(lname), trim(lname) ))
                call check( nf_put_att_text( ncido, varid6, "units"    , len_trim(units), trim(units) ))
                call check( nf_put_att_real( ncido, varid6, "_FillValue", nf_float, 1, missing ))
                lname  =  "Surface Latent Heat Flux"
                units  =  "W m-2"
                call check( nf_put_att_text( ncido, varid7, "long_name", len_trim(lname), trim(lname) ))
                call check( nf_put_att_text( ncido, varid7, "units"    , len_trim(units), trim(units) ))
                call check( nf_put_att_real( ncido, varid7, "_FillValue", nf_float, 1, missing ))
                call check( nf_enddef (ncido) )
           
                !*** Write data out
                call check( nf_put_var_real( ncido,      varid2, SH_ENT  ) )
                call check( nf_put_var_real( ncido,      varid3, LH_ENT  ) )
                call check( nf_put_var_real( ncido,      varid4, SH_TOT  ) )
                call check( nf_put_var_real( ncido,      varid5, LH_TOT  ) )
                call check( nf_put_var_real( ncido,      varid6, SH_SFC  ) )
                call check( nf_put_var_real( ncido,      varid7, LH_SFC  ) )
                call check( nf_close       ( ncido ) )

                !**** Deallocate variables that have a time dimension
                if(allocated(SH_ENT    )) deallocate(SH_ENT    )
                if(allocated(LH_ENT    )) deallocate(LH_ENT    )
                if(allocated(SH_SFC    )) deallocate(SH_SFC    )
                if(allocated(LH_SFC    )) deallocate(LH_SFC    )
                if(allocated(SH_TOT    )) deallocate(SH_TOT    )
                if(allocated(LH_TOT    )) deallocate(LH_TOT    )


            end if   !**** Turn on mixing diagrams calculated daily





            !**************************************
            !**************************************
            !**************************************
            !****                              ****
            !****      RH Tendency Section     ****
            !****                              ****
            !**************************************
            !**************************************
            !**************************************
            if( RH_OUT ) then 

                    write(*,*) 
                    write(*,*) 
                    write(*,*) 
                    write(*,*)  "  !!!!!!!!!!!      RH Tendency Begins      !!!!!!!!!!!!  "
                    !**********************************    
                    !*** Allocate output arrays
                    !**********************************    
                    if(.not.allocated( EF      )) allocate( EF      (ntim))
                    if(.not.allocated( NE      )) allocate( NE      (ntim))
                    if(.not.allocated( GROWTH  )) allocate( GROWTH  (ntim))
                    if(.not.allocated( HEATING )) allocate( HEATING (ntim))
                    if(.not.allocated( DRYING  )) allocate( DRYING  (ntim))
                    if(.not.allocated( dRH_dt  )) allocate( dRH_dt  (ntim))

                    !**********************************    
                    !*** Initialize output variables
                    !**********************************    
                    EF       =  missing
                    NE       =  missing
                    GROWTH   =  missing
                    HEATING  =  missing
                    DRYING   =  missing
                    dRH_dt   =  missing

                    plev  =  (/ 1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, &
                                 750, 725, 700, 675, 650, 625, 600, 575, 550, 525, &
                                 500, 475, 450, 425, 400, 375, 350, 325, 300, 275, &
                                 250, 225, 200, 175, 150, 125, 100 /)  * 1e2


                    !**********************************    
                    !*** Allocate input arrays
                    !**********************************    
                    if(.not.allocated(TK0  ))  allocate(TK0  (nlev,ntim))
                    if(.not.allocated(QK0  ))  allocate(QK0  (nlev,ntim))
                    if(.not.allocated(HGT0 ))  allocate(HGT0 (nlev,ntim))

                    if(.not.allocated(T2M0 ))  allocate(T2M0 (ntim))
                    if(.not.allocated(Q2M0 ))  allocate(Q2M0 (ntim))
                    if(.not.allocated(PSFC0))  allocate(PSFC0(ntim))
                    if(.not.allocated(SH0  ))  allocate(SH0  (ntim))
                    if(.not.allocated(LH0  ))  allocate(LH0  (ntim))
                    if(.not.allocated(PBLH0))  allocate(PBLH0(ntim))

                    if(.not.allocated( hgt2      ))  allocate( hgt2    (ntim,nlev1) )
                    if(.not.allocated( tmp2      ))  allocate( tmp2    (ntim,nlev1) )
                    if(.not.allocated( qhum2     ))  allocate( qhum2   (ntim,nlev1) )
                    if(.not.allocated( press2    ))  allocate( press2  (ntim,nlev1) )
                    if(.not.allocated( pbl2      ))  allocate( pbl2    (ntim) )
                    if(.not.allocated( shf2      ))  allocate( shf2    (ntim) )
                    if(.not.allocated( lhf2      ))  allocate( lhf2    (ntim) )
                    if(.not.allocated( ef0       ))  allocate( ef0     (ntim) )
                    if(.not.allocated( ne0       ))  allocate( ne0     (ntim) )
                    if(.not.allocated( heating0  ))  allocate( heating0(ntim) )
                    if(.not.allocated( growth0   ))  allocate( growth0 (ntim) )
                    if(.not.allocated( dry0      ))  allocate( dry0    (ntim) )
                    if(.not.allocated( drh_dt0   ))  allocate( drh_dt0 (ntim) )

                    !--------------------------
                    !--  3 dimensional read  
                    !-------------------------- 
                    !-- Open netcdf file
                    call check( nf_open(trim(in_files(ff)), nf_nowrite, ncid_in) )

                    !-- Temperature profile (K)
                    call check( nf_inq_varid    ( ncid_in, trim(vnames(1)), id_in) )
                    call check( nf_get_var_real ( ncid_in, id_in, TK0  ) )
 
                    !-- Specific Humidity profile (kg/kg)
                    call check( nf_inq_varid(ncid_in, trim(vnames(2)), id_in) ) 
                    call check( nf_get_var_real(ncid_in, id_in, QK0  ) )
   
                    !-- Geopotential Height (m)
                    call check( nf_inq_varid(ncid_in, trim(vnames(3)), id_in) ) 
                    call check( nf_get_var_real(ncid_in, id_in, HGT0  ) )

                    !--------------------------
                    !--  2 dimensional read  
                    !-------------------------- 
                    !-- 2-meter Temperature (K)
                    call check( nf_inq_varid(ncid_in, trim(vnames(5)), id_in) )   
                    call check( nf_get_var_real(ncid_in, id_in, T2M0  ) )
   
                    !-- 2-meter Specific Humidity (kg/kg)
                    call check( nf_inq_varid(ncid_in, trim(vnames(6)) , id_in) ) 
                    call check( nf_get_var_real(ncid_in, id_in, Q2M0  ) )
   
                    !-- Surface Pressure (Pa)
                    call check( nf_inq_varid(ncid_in, trim(vnames(7)), id_in) ) 
                    call check( nf_get_var_real(ncid_in, id_in, PSFC0  ) )
                    where( PSFC0.ne.missing )  PSFC0  =  PSFC0 * 1e2
                    do zz=1,nlev
                       press2(:,zz+1) =  plev(zz)
                    end do

                    !-- Surface Sensible Heat Flux (W/m2)
                    call check( nf_inq_varid(ncid_in, trim(vnames(8)) , id_in) ) 
                    call check( nf_get_var_real(ncid_in, id_in, SH0   ) )

                    !-- Surface Latent Heat Flux (W/m2)
                    call check( nf_inq_varid(ncid_in, trim(vnames(9)) , id_in) ) 
                    call check( nf_get_var_real(ncid_in, id_in, LH0   ) )

                    !-- PBL height (m)
                    call check( nf_inq_varid(ncid_in, trim(vnames(10)) , id_in) ) 
                    call check( nf_get_var_real(ncid_in, id_in, PBLH0   ) )
                    call check( nf_close(ncid_in) )





                    !---------------------------------------------
                    !---
                    !--- Call RH-Tendency routine
                    !---
                    !---------------------------------------------
                    ef0       =  missing
                    ne0       =  missing
                    heating0  =  missing
                    growth0   =  missing
                    dry0      =  missing
                    drh_dt0   =  missing

                    tmp2    (:,1)  =  T2M0 
                    qhum2   (:,1)  =  Q2M0 
                    hgt2    (:,1)  =  2.
                    press2  (:,1)  =  PSFC0

                    tmp2    (:,2:) =  transpose(TK0 ) 
                    qhum2   (:,2:) =  transpose(QK0 )
                    hgt2    (:,2:) =  transpose(HGT0)

                    pbl2           =  PBLH0 
                    shf2           =  SH0    
                    lhf2           =  LH0   

                    write(*,*) nlev1, ntim, "    LOOKY "
                    call rh_tend_calc ( nlev1    ,  ntim   ,                         &
                                        tmp2     ,  qhum2  ,  hgt2   ,  press2 ,     &
                                        pbl2     ,  shf2   ,  lhf2   ,  dt     ,     &
                                        ef0      ,  ne0    ,                         &
                                        heating0 ,  growth0,  dry0   ,  drh_dt0, missing  )
    
                    EF       =  ef0      
                    NE       =  ne0      
                    GROWTH   =  growth0 
                    HEATING  =  heating0  
                    DRYING   =  dry0     
                    dRH_dt   =  drh_dt0



                    !-- deallocate temporary arrays
                    if(allocated(TK0  ))  deallocate(TK0  )
                    if(allocated(QK0  ))  deallocate(QK0  )
                    if(allocated(HGT0 ))  deallocate(HGT0 )
                    if(allocated(T2M0 ))  deallocate(T2M0 )
                    if(allocated(Q2M0 ))  deallocate(Q2M0 )
                    if(allocated(PSFC0))  deallocate(PSFC0)
                    if(allocated(SH0  ))  deallocate(SH0  )
                    if(allocated(LH0  ))  deallocate(LH0  )
                    if(allocated(PBLH0))  deallocate(PBLH0)

                    if(allocated( hgt2    ))  deallocate( hgt2   )
                    if(allocated( tmp2    ))  deallocate( tmp2   )
                    if(allocated( qhum2   ))  deallocate( qhum2  )
                    if(allocated( press2  ))  deallocate( press2 )
                    if(allocated( pbl2    ))  deallocate( pbl2   )
                    if(allocated( shf2    ))  deallocate( shf2   )
                    if(allocated( lhf2    ))  deallocate( lhf2   )

                    if(allocated( ef0       ))  deallocate( ef0      )
                    if(allocated( ne0       ))  deallocate( ne0      )
                    if(allocated( heating0  ))  deallocate( heating0 )
                    if(allocated( growth0   ))  deallocate( growth0  )
                    if(allocated( dry0      ))  deallocate( dry0     )
                    if(allocated( drh_dt0   ))  deallocate( drh_dt0  )

                    !---------------------------------------------
                    !---                                       ---
                    !--- Begin the Output to NetCDF Section    ---
                    !---                                       ---
                    !---------------------------------------------   
               !*** Name the output file and open for writing
               write(*,*)" Beginning Output for file number "
               NAME1  = trim(out_dir2)//"/"//trim(out_name)//"."//trim(init_date)//".RH_TEND.nc"
               write(*,*) trim(NAME1)
               call check( nf_create(trim(NAME1), nf_64bit_offset, ncido) )
               
               !*** Define dimensions
               call check( nf_def_dim(ncido, "time", ntim, dimid3) )
               alldims = (/ dimid1 /)
      
               !*** Define variable
               call check( nf_def_var(ncido, "time"    , nf_double, 1, alldims(1), varid_time))

               call check( nf_def_var(ncido, "EF"      , nf_float, ndims, alldims, varid0))
               call check( nf_def_var(ncido, "NE"      , nf_float, ndims, alldims, varid1))
               call check( nf_def_var(ncido, "DRYING"  , nf_float, ndims, alldims, varid2))
               call check( nf_def_var(ncido, "GROWTH"  , nf_float, ndims, alldims, varid3))
               call check( nf_def_var(ncido, "HEATING" , nf_float, ndims, alldims, varid4))
               call check( nf_def_var(ncido, "dRH_dt"  , nf_float, ndims, alldims, varid5))
      
               !*** Assign long_name, missing values, and units
               lname  =  "Time offset from midnight"
               units  =  "seconds since 1998-1-1 0:00:00 0:00"
               call check( nf_put_att_text( ncido, varid_time, "long_name" , len_trim(lname), trim(lname) ))
               call check( nf_put_att_text( ncido, varid_time, "units"     , len_trim(units), trim(units) ))

               lname  =  "Evaporative Fraction"
               units  =  "unitless"
               call check( nf_put_att_text( ncido, varid0, "long_name" , len_trim(lname), trim(lname) ))
               call check( nf_put_att_text( ncido, varid0, "units"     , len_trim(units), trim(units) ))
               call check( nf_put_att_real( ncido, varid0, "_FillValue", nf_float, 1, missing ))
               lname  =  "Non-evaporative term"
               units  =  "unitless"
               call check( nf_put_att_text( ncido, varid1, "long_name", len_trim(lname), trim(lname) ))
               call check( nf_put_att_text( ncido, varid1, "units"    , len_trim(units), trim(units) ))
               call check( nf_put_att_real( ncido, varid1, "_FillValue", nf_float, 1, missing ))
               lname  =  "Change in RH at top of PBL due to drying air entrainment"
               units  =  "unitless"
               call check( nf_put_att_text( ncido, varid2, "long_name", len_trim(lname), trim(lname) ))
               call check( nf_put_att_text( ncido, varid2, "units"    , len_trim(units), trim(units) ))
               call check( nf_put_att_real( ncido, varid2, "_FillValue", nf_float, 1, missing ))
               lname  =  "Change in RH at top of PBL due to growth of PBL height"
               units  =  "unitless"
               call check( nf_put_att_text( ncido, varid3, "long_name", len_trim(lname), trim(lname) ))
               call check( nf_put_att_text( ncido, varid3, "units"    , len_trim(units), trim(units) ))
               call check( nf_put_att_real( ncido, varid3, "_FillValue", nf_float, 1, missing ))
               lname  =  "Change in RH at top of PBL due to heat entrainment"
               units  =  "unitless"
               call check( nf_put_att_text( ncido, varid4, "long_name", len_trim(lname), trim(lname) ))
               call check( nf_put_att_text( ncido, varid4, "units"    , len_trim(units), trim(units) ))
               call check( nf_put_att_real( ncido, varid4, "_FillValue", nf_float, 1, missing ))
               lname  =  "Estimated Change in Relative Humidity"
               units  =  "unitless"
               call check( nf_put_att_text( ncido, varid5, "long_name", len_trim(lname), trim(lname) ))
               call check( nf_put_att_text( ncido, varid5, "units"    , len_trim(units), trim(units) ))
               call check( nf_put_att_real( ncido, varid5, "_FillValue", nf_float, 1, missing ))
               call check( nf_enddef (ncido) )
           
               !*** Write data out
               call check( nf_put_var_double( ncido,  varid_time, time_stamp) )

               call check( nf_put_var_real( ncido,      varid0, EF      ) )  
               call check( nf_put_var_real( ncido,      varid1, NE      ) )
               call check( nf_put_var_real( ncido,      varid2, DRYING  ) )
               call check( nf_put_var_real( ncido,      varid3, GROWTH  ) )
               call check( nf_put_var_real( ncido,      varid4, HEATING ) )
               call check( nf_put_var_real( ncido,      varid5, dRH_dt  ) )
               call check( nf_close       ( ncido ) )

               !**** Deallocate variables that have a time dimension
               if(allocated(time_stamp)) deallocate(time_stamp)
               if(allocated(dRH_dt    )) deallocate(dRH_dt    )
               if(allocated(EF        )) deallocate(EF        )
               if(allocated(NE        )) deallocate(NE        )
               if(allocated(DRYING    )) deallocate(DRYING    )
               if(allocated(GROWTH    )) deallocate(GROWTH    )
               if(allocated(HEATING   )) deallocate(HEATING   )



            end if !*** 
                   !*** RH Tendency SECTION
                   !***



            if(allocated(time_stamp)) deallocate(time_stamp)
            write(*,*) " "
            write(*,*) " Finished  ---  On to the Next  ---"
            write(*,*) " "
            write(*,*) " "
            write(*,*) " "
            write(*,*) " "
            write(*,*) " "
            write(*,*) " "
         end do  !*** end of file loop
         write(*,*) " "
         write(*,*) " "
         write(*,*) " !!!!!! "
         write(*,*) " "
         write(*,*) " END!!! "
         write(*,*) " "
         write(*,*) " !!!!!! "
         write(*,*) " "
         write(*,*) " "
!-------------------------- End Interface ----------------------------------------

end Program Coupling_metrics






!------------------------------------------------------------------
     subroutine check(status)
        integer, intent ( in) :: status

        if(status /= 0) then
          ! print *, trim(nf_strerror(status))
          !print *, nf_strerror(status)
          write(*,*) " ERROR:   ",status
          stop "Stopped"
        end if
     end subroutine check
!------------------------------------------------------------------

