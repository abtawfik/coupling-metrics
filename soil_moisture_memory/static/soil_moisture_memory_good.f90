!---------------------------------------------------------------------------------
!
! Description:
! Calculates the lagged auto-correlation from daily soil moisture data
! This returns the soil moisture memory.  Different thresholds can be used for
! defining the "memory".  
!
! Pearson Correlation
! EQUATION:  covariance(Initial Soil Moisture, Subsequent Day Soil Moisture)
!            ---------------------------------------------------------------
!                      stddev(Initial SM)  *  stddev(Subsequent SM)
!
! Memory is said to be "lost" when the correlation drops below the the e-folding 
! time.  Other thresholds (e.g. correlation cutoffs) can be used if desired.  Just 
! change the "threshold" parameter variable below.
!
! Reference:
!
! Author and Revision History:
! A.B. Tawfik on May 2015
!
!---------------------------------------------------------------------------------
module Soil_Memory_Mod

     !
     ! subroutine name 
     !
     public soilm_memory

!---------------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------------



!---------------------------------------------------------------------------------
!
! subroutines:  calculates lagged autocorrelation of soil moisture.  
!               Assumes left most dimension is time and
!               the right most dimension is typically space (like lat/lon or site)
!               but it can also be depth OR both combined into a single dimension
!               Output is dimensioned:  
!               (days, spatial dimension like lat/lon or depth)
!---------------------------------------------------------------------------------
  subroutine soilm_memory( dim2,  ntim,  soilm,  smemory,  missing )
   implicit none

!
! Input/Output Variables
!
   integer, intent(in   )                         ::  dim2           ! *** missing value - useful for obs
   integer, intent(in   )                         ::  ntim           ! *** number of days

   real(4), intent(in   )                         ::  missing        ! *** missing value - useful for obs
   real(4), intent(in   ), dimension(ntim,dim2)   ::  soilm          ! *** soil moisture - (days,space)
   real(4), intent(out  ), dimension(     dim2)   ::  smemory        ! *** soil moisture memory [days]
!
! Local variables
!
   real(4), parameter   ::  threshold    = 1./exp(1.)   !*** threshold used to define where memory is lost
                                                        !*** can be changed depending on prefernce
   real(4), parameter   ::  sample_limit = 20.          !*** # of valid days needed to calculate memory

   integer                           ::  tt, n1
   real(4), dimension(     dim2)     ::  nsample
   real(4), dimension(     dim2)     ::  nsample_lagged
   real(4), dimension(ntim,dim2)     ::  soilm_mean                       
   real(4), dimension(ntim,dim2)     ::  soilm_diff     
   real(4), dimension(     dim2)     ::  soilm_std      
   real(4), dimension(ntim,dim2)     ::  soilm_lagged   
   real(4), dimension(ntim,dim2)     ::  lagged_mean    
   real(4), dimension(ntim,dim2)     ::  lagged_diff    
   real(4), dimension(     dim2)     ::  covar          
   real(4), dimension(     dim2)     ::  lagged_std     
   real(4), dimension(ntim,dim2)     ::  correlation
   real(4), dimension(ntim,dim2)     ::  correlation_diff
   real(4), dimension(ntim,dim2)     ::  time_ind

!-----------------------------------------------------------------------------


      !--------------------------------------------
      !--- Initialization and preliminary calculations
      !--------------------------------------------
      !----------------------------------------------------------------------------------
      !-- initialize output arrays
      !----------------------------------------------------------------------------------
      smemory    =  missing


      !-----------------------------------------------------------------------------
      !-- Check input array to make sure all are NOT missing
      !-----------------------------------------------------------------------------
      if( all(soilm.eq.missing) )  then
          return
      end if

      !----------------------------------------------------------------------------------
      !-- initialize working arrays
      !----------------------------------------------------------------------------------
      soilm_mean =  missing
      soilm_diff =  missing
      soilm_std  =  missing

    
      !----------------------------------------------------------------------------------
      !-- initialize correlation to 1 so values never reach e-folding timescale
      !----------------------------------------------------------------------------------
      correlation  =  1.
      

      !------------------------------------------------------------
      !-- A time index used to find when e-folding time is reached
      !------------------------------------------------------------
      do tt=1,ntim
         time_ind(tt,:)  =  tt * 1.
      end do

      !----------------------------------------------------------------------------------
      !-- Sample covariance equation without dividing through by sample size:
      !--               _         _
      !--     sum( [X - X] * [Y - Y] )
      !--     
      !--  
      !-- Sample standard deviation equation without dividing through by sample size:
      !--               _        
      !--     sum( [X - X]^2 )
      !--    
      !--  
      !-- Here X = the Initial Soil moisture and Y = Subsequent day Soil moisture
      !--      No need to divid by sample size because it will be eliminated in the 
      !--      in the correlation calculation anyway
      !----------------------------------------------------------------------------------
      !----------------------------------------------------------------------------------
      !--                                                       _          _
      !-- Get mean soil moisture and difference correspoding to X and [X - X]
      !-- expand to 2-dimensional to simplify subtraction and product 
      !--                                                                          
      !----------------------------------------------------------------------------------
      nsample  =  count(soilm.ne.missing, DIM=1) * 1.
      where( nsample.gt.0 )
          soilm_mean(1,:)  =  sum( soilm, DIM=1, MASK=soilm.ne.missing )  /  nsample 
      endwhere
      do n1=2,ntim
         soilm_mean(n1,:)  =  soilm_mean(1,:)
      end do
      where( soilm.ne.missing .and. soilm_mean.ne.missing )
         soilm_diff  =  soilm - soilm_mean
      endwhere

      !------------------------------------------------------------------
      !-- Calculate sample standard deviation for initial soil moisture
      !------------------------------------------------------------------
      soilm_std  =  sqrt( sum( soilm_diff**2, DIM=1, MASK=soilm_diff.ne.missing) )

      !----------------------------------------------------------------------------------
      !--                                                   
      !-- Loop over time to populate a tridiagonal matrix of covariances for all the 
      !-- subsequent (t>1) times
      !--                                                                          
      !----------------------------------------------------------------------------------
      do tt=2,ntim
         
         !------------------------------
         !-- initialize working arrays
         !------------------------------
         soilm_lagged    =  missing
         lagged_mean     =  missing
         lagged_diff     =  missing
         nsample_lagged  =  missing
         covar           =  missing
         lagged_std      =  missing

         !------------------------------
         !-- shift by one time step
         !------------------------------
         soilm_lagged(:ntim-tt+1,:)  =  soilm(tt:,:)

         !------------------------------
         !-- Get mean of lagged soil moisture
         !------------------------------
         nsample_lagged   =  count(soilm_lagged.ne.missing, DIM=1) * 1.
         where( nsample_lagged.gt.0 )
            lagged_mean(1,:)  =  sum( soilm_lagged, DIM=1, MASK=soilm_lagged.ne.missing )  /  nsample_lagged
         endwhere
         do n1=2,ntim 
            lagged_mean(n1,:)  =  lagged_mean(1,:)
         end do

         !------------------------------
         !-- Calculate difference
         !------------------------------
         where( soilm_lagged.ne.missing .and. lagged_mean.ne.missing )
            lagged_diff  =  soilm_lagged - lagged_mean
         endwhere
       
         !------------------------------------------------------------------
         !-- Calculate sample standard deviation for lagged soil moisture
         !------------------------------------------------------------------
         lagged_std  =  sqrt(  sum( lagged_diff**2, DIM=1, MASK=lagged_diff.ne.missing) )

         !------------------------------
         !-- Calculate covariance 
         !------------------------------
         covar  =   sum( soilm_diff * lagged_diff, DIM=1, MASK=lagged_diff.ne.missing  .and.  soilm_diff.ne.missing )

         !---------------------------------
         !-- Calculate Pearson Correlation
         !---------------------------------
         where( lagged_std.ne.missing  .and. soilm_std.ne.missing  .and.  covar.ne.missing  .and.  lagged_std*soilm_std.ne.0 )
            correlation(tt,:)  =  covar / (lagged_std * soilm_std)
         endwhere

do n1=1,1
   write(*,*) "------     ",n1,tt,correlation(tt,n1), covar(n1), lagged_std(n1), soilm_std(n1)
end do

         !----------------------------------------------------------------------------------------------------
         !-- if all of the correlations for each year are below the threshold then no need to iterate further
         !-- therefore exit for efficieny 
         !----------------------------------------------------------------------------------------------------
         if( all(correlation(tt,:)-threshold.le.0) ) exit


      end do !===  End of time loop      

      !-----------------------------------------------------------------------
      !-- Find when correlation falls below the e-folding time
      !-- This is defined as the Soil Moisture Memory (smemory)
      !-----------------------------------------------------------------------
      where( correlation.ne.missing ) 
         correlation_diff  =  correlation - threshold
      endwhere 
      smemory  =  minloc( time_ind, DIM=1, MASK=correlation_diff.le.0 )

stop
      return


end subroutine soilm_memory



end module Soil_Memory_Mod
