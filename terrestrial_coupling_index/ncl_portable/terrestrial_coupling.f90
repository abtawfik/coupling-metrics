!---------------------------------------------------------------------------------
!
! Description:
! Calculates the terrestrial coupling parameter which is simply the slope of 
! linear-fit between soil moisture and some surface fluxes (sensible or latent)
! multiplied by the sample standard deviation of soil moisture.  Physically
! this parameter quantifies to what extend variations in soil moisture
! correspond to change variations in surface fluxes.  If the coupling parameter
! is high then soil moisture is said to vary sufficiently to influence 
! variations in surface fluxes.
!
! Terrestrial Coupling Parameter [TCP]
! EQUATION:   std dev(soil moisture) * slope(soil moisture, sfc flux)
!
! Units returned correspond to the unit of the surface flux variable
!
! Reference:  
! Dirmeyer, The terrestrial segment of soil moisture–climate coupling (2011)
!
! Author and Revision History:
! A.B. Tawfik on May 2015
!
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
!
! subroutines:  calculates terrestrial coupling parameter (TCP)
!               Assumes left most dimension is time and
!               the right most dimension is typically space (like lat/lon or site)
!               but it can also be depth OR both combined into a single dimension
!               Therefore a TCP value will be returned across the other dimension
!               Output is dimensioned:  
!               (spatial dimension like lat/lon or depth)
!---------------------------------------------------------------------------------
  subroutine terra_coupling( dim2,  ntim,  soilm,  flux,  tcp,  missing )
   implicit none

!
! Input/Output Variables
!
   integer, intent(in   )                         ::  dim2           ! *** missing value - useful for obs
   integer, intent(in   )                         ::  ntim           ! *** size of time dimension

   real(4), intent(in   )                         ::  missing        ! *** missing value - useful for obs
   real(4), intent(in   ), dimension(ntim,dim2)   ::  soilm          ! *** soil moisture - (time,space)
   real(4), intent(in   ), dimension(ntim,dim2)   ::  flux           ! *** surface flux  - (time,space)
   real(4), intent(out  ), dimension(     dim2)   ::  tcp            ! *** terrestrial coupling parameter [flux units]
!
! Local variables
!
   integer, parameter                ::  sample_limit  =  100

   integer                           ::  n1
   integer, dimension(     dim2)     ::  nsample
   real(4), dimension(ntim,dim2)     ::  soilm_mean                       
   real(4), dimension(ntim,dim2)     ::  soilm_diff     
   real(4), dimension(     dim2)     ::  soilm_std      
   real(4), dimension(     dim2)     ::  flux_std
   real(4), dimension(     dim2)     ::  flux_std_with_n
   real(4), dimension(ntim,dim2)     ::  flux_mean                       
   real(4), dimension(ntim,dim2)     ::  flux_diff     
   real(4), dimension(     dim2)     ::  covar          

!-----------------------------------------------------------------------------


      !--------------------------------------------
      !--- Initialization and preliminary calculations
      !--------------------------------------------
      !----------------------------------------------------------------------------------
      !-- initialize output arrays
      !----------------------------------------------------------------------------------
      tcp   =  missing


      !-----------------------------------------------------------------------------
      !-- Check input array to make sure all are NOT missing
      !-----------------------------------------------------------------------------
      if( all(soilm.eq.missing) .or. all(flux.eq.missing) )  then
          return
      end if

      !-----------------------------------------------------------------------------
      !-- Check sample size is at least 100 times of data
      !-- !!!!! CURRENTLY THIS IS ARBITRARY ---  NEED BETTER ROBUSTNESS TEST  !!!!!!
      !-----------------------------------------------------------------------------
      nsample  =  count(soilm.ne.missing .and. flux.ne.missing, DIM=1)
      if( all( nsample.lt.sample_limit ) )  then
          write(*,*)  " ----  Sample size is too small for all points:  less than 20 non-missing values "
          return
      end if


      !----------------------------------------------------------------------------------
      !-- initialize working arrays
      !----------------------------------------------------------------------------------
      soilm_mean       =  missing
      soilm_diff       =  missing
      soilm_std        =  missing

      flux_std         =  missing
      flux_mean        =  missing
      flux_diff        =  missing
      covar            =  missing
      flux_std_with_n  =  missing    
      

      !----------------------------------------------------------------------------------
      !-- Slope equation is:
      !--                            
      !--                covar(X, Y)                 covar(X, Y)
      !--     slope  =  ------------  ===>  TCP  =  -------------
      !--                  var(X)                       std(X)
      !--     
      !--  
      !-- Sample covariance equation without dividing through by sample size:
      !--                               _         _
      !--     covar(X, Y)  =  sum( [X - X] * [Y - Y] )
      !--     
      !--  
      !-- Sample standard deviation equation without dividing through by sample size:
      !--               _        
      !--     sum( [X - X]^2 )
      !--    
      !--  
      !-- Here X = Soil moisture and Y = Surface Fluxes
      !----------------------------------------------------------------------------------

      !----------------------------------------------------------------------------------
      !--                                                       _          _
      !-- Get mean soil moisture and difference correspoding to X and [X - X]
      !-- expand to 2-dimensional to simplify subtraction and product 
      !--                                                                          
      !----------------------------------------------------------------------------------
      nsample  =  count(soilm.ne.missing, DIM=1)
      where( nsample.gt.0 )
          soilm_mean(1,:)  =  sum( soilm, DIM=1, MASK=soilm.ne.missing )  /  real(nsample)
      endwhere
      do n1=2,ntim
         soilm_mean(n1,:)  =  soilm_mean(1,:)
      end do
      where( soilm.ne.missing .and. soilm_mean.ne.missing )
         soilm_diff  =  soilm - soilm_mean
      endwhere


      !------------------------------------------------------------------
      !-- Calculate sample standard deviation soil moisture
      !------------------------------------------------------------------
      soilm_std  =  sqrt( sum( soilm_diff**2, DIM=1, MASK=soilm_diff.ne.missing) )


      !----------------------------------------------------------------------------------
      !--                                                        _          _
      !-- Get mean surface fluxes and difference correspoding to Y and [Y - Y]
      !-- expand to 2-dimensional to simplify subtraction and product 
      !--                                                                          
      !----------------------------------------------------------------------------------
      nsample  =  count(flux.ne.missing, DIM=1)
      where( nsample.gt.0 )
          flux_mean(1,:)  =  sum( flux, DIM=1, MASK=flux.ne.missing)  /  real(nsample)
      endwhere
      do n1=2,ntim
         flux_mean(n1,:)  =  flux_mean(1,:)
      end do
      where( flux.ne.missing .and. flux_mean.ne.missing )
         flux_diff  =  flux - flux_mean
      endwhere

      !------------------------------------------------------------------
      !-- Calculate sample standard deviation flux
      !------------------------------------------------------------------
      flux_std         =  sqrt( sum( flux_diff**2, DIM=1, MASK=flux_diff.ne.missing) )
      flux_std_with_n  =  sqrt( (1./real(nsample-1)) * sum( flux_diff**2, DIM=1, MASK=flux_diff.ne.missing) )

      !------------------------------
      !-- Calculate covariance 
      !------------------------------
      covar  =   sum( soilm_diff * flux_diff, DIM=1, MASK=flux_diff.ne.missing  .and.  soilm_diff.ne.missing )

      !---------------------------------------------
      !-- Calculate Terrestrial Coupling Parameter 
      !---------------------------------------------
      where( soilm_std.ne.missing  .and.  covar.ne.missing  .and.  soilm_std*flux_std.ne.0  .and. flux_std.ne.missing)
            tcp  =  flux_std_with_n * (covar / (soilm_std * flux_std))
      endwhere


      return


end subroutine terra_coupling
