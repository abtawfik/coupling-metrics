Contains several fortran 90 subroutines used to calculate various land-atmosphere coupling metrics.
For more information about each metric and use of the modules go to www.coupling-metrics.com

General use:
The modules available are meant to facilitate use of the most common land-atmosphere coupling metrics for research purposes. The modules are currently written in Fortran 90 and are individually packaged.  The most common use is to take a given module and create a library that is usable in the desired programming language.  For example if you are using the NCAR Common Langauge (NCL) then you could perform the following actions with the files in the */ncl_portable/ folder
> WRAPIT ctp_hilow.stub ctp_hilow.f90 <
This will create a file called ctp_hilow.so.  This shared object file can then be called from within your NCL script by adding a line like this:
> external CTP "ctp_hilow.so"
> CTP::ctp_hi_low( nlev, t, q, p, t2m , q2m , psfc, CTP, HILOW , missing ) <
Note that you have 
Below is a brief description of each metric and use of the modules in general:
1) Convective Triggering Potential (aka CTP-HIlow or just CTP)
Description - evaluates morning atmospheric profiles to determine whether dry or wet soils are more likely to trigger convection

