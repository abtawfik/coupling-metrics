#!/bin/bash


### Pick which technique
rm *.o *.mod
type="static"

### Compile
ifort -i4 -fPIC -O3 -m64 -ftz -fp-model precise -march=corei7 -axAVX -I/glade/apps/opt/netcdf/4.3.0/intel/12.1.5//include -L/glade/apps/opt/netcdf/4.3.0/intel/12.1.5/lib -lnetcdf -lnetcdff -c ../../heated_condensation/"$type"/hcfcalc.f90

ifort -i4 -fPIC -O3 -m64 -ftz -fp-model precise -march=corei7 -axAVX -I/glade/apps/opt/netcdf/4.3.0/intel/12.1.5//include -L/glade/apps/opt/netcdf/4.3.0/intel/12.1.5/lib -lnetcdf -lnetcdff -c  ../../rh_tendency/"$type"/rh_tendency_with_mixing.f90

ifort -i4 -fPIC -O3 -m64 -ftz -fp-model precise -march=corei7 -axAVX -I/glade/apps/opt/netcdf/4.3.0/intel/12.1.5//include -L/glade/apps/opt/netcdf/4.3.0/intel/12.1.5/lib -lnetcdf -lnetcdff -c ../../mixing_diagram/"$type"/mixing_diagram_daily.f90

ifort -i4 -fPIC -O3 -m64 -ftz -fp-model precise -march=corei7 -axAVX -I/glade/apps/opt/netcdf/4.3.0/intel/12.1.5//include -L/glade/apps/opt/netcdf/4.3.0/intel/12.1.5/lib -lnetcdf -lnetcdff MAIN_HCF_INTERFACE.f90 *.o -o ./RUN_NARR.exe
