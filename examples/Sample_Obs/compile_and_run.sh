#!/bin/bash

### Clean and select directory
rm *.o *.mod
type="static"

### Compile -- gfortran
gfortran -O3 -ffree-line-length-none -c ../../heated_condensation/"$type"/hcfcalc.f90
gfortran -O3 -ffree-line-length-none -c ../../rh_tendency/"$type"/rh_tendency.f90
gfortran -O3 -ffree-line-length-none -c ../../mixing_diagram/"$type"/mixing_diagram_daily.f90
gfortran -O3 -ffree-line-length-none Program_run_sample.f90 *.o -o ./sample_data_test.exe

### Compile
./sample_data_test.exe
