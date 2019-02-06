#!/bin/sh  

# To create Part-Func directory, copy CONTCAR and freq.dat in that and perform the code to get partition functions at different Temperature
# and store them in "lnq" file

mkdir Part-Func
cp freq.dat CONTCAR Part-Func
cd Part-Func/


       grep -a 'E0= ' OSZICAR | tail -1 > energy-line;
       awk '{print $5}' energy-line > SCF-Data;
       awk '{print $1}' freq.dat >> SCF-Data;
       rm energy-line
