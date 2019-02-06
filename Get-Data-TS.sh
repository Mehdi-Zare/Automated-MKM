#!/bin/sh  

# To get the energy and freqs of TS and remove imaginary frequency

       grep -a 'E0= ' OSZICAR | tail -1 > energy-line;
       awk '{print $5}' energy-line > SCF-Data;
       awk '{print $1}' freq.dat >> SCF-Data;
       sed -i '2d' SCF-Data
       rm energy-line
