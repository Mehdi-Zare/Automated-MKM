#!/bin/sh  

# To get the energy and freq for adsorbate

       grep -a 'E0= ' OSZICAR | tail -1 > energy-line;
       awk '{print $5}' energy-line > SCF-Data;
       awk '{print $1}' freq.dat >> SCF-Data;
       rm energy-line
