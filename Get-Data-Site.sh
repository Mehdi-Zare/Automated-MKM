#!/bin/sh  

# To Get the SCF-Data of clean surface (just energy from OSZICAR)

       grep -a 'E0= ' OSZICAR | tail -1 > energy-line;
       awk '{print $5}' energy-line > SCF-Data;
       rm energy-line
