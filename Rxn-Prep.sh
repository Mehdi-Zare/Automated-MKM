#!/bin/sh

# to cp scf-date of each element to reaction direcotrt and change names to Reactant-1 Ractant-2 TS-1 TS-2 Product-1 Procut-2
mkdir REACTIONS

line=`wc -l RXNs | awk '{print $1}'`
numRxn=$(($line-1))

i=2
while  [ "$i" -le "$line" ]; do
  dirname=`head -n$i  RXNs | tail -1 | awk '{print $1}'`
  R1=`head -n$i  RXNs | tail -1 | awk '{print $2}'`;R1t=Reactant-1
  R2=`head -n$i  RXNs | tail -1 | awk '{print $3}'`;R2t=Reactant-2
  TS1=`head -n$i  RXNs | tail -1 | awk '{print $4}'`;P1t=Product-1
  TS2=`head -n$i  RXNs | tail -1 | awk '{print $5}'`;P2t=Product-2
  P1=`head -n$i  RXNs | tail -1 | awk '{print $6}'`;TS1t=TS-1
  P2=`head -n$i  RXNs | tail -1 | awk '{print $7}'`;TS2t=TS-2
  mkdir REACTIONS/$dirname
  
  if [[ $R1 != none  ]]; then
  cp SCF-DATA/$R1/SCF-Data REACTIONS/$dirname/$R1t
  fi
  if [[ $R2 != none  ]]; then
  cp SCF-DATA/$R2/SCF-Data REACTIONS/$dirname/$R2t
  fi
  if [[ $P1 != none  ]]; then
  cp SCF-DATA/$P1/SCF-Data REACTIONS/$dirname/$P1t
  fi
  if [[ $P2 != none  ]]; then
  cp SCF-DATA/$P2/SCF-Data REACTIONS/$dirname/$P2t
  fi
  if [[ $TS1 != none  ]]; then
  cp SCF-DATA/$TS1/SCF-Data REACTIONS/$dirname/$TS1t
  fi
  if [[ $TS2 != none  ]]; then
  cp SCF-DATA/$TS2/SCF-Data REACTIONS/$dirname/$TS2t
  fi
i=$(($i+1))
done

i=2
while  [ "$i" -le "$line" ]; do
  dirname=`head -n$i  RXNs | tail -1 | awk '{print $1}'`;
  if [[ $dirname == *Gas  ]] || [[ $dirname == *H2  ]];then
  R1=`head -n$i  RXNs | tail -1 | awk '{print $2}'`
  cp SCF-DATA/$R1/Part-Func/lnq REACTIONS/$dirname/
  cp SCF-DATA/$R1/Gas-Data REACTIONS/$dirname/
  fi
  i=$(($i+1))
done

cd REACTIONS

for rxn in *

  do cd $rxn
  cp ~/MKM-Codes/Temp-File .
  if [[ $rxn == *Surf  ]]; then
     ~/MKM-Codes/Surface-Rxn
  elif [[ $rxn == *Gas  ]]; then
     ~/MKM-Codes/Adsorp-Rxn
  elif [[ $rxn == *H2  ]]; then
    ~/MKM-Codes/AdsH2-Rxn
  fi

cd ../  
done


