#!/bin/bash

# to prepare all properties for print

  mkdir MKM-DATA
  LT=`head -1 ~/MKM-Codes/Temp-File`
  HT=`head -2 ~/MKM-Codes/Temp-File | tail -1`
  IT=`head -3 ~/MKM-Codes/Temp-File | tail -1`
  nTemp=$(((($HT-$LT)/$IT)+1))

ii=1
Temp=$LT
linerxn=13
lineKin=10
linerxnG=12
while [ "$ii" -le "$nTemp" ]; do
  mkdir MKM-DATA/$Temp
  cd REACTIONS 
  
   
    for i in *Surf
      do cd $i
      value=`head -n$linerxn Rxn-Energy | tail -1 | awk '{print $5}'`; echo $i$'\t'$value >> ../../MKM-DATA/$Temp/Grxn
      value=`head -n$linerxn Rxn-Energy | tail -1 | awk '{print $6}'`; echo $i$'\t'$value >> ../../MKM-DATA/$Temp/Gact
      value=`head -n$lineKin Kinetics | tail -1 | awk '{print $2}'`; echo $i$'\t'$value >> ../../MKM-DATA/$Temp/Keq
      value=`head -n$lineKin Kinetics | tail -1 | awk '{print $3}'`; echo $i$'\t'$value >> ../../MKM-DATA/$Temp/kf
      value=`head -n$lineKin Kinetics | tail -1 | awk '{print $4}'`; echo $i$'\t'$value >> ../../MKM-DATA/$Temp/kr
      cd ../
    done
  
    for i in *Gas *H2
      do cd $i
      value=`head -n$linerxnG Rxn-Energy | tail -1 | awk '{print $3}'`; echo $i$'\t'$value >> ../../MKM-DATA/$Temp/Grxn
      value=`head -n$linerxnG Rxn-Energy | tail -1 | awk '{print $4}'`; echo $i$'\t'$value >> ../../MKM-DATA/$Temp/Keq
      value=`head -n$linerxnG Rxn-Energy | tail -1 | awk '{print $5}'`; echo $i$'\t'$value >> ../../MKM-DATA/$Temp/kf
      value=`head -n$linerxnG Rxn-Energy | tail -1 | awk '{print $6}'`; echo $i$'\t'$value >> ../../MKM-DATA/$Temp/kr    
      cd ../
    done

cd ../

ii=$(($ii+1))
Temp=$(($LT+($ii-1)*$IT))
linerxn=$(($linerxn+1))
lineKin=$(($lineKin+1))
linerxnG=$(($linerxnG+1))
done
