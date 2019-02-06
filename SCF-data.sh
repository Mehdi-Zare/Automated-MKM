#!/bin/bash

# this Code go to all Ads Site and TS directory to prepare SCF-Data file for each of them

echo " enter LT "
read LT
echo " eneter HT "
read HT
echo " eneter IT "
read IT

  echo $LT >  ~/MKM-Codes/Temp-File
  echo $HT >> ~/MKM-Codes/Temp-File
  echo $IT >> ~/MKM-Codes/Temp-File

for i in *

  do cd $i

    if [[ $i == Site ]]; then
       echo " We are in site directory"
       chmod +x ~/MKM-Codes/Get-Data-Site.sh
       ~/MKM-Codes/Get-Data-Site.sh 

    elif [[ $i == Ads* ]]; then
       echo " We are in Ads directory"
         chmod +x ~/MKM-Codes/Get-Data-Ads-Gas.sh
         ~/MKM-Codes/Get-Data-Ads-Gas.sh
    
    elif [[ $i == Gas* ]]; then
       echo " We are in Gas directory"
         chmod +x ~/MKM-Codes/Get-Data-Ads-Gas.sh
         ~/MKM-Codes/Get-Data-Ads-Gas.sh
         mkdir Part-Func
         cp CONTCAR freq.dat Gas-Data Part-Func/
           cd Part-Func
           cp ~/MKM-Codes/Temp-File .
           ~/MKM-Codes/partcal
           cd ../
           
    elif [[ $i == TS*  ]]; then
       echo "We are in TS directory"    
         chmod +x ~/MKM-Codes/Get-Data-TS.sh
         ~/MKM-Codes/Get-Data-TS.sh
    else
       echo "You have one directory more than Site, Ads and TS"

    fi

    cd ../
done  
