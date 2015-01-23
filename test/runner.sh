#!/bin/bash

directory=/eos/uscms/store/user/awhitbe1/SuSySubstructureAnalysisNtuples_V10
destination=/uscms_data/d2/fgior8/SubStructure/CMSSW_5_3_14_patch2/src/AWhitbeck/DataSetList
TEXT=$1.txt

for file in `ls -S $directory/$1*`
do
  a=`echo $file | awk -F_ '{print $6}'`
  copiedfiles[$a]=0
done  

for file in `ls -S $directory/$1*`
do
  b=`echo $file | awk -F_ '{print $6}'`
  echo $b
  if [ ${copiedfiles[$b]} -eq 1 ]
  then
    echo "Duplicate"
    continue
  fi
  echo "ok"
  echo "$file" >> $destination/$TEXT
#  echo "dcache:$directory/$file" >> $destination/$TEXT
#  dccp  $directory/$file $destination/$file
  copiedfiles[$b]=1
done                                  

exit 0
