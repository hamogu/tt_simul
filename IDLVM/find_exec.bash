#!/bin/bash
echo 'Searching for ' $1 
while read LINE
  do
  # echo $LINE
  grep -H -n $1 $LINE
done < /data/hspc62/steh305/idl/tt_simul/IDLVM/pro_filelist_tt_spectra_list.txt
