#!/bin/bash
echo 'Searching for ' $1 
while read LINE
  do
  # echo $LINE
  grep -H -n $1 $LINE
done < $2
