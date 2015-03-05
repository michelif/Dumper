#!/bin/bash
listfile=" "
for i in `cat lists/SingleMuPt100RelVal_$1.list`
do

listfile+=$i
listfile+=" "
done
echo $listfile

hadd  SMu$1.root $listfile



