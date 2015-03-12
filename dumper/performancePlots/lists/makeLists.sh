#!/bin/bash

version=$1

for dataset in DYToLLRelVal HggRelVal DYToEE GluGluToHToGG; do

    eos=/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select
    
    pathnoPU="root://eoscms//eos/cms/store/caf/user/micheli/ShashlikUpgrade/"$version"/"$dataset"/noPU/"
    pathPU="root://eoscms//eos/cms/store/caf/user/micheli/ShashlikUpgrade/"$version"/"$dataset"/PU140/"
    
    
    echo $pathPU
    
    $eos ls /eos/cms/store/caf/user/micheli/ShashlikUpgrade/$version/$dataset/noPU | awk -v pathnoPU="root://eoscms//eos/cms/store/caf/user/micheli/ShashlikUpgrade/"$version"/"$dataset"/noPU/" '{print pathnoPU $1}' > $dataset"_noPU_"$version.list
    $eos ls /eos/cms/store/caf/user/micheli/ShashlikUpgrade/$version/$dataset/PU140 | awk -v pathPU="root://eoscms//eos/cms/store/caf/user/micheli/ShashlikUpgrade/"$version"/"$dataset"/PU140/" '{print pathPU $1}' > $dataset"_"$version.list
done