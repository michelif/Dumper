#!/bin/bash

version=$1



eos=/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select

pathnoPU="root://eoscms//eos/cms/store/caf/user/micheli/ShashlikUpgrade/"$version"/DYToLLRelVal/noPU/"
pathPU="root://eoscms//eos/cms/store/caf/user/micheli/ShashlikUpgrade/"$version"/DYToLLRelVal/PU140/"


echo $pathPU

$eos ls /eos/cms/store/caf/user/micheli/ShashlikUpgrade/$version/DYToLLRelVal/noPU | awk -v pathnoPU="root://eoscms//eos/cms/store/caf/user/micheli/ShashlikUpgrade/"$version"/DYToLLRelVal/noPU/" '{print pathnoPU $1}' > DYToLLRelVal_noPU_$version.list
$eos ls /eos/cms/store/caf/user/micheli/ShashlikUpgrade/$version/DYToLLRelVal/PU140 | awk -v pathPU="root://eoscms//eos/cms/store/caf/user/micheli/ShashlikUpgrade/"$version"/DYToLLRelVal/PU140/" '{print pathPU $1}' > DYToLLRelVal_$version.list