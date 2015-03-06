#!/bin/bash

for version in V16 V17 V18;do
cd lists
./makeLists.sh $version
cd -
#./tmp/runCreateHistos2 lists/DYToLLRelVal_$version.list DYToLLRelVal_$version.root
#./tmp/runCreateHistos2 lists/DYToLLRelVal_noPU_$version.list DYToLLRelVal_noPU_$version.root
rm plots/h*
./tmp/plotHistos DYToLLRelVal_noPU_$version.root DYToLLRelVal_$version.root plots_DYToLLRelVal_$version.root
cp plots/h* ~/www/plots/shashlikUpgrade/2015-03-05/$version/

done