#!/bin/bash

versions=$1

for version in $versions;do
    cd lists
    ./makeLists.sh $version
    cd -
    dateDir=`date +"%Y-%m-%d"`
    mkdir ~/www/plots/shashlikUpgrade/$dateDir/
    mkdir ~/www/plots/shashlikUpgrade/$dateDir/$version/
    cp ~/www/plots/shashlikUpgrade/index.php  ~/www/plots/shashlikUpgrade/$dateDir/

for dataset in  HggRelVal DYToLLRelVal ZEERelVal; do
    if [ "$dataset" == "DYToLLRelVal" ]
    then
	./tmp/runCreateHistos2 lists/$dataset"_"$version.list $dataset"_"$version.root   
	./tmp/runCreateHistos2 lists/$dataset"_noPU_"$version.list $dataset"_noPU_"$version.root
    echo "DY ok"
    fi

    if [ "$dataset" == "ZEERelVal" ]
    then
	./tmp/runCreateHistos2 lists/$dataset"_"$version.list $dataset"_"$version.root   
	./tmp/runCreateHistos2 lists/$dataset"_noPU_"$version.list $dataset"_noPU_"$version.root
    echo "DY ok"
    fi


    if [ "$dataset" == "HggRelVal" ]    
    then
	./tmp/runCreateHistosPhotons lists/$dataset"_"$version.list $dataset"_"$version.root   
	./tmp/runCreateHistosPhotons lists/$dataset"_noPU_"$version.list $dataset"_noPU_"$version.root
    echo "Hgg ok" 
    fi

    rm plots/h*
   ./tmp/plotHistos $dataset"_noPU_"$version.root $dataset"_"$version.root plots_$dataset"_"$version.root
    mkdir ~/www/plots/shashlikUpgrade/$dateDir/$version/$dataset
    cp plots/h* ~/www/plots/shashlikUpgrade/$dateDir/$version/$dataset/
    cp ~/www/plots/shashlikUpgrade/index.php  ~/www/plots/shashlikUpgrade/$dateDir/$version/index.php
    cp ~/www/plots/shashlikUpgrade/index.php  ~/www/plots/shashlikUpgrade/$dateDir/$version/$dataset/index.php
done

done