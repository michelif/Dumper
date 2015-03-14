#!/bin/bash



for version in V22_b;do
    dateDir=`date +"%Y-%m-%d"`
    mkdir ~/www/plots/shashlikUpgrade/$dateDir/
    mkdir ~/www/plots/shashlikUpgrade/$dateDir/$version/
    cp ~/www/plots/shashlikUpgrade/index.php  ~/www/plots/shashlikUpgrade/$dateDir/

for dataset in  GluGluToHToGG DYToEE; do
    if [ "$dataset" == "DYToEE" ]
    then
	for i in $(seq 0 20); do
	    ./tmp/runCreateHistos2 lists/$dataset"_"$version"_"$i $dataset"_"$version"_"$i.root   
	    ./tmp/runCreateHistos2 lists/$dataset"_noPU_"$version"_"$i $dataset"_noPU_"$version"_"$i.root
	done
    echo "DY ok"
    fi

    if [ "$dataset" == "GluGluToHToGG" ]    
    then
	for i in $(seq 0 2); do
	    ./tmp/runCreateHistosPhotons lists/$dataset"_"$version$i $dataset"_"$version"_"$i.root   
	    ./tmp/runCreateHistosPhotons lists/$dataset"_noPU_"$version$i $dataset"_noPU_"$version"_"$i.root
	done
    echo "Hgg ok" 
    fi

    rm plots/h*
    hadd $dataset"_noPU_"$version.root $dataset"_noPU_"$version"_*".root
    hadd $dataset"_"$version.root $dataset"_"$version"_*".root
   ./tmp/plotHistos $dataset"_noPU_"$version.root $dataset"_"$version.root plots_$dataset"_"$version.root
    mkdir ~/www/plots/shashlikUpgrade/$dateDir/$version/$dataset
    cp plots/h* ~/www/plots/shashlikUpgrade/$dateDir/$version/$dataset/
    cp ~/www/plots/shashlikUpgrade/index.php  ~/www/plots/shashlikUpgrade/$dateDir/$version/index.php
    cp ~/www/plots/shashlikUpgrade/index.php  ~/www/plots/shashlikUpgrade/$dateDir/$version/$dataset/index.php
done

done