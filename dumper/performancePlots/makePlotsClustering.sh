#!/bin/bash

version1=$1
version2=$2

    dateDir=`date +"%Y-%m-%d"`
    mkdir ~/www/plots/shashlikUpgrade/$dateDir/
    mkdir ~/www/plots/shashlikUpgrade/$dateDir/$version1"-"$version2/
    mkdir ~/www/plots/shashlikUpgrade/$dateDir/$version1"-"$version2/PU140/
    mkdir ~/www/plots/shashlikUpgrade/$dateDir/$version1"-"$version2/noPU/
    cp ~/www/plots/shashlikUpgrade/index.php  ~/www/plots/shashlikUpgrade/$dateDir/

#for dataset in  HggRelVal DYToLLRelVal GluGluToHToGG DYToEE; do
for dataset in  DYToEE GluGluToHToGG; do

    n=`cat lists/$dataset"_"$version1.list |wc -l`
    echo "plotting"
    if [ $n -ne 0 ]
    then
	
	rm plots/h*
	./tmp/compareClusterings finalFiles/$dataset"_noPU_"$version1.root finalFiles/$dataset"_noPU_"$version2.root plots_$dataset"_noPU_"$version1"-"$version2.root
	mkdir ~/www/plots/shashlikUpgrade/$dateDir/$version1"-"$version2/noPU/$dataset
	cp plots/h* ~/www/plots/shashlikUpgrade/$dateDir/$version1"-"$version2/noPU/$dataset/
	cp ~/www/plots/shashlikUpgrade/index.php  ~/www/plots/shashlikUpgrade/$dateDir/$version1"-"$version2/index.php
	cp ~/www/plots/shashlikUpgrade/index.php  ~/www/plots/shashlikUpgrade/$dateDir/$version1"-"$version2/noPU/index.php
	cp ~/www/plots/shashlikUpgrade/index.php  ~/www/plots/shashlikUpgrade/$dateDir/$version1"-"$version2/noPU/$dataset/index.php
	
	rm plots/h*
	./tmp/compareClusterings  finalFiles/$dataset"_"$version1.root  finalFiles/$dataset"_"$version2.root plots_$dataset"_"$version1"-"$version2.root
	mkdir ~/www/plots/shashlikUpgrade/$dateDir/$version1"-"$version2/PU140/$dataset
	cp plots/h* ~/www/plots/shashlikUpgrade/$dateDir/$version1"-"$version2/PU140/$dataset/
	cp ~/www/plots/shashlikUpgrade/index.php  ~/www/plots/shashlikUpgrade/$dateDir/$version1"-"$version2/index.php
	cp ~/www/plots/shashlikUpgrade/index.php  ~/www/plots/shashlikUpgrade/$dateDir/$version1"-"$version2/PU140/index.php
	cp ~/www/plots/shashlikUpgrade/index.php  ~/www/plots/shashlikUpgrade/$dateDir/$version1"-"$version2/PU140/$dataset/index.php
    fi

done

