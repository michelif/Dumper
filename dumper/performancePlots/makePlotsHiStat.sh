#!/bin/bash



for version in V21_b_2 V20_b_2;do
    dateDir=`date +"%Y-%m-%d"`
    mkdir ~/www/plots/shashlikUpgrade/$dateDir/
    mkdir ~/www/plots/shashlikUpgrade/$dateDir/$version/
    cp ~/www/plots/shashlikUpgrade/index.php  ~/www/plots/shashlikUpgrade/$dateDir/

for dataset in  GluGluToHToGG DYToEE; do
#for dataset in  DYToEE; do
    if [ "$dataset" == "DYToEE" ]
    then
#	for i in $(seq 0 20); do
#	    ./tmp/runCreateHistos2 lists/$dataset"_"$version"_"$i $dataset"_"$version"_"$i.root   
#	    ./tmp/runCreateHistos2 lists/$dataset"_noPU_"$version"_"$i $dataset"_noPU_"$version"_"$i.root
#	done
    echo "DY ok"
    fi

    if [ "$dataset" == "GluGluToHToGG" ]    
    then
#	for i in $(seq 0 2); do
#	    ./tmp/runCreateHistosPhotons lists/$dataset"_"$version"_"$i $dataset"_"$version"_"$i.root   
#	    ./tmp/runCreateHistosPhotons lists/$dataset"_noPU_"$version"_"$i $dataset"_noPU_"$version"_"$i.root
#	done
    echo "Hgg ok" 
    fi

    n=`cat lists/$dataset"_"$version.list |wc -l`
    echo "plotting"
    if [ $n -ne 0 ]
    then
	rm plots/h*
#    hadd $dataset"_noPU_"$version.root $dataset"_noPU_"$version"_"*.root
#    hadd $dataset"_"$version.root $dataset"_"$version"_"*.root
	./tmp/plotHistos finalFiles/$dataset"_noPU_"$version.root finalFiles/$dataset"_"$version.root finalFiles/plots_$dataset"_"$version.root
	mkdir ~/www/plots/shashlikUpgrade/$dateDir/$version/$dataset
	cp plots/h* ~/www/plots/shashlikUpgrade/$dateDir/$version/$dataset/
	cp ~/www/plots/shashlikUpgrade/index.php  ~/www/plots/shashlikUpgrade/$dateDir/$version/index.php
	cp ~/www/plots/shashlikUpgrade/index.php  ~/www/plots/shashlikUpgrade/$dateDir/$version/$dataset/index.php

	rm plots_splitted/*
	./tmp/plotDistributions $dataset"_noPU_"$version.root plots_splitted_$dataset"_"$version
	mkdir ~/www/plots/shashlikUpgrade/$dateDir/$version/$dataset/
	mkdir ~/www/plots/shashlikUpgrade/$dateDir/$version/$dataset/noPU/
	cp plots_splitted/* ~/www/plots/shashlikUpgrade/$dateDir/$version/$dataset/noPU/
	cp ~/www/plots/shashlikUpgrade/index.php  ~/www/plots/shashlikUpgrade/$dateDir/$version/index.php
	cp ~/www/plots/shashlikUpgrade/index.php  ~/www/plots/shashlikUpgrade/$dateDir/$version/$dataset/index.php
	cp ~/www/plots/shashlikUpgrade/index.php  ~/www/plots/shashlikUpgrade/$dateDir/$version/$dataset/noPU/index.php

	rm plots_splitted/*
	./tmp/plotDistributions $dataset"_"$version.root plots_splitted_$dataset"_"$version
	mkdir ~/www/plots/shashlikUpgrade/$dateDir/$version/$dataset/
	mkdir ~/www/plots/shashlikUpgrade/$dateDir/$version/$dataset/PU140/
	cp plots_splitted/* ~/www/plots/shashlikUpgrade/$dateDir/$version/$dataset/PU140/
	cp ~/www/plots/shashlikUpgrade/index.php  ~/www/plots/shashlikUpgrade/$dateDir/$version/index.php
	cp ~/www/plots/shashlikUpgrade/index.php  ~/www/plots/shashlikUpgrade/$dateDir/$version/$dataset/index.php
	cp ~/www/plots/shashlikUpgrade/index.php  ~/www/plots/shashlikUpgrade/$dateDir/$version/$dataset/PU140/index.php


    fi
done

done