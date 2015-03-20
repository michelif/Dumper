#!/bin/bash


for version in V23 V22_b V21_b V20_b;do
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
#	./tmp/runCreateHistos2 lists/$dataset"_"$version.list $dataset"_"$version.root   
#	./tmp/runCreateHistos2 lists/$dataset"_noPU_"$version.list $dataset"_noPU_"$version.root
    echo "DY ok"
    fi

    if [ "$dataset" == "ZEERelVal" ]
    then
#	./tmp/runCreateHistos2 lists/$dataset"_"$version.list $dataset"_"$version.root   
#	./tmp/runCreateHistos2 lists/$dataset"_noPU_"$version.list $dataset"_noPU_"$version.root
    echo "ZEE ok"
    fi


    if [ "$dataset" == "HggRelVal" ]    
    then
#	./tmp/runCreateHistosPhotons lists/$dataset"_"$version.list $dataset"_"$version.root   
#	./tmp/runCreateHistosPhotons lists/$dataset"_noPU_"$version.list $dataset"_noPU_"$version.root
    echo "Hgg ok" 
    fi

    n=`cat lists/$dataset"_"$version.list |wc -l`
    echo "plotting"
    if [ $n -ne 0 ]
    then
	rm plots/h*
	./tmp/plotHistos $dataset"_noPU_"$version.root $dataset"_"$version.root plots_$dataset"_"$version.root
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