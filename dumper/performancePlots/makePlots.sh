#!/bin/bash



for version in V20 V21;do
    cd lists
    ./makeLists.sh $version
    cd -
    dateDir=`date +"%Y-%m-%d"`
    mkdir ~/www/plots/shashlikUpgrade/$dateDir/
    cp ~/www/plots/shashlikUpgrade/index.php  ~/www/plots/shashlikUpgrade/$dateDir/

for dataset in  HggRelVal DYToLLRelVal; do
    if [ "$dataset" == "DYToLLRelVal" ]
    then
	./tmp/runCreateHistos2 lists/$dataset"_"$version.list $dataset"_"$version.root   
	./tmp/runCreateHistos2 lists/$dataset"_noPU_"$version.list $dataset_noPU_$version.root
    fi

    if [ "$dataset" == "HggRelVal" ]    
    then
	./tmp/runCreateHistosPhotons lists/$dataset"_"$version.list $dataset"_"$version.root   
	./tmp/runCreateHistosPhotons lists/$dataset"_noPU_"$version.list $dataset_noPU_$version.root
    fi

    rm plots/h*
   ./tmp/plotHistos $dataset"_noPU_"$version.root $dataset"_"$version.root plots_$dataset"_"$version.root
   ./tmp/plotHistos $dataset"_noPU_"$version.root $dataset"_"$version.root plots_$dataset"_"$version.root
    mkdir ~/www/plots/shashlikUpgrade/$dateDir/plots_$dataset"_"$version
    cp plots/h* ~/www/plots/shashlikUpgrade/$dateDir/$version/
    cp ~/www/plots/shashlikUpgrade/index.php  ~/www/plots/shashlikUpgrade/$dateDir/$version/
done

done