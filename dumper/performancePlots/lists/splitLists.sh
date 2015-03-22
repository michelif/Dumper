#!/bin/bash

version=$1

split -a 1 -l 10 -d GluGluToHToGG_$version.list GluGluToHToGG_$version"_"
split -a 1 -l 10 -d DYToEE_$version.list DYToEE_$version"_"
split -l 10 -d DYToEE_$version.list DYToEE_$version"_"

split -a 1 -l 10 -d GluGluToHToGG_noPU_$version.list GluGluToHToGG_noPU_$version"_"
split -a 1 -l 10 -d DYToEE_noPU_$version.list DYToEE_noPU_$version"_"
split -l 10 -d DYToEE_noPU_$version.list DYToEE_noPU_$version"_"