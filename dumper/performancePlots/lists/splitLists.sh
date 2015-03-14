#!/bin/bash

version=$1

split -l 10 -d GluGluToHToGG_$version.list GluGluToHToGG_$version"_"
split -l 10 -d DYToEE_$version.list DYToEE_$version"_"