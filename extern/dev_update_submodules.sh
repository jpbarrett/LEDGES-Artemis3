#!/bin/sh

#simple script to update the canned copy of the git submodule (lunar-RFI)
#must be run from the <LEDGES-Artemis3>/extern directory
#make sure to clone the submodules first.

cd ./submodules/lunar-RFI
git archive --format=tar.gz HEAD > ../../lunar-RFI.tar.gz
cd ../../
mkdir -p lunar-RFI 
cd ./lunar-RFI
tar -xzvf ../lunar-RFI.tar.gz 
cd ../
rm ./lunar-RFI.tar.gz
