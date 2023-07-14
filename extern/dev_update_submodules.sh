#!/bin/sh

#simple script for developers to update the canned copy of the git submodules (date, json, pybind11)
#must be run from the <hops-git>/extern directory
#make sure to clone the submodules first!

cd ./submodules/lunar-RFI
git archive --format=tar.gz HEAD > ../../lunar-RFI.tar.gz
cd ../../
mkdir -p lunar-RFI 
cd ./lunar-RFI
tar -xzvf ../lunar-RFI.tar.gz 
cd ../
rm ./lunar-RFI.tar.gz

