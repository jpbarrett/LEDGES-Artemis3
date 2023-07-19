#!/bin/sh
#must run this from <LEDGES-Artemis3> repo dir
current_dir=$PWD
pythondir="$current_dir/extern/lunar-RFI/calc_width"
export PYTHONPATH=$pythondir:$PYTHONPATH
export LEDGES_DIR=$current_dir

#manually copy .pkl data file to this directory from submodule (submodule expects this to be in the root path)
echo "copying 2d_interp_h_vs_dB.pkl to $LEDGES_DIR"
cp "$current_dir/extern/lunar-RFI/calc_width/2d_interp_h_vs_dB.pkl" $current_dir
