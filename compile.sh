#!/bin/bash

rm -rf deploy
mkdir build
cd build 
cmake ../
make -j $NCORES
cd ..

mkdir deploy
cp ./build/strongmap deploy/
cp ./build/ecomap deploy/
cp ./build/fastmap deploy/
cp ./build/fastestmap deploy/
cp ./build/ssocialmap deploy/
cp ./build/esocialmap deploy/
cp ./build/fsocialmap deploy/
cp ./build/ffsocialmap deploy/

rm -rf build
