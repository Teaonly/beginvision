#!/bin/sh

rm -rf Eigen eigen-eigen-2249f9c22fe8 liblbfgs-master lbfgs
tar xfz eigen_3.1.3.tgz
mv eigen-eigen-2249f9c22fe8/Eigen ./
unzip iblbfgs-fb51be30d4.zip
mkdir lbfgs
cp liblbfgs-master/include/lbfgs.h ./lbfgs/
cp liblbfgs-master/lib/*.h ./lbfgs/
cp liblbfgs-master/lib/*.c ./lbfgs/                        
rm -rf eigen-eigen-2249f9c22fe8 liblbfgs-master
