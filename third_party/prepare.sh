#!/bin/sh

rm -rf Eigen eigen-eigen-2249f9c22fe8 liblbfgs-master liblbfgs
tar xfz eigen_3.1.3.tgz
mv eigen-eigen-2249f9c22fe8/Eigen ./
unzip liblbfgs-fb51be30d4.zip
mkdir liblbfgs
cp liblbfgs-master/include/lbfgs.h ./liblbfgs/
cp liblbfgs-master/lib/*.h ./liblbfgs/
cp liblbfgs-master/lib/*.c ./liblbfgs/                        
rm -rf eigen-eigen-2249f9c22fe8 liblbfgs-master
