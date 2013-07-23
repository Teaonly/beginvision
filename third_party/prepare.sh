#!/bin/sh

rm -rf Eigen eigen-eigen-2249f9c22fe8
tar xfz eigen_3.1.3.tgz
mv eigen-eigen-2249f9c22fe8/Eigen ./
rm -rf eigen-eigen-2249f9c22fe8
