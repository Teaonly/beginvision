#ifndef _BV_UTIL_H_
#define _BV_UTIL_H_

#include "beginvision.h"
#include <Eigen/Core>

namespace bv {

namespace Util {

template<typename T> 
T min(const T& v1,const T& v2) {
    if ( v1 < v2)
        return v1;
    return v2;
}

template<typename T> 
T max(const T& v1, const T& v2) {
    if ( v1 > v2)
        return v1;
    return v2;
}

double interp2(double x, double y, Eigen::MatrixXd& img) {
    int leftX = (int) (x);
    int rightX = leftX + 1;
    int topY = (int) (y);
    int bottomY = topY + 1;

    double top = (x - leftX) * img(rightX, topY) + (rightX - x) * img(leftX, topY);
    double bottom = (x - leftX) * img(rightX, bottomY) + (rightX - x) * img(leftX, bottomY);
                
    return (y - topY) * bottom + (bottomY - y) * top;
}

double mod2pi(double x) {
    while (x > (double)(2 * PI)) x -= (double) (2 * PI) ;
    while (x < 0.0) x += (double) (2 * PI);
    return x ;  
}

void saveAsImage(Eigen::MatrixXd& x, const std::string& fileName) {
    Image tmp(x.rows(), x.cols());
    ColorImage<3> tmp2(x.rows(), x.cols());
    Convert::matrixToGrayImage(x, tmp);
    Convert::grayImageToColorImage(tmp, tmp2); 
    tmp2.SaveImageToBMP(fileName);
}

}
}
#endif
