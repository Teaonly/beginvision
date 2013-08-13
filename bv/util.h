#ifndef _BV_UTIL_H_
#define _BV_UTIL_H_

#include <Eigen/Core>

namespace bv {

namespace Util {

template<typename T> 
T min(T& v1, T& v2) {
    if ( v1 < v2)
        return v1;
    return v2;
}

double interp2(double x, double y, Eigen::MatrixXd& img) {
    int leftX = (int) (x);
    int rightX = (int) (x + 1);
    int topY = (int) (y);
    int bottomY = (int) (y + 1);

    double top = (x - leftX) * img(leftX, topY) + (rightX - x) * img(rightX, topY);
    
    double bottom = (x - leftX) * img(leftX, bottomY) + (rightX - x) * img(rightX, bottomY);
                
    return (y - topY) * top + (bottomY - y) * bottom;
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
