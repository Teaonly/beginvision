#include <iostream>
#include <vector>

#include <bv/image.h>
#include <bv/image_convert.h>
#include <bv/filter.h>
#include <bv/detector_sift.h>
#include <bv/descriptor_sift.h>

int main(int argc, char *argv[]) {
    if ( argc < 3) {
        std::cout << "Please input two bmp file!" << std::endl;
        return -1;
    }
    
    int ret;  
    bv::ColorImage<3> c1( argv[1] );
    bv::ColorImage<3> c2( argv[2] );

    bv::Image i1( c1.color(0).width(), c1.color(0).height() );
    Eigen::MatrixXd g1(i1.width(), i1.height());
    bv::Image i2( c2.color(0).width(), c2.color(0).height() );
    Eigen::MatrixXd g2(i2.width(), i2.height());

    bv::Convert::colorImageToGrayImage(c1, i1);     //i1 = c1.color(1);
    bv::Convert::grayImageToMatrix( i1, g1);
    bv::Convert::colorImageToGrayImage(c2, i2);
    bv::Convert::grayImageToMatrix( i2, g2);

    bv::DT_Sift dt1;
    dt1.run(g1);
    bv::DS_Sift ds1(dt1); 
 
    bv::DT_Sift dt2;
    dt2.run(g2);
    bv::DS_Sift ds2(dt2); 
    
    std::cout << "Begin compute match...." << std::endl;
    std::vector<int> result;
    ds1.matchWith(ds2, result);
    ds1.showMatch(ds2, result, g1, g2);

    return 0;
}
