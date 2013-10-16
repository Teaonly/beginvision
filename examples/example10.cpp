#include <iostream>
#include <vector>

#include <bv/image.h>
#include <bv/image_convert.h>
#include <bv/filter.h>
#include <bv/detector_sift.h>
#include <bv/descriptor_sift.h>

int main(int argc, char *argv[]) {
    if ( argc < 2) {
        std::cout << "Please input bmp file!" << std::endl;
        return -1;
    }
    
    int ret;  
    bv::ColorImage<3> c1( argv[1] );

    bv::Image i1( c1.color(0).width(), c1.color(0).height() );
    Eigen::MatrixXd g1(i1.width(), i1.height());

    bv::Convert::colorImageToGrayImage(c1, i1);
    //i1 = c1.color(1);
    bv::Convert::grayImageToMatrix( i1, g1);

    bv::DT_Sift detector;
    detector.run(g1);
        
    
    return 0;
}
