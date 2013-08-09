#include <iostream>

#include <bv/image.h>
#include <bv/image_convert.h>
#include <bv/filter.h>

int main(int argc, char *argv[]) {
    if ( argc < 2) {
        std::cout << "Please input bmp file!" << std::endl;
        return -1;
    }
    
    int ret;  
    bv::ColorImage<3> colorImage( argv[1] );
    bv::Image img( colorImage.color(0).width(), colorImage.color(0).height() );
    Eigen::MatrixXd gray(img.width(), img.height());


    bv::Convert::colorImageToGrayImage(colorImage, img);
    bv::Convert::grayImageToMatrix(img, gray);

    Eigen::MatrixXd small(img.width()*2, img.height()*2);
    bv::Image imgSmall( small.rows(), small.cols() );
    bv::ColorImage<3> colorImageSmall(imgSmall.width(), imgSmall.height());

    bv::Convert::resizeImage(gray, small);
    
    bv::Convert::matrixToGrayImage(small, imgSmall);    
    bv::Convert::grayImageToColorImage(imgSmall, colorImageSmall);
    colorImageSmall.SaveImageToBMP("/tmp/xx.bmp");
}
