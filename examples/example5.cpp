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
    ret = bv::Convert::grayImageToMatrix( img, gray);   
    //ret = bv::Filter::average(gray, gray, 5);
    ret = bv::Filter::gaussianBlur(gray, gray, 9, 3.0);
    ret = bv::Convert::matrixToGrayImage(gray, img);    
    
    bv::Convert::grayImageToColorImage(img, colorImage);
    colorImage.SaveImageToBMP("/tmp/xx.bmp");
}
