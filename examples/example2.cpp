#include <iostream>

#include <bv/image.h>
#include <bv/image_convert.h>

int main(int argc, char *argv[]) {
    if ( argc < 2) {
        std::cout << "Please input bmp file!" << std::endl;
        return -1;
    }
    
    bv::ColorImage<3> colorImage( argv[1] );
    bv::Image img( colorImage.color(0).width(), colorImage.color(0).height() );
    
    bv::Convert::colorImageToGrayImage(colorImage, img);
    
    bv::Image img2( img.height(), img.width() );
    img2.data = img.data.transpose();

    bv::Convert::grayImageToColorImage(img2, colorImage);
    colorImage.SaveImageToBMP("/tmp/xx.bmp");
}
