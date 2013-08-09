#include <iostream>

#include <bv/image.h>

int main(int argc, char *argv[]) {
    if ( argc < 2) {
        std::cout << "Please input bmp file!" << std::endl;
        return -1;
    }
    
    bv::ColorImage<3> img( argv[1] );

    std::cout << "Width = " << img.color(0).width() << std::endl;
    std::cout << "Height = " << img.color(0).height() << std::endl;
    
    img.color(1).data *= 2;

    img.SaveImageToBMP("/tmp/xx.bmp");
}
