#include <iostream>

#include <bv/image.h>
#include <bv/image_io.h>

int main(int argc, char *argv[]) {
    if ( argc < 2) {
        std::cout << "Please input bmp file!" << std::endl;
        return -1;
    }
    
    bv::ColorImage<3>* img = bv::CreateImageFromBMP( argv[1] );

    std::cout << "Width = " << img->color(0).width() << std::endl;
    std::cout << "Height = " << img->color(0).data.cols() << std::endl;
    
    img->color(1).data *= 2;

    bv::SaveImageToBMP( img, "/tmp/xx.bmp");

    delete img;

    std::cout << "Destory image is also OK" << std::endl; 
}
