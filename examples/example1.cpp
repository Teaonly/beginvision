#include <iostream>

#include <src/image.h>
#include <src/image_io.h>

int main(int argc, char *argv[]) {
    if ( argc < 2) {
        std::cout << "Please input bmp file!" << std::endl;
        return -1;
    }

    bv::Image* img = bv::CreateImageFromBMP( argv[1] );
    
    std::cout << "Width = " << img->color(1).rows() << std::endl;
    std::cout << "Height = " << img->color(1).cols() << std::endl;
    delete img;
}
