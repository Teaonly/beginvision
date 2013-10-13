#ifndef _BV_IMAGE_H_
#define _BV_IMAGE_H_

#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <assert.h>

#include "beginvision.h"
#include "bmp.h"

namespace bv {

class Image {
public:
    Image(int width, int height, int maxValue = 255) : data(width, height), scale(maxValue) {
        
    }
    ~Image() {
    }
    int width() const {
        return data.rows();
    }
    int height() const {
        return data.cols();
    }

    Eigen::MatrixXi data;
    int scale; 
};

template<int ColorNumber>
class ColorImage {
public:
    ColorImage(int width, int height, int scale=255) {
        if ( ColorNumber <= 0) {
            assert(false);
        }

        for(int i = 0; i < ColorNumber; i++) {
            color_[i] = new Image(width, height, scale);
        }
    }
    ColorImage(const std::string& fileName) {
        bitmap_file_header bmfh;
        bitmap_info_header bmih;
        std::vector<unsigned char> data; 
        
        if ( ColorNumber != 3 ) {
            assert(false);
        } 
        for(int i = 0; i < ColorNumber; i++) {
            color_[i] = NULL;
        }

        // parse BMP file
        {
            std::ifstream in(fileName.c_str(),std::ios::binary);
            if (!in || !bmfh.read(in)) {
                return;
            }

            in.read((char*)&bmih, sizeof(bitmap_info_header));
            if (!in || bmih.biWidth <= 0 || bmih.biHeight <= 0 || bmih.biCompression != 0)
                return;
            if ( bmih.biBitCount != 24 || bmih.biSize != sizeof(bitmap_info_header) )
                return;
            
            unsigned int padding = (4-((bmih.biWidth*3) & 3)) &3;
            if ( (bmih.biWidth+padding) * bmih.biHeight >= bmih.biSizeImage )
                return;

            data.resize(bmih.biSizeImage);
            in.read((char*)&*data.begin(), data.size());
        }
        
        int width = bmih.biWidth;
        int height = bmih.biHeight; 

        for(int i = 0; i < ColorNumber; i++) {
            color_[i] = new Image(width, height);
        }
        
        unsigned int padding = (4-((bmih.biWidth*3) & 3)) &3; 
        int bmp_index = 0;
        for(int l = height-1; l >= 0; l -- ) { 
            for (int x = 0; x < width; x++) {
                unsigned r,g,b;
                b = data[bmp_index];
                bmp_index++;
                g = data[bmp_index];
                bmp_index++;
                r = data[bmp_index];
                bmp_index++;
                
                color(0).data(x, l) = r;
                color(1).data(x, l) = g;
                color(2).data(x, l) = b;
            }
            bmp_index += padding;
        }       
    }    

    ~ColorImage() {
        for(int i = 0; i < ColorNumber; i++) {
            if ( color_[i] != NULL) {
                delete color_[i];
            }    
        }
    }
 
    int width() const {
        return color_[0]->width();
    }
    int height() const {
        return color_[0]->height();
    }
    
    Image& color(int c) {
        return *(color_[c]);
    }

    int SaveImageToBMP(const std::string& fileName) {
        bitmap_file_header bmfh;
        bitmap_info_header bmih;
        std::vector<unsigned char> data; 

        int width = color(0).width();
        int height = color(0).height();
        
        if ( ColorNumber != 3)
            return BV_ERROR;


        // config the info header
        std::fill((char*)&bmih,(char*)&bmih+sizeof(bmih),0);
        bmih.biBitCount = 24;
        bmih.biPlanes = 1;
        bmih.biCompression = 0;;
        bmih.biWidth = width;
        bmih.biHeight = height;
        unsigned int padding = (4-((bmih.biWidth*3) & 3)) &3;
        bmih.biSizeImage = (width + padding) * height * 3;
        bmih.biSize = sizeof(bitmap_info_header);
        data.resize(bmih.biSizeImage);
        int bmp_index = 0;
        for(int l = height - 1; l >= 0; l -- ) {
            for (int x = 0; x < width; x++) {
                unsigned r,g,b;
               
                r = color(0).data(x,l);
                g = color(1).data(x,l);
                b = color(2).data(x,l);
                
                data[bmp_index] = b;
                bmp_index++;
                data[bmp_index] = g;
                bmp_index++;
                data[bmp_index] = r;
                bmp_index++;
            }
            bmp_index += padding;
        } 
       
        // write to file 
        std::ofstream out(fileName.c_str(),std::ios::binary);
        if (!bmfh.write(out))
            return BV_ERROR_IO;

        out.write((const char*)&bmih,sizeof(bitmap_info_header));
        out.write((const char*)&*data.begin(),data.size());

        return BV_OK;
    }

private:
    Image* color_[ColorNumber];
};

}
#endif
