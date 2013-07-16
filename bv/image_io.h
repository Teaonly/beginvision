#ifndef _BV_IMAGE_IO_H_
#define _BV_IMAGE_IO_H_

#ifdef _MSC_VER
typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
typedef __int16 int16_t;
typedef unsigned __int16 uint16_t;
typedef __int64 int64_t;
typedef unsigned __int64 uint64_t;
#else
#include <stdint.h>
#endif

#include <fstream>
#include <vector>
#include "image.h"

namespace bv {

struct bitmap_file_header {
public:    
    uint32_t   bfSize;
    bool read(std::istream& in) {
        uint16_t bfType;
        in.read((char*)&bfType,2);
        if (bfType !=  0x4D42) // BM
            return false;
        in.read((char*)&bfSize,4);
        uint32_t dummy;
        in.read((char*)&dummy,4);
        in.read((char*)&dummy,4);
        return in && dummy== 54;//bfOffBits
    }
    bool write(std::ostream& out) const {
        uint16_t bfType = 0x4D42;
        out.write((const char*)&bfType,2);
        out.write((const char*)&bfSize,4);
        uint32_t dummy = 0;
        out.write((const char*)&dummy,4);
        dummy = 54;
        out.write((const char*)&dummy,4);//bfOffBits
        return !(!out);
    }
};

struct bitmap_info_header {
    uint32_t      biSize;
    int32_t       biWidth;
    int32_t       biHeight;
    unsigned short       biPlanes;
    unsigned short       biBitCount;
    uint32_t      biCompression;
    uint32_t      biSizeImage;
    int32_t       biXPelsPerMeter;
    int32_t       biYPelsPerMeter;
    uint32_t      biClrUsed;
    uint32_t      biClrImportant;
};

Image* CreateImageFromBMP(const std::string& fileName) {
        bitmap_file_header bmfh;
        bitmap_info_header bmih;
        std::vector<unsigned char> data; 
        
        // parse BMP file
        {
            std::ifstream in(fileName.c_str(),std::ios::binary);
            if (!in || !bmfh.read(in)) {
                return NULL;
            }

            in.read((char*)&bmih, sizeof(bitmap_info_header));
            if (!in || bmih.biWidth <= 0 || bmih.biHeight <= 0 || bmih.biCompression != 0)
                return NULL;
            if ( bmih.biBitCount != 24 || bmih.biSize != sizeof(bitmap_info_header) )
                return NULL;
            
            unsigned int padding = (4-((bmih.biWidth*3) & 3)) &3;
            if ( (bmih.biWidth+padding) * bmih.biHeight >= bmih.biSizeImage )
                return NULL;

            data.resize(bmih.biSizeImage);
            in.read((char*)&*data.begin(), data.size());
        }
        
        int width = bmih.biWidth;
        int height = bmih.biHeight; 

        Image* img = new Image(width, height); 
        
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
                
                img->color(0)(x, l) = r;
                img->color(0)(x, l) = g;
                img->color(0)(x, l) = b;
            }
            bmp_index += padding;
        }           

        return img;
    }

}   //end of namespace

#endif
