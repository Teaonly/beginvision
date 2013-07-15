#ifndef _BV_BMP_H_
#define _BV_BMP_H_

#include <fstream>
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

namespace tea {
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
}   //end of namespace

#endif
