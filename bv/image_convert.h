#ifndef _BV_IMAGE_CONVERT_H_
#define _BV_IMAGE_CONVERT_H_

#include "image.h"

namespace bv {

class Convert {
public:
    static void colorImageToGrayImage(ColorImage<3>& in, Image& out) {
        for(int y = 0; y < in.height(); y++) {
            for ( int x = 0; x < in.width(); x++) {
                unsigned char r,g,b;
                r = in.color(0).data(x,y);
                g = in.color(1).data(x,y);
                b = in.color(2).data(x,y);

                unsigned char gray;
                rgbToGray(r,g,b, gray);
                out.data(x,y) = gray;
            }
        }
    }
    
    static void grayImageToColorImage(Image& in, ColorImage<3>& out) {
        for ( int i = 0; i < 3; i++) {
            out.color(i).data = in.data;
        }
    }

    template<typename T>
    static void rgbToGray(unsigned char r, unsigned char g, unsigned char b, T& output) {
        unsigned int y;
        y = (r*77)+(g*151)+(b*28);
        y = y >> 8;
        output = (T)y;
    }
    
    #define SATURATE(a,min,max) ((a)<(min)?(min):((a)>(max)?(max):(a)))
    static void yuvToRgb(unsigned char inY, unsigned char inU, unsigned char inV,
                unsigned char& R, unsigned char& G, unsigned char& B) {
        int Gi = 0, Bi = 0;
        int Y = 9535*(inY-16);
        int U = inU - 128;
        int V = inV - 128;
        int Ri = (Y + 13074*V) >> 13;
        Ri = SATURATE(Ri, 0, 255);
        Gi = (Y - 6660*V - 3203*U) >> 13;
        Gi = SATURATE(Gi, 0, 255);
        Bi = (Y + 16531*U) >> 13;
        Bi = SATURATE(Bi, 0, 255);
        R = Ri;
        G = Gi;
        B = Bi;
    }  

};

}
#endif
