#ifndef _TEA_HISTOGRAM_H_
#define _TEA_HISTOGRAM_H_

#include <vector>
#include "image.h"

namespace bv {

class Histogram {
public:    
    Histogram(Image& img) {
        float tv = 0.0;
        hist_.resize(255, tv);

        for(int y = 0; y < img.width(); y++) {
            for (int x = 0; x < img.height(); x++) {
                unsigned char charValue = img.data(y,x) % 256;
                hist_[charValue] = hist_[charValue] + 1;
            }
        }

    }
    
    float& value(unsigned char v) {
        return hist_[v];
    }

    float value(unsigned char v) const {
        return hist_[v];
    }
    
    // if float is float 
    void normalization() {
        float s = sum();
        for(int i = 0; i < 255; i++) {
            hist_[i] = hist_[i] / s;
        }
    }

    float sum() {
        float s = 0;
        for(int i = 0; i < 255; i++) {
            s = s + hist_[i];
        }
        return s;
    }

    void equaliseMatrix(Image& img) {
        float sumh = 0.0;
        unsigned char table[256];
        for(int i = 0; i <= 255; i++) {
            if ( hist_[i] < 1) {    // is normalized
                sumh = sumh + hist_[i] * (img.width()*img.height());       
            } else {
                sumh = sumh + hist_[i];
            }
            table[i] = (unsigned char)( 255.0/(img.width()*img.height()) * sumh + 0.0001 );
        }

        for(int y = 0; y < img.height(); y++) {
            for (int x = 0; x < img.width(); x++) {
                 img.data(x,y) = table[ img.data(x,y) ];
            }
        }
    }

    unsigned char autoThreshold() {
        normalization();

        // get the expection 
        float e = 0.0;
        for(int i = 0; i < 256; i++) {
            e = e + i * hist_[i];
        }

        // find the left/right expection
        float w1, w2, u1, u2;
        unsigned char thres = 0;
        float thresValue = 0;
        w1 = hist_[0];
        w2 = 1 - w1;
        u1 = 0;
        u2 = e;
        for(int i = 1; i < 256; i++) {
            // w1/w2 left/right sum
            w1 = w1 + hist_[i];
            w2 = 1 - w1;

            // w1/w2 left/right expection
            u1 = u1 + i * hist_[i];
            u2 = e - u1;

            // varB is optimazition value
            //float varB = (w1*e-u1) * (w1*e-u1) / (w1*w2) ;
            float varB = (u1-u2) * (u1-u2) * (w1*w2); 
            if ( varB > thresValue ) {
                thresValue = varB;
                thres = i;
            } else if ( varB < thresValue ) {
                break;
            }
        }    
        
        return thres;
    }

private:
    std::vector<float> hist_;
};

}
#endif
