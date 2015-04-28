#ifndef _TEA_HISTOGRAM_H_
#define _TEA_HISTOGRAM_H_

#include <vector>
#include "image.h"

namespace bv {

class Histogram {
public:    
    Histogram(Image& img) {
        float tv = 0.0;
        hist_.resize(img.scale + 1, tv);

        for(int x = 0; x < img.width(); x++) {
            for (int y = 0; y < img.height(); y++) {
                unsigned int charValue = img.data(x,y) % (img.scale+1);
                hist_[charValue] = hist_[charValue] + 1;
            }
        }
    }
    
    Histogram(Eigen::MatrixXd& mat, double scale, unsigned int bin) {
        float tv = 0.0;
        hist_.resize(bin, tv);
        
        for(int x = 0; x < mat.rows(); x++) {
            for(int y = 0; y < mat.cols(); y++) {
                unsigned int charValue = (int)(mat(x,y)/scale * bin) % bin;
                hist_[charValue] = hist_[charValue] + 1;
            }
        }
    }
   
    unsigned int size() {
        return hist_.size();
    } 

    float& value(unsigned int v) {
        return hist_[v];
    }

    float value(unsigned int v) const {
        return hist_[v];
    }
    
    // if float is float 
    void normalization() {
        float s = sum();
        for(int i = 0; i < hist_.size(); i++) {
            hist_[i] = hist_[i] / s;
        }
    }

    float sum() {
        float s = 0;
        for(int i = 0; i < hist_.size(); i++) {
            s = s + hist_[i];
        }
        return s;
    }

    void equaliseImage(Image& img) {
        float sumh = 0.0;
        unsigned int table[img.scale+1];

        for(int i = 0; i <= img.scale; i++) {
            if ( hist_[i] < 1) {    // is normalized
                sumh = sumh + hist_[i] * (img.width()*img.height());       
            } else {
                sumh = sumh + hist_[i];
            }
            table[i] = (unsigned int)( img.scale * 1.0 /(img.width()*img.height()) * sumh + 0.0001 );
        }

        for(int y = 0; y < img.height(); y++) {
            for (int x = 0; x < img.width(); x++) {
                 img.data(x,y) = table[ img.data(x,y) ];
            }
        }
    }

    unsigned int autoThreshold() {
        normalization();

        // get the expection 
        float e = 0.0;
        for(int i = 0; i < hist_.size(); i++) {
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
        for(int i = 1; i < hist_.size(); i++) {
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
