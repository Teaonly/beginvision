#ifndef _BV_IMAGE_H_
#define _BV_IMAGE_H_

#include <Eigen/Core>

namespace bv {

class Image {
public:
    Image(int width, int height) : data(width, height) {
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
};

template<int ColorNumber>
class ColorImage {
public:
    ColorImage(int width, int height) {
        for(int i = 0; i < ColorNumber; i++) {
            color_[i] = new Image(width, height);
        }
    }
    ColorImage(const int w[], const int h[]) {
        for (int i = 0; i < ColorNumber; i++)    {
            color_[i] = new Image(w, h);
        }
    }
    

    ~ColorImage() {
        for(int i = 0; i < ColorNumber; i++) {
            if ( color_[i] != NULL) {
                delete color_[i];
            }    
        }
    }
    
    Image& color(int c) {
        return *(color_[c]);
    }

private:
    Image* color_[ColorNumber];
};

}
#endif
