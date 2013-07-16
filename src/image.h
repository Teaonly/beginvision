#ifndef _BV_IMAGE_H_
#define _BV_IMAGE_H_

#include <Eigen/Core>

namespace bv {

class Image {
public:
    Image(int width, int height) {
        for(int i = 0; i < 3; i++) {
            color_[i] = new Eigen::MatrixXi(width, height);
        }
    }
    Image(int w1, int h1, int w2, int h2, int w3, int h3) {
        color_[0] = new Eigen::MatrixXi(w1,h1);
        color_[1] = new Eigen::MatrixXi(w2,h2);
        color_[2] = new Eigen::MatrixXi(w3,h3);
    }
    
    ~Image() {
        for(int i = 0; i < 3; i++) {
            if ( color_[i] != NULL) {
                delete color_[i];
            }    
        }
    }

    Eigen::MatrixXi& color(int i) {
        return *(color_[i]);
    }

private:
    Eigen::MatrixXi* color_[3];
};

}
#endif
