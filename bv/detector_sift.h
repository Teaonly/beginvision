#ifndef _BV_DETECTOR_SIFT_H_
#define _BV_DETECTOR_SIFT_H_ 

#include <vector>
#include <Eigen/Core>
#include <Eigen/LU> 
#include <cmath>
#include "beginvision.h"
#include "filter.h"
#include "image_convert.h"
#include "util.h"


namespace bv {

class DT_Sift {
public:
    DT_Sift(bool initDoubled = false): 
        initDoubled_(initDoubled)  {
    }
    
public:
    int run(Eigen::MatrixXd& img) {
        // 0. normlized image
        Eigen::MatrixXd I = img;
        I = I - Eigen::MatrixXd::Ones(img.rows(), img.cols()) * I.minCoeff();
        I = I / I.maxCoeff();    
            
        // 1. double original image
        if ( initDoubled_ ) {
            Eigen::MatrixXd doubleI ( img.rows(), img.cols() );
            Convert::resizeImage(I, doubleI);
            I = doubleI;
        }

        // 2. building scale-space image
         

        
        return BV_OK;
    }

protected:
    bool initDoubled_;
};

}


#endif
