#ifndef _BV_DETECTOR_SIFT_H_
#define _BV_DETECTOR_SIFT_H_ 

#include <vector>
#include <Eigen/Core>
#include <Eigen/LU> 
#include <cmath>
#include "beginvision.h"
#include "filter.h"
#include "util.h"

namespace bv {

class DT_Sift {
public:
    DT_Sift( ) {
    }
    
public:
    int run(Eigen::MatrixXd& img) {
        // normlized imag
        Eigen::MatrixXd I = img;
        //I = I - I.minCoeff();
        //I = I / I.maxCoeff();    

        return BV_OK;
    }

protected:
    
};

}


#endif
