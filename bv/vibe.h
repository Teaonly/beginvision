#ifndef _BV_DETECTOR_HARRIS_H_
#define _BV_DETECTOR_HARRIS_H_ 

#include <vector>
#include <Eigen/Core>
#include <Eigen/LU> 
#include <cmath>
#include "beginvision.h"
#include "util.h"

namespace bv {

class MD_ViBE {
public:
    MD_ViBE(int wid, int hei):
        sampleSize_(16) {
        
        
    }
    
public:
    int run(Image& in, Image& out) {
        
        return BV_OK;
    }

protected:
    std::vector<Eigen::MatrixXi> bgModle_;     //Backgournd data modle
    const unsigned int sampleSize_;
};


}


#endif
