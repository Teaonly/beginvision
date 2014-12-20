#ifndef _BV_DESCRIPTOR_HOG_H_
#define _BV_DESCRIPTOR_HOG_H_ 

#include <vector>
#include <utility>
#include <iterator>
#include <Eigen/Core>
#include <Eigen/LU> 
#include <cmath>
#include "beginvision.h"
#include "filter.h"
#include "image_convert.h"
#include "util.h"

namespace bv {

class DS_Hog {
public:   
    DS_Hog(int bin=36):_num_bin_(bin) {
        
    }

public:
    int run(Eigen::MatrixXd& img) {
                 
    }

private:
    unsigned int _num_bin;
};

}
#endif
