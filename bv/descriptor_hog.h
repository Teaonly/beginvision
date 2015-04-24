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
    DS_Hog(unsigned int bin=36):_num_bin(bin) {
        for(int i = 0; i < _num_bin; i++) {
            hist_.push_back(0.0);
        }
    }
    
    unsigned int size() {
        return (unsigned int)hist_.size();
    }
    
    float value(unsigned int v) const {
        return hist_[v];
    }
    
    float distWith(DS_Hog& that) {
        if ( that._num_bin != _num_bin) {
            return 0.0;
        }
        
        double dist = 0.0;
        for(unsigned int i = 0; i < (unsigned int)hist_.size(); i++) {
            dist += pow(hist_[i] - that.hist_[i], 2);
        }
        dist = sqrt(dist);
        return (float)dist;
    }
    
    int extrac(Eigen::MatrixXd& img) {
        for(int i = 0; i < _num_bin; i++) {
            hist_[i] = 0.0;
        }
        
        // sobel operator
        Eigen::MatrixXd h(3,3);
        h(0,0) = -1;
        h(1,0) = 0;
        h(2,0) = 1;
        
        h(0,1) = -1;
        h(1,1) = 0;
        h(2,1) = 1;
        
        h(0,2) = -1;
        h(1,2) = 0;
        h(2,2) = 1;
        Eigen::MatrixXd ht = h.transpose();

        Eigen::MatrixXd imgXd(img.rows(), img.cols());
        Eigen::MatrixXd imgYd(img.rows(), img.cols());
        Filter::withTemplate(img, imgXd, h);
        Filter::withTemplate(img, imgYd, ht);

        for (int x = 1; x < img.rows() - 1; x++) {
            for ( int y = 1; y < img.cols() - 1; y++) {
                double mag = sqrt(pow(imgXd(x,y), 2) + pow(imgYd(x,y), 2));
                double orient = atan2(imgYd(x,y), imgXd(x,y) );
                unsigned int bin = _binFor(orient);
                hist_[bin] += mag;
            }
        }

        // L1 normlize
        double sum = 0.00000001;
        for(unsigned int i = 0; i < _num_bin; i++) {
            sum += hist_[i];
        }
        for(unsigned int i = 0; i < _num_bin; i++) {
            hist_[i] /= sum;
        }

        return BV_OK;
    }

private:
    unsigned int _binFor(double radians) {
        radians += bv::PI;
        unsigned int bin = floor(_num_bin * radians / (2*bv::PI));
        return bin;
    }

private:
    unsigned int _num_bin;
    std::vector<double> hist_;

};

}
#endif
