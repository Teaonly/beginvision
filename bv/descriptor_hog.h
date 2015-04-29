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
    
    double value(unsigned int v) const {
        return hist_[v];
    }
    
    std::vector<double> getHitogram() {
        return hist_;
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
    
    void prepareImage(Eigen::MatrixXd& img) {
        // sobel operator
        Eigen::MatrixXd h(3,3);
        h(0,0) = -1;
        h(1,0) = 0;
        h(2,0) = 1;
        
        h(0,1) = -2;
        h(1,1) = 0;
        h(2,1) = 2;
        
        h(0,2) = -1;
        h(1,2) = 0;
        h(2,2) = 1;
        Eigen::MatrixXd ht = h.transpose();
        
        imgXd_.resize(img.rows(), img.cols());
        imgYd_.resize(img.rows(), img.cols());
        Filter::withTemplate(img, imgXd_, h);
        Filter::withTemplate(img, imgYd_, ht);
    }
    
    int extrac(int bx, int by, int wid, int hei, bool normlized = true) {
        if ( imgXd_.rows() == 0 || imgXd_.cols() == 0) {
            return BV_ERROR;
        }
        
        for(int i = 0; i < _num_bin; i++) {
            hist_[i] = 0.0;
        }
        
        for (int x = bx+1; x < bx+wid - 1; x++) {
            for ( int y = by+1; y < by+hei - 1; y++) {
                double mag = sqrt(pow(imgXd_(x,y), 2) + pow(imgYd_(x,y), 2));
                double orient = atan2(imgYd_(x,y), imgXd_(x,y) );
                unsigned int bin = _binFor(orient);
                hist_[bin] += mag;
            }
        }
        
        if ( normlized == true) {
            // L1 normlize
            double sum = 0.00000001;
            for(unsigned int i = 0; i < _num_bin; i++) {
                sum += hist_[i];
            }
            for(unsigned int i = 0; i < _num_bin; i++) {
                hist_[i] /= sum;
            }
        }

        return BV_OK;
    }
    
    int extrac(Eigen::MatrixXd& img, bool normlized = true) {
        prepareImage(img);
        return extrac(0, 0, (int)img.rows(), (int)img.cols(), normlized);
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
    Eigen::MatrixXd imgXd_;
    Eigen::MatrixXd imgYd_;

};

}
#endif
