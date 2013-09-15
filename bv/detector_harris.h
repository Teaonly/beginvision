#ifndef _BV_DETECTOR_HARRIS_H_
#define _BV_DETECTOR_HARRIS_H_ 

#include <vector>
#include <Eigen/Core>
#include <Eigen/LU> 
#include <cmath>
#include "beginvision.h"
#include "filter.h"
#include "util.h"

namespace bv {

class DT_Harris {
public:
    DT_Harris(double threshold = 0.0001,
              double alpha = 0.04, 
              int winR = 8):
        threshold_(threshold),
        alpha_(alpha),
        winR_(winR) {
         
    }
    
public:
    int run(Eigen::MatrixXd& img, std::vector<int>& resultX, std::vector<int>& resultY) {
        resultX.clear();
        resultY.clear();
        
        // Step.1 smooth image by Gaussian 
        Eigen::MatrixXd ker = Kernel::gaussian(4, 1.5);   
        Eigen::MatrixXd sImg(img.rows(), img.cols());
        Filter::withTemplate(img, sImg, ker);

        // Step.2 Computing second moment matrix 
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

        Eigen::MatrixXd imgXd(sImg.rows(), sImg.cols());
        Eigen::MatrixXd imgYd(sImg.rows(), sImg.cols());
        Filter::withTemplate(sImg, imgXd, h);
        Filter::withTemplate(sImg, imgYd, ht);
        
        // Step.3 Computing the response for every pixel 
        Eigen::MatrixXd win = Kernel::gaussian(winR_, winR_/2);
        Eigen::MatrixXd rep = Eigen::MatrixXd::Zero(sImg.rows(), sImg.cols());
        for (int x = winR_; x < sImg.rows() - winR_; x++) {
            for ( int y = winR_; y < sImg.cols() - winR_; y++) {

                Eigen::MatrixXd Hs = Eigen::MatrixXd::Zero(2,2);
                for(int xx = x - winR_; xx <=  x+winR_; xx++) {
                    for (int yy = y - winR_; yy <= y+winR_; yy++) {
                        /*
                        Hs(0,0) = imgXd(xx,yy) * imgXd(xx,yy) + Hs(0,0);
                        Hs(1,1) = imgYd(xx,yy) * imgYd(xx,yy) + Hs(1,1);
                        Hs(0,1) = imgXd(xx,yy) * imgYd(xx,yy) + Hs(0,1);
                        */
                        Hs(0,0) = imgXd(xx,yy) * imgXd(xx,yy)*win(xx-x+winR_,yy-y+winR_) + Hs(0,0);
                        Hs(1,1) = imgYd(xx,yy) * imgYd(xx,yy)*win(xx-x+winR_,yy-y+winR_) + Hs(1,1);
                        Hs(0,1) = imgXd(xx,yy) * imgYd(xx,yy)*win(xx-x+winR_,yy-y+winR_) + Hs(0,1);
                    }
                }
                Hs(1,0) = Hs(0,1);
                
                rep(x,y) = Hs.determinant() - alpha_ * Hs.trace() * Hs.trace();
            }
        }
        
        for (int x = 0; x < sImg.rows(); x++) {
            for ( int y = 0; y < sImg.cols(); y++) {
                if ( fabs(rep(x,y)) > threshold_  ) {
                    if ( fabs(rep(x,y)) > fabs(rep(x+1,y)) &&
                         fabs(rep(x,y)) > fabs(rep(x-1,y)) &&  
                         fabs(rep(x,y)) > fabs(rep(x,y+1)) &&
                         fabs(rep(x,y)) > fabs(rep(x,y-1)) &&
                         fabs(rep(x,y)) > fabs(rep(x+1,y-1)) &&
                         fabs(rep(x,y)) > fabs(rep(x-1,y+1)) &&
                         fabs(rep(x,y)) > fabs(rep(x+1,y+1)) &&  
                         fabs(rep(x,y)) > fabs(rep(x-1,y-1)) )
                         {
                        
                        resultX.push_back(x);
                        resultY.push_back(y);                        
                    }  
                }
            }
        }
        
        return BV_OK;
    }
protected:
    double threshold_;
    double alpha_;
    int winR_;
};


}


#endif
