#ifndef _BV_DETECTOR_SIFT_H_
#define _BV_DETECTOR_SIFT_H_ 

#include <vector>
#include <utility>
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
    DT_Sift(int numOctaves = 4, bool initDoubled = false): 
        numOctaves_(numOctaves),
        initDoubled_(initDoubled)  {
        
        K = sqrt(2);       
        startSigma_ = sqrt(2.0)/2;
        levels_ = 5; 
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
        buildOctaves(I); 
        
        return BV_OK;
    }

private:
    void buildOctaves( Eigen::MatrixXd& I) {
        // first Octave
        Eigen::MatrixXd I0 = I;
        Filter::gaussianBlur(I, I0, (int)(startSigma_*3+0.5)*2 + 1, startSigma_);  
        
        SiftImageOctave oct;
        oct.wid_ = I.rows();
        oct.hei_ = I.cols();
        oct.images_.push_back( std::pair<double, Eigen::MatrixXd> (startSigma_, I0) );
        
        Eigen::MatrixXd Ix = I0; 
        double sigma = startSigma_;
        for (int i = 0; i < levels_; i++) {            
            sigma = sigma * K;
            int kerSize = (int)(startSigma_*(K-1)*3+0.5) + 1;
            Filter::gaussianBlur(I0, Ix, kerSize, startSigma_*(K-1)); 
            oct.images_.push_back( std::pair<double, Eigen::MatrixXd> ( sigma, Ix) );
            Ix = I0;            
        } 
        octaves_.push_back(oct);
                
    }

protected:

    typedef struct {
        unsigned int wid_;
        unsigned int hei_;
        std::vector<std::pair<double, Eigen::MatrixXd> > images_;
    } SiftImageOctave;

    int numOctaves_;
    bool initDoubled_;
    double K;
    double startSigma_;
    int    levels_;

    std::vector<SiftImageOctave> octaves_;
};

}


#endif
