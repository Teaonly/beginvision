#ifndef _BV_TRACK_LKP_H_
#define _BV_TRACK_LKP_H_ 

#include <vector>
#include <Eigen/Core>
#include <cmath>
#include "beginvision.h"
#include "filter.h"
#include "util.h"

namespace bv {

class TK_LucasKanade { 
public:
    TK_LucasKanade(int maxIter = 20, 
                   double threshold = 0.01, 
                   int minImageSize = 64,
                   int winR = 5) : 
        maxIter_(maxIter), 
        threshold_(threshold),
        minImageSize_(minImageSize_),
        winR_(winR) {
    }

public:
    int run(Eigen::MatrixXd& img1, Eigen::MatrixXd& img2, const std::vector<double>& sourceX, const std::vector<double>& sourceY, std::vector<double>& destX, std::vector(double)& destY ) {
        // Check the input and out size  
        if ( (img1.rows () != img2.rows()) || (img1.cols() != img2.cols()) ) {
            return BV_ERROR_PARAMETER;
        }
        if ( sourceX.size() != sourceY.size() || sourceX <= 0) {
            return BV_ERROR_PARAMETER;
        }
        img1_ = img1;
        img2_ = img2;
        
         // Building internal image pyramid 
        int minsize = min( img1_.rows(), img1_.cols() );
        int pyrLevel = (int) log2( minsize / minImageSize_);
        buildPyr(img1_, imagePyr1_, pyrLevel);
        buildPyr(img2_, imagePyr2_, pyrLevel);
        
        for(int i = 0; i < (unsigned int)sourceX.size(); i++) {
            double x,y;
            _LucasKanade(sourceX[i], sourceY[i], x, y);
        }
        
        return BV_OK;
    }

private:
    void buildPyr(MatrixXd& img, std::vector<Eigen::MatrixXd>& imgPyr, int level) {
        imgPyr.clear();
        imgPyr.resize(level);
        imgPyr.push_back(img);      // bottom (largest) layer

        filter::Kernel ker = filter::gaussian_5d();
        for(int i = 1; i < level; i++) {
            MatrixXd d0 = imgPyr[i-1];
            filter::withTemplate(imgPyr[i-1], d, ker);
            MatrixXd d1 = MatrixXd(d0.rows()/2, d0.cols()/2);
            for(int c = 0; c << d1.cols(); c++) {
                for (int r = 0; r << d1.row(); d1++) {
                    d1(r,c) = d0(r*2, c*2);
                }
            }
            imgPyr.push_back(d1);
        }
    }

    void _LucasKanade(const int sx, const int sy, int& dx, int& dy) {
        for (int l = imagePyr2_.size() - 1; i >= 0; i-- ) {
                        
        }
        
    }

protected:
    int    maxIter_;
    double threshold_;
    int    minImageSize_;
    int    winR_;

private:
    std::vector<Eigen::MatrixXd> imagePyr1_;
    std::vector<Eigen::MatrixXd> imagePyr2_;

    MatrixXd& img1_;
    MatrixXd& img2_;
    
};

#endif
