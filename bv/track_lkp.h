#ifndef _BV_TRACK_LKP_H_
#define _BV_TRACK_LKP_H_ 

#include <vector>
#include <Eigen/Core>
#include <Eigen/LU> 
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
        minImageSize_(minImageSize),
        winR_(winR) {
    }

public:
    int run(Eigen::MatrixXd& img1, Eigen::MatrixXd& img2, const std::vector<double>& sourceX, const std::vector<double>& sourceY, std::vector<double>& destX, std::vector<double>& destY ) {
        // Check the input and out size  
        if ( (img1.rows () != img2.rows()) || (img1.cols() != img2.cols()) ) {
            return BV_ERROR_PARAMETER;
        }
        if ( sourceX.size() != sourceY.size() || sourceX.size() <= 0) {
            return BV_ERROR_PARAMETER;
        }
        destX.resize( sourceX.size() );
        destY.resize( sourceY.size() );
        
         // Building internal image pyramid 
        int wid = img1.rows();
        int hei = img1.cols();
        int minsize = Util::min( wid, hei );
        int pyrLevel = (int) (log(minsize * 1.0 / minImageSize_) / log(2.0));
        buildPyr(img1, imagePyr1_, pyrLevel);
        buildPyr(img2, imagePyr2_, pyrLevel);
        
        for(int i = 0; i < (int)sourceX.size(); i++) {
            double x,y;
            _LucasKanade(sourceX[i], sourceY[i], x, y);
            destX[i] = x;
            destY[i] = y;            
        }
        
        return BV_OK;
    }

private:
    void buildPyr(Eigen::MatrixXd& img, std::vector<Eigen::MatrixXd>& imgPyr, int level) {
        imgPyr.clear();
        imgPyr.resize(level);
        imgPyr[0] = img;      // bottom (largest) layer
        
        Eigen::MatrixXd ker = Kernel::gaussian(2, 1);        
        //Filter::withTemplate(img, imgPyr[0], ker);

        for(int i = 1; i < level; i++) {
            Eigen::MatrixXd d0 = imgPyr[i-1];
            Filter::withTemplate(imgPyr[i-1], d0, ker);
            Eigen::MatrixXd d1 = Eigen::MatrixXd(d0.rows()/2, d0.cols()/2);
            for(int c = 0; c < d1.cols(); c++) {
                for (int r = 0; r < d1.rows(); r++) {
                    d1(r,c) = d0(r*2, c*2);
                }
            }
            
            imgPyr.push_back(d1);
        }
    }

    void _LucasKanade(const double sx, const double sy, double& dx, double& dy) {
        double gx = 0.0;
        double gy = 0.0;
        
        double xt = sx;
        double yt = sy;
        for(int i = 0; i < (int)imagePyr1_.size(); i++) {
            xt = xt / 2;
            yt = yt / 2;
        }
        
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
        
        for (int l = imagePyr1_.size() - 1; l >= 0; l-- ) {
            xt = xt * 2;
            yt = yt * 2;
 
            Eigen::MatrixXd& img1(imagePyr1_[l]);
            Eigen::MatrixXd& img2(imagePyr2_[l]);
            Eigen::MatrixXd img1xd(img1.rows(), img1.cols());
            Eigen::MatrixXd img1yd(img1.rows(), img2.cols());
            Filter::withTemplate(img1, img1xd, h);
            Filter::withTemplate(img1, img1yd, ht);
            
            // fetch source image patch from img1 
            Eigen::MatrixXd ps( 2*winR_ + 1, 2*winR_ + 1);
            Eigen::MatrixXd px( 2*winR_ + 1, 2*winR_ + 1);
            Eigen::MatrixXd py( 2*winR_ + 1, 2*winR_ + 1);
            fetchImage(xt, yt, img1, ps);
            fetchImage(xt, yt, img1xd, px);
            fetchImage(xt, yt, img1yd, py);
                        
            Eigen::MatrixXd G = Eigen::MatrixXd::Zero(2,2);
            for(int i = 0; i < ps.cols(); i++) {
                for ( int j = 0; j < ps.rows(); j++) {
                    G(0,0) = px(j,i) * px(j,i) + G(0,0);
                    G(0,1) = px(j,i) * py(j, i) + G(0,1);
                    G(1,1) = py(j,i) * py(j,i) + G(1,1);
                }
            }
            G(1,0) = G(0,1);

            Eigen::MatrixXd pd( 2*winR_ + 1, 2*winR_ + 1);
            double vx = 0.0;
            double vy = 0.0;
            for ( int i = 0; i < maxIter_; i++) {
                fetchImage(xt + gx + vx, yt + gy + vy, img2, pd);
                pd = ps - pd;

                Eigen::VectorXd b(2); 
                b(0) = 0.0;
                b(1) = 0.0;
                
                double dsum = 0.0;
                for(int i = 0; i < ps.cols(); i++) {
                    for ( int j = 0; j < ps.rows(); j++) {
                        b(0) = pd(j,i) * px(j,i) + b(0);
                        b(1) = pd(j,i) * py(j,i) + b(1);
                        dsum = dsum + pd(j,i)*pd(j,i);
                    }   
                }
                
                Eigen::VectorXd delta(2);
                delta = G.inverse() * b;
                
                vx = vx + delta(0);
                vy = vy + delta(1);
                
                if ( (delta(0)*delta(0) + delta(1)*delta(1)) <= (threshold_ * threshold_) ) {
                    break;    
                }
            }

            gx = gx + vx;
            gy = gy + vy;
        }
    
        dx = gx;
        dy = gy;
    }
    
    bool fetchImage(double x, double y, Eigen::MatrixXd& img, Eigen::MatrixXd& patch) {
        if ( x < winR_ || x >= img.rows() - winR_ || y < winR_ || y >= img.cols() - winR_ ) {
            return false;
        }
        for ( int yy = -1*winR_; yy <= winR_; yy++) {
            for ( int xx = -1*winR_; xx <= winR_; xx++) {
                double r = Util::interp2(x + xx, y + yy, img);
                patch(xx + winR_, yy + winR_) = r;
            }
        }
        
        return true; 
    }

protected:
    int    maxIter_;
    double threshold_;
    int    minImageSize_;
    int    winR_;

private:
    std::vector<Eigen::MatrixXd> imagePyr1_;
    std::vector<Eigen::MatrixXd> imagePyr2_;
};

}
#endif
