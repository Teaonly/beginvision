#ifndef _BV_DESCRIPTOR_SIFT_H_
#define _BV_DESCRIPTOR_SIFT_H_ 

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

#include "detector_sift.h"

namespace bv {

class DS_Sift {
public:   
    class SiftDescriptor {
    public:    
        SiftDescriptor(int nbp, int nbo) {
            for ( int i = 0; i < nbo; i++) {
                Eigen::MatrixXd empty = Eigen::MatrixXd::Zero(nbp, nbp);
                values_.push_back(empty);    
            }    
        }

        std::vector<Eigen::MatrixXd> values_;
    }; 

   
    DS_Sift(DT_Sift& dt):detector_(dt) {
        NBP_ = 4;
        NBO_ = 8;
        winFactor_ = 3.0;
        run();
    }

private:
    void run() {
        std::vector<DT_Sift::SiftImageOctave>& octaves(detector_.octaves_);
        double weightSigma = NBP_ / 2.0;
        
        for(int ki = 0; ki < detector_.keyPoints_.size(); ki++) {  
            SiftDescriptor desc(NBP_, NBO_);                    
            DT_Sift::SiftKeyPoint n = detector_.keyPoints_[ki]; 
            
            int cx = floor( n.xx_ + 0.5);
            int cy = floor( n.yy_ + 0.5);
            double s = n.ss_;
            int si = floor( n.ss_ + 0.5);
            if ( si > detector_.numLevels_ - 1) {
                si = detector_.numLevels_ - 1;
            }
            int oi = n.octaveIndex_;
            double mainAngle = n.angle_; 

            int width = octaves[oi].images_[si].rows();
            int height = octaves[oi].images_[si].cols();
            double binSize = winFactor_ * detector_.sigma0_ * powf(2, s/detector_.S_ );
            int windowSize = floor( binSize * sqrt(2) * (NBP_+1) / 2.0 + 0.5);
            
            // Define the are of descriptor's window
            int leftX = Util::max(cx-windowSize, 0);
            int rightX = Util::min(cx+windowSize, width-1); 
            int topY = Util::max(cy-windowSize, 0);
            int bottomY = Util::min(cy+windowSize, height-1);
            for (int x = leftX; x <= rightX; x++) {
                for ( int y = topY; y <= bottomY; y++) {

                    double xx = x - n.xx_;
                    double yy = y - n.yy_;
                    double dx = (cos(mainAngle) * xx + sin(mainAngle) * yy) / binSize;
                    double dy = (-1*sin(mainAngle) * xx + cos(mainAngle) * yy) / binSize;
                    double weight = exp( (dx*dx+dy*dy) / (2.0*weightSigma*weightSigma) );
                    double mag =   octaves[oi].gradsX_[si](x,y) * octaves[oi].gradsX_[si](x,y) 
                                 + octaves[oi].gradsY_[si](x,y) * octaves[oi].gradsY_[si](x,y);
                    mag = sqrt(mag);
                    double angle = atan2( octaves[oi].gradsY_[si](x,y), octaves[oi].gradsX_[si](x,y) );
                    if ( angle < 0) {  
                        angle = angle + 2*PI;
                    }
                    angle = 1.0 * NBO_ * angle / (2.0*PI) ; 
                    int xbin = floor(dx - 0.5);
                    int ybin = floor(dy - 0.5);
                    int abin = floor(angle);  
                    
                    for ( int ibin = 0; ibin < 2; ibin ++) {
                        for ( int jbin = 0; jbin < 2; jbin++) {
                            for ( int kbin = 0; kbin < 2; kbin++) {
                                // check int descripting area sround by window
                                if ( (xbin + ibin) >= -(NBP_/2) &&
                                     (xbin + ibin) < (NBP_/2) &&
                                     (ybin + jbin) >= -(NBP_/2) &&
                                     (ybin + jbin) < (NBP_/2) ) {
                                    desc.values_[ (kbin + abin) % NBO_](xbin + ibin+NBP_/2, ybin+jbin+NBP_/2) += 
                                            mag * weight *
                                            (1-abs((xbin+ibin+0.5) - dx)) *
                                            (1-abs((ybin+jbin+0.5) - dy)) *
                                            (1-abs(kbin+abin-angle)) ;
                                }
                            }
                        }
                    }
               }
            }
            
            // Normalizes in norm L_2 a descriptor
            double l2sum = 0.0;
            for(int i = 0; i < NBO_; i++) {
                for ( int j = 0; j < NBP_; j++) {
                    for ( int k = 0; k < NBP_; k++) {
                        if ( desc.values_[i](j,k) > 0.2) {
                            desc.values_[i](j,k) = 0.2;
                        }
                        l2sum += desc.values_[i](j,k) * desc.values_[i](j,k);
                    } 
                }
            }
            l2sum = sqrt(l2sum);
            for(int i = 0; i < NBO_; i++) {
                for ( int j = 0; j < NBP_; j++) {
                    for ( int k = 0; k < NBP_; k++) {
                        desc.values_[i](j,k) = desc.values_[i](j,k) / l2sum;
                    } 
                }
            }
                       
            descs_.push_back( desc);
        }
    }    
    
private:
    DT_Sift& detector_;   
    double winFactor_; 
    int NBP_;
    int NBO_;

    std::vector<SiftDescriptor> descs_;
};

}
#endif
