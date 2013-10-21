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
    DS_Sift(DT_Sift& dt):detector_(dt) {
        NBP_ = 4;
        winFactor_ = 3.0;
        run();
    }

private:
    void run() {
        std::vector<DT_Sift::SiftImageOctave>& octaves(detector_.octaves_);
        
        for(int i = 0; i < detector_.keyPoints_.size(); i++) {  
            DT_Sift::SiftKeyPoint n = detector_.keyPoints_[i]; 

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
            double weightSigma = winFactor_ * (NBP_/2) * detector_.sigma0_ * powf(2, s/detector_.S_ );
            int windowSize = floor( weightSigma * sqrt(2) + 0.5);
 
                       
            // Define the are of descriptor's window
            int leftX = Util::max(cx-windowSize, 0);
            int rightX = Util::min(cx+windowSize, width-1); 
            int topY = Util::max(cy-windowSize, 0);
            int bottomY = Util::min(cy+windowSize, height-1);
            for (int x = leftX; x <= rightX; x++) {
                for ( int y = topY; y <= bottomY; y++) {
                    double xx = x - n.xx_;
                    double yy = y - n.yy_;
                    double dx = cos(mainAngle) * xx + sin(mainAngle) * yy;
                    double dy = -1*sin(mainAngle) * xx + cos(mainAngle) * yy;     
                      
                }
            }
        }
    }    
    
private:
    typedef struct {
                
    } SiftDescriptor; 

    DT_Sift& detector_;   
    double winFactor_; 
    int NBP_;

    std::vector<SiftDescriptor> descs_;
};

}
#endif
