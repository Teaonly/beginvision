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
        SiftDescriptor(int nbp, int nbo):nbp_(nbp), nbo_(nbo) {
            for ( int i = 0; i < nbo_*nbp_*nbp_; i++) {
                values_.push_back(0.0);    
            }    
        }
        double& val(int x, int y, int o) {
            return values_[x*nbo_ + y*nbp_*nbo_ + o]; 
        }
        
        double distanceWith(const SiftDescriptor& b) {
            if ( nbp_ != b.nbp_  ||
                 nbo_ != b.nbo_ )  {
                return -1;
            }
            double l2dist = 0.0;
            for (int i = 0; i < values_.size(); i++) {
                l2dist += (values_[i] - b.values_[i] ) * (values_[i] - b.values_[i] );
            }
            return sqrt(l2dist);
        }

        std::vector<double> values_;
        int nbp_;
        int nbo_;
    }; 

    DS_Sift(DT_Sift& dt):detector_(dt) {
        NBP_ = 4;
        NBO_ = 8;
        winFactor_ = 3.0;
        threshold_ = 0.8;
        run();
    }

    int matchWith(DS_Sift& other, std::vector<int>& results) {
        if ( other.NBP_ != NBP_ ||
             other.NBO_ != NBO_ ) {
            return BV_ERROR_PARAMETER;
        }
        
        int count = 0;
        results.clear();
        for(int i = 0; i < descs_.size(); i++) {
            double minValue[2];
            int minIndex;  

            minValue[0] = 1 / threshold_;
            minValue[1] = 1 / threshold_;
            for ( int j = 0; j < other.descs_.size(); j++) {
                double dist = descs_[i].distanceWith(other.descs_[j] );
                if ( dist <  minValue[0] ) {
                    minValue[1] = minValue[0];

                    minIndex = j;
                    minValue[0] = dist;
                } else if ( dist < minValue[1] ) {
                    minValue[1] = dist;
                }
            }
             
            if ( minValue[0] / minValue[1] < threshold_ ) {
                count ++;
                results.push_back(minIndex);
            } else {
                results.push_back(-1);
            }
        }   
        
        std::cout << "Matched points = " << count << std::endl;

        return BV_OK;
    }

    int saveMatch(DS_Sift& other, std::vector<int>& results , const std::string& fileName) {
        FILE* fp = fopen(fileName.c_str(), "wt");

        for (int i = 0; i < results.size(); i++) {
            if ( results[i] != -1) {
                DT_Sift::SiftKeyPoint a = detector_.keyPoints_[i];
                DT_Sift::SiftKeyPoint b = other.detector_.keyPoints_[ results[i] ];
                
                double xa, ya, xb, yb;
                xa = a.xx_ * powf(2, a.octaveIndex_ + detector_.minOctave_); 
                ya = a.yy_ * powf(2, a.octaveIndex_ + detector_.minOctave_);

                xb = b.xx_ * powf(2, b.octaveIndex_ + other.detector_.minOctave_); 
                yb = b.yy_ * powf(2, b.octaveIndex_ + other.detector_.minOctave_);
                
                fprintf(fp, "%f %f %f %f\n", xa,ya,xb,yb);
            }
        }
        
        fclose(fp); 
    } 

    // Just for debug
    int showMatch(DS_Sift& other, std::vector<int>& results, Eigen::MatrixXd& img1, Eigen::MatrixXd& img2) {
        std::cout << "Show result..." << std::endl;
        for (int i = 0; i < results.size(); i++) {
            if ( results[i] != -1) {
                detector_.showKeypoint(img1, i);
                int j = results[i];
                other.detector_.showKeypoint(img2, j);
            }
        }
        Util::saveAsImage(img1, "/tmp/xxx1.bmp");
        Util::saveAsImage(img2, "/tmp/xxx2.bmp");
        
        return BV_OK;
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
            int leftX = Util::max(cx-windowSize, 1);
            int rightX = Util::min(cx+windowSize, width-2); 
            int topY = Util::max(cy-windowSize, 1);
            int bottomY = Util::min(cy+windowSize, height-2);
            for ( int y = topY; y <= bottomY; y++) {
                for (int x = leftX; x <= rightX; x++) {

                    double xx = x - n.xx_;
                    double yy = y - n.yy_;
                    double dx = (cos(mainAngle) * xx + sin(mainAngle) * yy) / binSize;
                    double dy = (-1*sin(mainAngle) * xx + cos(mainAngle) * yy) / binSize;
                    double weight = exp( (dx*dx+dy*dy) / (2.0*weightSigma*weightSigma) );
                    double mag =   octaves[oi].gradsX_[si](x,y) * octaves[oi].gradsX_[si](x,y) 
                                 + octaves[oi].gradsY_[si](x,y) * octaves[oi].gradsY_[si](x,y);
                    mag = sqrt(mag);
                    double angle = atan2( octaves[oi].gradsY_[si](x,y), octaves[oi].gradsX_[si](x,y) ) - mainAngle;
                    angle = Util::mod2pi(angle);

                    angle = 1.0 * NBO_ * angle / (2.0*PI) ; 
                    int xbin = floor(dx - 0.5);
                    int ybin = floor(dy - 0.5);
                    int abin = floor(angle);  
                    
                    double rxbin = dx - (xbin + 0.5) ; 
                    double rybin = dy - (ybin + 0.5) ;
                    double rabin = angle - abin ;            

                    for ( int ibin = 0; ibin < 2; ibin ++) {
                        for ( int jbin = 0; jbin < 2; jbin++) {
                            for ( int kbin = 0; kbin < 2; kbin++) {
                                // check int descripting area sround by window
                                if ( (xbin + ibin) >= -(NBP_/2) &&
                                     (xbin + ibin) < (NBP_/2) &&
                                     (ybin + jbin) >= -(NBP_/2) &&
                                     (ybin + jbin) < (NBP_/2) ) {
                                    int o = (kbin + abin)%NBO_;
                                    int nx = xbin + ibin + NBP_/2;
                                    int ny = ybin + jbin + NBP_/2;
                                    desc.val(nx, ny, o) += 
                                            mag * weight *
                                            fabs(1.0 - 1.0*ibin - rxbin) *
                                            fabs(1.0 - 1.0*jbin - rybin) *
                                            fabs(1.0 - 1.0*kbin - rabin);
                                }
                            }
                        }
                    }
               }
            }
        
            // Normalizes in norm L_2 a descriptor
            double l2sum = 0.0;
            for(int i = 0; i < NBO_* NBP_ * NBP_; i++) {
                l2sum += desc.values_[i] * desc.values_[i];
            }
            l2sum = sqrt(l2sum);
            for(int i = 0; i < NBO_* NBP_ * NBP_; i++) {
                desc.values_[i] /= l2sum;
            }
            l2sum = 0.0;
            for(int i = 0; i < NBO_* NBP_ * NBP_; i++) {
                if ( desc.values_[i] > 0.2) {
                    desc.values_[i] = 0.2;
                }
                l2sum += desc.values_[i] * desc.values_[i];
            }
            l2sum = sqrt(l2sum);
            for(int i = 0; i < NBO_* NBP_ * NBP_; i++) {
                desc.values_[i] /= l2sum;
            }

            descs_.push_back( desc);
/*
            for(int i = 0; i < NBO_* NBP_ * NBP_; i++) {
                std::cout << desc.values_[i] << " " ;
            }
            std::cout << std::endl;
            exit(0);
*/
        }
    }    
    

private:
    DT_Sift& detector_;   
    double winFactor_; 
    int NBP_;
    int NBO_;
    double threshold_;

    std::vector<SiftDescriptor> descs_;
};

}
#endif
