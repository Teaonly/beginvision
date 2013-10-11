#ifndef _BV_DETECTOR_SIFT_H_
#define _BV_DETECTOR_SIFT_H_ 

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

class DT_Sift {
public:
    DT_Sift(int numOctaves = 4, int S = 3, int minOctave = 0): 
            numOctaves_(numOctaves), minOctave_(minOctave), S_(S) { 
        
        numLevels_ = S_ + 3;
        k = powf(2, 1/S_);
        sigmaNominal_ = 0.5;
        sigma0_ = 1.6; 
        dsigma_ = sqrt(powf(2, 2.0/S_) - 1); 

        peakThreshold_ = 0.003;
        diffThreshold_ = 0.0001;
        edgeThreshold_ = (10 + 1)*(10 + 1) / 10.0;
    }
    
public:
    int run(Eigen::MatrixXd& img) {
        
        // 0. normlized image
        Eigen::MatrixXd I = img;

        I = I - Eigen::MatrixXd::Ones(img.rows(), img.cols()) * I.minCoeff();
        I = I / I.maxCoeff();    
            
        // 1. upsamle/downsample the original image
        if ( minOctave_ != 0) {
            // TODO 
        }     
    

        // 2. building scale-space image
        buildOctaves(I);
        buildDoG();        
        
        // 3. detect maxima and minima of difference-of-Gaussian in scale space
        doDetect();
        //refineDetect();
        
        // 4. show detecing result , just for debug
        I = img;
        showDetect(I);
          
        return BV_OK;
    }
    
private:
    /*
       The scale space is defined as:
        sigma(o,s) = sigma0 * pow(2, o + s / S) = sigma0 * 2^o * 2^(s/S) 
        s is from smin to smin+numLevle

       I_sigman = g_sqrt( sigman^2 - sigma0^2) * I_sigma0  ( sigman > sigma0) 

       So the diffrent scale step is computed as follow:
       
            sqrt( sigma(o,s+1)^2 - sigma(o,s)^2) =
       =    sqrt( sigma0^2 * 2^(2o) * ( 2^((2s+2)/S) - 2^(2s/S) ))
       =    sigma0 * 2^o * sqrt ( 2^(2s/S) ( 2^(2/S) - 1) )
       =    sigma0 * 2^(o+s/S) * sqrt( 2^(2/S) - 1)

       we can ommit the 2^o in the same octvae.  
    */
    void buildOctaves( Eigen::MatrixXd& I) {
        Eigen::MatrixXd bottomLevel = I;
        double bottomSigma = sigma0_ * powf(2, minOctave_);
        octaves_.clear();

        if ( bottomSigma > sigmaNominal_ * powf(2, minOctave_)) {
            double sa = bottomSigma;
            double sb = sigmaNominal_ * powf(2, minOctave_); 
            double sigma = sqrt( sa*sa - sb*sb);
            
            Eigen::MatrixXd temp = bottomLevel;
            std::cout << " bottom sigma = " << bottomSigma << std::endl;
            std::cout << " first  sigma = " << sigma << std::endl;
            siftSmooth(bottomLevel, temp, sigma);
            bottomLevel = temp;
        }

        for(int oi = 0; oi < numOctaves_; oi++) {
            SiftImageOctave octave;
            octave.width_ = bottomLevel.rows();
            octave.height_ = bottomLevel.cols();
            octave.images_.push_back( bottomLevel );
            Eigen::MatrixXd lastLevel = bottomLevel;
            
            for ( int li = 1; li < numLevels_; li++) {

                Eigen::MatrixXd temp = lastLevel;
                
                double diffSigma = sigma0_ * dsigma_ * powf(2, (li-1)*1.0/S_);
                std::cout << " sigma = " << diffSigma << std::endl;
                siftSmooth(lastLevel, temp, diffSigma);
                octave.images_.push_back(temp);
                
                lastLevel = temp;
            }
            if ( oi != (numOctaves_ - 1) ) {
                // prepare for next octave
                bottomLevel.resize( bottomLevel.rows()/2, bottomLevel.cols()/2);
                // TODO using downsample replacing resize.
                Convert::resizeImage( octave.images_[S_], bottomLevel);
            }
            octaves_.push_back(octave);
        }    
    }

    void buildDoG() {
         for(int oi = 0; oi < numOctaves_; oi++) {
            for ( int li = 0; li < numLevels_ - 1; li++) {
                Eigen::MatrixXd dog = octaves_[oi].images_[li+1] - octaves_[oi].images_[li];
                octaves_[oi].dogs_.push_back(dog);
            }
         }
    }

    void doDetect() {
        keyPoints_.clear();

        for(int oi = 0; oi < numOctaves_; oi++) {
            for ( int si = 1; si < octaves_[oi].dogs_.size() - 1; si++) {
                int up = si+1; 
                int middle = si;
                int down = si-1;

                int boundary = (int)(sigma0_ * 6);

                for(int x = boundary; x < octaves_[oi].width_ - boundary; x++) {
                    for(int y = boundary; y < octaves_[oi].height_ - boundary; y++) {
                        double centerValue = octaves_[oi].dogs_[middle](x, y);
                        
                        /*
                        if ( fabs(centerValue) < peakThreshold_ ) {
                            continue;
                        }
                        */

                        bool isMin = true;
                        bool isMax = true;
                        for(int xx = x-1; xx <= x+1; xx++) {
                            for(int yy = y-1; yy <= y+1; yy++) {
                                if ( xx != x || yy != y) {
                                    if (    octaves_[oi].dogs_[up](xx,yy) - diffThreshold_ <= centerValue 
                                         || octaves_[oi].dogs_[middle](xx,yy) - diffThreshold_ <= centerValue 
                                         || octaves_[oi].dogs_[down](xx,yy) - diffThreshold_ <= centerValue ) {
                                        isMin = false;    
                                    } 

                                    if (    octaves_[oi].dogs_[up](xx,yy) + diffThreshold_ >= centerValue 
                                         || octaves_[oi].dogs_[middle](xx,yy) + diffThreshold_ >= centerValue 
                                         || octaves_[oi].dogs_[down](xx,yy) + diffThreshold_ >= centerValue ) {
                                        isMax = false;    
                                    } 

                                } else {
                                    if (    octaves_[oi].dogs_[up](xx,yy) - diffThreshold_ <= centerValue 
                                         || octaves_[oi].dogs_[down](xx,yy) - diffThreshold_ <= centerValue ) {
                                        isMin = false;    
                                    } 

                                    if (    octaves_[oi].dogs_[up](xx,yy) + diffThreshold_ >= centerValue 
                                         || octaves_[oi].dogs_[down](xx,yy) + diffThreshold_ >= centerValue ) {
                                        isMax = false;    
                                    }                                            
                                }

                                if ( isMin == false && isMax == false) {
                                    goto _detect_done;
                                }
                            }
                        }
_detect_done:
                        if ( isMax || isMin) {
                            SiftKeyPoint newKey;
                            newKey.x_ = x;
                            newKey.y_ = y;
                            newKey.levelIndex_ = si;
                            newKey.octaveIndex_ = oi;
                            newKey.xx_ = x;
                            newKey.yy_ = y;
                            newKey.ss_ = si;
                            keyPoints_.push_back(newKey);                            
                        }
                    }
                }
            }
        }
    
        std::cout << "Local Extream: Key points number is " << keyPoints_.size() << std::endl;
    }

#define AT(o,x,y,s) octaves_[(o)].dogs_[(s)]((x),(y))     
    void refineDetect() {
        for ( std::vector<SiftKeyPoint>::iterator n = keyPoints_.begin(); n != keyPoints_.end();  ) {
            // 3D quadratic refine the keypoint's location 
            int o = (*n).octaveIndex_;
            int x = (*n).x_;
            int y = (*n).y_;
            int s = (*n).levelIndex_;
            
            int dx, dy;
            double Dx, Dy, Ds, Dxx, Dyy, Dss, Dxy, Dxs, Dys;
            Eigen::MatrixXd A(3,3);
            Eigen::MatrixXd b(1,3);
            Eigen::MatrixXd c(1,3);
            
            for ( int i = 0; i < 3; i++) { 
                Dx = 0.5 * ( AT(o, x+1, y,   s  ) - AT(o, x-1, y,   s  ) );
                Dy = 0.5 * ( AT(o, x,   y+1, s  ) - AT(o, x,   y-1, s  ) );
                Ds = 0.5 * ( AT(o, x,   y,   s+1) - AT(o, x,   y,   s-1) );
                
                Dxx =  AT(o, x+1, y,   s  ) + AT(o, x-1, y,   s  ) - 2*AT(o, x, y, s);
                Dyy =  AT(o, x,   y+1, s  ) + AT(o, x,   y-1, s  ) - 2*AT(o, x, y, s);
                Dss =  AT(o, x,   y,   s+1) + AT(o, x,   y,   s-1) - 2*AT(o, x, y, s);
                
                Dxy = 0.25 * (   AT(o, x+1, y+1, s  ) + AT(o, x-1, y-1, s  ) 
                               - AT(o, x-1, y+1, s  ) - AT(o, x+1, y-1, s) ) ;            
                Dxs = 0.25 * (   AT(o, x+1, y,   s+1) + AT(o, x-1, y,   s-1) 
                               - AT(o, x+1, y,   s-1) - AT(o, x-1, y,   s+1) ) ;
                Dys = 0.25 * (   AT(o, x,   y+1, s+1) + AT(o, x,   y-1, s-1) 
                               - AT(o, x,   y+1, s-1) - AT(o, x,   y-1, s+1) ) ;

                A(0,0) = Dxx;
                A(1,1) = Dyy;
                A(2,2) = Dss;
                A(0,1) = A(1,0) = Dxy;
                A(0,2) = A(2,0) = Dxs;
                A(1,2) = A(2,1) = Dys;
                b(0,0) = -1 * Dx;    
                b(0,1) = -1 * Dy;
                b(0,2) = -1 * Ds;

                c = b * A.inverse();
                
                dx = dy = 0; 
                if ( c(0,0) > 0.6) {
                    dx = 1;
                } else if ( c(0,0) < -0.6) {
                    dx = -1;
                }
                if ( c(0,1) > 0.6) {
                    dy = 1;
                } else if ( c(0,1) < 0.6) {
                    dy = -1;
                }
                
                if ( dx == 0 && dy == 0) {
                    break;
                } else {       
                    x += dx;
                    y += dy;
                }
            }
            
            double score = (Dxx+Dyy)*(Dxx+Dyy) / (Dxx*Dyy - Dxy*Dxy) ; 
            double refineValue = AT(o,x,y,s) + 0.5 * (Dx * c(0, 0) + Dy * c(0, 1) + Ds * c(0, 2)) ;   
            
            if (    (fabs(refineValue) > peakThreshold_) 
                 && (score < edgeThreshold_)
                 && (score > 0)
                 && (fabs(c(0,0) ) < 1.5)
                 && (fabs(c(0,1) ) < 1.5)
                 && (fabs(c(0,2) ) < 1.5) ) {

                (*n).xx_ = x + c(0,0);
                (*n).yy_ = y + c(0,1);
                (*n).ss_ = s + c(0,2);
                
                n++;
            } else {
                // throw it out
                n = keyPoints_.erase(n);
            }
        }
        std::cout << "After refine: Key points number is " << keyPoints_.size() << std::endl;
    }

    void siftSmooth(Eigen::MatrixXd& in, Eigen::MatrixXd& out, double sigma) {
        //Eigen::MatrixXd ker = Kernel::gaussian((int)(sigma*2+0.5), sigma);   
        //Filter::withTemplate(in, out, ker, Filter::EXTENTION_ZERO);
        Filter::gaussianBlur(in, out, (int)(sigma*2+0.5)*2+1, sigma);
    }
    
    void showDetect(Eigen::MatrixXd& img) {
        for (int i = 0; i < keyPoints_.size(); i++) {
            if ( keyPoints_[i].octaveIndex_ >= 2) {
                int centerx = (int)keyPoints_[i].xx_ * powf(2, keyPoints_[i].octaveIndex_ + minOctave_); 
                int centery = (int)keyPoints_[i].yy_ * powf(2, keyPoints_[i].octaveIndex_ + minOctave_);
                
                double scale = powf(2, keyPoints_[i].octaveIndex_ + minOctave_) * sigma0_;
                scale = scale * powf(2, keyPoints_[i].levelIndex_/S_); 

                for ( double r = 0.0; r <= 2*pi ; r += pi/40) {
                    int xx = (int)( sinf(r) * scale + centerx);
                    int yy = (int)( cosf(r) * scale + centery);
                    if (    xx >= 0 
                         && xx < img.rows() 
                         && yy >= 0
                         && yy < img.cols() ) {
                        img(xx,yy) = 1;
                    }
                }
            }
        }
        Util::saveAsImage(img, "/tmp/xxx.bmp");
    }

public:
    typedef struct {
        unsigned int x_;
        unsigned int y_;
        unsigned int levelIndex_;
        unsigned int octaveIndex_;

        // refined value;
        double xx_;
        double yy_;
        double ss_;
    } SiftKeyPoint;

    typedef struct {
        unsigned int width_;
        unsigned int height_;
        std::vector<Eigen::MatrixXd> images_;
        std::vector<Eigen::MatrixXd> dogs_;
    } SiftImageOctave;

    // passed from parameters
    int numOctaves_;
    int numLevels_;
    int minOctave_;
    int S_;
    
    // initialed in creator
    double k;
    double sigmaNominal_;
    double sigma0_;
    double dsigma_;
    double peakThreshold_;
    double diffThreshold_;
    double edgeThreshold_;

    std::vector<SiftImageOctave> octaves_;
    std::vector<SiftKeyPoint> keyPoints_;
};

}


#endif
