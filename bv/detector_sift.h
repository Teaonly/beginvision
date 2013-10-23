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
    DT_Sift(int numOctaves = 5, int S = 3, int minOctave = 0): 
            numOctaves_(numOctaves), minOctave_(minOctave), S_(S) { 
        
        numLevels_ = S_ + 3;
        k = powf(2, 1/S_);
        sigmaNominal_ = 0.5;
        sigma0_ = 1.6; 
        dsigma_ = sqrt(powf(2, 2.0/S_) - 1); 

        peakThreshold_ = 0.003;
        edgeThreshold_ = (10 + 1)*(10 + 1) / 10.0;
        winFactor_ = 1.5;
        orientHistNumber_ = 36;
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
        doRefine();
        
        // 4. Caculating the orientation 
        buildGrad();
        getOrientation();

        // 5. show detecing result , just for debug
        //I = img;
        //showDetect(I);
          
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
                siftSmooth(lastLevel, temp, diffSigma);
                octave.images_.push_back(temp);
                
                lastLevel = temp;
            }
            if ( oi != (numOctaves_ - 1) ) {
                // prepare for next octave
                bottomLevel.resize( bottomLevel.rows()/2, bottomLevel.cols()/2);
                
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

                //int boundary = (int)(sigma0_*4.0) * si;
                int boundary = 1;
    
                for(int y = boundary; y < octaves_[oi].height_ - boundary; y++) {
                    for(int x = boundary; x < octaves_[oi].width_ - boundary; x++) {
                        double centerValue = octaves_[oi].dogs_[middle](x, y);
                        
                        if ( fabs(centerValue) < 0.8 * peakThreshold_ ) {
                            continue;
                        }

                        bool isMin = true;
                        bool isMax = true;
                        for(int xx = x-1; xx <= x+1; xx++) {
                            for(int yy = y-1; yy <= y+1; yy++) {
                                if ( xx != x || yy != y) {
                                    if (    octaves_[oi].dogs_[up](xx,yy) <= centerValue 
                                         || octaves_[oi].dogs_[middle](xx,yy) <= centerValue 
                                         || octaves_[oi].dogs_[down](xx,yy) <= centerValue ) {
                                        isMin = false;    
                                    } 

                                    if (    octaves_[oi].dogs_[up](xx,yy) >= centerValue 
                                         || octaves_[oi].dogs_[middle](xx,yy) >= centerValue 
                                         || octaves_[oi].dogs_[down](xx,yy) >= centerValue ) {
                                        isMax = false;    
                                    } 

                                } else {
                                    if (    octaves_[oi].dogs_[up](xx,yy) <= centerValue 
                                         || octaves_[oi].dogs_[down](xx,yy) <= centerValue ) {
                                        isMin = false;    
                                    } 

                                    if (    octaves_[oi].dogs_[up](xx,yy) >= centerValue 
                                         || octaves_[oi].dogs_[down](xx,yy) >= centerValue ) {
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
            std::cout << "Local Extream(" << oi << "): Key points number is " << keyPoints_.size() << std::endl;
        }
    }

#define AT(o,x,y,s) octaves_[(o)].dogs_[(s)]((x),(y))     
    void doRefine() {
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
            
            for ( int i = 0; i < 5; i++) { 
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
                
                if ( fabs(A.determinant()) > 0) {
                    c = b * A.inverse(); 
                } else {
                    c(0,0) = 0;
                    c(0,1) = 0;
                    c(0,2) = 0;
                }

                dx = dy = 0; 
                if ( c(0,0) > 0.6 && x < octaves_[o].width_ - 2) {
                    dx = 1;
                } else if ( c(0,0) < -0.6 && x > 1) {
                    dx = -1;
                }
                if ( c(0,1) > 0.6 && y < octaves_[o].height_ - 2) {
                    dy = 1;
                } else if ( c(0,1) < -0.6 && y > 1) {
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
            
            double xx = x + c(0,0);
            double yy = y + c(0,1);
            double ss = s + c(0,2);            
            if (    (fabs(refineValue) > peakThreshold_) 
                 && (score < edgeThreshold_)
                 && (score >= 0)
                 && (fabs(c(0,0) ) < 1.5)
                 && (fabs(c(0,1) ) < 1.5)
                 && (fabs(c(0,2) ) < 1.5) 
                 && (xx > 0)
                 && (xx < octaves_[o].width_-1)
                 && (yy > 0)
                 && (yy < octaves_[o].height_-1) 
                 && (ss >= 0) 
                 && (ss < numLevels_ )  ) {

                (*n).xx_ = xx;
                (*n).yy_ = yy;
                (*n).ss_ = ss;
                n++;
            } else {
                // throw it out
                n = keyPoints_.erase(n);
            }
        }
        std::cout << "After refine: Key points number is " << keyPoints_.size() << std::endl;
    }
    
    void buildGrad() {
        Eigen::MatrixXd hx(3,1);
        hx(0,0) = -0.5;
        hx(1,0) = 0.0;
        hx(2,0) = 0.5;
        Eigen::MatrixXd hy = hx.transpose();
        for(int oi = 0; oi < numOctaves_; oi++) {
            for ( int li = 0; li < numLevels_; li++) {
                int wid = octaves_[oi].images_[li].rows();
                int hei = octaves_[oi].images_[li].cols();

                Eigen::MatrixXd gradX(wid, hei);
                Eigen::MatrixXd gradY(wid, hei);

                Filter::withTemplate(octaves_[oi].images_[li], gradX, hx, Filter::EXTENTION_NEAR);
                Filter::withTemplate(octaves_[oi].images_[li], gradY, hy, Filter::EXTENTION_NEAR);
                
                octaves_[oi].gradsX_.push_back(gradX);
                octaves_[oi].gradsY_.push_back(gradY);   
            }
        }
    }   

    void getOrientation() {
        std::vector<SiftKeyPoint> newPoints;
        for ( std::vector<SiftKeyPoint>::iterator n = keyPoints_.begin(); n != keyPoints_.end(); n++ ) {
            int cx = floor( (*n).xx_ + 0.5);
            int cy = floor( (*n).yy_ + 0.5);
            double s = (*n).ss_;
            int oi = (*n).octaveIndex_;
            
            int si = floor( (*n).ss_ + 0.5);
            if ( si > numLevels_ - 1) {
                si = numLevels_ - 1;
            }
            int width = octaves_[oi].images_[si].rows();
            int height = octaves_[oi].images_[si].cols();

            double weightSigma = winFactor_ * sigma0_ * powf(2, s/S_ );
            int windowSize = floor( 3.0 * weightSigma ); 
 
          
            // caculating the histogram of grad's angle 
            Eigen::ArrayXd  hist(orientHistNumber_);
            for(int i = 0; i < orientHistNumber_; i++) {
                hist(i) = 0.0;
            }
            int leftX = Util::max(cx-windowSize, 0);
            int rightX = Util::min(cx+windowSize, width-1); 
            int topY = Util::max(cy-windowSize, 0);
            int bottomY = Util::min(cy+windowSize, height-1);
            
            for ( int y = topY; y <= bottomY; y++) {
                for (int x = leftX; x <= rightX; x++) {
                    /* limit to a circular window */ 
                    double r2 = (1.0*x-(*n).xx_)*(1.0*x-(*n).xx_) + (1.0*y-(*n).yy_)*(1.0*y-(*n).yy_);
                    if (r2 >= windowSize*windowSize + 0.6) continue ; 
                    
                    double mag =   octaves_[oi].gradsX_[si](x,y) * octaves_[oi].gradsX_[si](x,y) 
                                 + octaves_[oi].gradsY_[si](x,y) * octaves_[oi].gradsY_[si](x,y);
                    mag = sqrt(mag);

                    double weight = exp( 1.0 * r2 / (2*weightSigma*weightSigma) );
                    double angle = atan2(octaves_[oi].gradsY_[si](x,y), octaves_[oi].gradsX_[si](x,y));
                    if ( angle < 0) {
                        angle = angle + 2*PI;
                    }
                    double fbin = orientHistNumber_*angle/(2*PI);
                    int bin = (int) floor(fbin - 0.5) ; 
                    double rbin = fbin - bin - 0.5 ;
                    hist((bin + orientHistNumber_) % orientHistNumber_) += (1 - rbin) * mag * weight ;
                    hist((bin +                 1) % orientHistNumber_) += (    rbin) * mag * weight ;     
                }
            }
            

            /*  smooth histogram */ 
            for (int iter = 0; iter < 6; iter ++) {
                double prev  = hist (orientHistNumber_ - 1) ; 
                double first = hist (0) ;
                int i ;
                for (i = 0; i < orientHistNumber_ - 1; i++) {
                    double newh = (prev + hist(i) + hist((i+1) % orientHistNumber_)) / 3.0; 
                    prev = hist(i) ;
                    hist(i) = newh ;
                }    
                hist(i) = (prev + hist(i) + first) / 3.0 ;
            }
            
           // find the maxvalue 
            int maxIndex;
            double maxV = hist.maxCoeff(&maxIndex);
            (*n).angle_ = refineAngle(hist, maxIndex);

            for(int i = 0; i < orientHistNumber_; i++) {
                if ( i == maxIndex ) {
                    continue;
                }           
                int left = (i-1+orientHistNumber_) % orientHistNumber_;      
                int right = (i+1) % orientHistNumber_;
                if ( hist[i] >= 0.8 * maxV && hist[i] > hist[left] && hist[i] > hist[right] ) {
                    SiftKeyPoint np = *n;
                    np.angle_ = refineAngle(hist, i);            
                    newPoints.push_back(np);
                }
            }
        }
        keyPoints_.insert(keyPoints_.begin(), newPoints.begin(), newPoints.end());
        
        std::cout << "After orientation : " << keyPoints_.size() << std::endl;
    }
    
    void showKeypoint(Eigen::MatrixXd& img, int i) {
        SiftKeyPoint key = keyPoints_[i];

        int centerx = (int)key.xx_ * powf(2, key.octaveIndex_ + minOctave_); 
        int centery = (int)key.yy_ * powf(2, key.octaveIndex_ + minOctave_);
        
        double scale = powf(2, key.octaveIndex_ + minOctave_) * sigma0_;
        scale = scale * powf(2, key.ss_/S_); 
        //scale = scale * 2;

        for ( double r = 0.0; r <= 2*PI ; r += PI/40) {
            int xx = (int)( sinf(r) * scale + centerx);
            int yy = (int)( cosf(r) * scale + centery);
            if (    xx >= 0 
                 && xx < img.rows() 
                 && yy >= 0
                 && yy < img.cols() ) {
                img(xx,yy) = 1;
            }
        }
        double angle_ = key.angle_;
        for (int d = 0; d < scale; d++) {
            int xx = (int)( sinf(angle_) * d + centerx);
            int yy = (int)( cosf(angle_) * d + centery);
            if (    xx >= 0 
                 && xx < img.rows() 
                 && yy >= 0
                 && yy < img.cols() ) {
                img(xx,yy) = 1;
            }
        }
    }

    void showDetect(Eigen::MatrixXd& img) {
        for (int i = 0; i < keyPoints_.size(); i++) {
            if ( keyPoints_[i].octaveIndex_ >= 1) {
                showKeypoint(img, i);
            }
        }
        Util::saveAsImage(img, "/tmp/xxx.bmp");
    }

    double refineAngle(const Eigen::ArrayXd& hist, const int& peakIndex) {
        int leftIndex = (peakIndex - 1 + hist.size()) % hist.size();
        int rightIndex = (peakIndex + 1) % hist.size();
        double leftValue = hist[leftIndex];
        double rightValue = hist[rightIndex];
        double peakValue = hist[peakIndex];
        
#if 0        
        Eigen::MatrixXd A(3,3);
        Eigen::MatrixXd C(3,1);
        Eigen::MatrixXd B(3,1);
        A(0,0) = 1; A(1,0) = -1; A(2,0) = 1;
        A(0,1) = 0; A(1,1) = 0;  A(2,1) = 1;
        A(0,2) = 1; A(1,0) = 1; A(2,2) = 1;

        C(0,0) = leftValue;
        C(1,0) = peakValue;
        C(2,0) = rightValue;
       
        B = A.inverse() * C;
        double diff = -1 * B(1,0) / (2 * B(0,0));
        double refineIndex = peakIndex + diff;
        if ( refineIndex < 0) {
            refineIndex += hist.size();
        } else if ( refineIndex >= hist.size() ) {
            refineIndex -= hist.size(); 
        }
        
        return 2.0 * PI * refineIndex / hist.size();
#else
        double diff = - 0.5 * (rightValue - leftValue) / ( rightValue + leftValue - 2 * peakValue) ;
        double refineIndex = peakIndex + diff;
        double angle = 2.0 * PI * (refineIndex + 0.5) / hist.size();
        if ( angle > 2*PI) {
            angle = angle - 2*PI;
        }
        return angle;
#endif

            
    }

    void siftSmooth(Eigen::MatrixXd& in, Eigen::MatrixXd& out, double sigma) {
        int kerWidth = ceil(sigma*4.0);
        kerWidth = Util::max( kerWidth, 1); 
        
        //Eigen::MatrixXd ker = Kernel::gaussian(kerWidth, sigma);   
        //Filter::withTemplate(in, out, ker, Filter::EXTENTION_ZERO);
        Filter::gaussianBlur(in, out, kerWidth*2+1, sigma, Filter::EXTENTION_ZERO);
    }


private:
    typedef struct {
        unsigned int x_;
        unsigned int y_;
        unsigned int levelIndex_;
        unsigned int octaveIndex_;

        // refined value;
        double xx_;
        double yy_;
        double ss_;
        double angle_;
    } SiftKeyPoint;

    typedef struct {
        unsigned int width_;
        unsigned int height_;
        std::vector<Eigen::MatrixXd> images_;
        std::vector<Eigen::MatrixXd> dogs_;
        std::vector<Eigen::MatrixXd> gradsX_;
        std::vector<Eigen::MatrixXd> gradsY_;
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
    double edgeThreshold_;
    double winFactor_;
    int orientHistNumber_;    

    std::vector<SiftImageOctave> octaves_;
    std::vector<SiftKeyPoint> keyPoints_;
    friend class DS_Sift;
};

}


#endif
