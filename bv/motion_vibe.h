#ifndef _BV_MD_VIBE_H_
#define _BV_MD_VIBE_H_ 

#include <vector>
#include <Eigen/Core>
#include <Eigen/LU> 
#include <cmath>
#include "beginvision.h"
#include "util.h"

namespace bv {

class MD_ViBE {
public:
    MD_ViBE(int wid, int hei):
        isInited_(false),
        sampleNumber_(20),       
        grayDistance_(12),
        minSelected_(2),
        randomSub_(16) {
        Eigen::MatrixXi* onePicture;
        for(unsigned int i = 0; i < sampleNumber_; i++) {
            onePicture = new Eigen::MatrixXi(wid, hei);
            bgModle_.push_back(onePicture);    
        }
    }

    ~MD_ViBE() {
        for(unsigned int i = 0; i < sampleNumber_; i++) {
            delete bgModle_[i];
        }
    }
    
public:
    int run(Image& in, Image& out) {
        if ( isInited_ == false) {
            initWithFirstImage(in); 
            isInited_ = true;
        }
        
        int wid = bgModle_[0]->rows();
        int hei = bgModle_[0]->cols();
        int count = 0;
        for (int y = 0; y < hei; y++) {
            for (int x = 0; x < wid; x++) {
                int bg = checkPixel(x, y, in.data(x,y));
                updateModle(x, y, in.data(x,y), bg);
                out.data(x,y) = 1 - bg;
                count += bg;
            }
        } 

        return count;
    }

private:
    void initWithFirstImage(Image& firstPicture) {
        int wid = bgModle_[0]->rows();
        int hei = bgModle_[0]->cols();

        for (int y = 0; y < hei; y++) {
            for (int x = 0; x < wid; x++) {
                for(int n = 0; n < (int)sampleNumber_; n++) {
                    int x2 = x;
                    int y2 = y;
                    getRandomNeighbr(x2, y2);
                    (*bgModle_[n])(x,y) = firstPicture.data(x2,y2);
                }
            }
        }
    }
    
    void getRandomNeighbr(int& x, int& y) {
        if ( isInited_ == false)
            return;

        int r = 0;
        int xr, yr;
        do {
            r = random() % 3 - 1;        
            xr = x + r;
            r = random() % 3 - 1;
            yr = y + r;
        } while(   xr < 0 || xr >= bgModle_[0]->rows() 
                || yr < 0 || yr >= bgModle_[0]->cols() );
        x = xr;
        y = yr;
    }

    int checkPixel(int x, int y, int myValue) {
        unsigned int count = 0;
        for (int n = 0; n < (int)sampleNumber_; n++) {
            if ( abs(myValue - (*bgModle_[n])(x,y)) <= (int)grayDistance_) {
                count++;
            }
            if ( count >= minSelected_) {
                break;
            }
        }
        return count >= minSelected_;
    }
    
    void updateModle(int x, int y, int newData, int bg) {
        if (bg == 0) {
            return;
        }
        if ( random() % randomSub_ == 0) {
            int r = random() % sampleNumber_;
            (*bgModle_[r])(x,y) = newData;
        }
        if ( random() % randomSub_ == 0) {
            getRandomNeighbr(x, y);
            int r = random() % sampleNumber_;
            (*bgModle_[r])(x,y) = newData;
        }
    }

protected:
    std::vector<Eigen::MatrixXi*> bgModle_;     //Backgournd data modle
    bool isInited_;
    const unsigned int sampleNumber_;
    const unsigned int grayDistance_;          
    const unsigned int minSelected_;
    const unsigned int randomSub_;
};

}
#endif
