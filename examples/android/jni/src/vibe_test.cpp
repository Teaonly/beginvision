#include <bv/motion_vibe.h>
#include "helper.h"
#include "vibe_test.h"

static bv::MD_ViBE*    detector_ = NULL;

int VibeInit() {
    
    return 0;    
}

int VibeUpdateForResult(const unsigned char* frameIn, 
                        unsigned char* frameOut,
                        unsigned int wid, unsigned int hei ) {
    if ( detector_ == NULL) {
        detector_ = new bv::MD_ViBE(wid, hei);
    }
    if ( detector_ == NULL) {
        return -1;
    }

    bv::Image inImage(wid, hei);
    bv::Image outImage(wid, hei);
    for(int y = 0; y < (int)hei; y++) {
        for(int x = 0; x < (int)wid; x++) {
            inImage.data(x,y) = frameIn[wid*y + x];
        }
    }

    int ret = detector_->run(inImage, outImage);

    for(int y = 0; y < (int)hei; y++) {
        for(int x = 0; x < (int)wid; x++) {
            frameOut[wid*y + x] = outImage.data(x,y);
        }
    }
    
    LOGD(">>>>>>>>>>>>%d", ret); 

    return 0;
}

