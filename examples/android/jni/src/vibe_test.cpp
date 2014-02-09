#include <android/bitmap.h>
#include <bv/motion_vibe.h>
#include "helper.h"
#include "vibe_test.h"

static bv::MD_ViBE*    detector_ = NULL;

int VibeInit() {
    
    return 0;    
}

int VibeUpdateForResult(JNIEnv* env,
                        const unsigned char* frameIn, 
                        jobject bitmap,
                        unsigned int wid, unsigned int hei ) {
    if ( detector_ == NULL) {
        detector_ = new bv::MD_ViBE(wid/2, hei/2);
    }
    if ( detector_ == NULL) {
        return -1;
    }

    bv::Image inImage(wid/2, hei/2);
    bv::Image outImage(wid/2, hei/2);

    for(int y = 0; y < inImage.height(); y++) {
        for(int x = 0; x < inImage.width(); x++) {
            inImage.data(x,y) = 0;
            for(int yy = y*2; yy <= y*2 + 1; yy++) {
                for(int xx = x*2; xx <= x*2 + 1; xx++) {
                    inImage.data(x,y) += frameIn[wid*yy + xx];
                }
            }
            inImage.data(x,y) = inImage.data(x,y) / 4;
        }
    }
    
    int ret = detector_->run(inImage, outImage);
    LOGD(">>>>>>>>> Result = %d", ret); 

    
    AndroidBitmapInfo  info; 
    unsigned int*              pixels;
    if ((ret = AndroidBitmap_getInfo(env, bitmap, &info)) < 0) {
        LOGD("AndroidBitmap_getInfo() failed ! error=%d", ret);
        return -1;
    }
    if ((ret = AndroidBitmap_lockPixels(env, bitmap, (void**)&pixels)) < 0) {
        LOGD("AndroidBitmap_lockPixels() failed ! error=%d", ret);
        return -1;
    }
   
    int lineStride = info.stride / 4;
    for(int y = 0; y < (int)hei; y++) {
        for(int x = 0; x < (int)wid; x++) {
            if ( outImage.data(x/2, y/2) ) {
                pixels[y*lineStride+x] = 0xFFFFFFFF;
            } else {
                pixels[y*lineStride+x] = 0x00000000;
            }
        }
    }
 
    AndroidBitmap_unlockPixels(env, bitmap);
    return 0;
}

