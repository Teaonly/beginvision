#include <android/bitmap.h>
#include <bv/motion_vibe.h>
#include "helper.h"
#include "vibe_test.h"

static bv::MD_ViBE*    detector_ = NULL;

static const int SCALE = 3;

int VibeInit() {
    
    return 0;    
}

int VibeUpdateForResult(JNIEnv* env,
                        const unsigned char* frameIn, 
                        jobject bitmap,
                        unsigned int wid, unsigned int hei ) {
    if ( detector_ == NULL) {
        detector_ = new bv::MD_ViBE(wid/SCALE, hei/SCALE);
    }
    if ( detector_ == NULL) {
        return -1;
    }

    bv::Image inImage(wid/SCALE + 1, hei/SCALE + 1);
    bv::Image outImage(wid/SCALE + 1, hei/SCALE + 1);
    
    inImage.data *= 0;
    for(int y = 0; y < (int)hei; y++) {
        for(int x = 0; x < (int)wid; x++) {
            int xx = x/SCALE;
            int yy = y/SCALE; 
            inImage.data(xx, yy) = frameIn[y*wid+x] + inImage.data(xx, yy);
        }
    }
    inImage.data /= SCALE*SCALE;

    int ret;
    ret = detector_->run(inImage, outImage);

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
            int xx = x/SCALE;
            int yy = y/SCALE; 
            if ( outImage.data(xx, yy) ) {
                pixels[y*lineStride+x] = 0xFFFFFFFF;
            } else {
                pixels[y*lineStride+x] = 0x00000000;
            }
        }
    }
 
    AndroidBitmap_unlockPixels(env, bitmap);
    return 0;
}

