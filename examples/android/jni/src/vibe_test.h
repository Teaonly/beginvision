#ifndef _VIBE_TEST_H_
#define _VIBE_TEST_H_

#include <jni.h>

int VibeInit();

int VibeUpdateForResult(JNIEnv* env,
                        const unsigned char* frameIn, 
                        jobject result,
                        unsigned int wid, unsigned int hei );


#endif
