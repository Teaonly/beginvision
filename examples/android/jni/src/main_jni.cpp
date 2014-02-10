#include <string>
#include <jni.h>

#include "vibe_test.h"

#undef JNIEXPORT
#define JNIEXPORT __attribute__((visibility("default")))
#define JOW(rettype, name) extern "C" rettype JNIEXPORT JNICALL \
      Java_com_beginvision_demo_NativeAgent_##name

//
//  Global variables
//


//
//  Internal helper functions
//
static std::string convert_jstring(JNIEnv* env, const jstring &js) {
    static char outbuf[1024];
    int len = env->GetStringLength(js);
    env->GetStringUTFRegion(js, 0, len, outbuf);
    std::string str = outbuf;
    return str;
}

/*
static jint get_native_fd(JNIEnv* env, jobject fdesc) {
    jclass clazz;
    jfieldID fid;

    if (!(clazz = env->GetObjectClass(fdesc)) ||
            !(fid = env->GetFieldID(clazz,"descriptor","I"))) return -1;
    return env->GetIntField(fdesc,fid);
}
*/

//
//  Global functions called from Java side 
//
JOW(int, init)(JNIEnv* env, jclass, jstring target) {
    std::string objTarget = convert_jstring(env, target);
    
    if ( objTarget == "VIBE" ) {
        VibeInit();        
    }
    
    return 0;
}

JOW(int, updatePictureForResult)(JNIEnv* env, jclass, jstring target, jbyteArray yuvData, jobject result, jint wid, jint hei) {
    std::string objTarget = convert_jstring(env, target);

    jbyte* cameraFrame = env->GetByteArrayElements(yuvData, NULL);
    if ( objTarget == "VIBE") {
        VibeUpdateForResult(env, (unsigned char*)cameraFrame, result, wid, hei);
    }
    env->ReleaseByteArrayElements(yuvData, cameraFrame, 0);
    return 0;
}

JOW(int, updatePicture)(JNIEnv* env, jclass, jstring target, jbyteArray yuvData, jint wid, jint hei) {
    std::string objTarget = convert_jstring(env, target);

    jbyte* cameraFrame= env->GetByteArrayElements(yuvData, NULL);
    
    env->ReleaseByteArrayElements(yuvData, cameraFrame, 0);
    return 0;
}


//
//  Top level library load/unload 
//
extern "C" jint JNIEXPORT JNICALL JNI_OnLoad(JavaVM *jvm, void *reserved) {
    JNIEnv* jni;
    if (jvm->GetEnv(reinterpret_cast<void**>(&jni), JNI_VERSION_1_6) != JNI_OK)
        return -1;
    return JNI_VERSION_1_6;
}

extern "C" jint JNIEXPORT JNICALL JNI_OnUnLoad(JavaVM *jvm, void *reserved) {
    return 0;
}

