LOCAL_PATH:= $(call my-dir)
MY_LOCAL_PATH = $(LOCAL_PATH)

###########################################################
# building application library 
#
include $(CLEAR_VARS)
LOCAL_MODULE := libbeginvision
LOCAL_CPP_EXTENSION := .cc .cpp
LOCAL_CPPFLAGS := -O2 -Werror -Wall -Wno-switch -Wno-non-virtual-dtor -Wno-ctor-dtor-privacy -fno-rtti -fpic -fno-exceptions 
LOCAL_CPPFLAGS += -DPOSIX -DLINUX -DANDROID -DARCH_CPU_LITTLE_ENDIAN 
LOCAL_C_INCLUDES :=  $(MY_LOCAL_PATH)
include $(MY_LOCAL_PATH)/buildme.mk


LOCAL_SHARED_LIBRARIES := libcutils\
                          libgnustl\
                          libdl 

LOCAL_LDLIBS := -llog

include $(BUILD_SHARED_LIBRARY)

