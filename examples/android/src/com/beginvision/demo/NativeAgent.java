package com.beginvision.demo;

import java.io.*; 
import java.net.*;

import android.net.*;
import android.util.Log;

public class NativeAgent{
    public static native int updatePicture(String target, byte[]frame, int wid, int hei);

    public static native int updatePictureForResult(String target, 
                                                    byte[]frame, Object obj, int wid, int hei);
                                                    
}
