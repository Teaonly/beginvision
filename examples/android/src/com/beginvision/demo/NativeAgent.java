package com.beginvision.demo;

import java.io.*; 
import java.net.*;

import android.net.*;
import android.util.Log;

public class NativeAgent{
    public static native int updatePicture(byte[]frame, int wid, int hei);
}
