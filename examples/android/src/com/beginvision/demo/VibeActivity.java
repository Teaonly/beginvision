package com.beginvision.demo;

import java.util.*;
import java.util.concurrent.locks.ReentrantLock;
import org.apache.http.conn.util.InetAddressUtils;
import android.app.Activity;
import android.content.Context;
import android.content.Intent;
import android.content.SharedPreferences;
import android.hardware.Camera;
import android.hardware.Camera.PreviewCallback;
import android.hardware.Camera.PictureCallback;
import android.graphics.Bitmap;
import android.media.AudioFormat;
import android.media.MediaRecorder;
import android.media.AudioRecord;
import android.os.Bundle;
import android.os.Looper;
import android.os.Handler;
import android.util.*;
import android.view.SurfaceView;
import android.view.View;
import android.view.View.OnClickListener;
import android.view.Window;
import android.view.WindowManager;
import android.widget.Button;
import android.widget.TextView;

public class VibeActivity extends Activity
        implements CameraView.CameraReadyCallback, OverlayView.UpdateDoneCallback
{
    public static String TAG="BV";
    
    private ReentrantLock previewLock = new ReentrantLock();
    private CameraView cameraView = null;
    private OverlayView overlayView = null;
    private byte[]  resultFrame = null;
    private Bitmap  resultBitmap = null;    

    //
    //  Activiity's event handler
    //
    @Override
    public void onCreate(Bundle savedInstanceState) {
        // application setting
        requestWindowFeature(Window.FEATURE_NO_TITLE);
        Window win = getWindow();
        win.addFlags(WindowManager.LayoutParams.FLAG_KEEP_SCREEN_ON);    

        // load and setup GUI
        super.onCreate(savedInstanceState);
        setContentView(R.layout.vibe);

        TextView tv = (TextView)findViewById(R.id.tv_message);
        tv.setText("This a VIBE (VIsual Background Extractor) demo");
        Button btn = (Button)findViewById(R.id.btn_control);
        btn.setOnClickListener(controlAction);

        // init camera
        initCamera();
   }
    @Override
    public void onDestroy() {
        super.onDestroy();
    }

    @Override
    public void onPause() {      
        super.onPause();

        if ( cameraView != null) {
            previewLock.lock();
            cameraView.StopPreview();
            previewLock.unlock();
        }

        //finish();
        System.exit(0);
    }

    @Override
    public void onBackPressed() {
        super.onBackPressed();
    }

    //
    //  Interface implementation
    //
    public void onCameraReady() {
        cameraView.StopPreview();
        cameraView.setupCamera(640, 480, 4, previewCb);
        resultFrame = new byte[cameraView.Width() * cameraView.Height()]; 
        resultBitmap = Bitmap.createBitmap(640, 480, Bitmap.Config.ARGB_8888);
        cameraView.StartPreview();
    }

    public void onUpdateDone() {
                
    }

    //
    //  Internal help functions
    //
    private void initCamera() {
        SurfaceView cameraSurface = (SurfaceView)findViewById(R.id.surface_camera);
        cameraView = new CameraView(cameraSurface);        
        cameraView.setCameraReadyCallback(this);

        overlayView = (OverlayView)findViewById(R.id.surface_overlay);
        //overlayView_.setOnTouchListener(this);
        overlayView.setUpdateDoneCallback(this);
    }
     
    //
    //  Internal help class and object definment
    //
    private PreviewCallback previewCb = new PreviewCallback() {
        public void onPreviewFrame(byte[] yuvFrame, Camera c) {
            processNewFrame(yuvFrame, c);
        }
    };

    private void processNewFrame(final byte[] yuvFrame, final Camera c) {
        if ( previewLock.isLocked() ) {
            c.addCallbackBuffer(yuvFrame);
        }
        
        new Thread(new Runnable() {
                    public void run() {
                        previewLock.lock(); 
                        NativeAgent.updatePictureForResult("VIBE", yuvFrame, resultFrame, cameraView.Width(), cameraView.Height());
                        c.addCallbackBuffer(yuvFrame);
                        new Handler(Looper.getMainLooper()).post( resultAction );
                        previewLock.unlock();
                    }
                }).start();
    }

    private OnClickListener controlAction = new OnClickListener() {
        @Override
        public void onClick(View v) {
                
        }   
    };
    
    private Runnable resultAction = new Runnable() {
        @Override 
        public void run() {
            //overlayView.DrawResult(resultBitmap);
        }
    };

}
