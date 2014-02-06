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
        implements CameraView.CameraReadyCallback
{
    public static String TAG="BV";
    
    private ReentrantLock previewLock = new ReentrantLock();
    private CameraView cameraView = null;
    private OverlayView overlayView = null;
    private AudioRecord audioCapture = null;

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

        // init audio and camera
        initAudio();
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
        cameraView.StartPreview();
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
        //overlayView_.setUpdateDoneCallback(this);
    }
    private void initAudio() {
        int minBufferSize = AudioRecord.getMinBufferSize(16000, AudioFormat.CHANNEL_IN_MONO, AudioFormat.ENCODING_PCM_16BIT);
        int minTargetSize = 1600 * 2;      // 0.2 seconds buffer size
        if (minTargetSize < minBufferSize) {
            minTargetSize = minBufferSize;
        }
        if (audioCapture == null) {
            audioCapture = new AudioRecord(MediaRecorder.AudioSource.MIC,
                    16000,
                    AudioFormat.CHANNEL_IN_MONO,
                    AudioFormat.ENCODING_PCM_16BIT,
                    minTargetSize);
        }
        
        audioCapture.startRecording();
        AudioProcessor audioProcessor = new AudioProcessor();
        audioProcessor.start();  
    }   
     
    //
    //  Internal help class and object definment
    //
    private PreviewCallback previewCb = new PreviewCallback() {
        public void onPreviewFrame(byte[] frame, Camera c) {
            previewLock.lock(); 
            //nativeAgent.updatePicture(frame, cameraView.Width(), cameraView.Height(), is420);
            c.addCallbackBuffer(frame);
            previewLock.unlock();
        }
    };

    private OnClickListener controlAction = new OnClickListener() {
        @Override
        public void onClick(View v) {
                
        }   
    };

    private class AudioProcessor extends Thread {
        byte[] audioPackage = new byte[1024*16];
        int packageSize = 320*6;  // 60ms

        @Override
        public void run() {
            while(true) {
                int ret = audioCapture.read(audioPackage, 0, packageSize);
                if ( ret == AudioRecord.ERROR_INVALID_OPERATION ||
                     ret == AudioRecord.ERROR_BAD_VALUE) {
                    break; 
                }
                
                //nativeAgent.updatePCM(audioPackage, ret);  
            }
        }   
    }
    
}
