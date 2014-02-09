package com.beginvision.demo;

import android.content.Context;
import android.graphics.Bitmap;
import android.graphics.BitmapFactory;
import android.graphics.Canvas;
import android.graphics.Color;
import android.graphics.Paint;
import android.graphics.PorterDuff;
import android.graphics.PorterDuffXfermode;
import android.graphics.Rect;
import android.graphics.RectF;
import android.util.AttributeSet;
import android.util.Log;
import android.view.View;

public class OverlayView extends View {
    public static interface UpdateDoneCallback { 
        public void onUpdateDone(); 
    }  
   
    private UpdateDoneCallback updateDoneCb = null; 
    private Bitmap targetBMP = null;
    
    public OverlayView(Context c, AttributeSet attr) {
        super(c, attr); 
    }

    public void DrawResult(Bitmap bmp) {
        targetBMP = bmp;
        postInvalidate(); 
    }

    public void setUpdateDoneCallback(UpdateDoneCallback cb) {
        updateDoneCb = cb;
    }

    @Override
    protected void onDraw(Canvas canvas) {
        if ( targetBMP != null ) {            
            RectF targetRect = new RectF(0, 0, canvas.getWidth(), canvas.getHeight());
            canvas.drawBitmap(targetBMP, null, targetRect, null);
                        
            if ( updateDoneCb != null)
                updateDoneCb.onUpdateDone();
        }
    }

}
