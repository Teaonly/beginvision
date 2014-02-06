package com.beginvision.demo;

import android.app.Activity;
import android.content.Intent;
import android.os.Bundle;
import android.view.View;
import android.view.View.OnClickListener;
import android.widget.Button;
import android.widget.ImageButton;

public class DemoActivity extends Activity
{
    @Override
    public void onCreate(Bundle savedInstanceState)
    {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.main);
    
        Button btn = null;
        btn = (Button)findViewById(R.id.btn_vibe);
        btn.setOnClickListener(vibeAction);

        //Load jni libraries
        System.loadLibrary("beginvision");
    }

    private OnClickListener vibeAction = new OnClickListener() {
        @Override
        public void onClick(View v) {
            Intent intent = new Intent(DemoActivity.this, VibeActivity.class); 
            startActivity(intent);            
        }   
    };
}
