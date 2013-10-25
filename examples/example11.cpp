#include <iostream>
#include <vector>

#include <Eigen/Core>
#include <Eigen/LU> 
#include <bv/image.h>
#include <bv/image_convert.h>
#include <bv/util.h>

int doRANSAC(const std::vector<double>& data, std::vector<double>& result) {
    const int MAX_ITERATE = 32;
    srand48( stime(NULL) ); 
    
    Eigen::MatrixXd bestMode(6,1);
    int bestValue = 0;
    std::vector<int> bestItem;

    int matchedNumber = data.size() / 4;
    for(int i = 0; i < MAX_ITERATE; i++) {
        int r[3];
        r[0] = (int)( matchedNumber * drand48() );
        r[1] = (int)( matchedNumber * drand48() );
        r[2] = (int)( matchedNumber * drand48() );

        if ( r[0] == r[1] && r[1] == r[2] && r[0] == r[2]) {
            continue; 
        } 

        Eigen::MatrixXd A(6,6);
        Eigen::MatrixXd b(6,1);
        Eigen::MatrixXd m(6,1);
        
        for (int j = 0; j < 3; j++)  {
            double ax, ay, bx, by;
            ax = data[ r[j]*4 + 0];
            ay = data[ r[j]*4 + 1];
            bx = data[ r[j]*4 + 2];
            by = data[ r[j]*4 + 3];
            
            A(j*2, 0) = ax;
            A(j*2, 1) = ay;
            A(j*2, 2) = 0;
            A(j*2, 3) = 0;
            A(j*2, 4) = 1;
            A(j*2, 5) = 0;
            A(j*2+1, 0) = 0;
            A(j*2+1, 1) = 0;
            A(j*2+1, 2) = ax;
            A(j*2+1, 3) = ay;
            A(j*2+1, 4) = 0;
            A(j*2+1, 5) = 1;
            
            b(j*2, 0) = bx;
            b(j*2+1, 0) = by;
        }
        m = A.inverse() * b;
        
        Eigen::MatrixXd Xi(2,6);
        Eigen::MatrixXd xo(2,1);
        
        std::vector<int> inlierItem;
        for(int j = 0; j < matchedNumber; j++ ) {
            double ax, ay, bx, by;
            ax = data[ j*4 + 0];
            ay = data[ j*4 + 1];
            bx = data[ j*4 + 2];
            by = data[ j*4 + 3];
 
            Xi(0, 0) = ax;
            Xi(0, 1) = ay;
            Xi(0, 2) = 0;
            Xi(0, 3) = 0;
            Xi(0, 4) = 1;
            Xi(0, 5) = 0;
            Xi(1, 0) = 0;
            Xi(1, 1) = 0;
            Xi(1, 2) = ax;
            Xi(1, 3) = ay;
            Xi(1, 4) = 0;
            Xi(1, 5) = 1;
            
            xo = Xi * m;
            double offset = (xo(0,0) - bx) * (xo(0,0) - bx) + 
                            (xo(1,0) - by) * (xo(1,0) - by);
            offset = sqrt(offset);

            if ( offset < 10) {
                inlierItem.push_back(j);
            }
        } 
        if ( inlierItem.size() > bestValue) {
            bestValue = inlierItem.size();
            bestMode = m;
            bestItem = inlierItem;
        }
    }    
    
    std::cout << " >>>> Best mode from RANSAC: " << std::endl;
    std::cout << bestMode << std::endl;

    // We got the bsetMode ant it's inlier items.
    // Then we should do least squares refine from all matched item
    Eigen::MatrixXd A(bestItem.size()*2, 6);
    Eigen::MatrixXd B(bestItem.size()*2, 1);
    for ( int j = 0; j < bestItem.size(); j++) {
        double ax, ay, bx, by;
        ax = data[ bestItem[j]*4 + 0];
        ay = data[ bestItem[j]*4 + 1];
        bx = data[ bestItem[j]*4 + 2];
        by = data[ bestItem[j]*4 + 3];
        
        A(j*2, 0) = ax;
        A(j*2, 1) = ay;
        A(j*2, 2) = 0;
        A(j*2, 3) = 0;
        A(j*2, 4) = 1;
        A(j*2, 5) = 0;
        A(j*2+1, 0) = 0;
        A(j*2+1, 1) = 0;
        A(j*2+1, 2) = ax;
        A(j*2+1, 3) = ay;
        A(j*2+1, 4) = 0;
        A(j*2+1, 5) = 1;
        
        B(j*2, 0) = bx;
        B(j*2+1, 0) = by;
    }

    Eigen::MatrixXd ret(6,1);

    ret = (A.transpose() * A).inverse() * A.transpose() * B;
    
    std::cout << " >>> Refined value is " << std::endl;
    std::cout << ret << std::endl;
    
    result.clear();
    for ( int i = 0; i < 6; i++) {
        result.push_back(ret(i,0) );    
    }
    return 0;     
}

void doAffine(const char *fileName, const std::vector<double> mode) {
    if ( mode.size() != 6)
        return;
   
    bv::ColorImage<3> c1( fileName );
    bv::Image i1( c1.color(0).width(), c1.color(0).height() );
    Eigen::MatrixXd g1(i1.width(), i1.height());
    Eigen::MatrixXd g2(512, 384);
    bv::Convert::colorImageToGrayImage(c1, i1);
    bv::Convert::grayImageToMatrix( i1, g1);


    Eigen::MatrixXd t(2,1);
    Eigen::MatrixXd m(2,2);
    m(0,0) = mode[0];
    m(0,1) = mode[1];
    m(1,0) = mode[2];
    m(1,1) = mode[3];
    t(0,0) = mode[4];
    t(1,0) = mode[5];

    m = m.inverse();
    int height = g2.cols();
    int width = g2.rows();
    for ( int y = 0; y < height; y++) {
        for ( int x = 0; x < width; x++) {
            Eigen::MatrixXd pbar(2,1);
            Eigen::MatrixXd p(2,1);
            pbar(0,0) = x;
            pbar(1,0) = y;
            p = m * (pbar - t);
            if (    p(0,0) < g1.rows() && p(0,0) > 0 
                 && p(1,0) < g1.cols() && p(1,0) > 0) {
                int xx = p(0, 0);
                int yy = p(1, 0);
                g2(x,y) = g1(xx,yy);
            } else {
                g2(x,y) = 1;
            }
        }
    } 
    
    bv::Util::saveAsImage(g2, "/tmp/xxx.bmp");
}


int main(int argc, char *argv[]) {
    if ( argc < 2) {
        std::cout << "Please input one bmp file!" << std::endl;
        return -1;
    }
    
    std::vector<double> data;
    FILE* fp = fopen("./match.txt", "rb");
    while( !feof(fp) ) {
        double xa,ya,xb,yb;
        if ( fscanf(fp, "%lf %lf %lf %lf", &xa, &ya, &xb, &yb) == 4) {
            data.push_back(xa);
            data.push_back(ya);
            data.push_back(xb);
            data.push_back(yb);           
        }
    }    
    fclose(fp);

    if ( data.size() < 3*4) {
        std::cout << " Match is not enough for affine" << std::endl;
        return 0;    
    }
    
    std::vector<double> result;
    int ret = doRANSAC(data, result);
    if ( ret == 0) {
        doAffine(argv[1], result);  
    }
    return 0;
}




