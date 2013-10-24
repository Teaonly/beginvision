#include <iostream>
#include <vector>

#include <bv/image.h>
#include <bv/image_convert.h>

int doRANSAC(const std::vector<double>& data, std::vector<double>& result) {
    
    int matchedNumber = data.size() / 4;

    Eigen::MatrixXd A(matchedNumber*2, 6);
    Eigen::MatrixXd mode(6,1);
    Eigen::MatrixXd B(matchedNumber*2, 1); 
    
    for(int i = 0; i < matchedNumber; i++) {
        double ax, ay, bx, by;
        ax = data[i*4];
        ay = data[i*4+1];
        bx = data[i*4+2];
        by = data[i*4+3];
        
        A(i*2, 0) = ax;
        A(i*2, 1) = ay;
        A(i*2, 2) = 0;
        A(i*2, 3) = 0;
        A(i*2, 4) = 1;
        A(i*2, 5) = 0;
        A(i*2+1, 0) = 0;
        A(i*2+1, 1) = 0;
        A(i*2+1, 2) = ax;
        A(i*2+1, 3) = ay;
        A(i*2+1, 4) = 0;
        A(i*2+1, 5) = 1;

        B(i*2, 0) = bx;
        B(i*2+1, 0) = by;
    }    
    
    std::cout << A.transpose() * A << std::endl;

    return 0;     
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
                 
    }

    return 0;
}




