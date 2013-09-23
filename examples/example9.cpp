#include <iostream>
#include <vector>

#include <bv/image.h>
#include <bv/image_convert.h>
#include <bv/filter.h>
#include <bv/detector_sift.h>

int run(Eigen::MatrixXd& img, Eigen::MatrixXd& img2, const double sigma) {
    int hfSize = (int)(sigma*3.0);
    if (hfSize < 1) {
        hfSize = 1;
    }
    hfSize += 1;
    Eigen::MatrixXd ker = bv::Kernel::laplacianGaussian( hfSize, sigma); 
    bv::Filter::withTemplate(img, img2, ker);
    return bv::BV_OK;
}

int main(int argc, char *argv[]) {
    if ( argc < 2) {
        std::cout << "Please input bmp file!" << std::endl;
        return -1;
    }
    
    int ret;  
    bv::ColorImage<3> c1( argv[1] );

    bv::Image i1( c1.color(0).width(), c1.color(0).height() );
    Eigen::MatrixXd g1(i1.width(), i1.height());
    Eigen::MatrixXd g2(i1.width(), i1.height());

    bv::Convert::colorImageToGrayImage(c1, i1);
    bv::Convert::grayImageToMatrix( i1, g1);
    
    double maxValue = 0.0;
    double optSigma = 0.0;
    int xx = 287;
    int yy = 78;
    for(double sigma = 3.0; sigma < 6; sigma += 0.3) {
        run(g1, g2, sigma);
        if ( fabs(g2(xx, yy)) > maxValue ) {
            maxValue = fabs(g2(xx, yy));
            optSigma = sigma;
        }
        
        std::cout << "sigma = " << sigma << " LOG=" << fabs(g2(xx,yy)) << std::endl;
    }
    std::cout << "Result sigma = " << optSigma << " MaxValue = " << maxValue << std::endl;
    std::cout << "Result circle = " << optSigma*sqrt(2) << std::endl;

    int r = optSigma*sqrt(2);
    for ( int y = yy - r; y <= yy + r; y++) {
        c1.color(0).data(xx-r, y) = c1.color(2).data(xx-r, y) = 0;
        c1.color(1).data(xx-r, y) = 255;

        c1.color(0).data(xx+r, y) = c1.color(2).data(xx+r, y) = 0;
        c1.color(1).data(xx+r, y) = 255;
    }
    c1.SaveImageToBMP("/tmp/x.bmp");

    return 0;
}
