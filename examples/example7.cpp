#include <iostream>
#include <vector>

#include <bv/image.h>
#include <bv/image_convert.h>
#include <bv/filter.h>
#include <bv/track_lkp.h>

int main(int argc, char *argv[]) {
    if ( argc < 3) {
        std::cout << "Please input bmp files!" << std::endl;
        return -1;
    }
    
    int ret;  
    bv::ColorImage<3> c1( argv[1] );
    bv::ColorImage<3> c2( argv[2] );

    bv::Image i1( c1.color(0).width(), c1.color(0).height() );
    bv::Image i2( c1.color(0).width(), c1.color(0).height() );
    Eigen::MatrixXd g1(i1.width(), i1.height());
    Eigen::MatrixXd g2(i1.width(), i1.height());
    
    bv::Convert::colorImageToGrayImage(c1, i1);
    ret = bv::Convert::grayImageToMatrix( i1, g1);
    bv::Convert::colorImageToGrayImage(c2, i2);
    ret = bv::Convert::grayImageToMatrix( i2, g2);

    bv::TK_LucasKanade tracker;
    std::vector<double> xsource;
    std::vector<double> ysource;
    std::vector<double> xresult;
    std::vector<double> yresult;

    int xx = 182;
    int yy = 186;

    xsource.push_back(xx);
    ysource.push_back(yy);

    tracker.run(g1, g2, xsource, ysource, xresult, yresult);
    std::cout << "<vx, vy> = <" << xresult[0] << "," << yresult[0] << ">" << std::endl;
    
    for ( int y = yy - 5; y <= yy + 5; y++) {
        c1.color(0).data(xx, y) = c1.color(2).data(xx, y) = 0;
        c1.color(1).data(xx, y) = 255;

        c2.color(0).data(xx+xresult[0], y+yresult[0]) = c2.color(2).data(xx+xresult[0], y+yresult[0]) = 0;
        c2.color(1).data(xx+xresult[0], y+yresult[0]) = 255;
    }
    c1.SaveImageToBMP("/tmp/1x.bmp");
    c2.SaveImageToBMP("/tmp/2x.bmp");
}

