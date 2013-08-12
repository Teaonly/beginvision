#ifndef _BV_FILTER_H_
#define _BV_FILTER_H_

#include <math.h>
#include <assert.h>
#include <vector>
#include <Eigen/Core>
#include "beginvision.h"

namespace bv {

class Kernel {
public: 
    static Eigen::MatrixXd gaussian_5d() {
        double cw = 0.375;
        Eigen::VectorXd ker1d(5);
        ker1d << (0.25 - cw/2) , 0.25 , cw , 0.25 , (0.25 - cw/2) ;
        Eigen::MatrixXd ker = ker1d * ker1d.transpose();
        return ker;
    }

    static Eigen::MatrixXd gaussian_5i() {
        Eigen::VectorXd ker1d(5);
        ker1d << 1 , 4 , 6 , 4 , 1;
        Eigen::MatrixXd ker = ker1d * ker1d.transpose();
        return ker;
    }

};

class Filter {
public:    
    enum EXTENTION_MODE {
        EXTENTION_REPEAT = 1,
        EXTENTION_ZERO = 2,
        EXTENTION_SKIP = 3,
    };

    static int gaussianBlur(Eigen::MatrixXd& in, Eigen::MatrixXd& out, int size, double sigma, EXTENTION_MODE mode = EXTENTION_REPEAT) {
        if ( size%2 == 0) {
            return BV_ERROR_PARAMETER;
        }
        
        int hf_size = size >> 1;
        // Check the input and out size  
        if ( (in.rows () != out.rows()) || (in.cols() != out.cols()) ) {
            return BV_ERROR_PARAMETER;
        }

        // Create 1D filter
        std::vector<double> gaussian1D;
        for(int i = -1*hf_size; i <= hf_size; i++) {
            double v = exp(-1*i*i/(2*sigma*sigma));
            gaussian1D.push_back(v);
        } 
        
        // Compute the template sum 
        double tplSum = 0;
        for(int i = 0; i < size; i++) {
            for(int j = 0; j < size; j++) {
                tplSum = tplSum + gaussian1D[i]*gaussian1D[j]; 
            }
        }
        
        int cols = in.cols();
        int rows = in.rows();
        
        int ext_size = hf_size*2;
        Eigen::MatrixXd source(rows + ext_size*2 , cols + ext_size*2 );
        Eigen::MatrixXd source2(rows + ext_size*2 , cols + ext_size*2 );
        if ( mode == EXTENTION_REPEAT ) {
            for (int c = 0; c < cols + ext_size*2; c++) {
                int c_source = ( c + 2*cols - ext_size ) % cols;
                for(int r = 0; r < rows + ext_size*2; r++) {
                    int r_source = ( r + 2*rows - ext_size ) % rows;
                    source(r, c) = in(r_source, c_source);
                }
            } 
        } else {
            return BV_ERROR_PARAMETER;
        }

        for (int c = -1 * hf_size; c < cols + hf_size; c++) {
            for ( int r = -1*hf_size; r < rows + hf_size; r++) {
                double sum = 0;
                for ( int ri = -1 * hf_size; ri <= hf_size; ri++) {
                    sum = sum + gaussian1D[ri + hf_size] * source(r + ri + ext_size, c + ext_size);
                }
                source2(r + ext_size, c + ext_size) = sum;
            }
        }
        
        for(int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                double sum  = 0.0;
                for ( int ci = -1 * hf_size; ci <= hf_size; ci++) {
                    sum = sum + gaussian1D[ci + hf_size] * source2(r + ext_size, c + ci + ext_size);
                }
                out(r, c) = sum / tplSum;
            }
        }

        return BV_OK;
    }

    static int average(Eigen::MatrixXd& in, Eigen::MatrixXd& out, int size, EXTENTION_MODE mode = EXTENTION_REPEAT) { 
        // Change size to a odd value
        if ( size%2 == 0 ) {
            return BV_ERROR_PARAMETER;
        }
        
        // Check the input and out size  
        if ( (in.rows () != out.rows()) || (in.cols() != out.cols()) ) {
            return BV_ERROR_PARAMETER;
        }
        
        int hf_size = size >> 1; 
        int cols = in.cols();
        int rows = in.rows();

        Eigen::MatrixXd source(rows + hf_size*2 , cols + hf_size*2 );
        if ( mode == EXTENTION_REPEAT ) {
            for (int c = 0; c < cols + hf_size*2; c++) {
                int c_source = ( c + 2*cols - hf_size ) % cols;
                for(int r = 0; r < rows + hf_size*2; r++) {
                    int r_source = ( r + 2*rows - hf_size ) % rows;
                    source(r, c) = in(r_source, c_source);
                }
            } 
        } else {
            return BV_ERROR_PARAMETER;
        }
        
        for (int c = 0; c < cols; c++) {
            for(int r = 0; r < rows; r++) {
                double sum  = 0.0;
                for ( int cs = c - hf_size; cs <= c + hf_size; cs++) {
                    for ( int rs = r - hf_size; rs <= r + hf_size; rs++) {
                        sum += source(rs + hf_size, cs + hf_size); 
                    }
                }
                sum /= (size*size);
                out(r, c) = sum;
            }
        }
        return BV_OK;
    }


    static int withTemplate(Eigen::MatrixXd& in, Eigen::MatrixXd& out, 
                       Eigen::MatrixXd& t, EXTENTION_MODE mode = EXTENTION_REPEAT) { 
        if ( t.rows()%2 == 0 || t.cols()%2 == 0 ) {
            return BV_ERROR_PARAMETER;
        }
        int hf_rows = t.rows() >> 1;
        int hf_cols = t.cols() >> 1;
        
        // Check the input and out size  
        if ( (in.rows () != out.rows()) || (in.cols() != out.cols()) ) {
            return BV_ERROR_PARAMETER;
        }
        
        int cols = in.cols();
        int rows = in.rows();

        Eigen::MatrixXd source(rows + hf_rows*2, cols + hf_cols*2);

        if ( mode == EXTENTION_REPEAT ) {
            for (int c = 0; c < cols + hf_cols*2; c++) {
                int c_source = ( c + 2*cols - hf_cols ) % cols;
                for(int r = 0; r < rows + hf_rows*2; r++) {
                    int r_source = ( r + 2*rows - hf_rows ) % rows;
                    source(r, c) = in(r_source, c_source);
                }
            } 
        } else {
            return BV_ERROR_PARAMETER;
        }
        
        double tsum = 0;
        for ( int c = 0; c < t.cols(); c++) {
            for ( int r = 0; r < t.rows(); r++) {
                tsum += t(r,c);
            }
        }
                
        for (int c = 0; c < cols; c++) {
            for(int r = 0; r < rows; r++) {

                double sum  = 0.0;
                for ( int cs = c - hf_cols; cs <= c + hf_cols; cs++) {
                    for ( int rs = r - hf_rows; rs <= r + hf_rows; rs++) {
                        int ct = cs - c + hf_cols;
                        int rt = rs - r + hf_rows;
                        sum += source( rs + hf_rows, cs + hf_cols) * t(ct, rt);
                    }
                }
                sum /= tsum;
                out(r, c) = sum;
            }
        }

        return BV_OK;
    }
    
};    

}
#endif
