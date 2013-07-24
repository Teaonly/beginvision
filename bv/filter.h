#ifndef _BV_FILTER_H_
#define _BV_FILTER_H_

#include <assert.h>
#include <Eigen/Core>
#include "beginvision.h"

namespace bv {

class Filter {
public:    
    enum EXTENTION_MODE {
        EXTENTION_REPEAT = 1,
        EXTENTION_ZERO = 2,
        EXTENTION_SKIP = 3,
    };

    static int average(Eigen::MatrixXd& in, Eigen::EigenMatrixXd& out, int size, EXTENTION_MODE mode = EXTENTION_REPEAT) { 
        // Change size to a odd value
        if ( size%2 == 0 ) {
            return BV_ERROR_PARAMETER;
        }
        int hf_size = size >> 1; 
        
        // Check the input and out size  
        if ( (in.rows () != out.rows()) || (in.cols() != out.cols()) ) {
            return BV_ERROR_PARAMETER;
        }
        
        int cols = in.cols();
        int rows = in.rows();

        
        Eigen::EigenMatrixXd source(rows + hf_size*2, cols + hf_size*2);

        if ( mode == EXTENTION_REPEAT ) {
            for (int c = 1; c <= cols + hf_size*2; c++) {
                int c_source = ( c + 2*cols - hf_size ) % cols;
                if ( c_source == 0) c_source = cols;
                for(int r = 1: r <= rows + hf_size*2; r++) {
                    int r_source = ( r + 2*rows - hf_size ) % rows;
                    if ( r_source == 0) r_source = rows;
                    source(r, c) = in(r_source, c_source);
                }
            } 
        } else {
            return BV_ERROR_PARAMETER;
        }
        
        for (int c = 1; c <= cols; c++) {
            for(int r = 1: r <= rows; r++) {
                double sum  = 0.0;
                for ( int cs = c - hf_size; cs <= c + hf_size; cs++) {
                    for ( int rs = r - hf_size; rs <= r + hf_size; cr++) {
                        sum += source(rs + hf_size; cs + hf_size); 
                    }
                }
                sum /= (size*size);
                out(r, c) = sum;
            }
        }
        return BV_OK;
    }

    static int withTemplate(Eigen::MatrixXd& in, Eigen::EigenMatrixXd& out, 
                       Eigne::MatrixXd& t, EXTENTION_MODE mode = EXTENTION_REPEAT) { 
        if ( t.rows()%2 == 0 || t.cols()%2 == 0)
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

        Eigen::EigenMatrixXd source(rows + hf_rows*2, cols + hf_cols*2);

        if ( mode == EXTENTION_REPEAT ) {
            for (int c = 1; c <= cols + hf_cols*2; c++) {
                int c_source = ( c + 2*cols - hf_cols ) % cols;
                if ( c_source == 0) c_source = cols;
                for(int r = 1: r <= rows + hf_rows*2; r++) {
                    int r_source = ( r + 2*rows - hf_rows ) % rows;
                    if ( r_source == 0) r_source = rows;
                    source(r, c) = in(r_source, c_source);
                }
            } 
        } else {
            return BV_ERROR_PARAMETER;
        }
        
        double tsum = t.cols() * t.rows();

        for (int c = 1; c <= cols; c++) {
            for(int r = 1: r <= rows; r++) {
                double sum  = 0.0;
                for ( int cs = c - hf_cols; cs <= c + hf_cols; cs++) {
                    for ( int rs = r - hf_rows; rs <= r + hf_rows; cr++) {
                        int ct = cs - c + hf_cols + 1;
                        int rt = rs - r + hf_rows + 1;
                        sum += source( cs + hf_cols, rs + hf_rows) * t(ct, rt);
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
