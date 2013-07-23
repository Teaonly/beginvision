// Example of LBFGS 

#include <iostream>
#include <bv/minfunc_lbfgs.h>

#define N 100

class SimpleFunc : public bv::MF_LBFGS::IFunction {
protected:
    virtual int dim() const {
        return N;
    }

    virtual double operator()(const double * x, double * g) const {
        int i;
        double fx;
        for (i = 0;i < N;i += 2) {
            double t1 = 1.0 - x[i];
            double t2 = 10.0 * (x[i+1] - x[i] * x[i]);
            g[i+1] = 20.0 * t2;
            g[i] = -2.0 * (x[i] * g[i+1] + t1);
            fx += t1 * t1 + t2 * t2;
        }

        return fx;
    }
};

int main(int argc, char *argv[]) {
    SimpleFunc *func = new SimpleFunc();
    bv::MF_LBFGS lbfgs(func);

    double *x = (double *)malloc(N*sizeof(double) );
    for (int i = 0;i < N;i += 2) {
        x[i] = -1.2;
        x[i+1] = 1.0;
    }      
            
    double ret = lbfgs.run(&x[0]);

    std::cout << "Min value = " << ret << std::endl;
    std::cout << x[0] << " " << x[1] << " " << x[2] << std::endl;
}
