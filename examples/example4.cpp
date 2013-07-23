// Example of LBFGS 

#include <iostream>
#include <bv/minfunc_lbfgs.h>

class SimpleFunc : public bv::MF_LBFGS::IFunction {
protected:
    virtual int dim() const {
        return 3;
    }

    virtual double operator()(const double * x, double * g) const {
        double ret = x[0]*x[0]*x[0] + x[1]*x[1]*x[2] + x[2];
        g[0] = 3 * x[0]*x[0];
        g[1] = 2 * x[1]*x[2];
        g[2] = 1 + x[1]*x[1];
        return ret;
    }
};

int main(int argc, char *argv[]) {
    SimpleFunc *func = new SimpleFunc();
    bv::MF_LBFGS lbfgs(func);

    double x[9] = {10,100,10};
    double ret = lbfgs.run(&x[0]);

    std::cout << "Min value = " << ret << std::endl;
}
