#include "types.hpp"
#include "measures.hpp"
#include "iterative_methods.hpp"
#include <cblas.h>

#ifndef MLM_HPP
#define MLM_HPP

//#define MVR_ITER_MAX 100000

class MinLagMethod : public IterativeMethod {
public:
    MinLagMethod(fp eps_required, uint iter_max) :
                    IterativeMethod(eps_required, iter_max)
                    { }

    void run(fp* x, const fp* b, const fp* A, uint32_t n) override;
};


inline void MinLagMethod::run(
    fp* v,
    const fp* b,
    const fp* A,
    uint32_t n
) {
    uint iter = 0;
    fp eps = 0.;
    fp* r =(fp*)malloc(sizeof(fp) * (n));
    fp* Ar = (fp*)malloc(sizeof(fp) * (n));
    //inf_lag = 0.0;
    for(; iter < ITER_MAX; iter++) {
        memcpy(r, b, sizeof(fp) * (n));
        cblas_dgbmv(CblasRowMajor, CblasNoTrans, n, n, 1, 1, 1.0, A, n, v, 1, -1.0, r, 1);
        cblas_dgbmv(CblasRowMajor, CblasNoTrans, n, n, 1, 1, 1.0, A, n, r, 1, 0.0, Ar, 1); //matrix Ar = A_mult(r, hsq_inverse, ksq_inverse, mA, S);
        fp Ar_sq = cblas_ddot(n, Ar, 1, Ar, 1);
        fp Ar_r = cblas_ddot(n, r, 1, Ar, 1);
        fp tau = Ar_r/Ar_sq;

        eps = (cblas_damax(n, r, 1));

        cblas_daxpy(n, -tau, r, 1, v, 1);

        //for(int i=0;i<n;++i) printf("%lf ", r[i]);
        //printf("\n"); 
        //printf("eps = %lf, tau = %lf\n", eps, tau);

        if (fabs(eps) < EPS_REQUIRED) {
            iter++;
            break;
        }
    }
    //printf("FIN\n");
    free(Ar);
    free(r);
    setErrorEstimates(eps);
    setCountIter(iter);
}

#endif // MLM_HPP