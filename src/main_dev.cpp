#include "problems.hpp"
#include "min_lag_method.hpp"
#include <cblas.h>
uint n = 1024, m = 1024;

void minl() {
    fp eps = 1e-8;
    uint Nmax = 10000;
    IterativeMethod* meth=  new MinLagMethod(eps, Nmax);

    

    //test_problem(n, m, &method, false);
    problem1(n, m, 1.0, 1.0, meth);
}


void system_check(){

    fp eps = 1e-15;
    uint Nmax = 10000;

    IterativeMethod* meth=  new MinLagMethod(eps, Nmax);
    fp A[] = {3,-2,0,0,
              -2,3,-2,0,
              0,-2,3,-2,
              0,0,-2,3};
    fp a[16];
    fp b[4] = {5,4,3,2};
    fp* x = (fp*)malloc(4*sizeof(fp));

    int kl = 1;
    int ku = 1;
    int m =4;
    int n = 4;
    int lda = 4;
    int ldm = n;
    for (int i = 0; i < m; i++) {
        int k = kl - i;
        for (int j = std::max(0, i - kl); j < std::min(n, i + ku + 1); j++) {
            a[(k + j) + i * lda] = A[j + i * ldm];
        }
    }

    memset(x, 0, 4*sizeof(fp));


    //cblas_dgbmv(CblasRowMajor, CblasNoTrans, 4, 4, 1, 1, 1.0, a, 4, b, 1, 0.0, x, 1);

    meth->run(x,b,a,4);
    printf("===\n");
    for (int i=0;i<4;++i) printf("%lf ", x[i]);

    free(x);


}


int main() 
{

    std::cin>>n>>m;
    minl();
    //system_check();
    return 0;
}
