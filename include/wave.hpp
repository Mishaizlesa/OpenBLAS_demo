#ifndef POISSON_HPP
#define POISSON_HPP

#include "types.hpp"
#include "measures.hpp"
#include "iterative_methods.hpp"

void solve_wave(
    matrix& V,
    fp T,
    uint n, uint m, fp a_sq,
    FUNCTION_XY_TYPE f,
    FUNCTION_X_TYPE phi1, FUNCTION_X_TYPE phi2,
    FUNCTION_X_TYPE xi1, FUNCTION_X_TYPE xi2,
    IterativeMethod* method
) {

    V = matrix(m+1, vector(n+1));
    fp h = 1.0 / (fp)n;
    fp k = T / (fp)m;

    fp h_sq_inv = 1/(h*h);
    fp k_sq_inv = 1/(k*k);

    fp* A = (fp*)(malloc(sizeof(fp)*(n-1)*(n-1)));
    fp* b = (fp*)(malloc(sizeof(fp)*(n-1)));


    fp mA = k_sq_inv + 2.0 * a_sq * h_sq_inv / 3;
    fp left = -a_sq * h_sq_inv /3;
    fp right = left;  
    A[1] = right;
    A[0] = mA;
    for(int i=1;i<n - 2;++i){
        A[i*(n-1) + i-1] = left;
        A[i*(n-1) + i] = mA;
        A[i*(n-1) + i+1] = right;
    }
    A[(n-1)*(n-1) - 1] = mA;
    A[(n-1)*(n-1) - 2] = left;


    fp* a = (fp*)(malloc(sizeof(fp)*(n-1)*(n-1)));
    int kl = 1;
    int ku = 1;
    int n_ = n-1;
    int m_ = m - 1;
    int lda = n_;
    int ldm =lda;
    for (int i = 0; i < m_; i++) {
        int k = kl - i;
        for (int j = std::max(0, i - kl); j < std::min(int(n_), i + ku + 1); j++) {
            a[(k + j) + i * lda] = A[j + i * ldm];
        }
    }


    for(uint j = 0; j <= m; j++) {
        V[j][0] = xi1(j * k);
        V[j][n]= xi2(j * k);

    }


   for (uint i = 0; i <= n; i++) {
        V[0][i] = phi1(i * h);    //нулевой слой
    }

    for (int i = 1; i < n; ++i){
        fp tmp = phi2 (i * h)
        + k * 0.5 * (a_sq * h_sq_inv * (phi1((i+1)*h) - 2 * phi1(i * h) + phi1((i-1)*h)) + f(i*h,0)); //отыскание 1-го слоя со вторым порядком погрешности по времени
        V [1][i] = k * tmp + V[0][i];
    }

    for (int j = 2; j <=m; ++j) 
    {
        for (int i=1; i < n ;++i)
        {
           b[i-1] = f(h * i, k * (j-1)) + a_sq * h_sq_inv / 3 * 
           (V[j-2][i-1] - 2.0 * V[(j-2)][i] + V[(j-2)][i+1] + 
            V[j-1][i-1] - 2.0 * V[(j-1)][i] + V[(j-1)][i+1]) 
            + k_sq_inv * (V[(j-1)][i]* 2 - V[(j-2)][i]);

           //V[j][i] = (f(i*h, (j-1)*k) + a_sq * h_sq_inv* (V[j-1][i-1]- 2 * V[j-1][i] + V[j-1][i+1]) + k_sq_inv * (2 * V[j-1][i] - V[j-2][i])) * (k*k);
        }

        method->run(&V[j][1], b, a, n - 1);

    }


    free(A);
    free(a);
    free(b);
}

#endif // POISSON_HPP