#include "problems.hpp"
#include "wave.hpp"

FILE* output1 = fopen("tmp1.txt", "w+");
FILE* output2 = fopen("tmp2.txt", "w+");

void problem1(uint n, uint m, fp T, fp a_sq, IterativeMethod* method) {

    const auto phi1 = [&](fp x) -> fp { return -x*x + fp(1); };
    const auto phi2 = [&](fp x) -> fp { return  std::sin(x); };
    const auto xi1 = [&](fp t) -> fp { return std::cos(t); };
    const auto xi2 = [&](fp t) -> fp { return 0; };

    const auto f = [](fp x, fp t) -> fp { return x*x + t; };

    /*const auto phi1 = [&](fp x) -> fp { return 0; };
    const auto phi2 = [&](fp x) -> fp { return 0; };
    const auto xi1 = [&](fp t) -> fp { return 0; };
    const auto xi2 = [&](fp t) -> fp { return 0; };

    const auto f = [](fp x, fp t) -> fp { return x + t; };*/


    matrix V; 
    matrix V2;

    

    //V = (fp*)malloc(sizeof(fp)*(m + 1)*(n + 1));
    //memset(V, 0, sizeof(fp)*(m + 1)*(n + 1));

    solve_wave(V, T, 
                  n, m, a_sq,
                  f, 
                  phi1, phi2, xi1, xi2, method);
    //printf("%lf\n", V[2]);
    fprintf(output1, "%d ", method->count_iter);
    fprintf(output1, "%.17lf ", method->error_estimates);
    method->count_iter = 0;
    method->error_estimates = 0;


    //V2 = (fp*)malloc(sizeof(fp)*(2 *m + 1)*(2*n + 1));


   // memset(V2, 0, sizeof(fp)*(2*m + 1)*(2*n + 1));


    solve_wave(V2, T, 
                  2 * n, 2 * m, a_sq,
                  f, 
                  phi1, phi2, xi1, xi2, method);
    fp max_diff1 = 0.;
    fp h = 1.0 / (fp)n;
    fp k = T / (fp)m;
    fp mx, mt;
    for (uint j = 1; j <= m; j++) {
        for (uint i = 1; i < n; i++) {
            fp _x = i * h;
            fp _t = j * k;
            //printf("))\n");
            fp diff = std::abs(V2[j*2][i*2] - V[j][i]);
           // printf("%lf %lf\n", V2[j*2][i*2], V[i][j]);
            if (diff > max_diff1) {
                max_diff1 = diff;
                mx = _x;
                mt = _t;
            }
        }
    }

    printf("err = %.17lf\n", max_diff1);
    printf("%lf %lf\n", mx, mt);
    fprintf(output1, "%d ", method->count_iter);
    fprintf(output1, "%.17lf ", method->error_estimates);
    fprintf(output1, "%.17lf ", max_diff1); //"Погрешность решения СЛАУ по норме1 
    fprintf(output1, "%.17lf %.17lf ", mx, mt);//"Максимальная погрешность достигается в x = %.17lf, y = %.17lf\n"
    fprintf(output1, "%.17lf ", method->inf_lag);


   fprintf(output2, "i j x t V(x,t) V2(x,t) |V(x,t)-V2(x,t)|\n");
    for (uint i = 0; i <= n; i++) {
        for (uint j = 0; j <= m; j++) {
            fp _x = i * h;
            fp _t = j * k;
            fp diff = std::abs(V2[j*2][i*2] - V[j][i]);
            fprintf(output2, "%d %d %.10lf %.10lf %.10lf %.10lf %.10lf\n", i, j, _x, _t, V[j][i], V2[j*2][i*2], diff);
        }
    }
    //free(V);
    //free(V2);
}


