#include "matrix.h"
#include <time.h>

using namespace std;
int main(int argc, char **argv)
{
    unsigned long long time = 0;

    if (argc < 4 || argc > FIVE) {
        cout << "Incorrect number of arguments of comand line \n";
        return -1;
    }
    int n = atoi(argv[1]);
    int m = atoi(argv[2]);
    int k = atoi(argv[3]);
    if (n <= 0 || m <= 0 || n < m) {
        return -1;
    }
    if (n > TEN_THOUSAND) {
        return -2;
    }
    double *mat;
    mat = new double[n * n];
    if (k == 0) {
        if (argc == 4) {
            delete[] mat;
            return -1;
        }
        bool error3 = kfile(n, mat, argv[4]);
        if (!error3) {
            delete[] mat;
            return -3;
        }
    } else if (k > 0 && k <= 4) {
        kformul(n, mat, k);
    } else {
        delete[] mat;
        return -1;
    }
    printMatrix(n, mat, m);

    bool error4 = is_zero(mat, n);
    if (!error4) {
        delete[] mat;
        return -4;
    }

    double *mat_tmp;
    mat_tmp = new double[n * n];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mat_tmp[i * n + j] = mat[i * n + j];
        }
    }

    double *eigen;
    eigen = new double[n];
    for (int i = 0; i < n; i++) {
        eigen[i] = 0;
    }
    time = currentTimeNano();
    eigenvalues(n, mat_tmp, eigen);
    time = currentTimeNano() - time;

    cout << "Eigenvalues:\n  ";
    for (int i = 0; i < m; i++) {
        cout << eigen[i] << " ";
    }
    cout << "\n";
    printf("Time: %llu ns\n", time);
    double res1 = residual1(n, eigen, mat);
    double res2 = residual2(n, eigen, mat);
    cout << "Residual1: " << res1 << "\n";
    cout << "Residual2: " << res2 << "\n";

    delete[] mat_tmp;
    delete[] mat;
    return 0;
}
