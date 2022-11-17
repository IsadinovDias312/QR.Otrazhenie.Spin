#include "matrix.h"
#include <time.h>

using namespace std;
int main(int argc, char **argv)
{
    unsigned long long time = 0;

    if (argc < FIVE || argc > FIVE + 1) {
        cout << "Incorrect number of arguments of comand line \n";
        return -1;
    }
    int n = atoi(argv[1]);
    int m = atoi(argv[2]);

    double eps;
    int ps;

    if (int(strlen(argv[3])) < 4 || int(strlen(argv[3])) > FIVE) {
        return -1;
    }
    if (int(argv[3][0]) < '1' || int(argv[3][0]) > '9') {
        return -1;
    }
    int e = int(argv[3][0]) - '0';
    if (argv[3][1] != 'e' || argv[3][2] != '-') {
        return -1;
    }
    if (int(argv[3][3]) < '1' || int(argv[3][3]) > '9') {
        return -1;
    }
    if (int(strlen(argv[3])) == 4) {
        ps = int(argv[3][3]) - '0';
    } else {
        if (int(argv[3][4]) < '0' || int(argv[3][4]) > '9') {
            return -1;
        }
        ps = (int(argv[3][3]) - '0') * (FIVE + FIVE) + (int(argv[3][4]) - '0');
    }
    eps = e / pow((FIVE + FIVE), ps);

    int k = atoi(argv[4]);
    if (n <= 0 || m <= 0 || n < m) {
        return -1;
    }
    if (n > TEN_THOUSAND) {
        return -2;
    }
    double *mat;
    mat = new double[n * n];
    if (k == 0) {
        if (argc == FIVE) {
            delete[] mat;
            return -1;
        }
        bool error3 = kfile(n, mat, argv[FIVE]);
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
    int iter = eigenvalues(n, mat_tmp, eigen, eps);
    time = currentTimeNano() - time;

    cout << "Eigenvalues:\n  ";
    for (int i = 0; i < m; i++) {
        cout << eigen[i] << " ";
    }
    cout << "\n";
    printf("Time: %llu ns\n", time);
    double res1 = residual1(n, eigen, mat);
    double res2 = residual2(n, eigen, mat);
    cout << "Residual 1: " << res1 << "\n";
    cout << "Residual 2: " << res2 << "\n";
    cout << "Iterations: " << iter << "\n";

    delete[] eigen;
    delete[] mat_tmp;
    delete[] mat;
    return 0;
}
