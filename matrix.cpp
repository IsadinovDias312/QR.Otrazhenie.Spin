#include "matrix.h"
using namespace std;

unsigned long long currentTimeNano()
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (long long)(t.tv_sec * GIGA_MODIFIER + t.tv_nsec);
}

void GenerateId(int size, double *mas)
{
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            mas[i * size + j] = (double)(i == j);
        }
    }
}

bool kfile(int size, double *mas, char *argv)
{
    int d_size = size * size;
    ifstream input;
    input.open(argv);
    if (!input.is_open()) {
        return false;
    }
    cout << scientific;
    cout << setprecision(3);
    int count = 0;
    while (input >> mas[count]) {
        count++;
        if (count == d_size) {
            break;
        }
    }
    return count >= d_size;
}

void kformul(int size, double *mas, int c)
{
    cout << scientific;
    cout << setprecision(3);
    if (c == 1) {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                mas[i * size + j] = size - max(i, j);
            }
        }
    } else if (c == 2) {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                mas[i * size + j] = 0;
                if (abs(i - j) == 1) {
                    mas[i * size + j] = -1;
                }
            }
            mas[i * size + i] = 2;
        }
    } else if (c == 3) {
        for (int i = 0; i < size - 1; i++) {
            for (int j = 0; j < size - 1; j++) {
                mas[i * size + j] = 0;
            }
            mas[i * size + i] = 1;
            mas[(size - 1) * size + i] = i + 1;
            mas[i * size + size - 1] = i + 1;
        }
        mas[size * (size - 1) + size - 1] = size;
    } else {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                mas[i * size + j] = 1 / double(i + j + 1);
            }
        }
    }
}

void SpinTriangle(int n, double *mat)
{
    for (int i = 1; i < n; i++) {
        double *cos = new double[n - i - 1];
        double *sin = new double[n - i - 1];
        for (int k = 0; k < n - i - 1; k++) {
            cos[k] = 1;
        }
        for (int k = 0; k < n - i - 1; k++) {
            sin[k] = 0;
        }
        for (int k = i + 1; k < n; k++) {
            while (fabs(mat[i * n + i - 1]) < EPS &&
                   fabs(mat[k * n + i - 1]) < EPS && k != n) {
                k++;
            }
            if (k == n) {
                break;
            }
            cos[k - i - 1] = mat[i * n + i - 1] /
                             hypot(mat[i * n + i - 1], mat[k * n + i - 1]);
            sin[k - i - 1] = -(mat[k * n + i - 1] /
                               hypot(mat[i * n + i - 1], mat[k * n + i - 1]));
            double *t = new double[n - i + 1];
            double *q = new double[n - i + 1];
            for (int j = i - 1; j < n; j++) {
                t[j - i + 1] = mat[i * n + j] * cos[k - i - 1] -
                               mat[k * n + j] * sin[k - i - 1];
                q[j - i + 1] = mat[i * n + j] * sin[k - i - 1] +
                               mat[k * n + j] * cos[k - i - 1];
            }
            for (int j = i - 1; j < n; j++) {
                mat[i * n + j] = t[j - i + 1];
                mat[k * n + j] = q[j - i + 1];
            }
            delete[] t;
            delete[] q;
        }
        for (int k = i + 1; k < n; k++) {
            double *t = new double[n - i + 1];
            double *q = new double[n - i + 1];
            for (int j = i - 1; j < n; j++) {
                t[j - i + 1] = mat[j * n + i] * cos[k - i - 1] -
                               mat[j * n + k] * sin[k - i - 1];
                q[j - i + 1] = mat[j * n + k] * cos[k - i - 1] +
                               mat[j * n + i] * sin[k - i - 1];
            }
            for (int j = i - 1; j < n; j++) {
                mat[j * n + i] = t[j - i + 1];
                mat[j * n + k] = q[j - i + 1];
            }
            delete[] t;
            delete[] q;
        }
        delete[] cos;
        delete[] sin;
    }
}

void toQR(int n, double *mas, double *x1, double *x2)
{
    cout << " --- \n";
    cout << mas[(n - 2) * n + n - 2] << "\n";
    for (int i = 0; i < n - 1; i++) {
        double s = mas[(i + 1) * n + i] * mas[(i + 1) * n + i];
        double normA = sqrt(mas[i * n + i] * mas[i * n + i] + s);
        x1[i] = mas[i * n + i] - normA;
        x2[i] = mas[(i + 1) * n + i];
        double normX = sqrt(x1[i] * x1[i] + s);
        cout << normX << " -- normX\n"
             << x1[i] << " " << x2[i] << " " << mas[i * n + i] << "\n";
        x1[i] /= normX;
        x2[i] /= normX;
        double u11 = 1 - 2 * x1[i] * x1[i];
        double u12 = (-2) * x1[i] * x2[i];
        double u22 = 1 - 2 * x2[i] * x2[i];
        for (int j = i; j < n; j++) {
            double tmp1 = u11 * mas[i * n + j] + u12 * mas[(i + 1) * n + j];
            double tmp2 = u12 * mas[i * n + j] + u22 * mas[(i + 1) * n + j];
            mas[i * n + j] = tmp1;
            mas[(i + 1) * n + j] = tmp2;
        }
    }
}

void eigenvalues(int n, double *mas, double *eigen)
{
    SpinTriangle(n, mas);
    cout << "Triangled:\n";
    printMatrix(n, mas, n);
    double *x1;
    double *x2;
    x1 = new double[n - 1];
    x2 = new double[n - 1];
    for (int i = 0; i < n - 1; i++) {
        x1[i] = 0;
        x2[i] = 0;
    }
    for (int s = 0; s < FIVE; s++) {
        toQR(n, mas, x1, x2);
        for (int i = 0; i < n - 1; i++) {
            double u11 = 1 - 2 * x1[i] * x1[i];
            double u12 = (-2) * x1[i] * x2[i];
            double u22 = 1 - 2 * x2[i] * x2[i];
            for (int j = 0; j < i + 2; j++) {
                double tmp1 = mas[j * n + i] * u11 + u12 * mas[j * n + i + 1];
                double tmp2 = mas[j * n + i] * u12 + u22 * mas[j * n + i + 1];
                mas[j * n + i] = tmp1;
                mas[j * n + i + 1] = tmp2;
            }
        }
    }
    for (int s = 0; s < n; s++) {
        eigen[s] = mas[s * n + s];
    }
}

bool is_zero(const double *mas, int size)
{
    double *copy;
    copy = new double[size * size];
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            copy[i * size + j] = mas[i * size + j];
        }
    }
    for (int j = size - 1; j > 0; j--) {
        for (int ii = j; ii >= 0; ii--) {
            if (fabs(copy[ii * size + j]) > EPS) {
                if (ii != j) {
                    for (int k = 0; k < size; k++) {
                        double tmp = copy[ii * size + k];
                        copy[ii * size + k] = copy[j * size + k];
                        copy[j * size + k] = tmp;
                    }
                }
            }
        }
        if (fabs(copy[j * size + j]) < EPS) {
            delete[] copy;
            return false;
        }
        for (int i = 0; i < j; i++) {
            double tmp = copy[i * size + j] / copy[j * size + j];
            for (int k = 0; k < size; k++) {
                copy[i * size + k] -= copy[j * size + k] * tmp;
            }
        }
    }
    if (fabs(copy[0]) < EPS) {
        delete[] copy;
        return false;
    }
    delete[] copy;
    return true;
}

void printMatrix(int size, const double *mas, int m)
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            if (mas[i * size + j] < 0) {
                cout << " " << mas[i * size + j];
            } else {
                cout << "  " << mas[i * size + j];
            }
        }
        cout << "\n";
    }
}

double residual1(int n, double *eigen, double *mat)
{
    double res = 0;
    double eigensum = 0;
    double trace = 0;
    for (int i = 0; i < n; i++) {
        eigensum += eigen[i];
        trace += mat[i * n + i];
    }
    res = fabs(eigensum - trace);
    return res;
}

double residual2(int n, double *eigen, double *mat)
{
    double res = 0;
    double eigensum = 0;
    double trace = 0;
    int n2 = n * n;
    for (int i = 0; i < n2; i++) {
        trace += mat[i];
    }
    for (int i = 0; i < n; i++) {
        eigensum += eigen[i] * eigen[i];
    }
    res = fabs(sqrt(fabs(eigensum)) - sqrt(fabs(trace)));
    return res;
}

double *MatOnMat(int n, const double *mas, const double *mat)
{
    double *rez;
    rez = new double[n * n];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            rez[i * n + j] = 0.;
            for (int k = 0; k < n; k++) {
                rez[i * n + j] += mas[i * n + k] * mat[k * n + j];
            }
        }
    }
    return rez;
}

double *MatMinusMat(int size, const double *mas, const double *mat)
{
    double *rez;
    rez = new double[size * size];
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            rez[i * size + j] = mas[i * size + j] - mat[i * size + j];
        }
    }
    return rez;
}

double residual(int size, const double *mas, const double *mat)
{
    double *id;
    id = new double[size * size];
    GenerateId(size, id);
    double *ed;
    double *nd;
    ed = MatOnMat(size, mas, mat);
    nd = MatMinusMat(size, ed, id);
    delete[] id;
    delete[] ed;
    double maxi = 0;
    for (int j = 0; j < size; j++) {
        double maybemaxi = 0;
        for (int i = 0; i < size; i++) {
            maybemaxi += fabs(nd[i * size + j]);
        }
        if (maybemaxi > maxi) {
            maxi = maybemaxi;
        }
    }
    delete[] nd;
    return maxi;
}
