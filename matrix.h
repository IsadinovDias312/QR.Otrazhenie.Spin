#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <time.h>
#include <string.h>

#define EPS 1e-18
#define CHECK 1e-10
#define GIGA_MODIFIER 1e9
#define TEN_THOUSAND 1e4
#define KILO_MODIFIER 1e3
#define MILLION 1e6
#define MICRO_MODIFIER 1e-6
#define FIVE 5

unsigned long long currentTimeNano();
void GenerateId(int size, double *mas);
bool kfile(int size, double *mas, char *argv);
void kformul(int size, double *mas, int c);

void SpinTriangle(int n, double *mat);
void toQR(int n, int end, double *mas, double *x1, double *x2);
void Iteration(int start, int end, int n, double *mas, double *x1, double *x2, double strochA, int &iter);
int eigenvalues(int n, double *mas, double *eigenvalues, double eps);

void printMatrix(int size, const double *mas, int m);
bool is_zero(const double *mas, int size);
double residual1(int n, const double *eigen, const double *mat);
double residual2(int n, const double *eigen, const double *mat);

double residual(int size, const double *mas);
