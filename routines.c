/* routines.c */

#include <stdlib.h>
#include <math.h>

void spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
{
    double* u = (double*)malloc(sizeof(double) * n);
    double* p = (double*)malloc(sizeof(double) * n);
    double sig, p1, p2, u1, u2;
    int i, k;

    if (yp1 > 0.99e30)
        y2[0] = u[0] = 0.0;
    else {
        y2[0] = -0.5;
        p[0] = 0.0;
        u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
    }

    for (i = 1; i < n - 1; i++) {
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
        p[i] = sig * y2[i - 1] + 2.0;
        y2[i] = (sig - 1.0) / p[i];
        u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
        u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p[i];
    }

    if (ypn > 0.99e30)
        u1 = p1 = 0.0;
    else {
        p1 = 0.0;
        u1 = (3.0 / (x[n - 1] - x[n - 2])) * (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
    }

    y2[n - 1] = (u1 - p1 * u[n - 2]) / (p1 * y2[n - 2] + 1.0);
    for (k = n - 2; k >= 0; k--) {
        y2[k] = y2[k] * y2[k + 1] + u[k];
    }

    free(p);
    free(u);
}

#include <stdlib.h>
#include <math.h>

double gasdev()
{
    static int iset = 0;
    static double gset;
    double fac, rsq, v1, v2;

    if (iset == 0) {
        do {
            v1 = 2.0 * ((double)rand() / RAND_MAX) - 1.0;
            v2 = 2.0 * ((double)rand() / RAND_MAX) - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);

        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        iset = 1;

        return v2 * fac;
    } else {
        iset = 0;
        return gset;
    }
}


void splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
{
    int klo, khi, k;
    double h, b, a;

    klo = 0;
    khi = n - 1;

    while (khi - klo > 1) {
        k = (khi + klo) >> 1;
        if (xa[k] > x)
            khi = k;
        else
            klo = k;
    }

    h = xa[khi] - xa[klo];
    if (h == 0.0) {
        // Handle the case where xa[k] = xa[k+1]
        *y = ya[klo];
        return;
    }

    a = (xa[khi] - x) / h;
    b = (x - xa[klo]) / h;
    *y = a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * (h * h) / 6.0;
}
