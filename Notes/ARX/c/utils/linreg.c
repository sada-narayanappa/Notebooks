#include <stdio.h>
#include <stdlib.h>
#include <math.h>                           /* math functions */

#include "marray.h"
#include "linreg.h"
#include "CSV.h"

// Regress on first n columns on array y
// x = a0 * x0 + a1 x1 + a2 * x2 + ... + an * mn + bn 
double randd() { return ((double) rand() / (RAND_MAX)) + 1; }
double Pred1(CSV& csv, double *params, int n, int row, double *y){
    double yh = params[n+1];
    for (int k=1; k <= n; k++) {
        yh += csv.data[k][row] * params[k];
    }
    return yh;
}
double LinearRegression(CSV& csv, double *params, double *y, int n, int nIters = 128, int batch=10, double alpha=0.1){
    int i,j,k,b,m = 1;
    double bp[n+1];
    
    for (j=0; j < csv.nRows; j++) {
        printf("%lf ", y[j]);
    }
    printf("%d \n", n);
    
    for (i =0; i <= n; i++) {
        params[i] = randd();
        bp[i] = 0;
    }

    for (i=0; i < nIters; i++) {
        j = 0;
        for (j=0; j < csv.nRows; j+= batch) {
            double err = 0;
            for (int ii =0; ii <= n; ii++) {
                bp[ii] = 0;
            }
            for (m=0; m < csv.nRows && m < batch; m++) {
                double yh = Pred1(csv, params, n, j+m, y);
                err = yh - y[j+m];
                for (k =0; k < n; k++) {
                   bp[k] += err * csv.data[k][j+m];
                }
                bp[n] += err;
            }
            // Update Paramters now;
            params[n] = alpha * bp[n]/m;
            for (k =0; k < n; k++) {
                params[k] = bp[k]/m;
            }
       }
    }

    double mse = 0;
    for (j=0; j < csv.nRows; j+= batch) {
        double yh = Pred1(csv, params, n, j+m, y);
        mse = (yh - y[j]) * (yh - y[j]);
    }
    return mse/csv.nRows;
}

void test_lr1(){
    CSV csv = CSV();
    csv.Read("../data/test.csv");
    csv.Dump();
    
    double params[csv.nColumns-1];
    LinearRegression(csv, params, csv.data[csv.nColumns-1].a, csv.nColumns-2, 1024,10);
    
    for (int i=0; i < csv.nColumns; i++) {
        printf("%d : %lf \n", i+1, params[i]);
    }
}
