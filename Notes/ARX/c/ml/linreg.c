#include <stdio.h>
#include <stdlib.h>
#include <math.h>                           /* math functions */

#include "linreg.h"

int LinearRegression(int n, const REAL x[], const REAL y[], REAL* m, REAL* b, REAL* r){
    REAL   sumx = 0.0;                      /* sum of x     */
    REAL   sumx2 = 0.0;                     /* sum of x**2  */
    REAL   sumxy = 0.0;                     /* sum of x * y */
    REAL   sumy = 0.0;                      /* sum of y     */
    REAL   sumy2 = 0.0;                     /* sum of y**2  */

    for (int i=0;i<n;i++){ 
        sumx  += x[i];       
        sumx2 += Square(x[i]);  
        sumxy += x[i] * y[i];
        sumy  += y[i];      
        sumy2 += Square(y[i]); 
    } 

    REAL denom = (n * sumx2 - Square(sumx));
    if (denom == 0) {
        // singular matrix. can't solve the problem.
        *m = 0;
        *b = 0;
        if (r) *r = 0;
            return 1;
    }

    *m = (n * sumxy  -  sumx * sumy) / denom;
    *b = (sumy * sumx2  -  sumx * sumxy) / denom;
    if (r!=NULL) {
        *r = (sumxy - sumx * sumy / n) /    /* compute correlation coeff */
              sqrt((sumx2 - sqr(sumx)/n) *
              (sumy2 - sqr(sumy)/n));
    }
    return 0; 
}