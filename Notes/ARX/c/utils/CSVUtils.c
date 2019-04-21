#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "marray.h"
#include "CSV.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

void GetMatrix(const CSV& csv, Eigen::MatrixXd & x, int getall, int uniq){
    if (getall) {
        x.resize( csv.nRows, csv.nColumns );
        for (int i=0; i < csv.nColumns; i++)
            for (int j=0; j < csv.nRows; j++){
                x(j,i) = csv.data[i][j];
            }
        return;
    }
    int count =0;
    int willget[csv.nColumns];
    for (int i=0; i < csv.nColumns; i++) {
        if (!haveNUniqueVals(csv.data[i].a, csv.nRows, uniq)){
            //printf( "NUNIQUE %d IGNORING %s \n", uniq, csv.header[i]); 
            continue;
        }
        double std = mstd(csv.data[i].a, csv.nRows);
        if (std < 1e-17) {
            //printf( "STD == %.9f IGNORING %s \n", std, csv.header[i]); 
            continue;
        }
        willget[count++] = i;
    }
    //printf("#Columns: %d ",count); 
    for (int i =0; i <count; i++) {
        //printf("%s ", csv.header[ willget[i]] );
    }
    x.resize( csv.nRows, count );
    for (int i=0; i < count; i++)
        for (int j=0; j < csv.nRows; j++){            
            int c = willget[i];
            x(j,i) = csv.data[ c ][j];
        }
}
