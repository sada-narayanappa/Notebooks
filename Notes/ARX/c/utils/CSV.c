#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "marray.h"
#include "CSV.h"

#define LINESIZE   8 * 1024 * 256

//-------------------------------------------------------------------------
double mmean(double *a, int n) {
    double sum = 0;
    for (int i=0; i< n; i++)
        sum += a[i];
    return sum/n;
}
double mstd(double *a, int n) {
    double mean = mmean(a,n);
    double var =0;
    for (int i=0; i< n; i++)
        var += (a[i] - mean) * (a[i] - mean);
    return sqrt(var/(n-1));
}
int inSet(double *s, double d, int n) {
   for (int i=0; i < n; i++)
      if (s[i] == d)
          return 1;
   return 0;
}
int haveNUniqueVals(double *a, int n, int nu) {
    double u[nu];
    if ( n < nu )
        return 0;
    int uf = 0;    
    for (int i=0; i < n; i++){
        if (!inSet(u, a[i], uf)){
            u[uf++]=a[i];
        }
        if ( uf == nu )
            return 1;
    }
    return 0;
}

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

//-------------------------------------------------------------------------
void CSV::Dump(int numRows) {
    printf("%s, %d columns, %d Rows\n", fileName, nColumns, nRows);
    for (int i=0; i < nColumns; i++) {
        printf("%s,", header[i]);
    }
    printf("\n");

    // print all rows and be done!!
    if (nRows <= 2*numRows) {
        for (int i=0; i < numRows; i++) {
            printf("%d: ",i);
            for (int j=0; j < nColumns; j++) {
                printf("%3.2f,", data[j][i]);
            }
            printf("\n");
        }
        return;
    }
    
    for (int i=0; i < numRows; i++) {
        printf("%d: ",i);
        for (int j=0; j < nColumns; j++) {
            printf("%3.2f,", data[j][i]);
        }
        printf("\n");
    }
    printf("...\n");    
    for (int i=nRows-numRows; i < nRows ; i++) {
        printf("%d: ",i);
        for (int j=0; j < nColumns; j++) {
            printf("%3.2f,", data[j][i]);
        }
        printf("\n");
    }
}

const char* CSV::Read(const char *filename, int nrows, int *columns, const char *ignore){
    fileName = strdup(filename);
    FILE *file;
    file = fopen(filename, "r");
    if (file == NULL)  {
        error = "File Does not Exists";
        return error;
    }
    char line[LINESIZE];

    char* ret = fgets(line, LINESIZE, file);
    if (ret == NULL)  return "Could not read File contents";

    GetColumns(line);
    int i = 0;
    while (fgets(line, LINESIZE, file) && (i < nrows || nrows < 0)){
        char * line1 = Trim(line);
        if (strlen(line1) <= 0 || strncmp(line, ignore, strlen(ignore)) == 0 ){
            continue;
        }
        int j = 0;
        const char* tok;

        for (tok = strtok(line, ","); tok && *tok; j++, tok = strtok(NULL, ",\n")) {
            data[j][i] = atof(tok);
        }
        i++;
        nRows++;
    }
    fclose(file);
    return NULL;
}

int main_read_csv(int argc, char const *argv[]){
    printf("1. *******Reading CSV start\n");
    CSV csv;
    csv.Read("../data/test.csv");
    csv.Dump();
    return 0;
}
