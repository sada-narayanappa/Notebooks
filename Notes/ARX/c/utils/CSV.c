#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "marray.h"
#include "CSV.h"
#include <math.h>

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

char *Trim(char *str){
    char *end;
    while(isspace((unsigned char)*str)) str++;

    if(*str == 0) return str;

    end = str + strlen(str) - 1;                // Trim trailing space
    while(end > str && isspace((unsigned char)*end)) end--;

    if ( end > str)
        end[1] = '\0';                              // Write new null terminator character
    return str;
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
            printf("=>%d: ",i);
            for (int j=0; j < nColumns; j++) {
                printf("%3.2f,", data[j][i]);
            }
            printf("\n");
        }
        return;
    }
    int minColumns=10;
    for (int i=0; i < numRows; i++) {
        printf("=>%d: ",i);
        for (int j=0; j < minm(minColumns,nColumns); j++) {
            printf("%3.2f,", data[j][i]);
        }
        if ( nColumns > minColumns)
            printf( "... Omitted %d columns/%d...", nColumns-minColumns, nColumns);
        printf("\n");
    }
    printf("...(%d ... %d) not print ;)\n", numRows, nRows-numRows-1);    
    for (int i=nRows-numRows; i < nRows ; i++) {
        printf("=>%d: ",i);
        for (int j=0; j < minm(minColumns,nColumns); j++) {
            printf("%3.2f,", data[j][i]);
        }
        if ( nColumns > minColumns)
            printf( "... Omitted %d columns/%d...", nColumns-minColumns, nColumns);
        printf("\n");
    }
    printf("--END: %s, Shape= %d X %d \n", fileName, nColumns, nRows);
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
        char * sline = line1;
        
        if (strlen(line1) <= 0 || strncmp(line, ignore, strlen(ignore)) == 0 ){
            continue;
        }
        int j = 0;
        const char* tok;

        while ((tok = strsep(&sline, ",")) != NULL){
            data[j++][i] = strlen(tok)>0 ? atof(tok): atof("nan");
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
