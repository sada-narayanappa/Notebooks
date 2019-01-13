#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define LINESIZE   8 * 1024 * 256

struct CSV {
    char*    fileName;
    int      nColumns;
    int      nRows;
    char**   header;
    double** data;
    
    void dump() {
        printf("%s, %d: columns, %d: Rows\n", fileName, nColumns, nRows);
        for (int i=0; i < nColumns; i++) {
            printf("%d %s\n", i, header[i]);
        }
    }
    ~CSV() {
        free(fileName);
        for (int i=0; i < nColumns; i++) {
            free(header[i]);
        }
        //TODO: delete data
    }
};

void get_cols(char *iline, CSV* csv) {
    const char* tok;
    int j =0;
    printf("%s\n\n", iline);
    
    char* line = strdup(iline);
    for (tok = strtok(line, ","); tok && *tok; j++, tok = strtok(NULL, ",\n")) {
    }
    free(line);
    
    csv->nColumns = j;
    csv->header = (char **)malloc(j * sizeof(char *));
    
    //Run through it again
    line = strdup(iline);
    j = 0;
    for (tok = strtok(line, ","); tok && *tok; j++, tok = strtok(NULL, ",\n")) {
        csv->header[j] = strdup(tok); 
    }
    free(line);
}
void read_csv(int nrows, const char *filename, CSV* csv){
    csv->fileName = strdup(filename);
    FILE *file;
    file = fopen(filename, "r");

    char line[LINESIZE];
    int i = 0;
    while (fgets(line, LINESIZE, file) && (i < nrows || nrows < 0)){
        //char* tmp = line; // strdup(line);
        //printf("%s", line);
        
        if ( i == 0 ) {
            get_cols(line, csv);
            i++;
            continue;
        }
        int j = 0;
        const char* tok;
        
        for (tok = strtok(line, ","); tok && *tok; j++, tok = strtok(NULL, ",\n")) {
            //data[i][j] = atof(tok);
            //printf("%f\t", data[i][j]);
        }
        i++;
    }
    fclose(file);
}

int main(int argc, char const *argv[]){
    printf("1. *******Reading CSV start\n");
    CSV csv;
    read_csv(-1, "../data/test.csv", &csv);
    csv.dump();
}

int main1(int argc, char const *argv[]){
    if (argc < 3){
        printf("Please specify the CSV file as an input.\n");
        exit(0);
    }
    int row     = atoi(argv[1]);
    int col     = atoi(argv[2]);
    char fname[256];
    strcpy(fname, argv[3]);

    double **data;
    data = (double **)malloc(row * sizeof(double *));
    for (int i = 0; i < row; ++i){
        data[i] = (double *)malloc(col * sizeof(double));
    }
    //read_csv(row, col, fname, data);

    return 0;
}