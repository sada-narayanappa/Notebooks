#ifndef CSV_HPP
#define CSV_HPP

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#define LINESIZE   8 * 1024 * 256

int isspace ( int c );

extern char *Trim(char *str);
int haveNUniqueVals(double *a, int n, int nu);
double mmean(double *a, int n);
double mstd(double *a, int n);

struct CSV {
public:
    char*            fileName;
    int              nColumns;
    int              nRows;
    char**           header;
    const char*      msg;
    marray<double>  *data;
    void            *user;
    const char      *error;
    
    CSV() {
        fileName = NULL; header = NULL; msg = NULL; data = NULL; user = NULL;
        nColumns = 0, nRows = 0;
        error = NULL;
    }
    CSV(const char *f){
        fileName = NULL; header = NULL; msg = NULL; data = NULL; user = NULL;
        nColumns = 0, nRows = 0; error = NULL;
        Read(f);
    }
    ~CSV() {
        if ( fileName)
            delete(fileName);
        for (int i=0; i < nColumns; i++) {
            delete(header[i]);
        }
        if ( data )
            delete []data;
    }
    void Dump(int numRows = 4);
    
    void GetColumns(char *iline) {
        const char* tok;
        int j =0;

        char* line = strdup(iline);
        for (tok = strtok(line, ","); tok && *tok; j++, tok = strtok(NULL, ",\n")) {
        }
        free(line);

        nColumns = j;
        header = new char*[j];

        //Run through it again
        line = strdup(iline);
        j = 0;
        for (tok = strtok(line, ","); tok && *tok; j++, tok = strtok(NULL, ",\n")) {
            header[j] = strdup(tok); 
            Trim(header[j]);
        }
        free(line);
        
        data = new marray<double>[nColumns];   //Allocate Array for data
    }

    const char* Read(const char *filename, int nrows = -1, int *columns= NULL, const char *ignore="##");
};

#endif
