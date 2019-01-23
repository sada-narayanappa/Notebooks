#ifndef CSV_HPP
#define CSV_HPP

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#define LINESIZE   8 * 1024 * 256

int isspace ( int c );


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
    CSV(const char *f) {
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
    
    char *Trim(char *str) {
        char *end;
        while(isspace((unsigned char)*str)) str++;

        if(*str == 0) return str;
        
        end = str + strlen(str) - 1;                // Trim trailing space
        while(end > str && isspace((unsigned char)*end)) end--;
        
        end[1] = '\0';                              // Write new null terminator character
        return str;
    }
    
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

//Must have at least 3 unique values
void GetMatrix(const CSV& csv, Eigen::MatrixXd & x, int getall=0, int uniq=3);

#endif