#ifndef ncsv_HPP
#define ncsv_HPP
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "any.h"

#define LINESIZE   8 * 1024 * 256
int isspace ( int c );
extern char *Trim(char *str);
extern const char* move(const char *p);
extern const char* trimFromEnd(char *p, const char* beginning);

#define MAXCOLUMNS 32
struct ncsv {
    marray<char*>         lines;
    marray<any>           head;  //Support upto 32 columns now    
    marray<any>           data[MAXCOLUMNS];  //Support upto 32 columns now    
    int ncols; // n columns
    int nrows; // n Rows
    char * fileName;
    
    ncsv(): ncols(0), nrows(0),fileName(NULL) {    }
    ~ncsv() {
        if (fileName) delete fileName;        
    }
    void add(char * p)      {  lines[lines.n] = p;}
    
    int process(const char* iline, int header){
        char * line = strdup(iline);
        add(line);
        int i = 0;
        const char  *data[128];
        
        data[i++] = move(line);
        
        for (char * c=line; *c ; c++){
            if ( *c == '"'){
                c++;
                while (*c != '"') c++;
                c++;                
                *c = '\0';
                trimFromEnd(c, data[i]);
                data[i++] = move(c);
            } else if ( *c == ',' ) {
                *c = '\0';
                trimFromEnd(c, data[i]);
                data[i++] = move(c+1);
            }
        }
        //for (int j=0; j < i; j++) printf("%d %s|", j, data[j]); printf("\n");
        if ( header){
            for (int j=0; j < i; j++){                
                this->head[j].set(data[j]);
            }
            ncols = i;
            if (ncols > MAXCOLUMNS) {
                printf("too many columns cannot handle more than %d\n", MAXCOLUMNS);
                exit(1);
            }
        } else{
            for (int j=0; j < i; j++){
                this->data[j][this->data[j].n].set(data[j]);
            }
            nrows++;
        }
        return ncols;
    }
    void dump();
    const char * Read(const char *file1, int nrows=1024*1024, const char *ignore="#");
    
    void ToInt(int column){
        for (int j=0; j < nrows; j++){
            data[column][j].ToInt();
        }
    }
    void ToDouble(int column){
        for (int j=0; j < nrows; j++){
            data[column][j].ToDouble();
        }
    }
    void AddColumn(const char* c, void ** values = NULL);
};

#endif