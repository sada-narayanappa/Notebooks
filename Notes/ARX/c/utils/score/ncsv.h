#ifndef ncsv_HPP
#define ncsv_HPP
#include <stdio.h>
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
    
    int process(const char* iline, int header);

    void dumphead(FILE * fd=stderr);
    void dumprow(int i, FILE *fd=stderr);
    void dump(FILE* fd=stderr);
    const char * Read(const char *file1, int nrows=1024*1024, const char *ignore="#");
    
    void ToInt(int column){
        for (int j=0; j < nrows; j++)
            data[column][j].ToInt();
    }
    void ToDouble(int column){
        for (int j=0; j < nrows; j++)
            data[column][j].ToDouble();
    }
    void AddColumn(const char* c, void ** values = NULL);
};

#endif