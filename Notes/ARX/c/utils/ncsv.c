#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "marray.h"
#include "any.h"
#include "ncsv.h"

const char* move(const char *p){
    while(*p && isspace((unsigned char)*p)) {
        p++;
    }
    return p;
}
const char* trimFromEnd(char *p, const char* beginning){
    while(p != beginning && (isspace((unsigned char)*p) || *p =='\0') ) {
        *p = '\0';
        p--;
    }
    return p;    
}
const char * ncsv::Read(const char *file1, int nrows, const char *ignore){
    fileName = strdup(file1);
    FILE *file;
    file = fopen(fileName, "r");
    if (file == NULL)  {
        return "File Does not Exists";
    }
    char line[LINESIZE];
    char* ret;
    int i = 0;
    while ( (ret = fgets(line, LINESIZE, file)) && (i < nrows || nrows < 0)){
        char * line1 = Trim(line);
        if (strlen(line1) <= 0 || strncmp(line1, ignore, strlen(ignore)) == 0 ){
            continue;
        }
        process(line1, i == 0);
        i++;
    }
    fclose(file);
    return "";
}

void ncsv::dumphead(FILE* fd){
    for (int i=0; i < ncols; i++)  
        fprintf(fd,"%s%s", head[i].data.c,(i==ncols-1)?"":",");
    fprintf(fd,"\n");
}
void ncsv::dumprow(int row, FILE *fd) {
    int j = row;
    int nd = ncols;
    for (int i=0; i < nd; i++){
        if (data[i][j].type == TYPECHAR)
            fprintf(fd,"%s%s", data[i][j].data.c, (i == nd-1) ?"":"," );
        else if (data[i][j].type == TYPEINT)
            fprintf(fd,"%d%s", data[i][j].data.i, (i == nd-1) ?"":"," );
        else if (data[i][j].type == TYPEDOUBLE)
            fprintf(fd,"%f%s", data[i][j].data.d, (i == nd-1) ?"":"," );
        else
            fprintf(fd,"%s%s", "UKNOWN-TYPE", (i == nd-1) ?"":"," );
        //printf("%d%s", data[i][j].type, (i == ncols-1) ?"":"," );
    }
    fprintf(fd,"\n");
}
void ncsv::dump(FILE* fd){
    dumphead(fd);
    for (int j=0; j < nrows; j++){
        dumprow(j,fd);
    }
    fprintf(fd, "Size: %d x %d\n", nrows, ncols);
}

void ncsv::AddColumn(const char* c, void ** values){
    head[ncols++].set(strdup(c));
    if (ncols > MAXCOLUMNS) {
        printf("too many columns cannot handle more than %d\n", MAXCOLUMNS);
        exit(1);
    }
    for (int row=0; row < nrows; row++){
        void * v = (void *)"";
        data[ncols-1][row].set(0.);
        //int t= data[ncols-1][row].type;
        //printf("%d;",t);
    }
}