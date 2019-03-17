#ifndef any_HPP
#define any_HPP
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

//0 int;1 double; 2 const char*;
#define TYPEINT     0
#define TYPEDOUBLE  1
#define TYPECHAR    2
#define TYPEVOID    3
struct any {
    union Data {
        int     i;
        double  d;
        const char *c;
        void *  v;
    }data;
    int type;
    int deletePtr;
    
    any(): type(TYPECHAR), deletePtr(0){data.c=NULL;}  
    any(const char * c, int del=0) {
        //data.c = c; deletePtr = del;
    }
    void ToInt() {
        if ( type == TYPEINT)          return;
        else if( type == TYPEDOUBLE)   data.i = (int) data.d;
        else if( type == TYPECHAR  )   data.i = atoi(data.c);
        
        type = TYPEINT;
    }
    void ToDouble() {
        if ( type == TYPEDOUBLE)    return;
        else if( type == TYPEINT )  data.d = (double) data.i;
        else if( type == TYPECHAR ) data.d = atof(data.c);
        type = TYPEDOUBLE;
    }    
    void set(void *v){
        data.v = v;
        type = TYPEVOID;
    }
    void set(double d){
        data.d = d;
        type = TYPEDOUBLE;
    }
    const char * set(const char* d, int t= TYPECHAR){
        data.c = d;
             if (t == 0 )  { data.i = atof(d); }
        else if (t == 1 )  { data.d = atof(d);}
        type =t;
        return d;
    }
    const char * tostring(char *buff = NULL) {
        const char * ret = buff;
        switch (type) {
            case TYPEINT:    sprintf(buff,"%d", data.i);  break;
            case TYPEDOUBLE: sprintf(buff,"%lf", data.d); break;
            case TYPECHAR:   ret = data.c; break;
            case TYPEVOID:   sprintf(buff,"Void*: %p", data.v); break;
            default:
               sprintf(buff,"%s", "UKNOWN"); break;
        }
        return ret;
    }    
};
#endif