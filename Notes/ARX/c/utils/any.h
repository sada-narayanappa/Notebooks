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
        const char *  c;
        void *  v;
    }data;
    int type;
    char buff[32];
    const char* str;
    
    any();
        
    void ToInt() {
        if ( type == TYPEINT) {
            return;
        } else if( type == TYPEDOUBLE )
            data.i = (int) data.d;
        else if( type == TYPECHAR )
            data.i = atoi(data.c);
        
        type = TYPEINT;
        dump();
    }
    void ToDouble() {
        if ( type == TYPEDOUBLE) {
            return;
        } else if( type == TYPEINT )
            data.d = (double) data.i;
        else if( type == TYPECHAR )
            data.d = atof(data.c);
            
        type = TYPEDOUBLE;
        dump();
    }    
    void set(void *v){
        data.v = v;
        type = TYPEVOID;
        dump();
    }
    void set(double d){
        data.d = 0.0;
        type = TYPEDOUBLE;
        dump();
    }
    const char * set(const char* d, int t=2){
        data.c = d;
             if (t == 0 )  { data.i = atof(d); }
        else if (t == 1 )  { data.d = atof(d);}
        type =t;
        dump();
        return d;
    }
    const char * dump() {
        switch (type) {
            case TYPEINT:    sprintf(buff,"%d", data.i);  str=buff; break;
            case TYPEDOUBLE: sprintf(buff,"%lf", data.d); str=buff; break;
            case TYPECHAR:   str = data.c; break;
            case TYPEVOID:   sprintf(buff,"Void*: %p", data.v); str=buff; break;
            default:
               sprintf(buff,"%s", "UKNOWN"); str=buff; break;
        }
        return str;
    }    
};
#endif