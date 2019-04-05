#ifndef any_HPP
#define any_HPP
#include <stdio.h>
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
    short int type;
    short int deletePtr;
    
    any()         : type(TYPECHAR),   deletePtr(0){data.c=NULL;}  
    any(int i    ): type(TYPEINT),    deletePtr(0){data.i=i;}  
    any(double d ): type(TYPEDOUBLE), deletePtr(0){data.d=d;}  
    any(const any& o){
        if ( this == &o)  return;
        
        deletePtr = o.deletePtr;
        type = o.type;
        data = o.data;
        if ( !deletePtr )
            return;
        data.c = strdup(data.c);
    }
    any(const char * c, int makeCopyAndDelete=0):type(TYPECHAR),deletePtr(makeCopyAndDelete) {
        if(makeCopyAndDelete){
            char * cc = strdup(c);
            set(cc, TYPECHAR);
            if(type != TYPECHAR){
                printf("*** Programming error, cannot copy and type!=TYPECHAR\n");
                abort();
            }
        } else {
            set(c, TYPECHAR);
        }
    }
    ~any(){
        if(deletePtr && data.c)
            delete data.c;
    };
    //~~~~~~~~~ Operators
    operator double const ()        { return data.d; }
    operator int const ()           { return data.i; }
    operator const char * const ()  { return data.c; }
    int operator == (const any &o)  { return o.type == type && strcmp(o.data.c, data.c)==0; }
    int operator == (const char *c) { return strcmp(c, data.c)==0; }
    //~~~~~~~~~ Other funcs
    void ToInt() {
        switch (type){
            case TYPECHAR:   data.i = atoi(data.c); break;
            case TYPEDOUBLE: data.i = (int)(data.d); break;
        }
        type = TYPEINT;
    }
    void ToDouble() {
        switch (type){
            case TYPECHAR:   data.d = strlen(data.c) > 0 ? atof(data.c): atof("nan"); break;
            case TYPEINT:    data.d = (double)(data.i); break;
            case TYPEDOUBLE: break;
        }
        type = TYPEDOUBLE;
    }    
    void set(void *v ) {  data.v = v; type = TYPEVOID;   }
    void set(double d) {  data.d = d; type = TYPEDOUBLE; }
    
    const char * set(const char* d, int t= TYPECHAR){
        data.c = d; 
        type =t;
        switch (type){
            case TYPECHAR:   break;
            case TYPEINT:    data.i = atoi(data.c); break;
            case TYPEDOUBLE: data.d = atof(data.c); break;
        }
        return d;
    }
};
#endif