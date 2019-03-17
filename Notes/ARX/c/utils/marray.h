#ifndef MARRAY_H
#define MARRAY_H
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define min(x,y) ( ((x) < (y))? (x): (y))

template<class T> struct marray {
public:
    int n;
    int capacity;
    T*  a;
    marray(int len=1) {
        int cap = len > 1? len:1;
        n = capacity = 0;
        a = new T[n];
        ExtendTo(cap);
    }
    virtual ~marray() {
        n = 0;
        delete [] a;
    }
    void Dump() {
        printf("Length %d: capacity: %d, First few elements: ", n, capacity);
        for(int i=0 ; i < min(10,n); i++) {
            printf("%.2f ", a[i]);
        }
        printf("%s \n", ( n > 10)? "...":"");
    }
    void ExtendTo(int cap) {
        if ( capacity >= cap ) return;
        
        T * b  = new T[cap];
        memset(b,0, sizeof(T)*cap);
        memcpy(b, a, sizeof(T)*n);
        delete [] a;
        
        a = b;
        capacity = cap;
    }

    void Clear()           { n=0; a = (T*) realloc(a, sizeof(T)*0); }
    T& operator [] (int i) { 
        //printf("+%d assigning %d %d\n", i, n, capacity);
        if (i+1 > capacity ) {
            ExtendTo(i+1024);
        }
        if (i+1 > n){
            n = i+1;
        }
        //printf("++%d assigning %d %d\n", i, n, capacity);
        return a[i];  
    }
    const T& operator[](int i) const { return a[i]; }
    
    marray<T>& operator = (const marray &o) { 
        if (this == &o)
            return *this;

        Clear();
        Extend(o.capacity);
        memcpy(a, o.a, sizeof(T)*o.capacity);
        return *this;
    } 
};
#endif