#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <stdarg.h>

#include "array.h"
#include "Watch.h"

int test1(char *s, ...) {
    va_list ap;
    va_start(ap, s);
    char * n1 = (char *)va_arg(ap, char *); 
    printf("=>%s \n", s);
    char * last = NULL;
    while ( n1 != last && n1 != NULL ) {
        printf("=>%p %s \n", n1, n1 );
        last = n1;
        n1 = (char *)va_arg(ap, char *); 
    } 
    va_end(ap);
    return 0;
}
//-------------------------------------------------------------------------------
void testArrays2() {
    marray<double> a = marray<double>();
    a[1] = 1;
    a[3] = 3;
    a.Dump();
    a[7]=2;
    a.Dump();
}
//-------------------------------------------------------------------------------
int main(int argc, char const *argv[]){
    test1( (char *)"sada", "A1", "A2","A3", "A4");
    
    int x[] = {0,0,0,0,0,0,0,0};

    printf("%d %d %d\n", x[1], sizeof(x)); 
    
    test_lr1();
    return 0;
}