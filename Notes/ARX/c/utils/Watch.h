#ifndef WATCH_H
#define WATCH_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>

struct Watch {
    struct timeval t1, t2;
    char   msg[40];
    
    Watch() { Start(); };
    char* Start() { 
        gettimeofday(&t1, NULL); 
        struct tm* ptm = localtime (&t1.tv_sec);
        strftime (msg, sizeof (msg), "%Y-%m-%d %H:%M:%S", ptm);
        return (char*) msg;  
    }
    long Stop (const char *print=NULL) { 
        gettimeofday(&t2, NULL);
        int seconds = (t2.tv_sec - t1.tv_sec);
        //sprintf(msg, "%lu", ret); 
        char * now = Start(); 
        if (print) {
            int hr = seconds / 60 /60;
            int mn = (seconds / 60) % 60;
            int se = (seconds % 60);
            long micros = ((seconds * 1000000) + t2.tv_usec) - (t1.tv_usec);
            printf("%s %s: %d => %0d:%0d:%0d %ld micros\n", print, now,seconds, hr,mn,se, micros);
        }
        return seconds;
    }
};

#endif
