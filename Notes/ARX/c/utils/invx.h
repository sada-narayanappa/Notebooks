#ifndef INVX_H
#define INVX_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include "Common.h"
#include "Watch.h"
#include "marray.h"
#include "LR.h"
#include <unistd.h> 
#include <pthread.h> 

typedef const Eigen::VectorXd& CONSTEV;

struct invx{
    
};


struct MParams {
    const char      *file;
    const char      *modelFile;
    int              from;
    int              to;
    char             out[2*1024];
    const CSV        &df;
    Eigen::MatrixXd  &xx;
    MParams(const MParams& c):  df(c.df), xx(c.xx) {
        file = c.file;
        from = c.from; to = c.to;
        strcpy(out,c.out);
    }
    MParams(const char* fi, int f, int t,const CSV &df1, Eigen::MatrixXd  &xx1): 
                df(df1), xx(xx1) {
        file = fi;
        from = f; to = t;
        out[0]=0;
    }
    MParams(const char* fi, const char* mfile, int f, int t,const CSV &df1, Eigen::MatrixXd  &xx1): 
                df(df1), xx(xx1) {
        file = fi;
        modelFile = mfile;
        from = f; to = t;
        out[0]=0;
    }
};


struct Best{
    const char * u;
    const char * v;
    double           fitscore;
    int              n,m,k;
    Eigen::VectorXd  theta;
    double           threshold;
    double           corr;
    char             buff[1024];
    Best() {
        corr = 1;
    }
    static const char* Head() {
        const char* head ="uName,yName,fitness,correlation,n,m,k,threshold,theta\n";
        return head;
    }
    const char* Dump() {
        char t[1024];
        char t1[64];
        t[0] = 0;
        for (int i=0; i < theta.rows(); i++) {
            sprintf(t1,"%+1.16f,", theta[i]);
            strcat(t, t1);
        }
        t[strlen(t)-2] = 0;
        //cols = 'uName,yName,fitness,correlation,theta,n,m,k,threshold'.split(',')
        sprintf(buff,"%16s,%16s, %+.16f, %+.4f, %d,%d,%d,  %24.16f,\"%s\"\n", 
                            u,v,fitscore, corr, n, m, k, threshold,t);
        
        return buff;
    }
};

int Regressors (CONSTEV y, CONSTEV x, int n, int m, int k, Eigen::MatrixXd& r);
void ARXModelLR(CONSTEV y, CONSTEV x, int n, int m, int k, Eigen::VectorXd& lr);
double FitnessScore(CONSTEV x, CONSTEV y, 
                  int n, int m, int k, Eigen::VectorXd& theta);
double ComputeResid(CONSTEV y, CONSTEV x, 
                    Eigen::VectorXd& theta, int n, int m, int k);
double StdDev(CONSTEV v);
void findBest(CONSTEV y, CONSTEV x, Best& best );
double mcorr(CONSTEV x, CONSTEV y);
void CreateInvariants(const char * file, const CSV& df,const Eigen::MatrixXd& xx, 
                      const char *out=NULL, int from=1, int to=200000);

void* runINVX(void *inp);
void SplitRun(int n, MParams& p);
int main_invx(int argc, char **argv);

#endif
