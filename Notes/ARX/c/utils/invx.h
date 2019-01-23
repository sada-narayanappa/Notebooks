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
    int              from;
    int              to;
    char             out[1024];
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
        //printf("BEST=>(%s) Fitness: %.4f (n,m,k): %d %d %d %.5f ", str, fitscore, n,m,k,threshold);
        char t[1024];
        char t1[64];
        t[0] = 0;
        for (int i=0; i < theta.rows(); i++) {
            sprintf(t1,"%.10f,", theta[i]);
            strcat(t, t1);
        }
        t[strlen(t)-2] = 0;
        //cols = 'uName,yName,fitness,correlation,theta,n,m,k,threshold'.split(',')
        sprintf(buff,"%16s,%16s, %.10f, %+.4f, %d,%d,%d,  %13.9f,  \"[%s]\"\n", 
                            u,v,fitscore, corr, n, m, k, threshold,t);
        
        return buff;
    }
};

int Regressors(CONSTEV y, const Eigen::VectorXd& x, int n, int m, int k, Eigen::MatrixXd& r);
void ARXModelLR(const Eigen::VectorXd& y, const Eigen::VectorXd& x, int n, int m, int k, LinearRegression& lr);
double predict3(const Eigen::VectorXd& x, const Eigen::VectorXd& y, 
                int n, int m, int k,Eigen::VectorXd& theta, double t, double& yh, double& rs);
double FitnessScore(const Eigen::VectorXd& x, const Eigen::VectorXd& y, 
                  int n, int m, int k, Eigen::VectorXd& theta);
double ComputeResid(const Eigen::VectorXd& y, const Eigen::VectorXd& x, 
                    Eigen::VectorXd& theta, int n, int m, int k);
double StdDev(const Eigen::VectorXd& v);
void findBest(const Eigen::VectorXd& y, const Eigen::VectorXd& x, Best& best );
double mcorr(const Eigen::VectorXd& x, const Eigen::VectorXd& y);
void CreateInvariants(const char * file, const CSV& df,const Eigen::MatrixXd& xx, 
                      const char *out=NULL, int from=1, int to=200000);

void* runINVX(void *inp);
void SplitRun(int n, MParams& p);
int main_invx(int argc, char const *argv[]);

#endif