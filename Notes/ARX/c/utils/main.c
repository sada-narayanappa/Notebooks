#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

#include <Eigen/Dense>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include "Common.h"
#include "Watch.h"
#include "marray.h"
#include "LR.h"
#include <unistd.h> 
#include <pthread.h> 

int Regressors(const Eigen::VectorXd& y, const Eigen::VectorXd& x, int n, int m, int k, Eigen::MatrixXd& r){
    int offset = MAX(n, m + k);
    int llen  = y.rows() - offset;
    int ncols = n + m +1;
    
    if( ncols <= 1){
        r.resize(llen,2);
        r.col(0) = x.segment(0,llen);
        SetVal(r, -1, 1);
        return offset;
    }
    r.resize(llen, ncols+1);
    int i = 0;
    for (; i < n; i++) {
        r.col(i) = y.segment(offset-i-1,llen);
    }
    for (int j=0; j < m+1; j++) {
        int upto = x.rows() - k;
        r.col(i+j) = x.segment(upto-llen-j, llen);
    }
    SetVal(r, -1, 1);
    return offset;
}

void ARXModelLR(const Eigen::VectorXd& y, const Eigen::VectorXd& x, int n, int m, int k, LinearRegression& lr){
    Eigen::MatrixXd r;
    int offset = Regressors(y, x, n,m,k, r);
    Eigen::VectorXd z = y.segment(offset, y.rows()-offset);                            
    lr.QR(r, z);
}

double predict3(const Eigen::VectorXd& x, const Eigen::VectorXd& y, 
                int n, int m, int k,Eigen::VectorXd& theta, double t, double& yh, double& rs) {
    yh = theta[theta.rows()-1];
    for(int i=0; i <n; i++){
        yh += y[t-n+i] * theta[n-1-i];
    }
    for(int i=0; i < m+1; i++){
        yh += x[t-k-i] * theta[i+n];
    }    
    rs = (y[t] - yh);
    return yh;
}

double FitnessScore(const Eigen::VectorXd& x, const Eigen::VectorXd& y, 
                  int n, int m, int k, Eigen::VectorXd& theta) {
    
    int s= MAX(n,m+k);
    
    Eigen::VectorXd ys = y.segment(s, y.rows()-s);
    double denom = (ys.array() - ys.mean()).array().pow(2).sum();
    
    Eigen::VectorXd yhat= y;
    double sumResidue = 0;
    double yyFit, rrFit;
    
    for (int t=s; t < y.rows(); t++) { //  # <= predict all possible candidates
        predict3(x, yhat, n,m,k, theta, t, yyFit,rrFit);
        yhat[t] = yyFit;
        sumResidue += rrFit*rrFit;
    }
    double fitness = 1 - sqrt(sumResidue/denom);
    return fitness;
}

//# Theta is the parameter vector with last index is the constant of ARX model
//#
double ComputeResid(const Eigen::VectorXd& y, const Eigen::VectorXd& x, Eigen::VectorXd& theta, int n, int m, int k){
    Eigen::MatrixXd r;
    int offset = Regressors(y, x, n,m,k, r);
    Eigen::VectorXd rs1 = (r * theta);
    Eigen::VectorXd yh1= y.segment(y.rows()-r.rows(), r.rows()) - rs1;
    
    double threshold= (yh1.array().abs()).maxCoeff() * 1.05;
    return threshold; 
}

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
        const char* head ="U,V,fitness,corr,n,m,k,threshold,theta\n";
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
        sprintf(buff,"%s,%s,%.9f,%.1f,%d,%d,%d,%.9f,[%s]\n", u,v,fitscore, corr, n, m, k, threshold,t);
        
        return buff;
    }
};

void findBest(const Eigen::VectorXd& y, const Eigen::VectorXd& x, Best& best ) {
    double fitscore, yh, rs; fitscore =  yh = rs = 0;
    
    best.fitscore = -1;    
    int n,m,k, bestn, bestm, bestk;
    LinearRegression lr;
    
    //theta, arx,theta1 = None, None, None
    for(n=0; n < 3; n++)
        for (m=0; m < 2; m++)
            for (k=0; k  < 3; k++) {                
                ARXModelLR(y, x,n,m,k, lr);
                double fitscore = FitnessScore(x,y,n,m,k, lr.theta);
                //lr.Dump();
                if ( fitscore > best.fitscore ) { 
                    best.theta = lr.theta; 
                    best.fitscore = fitscore; 
                    best.n=n; best.m=m; best.k=k;
                }
            }
    best.threshold =ComputeResid(y, x, best.theta, best.n, best.m, best.k);
}
double StdDev(const Eigen::VectorXd& v) {
    int n = v.rows()-1;
    if (v.rows() <= 0) return 0;
    
    double std1 = std::sqrt((v.array() - v.array().mean()).array().square().sum()/n);
    return std1;
}

void CreateInvariants(const char * file, const CSV& df,const Eigen::MatrixXd& xx, 
                      const char *out=NULL, int from=1, int to=200000){    
    
    Watch w;
    
    Best best;
    char temp[256];
    FILE *ofile = (out)? fopen(out, "w") : stdout;
    
    if ( from ==1)
        fprintf(ofile, "%s##File: %s %ld %ld %d %d OUT: %s\n",best.Head(),file,xx.rows(),xx.cols(),from, to, out);

    int nvars = 0;
    for (int i = from; i < xx.cols(); i++) {
        if ( i > to )
            break;
        char * u = df.header[i];
        
        const Eigen::VectorXd& x = xx.col(i);
        if ( StdDev(xx.col(i)) == 0){
            //printf("===> U STDDEV ==0 enough unique values in Y: {%s}\n",u);
            continue;
        }
        for (int j=1; j < xx.cols(); j++) {
            char *v = df.header[j];
            if (i == j )
                continue;
            
            const Eigen::VectorXd& y = xx.col(j);   
            
            if ( StdDev(xx.col(j)) == 0){
                printf("===> V STDDEV ==0 enough unique values in Y: {%s}\r",v);
                continue;
            }
            best.u = u;
            best.v = v;
            findBest(y, x, best);
            const char * o = best.Dump();
            if (out) { fprintf(ofile, "%s", o);} 
            else printf("%s",o);
            nvars++;
            //best.corr = np.corrcoef(x,y)[0][1]
        }
    }
    (ofile !=stdout)? fclose(ofile) : NULL;
}

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

void* run(void *inp){
    MParams& p1 = *(MParams*)inp;
    
    CreateInvariants(p1.file, p1.df, p1.xx, p1.out, p1.from, p1.to);
    
    delete (MParams*)inp;
    return NULL;
}

void SplitRun(int n, MParams& p) {
    n = n <=0 ? 1:n;
    int        each = ceil(p.xx.cols()*1.0/n);
    pthread_t  ids[n+1];
    int  j = 0;
    
    for (int i=0; i < p.xx.cols()-1; i +=each,j++ ) {
            MParams *l1 = new MParams(p);
            MParams& l = *l1;
            l.from = i+1;
            l.to = MIN(i+each, p.xx.cols());
            sprintf(l.out, "OUT-%05d-%05d.out",l.from, l.to);
            printf("%d Run from %d - %d (each: %d/%ld) out: %s\n",
                                   j, l.from, l.to, each, p.xx.cols(), l.out);
            int ret = pthread_create(&ids[j],NULL,&run, (void*)&l);
            if(ret !=0){
                printf("Thread Creating Failed!!.\n");
                exit(1); 
            }

    }
    for (--j; j >=0; j--){
        pthread_join(ids[j], NULL);
    }
}
/*
void test1(){
    Eigen::MatrixXd m(3,3);
    m << 1,2,3,
       4,5,6,
       7,8,9;
    cout << "Here is the matrix m:" << endl << m << endl;
    cout << "2nd Row: " << m.row(1) << endl;
    m.col(2) += 3 * m.col(0);
    cout << "After adding 3 times the first column into the third column, the matrix m is:\n";
    cout << m << endl;
    
    const Eigen::VectorXd& v = m.col(1);
    
    Eigen::VectorXd *vv = (Eigen::VectorXd *) &v;
        
    (*vv)[1] = 890;
    cout << "\n\n" << v << endl;    
    cout << "\n\n" << m << endl;    
}
void test2(Eigen::MatrixXd& xx, int iters=100){
    Watch w;
    printf("%d %s \n", iters, w.Start());
    LinearRegression lr;
    Eigen::VectorXd x;
    Eigen::VectorXd y;
    
    for (int i=0; i < iters; i++) {
        x = xx.col(i % xx.cols());
        if ( StdDev(x) == 0){
            continue;
        }
            
        for (int j=0; j < iters; j++) {
            y = xx.col(j % xx.cols());   
            
            if ( StdDev(y) == 0){
                //printf("===> V STDDEV ==0 enough unique values in Y: {%s}\n",v);
            }
        }
    }
    w.Stop("## Time to Complete: ");
}
*/
//-------------------------------------------------------------------------------
int main(int argc, char const *argv[]){
        
    Watch w;
    const char * file = (argc > 1) ? argv[1]: "../data/test1.csv";
    
    int ret;
    CSV df(file);
    if ( df.error != NULL || df.nRows <= 2 || df.nColumns <= 1) {
        printf("**NOTE** File not found or Not enough data in %s ... ending \n", file);
        return 0;
    }
    Eigen::MatrixXd xx;
    GetMatrix(df, xx, 1, 2);

    printf("## %ld metrics pair-wise compared\n", xx.cols());
    
    MParams p(file, 1, 200000, df,xx);
    printf("##File: %s %ld %ld %d %d: %s\n",file,
           p.xx.rows(),p.xx.cols(),p.df.nRows,p.df.nColumns,w.Start());

    int nThreads = (argc > 2) ? atoi(argv[2]): 16;
    SplitRun(nThreads, p);    
    w.Stop("## Time to Complete: ");
    return 0;
    //test1();
}
