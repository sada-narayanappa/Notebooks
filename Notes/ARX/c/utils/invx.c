#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include "Common.h"
#include "Watch.h"
#include "marray.h"
#include "invx.h"
#include "LR.h"
#include <unistd.h> 
#include <pthread.h> 

//---------------------------------------------------------------------------------
struct ARMap{
    int isAR;
    int bestN;
    const char * u;
    float fit;
    ARMap():isAR(0), bestN(1), u(NULL), fit(0){}
};
ARMap *aRMap;
float ARModelThreshold = 0.7;

void CheckARModel(const CSV& df , const Eigen::MatrixXd& xx){
    double fitscore = 0., yh, rs; 
    LinearRegression lr;
    aRMap = new ARMap[ df.nColumns];
    
    for (int i=1; i < df.nColumns; i++ ) {
        const Eigen::VectorXd& x = xx.col(i);
        aRMap[i].u = df.header[i];
        if ( StdDev(x) == 0) {
            aRMap[i].isAR = 1; 
            continue;
        }
        for(int n=1; n < 4; n++) {
            ARXModelLR(x, x,n,-1,0, lr);
            double fitscore = FitnessScore(x,x,n,-1,0, lr.theta);
            //printf("*** %d %s %f \n", n, aRMap[i].u, fitscore);
            if ( fitscore > ARModelThreshold && fitscore > aRMap[i].fit ) { 
                aRMap[i].isAR = 1;
                aRMap[i].bestN = n;
                aRMap[i].fit = fitscore;
            }
        }
    }
    printf("AR Models: ");
    for (int i=0; i < df.nColumns; i++ ) {
        if(aRMap[i].isAR){
            printf("%s:%d; ", aRMap[i].u, aRMap[i].bestN);
        }
    }
}
//---------------------------------------------------------------------------------
int Regressors(const Eigen::VectorXd& y,const Eigen::VectorXd& x,int n,int m,int k,Eigen::MatrixXd& r){
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
                int n, int m, int k,const Eigen::VectorXd& theta, double t, double& rs) {
    //double yh = theta[theta.rows()-1];
    double yh = theta[n+m+1];
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
    double rrFit;
    
    for (int t=s; t < y.rows(); t++) { //  # <= predict all possible candidates
        yhat[t] = predict3(x, yhat, n,m,k, theta, t,rrFit);
        sumResidue += rrFit*rrFit;
    }
    double fitness = 1 - sqrt(sumResidue/denom);
    return fitness;
}

//# Theta is the parameter vector with last index is the constant of ARX model
//#
double ComputeResid(const Eigen::VectorXd& y, const Eigen::VectorXd& x, 
                    Eigen::VectorXd& theta, int n, int m, int k){
    Eigen::MatrixXd r;
    int offset = Regressors(y, x, n,m,k, r);
    Eigen::VectorXd rs1 = (r * theta);
    Eigen::VectorXd yh1= y.segment(y.rows()-r.rows(), r.rows()) - rs1;
    
    double threshold= (yh1.array().abs()).maxCoeff() * 1.05;
    return threshold; 
}

void findBest(const Eigen::VectorXd& y, const Eigen::VectorXd& x, 
              Best& best, int indexOfyColumn, const char* uName ) {
    double fitscore, yh, rs; fitscore =  yh = rs = 0;
    
    best.fitscore = -1;    
    int n,m,k, bestn, bestm, bestk;
    LinearRegression lr;
    
    int nMAX = aRMap[indexOfyColumn].isAR ? 1:3;
        
    for(n=0; n < nMAX; n++)
        for (m=0; m < 2; m++)
            for (k=0; k  < 3; k++) {                
                ARXModelLR(y, x,n,m,k, lr);
                double fitscore = FitnessScore(x,y,n,m,k, lr.theta);
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

double mcorr(const Eigen::VectorXd& x, const Eigen::VectorXd& y) { 
    double sum_X = 0, sum_Y = 0, sum_XY = 0; 
    double squareSum_X = 0, squareSum_Y = 0; 
  
    int n= x.rows();
    for (int i = 0; i < x.rows(); i++) { 
        sum_X       = sum_X + x[i]; 
        sum_Y       = sum_Y + y[i]; 
        sum_XY      = sum_XY + x[i] * y[i]; 
        squareSum_X = squareSum_X + x[i] * x[i]; 
        squareSum_Y = squareSum_Y + y[i] * y[i]; 
    } 
  
    // use formula for calculating correlation coefficient. 
    double d = sqrt((n * squareSum_X - sum_X * sum_X) * (n * squareSum_Y - sum_Y * sum_Y));
    double corr = (n * sum_XY - sum_X * sum_Y) / d;
    return corr; 
} 

double eps= 1e-13;

int allUnique(const Eigen::VectorXd& x){
    for (int i = 1; i < x.rows(); i++) {
        double d = abs(x[i-1] - x[i]);
        if (x[i-1] != x[i] && d > eps) {
            //printf("%f <=> %f  d:%E  %d ", x[i-1], x[i], d, (int)(d<eps) );
            return 0;
        }
    }
    return 1;
}

void CreateInvariants(const char * file, const CSV& df,const Eigen::MatrixXd& xx, 
                      const char *out, int from, int to){    
    Watch w;
    
    Best best;
    char temp[256];
    FILE *ofile = (out)? fopen(out, "w") : stdout;

    int ignore[xx.cols()];
    for (int i = 0; i < xx.cols(); i++) ignore[i] = 0;
        
    if ( from ==1)
        fprintf(ofile, "%s",best.Head());
        //fprintf(ofile, "#File: %s %ld %ld %d %d OUT: %s\n",
        //                    file,xx.rows(),xx.cols(),from, to, out);
    int nvars = 0;
    for (int i = from; i < xx.cols(); i++) {
        if ( i > to )
            break;
        char * u = df.header[i];
        if ( ignore[i] ){
            printf("===> U STDDEV ==0 NO unique values in Y: {%s}\r",u);
            continue;
        }
        const Eigen::VectorXd& x = xx.col(i);
        double sx = StdDev(x);
        if ( allUnique(x) || sx < eps ){
            ignore[i] = 1;
            printf("===> U STDDEV ==0 NO unique values in Y: {%s}\r",u);
            continue;
        }
        //printf("===> Doing {%s} %f %d\n",u, sx, sx < eps );
        for (int j=1; j < xx.cols(); j++) {
            if (i == j || ignore[j] )
                continue;
            
            char *v = df.header[j];
            const Eigen::VectorXd& y = xx.col(j);   
            if ( ignore[j] || allUnique(y) || StdDev(y) < eps ){
                ignore[j] = 1;
                printf("===> V STDDEV ==0 NO unique values in Y: {%s}\r",v);
                continue;
            }
            best.u = u;
            best.v = v;
            findBest(y, x, best, j, v);
            const char * o = best.Dump();
            if (out) { fprintf(ofile, "%s", o);} 
            else printf("%s",o);
            best.corr = mcorr(x,y);
            nvars++;
        }
    }
    (ofile !=stdout)? fclose(ofile) : NULL;
}
//-------------------------------------------------------------------------------
void Test1(const char * file, const CSV& df,const Eigen::MatrixXd& xx, 
                      const char *out, int from, int to){    
    int c = 0;
    int n = xx.cols();
    printf("**********#===> Total Number of Columns: %d \n", (int)(xx.cols()) );
    
    for (int i = 1; i < n; i++) {
        const char * u = df.header[i];
        //if ( strcmp(u, "PIC401.SV") != 0 ) continue;
     
        const Eigen::VectorXd& x = xx.col(i);
        double sx = StdDev(x);
        double eps= 1e-14;
        int au = allUnique(x);
        if ( au || sx < eps) {
            c++;
            printf("+grep %s normal1.inv.xml.csv #%d: %d %f\n",u, i, au, sx);
        } else {
            printf("-grep %s normal1.inv.xml.csv #%d: %d %f\n",u, i, au, sx);
        }
    }
    printf("Total %d\n",c);
}

void* runINVX(void *inp){
    MParams& p1 = *(MParams*)inp;
    CreateInvariants(p1.file, p1.df, p1.xx, p1.out, p1.from, p1.to);
    delete (MParams*)inp;
    return NULL;
}
void append(FILE *head, const char* tailFile) {
    FILE *tail = fopen(tailFile, "rb");

    char buf[BUFSIZ];
    size_t n;
    while ((n = fread(buf, 1, sizeof buf, tail)) > 0)
        if (fwrite(buf, 1, n, head) != n)
            abort();
    if (ferror(tail))
        abort();
    fclose(tail);
}
//---------------------------------------------------------------------------------
// From column - start at 0
void SplitRun(int n, MParams& p) {
    n = n <=0 ? 1:n;
    int each = ceil(p.xx.cols()*1.0/n);
    pthread_t  ids[n+1];
    int  j = 0;
    
    CheckARModel(p.df, p.xx);
        
    int to = MIN(p.xx.cols(), p.to);
    
    for (int i=p.from; i < to-1; i +=each,j++ ) {
        MParams *l1 = new MParams(p);
        MParams& l = *l1;
        l.from = i+1;
        l.to = MIN(i+each, p.xx.cols());
        sprintf(l.out, "%s-OUT-%05d-%05d.csv", l.file,l.from, l.to);
        printf("%d Run from %d - %d (each: %d/%ld) out: %s\n",
                               j, l.from, l.to, each, p.xx.cols(), l.out);
        int ret = pthread_create(&ids[j], NULL, &runINVX, (void*)&l);
        if(ret !=0){
            printf("Thread Creating Failed!!.\n");
            exit(1); 
        }
    }
    for (--j; j >=0; j--){
        pthread_join(ids[j], NULL);
    }
    // COMBINE All files ....
    char  invFile[2*1024];
    sprintf(invFile, "%s.model.csv", p.file);
    printf("\n**Combining all files** into %s\n", invFile);
    
    FILE *file;
    file = fopen(invFile, "wb");
    if (file == NULL)  {
        printf("Cannot combine Files into %s", invFile);
        return ;
    }

    for (int i=0; i < p.xx.cols()-1; i +=each,j++ ) {
        int from = i+1;
        int to = MIN(i+each, p.xx.cols());
        char  tmpFile[2*1024];
        sprintf(tmpFile, "%s-OUT-%05d-%05d.csv", p.file, from, to);
        append(file, tmpFile);
        remove(tmpFile);
    }    
    fclose(file);
}
//-------------------------------------------------------------------------------
void GetMatrix(const CSV& csv, Eigen::MatrixXd & x, int getall=0, int uniq=3);

int main_invx(int argc, char const *argv[]){
    Watch w;
    if (argc <= 1){
        printf("Ex: INVX.exe <csv-file> <#threads> <from-column> <to-column> MAX_N MAX_M MAX_K");
        return 0;
    }
    const char * file = (argc > 1) ? argv[1]: "../data/test1.csv";
    int nThreads = (argc > 2) ? atoi(argv[2]): 16;
    
    int from = (argc > 3) ? atoi(argv[3]): 0;
    int upto = (argc > 4) ? atoi(argv[4]): 1024 * 200;
    
    int ret;
    CSV df(file);
    if ( df.error != NULL || df.nRows <= 2 || df.nColumns <= 1) {
        printf("**NOTE** File not found or Not enough data in %s ... CHECK File contents \n", file);
        return 0;
    }
    Eigen::MatrixXd xx;
    GetMatrix(df, xx, 1, 2);

    printf("## %ld metrics pair-wise compared\n", xx.cols());
    
    MParams p(file, from, upto, df,xx);
    printf("##File: %s %ld %ld %d %d: %s\n",file,
           p.xx.rows(),p.xx.cols(),p.df.nRows,p.df.nColumns,w.Start());

    //lets not create threads if it does not makes sense
    nThreads = (df.nColumns/nThreads) > 1 ? nThreads: 1;
    
    SplitRun(nThreads, p);    
    w.Stop("## Time to Complete: ");
    return 0;
}
