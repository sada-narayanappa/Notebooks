#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <map>
#include <list>
#include <vector>

#include "Common.h"
#include "Watch.h"
#include "marray.h"
#include "invx.h"
#include "CSV.h"
#include "ncsv.h"
#include "LR.h"
#include <unistd.h> 
#include <pthread.h> 
using namespace std;

float ARModelThreshold = 0.7;

typedef vector< vector<double> > ThetaType;
//---------------------------------------------------------------------------------
template<typename A, typename B>
std::pair<B,A> flip_pair(const std::pair<A,B> &p){
    return std::pair<B,A>(p.second, p.first);
}

template<typename A, typename B>
std::multimap<B,A> flip_map(const std::map<A,B> &src){
    std::multimap<B,A> dst;
    std::transform(src.begin(), src.end(), std::inserter(dst, dst.begin()), 
                   flip_pair<A,B>);
    return dst;
}

double predict3(const double *x, const double *y, 
                int n, int m, int k, const vector<double> &theta, int t, double& rs) {
    double yh = theta[n+m+1];
    
    for(int i=0; i <n; i++){
        yh += y[t-n+i] * theta[n-1-i];
    }
    for(int i=0; i < m+1; i++){
        yh += x[t-k-i] * theta[i+n];
    }    
    rs = (y[t] - yh);
    return yh;
}//---------------------------------------------------------------------------------
void createTheta(const ncsv &csv, Eigen::MatrixXd & theta, vector< vector<double> > &mtheta){
    int nrows = csv.nrows;
    theta.resize( 8, csv.nrows);
    mtheta.resize(csv.nrows, vector<double>(8, 0.));
    
    int column = csv.ncols-1;
    for (int r=0; r < nrows; r++) {
        const char * str = csv.data[column][r].data.c;
        if ( *str == '"')
            str++;
        int n = 0;
        while(*str) {
            theta(n, r) = atof(str);
            mtheta[r][n]= atof(str);
            n++;
            while (*str && *str != ',') str++;
            if (*str == ',') str++;
        }
    }
    //std::cout << theta.block(0,0,8,4) << endl;
}

#define IuName       0
#define IyName       1
#define Ifitness     2
#define Icorrelation 3
#define In           4
#define Im           5
#define Ik           6
#define Ithreshold   7
#define Itheta       8
void ModelMatrix(const ncsv& csv, Eigen::MatrixXd & x){
    x.resize( csv.nrows, csv.ncols );
    for (int i=0; i < 8; i++)
        for (int j=0; j < csv.nrows; j++){
            switch(i){
                case 0: 
                case 1: x(j,i) = (long)csv.data[i][j].data.c; break;
                case 2: 
                case 7: x(j,i) = csv.data[i][j].data.d; break;
                default:x(j,i) = csv.data[i][j].data.i; break;
            }
        }
}
struct BrokenPair{
    const char * u;
    const char * v;
    double res, sig;
    BrokenPair():u(NULL), v(NULL), res(0), sig(0) {}
    BrokenPair(const char *u1, const char * v1, double r1, double s1):
        u(u1), v(v1), res(r1), sig(s1) {
    }
    const char* dump(char * buff = NULL){
        char b[128];
        char * lb = (buff == NULL) ? b:buff;
        sprintf(lb, "%s,%s,%lf,%lf",u,v,res,sig);
        //printf("%s\n",lb);
        return lb;
    }
};
#define MAXBROKENS_TO_PRINT 5000
/*
void FindScore(int t, double sumFit,
               const ncsv &model, Eigen::MatrixXd & xx, 
               const CSV  &adf,   Eigen::MatrixXd & cxx, 
               map<string, int> &edgeCounts, Eigen::MatrixXd& theta,
               map<string, int> &colIdx) {
    
    int n,m,k; 
    double fitness, threshold; 
    const char *uName, *yName;
    
    const Eigen::VectorXd& fitnessA = xx.col(Ifitness);
    const Eigen::VectorXd& na = xx.col(In);
    const Eigen::VectorXd& ma = xx.col(Im);
    const Eigen::VectorXd& ka = xx.col(Ik);
    const Eigen::VectorXd& thresh = xx.col(Ithreshold);
    
    double anomScore = 0;
    int    brknCount = 0;
    
    map<string, double>  brknPairsSig;
    map<string, double>  brknPairsRes;
    map<string, double>  brknSensors;
    map<string, int>     brknSensorsCount;
    
    int debug = 0;
    for(int i=0; i < model.nrows; i++) {
        n = na[i]; m = ma[i]; k = ka[i];
        uName = model.data[0][i].data.c;
        yName = model.data[1][i].data.c;
        threshold = thresh[i];
        fitness   = fitnessA[i];
        const Eigen::VectorXd& theta1 = theta.col(i);
        
        int ui = colIdx[uName];
        int vi = colIdx[yName];
        if (!ui || !vi){
            printf("Column Index not valid %s:%d %s:%d \n",uName,ui, yName,vi);
            continue;
        }
        const Eigen::VectorXd& x = cxx.col(ui-1);
        const Eigen::VectorXd& y = cxx.col(vi-1);
        
        double res;
        double yh = predict3(x, y, n,m,k, theta1, t,res);
        double sigf = (threshold) ? abs(res)/threshold : INFINITY;
        
        if ( debug) {
            printf("%d %s,%s, (%d,%d,%d), fit: %lf, thr: %lf\n",
                       i, uName, yName, n,m,k, fitness, threshold);
            cout << theta1 << endl;
            printf("%lf %lf %lf %lf \n", x[t], y[t], yh, res);
        }
        if ( abs(res) > 1.1 * threshold ){
            anomScore += fitness;
            brknCount++;
            string ks = string(uName)+","+yName;
            brknPairsSig[ks]= sigf;
            brknPairsRes[ks]= res;
                
            brknSensors[uName] += fitness;
            brknSensors[yName] += fitness;
            brknSensorsCount[uName] += 1;
            brknSensorsCount[yName] += 1;
            
            if ( debug) {
                printf("%d %s,%s, (%d,%d,%d), fit: %lf, thr: %lf\n",
                           brknCount, uName, yName, n,m,k, fitness, anomScore);

                printf("%d %s,%s, (%d,%d,%d), fit: %lf, thr: %lf\n",
                           i, uName, yName, n,m,k, fitness, threshold);
                cout << theta1 << endl;
                printf("%lf %lf %lf %lf \n", x[t], y[t], yh, res);
            }
        }
    }
    printf("%d,%d,%lf,%d,\"[",t,(int)cxx(t,0),anomScore,brknCount);
    char buff[128];
    int i =0;
    //Sort Pairs & sensors
    multimap<double, string> dstp = flip_map(brknPairsSig);
    for (multimap<double, string>::reverse_iterator it=dstp.rbegin(); it != dstp.rend(); it++){
        const char * pr = it->second.c_str();
        double res = brknPairsRes[it->second];
        printf("[%s,%f,%f],", pr, res, it->first);
        if ( i++ > MAXBROKENS_TO_PRINT)
            break;
    }
    printf("]\",\"[");
    i=0;
    for ( map<string, double>::iterator it = brknSensors.begin(); it != brknSensors.end(); it++ ){
        brknSensors[it->first] = brknSensors[it->first]/edgeCounts[it->first];
        //brknSensors[it->first] /= edgeCounts[it->first];
        double v = it->second / edgeCounts[it->first];
        //printf("[%s,%lf,%d],",it->first.c_str(), v, brknSensorsCount[it->first]);
    }
    multimap<double, string> dst = flip_map(brknSensors);
    for (multimap<double, string>::reverse_iterator it=dst.rbegin(); it != dst.rend(); it++){
        printf("[%s,%lf,%d],",it->second.c_str(), it->first, brknSensorsCount[it->second]);
        if ( i++ > MAXBROKENS_TO_PRINT)
            break;
    }
    printf("]\"\n");
}
*/
void FindScore1(int t, double sumFit,
               const ncsv &model,
               const CSV  &adf,   
               map<string, int> &edgeCounts, 
               vector< vector<double> > &theta,
               map<string, int> &colIdx,
               marray<int>& uIdx, marray<int>& yIdx ) {
    
    int n,m,k; 
    double fitness, threshold; 
    const char *uName, *yName;
        
    double anomScore = 0;
    int    brknCount = 0;
    
    map<string, double>  brknPairsSig;
    map<string, double>  brknPairsRes;
    map<string, double>  brknSensors;
    map<string, int>     brknSensorsCount;
    
    int debug = 0;
    for(int i=0; i < model.nrows; i++) {
        n = model.data[In][i].data.i; 
        m = model.data[Im][i].data.i; 
        k = model.data[Ik][i].data.i;
        uName = model.data[0][i].data.c;
        yName = model.data[1][i].data.c;
        threshold = model.data[Ithreshold][i].data.d;
        fitness   = model.data[Ifitness][i].data.d;
        const vector<double> & theta1 = theta[i];
        
        //int ui = colIdx[uName];
        //int vi = colIdx[yName];
        int ui = uIdx[i];
        int vi = yIdx[i];
        
        if (!ui || !vi){
            printf("Column Index not valid %s:%d %s:%d \n",uName,ui, yName,vi);
            continue;
        }
        const double * x = adf.data[ui-1].a;
        const double * y = adf.data[vi-1].a;
        double res;
        double yh = predict3(x, y, n,m,k, theta1, t,res);
        double sigf = (threshold) ? abs(res)/threshold : INFINITY; 
        
        if ( abs(res) > 1.1 * threshold ){
            anomScore += fitness;
            brknCount++;
            string ks = string(uName)+","+yName;
            brknPairsSig[ks]= sigf;
            brknPairsRes[ks]= res;
                
            brknSensors[uName] += fitness;
            brknSensors[yName] += fitness;
            brknSensorsCount[uName] += 1;
            brknSensorsCount[yName] += 1;
        }
    }
    printf("=>%d,%d,%lf,%d,\"[",t,(int)adf.data[0][t],anomScore,brknCount);
    char buff[128];
    int i =0;
    //Sort Pairs & sensors
    /*
    for(list<BrokenPair>::iterator it = brknPairs.begin(); it != brknPairs.end(); it++,i++){
        printf("[%s],",(it)->dump(buff));
        if ( i++ > MAXBROKENS_TO_PRINT) break;
    }
    */
    multimap<double, string> dstp = flip_map(brknPairsSig);
    for (multimap<double, string>::reverse_iterator it=dstp.rbegin(); it != dstp.rend(); it++){
        const char * pr = it->second.c_str();
        double res = brknPairsRes[it->second];
        printf("[%s,%f,%f],", pr, res, it->first);
        if ( i++ > MAXBROKENS_TO_PRINT)
            break;
    }
    printf("]\",\"[");
    i=0;
    for ( map<string, double>::iterator it = brknSensors.begin(); it != brknSensors.end(); it++ ){
        brknSensors[it->first] = brknSensors[it->first]/edgeCounts[it->first];
        //brknSensors[it->first] /= edgeCounts[it->first];
        double v = it->second / edgeCounts[it->first];
        //printf("[%s,%lf,%d],",it->first.c_str(), v, brknSensorsCount[it->first]);
    }
    multimap<double, string> dst = flip_map(brknSensors);
    for (multimap<double, string>::reverse_iterator it=dst.rbegin(); it != dst.rend(); it++){
        printf("[%s,%lf,%d],",it->second.c_str(), it->first, brknSensorsCount[it->second]);
        if ( i++ > MAXBROKENS_TO_PRINT)
            break;
    }
    printf("]\"\n");
}

void FindResiduals(const ncsv &model, const CSV &adf, int start=0, int stop=1024*1024){
    Eigen::MatrixXd xx;
    ModelMatrix(model, xx);
    const Eigen::VectorXd& fitness = xx.col(Ifitness);
    const Eigen::VectorXd& na = xx.col(In);
    const Eigen::VectorXd& ma = xx.col(Im);
    const Eigen::VectorXd& ka = xx.col(Ik);
    const Eigen::VectorXd& mka = ma+ka;
    
    int maxLag = max(na.array().maxCoeff(), mka.array().maxCoeff());
    start  = start < maxLag? maxLag: start;
    stop   = min(stop, adf.nRows);
    
    double sumFit = fitness.array().sum();
        
    fprintf(stderr, "==>start: %d, stop: %d sumFit: %lf\n", start, stop, sumFit);
    //start = min(start, na.array().max(),  
    //cout << xx << endl;
    map<string, int> colIdx;
    for (int i=0; i < adf.nColumns; i++) {
        colIdx[adf.header[i]] = i+1;
    }

    marray<int> uIdx;
    marray<int> yIdx;
    map<string, int> edgeCounts;
    for (int i =0; i < model.nrows; i++){
        string s1 = model.data[0][i].data.c;
        string s2 = model.data[1][i].data.c;
        //string ss = s1 + ":"+ s2;
        edgeCounts[s1] += 1;
        edgeCounts[s2] += 1;
        uIdx[i] = colIdx[s1];
        yIdx[i] = colIdx[s2];
    }
    /*
    map<string, int>::iterator it;
    for ( it = edgeCounts.begin(); it != edgeCounts.end(); it++ ){
        std::cout << it->first << ":"<< it->second << endl;
    }
    */
    vector<vector<double> > mtheta;
    Eigen::MatrixXd theta;
    createTheta(model, theta, mtheta);
    //Eigen::MatrixXd cxx;
    //GetMatrix(adf, cxx, 1, 0);
    
    for (int t=start; t < stop; t++){
        //FindScore(t, sumFit, model, xx, adf, cxx, edgeCounts, theta,colIdx );
        FindScore1(t, sumFit, model, adf, edgeCounts, mtheta,colIdx,uIdx, yIdx );
        if (t>1024*1024) 
            break;
    }
}

void test(){
    vector< vector<double> > t1;
    t1.resize(1024*1024, vector<double>(8, 0.));
    for (int i=0; i< 8; i++){
        printf("%f", t1[0][i]);
    }
    printf("\n");
}
int main(int argc, char const *argv[]){
    //test(); return 0;
    Watch w;
    //createTheta();  return 0;
    //const char * mfile = "/NEC/SIAT-OLD/SIAT-OLD/benchmarks/normal1.inv.xml.csv";
    const char * mfile = "/NEC/SIAT-OLD/SIAT-OLD/benchmarks/normal1.csv.model.csv";
    //const char * cfile = "/NEC/SIAT-OLD/SIAT-OLD/benchmarks/ab1.csv";
    const char * cfile = "/NEC/SIAT-OLD/SIAT-OLD/benchmarks/abnormal1.csv";

    fprintf(stderr, "SCORE.EXE %s %s\n", mfile, cfile);    
    
    ncsv model;  
    model.Read(mfile);
    model.ToDouble(2);
    model.ToInt(4);
    model.ToInt(5);
    model.ToInt(6);
    model.ToDouble(7);
    
    //model.dump();
    CSV df(cfile);
    FindResiduals(model, df);
    w.Stop("## Time to Complete: ");
    return 0;
    df.Dump();
}
