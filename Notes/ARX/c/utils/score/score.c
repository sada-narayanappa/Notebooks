#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <map>
#include <list>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>

#include "Common.h"
#include "Watch.h"
#include "marray.h"
#include "CSV.h"
#include "ncsv.h"

using namespace std;

float ARModelThreshold = 0.7;

typedef vector< vector<double> > ThetaType;
map<char, any> opts;

any *getfopt(char o){
    if (opts.find(o) != opts.end()){
        return &opts[o];
    }
    return NULL;
}
const char* getfopt(char o, const char* def){
    return opts.find(o) != opts.end() ? opts[o].data.c:def;
}
double getfopt(char o, double def){
    return opts.find(o) != opts.end() ? atof(opts[o].data.c):def;
}
int getfopt(char o, int def){
    return opts.find(o) != opts.end() ? atoi(opts[o].data.c):def;
}
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
    //long double yh = theta[n+m+1];
    double yh = theta[n+m+1];
    
    for(int i=0; i <n; i++){
        yh += y[t-n+i] * theta[n-1-i];
        //printf("predict3 y++ %lf %lf %lf\n", y[t-n+i], theta[n-1-i], yh); 
    }
    for(int i=0; i < m+1; i++){
        yh += x[t-k-i] * theta[i+n];
        //printf("predict3 x++ %lf %lf %lf\n" , x[t-k-i], theta[i+n],yh);
    }    
    //long double rs1 = (y[t] - yh);
    rs = (y[t] - yh);
    //printf("res: %e \n", rs);
    return yh;
}

double predict3(const double *x, const double *y, 
                int n, int m, int k, const double *theta, int t, double& rs) {
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

//---------------------------------------------------------------------------------
void createTheta(const ncsv &csv, vector< vector<double> > &mtheta){
    int nrows = csv.nrows;
    mtheta.resize(csv.nrows, vector<double>(8, 0.));
    
    int column = csv.ncols-1;
    for (int r=0; r < nrows; r++) {
        const char * str = csv.data[column][r].data.c;
        if ( *str == '"')
            str++;
        int n = 0;
        while(*str) {
            mtheta[r][n]= atof(str);
            n++;
            while (*str && *str != ',') str++;
            if (*str == ',') str++;
        }
    }
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
int MAXBROKENS_TO_PRINT = 500;

void FindScore1(int t, double sumFit,
               const ncsv &model,
               const CSV  &adf,   
               map<string, int> &edgeCounts, 
               vector< vector<double> > &theta,
               map<string, int> &colIdx,
               marray<int>& uIdx, marray<int>& yIdx ) {
    
    int n,m,k; 
    double fitness, threshold, corr; 
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
        uName = model.data[IuName][i].data.c;
        yName = model.data[IyName][i].data.c;
        threshold = model.data[Ithreshold][i].data.d;
        fitness   = model.data[Ifitness][i].data.d;
        corr   = model.data[Icorrelation][i].data.d;
        if ( strcmp(yName, "FIC301.SV") !=0 || strcmp(uName, "QI301.Q") !=0 ) {
            //continue;
        }
        if (corr <= -2.0)
            continue;

        const vector<double> & theta1 = theta[i];
        
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
        double sigf = (threshold) ? std::abs(res)/threshold : std::numeric_limits<int>::max(); 
        
        if ( std::abs(res) > 1.1 * threshold ){
            anomScore += fitness;
            brknCount++;
            string ks = string(uName)+","+yName;
            brknPairsSig[ks]= sigf;
            brknPairsRes[ks]= res;
                
            brknSensors[uName] += fitness;
            brknSensors[yName] += fitness;
            brknSensorsCount[uName] += 1;
            brknSensorsCount[yName] += 1;
            
            int debug=0;
            if ( debug) {
                printf("-->%d %d %s,%s, (%d,%d,%d), fit: %lf, thr: %e, cumscore: %lf\n",
                           i, brknCount, uName, yName, n,m,k, fitness, threshold, anomScore);

                for (int tt=0; tt < 6; tt++){
                    printf("%f;", theta1[tt]);
                }
                printf("%lf %lf %f %e %e \n\n", x[t], y[t], yh, res, std::abs(res) - 1.1 * threshold);
            }
        }
    }
    printf("=>%d,%d,%lf,%d,\"[",t,(int)adf.data[0][t],anomScore,brknCount);
    if (getfopt('a', (const char*) (NULL))) {
        printf("]\"\n");
        return;
    }
    
    char buff[128];
    int i =0;
    
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
//------------------------------------------------------------------------------
int dmax(const ncsv &model, int idx1, int idx2){
    int ret = -1;
    for (int i=0; i < model.nrows; i++){
        int s1 = model.data[idx1][i].data.i;
        if ( idx2 >= 0)
            s1 += model.data[idx2][i].data.i;
        ret = MAX(ret, s1);
    }
    return ret;
}
double dsum(const ncsv &model, int idx){
    double ret = 0.0;
    for (int i=0; i < model.nrows; i++){
        ret +=  model.data[idx][i].data.d;
    }
    return ret;
}

void FindResiduals(const ncsv &model, const CSV &adf, int start=0, int stop=1024*1024){
    
    int nmax = dmax(model, In, -1); 
    int mkmax= dmax(model, Im, Ik);  
    double sumFit  = dsum (model, Ifitness);
   
    int maxLag = max(nmax, mkmax);
    
    start  = start < maxLag? maxLag: start;
    stop   = minm(stop, adf.nRows);
    
    fprintf(stderr, "==>start: %d, stop: %d Normalize to Fit: %lf\n", start, stop, sumFit);
    //start = minm(start, na.array().max(),  
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
    createTheta(model, mtheta);
    
    for (int t=start; t < stop; t++){
        FindScore1(t, sumFit, model, adf, edgeCounts, mtheta,colIdx,uIdx, yIdx );
        if (t> 1024*1024) 
            break;
    }
}
/*-----------------------------------------------------------------------------------
Parameters: 
 1. model file
 2. csv file upon which you need to do inferences
 3. filter the model -f -e 0.7 -b 0.5 - parameter: ;) you know what I means
 4. start row for inference file
 5. End row for inference file
*/
#include <getopt.h>
void getopts(int argc, char **argv) {
    fprintf(stderr, 
            "SCORE.EXE <options> model-file inference-file \n"
            "     -i *REQUIRED* model-file \n"
            "     -t *REQUIRED* time-series file \n"
            "     -E filter threshold \n"
            "     -B filter both(any) fitness score must be atleast this much \n"
            "     -O output filtered model \n"
            "     -s start row in inference file\n"
            "     -e end row in inference file \n"
            "     -a Print only Anomaly Score \n"
            "     -A Amount of Broken Pairs to Return \n"
            "     -W *IGNORED* window size \n"
           );    

    int c;
    while ((c = getopt (argc, argv, "i:t:E:B:OseaAW:")) != -1){
        opts[c] = optarg ? optarg : "-";
    }
    for (int i = optind; i < argc; i++)
        opts[128+i] = argv[i];
    
    // Must have these options
    if ( opts.find('i') == opts.end() || opts.find('t') == opts.end() ) {
        fprintf(stderr, "*ERROR: Must provide model file and time series file! \n");
        abort();
    }
    
    //Debug - print options
    fprintf(stderr, "WILL USE OPTIONS:\n");
    for ( map<char, any>::iterator it = opts.begin(); it != opts.end(); it++ ){
        int hj=(it->first);
        fprintf(stderr,"  -[%c] (%5d) [%s]\n", it->first, hj, it->second.data.c);
    }
}
/*-----------------------------------------------------------------------------------
 MODEL FILTERS
*/
struct filter_d2 {double f; int i;};
void filter(ncsv &model){
    double ei = getfopt('E', 1.0); 
    double bo = getfopt('B', 1.0);
    
    if ( ei >= 1 || bo <= 0 || bo >= ei) {
        fprintf(stderr, "*ERROR: opt[E] must be >= opt[B] ! \n");
        abort();
    }
    map<string, double> keep;
    map<string, filter_d2> ques;
    double k[2];
    for(int i=0; i < model.nrows; i++) {
        const char * uName = model.data[0][i].data.c;
        const char * yName = model.data[1][i].data.c;
        double f1 = model.data[Ifitness][i].data.d;

        string uv = string(uName) + ","+yName;
        string vu = string(yName) + ","+uName;

        if (ques.find(vu) == ques.end()){
            filter_d2 d; d.f=f1; d.i = i;
            ques[uv]=d;
            continue;
        }
        double f2 = ques[vu].f;
        double ii = ques[vu].i;
        ques.erase (vu);
        
        if ( f1 < bo or f2 < bo) {
            model.data[Icorrelation][ii].data.d = -2;   //# pass # Remove uv, vu
            model.data[Icorrelation][i].data.d  = -2; 
        }
        else if ( f1 < ei and f2 < ei) {
            model.data[Icorrelation][ii].data.d = -2;   //# pass # Remove uv, vu
            model.data[Icorrelation][i].data.d  = -2; 
            //model.data[Ifitness][ii] = 0;   //# pass # Remove uv, vu
            //model.data[Ifitness][i] = 0; 
        }
        else if ( f1 > f2) {
            //keep[uv]=f1;
            model.data[Icorrelation][ii].data.d = -2;   //# pass # Remove uv, vu
        }
        else {
            //keep[vu]=f2;
            model.data[Icorrelation][i].data.d = -2;   //# pass # Remove uv, vu
        }
    }
    for ( map<string, filter_d2>::iterator it = ques.begin(); it != ques.end(); it++ ){
        if (it->second.f < ei) {
            int idx = it->second.i;
            model.data[Icorrelation][idx].data.d = -2;   //# pass # Remove uv, vu
        }
    }
    int fil=0;
    
    if (!getfopt('O', (const char*) (NULL))) {
        return;
    }
    FILE *file;
    char ffile[4*1024];
    sprintf(ffile, "%s_f%f_%f_filter.csv",model.fileName,bo,ei);
    file = fopen(ffile, "r");
    if (file != NULL)  {
        printf("Filter file Already Exists not doing anything...");
        fclose(file);
        return;
    }
    fclose(file);
    file = fopen(ffile, "w");
    
    model.dumphead(file);
    for(int j=0; j < model.nrows; j++) {
        if ( model.data[Icorrelation][j].data.d <= -2) {
            fprintf(file, "#");
            //fil++;
        }
        model.dumprow(j,file);                  
    }
    fclose(file);
    fprintf(file,"#%d x %d filtered: %d%%(%f), : %d\n", model.nrows, model.ncols, fil, 
                   (1.0 * fil/model.nrows), (model.nrows-fil));
}


void test() {
    any a1("Sada");
    any a2(10);
    any a3(11.0);
    print("%s %d %f \n",(const char*)a1, (int)a2, (double)a3);
}
int main(int argc, char **argv){
    //test(); return 0;
    getopts(argc, argv); 
    
    Watch w;
    const char * mfile;
    const char * cfile;
    mfile = "/NEC/SIAT-OLD/SIAT-OLD/benchmarks/normal1.csv.model.csv";
    cfile= "/NEC/SIAT-OLD/SIAT-OLD/benchmarks/ab1.csv";
    mfile = opts['i'];
    cfile = opts['t'];
    

    ncsv model;  

    model.Read(mfile);

    model.ToDouble(Ifitness);
    model.ToDouble(Icorrelation);
    model.ToInt(In);
    model.ToInt(Im);
    model.ToInt(Ik);
    model.ToDouble(Ithreshold);
       
    if ( getfopt('E') || getfopt('B') ){
        printf ("#Filtering the model ... \n");
        filter(model);
    }
    //model.dump(); return 0;
    CSV df(cfile);
    
    MAXBROKENS_TO_PRINT  = getfopt('A', 5000);
    FindResiduals(model, df, getfopt('s',0), getfopt('e',1024*1024));
    w.Stop("## Time to Complete: ");
    return 0;
    df.Dump();
}
