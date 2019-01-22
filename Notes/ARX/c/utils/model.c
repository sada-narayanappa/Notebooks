#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#include "linreg.h"
#include "Common.h"

int __DEBUG__ = 1;


void logd(char const *argv[] = NULL){
    if (not __DEBUG__) return;
    const char *a;
    a = argv[0];
    
    for (int i=0; a;  i++ ) {
        a = argv[i];
        printf("%s", a);
    }
}
/*
// We use this function to computer ARX b/w u & v
// For creating model using ARX model from ARCH modellers - see ARXClass.ipynb
// than  the LR model hand coded from scratch.
//
void* Regressors(double * y, double * x, int len, int n, int m, int k){
    logd('Regressors {n, m, k}')
    int offset = MAX(n, m + k);
    int llen = len(y) - offset;
    logd( "Length of each array: {llen}, n+m+k={offset}, Length of Original: {len(y)}\n" )
    
    int xxx[] = {0,0,0,0,0,0};
    
    for (int i=0 i <  n+1; i++){
        logd(debug, f'==>{offset-i} - {offset-i+llen} {llen}\n' )
        xxx.append(y[offset-i:offset-i+llen].reshape((llen,1)));
    }
    for i in range(0, m+1):
        upto = len(x) - k
        x1 = x[upto-llen-i: upto-i]
        #x1 = x[i: i+llen]
        xxx.append(x1);
    
    for i in range(len(xxx)):
        logd(debug, "+==>" , len(xxx[i]), "\n")
        
    if (len(xxx) > 0):
        xx = np.hstack(xxx)
    else:
        xx = x
        
    return xx,offset;
}

// Theta is the parameter vector with last index is the constant of ARX model
def ComputeResid(y, x, theta, n,m,k):
    xx, offset= Regressors(y,x, n,m,k)
    rs1=np.dot(xx, theta[:-1]) + theta[-1]
    return rs1

    
def ARXModelLR(y, x, n,m,k, debug=False):
    xx, offset = Regressors(y, x, n,m,k, debug)
    arx = LinearRegression(fit_intercept=True).fit(xx, y[offset:])
    ret = np.append(arx.coef_, arx.intercept_)
    
    return ret, arx, None;


def predict3(x, y, n,m,k, theta, t):
    yh = theta[-1];
    for i in range(n):
        yh += y[t-n+i] * theta[n-1-i]
    for i in range(m+1):
        yh += x[t-k-i] * theta[i+n]
    
    rs = (y[t] - yh)
    return yh ,rs


# Compute the Fitness Score
#
def FitnessScore(x, y, n,m,k,theta, needArrays=True):
    s=max(n,m+k)
    denom = np.sum((y[s:]- np.mean(y[s:]))**2)
    yhat=np.array(y.copy())
    residueFit=[]
    sumResidue = 0;
    yyFit, rrFit =0.5,0.5;
    for t in range(s,len(y)):  # <= predict all possible candidates
        yyFit,rrFit = predict3(x, yhat, n,m,k, theta, t)
        yhat[t] = yyFit
        sumResidue += rrFit ** 2
        
        if (needArrays):
            residueFit.append(rrFit)
   
    fitness = 1- np.sqrt(sumResidue/denom)
    return fitness, yhat,residueFit;

# Compute best fitness score and return the results
#
def findBest(y, x, dfi):
    best=Map({});
    #x=x.reshape((len(x),1))
    fitscore, yh, rs = 0,0,0
    best.fitscore = -1;
    
    theta, arx,theta1 = None, None, None
    for n in range(3):
        for m in range(2):
            for k in range(3):
                theta, arx,theta1 = ARXModelLR(y, x,n,m,k)
                fitscore, yh, rs = FitnessScore(x,y,n,m,k, theta, False)
                    
                #print(f'{(n,m,k)}, {theta}')
                if (best.res is None or fitscore > best.fitscore ): 
                    best.res = theta; best.rs=rs; best.nmk= (n,m,k); best.fitscore = fitscore; 
                    best.n=n; best.m=m; best.k=k;
    
    rs1=ComputeResid(y, x, best.res, best.n, best.m, best.k)
    yh1=y[-len(rs1):]-rs1
    best.threshold=max(abs(yh1)) * 1.05
        
    return best;

            
//#This will create a Invariant file. Note the CSV file has the following format
//# Time, A, B, C, D => first columns in time and time series for subsequent columns
//#
CSV csv;
void* CreateInvariants(file, outFileName=None, columns_from=0, columns_to=100000){
    //csv.Read("../data/test.csv");
    csv.Read(file);
    CSV & df=csv;
    if ( df.nRows <= 2 or df.nColumns <= 1){
        logd("Not enough data in {file} ... ending")
        throw Exception("Not enough data in the dataframe {file}")
    }
    logd(f"Creating invariants using {df.columns[1:]}" )
    
    cols = 'uName,yName,fitness,correlation,theta,n,m,k,threshold'.split(',')
    dfi1 = pd.DataFrame(columns=cols);
    for i,u in enumerate(df.columns[columns_from:]):
        if ( i > columns_to):
            break;
        if ( i == 0 or len(df[u].unique()) <=1 ):
            log.debug(f'either index is 0 : index={i} or not enough unique values in {u}')
            continue;
                    
        for v in df.columns[1:]:
            if (u == v ):
                continue;
            if (len(df[v].unique()) <= 2 ):
                log.debug(f'not enough unique values in {u}')
                continue;
            x=df[u].values    
            y=df[v].values
            print(f"Finding Best of {i}/{len(df.columns)} {u} and {v} \r", end='')
            
            ret = findBest(y, x, dfi1);
            theta = ",".join([str(c) for c in ret.res])
            corr = np.corrcoef(x,y)[0][1]
            
            inv1 = [u,v,ret.fitscore, corr, theta, ret.n, ret.m, ret.k, ret.threshold]
            log.debug(f"{inv1}" )
            dfi1.loc[len(dfi1)] = inv1
    
    dfi1.sort_values(['uName', 'yName'], inplace=True)
    if ( outFileName is not None):
        dfi1.to_csv(outFileName, index=False)
    
    return dfi1
}

*/         
