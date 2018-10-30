import logging as log
import sys
import os
import pandas as pd
import numpy as np;
import datetime
import getopt
from collections import defaultdict
import math
import re;
import collections
import gc
import json
from Jupytils import Map
from sklearn.linear_model import LinearRegression
from numba import njit

def logd(debug = True, *args):
    if not debug: return
    for a in args:
        print(a, sep='', end=' ')

log.getLogger().setLevel(log.INFO)

# We use this function to computer ARX b/w u & v- ARXModelLR is 100 times faster
# For creating model using ARX model from ARCH modellers - see ARXClass.ipynb
# than  the LR model hand coded from scratch.
#
# 
def Regressors(y, x, n,m,k, debug=False):
    x = x.reshape(len(x),1)
    logd(f'Regressors {n, m, k}')
    offset = max(n, m + k);
    llen = len(y) - offset;
    logd(debug, f"Length of each array: {llen}, n+m+k={offset}, Length of Original: {len(y)}\n" )
    
    xxx=[]
    for i in range(1, n+1):
        logd(debug, f'==>{offset-i} - {offset-i+llen} {llen}\n' )
        xxx.append(y[offset-i:offset-i+llen].reshape((llen,1)));
        
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

# Theta is the parameter vector with last index is the constant of ARX model
#
def ComputeResid(y, x, theta, n,m,k):
    xx, offset= Regressors(y,x, n,m,k)
    rs1=np.dot(xx, theta[:-1]) + theta[-1]
    return rs1

    
def ARXModelLR(y, x, n,m,k, debug=False):
    xx, offset = Regressors(y, x, n,m,k, debug)
    arx = LinearRegression(fit_intercept=True).fit(xx, y[offset:])
    ret = np.append(arx.coef_, arx.intercept_)
    
    # The following is used until we stabilize the results and then can be removed
    #ret1 = ret.copy();
    #ret1[0:n] *= -1;
    #del xx
    return ret, arx, None;


#Theta is a parameter matrix [n-coeffs, m-coffients, constant]
def predict(x, y, n,m,k, theta, t):
    s = max(n, (m+k))
    if( t < s):
        print("Hmmm Passing an index i:{} resetting to {}".format(i,s))
        t=s
    if (x is None or len(x) < (m+k+1) ):
        #p= list(reversed(y[t-n:t] * -1)) + [0] * (1+m)  + [1]
        p= list(reversed(y[t-n:t])) + [0] * (1+m)  + [1]
    else:
        p= list(reversed(y[t-n:t])) + list(reversed(x[t-m-k:t-k+1])) + [1]

    yh = np.sum(np.array(p).dot(theta))
    rs = (y[t] - yh)
    return yh,rs

# Compute the Fitness Score
#
def FitnessScore(x, y, n,m,k,theta, needArrays=True):
    s=max(n,m+k)
    denom = np.sum((y[s:]- np.mean(y[s:]))**2)
    yhat=np.array(y.copy())
    residueFit=[]
    sumResidue = 0;

    for t in range(s,len(y)):  # <= predict all possible candidates
        yyFit,rrFit = predict(x, yhat, n,m,k, theta, t)
            
        yhat[t] = yyFit
        sumResidue += rrFit ** 2
        
        if (needArrays):
            residueFit.append(rrFit)
   
    fitness = 1- np.sqrt(sumResidue/denom)
    return fitness, yhat,residueFit;

# Compute best fitness score and return the results
#
def findBest(y, x):
    best=Map({});
    x=x.reshape((len(x),1))
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

#This will create a Invariant file. Note the CSV file has the following format
# Time, A, B, C, D => first columns in time and time series for subsequent columns
#
def CreateInvariants(file, outFileName=None, columns_from=0, columns_to=100000):
    df=pd.read_csv(file)
    if ( len(df) <= 2 or len(df.columns) <= 1):
        log.info(f"Not enough data in {file} ... ending")
        raise Exception(f"Not enough data in the dataframe {file}")

    log.info(f"Creating invariants using {df.columns[1:]}" )
    
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
            
            ret = findBest(y, x);
            theta = ",".join([str(c) for c in ret.res])
            corr = np.corrcoef(x,y)[0][1]
            
            inv1 = [u,v,ret.fitscore, corr, theta, ret.n, ret.m, ret.k, ret.threshold]
            log.debug(f"{inv1}" )
            dfi1.loc[len(dfi1)] = inv1
    
    dfi1.sort_values(['uName', 'yName'], inplace=True)
    if ( outFileName is not None):
        dfi1.to_csv(outFileName, index=False)
    
    return dfi1
#
GLOBAL_ARGS=defaultdict(int)
def Usage():
    print('''Usage: sys.argv[0]} csvfile <output file> [from -f columnnumber] [to -t columnnumner]
          Ex: sys.argv[0]}  -f 0 -t 10 test.csv test.inv.0.csv
          ''')
def getargs(opts="hf:t:"): 
    try:
        opts, args = getopt.getopt(sys.argv[1:],opts)
    except getopt.GetoptError:
        Usage("Exception~~")
        
    for opt, arg in opts:
        if opt == '-h': 
            Usage();
        GLOBAL_ARGS[opt] = 1 if not arg else arg;
    GLOBAL_ARGS['__ARGS__'] = args
    
def inJupyter():
    try:
        get_ipython
        return True
    except:
        return False

def main():
    
    global GLOBAL_ARGS
    args = GLOBAL_ARGS['__ARGS__']
    if (len(args) < 2 ): 
        Usage();
        return;
    csvp = args[0]
    outp = args[1]
    cFrom = int(GLOBAL_ARGS['-f']) if ('-f' in GLOBAL_ARGS) else 0
    cTo = int(GLOBAL_ARGS['-t']) if ('-t' in GLOBAL_ARGS) else 100000
    CreateInvariants( csvp, outp, cFrom, cTo)

if __name__ == '__main__':
    if (not inJupyter()):
        t1 = datetime.datetime.now()
        getargs("-hf:t:")
        main()
        t2 = datetime.datetime.now()
        print(f"All Done in {str(t2-t1)} ***")


        
