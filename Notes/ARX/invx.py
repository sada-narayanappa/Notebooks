
from sklearn.linear_model import LinearRegression
import numpy as np;

def logd(debug = True, *args):
    if not debug: return
    for a in args:
        print(a, sep='', end=' ')

# We do not use this function - ARXModelLR is 100 times faster
# Just incase if we need to debug ARX implementation, then we can use this
# 
def ARXModel(y, x, n,m,k):
    if (x is None):
        xx = None
    else:
        xxx=[]
        for i in range(m, -1, -1):
            to = -(m+k-i) if ((m+k-i) > 0) else None
            xxx.append(x[i:to]);
            #print(i, -(m+k-i))

        xx = np.hstack(xxx)
        
    #print ( "===>", len(xx))
    arx = ARX(y[m+k:],xx, lags=n)
    res = arx.fit()
    return res, arx, res;
# We use this function to computer ARX b/w u & v- ARXModelLR is 100 times faster
#
# 
def ARXModelLR(y, x, n,m,k, debug=False):
    logd(f'ARXModelLR {n, m, k}')
    l = max(n, m + k);
    #if (n > 0): l -= 1;
    llen = len(y) - l;
    logd(debug, f"Length of each array: {llen}, n+m+k={l}, Length of Original: {len(y)}\n" )
    
    xxx=[]
    for i in range(1, n+1):
        logd(debug, f'==>{l-i} - {l-i+llen} {llen}\n' )
        xxx.append(y[l-i:l-i+llen].reshape((llen,1)));
        
    for i in range(m, -1, -1):
        upto = len(x) - k
        x1 = x[upto-llen-i: upto-i]
        #x1 = x[i: i+llen]
        #xxx.append(x1);
        
    for i in range(0, m+1):
        upto = len(x) - k
        x1 = x[upto-llen-i: upto-i]
        #x1 = x[i: i+llen]
        xxx.append(x1);
    
    for i in range(len(xxx)):
        logd(debug, "+==>" , len(xxx[i]), "\n")
        
    if (len(xxx) > 0):
        xx = np.hstack(xxx)
        
    arx = LinearRegression(fit_intercept=True).fit(xx, y[l:])
    ret = np.append(arx.coef_, arx.intercept_)
    
    # THe following is used until we stabilize the results and then can be removed
    ret1 = ret.copy();
    ret1[0:n] *= -1;
    return ret, arx, ret1, xx;

u,v,f,c,_,t,n,m,k,_,_,_,etheta=dfi.loc[13]
x=df[u].values    
y=df[v].values
r,_,r1,xx=ARXModelLR (y, x.reshape(len(x),1), n, m, k, True)
print(f'{(n,m,k)}, {r}\n{r1}\n{etheta}')
# We use this function to computer ARX b/w u & v- ARXModelLR is 100 times faster
#
# 
def Regressors(y, x, n,m,k, debug=False):
    logd(f'ARXModelLR {n, m, k}')
    l = max(n, m + k);
    llen = len(y) - l;
    logd(debug, f"Length of each array: {llen}, n+m+k={l}, Length of Original: {len(y)}\n" )
    
    xxx=[]
    for i in range(1, n+1):
        logd(debug, f'==>{l-i} - {l-i+llen} {llen}\n' )
        xxx.append(y[l-i:l-i+llen].reshape((llen,1)));
        
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
        
    return xx;
    
def ARXModelLR(y, x, n,m,k, debug=False):
    Regressors(y, x)
    arx = LinearRegression(fit_intercept=True).fit(xx, y[l:])
    ret = np.append(arx.coef_, arx.intercept_)
    
    # THe following is used until we stabilize the results and then can be removed
    ret1 = ret.copy();
    ret1[0:n] *= -1;
    return ret, arx, ret1, xx;

u,v,f,c,_,t,n,m,k,_,_,_,etheta=dfi.loc[13]
x=df[u].values    
y=df[v].values
r,_,r1,xx=ARXModelLR (y, x.reshape(len(x),1), n, m, k, True)
print(f'{(n,m,k)}, {r}\n{r1}\n{etheta}')
# We use this function to computer ARX b/w u & v- ARXModelLR is 100 times faster
#
# 
def Regressors(y, x, n,m,k, debug=False):
    x = x.reshape(len(x),1)
    logd(f'ARXModelLR {n, m, k}')
    l = max(n, m + k);
    llen = len(y) - l;
    logd(debug, f"Length of each array: {llen}, n+m+k={l}, Length of Original: {len(y)}\n" )
    
    xxx=[]
    for i in range(1, n+1):
        logd(debug, f'==>{l-i} - {l-i+llen} {llen}\n' )
        xxx.append(y[l-i:l-i+llen].reshape((llen,1)));
        
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
        
    return xx;
    
def ARXModelLR(y, x, n,m,k, debug=False):
    Regressors(y, x)
    arx = LinearRegression(fit_intercept=True).fit(xx, y[l:])
    ret = np.append(arx.coef_, arx.intercept_)
    
    # THe following is used until we stabilize the results and then can be removed
    ret1 = ret.copy();
    ret1[0:n] *= -1;
    return ret, arx, ret1, xx;

u,v,f,c,_,t,n,m,k,_,_,_,etheta=dfi.loc[13]
x=df[u].values    
y=df[v].values
r,_,r1,xx=ARXModelLR (y, x, n, m, k, True)
print(f'{(n,m,k)}, {r}\n{r1}\n{etheta}')
# We use this function to computer ARX b/w u & v- ARXModelLR is 100 times faster
#
# 
def Regressors(y, x, n,m,k, debug=False):
    x = x.reshape(len(x),1)
    logd(f'ARXModelLR {n, m, k}')
    l = max(n, m + k);
    llen = len(y) - l;
    logd(debug, f"Length of each array: {llen}, n+m+k={l}, Length of Original: {len(y)}\n" )
    
    xxx=[]
    for i in range(1, n+1):
        logd(debug, f'==>{l-i} - {l-i+llen} {llen}\n' )
        xxx.append(y[l-i:l-i+llen].reshape((llen,1)));
        
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
        
    return xx;
    
def ARXModelLR(y, x, n,m,k, debug=False):
    Regressors(y, x)
    arx = LinearRegression(fit_intercept=True).fit(xx, y[l:])
    ret = np.append(arx.coef_, arx.intercept_)
    
    # THe following is used until we stabilize the results and then can be removed
    ret1 = ret.copy();
    ret1[0:n] *= -1;
    return ret, arx, ret1, xx;

u,v,f,c,_,t,n,m,k,_,_,_,etheta=dfi.loc[13]
x=df[u].values    
y=df[v].values
r,_,r1,xx=ARXModelLR (y, x, n, m, k, True)
print(f'{(n,m,k)}, {r}\n{r1}\n{etheta}')
# We use this function to computer ARX b/w u & v- ARXModelLR is 100 times faster
#
# 
def Regressors(y, x, n,m,k, debug=False):
    x = x.reshape(len(x),1)
    logd(f'ARXModelLR {n, m, k}')
    l = max(n, m + k);
    llen = len(y) - l;
    logd(debug, f"Length of each array: {llen}, n+m+k={l}, Length of Original: {len(y)}\n" )
    
    xxx=[]
    for i in range(1, n+1):
        logd(debug, f'==>{l-i} - {l-i+llen} {llen}\n' )
        xxx.append(y[l-i:l-i+llen].reshape((llen,1)));
        
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
        
    return xx;
    
def ARXModelLR(y, x, n,m,k, debug=False):
    Regressors(y, x)
    arx = LinearRegression(fit_intercept=True).fit(xx, y[l:])
    ret = np.append(arx.coef_, arx.intercept_)
    
    # THe following is used until we stabilize the results and then can be removed
    ret1 = ret.copy();
    ret1[0:n] *= -1;
    return ret, arx, ret1, xx;

u,v,f,c,_,t,n,m,k,_,_,_,etheta=dfi.loc[13]
x=df[u].values    
y=df[v].values
r,_,r1,xx=ARXModelLR (y, x, n, m, k, True)
print(f'{(n,m,k)}, {r}\n{r1}\n{etheta}')
# We use this function to computer ARX b/w u & v- ARXModelLR is 100 times faster
#
# 
def Regressors(y, x, n,m,k, debug=False):
    x = x.reshape(len(x),1)
    logd(f'ARXModelLR {n, m, k}')
    l = max(n, m + k);
    llen = len(y) - l;
    logd(debug, f"Length of each array: {llen}, n+m+k={l}, Length of Original: {len(y)}\n" )
    
    xxx=[]
    for i in range(1, n+1):
        logd(debug, f'==>{l-i} - {l-i+llen} {llen}\n' )
        xxx.append(y[l-i:l-i+llen].reshape((llen,1)));
        
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
        
    return xx;
    
def ARXModelLR(y, x, n,m,k, debug=False):
    Regressors(y, x)
    arx = LinearRegression(fit_intercept=True).fit(xx, y[l:])
    ret = np.append(arx.coef_, arx.intercept_)
    
    # THe following is used until we stabilize the results and then can be removed
    ret1 = ret.copy();
    ret1[0:n] *= -1;
    return ret, arx, ret1, xx;

u,v,f,c,_,t,n,m,k,_,_,_,etheta=dfi.loc[13]
x=df[u].values    
y=df[v].values
r,_,r1,xx=ARXModelLR (y, x, n, m, k, True)
print(f'{(n,m,k)}, {r}\n{r1}\n{etheta}')
# We use this function to computer ARX b/w u & v- ARXModelLR is 100 times faster
#
# 
def Regressors(y, x, n,m,k, debug=False):
    x = x.reshape(len(x),1)
    logd(f'ARXModelLR {n, m, k}')
    l = max(n, m + k);
    llen = len(y) - l;
    logd(debug, f"Length of each array: {llen}, n+m+k={l}, Length of Original: {len(y)}\n" )
    
    xxx=[]
    for i in range(1, n+1):
        logd(debug, f'==>{l-i} - {l-i+llen} {llen}\n' )
        xxx.append(y[l-i:l-i+llen].reshape((llen,1)));
        
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
        
    return xx;
    
def ARXModelLR(y, x, n,m,k, debug=False):
    Regressors(y, x)
    arx = LinearRegression(fit_intercept=True).fit(xx, y[l:])
    ret = np.append(arx.coef_, arx.intercept_)
    
    # THe following is used until we stabilize the results and then can be removed
    ret1 = ret.copy();
    ret1[0:n] *= -1;
    return ret, arx, ret1, xx;

u,v,f,c,_,t,n,m,k,_,_,_,etheta=dfi.loc[13]
x=df[u].values    
y=df[v].values
r,_,r1,xx=ARXModelLR (y, x, n, m, k, True)
print(f'{(n,m,k)}, {r}\n{r1}\n{etheta}')
# We use this function to computer ARX b/w u & v- ARXModelLR is 100 times faster
#
# 
def Regressors(y, x, n,m,k, debug=False):
    x = x.reshape(len(x),1)
    logd(f'ARXModelLR {n, m, k}')
    l = max(n, m + k);
    llen = len(y) - l;
    logd(debug, f"Length of each array: {llen}, n+m+k={l}, Length of Original: {len(y)}\n" )
    
    xxx=[]
    for i in range(1, n+1):
        logd(debug, f'==>{l-i} - {l-i+llen} {llen}\n' )
        xxx.append(y[l-i:l-i+llen].reshape((llen,1)));
        
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
        
    return xx;
    
def ARXModelLR(y, x, n,m,k, debug=False):
    Regressors(y, x)
    arx = LinearRegression(fit_intercept=True).fit(xx, y[l:])
    ret = np.append(arx.coef_, arx.intercept_)
    
    # THe following is used until we stabilize the results and then can be removed
    ret1 = ret.copy();
    ret1[0:n] *= -1;
    return ret, arx, ret1, xx;

u,v,f,c,_,t,n,m,k,_,_,_,etheta=dfi.loc[13]
x=df[u].values    
y=df[v].values
r,_,r1,xx=ARXModelLR (y, x, n, m, k, True)
print(f'{(n,m,k)}, {r}\n{r1}\n{etheta}')
# We use this function to computer ARX b/w u & v- ARXModelLR is 100 times faster
#
# 
def Regressors(y, x, n,m,k, debug=False):
    x = x.reshape(len(x),1)
    logd(f'ARXModelLR {n, m, k}')
    l = max(n, m + k);
    llen = len(y) - l;
    logd(debug, f"Length of each array: {llen}, n+m+k={l}, Length of Original: {len(y)}\n" )
    
    xxx=[]
    for i in range(1, n+1):
        logd(debug, f'==>{l-i} - {l-i+llen} {llen}\n' )
        xxx.append(y[l-i:l-i+llen].reshape((llen,1)));
        
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
        
    return xx;
    
def ARXModelLR(y, x, n,m,k, debug=False):
    Regressors(y, x)
    arx = LinearRegression(fit_intercept=True).fit(xx, y[l:])
    ret = np.append(arx.coef_, arx.intercept_)
    
    # THe following is used until we stabilize the results and then can be removed
    ret1 = ret.copy();
    ret1[0:n] *= -1;
    return ret, arx, ret1, xx;

u,v,f,c,_,t,n,m,k,_,_,_,etheta=dfi.loc[13]
x=df[u].values    
y=df[v].values
r,_,r1,xx=ARXModelLR (y, x, n, m, k, True)
print(f'{(n,m,k)}, {r}\n{r1}\n{etheta}')