import matplotlib.pyplot as plt
from numpy import *
from collections import Counter
import numpy as np
import pylab as pl
import pandas as pd
from pandas.tools.plotting import scatter_matrix
import numpy.random as random
from mpl_toolkits.mplot3d import Axes3D
from sklearn.feature_extraction import DictVectorizer
import time
import datetime
from sklearn.preprocessing import StandardScaler
from sklearn import neighbors, datasets, cluster, preprocessing, decomposition
  
def ConvertCaegoricalToBool(df, n):
    p = df[n];
    
    yy=p.astype(str)
    unique=np.unique(yy)

    if ( len(unique) > 2 ):
        raise Exception("supports only binary - Found (%d) categories!" % len(unique))
    unique.sort();
    y= np.where(yy==unique[0], int(0), int(1) )
    return y;
#
# This function will take the categorical value or any column and
# either replaces inline or create a new column with enumeration of 
# Values for a left column will be changed as shown to the right array:
# For ex. [a,a,a,b,b,c,c,c,c,d, d] => [0,0,0,1,1,2,2,2,2,3, 3]
#
# pd.Categorical(['a', 'a', 'a', 'b', 'b', 'c']).codes
# pd.dummies    

#def encodeCategorical(df, columnName, newColumnName = None, makeCopy = False):
#    print "====> encodeCategorical ", makeCopy;
#    df_mod = df.copy() if makeCopy else df;
#    targets = df_mod[columnName].unique()
#    targets.sort();
#    mapToInt = {name: n for n, name in enumerate(targets)}
#    newColumnName = newColumnName if (newColumnName!=None) else columnName;
#    df_mod.loc[:,newColumnName] = df_mod[columnName].replace(mapToInt)#
#
#   df_mod[columnName].replace(mapToInt, inplace=True)
#
#    return (df_mod, targets, mapToInt)
#
    
#def encodeCategorical(df, columnName = None):
#    l = preprocessing.LabelEncoder()
#    y = df[columnName] if columnName else df;
#    l.fit_transform(dfL[k])
#    return y;

#
# This calls for get_dummies
#
def vectorize(df, cols):
    vec = DictVectorizer(dtype=float32)
    mkdict = lambda row: dict((col, row[col]) for col in cols)
    vecData = pd.DataFrame(vec.fit_transform(df[cols].apply(mkdict, axis=1)).toarray())
    vecData.columns = vec.get_feature_names()
    vecData.index = df.index

    df = df.drop(cols, axis=1)
    df = df.join(vecData)

    return (df, vecData, vec)

def ToUNIXTime(dfi, columnName = None):
    y = dfi if columnName is None else dfi[columnName];
    Lam=lambda k: int((k-datetime.datetime(1970,1,1)).total_seconds())   
    r = y.apply(Lam)

    return r;

def ToYYYYTime(dfi, columnName = None):
    Lam=lambda k: int(datetime64(k).item().strftime("%Y%m%d%H%M%S"))  
    df = dfi if type(dfi) == pd.core.series.Series else dfi[columnName];
    df.fillna(0, inplace=True)
    r = df.apply(Lam)

    return r;

def FromUNIXTime(dfi, columnName):
    l=lambda k:  np.datetime64(int(k), "ms") if (k > 1000000000000) else np.datetime64(int(k), "s")
    #newColumnName = newColumnName if (newColumnName!=None) else columnName;
    #dfi.loc[:,newColumnName] = dfi[columnName].apply(l);
    r = dfi[columnName].apply(l);

    return r

#this counts upto count unique values, if it exceeds the count, it returns None
def myUnique(y, threshhold =50):
    m=set()
    for i in y:
        m.add(i)
        if (len(m) > threshhold):
            return None;
    return m
    
# 
# Category = True - instead suse Vectoriser 
#
def prepareDF(dfi, makeCopy = False, 
              threshHold = 50, fillna=0,
              category=True):
    df = dfi.copy() if makeCopy else dfi;

    df.fillna(fillna, inplace=True)
    t = df.select_dtypes(exclude=[np.number])
    
    if ( len(t.columns) <= 0) :
        return df;

    l = preprocessing.LabelEncoder()
    
    n = df.shape[0]    # number of rows
    for k in t.columns:
        if (len(k) <= 0):
            continue;
        if ( df[k].dtype == 'bool'):
            df[k] = df[k].astype('int');
            continue;
        if (df[k].dtype == np.dtype('datetime64[ns]')):
            #print "YEY";
            #df[k] = ToUNIXTime(df,k);
            df[k] = ToYYYYTime(df,k);
            continue;
        
        
        targets = df[k].unique()
        #targets.sort()
        #mapToInt = {name: n for n, name in enumerate(targets)}
        o = len(targets)
        r = float(o)/n*100;
        if (threshHold == None or threshHold ==0 or (o < threshHold or r < 10)) :
            if ( o == 2 or category ):
                if ( o < threshHold or r < 10 ):
                   try:
                      #print ("processing ", k);
                      df[k] = l.fit_transform(df[k])
                   except:
                      print ("Error while processing {} {} ".format(k,df.columns[k]))
            else:
                (df,v1,v2) = vectorize(df,[k])
    
    t = df.select_dtypes(exclude=[np.number])
    df = df.drop(t, axis=1)
    print ("Dropping in prepareDF - ", t.columns);

    return df

