# matplotlib inline
#import matplotlib.pyplot as plt
#from numpy import *
#from collections import Counter
import numpy as np
#import pylab as pl
#from matplotlib.colors import ListedColormap
#from sklearn import neighbors, datasets, cluster, preprocessing, decomposition
from sklearn import preprocessing
from sklearn.decomposition import PCA
import pandas as pd
#from pandas.tools.plotting import scatter_matrix
#import numpy.random as random
#from mpl_toolkits.mplot3d import Axes3D
#import glob
import os
#import time
from IPython.display import display
#from IPython.display import Image
#from StringIO import StringIO

#import scipy;
#from scipy	import stats;

#import sklearn
from IPython.display import HTML
#import glob;
#import traceback;
import dateutil;
import json;
import urllib;

import matplotlib
matplotlib.style.use('ggplot')


np.set_printoptions(precision=2, linewidth=100)
pd.set_option('display.width', 1000)
pd.set_option('precision',5);
pd.set_option('display.precision', 5)
pd.set_option('mode.sim_interactive', True);
pd.set_option('display.max_rows', 9)
pd.set_option('display.float_format', lambda x: '%.3f' % x)

'''
Loads data set for person comparative analysis
Problems and ISSUES - Check the following:
==========================================
** Does your data have headers! If not you need more complex call
If so, Pass "headers=None" and pass names

'''

def DetermineSeperator(line):
    sep = ","    
    split2  = line.split("\t");
    if (len(split2) > 1 ):
        sep = "\t";   
    return sep;
    
def getAuraDF(link):
    f = urllib.urlopen(link)
    js = f.read()
    if (js.find("$rs=") > 0):
        js = js[js.find("$rs=")+5:]
        #print js[0:100]
        data = json.loads(js)
        df=pd.DataFrame(data['rows'],columns=data['colnames'])
        return df
    else:
        return js
        
def getDF(fileName, debug=False, headers=0, names=None):
    if (    not (fileName.startswith("http://"))  and
            not (fileName.startswith("https://")) and
            not os.path.exists(fileName)):
        #raise Exception( fileName + " does not exist")
        print ("ERROR: *** " +fileName + " does not exist");
        return None;
 
    if fileName.endswith(".xlsx"): 
       df1 = pd.read_excel(fileName)
    elif ("/aura/" in fileName):
        df1 = getAuraDF(fileName);
        return df1
    else:
        sep = ","
        if not fileName.endswith(".csv"):
            with open(fileName, 'r') as f:
                line    = f.readline();    
                #split1  = line.split(",");
                sep = DetermineSeperator(line);
                
        df1 = pd.read_csv(fileName, sep=sep, header=headers, low_memory=False,
                          skipinitialspace =True, names=names, comment='#')
    return df1;

    
def LoadDataSet(fileOrString, columns=None, 
                debug=False, headers=0, names=None, checkForDateTime=False):
    if (fileOrString.find("\n") >=0 ):
        ps = [line.strip() for line in fileOrString.split('\n')
                if line.strip() != '' and not line.startswith("#") ];
        sep = DetermineSeperator(ps[0]);
        ns = [p.split(sep) for p in ps]
        df1 = pd.DataFrame(ns[1:], columns=ns[0]);
    else:               
        df1 = getDF(fileOrString, debug=False, headers=0, names=None)     

    if ( df1 is None or str(type(df1)).find("DataFrame") < 0):
        return df1;
    #df1=df1.convert_objects(convert_numeric=False)
    df2 = df1[columns] if (columns != None ) else df1;

    if (checkForDateTime):
        for i, c in enumerate(df1.columns):
            if (df2.dtypes[i] != object ):
                continue;
            s = df2[c][0]
            if ( len(s) < 8): continue;
            try:
                dateutil.parser.parse(s);
            except:
                continue; 
            print "Trying to convert to datetime:"+ c;
            df2[c] =  pd.to_datetime(df2[c])  
            
    if debug:
        print ("Printing 5 of %d rows", df2.shape[0]);
        print (df2[:5]);
        
    return df2;


def normalizeData(df):
    df1 = df.select_dtypes(exclude=[object])
    vals = df1.values
    cols = df1.columns

    d = preprocessing.scale(vals)
    df2  = pd.DataFrame(d, columns=cols)
    return df2;

def computePCA(d, components = 2, columnNames=None):

    if (d.shape[1] < 2) :
        print ("Number of components are already less than 2")

    pca = PCA(n_components= components)
    pca.fit(d)

    #usorted=  zip(df.columns.values, pca.components_[0])
    #sort = sorted (usorted, key = lambda x:x[1])
    #print (sort)

    return pca
#======================================================================
'''
Evaluating the cost of the clusters
'''
'''
The costOfCluster is the cost assigned to a k-means cluster.
K-means algorithm optimizes on the cost function SUM( ||xi - Mj||^2) - where Mj is the
centroid of of cluster. In this case, if xi ( i is from 1..m) is assigned to cluster j ( j from 1..k)
then the ||xi - Mj||^2) is the square of the distaance
'''
#from scipy.spatial.distance import *
#from sklearn import metrics
#from sklearn.metrics import pairwise_distances

def costOfCluster(kmeans, data) :
    l = kmeans.labels_
    c = kmeans.cluster_centers_
    s = 0.0
    for i, x in enumerate(data):
        c1 = euclidean (x, c[l[i]]);
        s =  s + c1 ** 2
    return s;

# x must be an np array. Thi
def elbowIndex(x):
    if ( x.shape[0] <= 2):
        return 0
    d = 0
    d = x[:-1] - x[1:]
    second_d = np.abs(d[:-1] - d[1:])
    if (second_d.shape[0] <=0):
        return 0;
    print "+++", second_d.shape
    return 1 + np.argmax(second_d)

def DisplayInit():
    style= '''   
        <style>   
        .container { width: 100% !important; }
            div.cell{
                width:100%;
                margin-left:0%;
                margin-right:auto;
            }
        </style>    ''';
    display(HTML(style))
    
def pcolumns(df):
    c1 = [ k for k in df.columns]
    print c1
    return c1
    
DisplayInit()

